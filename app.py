from __future__ import annotations

import io
from datetime import datetime
from html import escape

import pandas as pd
import streamlit as st

from media_mixer.calculations import (
    Component,
    MixResult,
    density_correction_add_stock,
    fast_dirty_percoll,
    ips74_osmotic,
    ips_basic,
    optiprep_endpoint_fibrinogen,
    optiprep_table_density,
    optiprep_working_solution,
    rbc_density_gradient_targets,
    sdm_ips_pbs,
    sdm_with_whole_blood,
    ternary_ips_optiprep_pbs,
    two_endpoint_density_mix,
)

st.set_page_config(page_title="Percoll / OptiPrep Density Media Mixer", layout="wide")

DEFAULTS = {
    "rho_percoll": 1.1300,
    "rho_10xpbs": 1.0674,
    "rho_pbs": 1.0047,
    "rho_ips": 1.1116,
    "rho_ips74": 1.1189,
    "rho_plasma": 1.0255,
    "rho_optiprep_raw": 1.320,
    "rho_optiprep_40": 1.215,
    "rho_hcl": 1.000,
    "rho_fib_stock": 1.020,
}


def fmt(x, nd=6):
    if x is None:
        return ""
    try:
        return f"{x:.{nd}f}"
    except Exception:
        return str(x)


def result_to_df(result: MixResult) -> pd.DataFrame:
    df = pd.DataFrame(result.as_rows())
    for col in ["Volume [ml]", "Density [g/ml]", "Mass [g]", "v/v fraction"]:
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def nice_tick_step(volume_ml: float) -> float:
    """Choose a readable major tick spacing in ml."""
    candidates = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 25, 50, 100, 200, 500, 1000]
    if volume_ml <= 0:
        return 1.0
    target_marks = 6
    ideal = volume_ml / target_marks
    for step in candidates:
        if step >= ideal:
            return step
    return candidates[-1]


def tube_svg(result: MixResult, width: int = 440, height: int = 640) -> str:
    """Create a simple SVG tube with stacked component fills and volume marks."""
    total = max(result.target_volume_ml, 1e-9)
    colors = [
        '#c94f4f', '#4f8cc9', '#67b567', '#d4a64a', '#9b6fcf', '#53b7b0',
        '#e07f39', '#8d9aa5', '#de5d9e', '#62a04a'
    ]
    left_margin = 86
    tube_x = 115
    tube_y = 30
    tube_w = 120
    tube_h = 500
    bottom_r = tube_w / 2
    body_h = tube_h - bottom_r
    right_legend_x = 275
    tick_len_major = 14
    tick_len_minor = 8
    label_fs = 15
    title_fs = 20

    def y_for_volume(v: float) -> float:
        frac = min(max(v / total, 0.0), 1.0)
        return tube_y + tube_h - frac * tube_h

    clip_id = f"tubeclip_{abs(hash((result.title, result.target_volume_ml, len(result.components)))) % 10**8}"

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append(f'<text x="{width/2}" y="22" text-anchor="middle" font-size="{title_fs}" font-family="Arial, sans-serif" font-weight="bold">{escape(result.title)}</text>')
    svg.append(f'<defs><clipPath id="{clip_id}"><rect x="{tube_x}" y="{tube_y}" width="{tube_w}" height="{tube_h}" rx="{tube_w/2}" ry="{tube_w/2}"/></clipPath></defs>')

    # Tube background
    svg.append(f'<rect x="{tube_x}" y="{tube_y}" width="{tube_w}" height="{tube_h}" rx="{tube_w/2}" ry="{tube_w/2}" fill="#f8fafc" stroke="#455a64" stroke-width="4"/>')

    # Fills
    running = 0.0
    for idx, comp in enumerate(result.components):
        if comp.volume_ml <= 0:
            continue
        y0 = y_for_volume(running)
        running += comp.volume_ml
        y1 = y_for_volume(running)
        section_y = y1
        section_h = max(y0 - y1, 0)
        color = colors[idx % len(colors)]
        svg.append(f'<rect x="{tube_x}" y="{section_y}" width="{tube_w}" height="{section_h}" fill="{color}" clip-path="url(#{clip_id})"/>')

    # Top meniscus shine
    svg.append(f'<rect x="{tube_x + 14}" y="{tube_y + 18}" width="18" height="{tube_h - 36}" rx="9" fill="rgba(255,255,255,0.32)" clip-path="url(#{clip_id})"/>')

    # Component boundary lines
    running = 0.0
    for idx, comp in enumerate(result.components[:-1]):
        running += comp.volume_ml
        y = y_for_volume(running)
        if tube_y < y < tube_y + tube_h:
            svg.append(f'<line x1="{tube_x+6}" x2="{tube_x+tube_w-6}" y1="{y}" y2="{y}" stroke="rgba(255,255,255,0.7)" stroke-width="2" clip-path="url(#{clip_id})"/>')

    # Tube outline again for crisp edge
    svg.append(f'<rect x="{tube_x}" y="{tube_y}" width="{tube_w}" height="{tube_h}" rx="{tube_w/2}" ry="{tube_w/2}" fill="none" stroke="#37474f" stroke-width="4"/>')

    # Major and minor volume marks
    major = nice_tick_step(total)
    minor = major / 2 if major >= 0.2 else major
    n_minor = int(total / minor)
    for i in range(n_minor + 1):
        v = i * minor
        if v > total + 1e-9:
            break
        y = y_for_volume(v)
        is_major = abs((v / major) - round(v / major)) < 1e-9
        tick_len = tick_len_major if is_major else tick_len_minor
        stroke = '#607d8b' if is_major else '#b0bec5'
        svg.append(f'<line x1="{tube_x - tick_len}" x2="{tube_x}" y1="{y}" y2="{y}" stroke="{stroke}" stroke-width="2"/>')
        if is_major:
            label = f"{v:.0f}" if major >= 1 else f"{v:.1f}"
            svg.append(f'<text x="{tube_x - tick_len - 6}" y="{y + 5}" text-anchor="end" font-size="{label_fs}" font-family="Arial, sans-serif" fill="#37474f">{label}</text>')
    svg.append(f'<text x="{left_margin/2 + 10}" y="{tube_y + tube_h/2}" transform="rotate(-90 {left_margin/2 + 10} {tube_y + tube_h/2})" text-anchor="middle" font-size="16" font-family="Arial, sans-serif" fill="#37474f">volume [ml]</text>')

    # Total label
    svg.append(f'<text x="{tube_x + tube_w/2}" y="{tube_y + tube_h + 32}" text-anchor="middle" font-size="16" font-family="Arial, sans-serif" fill="#37474f">total = {total:.4g} ml</text>')

    # Legend
    svg.append(f'<text x="{right_legend_x}" y="{tube_y + 8}" font-size="17" font-family="Arial, sans-serif" font-weight="bold" fill="#263238">Mixture</text>')
    legend_y = tube_y + 34
    for idx, comp in enumerate(result.components):
        color = colors[idx % len(colors)]
        frac_pct = 100 * comp.volume_ml / total if total > 0 else 0
        label = f"{comp.name}: {comp.volume_ml:.4g} ml ({frac_pct:.1f}%)"
        svg.append(f'<rect x="{right_legend_x}" y="{legend_y - 12}" width="16" height="16" fill="{color}" stroke="#455a64" stroke-width="1"/>')
        svg.append(f'<text x="{right_legend_x + 24}" y="{legend_y}" font-size="14" font-family="Arial, sans-serif" fill="#263238">{escape(label)}</text>')
        legend_y += 24
    if result.apparent_density_g_ml is not None:
        density_txt = f"density = {result.apparent_density_g_ml:.6f} g/ml"
        svg.append(f'<text x="{right_legend_x}" y="{legend_y + 8}" font-size="14" font-family="Arial, sans-serif" fill="#263238">{escape(density_txt)}</text>')

    svg.append('</svg>')
    return ''.join(svg)


def render_tube_visualization(result: MixResult):
    st.subheader("Tube visualization")
    st.caption("Stacked by component volume. Volume marks use the final/target total volume.")
    svg = tube_svg(result)
    st.image(svg.encode('utf-8'))
    st.download_button(
        "Download tube SVG",
        svg.encode('utf-8'),
        file_name=f"{safe_name(result.title)}_tube_{datetime.now().strftime('%Y%m%d_%H%M')}.svg",
        mime='image/svg+xml',
    )


def render_result(result: MixResult):
    st.subheader(result.title)
    c1, c2, c3 = st.columns(3)
    c1.metric("Final / target volume", f"{result.target_volume_ml:.6g} ml")
    if result.apparent_density_g_ml is not None:
        c2.metric("Calculated density", f"{result.apparent_density_g_ml:.6f} g/ml")
    c3.metric("Number of components", str(len(result.components)))

    viz_col, table_col = st.columns([1, 1.3])
    with viz_col:
        render_tube_visualization(result)

    df = result_to_df(result)
    with table_col:
        st.subheader("Component table")
        st.dataframe(
            df.style.format({
                "Volume [ml]": "{:.6f}",
                "Density [g/ml]": "{:.6f}",
                "Mass [g]": "{:.6f}",
                "v/v fraction": "{:.6f}",
            }, na_rep=""),
            use_container_width=True,
        )

    with st.expander("Calculation checks / internal fractions", expanded=False):
        if result.checks:
            check_df = pd.DataFrame([{"Quantity": k, "Value": v} for k, v in result.checks.items()])
            st.dataframe(check_df, use_container_width=True)
        else:
            st.write("No additional checks.")

    if result.notes:
        st.info("\n".join(f"• {n}" for n in result.notes))

    protocol = make_protocol_text(result)
    st.download_button(
        "Download protocol as .txt",
        protocol,
        file_name=f"{safe_name(result.title)}_{datetime.now().strftime('%Y%m%d_%H%M')}.txt",
        mime="text/plain",
    )
    st.download_button(
        "Download table as CSV",
        df.to_csv(index=False).encode("utf-8"),
        file_name=f"{safe_name(result.title)}_{datetime.now().strftime('%Y%m%d_%H%M')}.csv",
        mime="text/csv",
    )


def safe_name(name: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in name).strip("_")


def make_protocol_text(result: MixResult) -> str:
    lines = []
    lines.append(result.title)
    lines.append("=" * len(result.title))
    lines.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Target/final volume: {result.target_volume_ml:.10g} ml")
    if result.apparent_density_g_ml is not None:
        lines.append(f"Calculated density: {result.apparent_density_g_ml:.10g} g/ml")
    lines.append("")
    lines.append("Components")
    lines.append("----------")
    for c in result.components:
        m = "" if c.mass_g is None else f", mass {c.mass_g:.10g} g"
        d = "" if c.density_g_ml is None else f", density {c.density_g_ml:.10g} g/ml"
        f = "" if c.fraction_vv is None else f", v/v fraction {c.fraction_vv:.10g}"
        note = "" if not c.note else f" [{c.note}]"
        lines.append(f"- {c.name}: {c.volume_ml:.10g} ml{d}{m}{f}{note}")
    if result.checks:
        lines.append("")
        lines.append("Checks")
        lines.append("------")
        for k, v in result.checks.items():
            lines.append(f"- {k}: {v:.10g}")
    if result.notes:
        lines.append("")
        lines.append("Notes")
        lines.append("-----")
        lines.extend(f"- {n}" for n in result.notes)
    return "\n".join(lines) + "\n"


def common_density_inputs():
    st.caption("Use measured densities at the working temperature whenever possible.")
    c1, c2, c3 = st.columns(3)
    with c1:
        rho_ips = st.number_input("ρ IPS / dense Percoll stock [g/ml]", 0.5, 2.0, DEFAULTS["rho_ips"], 0.0001, format="%.6f")
    with c2:
        rho_pbs = st.number_input("ρ 1x PBS [g/ml]", 0.5, 2.0, DEFAULTS["rho_pbs"], 0.0001, format="%.6f")
    with c3:
        rho_opt = st.number_input("ρ OptiPrep working solution [g/ml]", 0.5, 2.0, DEFAULTS["rho_optiprep_40"], 0.0001, format="%.6f")
    return rho_ips, rho_pbs, rho_opt


st.title("Percoll / OptiPrep Density Media Mixer")
st.write(
    "A Streamlit calculator for RBC density media: IPS, IPS74, SDM, blood-aware SDM, "
    "OptiPrep/fibrinogen endpoints, ternary Percoll–OptiPrep–PBS mixes, and density corrections."
)
st.warning(
    "Bench-use warning: calculations assume additive volumes and linear density mixing. "
    "For RBC work, verify final pH, osmolarity, and density experimentally when precision matters.",
    icon="⚠️",
)

with st.sidebar:
    st.header("Global defaults")
    st.write("These are starting values from your protocol collection and data sheet. Override per tab when needed.")
    st.code(
        "\n".join(f"{k} = {v}" for k, v in DEFAULTS.items()),
        language="text",
    )
    st.markdown("---")
    st.markdown("**Assumptions**")
    st.markdown("- Volume-weighted density mixing")
    st.markdown("- Mass = density × volume")
    st.markdown("- Blood mode targets suspending liquid density, not whole suspension density")


tabs = st.tabs([
    "Fast Percoll",
    "IPS / IPS74",
    "SDM IPS+PBS",
    "SDM + blood",
    "OptiPrep",
    "OptiPrep + fibrinogen",
    "Percoll + OptiPrep",
    "Density correction",
    "Gradient helper",
    "Formula reference",
])

with tabs[0]:
    st.header("Fast approximate Percoll protocol")
    st.write("Legacy formula: 0.9 V Percoll + 0.1 V 10×PBS + 0.11 V 1×PBS.")
    exact = st.radio("Interpret entered volume as", ["Exact final volume", "Legacy base V"], horizontal=True) == "Exact final volume"
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        V = st.number_input("Volume [ml]", 0.001, 10000.0, 50.0, 1.0, key="fast_V")
    with c2:
        rho_pc = st.number_input("ρ Percoll [g/ml]", 0.5, 2.0, DEFAULTS["rho_percoll"], 0.0001, format="%.6f", key="fast_pc")
    with c3:
        rho_10 = st.number_input("ρ 10×PBS [g/ml]", 0.5, 2.0, 1.0400, 0.0001, format="%.6f", key="fast_10")
    with c4:
        rho_1 = st.number_input("ρ 1×PBS [g/ml]", 0.5, 2.0, DEFAULTS["rho_pbs"], 0.0001, format="%.6f", key="fast_1")
    try:
        render_result(fast_dirty_percoll(V, exact_final_volume=exact, rho_percoll=rho_pc, rho_pbs10=rho_10, rho_pbs1=rho_1))
    except Exception as e:
        st.error(str(e))

with tabs[1]:
    st.header("IPS and IPS74 preparation")
    mode = st.radio("Recipe", ["IPS empirical", "IPS 90:10", "IPS74 empirical", "IPS74 osmotic calculation"], horizontal=False)
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        V = st.number_input("Target volume [ml]", 0.001, 10000.0, 100.0, 1.0, key="ips_V")
    with c2:
        rho_pc = st.number_input("ρ Percoll [g/ml]", 0.5, 2.0, DEFAULTS["rho_percoll"], 0.0001, format="%.6f", key="ips_pc")
    with c3:
        rho_10 = st.number_input("ρ 10×PBS [g/ml]", 0.5, 2.0, DEFAULTS["rho_10xpbs"], 0.0001, format="%.6f", key="ips_10")
    with c4:
        rho_hcl = st.number_input("ρ HCl solution [g/ml]", 0.5, 2.0, DEFAULTS["rho_hcl"], 0.0001, format="%.6f", key="ips_hcl_rho")

    try:
        if mode == "IPS empirical":
            render_result(ips_basic(V, mode="empirical", rho_percoll=rho_pc, rho_pbs10=rho_10))
        elif mode == "IPS 90:10":
            render_result(ips_basic(V, mode="simple_90_10", rho_percoll=rho_pc, rho_pbs10=rho_10))
        elif mode == "IPS74 empirical":
            render_result(ips74_osmotic(V, rho_percoll=rho_pc, rho_pbs10=rho_10, rho_hcl=rho_hcl, empirical_fixed=True))
        else:
            st.subheader("Osmolarity inputs")
            c1, c2, c3, c4, c5 = st.columns(5)
            with c1:
                osm_target = st.number_input("target osm [mOsm/kg]", 0.0, 5000.0, 300.0, 1.0)
            with c2:
                osm_pc = st.number_input("Percoll osm", 0.0, 5000.0, 13.0, 1.0)
            with c3:
                osm_10 = st.number_input("10×PBS osm", 0.0, 10000.0, 3032.0, 1.0)
            with c4:
                osm_hcl = st.number_input("HCl osm", 0.0, 10000.0, 355.0, 1.0)
            with c5:
                p_hcl = st.number_input("HCl fraction v/v", 0.0, 1.0, 0.032, 0.001, format="%.6f")
            render_result(ips74_osmotic(V, osm_target=osm_target, osm_pbs10=osm_10, osm_hcl=osm_hcl,
                                        osm_percoll=osm_pc, p_hcl=p_hcl, rho_percoll=rho_pc,
                                        rho_pbs10=rho_10, rho_hcl=rho_hcl, empirical_fixed=False))
    except Exception as e:
        st.error(str(e))

with tabs[2]:
    st.header("Specific density medium: IPS + PBS")
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        V = st.number_input("Final volume [ml]", 0.001, 10000.0, 20.0, 1.0, key="sdm_V")
    with c2:
        rho_target = st.number_input("Target density [g/ml]", 0.5, 2.0, 1.1000, 0.0001, format="%.6f", key="sdm_rho")
    with c3:
        rho_ips = st.number_input("ρ IPS [g/ml]", 0.5, 2.0, DEFAULTS["rho_ips"], 0.0001, format="%.6f", key="sdm_ips")
    with c4:
        rho_pbs = st.number_input("ρ PBS [g/ml]", 0.5, 2.0, DEFAULTS["rho_pbs"], 0.0001, format="%.6f", key="sdm_pbs")
    try:
        render_result(sdm_ips_pbs(V, rho_target, rho_ips, rho_pbs))
    except Exception as e:
        st.error(str(e))

with tabs[3]:
    st.header("SDM with whole blood addition")
    st.write("Calculates IPS/PBS carrier so the final suspending liquid phase reaches the target density after adding whole blood.")
    c1, c2, c3 = st.columns(3)
    with c1:
        V = st.number_input("Final total volume incl. RBCs [ml]", 0.001, 10000.0, 50.0, 1.0, key="blood_V")
        target_rbc = st.number_input("Final RBC volume fraction", 0.0, 0.95, 0.05, 0.005, format="%.6f", key="blood_rbc")
    with c2:
        rho_target = st.number_input("Target liquid density [g/ml]", 0.5, 2.0, 1.1000, 0.0001, format="%.6f", key="blood_rho")
        hct = st.number_input("Sample hematocrit", 0.001, 0.99, 0.45, 0.01, format="%.6f", key="blood_hct")
    with c3:
        rho_ips = st.number_input("ρ IPS [g/ml]", 0.5, 2.0, DEFAULTS["rho_ips"], 0.0001, format="%.6f", key="blood_ips")
        rho_pbs = st.number_input("ρ PBS [g/ml]", 0.5, 2.0, DEFAULTS["rho_pbs"], 0.0001, format="%.6f", key="blood_pbs")
        rho_plasma = st.number_input("ρ plasma / sample liquid [g/ml]", 0.5, 2.0, DEFAULTS["rho_plasma"], 0.0001, format="%.6f", key="blood_plasma")
    try:
        render_result(sdm_with_whole_blood(V, rho_target, rho_ips, rho_pbs, target_rbc, hct, rho_sample_liquid=rho_plasma))
    except Exception as e:
        st.error(str(e))

with tabs[4]:
    st.header("OptiPrep working solution and table density")
    sub = st.radio("Task", ["Make working solution", "Look up/interpolate data-sheet density"], horizontal=True)
    if sub == "Make working solution":
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            V = st.number_input("Final WS volume [ml]", 0.001, 10000.0, 30.0, 1.0, key="op_ws_V")
        with c2:
            target_pct = st.number_input("Target iodixanol [% w/v]", 0.001, 60.0, 40.0, 1.0, key="op_ws_pct")
        with c3:
            stock_pct = st.number_input("Raw stock iodixanol [% w/v]", 0.001, 100.0, 60.0, 1.0, key="op_ws_stockpct")
            rho_stock = st.number_input("ρ raw OptiPrep [g/ml]", 0.5, 2.0, DEFAULTS["rho_optiprep_raw"], 0.0001, format="%.6f", key="op_ws_rhostock")
        with c4:
            rho_diluent = st.number_input("ρ diluent [g/ml]", 0.5, 2.0, 1.006, 0.0001, format="%.6f", key="op_ws_rhodil")
        try:
            render_result(optiprep_working_solution(V, target_pct, stock_pct, rho_stock, rho_diluent))
        except Exception as e:
            st.error(str(e))
    else:
        c1, c2 = st.columns(2)
        with c1:
            pct = st.number_input("Iodixanol [% w/v]", 6.0, 40.0, 18.0, 0.5, key="op_table_pct")
        with c2:
            medium = st.selectbox("Data-sheet medium", ["NaCl/Tricine", "DMEM"])
        try:
            rho, method = optiprep_table_density(pct, medium)
            st.metric("Density", f"{rho:.6f} g/ml")
            st.caption(method)
            v_ws = pct / 40.0 * 4.0
            v_medium = 4.0 - v_ws
            render_result(MixResult(
                title=f"Data-sheet OptiPrep table mix ({medium})",
                target_volume_ml=4.0,
                apparent_density_g_ml=rho,
                components=[
                    Component.from_volume("40% iodixanol working solution", v_ws, DEFAULTS["rho_optiprep_40"], 4.0),
                    Component.from_volume(medium, v_medium, 1.006 if medium == "NaCl/Tricine" else 1.0065, 4.0),
                ],
                notes=["Density shown is from data-sheet table/interpolation, not from linear density mixing."],
                checks={"iodixanol_percent_wv": pct, "table_density_g_ml": rho},
            ))
        except Exception as e:
            st.error(str(e))

with tabs[5]:
    st.header("OptiPrep + fibrinogen endpoint and endpoint-to-density mixing")
    sub = st.radio("Task", ["Prepare one fibrinogen endpoint", "Mix high/low endpoints to target density"], horizontal=True)
    if sub == "Prepare one fibrinogen endpoint":
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            V = st.number_input("Endpoint volume [ml]", 0.001, 10000.0, 10.0, 1.0, key="fib_endpoint_V")
        with c2:
            target_fib = st.number_input("Target fibrinogen [mg/ml]", 0.0, 200.0, 5.0, 0.5, key="fib_endpoint_target")
            stock_fib = st.number_input("Fibrinogen stock [mg/ml]", 0.001, 500.0, 50.0, 1.0, key="fib_endpoint_stock")
        with c3:
            base_name = st.selectbox("Base medium", ["OptiPrep 40% working solution", "PBS", "Custom medium"])
            rho_base_default = DEFAULTS["rho_optiprep_40"] if "OptiPrep" in base_name else DEFAULTS["rho_pbs"]
            rho_base = st.number_input("ρ base [g/ml]", 0.5, 2.0, rho_base_default, 0.0001, format="%.6f", key="fib_endpoint_rhobase")
        with c4:
            rho_fib = st.number_input("ρ fibrinogen stock [g/ml]", 0.5, 2.0, DEFAULTS["rho_fib_stock"], 0.0001, format="%.6f", key="fib_endpoint_rhofib")
        try:
            render_result(optiprep_endpoint_fibrinogen(V, target_fib, stock_fib, base_name, rho_base, rho_fib))
        except Exception as e:
            st.error(str(e))
    else:
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            V = st.number_input("Final volume [ml]", 0.001, 10000.0, 20.0, 1.0, key="endpoint_mix_V")
        with c2:
            rho_target = st.number_input("Target density [g/ml]", 0.5, 2.0, 1.1000, 0.0001, format="%.6f", key="endpoint_mix_target")
        with c3:
            rho_low = st.number_input("ρ low endpoint [g/ml]", 0.5, 2.0, 1.0100, 0.0001, format="%.6f", key="endpoint_mix_low")
        with c4:
            rho_high = st.number_input("ρ high endpoint [g/ml]", 0.5, 2.0, 1.2000, 0.0001, format="%.6f", key="endpoint_mix_high")
        try:
            render_result(two_endpoint_density_mix(V, rho_target, rho_low, rho_high,
                                                   low_name="PBS + fibrinogen endpoint",
                                                   high_name="OptiPrep + fibrinogen endpoint"))
        except Exception as e:
            st.error(str(e))

with tabs[6]:
    st.header("Ternary Percoll/IPS + OptiPrep + PBS")
    st.write("Fix the Percoll/IPS fraction, then solve OptiPrep and PBS fractions for the target density.")
    c1, c2, c3 = st.columns(3)
    with c1:
        V = st.number_input("Final volume [ml]", 0.001, 10000.0, 10.0, 1.0, key="ternary_V")
        rho_target = st.number_input("Target density [g/ml]", 0.5, 2.0, 1.1000, 0.0001, format="%.6f", key="ternary_rho")
    with c2:
        kind = st.selectbox("Fixed input means", ["IPS fraction", "Percoll-equivalent fraction"])
        fixed = st.number_input("Fixed fraction v/v", 0.0, 1.0, 0.45, 0.01, format="%.6f", key="ternary_fixed")
        percoll_in_ips = st.number_input("Percoll fraction inside IPS", 0.001, 1.0, 0.9, 0.01, format="%.6f", key="ternary_pcips")
    with c3:
        rho_ips = st.number_input("ρ IPS [g/ml]", 0.5, 2.0, DEFAULTS["rho_ips"], 0.0001, format="%.6f", key="ternary_ips")
        rho_opt = st.number_input("ρ OptiPrep WS [g/ml]", 0.5, 2.0, DEFAULTS["rho_optiprep_40"], 0.0001, format="%.6f", key="ternary_opt")
        rho_pbs = st.number_input("ρ PBS [g/ml]", 0.5, 2.0, DEFAULTS["rho_pbs"], 0.0001, format="%.6f", key="ternary_pbs")
    try:
        render_result(ternary_ips_optiprep_pbs(V, rho_target, fixed, rho_ips, rho_opt, rho_pbs, kind, percoll_in_ips))
    except Exception as e:
        st.error(str(e))

with tabs[7]:
    st.header("Density correction")
    st.write("Add one dense or light correction stock to an existing medium.")
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        V = st.number_input("Existing volume [ml]", 0.001, 10000.0, 10.0, 1.0, key="corr_V")
    with c2:
        rho_existing = st.number_input("Existing density [g/ml]", 0.5, 2.0, 1.0950, 0.0001, format="%.6f", key="corr_existing")
    with c3:
        rho_target = st.number_input("Target density [g/ml]", 0.5, 2.0, 1.1000, 0.0001, format="%.6f", key="corr_target")
    with c4:
        rho_stock = st.number_input("Correction stock density [g/ml]", 0.5, 2.0, DEFAULTS["rho_ips"], 0.0001, format="%.6f", key="corr_stock")
    try:
        render_result(density_correction_add_stock(V, rho_existing, rho_target, rho_stock, "Correction stock"))
    except Exception as e:
        st.error(str(e))

with tabs[8]:
    st.header("Gradient helper")
    st.write("Choose low/high endpoint densities around an average RBC density, then prepare each endpoint using one of the other tabs.")
    c1, c2 = st.columns(2)
    with c1:
        rho_rbc = st.number_input("Average RBC density [g/ml]", 0.5, 2.0, 1.0960, 0.0001, format="%.6f")
    with c2:
        half_width = st.number_input("Half-width Δρ [g/ml]", 0.0, 0.1, 0.0050, 0.0005, format="%.6f")
    targets = rbc_density_gradient_targets(rho_rbc, half_width)
    st.metric("Low endpoint density", f"{targets['rho_low']:.6f} g/ml")
    st.metric("High endpoint density", f"{targets['rho_high']:.6f} g/ml")
    st.info("For a prepared linear gradient, make low and high endpoint media with matched additive concentrations, then pump with complementary flow profiles.")

with tabs[9]:
    st.header("Formula reference")
    st.markdown(r"""
### Two-component density mixing

$$p_H = \frac{\rho_T-\rho_L}{\rho_H-\rho_L}, \quad p_L=1-p_H$$

### SDM from IPS and PBS

$$p_{IPS}=\frac{\rho_{SDM}-\rho_{PBS}}{\rho_{IPS}-\rho_{PBS}}, \quad p_{PBS}=1-p_{IPS}$$

### Whole blood addition

$$V_{RBC}=p_{RBC}V, \quad V_{blood}=\frac{V_{RBC}}{hct}, \quad V_{plasma}=V_{blood}-V_{RBC}$$

The target density is applied to the liquid phase:

$$f_{plasma}=\frac{V_{plasma}}{V-V_{RBC}}$$

$$p_{IPS}=\frac{\rho_T-\rho_{PBS}-f_{plasma}(\rho_{plasma}-\rho_{PBS})}{\rho_{IPS}-\rho_{PBS}}$$

$$p_{PBS}=1-p_{IPS}-f_{plasma}$$

### Ternary IPS + OptiPrep + PBS

For fixed IPS fraction \(P\):

$$O=\frac{\rho_T-P\rho_{IPS}-(1-P)\rho_{PBS}}{\rho_{OptiPrep}-\rho_{PBS}}, \quad B=1-P-O$$

### IPS74 osmolarity calculation

$$p_{PC}=\frac{osm_T-osm_{10PBS}-p_{HCl}(osm_{HCl}-osm_{10PBS})}{osm_{PC}-osm_{10PBS}}$$

$$p_{10PBS}=1-p_{PC}-p_{HCl}$$
""")
    st.caption("All equations use volume fractions unless stated otherwise.")
