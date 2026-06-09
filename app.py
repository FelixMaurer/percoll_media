from __future__ import annotations

import io
import re
from datetime import datetime
from html import escape

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

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
    paired_gradient_endpoint_batch,
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


def _svg_tspans(lines, x, y, *, font_size=14, line_height=17, fill="#263238", anchor="start"):
    """Return SVG <text> with multiple tspans for wrapped labels."""
    escaped_lines = [escape(line) for line in lines]
    tspans = []
    for i, line in enumerate(escaped_lines):
        dy = 0 if i == 0 else line_height
        tspans.append(f'<tspan x="{x}" dy="{dy}">{line}</tspan>')
    return (
        f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{font_size}" '
        f'font-family="Arial, sans-serif" fill="{fill}">' + "".join(tspans) + "</text>"
    )


def _wrap_label(label: str, max_chars: int = 48) -> list[str]:
    """Simple word wrap for SVG legend labels."""
    words = label.split()
    if not words:
        return [""]
    lines = []
    current = words[0]
    for word in words[1:]:
        if len(current) + 1 + len(word) <= max_chars:
            current += " " + word
        else:
            lines.append(current)
            current = word
    lines.append(current)
    return lines


def tube_svg(result: MixResult, width: int = 760, height: int | None = None) -> str:
    """Create a responsive SVG tube with stacked component fills and volume marks.

    The legend area is deliberately wide and wraps labels so component names do
    not get clipped in Streamlit.
    """
    total = max(result.target_volume_ml, 1e-9)
    colors = [
        '#c94f4f', '#4f8cc9', '#67b567', '#d4a64a', '#9b6fcf', '#53b7b0',
        '#e07f39', '#8d9aa5', '#de5d9e', '#62a04a'
    ]

    # Layout
    tube_x = 130
    tube_y = 46
    tube_w = 126
    tube_h = 500
    legend_x = 320
    legend_text_x = legend_x + 26
    legend_max_chars = 46
    tick_len_major = 14
    tick_len_minor = 8
    label_fs = 15
    title_fs = 20

    # Precompute wrapped legend lines and needed height.
    legend_blocks = []
    for comp in result.components:
        frac_pct = 100 * comp.volume_ml / total if total > 0 else 0
        label = f"{comp.name}: {comp.volume_ml:.4g} ml ({frac_pct:.1f}%)"
        legend_blocks.append(_wrap_label(label, legend_max_chars))

    legend_line_height = 17
    legend_block_gap = 9
    legend_height = 42 + sum(len(lines) * legend_line_height + legend_block_gap for lines in legend_blocks) + 28
    if result.apparent_density_g_ml is not None:
        legend_height += 22
    if height is None:
        height = int(max(640, tube_y + tube_h + 76, tube_y + legend_height + 40))

    def y_for_volume(v: float) -> float:
        frac = min(max(v / total, 0.0), 1.0)
        return tube_y + tube_h - frac * tube_h

    clip_id = f"tubeclip_{abs(hash((result.title, result.target_volume_ml, len(result.components)))) % 10**8}"

    svg = []
    # width="100%" makes it scale to the available Streamlit column, but the
    # viewBox preserves the wider legend canvas.
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="100%" height="{height}" viewBox="0 0 {width} {height}" preserveAspectRatio="xMidYMin meet">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append(f'<text x="{width/2}" y="24" text-anchor="middle" font-size="{title_fs}" font-family="Arial, sans-serif" font-weight="bold">{escape(result.title)}</text>')
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

    # Shine
    svg.append(f'<rect x="{tube_x + 15}" y="{tube_y + 18}" width="18" height="{tube_h - 36}" rx="9" fill="#ffffff" opacity="0.32" clip-path="url(#{clip_id})"/>')

    # Component boundary lines
    running = 0.0
    for comp in result.components[:-1]:
        running += comp.volume_ml
        y = y_for_volume(running)
        if tube_y < y < tube_y + tube_h:
            svg.append(f'<line x1="{tube_x+7}" x2="{tube_x+tube_w-7}" y1="{y}" y2="{y}" stroke="#ffffff" opacity="0.7" stroke-width="2" clip-path="url(#{clip_id})"/>')

    # Tube outline
    svg.append(f'<rect x="{tube_x}" y="{tube_y}" width="{tube_w}" height="{tube_h}" rx="{tube_w/2}" ry="{tube_w/2}" fill="none" stroke="#37474f" stroke-width="4"/>')

    # Volume marks
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

    axis_x = 54
    svg.append(f'<text x="{axis_x}" y="{tube_y + tube_h/2}" transform="rotate(-90 {axis_x} {tube_y + tube_h/2})" text-anchor="middle" font-size="16" font-family="Arial, sans-serif" fill="#37474f">volume [ml]</text>')
    svg.append(f'<text x="{tube_x + tube_w/2}" y="{tube_y + tube_h + 32}" text-anchor="middle" font-size="16" font-family="Arial, sans-serif" fill="#37474f">total = {total:.4g} ml</text>')

    # Legend with wrapping
    svg.append(f'<text x="{legend_x}" y="{tube_y + 8}" font-size="18" font-family="Arial, sans-serif" font-weight="bold" fill="#263238">Mixture</text>')
    legend_y = tube_y + 36
    for idx, lines in enumerate(legend_blocks):
        color = colors[idx % len(colors)]
        svg.append(f'<rect x="{legend_x}" y="{legend_y - 13}" width="16" height="16" fill="{color}" stroke="#455a64" stroke-width="1"/>')
        svg.append(_svg_tspans(lines, legend_text_x, legend_y, font_size=14, line_height=legend_line_height))
        legend_y += len(lines) * legend_line_height + legend_block_gap

    if result.apparent_density_g_ml is not None:
        density_txt = f"density = {result.apparent_density_g_ml:.6f} g/ml"
        svg.append(f'<text x="{legend_x}" y="{legend_y + 8}" font-size="14" font-family="Arial, sans-serif" fill="#263238">{escape(density_txt)}</text>')

    svg.append('</svg>')
    return ''.join(svg)

def render_tube_visualization(result: MixResult):
    st.subheader("Tube visualization")
    st.caption("Stacked by component volume. Volume marks use the final/target total volume.")
    svg = tube_svg(result)
    # Render SVG as HTML instead of st.image(). Some Streamlit/Pillow
    # versions try to decode SVG bytes as a raster image and fail with:
    # "cannot identify image file <_io.BytesIO ...>".
    components.html(f'<div style="width:100%; overflow-x:auto;">{svg}</div>', height=720, scrolling=True)
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

    viz_col, table_col = st.columns([1.35, 1.15])
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


def parse_fraction_list(raw: str) -> list[float]:
    """Parse comma/newline/semicolon-separated fractions.

    Values > 1 are interpreted as percentages, e.g. 50 -> 0.50.
    """
    tokens = [tok.strip() for tok in re.split(r"[,;\n\t ]+", raw) if tok.strip()]
    fractions: list[float] = []
    for tok in tokens:
        tok = tok.replace("%", "")
        value = float(tok)
        if value > 1:
            value = value / 100.0
        if not 0 <= value <= 1:
            raise ValueError(f"Fraction {tok!r} is outside 0–1 / 0–100%.")
        fractions.append(value)
    if not fractions:
        raise ValueError("Enter at least one IPS fraction.")
    return fractions


def make_paired_endpoint_protocol_text(recipe_df: pd.DataFrame, totals_df: pd.DataFrame, summary: dict, intervals: dict) -> str:
    lines = []
    lines.append("Paired low/high gradient endpoint batch")
    lines.append("=======================================")
    lines.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("Summary")
    lines.append("-------")
    for key, value in summary.items():
        if isinstance(value, float):
            lines.append(f"- {key}: {value:.10g}")
        else:
            lines.append(f"- {key}: {value}")
    lines.append("")
    lines.append("Feasible fixed-IPS intervals")
    lines.append("----------------------------")
    for endpoint, interval in intervals.items():
        lines.append(
            f"- {endpoint}: {100*interval['p_min']:.4f}% to {100*interval['p_max']:.4f}% IPS"
        )
    lines.append("")
    lines.append("Endpoint recipes")
    lines.append("----------------")
    for _, row in recipe_df.iterrows():
        lines.append(
            f"Condition {int(row['Condition'])}, {row['Endpoint']} endpoint, "
            f"target density {row['Target density [g/ml]']:.6f} g/ml, "
            f"IPS {row['IPS fraction [%]']:.3f}%: {row['Status']}"
        )
        if row["Status"] == "OK":
            lines.append(f"  - IPS:       {row['IPS volume [ml]']:.6f} ml = {row['IPS mass [g]']:.6f} g")
            lines.append(f"  - OptiPrep:  {row['OptiPrep volume [ml]']:.6f} ml = {row['OptiPrep mass [g]']:.6f} g")
            lines.append(f"  - PBS:       {row['PBS volume [ml]']:.6f} ml = {row['PBS mass [g]']:.6f} g")
        else:
            lines.append(f"  - Warning: {row['Warning']}")
    lines.append("")
    lines.append("Total original media required")
    lines.append("-----------------------------")
    for _, row in totals_df.iterrows():
        lines.append(
            f"- {row['Original medium']}: {row['Total volume needed [ml]']:.6f} ml "
            f"= {row['Total mass needed [g]']:.6f} g"
        )
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
    "Paired endpoints",
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
    st.header("Paired low/high endpoint batch")
    st.write(
        "Prepare paired low- and high-density endpoint media for several tubes. "
        "Each tube/condition keeps the same fixed IPS fraction in both endpoints; "
        "OptiPrep medium and PBS are solved separately for the low and high target densities."
    )

    c1, c2, c3 = st.columns([1.15, 1, 1])
    with c1:
        fractions_raw = st.text_area(
            "IPS fractions for conditions",
            "50, 60, 75, 87",
            help="Enter comma/newline-separated values. 50 means 50%; 0.5 also works.",
            key="paired_fractions",
        )
        endpoint_volume = st.number_input(
            "Nominal volume needed per endpoint [ml]",
            0.001, 10000.0, 6.0, 0.5, key="paired_endpoint_volume"
        )
        excess_percent = st.number_input(
            "Preparation excess [%]",
            0.0, 200.0, 10.0, 1.0,
            help="Scales each endpoint recipe. 10% means prepare 6.6 ml when 6 ml are needed.",
            key="paired_excess_percent",
        )
    with c2:
        rho_low = st.number_input(
            "Low endpoint density [g/ml]",
            0.5, 2.0, 1.0950, 0.0001, format="%.6f", key="paired_rho_low"
        )
        rho_high = st.number_input(
            "High endpoint density [g/ml]",
            0.5, 2.0, 1.1090, 0.0001, format="%.6f", key="paired_rho_high"
        )
        st.caption("Low and high density are applied to every condition.")
    with c3:
        rho_ips = st.number_input(
            "ρ IPS [g/ml]",
            0.5, 2.0, DEFAULTS["rho_ips74"], 0.0001, format="%.6f", key="paired_rho_ips"
        )
        rho_opt = st.number_input(
            "ρ OptiPrep medium [g/ml]",
            0.5, 2.0, DEFAULTS["rho_optiprep_40"], 0.0001, format="%.6f", key="paired_rho_opt"
        )
        rho_pbs = st.number_input(
            "ρ PBS [g/ml]",
            0.5, 2.0, DEFAULTS["rho_pbs"], 0.0001, format="%.6f", key="paired_rho_pbs"
        )

    st.subheader("Available stock volumes")
    a1, a2, a3 = st.columns(3)
    with a1:
        avail_ips = st.number_input("Available IPS [ml]", 0.0, 100000.0, 100.0, 1.0, key="paired_avail_ips")
    with a2:
        avail_opt = st.number_input("Available OptiPrep medium [ml]", 0.0, 100000.0, 100.0, 1.0, key="paired_avail_opt")
    with a3:
        avail_pbs = st.number_input("Available PBS [ml]", 0.0, 100000.0, 100.0, 1.0, key="paired_avail_pbs")

    try:
        ips_fractions = parse_fraction_list(fractions_raw)
        batch = paired_gradient_endpoint_batch(
            endpoint_volume_ml=endpoint_volume,
            ips_fractions=ips_fractions,
            rho_low_target=rho_low,
            rho_high_target=rho_high,
            rho_ips=rho_ips,
            rho_optiprep=rho_opt,
            rho_pbs=rho_pbs,
            excess_fraction=excess_percent / 100.0,
        )
        recipe_df = pd.DataFrame(batch["recipe_rows"])
        totals_df = pd.DataFrame(batch["total_rows"])
        summary = batch["summary"]
        intervals = batch["endpoint_feasible_intervals"]

        st.subheader("Feasibility")
        low_interval = intervals["Low"]
        high_interval = intervals["High"]
        combined_min = summary["combined_feasible_ips_min"]
        combined_max = summary["combined_feasible_ips_max"]
        f1, f2, f3 = st.columns(3)
        f1.metric("Low endpoint feasible IPS", f"{100*low_interval['p_min']:.2f}–{100*low_interval['p_max']:.2f}%")
        f2.metric("High endpoint feasible IPS", f"{100*high_interval['p_min']:.2f}–{100*high_interval['p_max']:.2f}%")
        f3.metric("Feasible for both endpoints", f"{100*combined_min:.2f}–{100*combined_max:.2f}%")

        warnings = recipe_df.loc[recipe_df["Status"] != "OK", ["Condition", "Endpoint", "IPS fraction [%]", "Warning"]]
        if not warnings.empty:
            st.error(
                "Some requested IPS fractions are not feasible with the chosen densities. "
                "The table gives the maximum/minimum feasible IPS fraction for each endpoint."
            )
            st.dataframe(warnings, use_container_width=True)
        else:
            st.success("All requested endpoint recipes are feasible with non-negative IPS, OptiPrep, and PBS volumes.")

        st.subheader("Endpoint recipe table")
        display_cols = [
            "Condition", "Endpoint", "Target density [g/ml]", "IPS fraction [%]",
            "Nominal endpoint volume [ml]", "Preparation volume incl. excess [ml]",
            "IPS volume [ml]", "OptiPrep volume [ml]", "PBS volume [ml]",
            "IPS mass [g]", "OptiPrep mass [g]", "PBS mass [g]",
            "Calculated density [g/ml]", "Status", "Warning",
        ]
        st.dataframe(
            recipe_df[display_cols].style.format({
                "Target density [g/ml]": "{:.6f}",
                "IPS fraction [%]": "{:.3f}",
                "Nominal endpoint volume [ml]": "{:.4f}",
                "Preparation volume incl. excess [ml]": "{:.4f}",
                "IPS volume [ml]": "{:.6f}",
                "OptiPrep volume [ml]": "{:.6f}",
                "PBS volume [ml]": "{:.6f}",
                "IPS mass [g]": "{:.6f}",
                "OptiPrep mass [g]": "{:.6f}",
                "PBS mass [g]": "{:.6f}",
                "Calculated density [g/ml]": "{:.6f}",
            }, na_rep="—"),
            use_container_width=True,
        )

        st.subheader("Total original media required")
        totals_df["Available [ml]"] = [
            avail_ips if name.startswith("IPS") else avail_opt if name.startswith("OptiPrep") else avail_pbs
            for name in totals_df["Original medium"]
        ]
        totals_df["Remaining [ml]"] = totals_df["Available [ml]"] - totals_df["Total volume needed [ml]"]
        totals_df["Enough stock?"] = totals_df["Remaining [ml]"] >= -1e-9

        st.dataframe(
            totals_df.style.format({
                "Density [g/ml]": "{:.6f}",
                "Total volume needed [ml]": "{:.6f}",
                "Total mass needed [g]": "{:.6f}",
                "Available [ml]": "{:.3f}",
                "Remaining [ml]": "{:.3f}",
            }),
            use_container_width=True,
        )

        if not totals_df["Enough stock?"].all():
            st.warning("At least one original medium is insufficient for the feasible recipes shown above.")
        else:
            st.info(
                f"Requested nominal endpoint volume: {summary['total_nominal_endpoint_volume_requested_ml']:.3f} ml. "
                f"Feasible prepared volume including excess: {summary['total_preparation_volume_feasible_ml']:.3f} ml "
                f"({int(summary['number_endpoints_feasible'])}/{int(summary['number_endpoints_requested'])} endpoints feasible)."
            )

        # Downloads
        protocol_text = make_paired_endpoint_protocol_text(recipe_df, totals_df, summary, intervals)
        csv_recipe = recipe_df.to_csv(index=False).encode("utf-8")
        csv_totals = totals_df.to_csv(index=False).encode("utf-8")
        combined_csv = (
            "# Endpoint recipes\n" + recipe_df.to_csv(index=False) +
            "\n# Total original media required\n" + totals_df.to_csv(index=False)
        ).encode("utf-8")

        d1, d2, d3 = st.columns(3)
        with d1:
            st.download_button(
                "Download full protocol .txt",
                protocol_text,
                file_name=f"paired_gradient_endpoint_protocol_{datetime.now().strftime('%Y%m%d_%H%M')}.txt",
                mime="text/plain",
            )
        with d2:
            st.download_button(
                "Download endpoint recipe CSV",
                csv_recipe,
                file_name=f"paired_gradient_endpoint_recipes_{datetime.now().strftime('%Y%m%d_%H%M')}.csv",
                mime="text/csv",
            )
        with d3:
            st.download_button(
                "Download totals CSV",
                csv_totals,
                file_name=f"paired_gradient_endpoint_totals_{datetime.now().strftime('%Y%m%d_%H%M')}.csv",
                mime="text/csv",
            )
        st.download_button(
            "Download combined CSV",
            combined_csv,
            file_name=f"paired_gradient_endpoint_full_{datetime.now().strftime('%Y%m%d_%H%M')}.csv",
            mime="text/csv",
        )

    except Exception as e:
        st.error(str(e))


with tabs[10]:
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
