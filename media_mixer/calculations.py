"""Density-media mixing calculations for Percoll, OptiPrep, PBS and blood additions.

The implemented equations assume additive volumes and linear density mixing. They are
intended as bench-calculation helpers; final experimental media should be checked with
an osmometer/pH meter/density meter when high precision is required.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

EPS = 1e-12


@dataclass(frozen=True)
class Component:
    name: str
    volume_ml: float
    density_g_ml: Optional[float] = None
    mass_g: Optional[float] = None
    fraction_vv: Optional[float] = None
    note: str = ""

    @staticmethod
    def from_volume(name: str, volume_ml: float, density_g_ml: Optional[float] = None,
                    total_volume_ml: Optional[float] = None, note: str = "") -> "Component":
        mass = None if density_g_ml is None else volume_ml * density_g_ml
        frac = None if total_volume_ml in (None, 0) else volume_ml / total_volume_ml
        return Component(name, volume_ml, density_g_ml, mass, frac, note)


@dataclass(frozen=True)
class MixResult:
    title: str
    target_volume_ml: float
    apparent_density_g_ml: Optional[float]
    components: List[Component]
    notes: List[str]
    checks: Dict[str, float]

    def as_rows(self) -> List[Dict[str, object]]:
        rows = []
        for c in self.components:
            rows.append({
                "Component": c.name,
                "Volume [ml]": c.volume_ml,
                "Density [g/ml]": c.density_g_ml,
                "Mass [g]": c.mass_g,
                "v/v fraction": c.fraction_vv,
                "Note": c.note,
            })
        return rows


def _validate_fraction(name: str, x: float, tol: float = 1e-9) -> None:
    if x < -tol or x > 1 + tol:
        raise ValueError(f"{name}={x:.6g} is outside [0, 1]. This target is not feasible with the chosen stocks.")


def _validate_nonnegative(name: str, x: float, tol: float = 1e-9) -> None:
    if x < -tol:
        raise ValueError(f"{name}={x:.6g} is negative. This target is not feasible with the chosen stocks.")


def linear_fraction_for_density(rho_target: float, rho_low: float, rho_high: float) -> float:
    """Fraction of high-density stock needed in a two-component mixture."""
    if abs(rho_high - rho_low) < EPS:
        raise ValueError("Stock densities are identical; cannot solve a density mixture.")
    f_high = (rho_target - rho_low) / (rho_high - rho_low)
    _validate_fraction("fraction_high", f_high)
    return f_high


def weighted_density(components: List[Component]) -> float:
    vol = sum(c.volume_ml for c in components)
    if vol <= 0:
        raise ValueError("Total volume must be positive.")
    mass = 0.0
    for c in components:
        if c.density_g_ml is None:
            raise ValueError(f"Missing density for {c.name}.")
        mass += c.volume_ml * c.density_g_ml
    return mass / vol


def fast_dirty_percoll(final_or_base_volume_ml: float, *, exact_final_volume: bool,
                       rho_percoll: float = 1.130, rho_pbs10: float = 1.040,
                       rho_pbs1: float = 1.0048) -> MixResult:
    """Fast approximate Percoll + 10xPBS + 1xPBS recipe.

    Legacy protocol: 0.9 V Percoll + 0.1 V 10xPBS + 0.11 V 1xPBS, so final is 1.11 V.
    If exact_final_volume is True, V is back-calculated as final/1.11.
    """
    if final_or_base_volume_ml <= 0:
        raise ValueError("Volume must be positive.")
    base = final_or_base_volume_ml / 1.11 if exact_final_volume else final_or_base_volume_ml
    final = 1.11 * base
    comps = [
        Component.from_volume("Percoll stock", 0.9 * base, rho_percoll, final),
        Component.from_volume("10x PBS", 0.1 * base, rho_pbs10, final),
        Component.from_volume("1x PBS", 0.11 * base, rho_pbs1, final),
    ]
    return MixResult(
        title="Fast approximate Percoll medium",
        target_volume_ml=final,
        apparent_density_g_ml=weighted_density(comps),
        components=comps,
        notes=[
            "Approximate protocol for quick biological testing and rough RBC-density matching.",
            "The legacy formula sums to 1.11 × V. Use exact-final mode when a precise final volume matters.",
        ],
        checks={"base_V_ml": base, "sum_volume_ml": final},
    )


def ips_basic(volume_ml: float, mode: str = "empirical", rho_percoll: float = 1.130,
              rho_pbs10: float = 1.0674) -> MixResult:
    if mode == "simple_90_10":
        p_pc, p_pbs10 = 0.9, 0.1
        title = "IPS basic 90:10"
    elif mode == "empirical":
        p_pc, p_pbs10 = 0.90493541, 0.09506459
        title = "IPS empirical osmolarity recipe"
    else:
        raise ValueError(f"Unknown IPS mode: {mode}")
    comps = [
        Component.from_volume("Percoll stock", p_pc * volume_ml, rho_percoll, volume_ml),
        Component.from_volume("10x PBS", p_pbs10 * volume_ml, rho_pbs10, volume_ml),
    ]
    return MixResult(
        title=title,
        target_volume_ml=volume_ml,
        apparent_density_g_ml=weighted_density(comps),
        components=comps,
        notes=["IPS density is calculated from supplied component densities; measure the final density for precision."],
        checks={"p_percoll": p_pc, "p_pbs10": p_pbs10},
    )


def ips74_osmotic(volume_ml: float, *, osm_target: float = 300.0, osm_pbs10: float = 3032.0,
                  osm_hcl: float = 355.0, osm_percoll: float = 13.0,
                  p_hcl: float = 0.032, rho_percoll: float = 1.130,
                  rho_pbs10: float = 1.0674, rho_hcl: float = 1.0,
                  empirical_fixed: bool = False) -> MixResult:
    if empirical_fixed:
        p_pc, p_pbs10, p_hcl = 0.87656045, 0.09143955, 0.032
        note = "Uses the empirical IPS74 recipe: 87.656045% Percoll, 9.143955% 10xPBS, 3.2% 0.178 M HCl."
    else:
        if abs(osm_percoll - osm_pbs10) < EPS:
            raise ValueError("osm_percoll and osm_pbs10 are identical; cannot solve osmotic balance.")
        p_pc = (osm_target - osm_pbs10 - p_hcl * (osm_hcl - osm_pbs10)) / (osm_percoll - osm_pbs10)
        p_pbs10 = 1.0 - p_pc - p_hcl
        _validate_fraction("p_percoll", p_pc)
        _validate_fraction("p_10xPBS", p_pbs10)
        _validate_fraction("p_HCl", p_hcl)
        note = "Calculated from target osmolarity and fixed HCl fraction. HCl fraction remains an empirical pH setting."
    comps = [
        Component.from_volume("Percoll stock", p_pc * volume_ml, rho_percoll, volume_ml),
        Component.from_volume("10x PBS", p_pbs10 * volume_ml, rho_pbs10, volume_ml),
        Component.from_volume("HCl solution", p_hcl * volume_ml, rho_hcl, volume_ml),
    ]
    return MixResult(
        title="IPS74 / pH-adjusted isotonic Percoll",
        target_volume_ml=volume_ml,
        apparent_density_g_ml=weighted_density(comps),
        components=comps,
        notes=[note, "Check pH and osmolarity experimentally for final use with RBCs."],
        checks={"p_percoll": p_pc, "p_pbs10": p_pbs10, "p_hcl": p_hcl},
    )


def sdm_ips_pbs(volume_ml: float, rho_target: float, rho_ips: float, rho_pbs: float) -> MixResult:
    p_ips = linear_fraction_for_density(rho_target, rho_pbs, rho_ips)
    p_pbs = 1 - p_ips
    comps = [
        Component.from_volume("IPS / dense Percoll stock", p_ips * volume_ml, rho_ips, volume_ml),
        Component.from_volume("1x PBS", p_pbs * volume_ml, rho_pbs, volume_ml),
    ]
    return MixResult(
        title="SDM from IPS + PBS",
        target_volume_ml=volume_ml,
        apparent_density_g_ml=weighted_density(comps),
        components=comps,
        notes=["Two-component volume-weighted density mixture."],
        checks={"p_ips": p_ips, "p_pbs": p_pbs},
    )


def sdm_with_whole_blood(final_volume_ml: float, rho_target_liquid: float, rho_ips: float,
                         rho_pbs: float, target_rbc_fraction: float, sample_hematocrit: float,
                         rho_sample_liquid: float = 1.0255, rho_rbc_for_mass: Optional[float] = None,
                         sample_liquid_name: str = "plasma / blood supernatant") -> MixResult:
    """Prepare IPS/PBS liquid so final liquid phase reaches rho_target after whole blood addition.

    The target density is the continuous/suspending liquid phase after plasma/buffer from the blood
    sample is included; the RBC volume is excluded from the liquid density balance.
    """
    _validate_fraction("target_rbc_fraction", target_rbc_fraction)
    _validate_fraction("sample_hematocrit", sample_hematocrit)
    if sample_hematocrit <= 0:
        raise ValueError("sample_hematocrit must be positive.")
    if target_rbc_fraction >= sample_hematocrit:
        raise ValueError("target RBC fraction must be lower than the sample hematocrit for dilution with medium.")

    v_rbc = final_volume_ml * target_rbc_fraction
    v_blood = v_rbc / sample_hematocrit
    v_sample_liquid = v_blood - v_rbc
    v_liquid = final_volume_ml - v_rbc
    if v_liquid <= 0:
        raise ValueError("No liquid volume remains after RBC volume; choose lower target RBC fraction.")
    f_sample_liquid = v_sample_liquid / v_liquid
    p_ips = (rho_target_liquid - rho_pbs - f_sample_liquid * (rho_sample_liquid - rho_pbs)) / (rho_ips - rho_pbs)
    p_pbs = 1 - p_ips - f_sample_liquid
    _validate_fraction("f_sample_liquid", f_sample_liquid)
    _validate_fraction("p_ips", p_ips)
    _validate_fraction("p_pbs", p_pbs)

    v_ips = p_ips * v_liquid
    v_pbs = p_pbs * v_liquid
    rho_rbc_for_mass = rho_rbc_for_mass if rho_rbc_for_mass is not None else rho_target_liquid
    est_blood_mass = v_sample_liquid * rho_sample_liquid + v_rbc * rho_rbc_for_mass
    comps = [
        Component.from_volume("IPS / dense Percoll stock", v_ips, rho_ips, final_volume_ml),
        Component.from_volume("1x PBS", v_pbs, rho_pbs, final_volume_ml),
        Component.from_volume("Whole blood sample", v_blood, None, final_volume_ml,
                              note=f"Contains {v_rbc:.6g} ml RBC + {v_sample_liquid:.6g} ml {sample_liquid_name}. Estimated mass: {est_blood_mass:.6g} g."),
    ]
    return MixResult(
        title="SDM + whole blood addition",
        target_volume_ml=final_volume_ml,
        apparent_density_g_ml=rho_target_liquid,
        components=comps,
        notes=[
            "Target density refers to the continuous liquid phase after the blood-sample liquid is included; RBC volume is not part of this density balance.",
            "Whole-blood mass is only an estimate unless RBC density is measured.",
        ],
        checks={
            "v_rbc_ml": v_rbc,
            "v_blood_ml": v_blood,
            "v_sample_liquid_ml": v_sample_liquid,
            "v_liquid_phase_ml": v_liquid,
            "f_sample_liquid_in_liquid_phase": f_sample_liquid,
            "p_ips_in_liquid_phase": p_ips,
            "p_pbs_in_liquid_phase": p_pbs,
            "estimated_blood_mass_g": est_blood_mass,
        },
    )


def density_correction_add_stock(existing_volume_ml: float, rho_existing: float,
                                 rho_target: float, rho_stock_to_add: float,
                                 stock_name: str = "correction stock") -> MixResult:
    """Add one stock to an existing medium to correct density."""
    denom = rho_stock_to_add - rho_target
    if abs(denom) < EPS:
        raise ValueError("Stock density equals target density; cannot solve added volume.")
    v_add = existing_volume_ml * (rho_target - rho_existing) / denom
    _validate_nonnegative("added volume", v_add)
    final_vol = existing_volume_ml + v_add
    comps = [
        Component.from_volume("Existing medium", existing_volume_ml, rho_existing, final_vol),
        Component.from_volume(stock_name, v_add, rho_stock_to_add, final_vol),
    ]
    return MixResult(
        title="Density correction by adding one stock",
        target_volume_ml=final_vol,
        apparent_density_g_ml=weighted_density(comps),
        components=comps,
        notes=["Target density must lie between existing medium and correction stock."],
        checks={"added_volume_ml": v_add, "final_volume_ml": final_vol},
    )


def optiprep_working_solution(final_volume_ml: float, target_iodixanol_percent_wv: float = 40.0,
                              stock_iodixanol_percent_wv: float = 60.0,
                              rho_optiprep_stock: float = 1.320,
                              rho_diluent: float = 1.006,
                              diluent_name: str = "concentrated diluent") -> MixResult:
    """Make an iodixanol working solution from raw OptiPrep stock and diluent by concentration balance."""
    if target_iodixanol_percent_wv <= 0 or stock_iodixanol_percent_wv <= 0:
        raise ValueError("Iodixanol percentages must be positive.")
    if target_iodixanol_percent_wv > stock_iodixanol_percent_wv:
        raise ValueError("Target iodixanol concentration cannot exceed raw stock concentration.")
    p_stock = target_iodixanol_percent_wv / stock_iodixanol_percent_wv
    p_diluent = 1 - p_stock
    comps = [
        Component.from_volume("Raw OptiPrep stock", p_stock * final_volume_ml, rho_optiprep_stock, final_volume_ml,
                              note=f"{stock_iodixanol_percent_wv:g}% w/v iodixanol"),
        Component.from_volume(diluent_name, p_diluent * final_volume_ml, rho_diluent, final_volume_ml),
    ]
    notes = [
        "For mammalian-cell 40% WS, the data sheet uses 2 vol OptiPrep + 1 vol diluent.",
        "Density from linear mixing is approximate; use data-sheet values or measured density when available.",
    ]
    return MixResult(
        title=f"OptiPrep {target_iodixanol_percent_wv:g}% w/v working solution",
        target_volume_ml=final_volume_ml,
        apparent_density_g_ml=weighted_density(comps),
        components=comps,
        notes=notes,
        checks={"p_raw_optiprep": p_stock, "p_diluent": p_diluent},
    )


def optiprep_table_density(target_iodixanol_percent_wv: float, medium: str = "NaCl/Tricine") -> Tuple[float, str]:
    """Piecewise-linear interpolation of data-sheet density tables for 40% WS diluted with CSM/DMEM."""
    tables = {
        "NaCl/Tricine": {
            6: 1.037, 8: 1.048, 10: 1.058, 12: 1.069, 14: 1.079, 16: 1.090,
            18: 1.100, 20: 1.111, 22: 1.121, 24: 1.132, 26: 1.142, 28: 1.153,
            30: 1.163, 32: 1.174, 34: 1.184, 36: 1.195, 38: 1.205, 40: 1.215,
        },
        "DMEM": {
            6: 1.038, 8: 1.049, 10: 1.059, 12: 1.070, 14: 1.080, 16: 1.090,
            18: 1.101, 20: 1.111, 22: 1.122, 24: 1.132, 26: 1.143, 28: 1.153,
            30: 1.163, 32: 1.174, 34: 1.184, 36: 1.195, 38: 1.205, 40: 1.216,
        },
    }
    if medium not in tables:
        raise ValueError(f"Unknown table medium: {medium}")
    table = tables[medium]
    xs = sorted(table.keys())
    x = target_iodixanol_percent_wv
    if x < xs[0] or x > xs[-1]:
        raise ValueError(f"Iodixanol percentage must be between {xs[0]} and {xs[-1]} for table interpolation.")
    if x in table:
        return table[x], "exact table value"
    for a, b in zip(xs[:-1], xs[1:]):
        if a <= x <= b:
            y = table[a] + (table[b] - table[a]) * (x - a) / (b - a)
            return y, f"linear interpolation between {a:g}% and {b:g}%"
    raise AssertionError("Interpolation failed.")


def optiprep_endpoint_fibrinogen(endpoint_volume_ml: float, target_fib_mg_ml: float = 5.0,
                                 fib_stock_mg_ml: float = 50.0,
                                 base_name: str = "OptiPrep 40% working solution",
                                 rho_base: float = 1.215,
                                 rho_fib_stock: float = 1.02) -> MixResult:
    """Create an endpoint medium with a desired fibrinogen concentration from fib stock + base."""
    if fib_stock_mg_ml <= 0 or target_fib_mg_ml < 0:
        raise ValueError("Fibrinogen concentrations must be valid.")
    if target_fib_mg_ml > fib_stock_mg_ml:
        raise ValueError("Target fibrinogen concentration cannot exceed fibrinogen stock concentration.")
    p_fib_stock = target_fib_mg_ml / fib_stock_mg_ml
    p_base = 1 - p_fib_stock
    comps = [
        Component.from_volume(base_name, p_base * endpoint_volume_ml, rho_base, endpoint_volume_ml),
        Component.from_volume(f"Fibrinogen stock ({fib_stock_mg_ml:g} mg/ml)", p_fib_stock * endpoint_volume_ml,
                              rho_fib_stock, endpoint_volume_ml),
    ]
    return MixResult(
        title=f"{base_name} + fibrinogen endpoint",
        target_volume_ml=endpoint_volume_ml,
        apparent_density_g_ml=weighted_density(comps),
        components=comps,
        notes=["Endpoint density should be measured after preparation for precise downstream mixing."],
        checks={"p_fib_stock": p_fib_stock, "p_base": p_base, "target_fib_mg_ml": target_fib_mg_ml},
    )


def two_endpoint_density_mix(final_volume_ml: float, rho_target: float, rho_low: float, rho_high: float,
                             low_name: str = "low-density endpoint",
                             high_name: str = "high-density endpoint") -> MixResult:
    f_high = linear_fraction_for_density(rho_target, rho_low, rho_high)
    f_low = 1 - f_high
    comps = [
        Component.from_volume(high_name, f_high * final_volume_ml, rho_high, final_volume_ml),
        Component.from_volume(low_name, f_low * final_volume_ml, rho_low, final_volume_ml),
    ]
    return MixResult(
        title="Two-endpoint density mixture",
        target_volume_ml=final_volume_ml,
        apparent_density_g_ml=weighted_density(comps),
        components=comps,
        notes=["Use when both endpoints contain the same additive concentration, e.g. same fibrinogen concentration."],
        checks={"p_high": f_high, "p_low": f_low},
    )


def ternary_ips_optiprep_pbs(final_volume_ml: float, rho_target: float, fixed_fraction: float,
                             rho_ips: float, rho_optiprep: float, rho_pbs: float,
                             fixed_fraction_kind: str = "IPS fraction",
                             percoll_fraction_in_ips: float = 0.9) -> MixResult:
    """Ternary mix with fixed IPS/Percoll-containing fraction and solve OptiPrep/PBS fractions.

    fixed_fraction_kind:
      - "IPS fraction": fixed_fraction is final v/v of IPS stock.
      - "Percoll-equivalent fraction": fixed_fraction is desired stock-Percoll-equivalent fraction;
        it is converted to IPS stock fraction by dividing by percoll_fraction_in_ips.
    """
    if fixed_fraction_kind == "IPS fraction":
        p_ips = fixed_fraction
    elif fixed_fraction_kind == "Percoll-equivalent fraction":
        if percoll_fraction_in_ips <= 0:
            raise ValueError("percoll_fraction_in_ips must be positive.")
        p_ips = fixed_fraction / percoll_fraction_in_ips
    else:
        raise ValueError(f"Unknown fixed fraction kind: {fixed_fraction_kind}")
    _validate_fraction("p_ips", p_ips)
    if abs(rho_optiprep - rho_pbs) < EPS:
        raise ValueError("OptiPrep and PBS densities are identical; cannot solve ternary mixture.")
    p_opt = (rho_target - p_ips * rho_ips - (1 - p_ips) * rho_pbs) / (rho_optiprep - rho_pbs)
    p_pbs = 1 - p_ips - p_opt
    _validate_fraction("p_optiprep", p_opt)
    _validate_fraction("p_pbs", p_pbs)
    comps = [
        Component.from_volume("IPS / Percoll-containing stock", p_ips * final_volume_ml, rho_ips, final_volume_ml),
        Component.from_volume("OptiPrep working solution", p_opt * final_volume_ml, rho_optiprep, final_volume_ml),
        Component.from_volume("1x PBS", p_pbs * final_volume_ml, rho_pbs, final_volume_ml),
    ]
    notes = ["Fixed IPS/Percoll fraction controls Percoll content; OptiPrep/PBS complete the density balance."]
    if fixed_fraction_kind == "Percoll-equivalent fraction":
        notes.append("Percoll-equivalent mode assumes the entered Percoll fraction is contained inside IPS, not pure Percoll added directly.")
    return MixResult(
        title="Ternary IPS/Percoll + OptiPrep + PBS mixture",
        target_volume_ml=final_volume_ml,
        apparent_density_g_ml=weighted_density(comps),
        components=comps,
        notes=notes,
        checks={"p_ips": p_ips, "p_optiprep": p_opt, "p_pbs": p_pbs},
    )


def rbc_density_gradient_targets(rho_rbc_mean: float, half_width: float = 0.005) -> Dict[str, float]:
    return {"rho_low": rho_rbc_mean - half_width, "rho_high": rho_rbc_mean + half_width}
