"""Density-media mixing calculations for Percoll, OptiPrep, PBS and blood additions.

The implemented equations assume additive volumes and linear density mixing. They are
intended as bench-calculation helpers; final experimental media should be checked with
an osmometer/pH meter/density meter when high precision is required.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any

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


def _linear_ge_interval(c: float, d: float) -> Tuple[Optional[float], Optional[float], bool]:
    """Return interval contribution for c + d*p >= 0.

    Returns (lower, upper, feasible). None means unbounded on that side.
    """
    if abs(d) < EPS:
        return (None, None, c >= -1e-12)
    root = -c / d
    if d > 0:
        return (root, None, True)
    return (None, root, True)


def feasible_fixed_ips_interval_for_density(rho_target: float, rho_ips: float,
                                            rho_optiprep: float, rho_pbs: float) -> Dict[str, float]:
    """Feasible fixed IPS-fraction interval for a ternary IPS/OptiPrep/PBS target.

    For a fixed IPS fraction p, the remaining 1-p is filled by OptiPrep and PBS.
    This function returns the p interval for which both required OptiPrep and PBS
    fractions are non-negative.
    """
    if abs(rho_optiprep - rho_pbs) < EPS:
        raise ValueError("OptiPrep and PBS densities are identical; cannot solve fixed-IPS feasibility.")

    # O(p) = a + b*p
    a = (rho_target - rho_pbs) / (rho_optiprep - rho_pbs)
    b = -(rho_ips - rho_pbs) / (rho_optiprep - rho_pbs)

    # PBS fraction B(p) = 1 - p - O(p) = c + d*p
    c = 1.0 - a
    d = -1.0 - b

    lower, upper = 0.0, 1.0
    for cc, dd in [(a, b), (c, d)]:
        lo, hi, ok = _linear_ge_interval(cc, dd)
        if not ok:
            return {
                "p_min": float("nan"),
                "p_max": float("nan"),
                "feasible": 0.0,
                "a_opt": a,
                "b_opt": b,
            }
        if lo is not None:
            lower = max(lower, lo)
        if hi is not None:
            upper = min(upper, hi)

    feasible = lower <= upper + 1e-9
    return {
        "p_min": max(0.0, lower),
        "p_max": min(1.0, upper),
        "feasible": 1.0 if feasible else 0.0,
        "a_opt": a,
        "b_opt": b,
    }


def solve_fixed_ips_ternary_fractions(rho_target: float, p_ips: float, rho_ips: float,
                                      rho_optiprep: float, rho_pbs: float) -> Dict[str, float]:
    """Solve OptiPrep/PBS fractions for a fixed IPS fraction and target density."""
    _validate_fraction("p_ips", p_ips)
    if abs(rho_optiprep - rho_pbs) < EPS:
        raise ValueError("OptiPrep and PBS densities are identical; cannot solve ternary mixture.")
    p_opt = (rho_target - p_ips * rho_ips - (1 - p_ips) * rho_pbs) / (rho_optiprep - rho_pbs)
    p_pbs = 1.0 - p_ips - p_opt
    rho_check = p_ips * rho_ips + p_opt * rho_optiprep + p_pbs * rho_pbs
    interval = feasible_fixed_ips_interval_for_density(rho_target, rho_ips, rho_optiprep, rho_pbs)
    feasible = (
        interval["feasible"] == 1.0
        and p_ips >= interval["p_min"] - 1e-9
        and p_ips <= interval["p_max"] + 1e-9
        and p_opt >= -1e-9
        and p_pbs >= -1e-9
    )
    return {
        "p_ips": p_ips,
        "p_optiprep": p_opt,
        "p_pbs": p_pbs,
        "rho_check": rho_check,
        "feasible": 1.0 if feasible else 0.0,
        "p_ips_min_feasible": interval["p_min"],
        "p_ips_max_feasible": interval["p_max"],
    }



def _blood_sample_fractions(target_rbc_fraction: float, sample_hematocrit: float) -> Dict[str, float]:
    """Return final-volume fractions for an RBC suspension addition.

    target_rbc_fraction is the desired final RBC volume fraction in the tube.
    sample_hematocrit is the RBC volume fraction inside the stock blood/RBC suspension.
    """
    _validate_fraction("target_rbc_fraction", target_rbc_fraction)
    _validate_fraction("sample_hematocrit", sample_hematocrit)
    if sample_hematocrit <= 0:
        raise ValueError("sample_hematocrit must be positive.")
    if target_rbc_fraction >= sample_hematocrit:
        raise ValueError(
            "target_rbc_fraction must be lower than sample_hematocrit; otherwise no carrier medium can be added."
        )
    p_rbc = target_rbc_fraction
    p_sample_total = p_rbc / sample_hematocrit if p_rbc > 0 else 0.0
    p_sample_liquid = p_sample_total - p_rbc
    p_carrier = 1.0 - p_sample_total
    p_liquid = 1.0 - p_rbc
    return {
        "p_rbc": p_rbc,
        "p_sample_total": p_sample_total,
        "p_sample_liquid": p_sample_liquid,
        "p_carrier": p_carrier,
        "p_liquid": p_liquid,
    }


def _normalise_ips_fraction_basis(ips_fraction_basis: str) -> str:
    key = (ips_fraction_basis or "").strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "carrier": "carrier",
        "carrier_excluding_sample": "carrier",
        "carrier_excluding_rbc_suspension": "carrier",
        "carrier_before_rbc_suspension": "carrier",
        "formulated_carrier": "carrier",
        "final": "final_total",
        "final_total": "final_total",
        "final_total_mixture": "final_total",
        "final_mixture": "final_total",
    }
    if key not in aliases:
        raise ValueError(
            "ips_fraction_basis must be 'carrier_excluding_sample' or 'final_total_mixture'."
        )
    return aliases[key]


def feasible_fixed_ips_interval_with_sample(
    rho_target: float,
    rho_ips: float,
    rho_optiprep: float,
    rho_pbs: float,
    *,
    target_rbc_fraction: float = 0.0,
    sample_hematocrit: float = 0.45,
    rho_sample_liquid: float = 1.0255,
    ips_fraction_basis: str = "carrier_excluding_sample",
) -> Dict[str, float]:
    """Feasible requested-IPS interval for a density endpoint with optional RBC suspension.

    If ips_fraction_basis is "carrier", the requested IPS fraction is the IPS fraction
    within the formulated carrier only, excluding the RBC suspension volume.

    If ips_fraction_basis is "final_total", the requested IPS fraction is the final
    tube-volume fraction including the RBC suspension volume.
    """
    if abs(rho_optiprep - rho_pbs) < EPS:
        raise ValueError("OptiPrep and PBS densities are identical; cannot solve fixed-IPS feasibility.")

    basis = _normalise_ips_fraction_basis(ips_fraction_basis)
    blood = _blood_sample_fractions(target_rbc_fraction, sample_hematocrit)
    p_rbc = blood["p_rbc"]
    p_sample_liquid = blood["p_sample_liquid"]
    p_carrier = blood["p_carrier"]
    p_liquid = blood["p_liquid"]

    if p_carrier <= EPS:
        return {
            "p_min": float("nan"),
            "p_max": float("nan"),
            "feasible": 0.0,
            "basis": 0.0 if basis == "carrier" else 1.0,
        }

    denom = rho_optiprep - rho_pbs
    adjusted_mass_density_term = rho_target * p_liquid - p_sample_liquid * rho_sample_liquid

    if basis == "carrier":
        # q_opt(x) = a + b*x, q_pbs(x) = 1 - x - q_opt(x)
        # where x is IPS fraction inside the formulated carrier.
        lower, upper = 0.0, 1.0
        a = (adjusted_mass_density_term - p_carrier * rho_pbs) / (p_carrier * denom)
        b = -(rho_ips - rho_pbs) / denom
    else:
        # p_opt(x) = a + b*x, p_pbs(x) = p_carrier - x - p_opt(x)
        # where x is final total-volume fraction of IPS.
        lower, upper = 0.0, p_carrier
        a = (adjusted_mass_density_term - p_carrier * rho_pbs) / denom
        b = -(rho_ips - rho_pbs) / denom

    # OptiPrep >= 0
    lo, hi, ok = _linear_ge_interval(a, b)
    if not ok:
        feasible = False
    else:
        feasible = True
        if lo is not None:
            lower = max(lower, lo)
        if hi is not None:
            upper = min(upper, hi)

    # PBS >= 0
    # carrier basis: q_pbs = 1 - x - q_opt = (1-a) + (-1-b)x
    # final basis:   p_pbs = p_carrier - x - p_opt = (p_carrier-a) + (-1-b)x
    c = (1.0 - a) if basis == "carrier" else (p_carrier - a)
    d = -1.0 - b
    lo, hi, ok = _linear_ge_interval(c, d)
    if not ok:
        feasible = False
    else:
        if lo is not None:
            lower = max(lower, lo)
        if hi is not None:
            upper = min(upper, hi)

    feasible = feasible and lower <= upper + 1e-9
    return {
        "p_min": max(0.0, lower) if feasible else float("nan"),
        "p_max": min(1.0 if basis == "carrier" else p_carrier, upper) if feasible else float("nan"),
        "feasible": 1.0 if feasible else 0.0,
        "basis": 0.0 if basis == "carrier" else 1.0,
        "p_carrier_fraction_final_volume": p_carrier,
        "p_sample_total": blood["p_sample_total"],
        "p_sample_liquid": p_sample_liquid,
        "p_rbc": p_rbc,
        "a_opt": a,
        "b_opt": b,
    }


def solve_fixed_ips_ternary_fractions_with_sample(
    rho_target: float,
    ips_fraction: float,
    rho_ips: float,
    rho_optiprep: float,
    rho_pbs: float,
    *,
    target_rbc_fraction: float = 0.0,
    sample_hematocrit: float = 0.45,
    rho_sample_liquid: float = 1.0255,
    ips_fraction_basis: str = "carrier_excluding_sample",
) -> Dict[str, float]:
    """Solve IPS/OptiPrep/PBS fractions with optional RBC suspension addition."""
    basis = _normalise_ips_fraction_basis(ips_fraction_basis)
    blood = _blood_sample_fractions(target_rbc_fraction, sample_hematocrit)
    p_rbc = blood["p_rbc"]
    p_sample_total = blood["p_sample_total"]
    p_sample_liquid = blood["p_sample_liquid"]
    p_carrier = blood["p_carrier"]
    p_liquid = blood["p_liquid"]

    interval = feasible_fixed_ips_interval_with_sample(
        rho_target,
        rho_ips,
        rho_optiprep,
        rho_pbs,
        target_rbc_fraction=target_rbc_fraction,
        sample_hematocrit=sample_hematocrit,
        rho_sample_liquid=rho_sample_liquid,
        ips_fraction_basis=ips_fraction_basis,
    )

    if basis == "carrier":
        _validate_fraction("ips_fraction_in_carrier", ips_fraction)
        q_ips = ips_fraction
        if p_carrier <= EPS:
            raise ValueError("No carrier volume remains after RBC suspension addition.")
        q_opt = (
            rho_target * p_liquid
            - p_sample_liquid * rho_sample_liquid
            - p_carrier * (q_ips * rho_ips + (1.0 - q_ips) * rho_pbs)
        ) / (p_carrier * (rho_optiprep - rho_pbs))
        q_pbs = 1.0 - q_ips - q_opt
        p_ips_final = p_carrier * q_ips
        p_opt_final = p_carrier * q_opt
        p_pbs_final = p_carrier * q_pbs
        requested_basis_fraction = q_ips
        carrier_ips_fraction = q_ips
    else:
        _validate_fraction("ips_fraction_final_total", ips_fraction)
        p_ips_final = ips_fraction
        p_opt_final = (
            rho_target * p_liquid
            - p_sample_liquid * rho_sample_liquid
            - p_ips_final * rho_ips
            - (p_carrier - p_ips_final) * rho_pbs
        ) / (rho_optiprep - rho_pbs)
        p_pbs_final = p_carrier - p_ips_final - p_opt_final
        carrier_ips_fraction = p_ips_final / p_carrier if p_carrier > EPS else float("nan")
        q_opt = p_opt_final / p_carrier if p_carrier > EPS else float("nan")
        q_pbs = p_pbs_final / p_carrier if p_carrier > EPS else float("nan")
        requested_basis_fraction = p_ips_final

    rho_check = (
        p_ips_final * rho_ips
        + p_opt_final * rho_optiprep
        + p_pbs_final * rho_pbs
        + p_sample_liquid * rho_sample_liquid
    ) / p_liquid

    feasible = (
        interval["feasible"] == 1.0
        and requested_basis_fraction >= interval["p_min"] - 1e-9
        and requested_basis_fraction <= interval["p_max"] + 1e-9
        and p_opt_final >= -1e-9
        and p_pbs_final >= -1e-9
        and p_ips_final >= -1e-9
        and p_carrier >= -1e-9
    )

    return {
        "requested_ips_fraction": requested_basis_fraction,
        "carrier_ips_fraction": carrier_ips_fraction,
        "carrier_optiprep_fraction": q_opt,
        "carrier_pbs_fraction": q_pbs,
        "p_ips": p_ips_final,
        "p_optiprep": p_opt_final,
        "p_pbs": p_pbs_final,
        "p_rbc": p_rbc,
        "p_sample_total": p_sample_total,
        "p_sample_liquid": p_sample_liquid,
        "p_carrier": p_carrier,
        "p_liquid": p_liquid,
        "rho_check": rho_check,
        "feasible": 1.0 if feasible else 0.0,
        "p_ips_min_feasible": interval["p_min"],
        "p_ips_max_feasible": interval["p_max"],
    }


def _expand_condition_values(values: Optional[List[float]], n: int, default: float, name: str) -> List[float]:
    if values is None or len(values) == 0:
        return [default] * n
    if len(values) == 1:
        return values * n
    if len(values) != n:
        raise ValueError(f"{name} must contain either one value or exactly {n} values.")
    return values


def paired_gradient_endpoint_batch(endpoint_volume_ml: float, ips_fractions: List[float],
                                   rho_low_target: float, rho_high_target: float,
                                   rho_ips: float, rho_optiprep: float, rho_pbs: float,
                                   excess_fraction: float = 0.10,
                                   *,
                                   include_rbc_suspension: bool = False,
                                   target_rbc_fractions: Optional[List[float]] = None,
                                   sample_hematocrit: float = 0.45,
                                   rho_sample_liquid: float = 1.0255,
                                   sample_liquid_name: str = "sample liquid",
                                   rho_rbc_for_mass: float = 1.100,
                                   ips_fraction_basis: str = "carrier_excluding_sample") -> Dict[str, Any]:
    """Batch recipe for paired low/high endpoints with fixed IPS fractions.

    Each condition receives one low-density endpoint and one high-density endpoint.
    The low and high endpoint of a condition share the same requested IPS fraction.
    OptiPrep and PBS are calculated independently for the low/high target densities.

    endpoint_volume_ml is the nominal final volume needed per endpoint. The recipe is
    scaled by (1 + excess_fraction) so users can prepare extra medium for pipette
    dead volume, density checks, or losses.

    If include_rbc_suspension is True, endpoint_volume_ml is interpreted as the final
    tube mixture volume including RBCs and the stock RBC suspension. The density target
    is applied to the continuous liquid phase after the sample liquid has been included.
    RBC volume is excluded from the density denominator, following the blood-aware SDM
    convention used elsewhere in this app.
    """
    if endpoint_volume_ml <= 0:
        raise ValueError("endpoint_volume_ml must be positive.")
    if excess_fraction < 0:
        raise ValueError("excess_fraction must be non-negative.")
    if not ips_fractions:
        raise ValueError("At least one IPS fraction is required.")

    basis = _normalise_ips_fraction_basis(ips_fraction_basis)
    rbc_values = _expand_condition_values(
        target_rbc_fractions if include_rbc_suspension else [0.0],
        len(ips_fractions),
        0.0,
        "target_rbc_fractions",
    )
    if not include_rbc_suspension:
        rbc_values = [0.0] * len(ips_fractions)

    prep_volume_ml = endpoint_volume_ml * (1.0 + excess_fraction)
    targets = [("Low", rho_low_target), ("High", rho_high_target)]
    rows: List[Dict[str, Any]] = []

    totals = {
        "IPS / Percoll-containing stock": {"volume_ml": 0.0, "density_g_ml": rho_ips, "mass_g": 0.0},
        "OptiPrep medium": {"volume_ml": 0.0, "density_g_ml": rho_optiprep, "mass_g": 0.0},
        "1x PBS": {"volume_ml": 0.0, "density_g_ml": rho_pbs, "mass_g": 0.0},
    }
    if include_rbc_suspension:
        totals["RBC suspension sample"] = {"volume_ml": 0.0, "density_g_ml": float("nan"), "mass_g": 0.0}

    interval_groups: Dict[str, List[Dict[str, float]]] = {"Low": [], "High": []}

    for endpoint_name, rho_target in targets:
        for rbc_fraction in rbc_values:
            interval_groups[endpoint_name].append(
                feasible_fixed_ips_interval_with_sample(
                    rho_target,
                    rho_ips,
                    rho_optiprep,
                    rho_pbs,
                    target_rbc_fraction=rbc_fraction,
                    sample_hematocrit=sample_hematocrit,
                    rho_sample_liquid=rho_sample_liquid,
                    ips_fraction_basis=basis,
                )
            )

    endpoint_maxima: Dict[str, Dict[str, float]] = {}
    for endpoint_name, intervals in interval_groups.items():
        feasible_intervals = [itv for itv in intervals if itv["feasible"] == 1.0]
        if feasible_intervals:
            endpoint_maxima[endpoint_name] = {
                "p_min": max(itv["p_min"] for itv in feasible_intervals),
                "p_max": min(itv["p_max"] for itv in feasible_intervals),
                "feasible": 1.0,
                "basis": 0.0 if basis == "carrier" else 1.0,
            }
            if endpoint_maxima[endpoint_name]["p_min"] > endpoint_maxima[endpoint_name]["p_max"] + 1e-9:
                endpoint_maxima[endpoint_name]["feasible"] = 0.0
        else:
            endpoint_maxima[endpoint_name] = {
                "p_min": float("nan"),
                "p_max": float("nan"),
                "feasible": 0.0,
                "basis": 0.0 if basis == "carrier" else 1.0,
            }

    combined_p_max = min(endpoint_maxima["Low"]["p_max"], endpoint_maxima["High"]["p_max"])
    combined_p_min = max(endpoint_maxima["Low"]["p_min"], endpoint_maxima["High"]["p_min"])

    for condition_index, (p_ips_request, rbc_fraction) in enumerate(zip(ips_fractions, rbc_values), start=1):
        _validate_fraction(f"IPS fraction #{condition_index}", p_ips_request)
        _validate_fraction(f"RBC target fraction #{condition_index}", rbc_fraction)

        for endpoint_name, rho_target in targets:
            solved = solve_fixed_ips_ternary_fractions_with_sample(
                rho_target,
                p_ips_request,
                rho_ips,
                rho_optiprep,
                rho_pbs,
                target_rbc_fraction=rbc_fraction,
                sample_hematocrit=sample_hematocrit,
                rho_sample_liquid=rho_sample_liquid,
                ips_fraction_basis=basis,
            )
            feasible = solved["feasible"] == 1.0
            p_opt = solved["p_optiprep"]
            p_pbs = solved["p_pbs"]

            if feasible:
                v_ips = solved["p_ips"] * prep_volume_ml
                v_opt = p_opt * prep_volume_ml
                v_pbs = p_pbs * prep_volume_ml
                v_rbc = solved["p_rbc"] * prep_volume_ml
                v_sample = solved["p_sample_total"] * prep_volume_ml
                v_sample_liquid = solved["p_sample_liquid"] * prep_volume_ml
                v_carrier = solved["p_carrier"] * prep_volume_ml

                m_ips = v_ips * rho_ips
                m_opt = v_opt * rho_optiprep
                m_pbs = v_pbs * rho_pbs
                m_sample_liquid = v_sample_liquid * rho_sample_liquid
                m_rbc = v_rbc * rho_rbc_for_mass
                m_sample = m_sample_liquid + m_rbc if include_rbc_suspension else 0.0

                totals["IPS / Percoll-containing stock"]["volume_ml"] += v_ips
                totals["IPS / Percoll-containing stock"]["mass_g"] += m_ips
                totals["OptiPrep medium"]["volume_ml"] += v_opt
                totals["OptiPrep medium"]["mass_g"] += m_opt
                totals["1x PBS"]["volume_ml"] += v_pbs
                totals["1x PBS"]["mass_g"] += m_pbs
                if include_rbc_suspension:
                    totals["RBC suspension sample"]["volume_ml"] += v_sample
                    totals["RBC suspension sample"]["mass_g"] += m_sample
            else:
                v_ips = v_opt = v_pbs = v_rbc = v_sample = v_sample_liquid = v_carrier = float("nan")
                m_ips = m_opt = m_pbs = m_sample_liquid = m_rbc = m_sample = float("nan")

            max_for_endpoint = solved["p_ips_max_feasible"]
            min_for_endpoint = solved["p_ips_min_feasible"]
            if feasible:
                status = "OK"
                warning = ""
            elif p_ips_request > max_for_endpoint:
                status = "Not feasible"
                warning = (
                    f"Requested IPS fraction too high for {endpoint_name.lower()} endpoint; "
                    f"max is {max_for_endpoint:.6f} ({100*max_for_endpoint:.2f}%) in the selected basis."
                )
            elif p_ips_request < min_for_endpoint:
                status = "Not feasible"
                warning = (
                    f"Requested IPS fraction too low for {endpoint_name.lower()} endpoint; "
                    f"min is {min_for_endpoint:.6f} ({100*min_for_endpoint:.2f}%) in the selected basis."
                )
            else:
                status = "Not feasible"
                warning = "Target cannot be reached with non-negative OptiPrep and PBS fractions."

            rows.append({
                "Condition": condition_index,
                "Endpoint": endpoint_name,
                "Target density [g/ml]": rho_target,
                "Requested IPS fraction": p_ips_request,
                "IPS fraction [%]": 100 * p_ips_request,
                "IPS fraction basis": "carrier excluding RBC suspension" if basis == "carrier" else "final total mixture",
                "Carrier IPS fraction [%]": 100 * solved["carrier_ips_fraction"],
                "Actual final IPS fraction [%]": 100 * solved["p_ips"],
                "Target RBC fraction": rbc_fraction,
                "Target RBC fraction [%]": 100 * rbc_fraction,
                "Sample hematocrit": sample_hematocrit if include_rbc_suspension else float("nan"),
                "Nominal endpoint volume [ml]": endpoint_volume_ml,
                "Preparation volume incl. excess [ml]": prep_volume_ml,
                "Carrier volume [ml]": v_carrier,
                "IPS volume [ml]": v_ips,
                "OptiPrep volume [ml]": v_opt,
                "PBS volume [ml]": v_pbs,
                "RBC suspension volume [ml]": v_sample,
                "RBC volume [ml]": v_rbc,
                f"{sample_liquid_name} volume [ml]": v_sample_liquid,
                "IPS mass [g]": m_ips,
                "OptiPrep mass [g]": m_opt,
                "PBS mass [g]": m_pbs,
                "RBC suspension mass estimate [g]": m_sample,
                f"{sample_liquid_name} mass [g]": m_sample_liquid,
                "RBC mass estimate [g]": m_rbc,
                "Calculated density [g/ml]": solved["rho_check"] if feasible else float("nan"),
                "Feasible IPS min [%]": 100 * min_for_endpoint,
                "Feasible IPS max [%]": 100 * max_for_endpoint,
                "Status": status,
                "Warning": warning,
            })

    total_rows = []
    for name, data in totals.items():
        density = data["density_g_ml"]
        if name == "RBC suspension sample" and data["volume_ml"] > EPS:
            density = data["mass_g"] / data["volume_ml"]
        total_rows.append({
            "Original medium": name,
            "Density [g/ml]": density,
            "Total volume needed [ml]": data["volume_ml"],
            "Total mass needed [g]": data["mass_g"],
        })

    total_nominal = endpoint_volume_ml * 2 * len(ips_fractions)
    total_prepared_requested = prep_volume_ml * 2 * len(ips_fractions)
    feasible_endpoint_count = sum(1 for row in rows if row["Status"] == "OK")
    total_prepared_feasible = prep_volume_ml * feasible_endpoint_count

    summary = {
        "nominal_endpoint_volume_ml": endpoint_volume_ml,
        "preparation_endpoint_volume_ml": prep_volume_ml,
        "excess_fraction": excess_fraction,
        "number_conditions": len(ips_fractions),
        "number_endpoints_requested": 2 * len(ips_fractions),
        "number_endpoints_feasible": feasible_endpoint_count,
        "total_nominal_endpoint_volume_requested_ml": total_nominal,
        "total_preparation_volume_requested_ml": total_prepared_requested,
        "total_preparation_volume_feasible_ml": total_prepared_feasible,
        "combined_feasible_ips_min": combined_p_min,
        "combined_feasible_ips_max": combined_p_max,
        "include_rbc_suspension": 1.0 if include_rbc_suspension else 0.0,
        "sample_hematocrit": sample_hematocrit if include_rbc_suspension else 0.0,
        "rho_sample_liquid_g_ml": rho_sample_liquid if include_rbc_suspension else 0.0,
        "rho_rbc_for_mass_g_ml": rho_rbc_for_mass if include_rbc_suspension else 0.0,
        "ips_fraction_basis_code": 0.0 if basis == "carrier" else 1.0,
    }

    return {
        "recipe_rows": rows,
        "total_rows": total_rows,
        "summary": summary,
        "endpoint_feasible_intervals": endpoint_maxima,
    }