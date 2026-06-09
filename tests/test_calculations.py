import math

from media_mixer.calculations import (
    density_correction_add_stock,
    ips74_osmotic,
    sdm_ips_pbs,
    sdm_with_whole_blood,
    ternary_ips_optiprep_pbs,
    two_endpoint_density_mix,
    paired_gradient_endpoint_batch,
    feasible_fixed_ips_interval_for_density,
    solve_fixed_ips_ternary_fractions_with_sample,
)


def test_sdm_ips_pbs_example():
    r = sdm_ips_pbs(100, 1.100, 1.1116, 1.0047)
    p_ips = r.checks["p_ips"]
    assert math.isclose(p_ips, (1.100 - 1.0047) / (1.1116 - 1.0047), rel_tol=1e-12)
    assert math.isclose(r.apparent_density_g_ml, 1.100, abs_tol=1e-12)


def test_ips74_empirical_sums_to_one():
    r = ips74_osmotic(100, empirical_fixed=True)
    assert math.isclose(sum(c.volume_ml for c in r.components), 100, abs_tol=1e-12)
    assert math.isclose(r.checks["p_percoll"] + r.checks["p_pbs10"] + r.checks["p_hcl"], 1.0, abs_tol=1e-12)


def test_blood_mode_volume_balance():
    r = sdm_with_whole_blood(50, 1.100, 1.1208, 1.0047, 0.05, 0.45, 1.0255)
    assert math.isclose(sum(c.volume_ml for c in r.components), 50, abs_tol=1e-12)
    assert math.isclose(r.checks["v_blood_ml"], (50 * 0.05) / 0.45, rel_tol=1e-12)


def test_two_endpoint_mix():
    r = two_endpoint_density_mix(10, 1.100, 1.010, 1.210)
    assert math.isclose(r.apparent_density_g_ml, 1.100, abs_tol=1e-12)


def test_density_correction():
    r = density_correction_add_stock(10, 1.095, 1.100, 1.120)
    assert r.checks["added_volume_ml"] > 0
    assert math.isclose(r.apparent_density_g_ml, 1.100, abs_tol=1e-12)


def test_ternary():
    r = ternary_ips_optiprep_pbs(10, 1.100, 0.45, 1.1116, 1.215, 1.0047)
    assert math.isclose(r.apparent_density_g_ml, 1.100, abs_tol=1e-12)
    assert math.isclose(r.checks["p_ips"] + r.checks["p_optiprep"] + r.checks["p_pbs"], 1.0, abs_tol=1e-12)


def test_paired_endpoint_batch_flags_too_high_ips():
    batch = paired_gradient_endpoint_batch(
        endpoint_volume_ml=6,
        ips_fractions=[0.50, 0.60, 0.75, 0.87],
        rho_low_target=1.095,
        rho_high_target=1.109,
        rho_ips=1.1189,
        rho_optiprep=1.215,
        rho_pbs=1.0047,
        excess_fraction=0.0,
    )
    rows = batch["recipe_rows"]
    assert len(rows) == 8
    bad = [row for row in rows if row["Status"] != "OK"]
    assert len(bad) == 1
    assert bad[0]["Condition"] == 4
    assert bad[0]["Endpoint"] == "Low"
    assert "too high" in bad[0]["Warning"]
    assert math.isclose(batch["summary"]["combined_feasible_ips_max"], (1.095 - 1.0047) / (1.1189 - 1.0047), rel_tol=1e-10)


def test_fixed_ips_feasible_interval_low_density():
    interval = feasible_fixed_ips_interval_for_density(1.095, 1.1189, 1.215, 1.0047)
    assert interval["feasible"] == 1.0
    assert interval["p_max"] < 0.80


def test_paired_endpoint_batch_with_rbc_suspension_balances_liquid_density():
    batch = paired_gradient_endpoint_batch(
        endpoint_volume_ml=6,
        ips_fractions=[0.50, 0.60],
        rho_low_target=1.095,
        rho_high_target=1.109,
        rho_ips=1.1189,
        rho_optiprep=1.215,
        rho_pbs=1.0047,
        excess_fraction=0.0,
        include_rbc_suspension=True,
        target_rbc_fractions=[0.03],
        sample_hematocrit=0.45,
        rho_sample_liquid=1.0047,
        sample_liquid_name="sample PBS",
        rho_rbc_for_mass=1.100,
        ips_fraction_basis="carrier_excluding_sample",
    )
    rows = batch["recipe_rows"]
    assert len(rows) == 4
    assert all(row["Status"] == "OK" for row in rows)
    for row in rows:
        assert math.isclose(row["RBC suspension volume [ml]"], (0.03 / 0.45) * 6, rel_tol=1e-12)
        assert math.isclose(row["RBC volume [ml]"], 0.03 * 6, rel_tol=1e-12)
        assert math.isclose(row["Calculated density [g/ml]"], row["Target density [g/ml]"], abs_tol=1e-12)


def test_solve_fixed_ips_with_rbc_suspension_can_flag_too_high_final_fraction():
    solved = solve_fixed_ips_ternary_fractions_with_sample(
        rho_target=1.095,
        ips_fraction=0.87,
        rho_ips=1.1189,
        rho_optiprep=1.215,
        rho_pbs=1.0047,
        target_rbc_fraction=0.03,
        sample_hematocrit=0.45,
        rho_sample_liquid=1.0047,
        ips_fraction_basis="final_total_mixture",
    )
    assert solved["feasible"] == 0.0
    assert solved["p_ips_max_feasible"] < 0.87
