#!/usr/bin/env python3
"""Recompute all numerical outputs for the laser-characterization notebooks.

The script is intentionally separate from notebook presentation.  It writes tidy
CSV/JSON outputs that the notebooks read, while every algorithm remains available
in ``laser_characterization.py`` for interactive use.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from laser_characterization import (
    CameraConfig,
    SpectralConfig,
    TemporalConfig,
    aggregate_spatial_metrics,
    analyze_scope_directory,
    analyze_spatial_directory,
    analyze_spatial_png_directory,
    analyze_spectral_directory,
    find_data_root,
    fit_beam_propagation,
    summarize_replicate_spectra,
    transform_limited_duration,
)

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "outputs"


def _write(frame: pd.DataFrame, name: str) -> Path:
    path = OUTPUT / name
    frame.to_csv(path, index=False)
    print(f"wrote {path.name}: {frame.shape}", flush=True)
    return path


def _fit_row(
    groups: pd.DataFrame,
    measurement: str,
    display_name: str,
    category: str,
    weighting: str,
) -> dict[str, Any]:
    divisor = 2.0 if measurement.startswith("d4sigma") else 1.0
    radius = groups[f"{measurement}_mean"].to_numpy(float) * 1e-6 / divisor
    std = groups[f"{measurement}_std"].to_numpy(float) * 1e-6 / divisor
    use_std = std if weighting == "repeatability_sd" else None
    fit = fit_beam_propagation(
        groups["position_m"].to_numpy(float),
        radius,
        wavelength_nm=532.0,
        radius_std_m=use_std,
    )
    return {
        "measurement": measurement,
        "display_name": display_name,
        "category": category,
        "fit_weighting": weighting,
        **{k: v for k, v in fit.items() if k != "covariance"},
        "strict_unsaturated_z_positions": int((groups["n_strict_unsaturated"] >= 2).sum()),
        "waist_region_contains_primary_saturation": bool(
            groups.loc[groups[measurement + "_mean"].idxmin(), "primary_saturation_any"]
        ),
        "reported_status": "provisional sensitivity analysis; not a standards-grade M2 result",
    }


def spatial_derived(groups: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    focused = groups[groups["focal"].eq("yes")].sort_values("position_m")
    methods = [
        ("d4sigma_x_um", "D4sigma x radius", "second_moment_sensor_axis"),
        ("d4sigma_y_um", "D4sigma y radius", "second_moment_sensor_axis"),
        ("d4sigma_major_um", "D4sigma principal-major radius", "second_moment_principal_axis"),
        ("d4sigma_minor_um", "D4sigma principal-minor radius", "second_moment_principal_axis"),
        ("gaussian_w_major_um", "Gaussian 1/e2 major radius", "model_dependent_gaussian"),
        ("gaussian_w_minor_um", "Gaussian 1/e2 minor radius", "model_dependent_gaussian"),
    ]
    fits: list[dict[str, Any]] = []
    for measurement, display, category in methods:
        for weighting in ("unweighted", "repeatability_sd"):
            try:
                fits.append(_fit_row(focused, measurement, display, category, weighting))
            except Exception as exc:
                fits.append(
                    {
                        "measurement": measurement,
                        "display_name": display,
                        "category": category,
                        "fit_weighting": weighting,
                        "reported_status": f"fit failed: {exc}",
                    }
                )
    fit_table = pd.DataFrame(fits)

    no_lens = groups[groups["focal"].eq("no")].sort_values("position_m")
    slope_rows: list[dict[str, Any]] = []
    for measurement, display, category in methods:
        divisor = 2.0 if measurement.startswith("d4sigma") else 1.0
        z = no_lens["position_m"].to_numpy(float)
        radius = no_lens[f"{measurement}_mean"].to_numpy(float) * 1e-6 / divisor
        slope, intercept = np.polyfit(z, radius, 1)
        predicted = slope * z + intercept
        r2 = 1 - np.sum((radius - predicted) ** 2) / np.sum((radius - np.mean(radius)) ** 2)
        slope_rows.append(
            {
                "measurement": measurement,
                "display_name": display,
                "category": category,
                "apparent_radius_slope_mrad": slope * 1e3,
                "linear_intercept_um": intercept * 1e6,
                "linear_r2": r2,
                "z_min_cm": no_lens["position_cm"].min(),
                "z_max_cm": no_lens["position_cm"].max(),
                "interpretation": "local expansion slope over measured range; not a far-field divergence or M2",
            }
        )
    slopes = pd.DataFrame(slope_rows)

    pointing_rows: list[dict[str, Any]] = []
    for focal, group in groups.groupby("focal"):
        group = group.sort_values("position_m")
        subsets = [("all_positions", group)]
        if focal == "no" and (group["position_cm"] >= 10).sum() >= 3:
            subsets.append(("z_ge_10_cm_diagnostic", group[group["position_cm"] >= 10]))
        for subset_name, subset in subsets:
            row: dict[str, Any] = {
                "focal": focal,
                "subset": subset_name,
                "n_positions": len(subset),
                "within_position_x_sd_median_um": float(np.nanmedian(subset["primary_center_x_mm_std"]) * 1000),
                "within_position_y_sd_median_um": float(np.nanmedian(subset["primary_center_y_mm_std"]) * 1000),
                "between_position_x_peak_to_peak_mm": float(np.ptp(subset["primary_center_x_mm_mean"])),
                "between_position_y_peak_to_peak_mm": float(np.ptp(subset["primary_center_y_mm_mean"])),
                "interpretation": "centroid trajectory across separate z acquisitions; not time-domain pointing jitter",
            }
            for axis in ("x", "y"):
                z = subset["position_m"].to_numpy(float)
                coordinate = subset[f"primary_center_{axis}_mm_mean"].to_numpy(float) * 1e-3
                slope, intercept = np.polyfit(z, coordinate, 1)
                predicted = slope * z + intercept
                denom = np.sum((coordinate - np.mean(coordinate)) ** 2)
                r2 = 1 - np.sum((coordinate - predicted) ** 2) / denom if denom > 0 else np.nan
                row[f"apparent_{axis}_axis_slope_mrad"] = slope * 1e3
                row[f"apparent_{axis}_axis_linear_r2"] = r2
            pointing_rows.append(row)
    pointing = pd.DataFrame(pointing_rows)
    return fit_table, slopes, pointing


def temporal_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    proper = metrics[metrics["termination"].eq("50ohm")].copy()
    summary = (
        proper.sort_values("record_duration_s")
        .groupby("mode", as_index=False, group_keys=False)
        .tail(1)
        .sort_values("mode")
        .copy()
    )
    keep = [
        "mode", "file", "record_duration_us", "sample_interval_ns", "n_cycles",
        "amplitude_v", "fwhm_ns", "rise_10_90_ns", "fall_90_10_ns",
        "equivalent_width_ns", "amplitude_cv_percent", "area_cv_percent",
        "rep_rate_measured_hz", "trigger_edge_count", "trigger_edge_residual_rms_ps",
        "trigger_to_electrical_peak_median_ns", "trigger_to_electrical_peak_std_ns",
        "trigger_delay_repeatability_usable", "pulse_resolution_class",
        "measured_fwhm_interpretation", "absolute_optical_pulse_width_reportable",
        "intrinsic_laser_timing_jitter_reportable", "timing_metric_scope", "absolute_power_calibrated",
        "ch1_imp", "ch1_bw", "ch3_imp", "trigger_source", "trigger_coup",
        "acquisition_mode", "acquisition_realtime", "scope_model", "scope_serial",
    ]
    return summary[[c for c in keep if c in summary.columns]].reset_index(drop=True)


def spectral_derived(metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    repeats = summarize_replicate_spectra(metrics, exclude_near_saturation=True)
    ps = metrics[(metrics["wavelength_family_nm"].eq(1064)) & metrics["mode"].eq("picosecond")].copy()
    ps["transform_limit_gaussian_ps"] = ps["fwhm_frequency_hz"].map(
        lambda value: transform_limited_duration(value, "gaussian") * 1e12
    )
    ps["transform_limit_sech2_ps"] = ps["fwhm_frequency_hz"].map(
        lambda value: transform_limited_duration(value, "sech2") * 1e12
    )
    ps["transform_limit_interpretation"] = (
        "ideal lower bound from spectral intensity FWHM; actual duration requires spectral phase or autocorrelation/FROG"
    )
    return repeats, ps


def integrated_results(
    spatial_frames: pd.DataFrame,
    spatial_groups: pd.DataFrame,
    caustic: pd.DataFrame,
    temporal: pd.DataFrame,
    temporal_modes: pd.DataFrame,
    spectral: pd.DataFrame,
    ps: pd.DataFrame,
) -> dict[str, Any]:
    principal = caustic[
        caustic["category"].eq("second_moment_principal_axis")
        & caustic["fit_weighting"].eq("repeatability_sd")
    ]
    m2_values = principal["m2"].dropna().to_numpy(float)
    narrow = spectral[
        spectral["wavelength_family_nm"].eq(532)
        | spectral["mode"].isin(["cw", "65_ns"])
    ]
    return {
        "spatial": {
            "native_raw_frames": int(len(spatial_frames)),
            "focused_z_positions": int((spatial_groups["focal"] == "yes").sum()),
            "strict_unsaturated_focused_z_positions": int(
                ((spatial_groups["focal"] == "yes") & (spatial_groups["n_strict_unsaturated"] >= 2)).sum()
            ),
            "focused_positions_with_primary_saturation": int(
                ((spatial_groups["focal"] == "yes") & spatial_groups["primary_saturation_any"]).sum()
            ),
            "provisional_principal_axis_m2_range": [float(np.min(m2_values)), float(np.max(m2_values))] if m2_values.size else None,
            "reportable_m2": False,
            "reason": "only eight focused z positions; waist-region saturation; both measured sides extend less than about two fitted Rayleigh ranges; no z-position uncertainty",
        },
        "temporal": {
            "records": int(len(temporal)),
            "proper_50ohm_records": int((temporal["termination"] == "50ohm").sum()),
            "measured_repetition_rate_hz_range_long_records": [
                float(temporal_modes["rep_rate_measured_hz"].min()),
                float(temporal_modes["rep_rate_measured_hz"].max()),
            ],
            "electrical_fwhm_ns_by_mode": {
                row["mode"]: float(row["fwhm_ns"]) for _, row in temporal_modes.iterrows()
            },
            "optical_ps_and_1ns_width_reportable": False,
            "optical_65ns_width_reportable": False,
            "intrinsic_timing_jitter_reportable": False,
            "pulse_energy_ready": False,
        },
        "spectral": {
            "spectra": int(len(spectral)),
            "high_count_warning_spectra": int(spectral["near_saturation"].sum()),
            "narrow_spectra_with_intrinsic_linewidth_reportable": int(narrow["intrinsic_narrow_linewidth_reportable"].astype("boolean").fillna(False).sum()),
            "picosecond_fwhm_nm_range": [float(ps["fwhm_nm"].min()), float(ps["fwhm_nm"].max())],
            "picosecond_gaussian_transform_limit_ps_range": [
                float(ps["transform_limit_gaussian_ps"].min()),
                float(ps["transform_limit_gaussian_ps"].max()),
            ],
            "absolute_1064_wavelength_validated": False,
            "narrow_1059_1061_nm_consistent_with_gsl_nominal_band": bool(
                spectral[
                    spectral["wavelength_family_nm"].eq(1064)
                    & spectral["mode"].isin(["cw", "65_ns"])
                ]["peak_within_gsl_1064_plusminus10_nm"].astype("boolean").fillna(False).all()
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=None)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-legacy-png", action="store_true")
    args = parser.parse_args()

    OUTPUT.mkdir(exist_ok=True)
    data_root = args.data_root.resolve() if args.data_root else find_data_root(ROOT)
    print(f"data root: {data_root}", flush=True)

    inventory_sets = [
        ("native spatial RAW", sorted((data_root / "532nm_Spatial" / "Beam").glob("*.raw"))),
        ("supplemental spatial PNG", sorted(data_root.glob("532nm_*Focal*/*.png"))),
        ("oscilloscope CSV", sorted((data_root / "Pulse_Data").rglob("*.csv"))),
        ("HORIBA spectrum", sorted((data_root / "Spectral_Data").rglob("*.txt"))),
    ]
    inventory_rows = []
    for family, files in inventory_sets:
        inventory_rows.append(
            {
                "family": family,
                "count": len(files),
                "total_bytes": sum(path.stat().st_size for path in files),
                "example": files[0].name if files else "",
            }
        )
    _write(pd.DataFrame(inventory_rows), "data_inventory.csv")

    camera = CameraConfig(sensor_pixel_pitch_um=2.4, magnification=1.0, wavelength_nm=532.0)
    spatial_file = OUTPUT / "spatial_frame_metrics.csv"
    if args.force or not spatial_file.exists():
        spatial_frames, _ = analyze_spatial_directory(data_root / "532nm_Spatial" / "Beam", camera)
        _write(spatial_frames, spatial_file.name)
    else:
        spatial_frames = pd.read_csv(spatial_file)
        print(f"loaded {spatial_file.name}", flush=True)
    spatial_groups = aggregate_spatial_metrics(spatial_frames)
    _write(spatial_groups, "spatial_group_metrics.csv")

    if not args.skip_legacy_png:
        legacy_frames: list[pd.DataFrame] = []
        for folder_name in (
            "532nm_NoFocal_Gain0_Exposure10k_BayerRG12",
            "532nm_YesFocal_Gain0_Exposure10k_BayerRG12",
        ):
            output_name = "legacy_png_" + ("no_focal" if "NoFocal" in folder_name else "focused") + "_metrics.csv"
            output_path = OUTPUT / output_name
            if args.force or not output_path.exists():
                frame, _ = analyze_spatial_png_directory(data_root / folder_name, camera)
                _write(frame, output_name)
            else:
                frame = pd.read_csv(output_path)
                print(f"loaded {output_name}", flush=True)
            legacy_frames.append(frame)
        _write(pd.concat(legacy_frames, ignore_index=True), "legacy_png_spatial_metrics.csv")

    caustic, slopes, pointing = spatial_derived(spatial_groups)
    _write(caustic, "spatial_caustic_fit_sensitivity.csv")
    _write(slopes, "spatial_no_lens_expansion.csv")
    _write(pointing, "spatial_pointing_summary.csv")

    temporal_file = OUTPUT / "temporal_record_metrics.csv"
    if args.force or not temporal_file.exists():
        temporal, _ = analyze_scope_directory(data_root / "Pulse_Data", TemporalConfig())
        _write(temporal, temporal_file.name)
    else:
        temporal = pd.read_csv(temporal_file)
        print(f"loaded {temporal_file.name}", flush=True)
    temporal_modes = temporal_summary(temporal)
    _write(temporal_modes, "temporal_mode_summary.csv")

    power_template = OUTPUT / "average_power_inputs.csv"
    if not power_template.exists():
        _write(
            pd.DataFrame(
                {
                    "mode": ["picosecond", "1_ns", "65_ns"],
                    "average_power_mW": [np.nan, np.nan, np.nan],
                    "average_power_standard_uncertainty_mW_optional": [np.nan, np.nan, np.nan],
                    "notes": ["enter measured value", "enter measured value", "enter measured value"],
                }
            ),
            power_template.name,
        )

    spectral_file = OUTPUT / "spectral_metrics.csv"
    if args.force or not spectral_file.exists():
        spectral, _ = analyze_spectral_directory(data_root / "Spectral_Data", SpectralConfig())
        _write(spectral, spectral_file.name)
    else:
        spectral = pd.read_csv(spectral_file, parse_dates=["acquired_datetime"])
        print(f"loaded {spectral_file.name}", flush=True)
    repeats, ps = spectral_derived(spectral)
    _write(repeats, "spectral_repeatability_sequences.csv")
    _write(ps, "spectral_picosecond_transform_limits.csv")

    results = integrated_results(spatial_frames, spatial_groups, caustic, temporal, temporal_modes, spectral, ps)
    with (OUTPUT / "integrated_results.json").open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    print("wrote integrated_results.json", flush=True)


if __name__ == "__main__":
    main()
