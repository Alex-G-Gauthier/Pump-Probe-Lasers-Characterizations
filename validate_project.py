"""Validate parsers, derived outputs, numerical regressions, figures, and notebooks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import nbformat
import numpy as np
import pandas as pd

from laser_characterization import (
    CameraConfig,
    SpectralConfig,
    TemporalConfig,
    analyze_scope_record,
    analyze_spectrum,
    find_data_root,
    green_superpixel,
    load_bayer_raw,
)


EXPECTED_NOTEBOOKS = [
    "00_Start_Here",
    "01_Spatial_532nm",
    "02_Temporal_1064nm",
    "03_Spectral_532_1064nm",
    "04_Integrated_Summary",
]

EXPECTED_FIGURES = [
    "spatial_native_contact_sheet_v2.png",
    "spatial_size_ellipticity_v2.png",
    "spatial_focused_caustic_v2.png",
    "spatial_provisional_m2_v2.png",
    "spatial_centroid_trajectory_v2.png",
    "spatial_legacy_png_v2.png",
    "temporal_long_records_v2.png",
    "temporal_coherent_average_v2.png",
    "temporal_termination_v2.png",
    "spectral_532_gratings_v2.png",
    "spectral_1064_slit_broadening_v2.png",
    "spectral_65ns_sequences_v2.png",
    "spectral_1064_modes_v2.png",
    "spectral_transform_limits_v2.png",
]

EXPECTED_OUTPUT_ROWS = {
    "data_inventory.csv": 4,
    "spatial_frame_metrics.csv": 36,
    "spatial_group_metrics.csv": 12,
    "legacy_png_no_focal_metrics.csv": 12,
    "legacy_png_focused_metrics.csv": 8,
    "legacy_png_spatial_metrics.csv": 20,
    "spatial_caustic_fit_sensitivity.csv": 12,
    "spatial_no_lens_expansion.csv": 6,
    "spatial_pointing_summary.csv": 3,
    "temporal_record_metrics.csv": 19,
    "temporal_mode_summary.csv": 3,
    "average_power_inputs.csv": 3,
    "spectral_metrics.csv": 40,
    "spectral_repeatability_sequences.csv": 21,
    "spectral_picosecond_transform_limits.csv": 10,
}


def _assert_close(actual: float, expected: float, atol: float, label: str) -> None:
    assert np.isfinite(actual), f"{label} is not finite: {actual}"
    assert abs(actual - expected) <= atol, (
        f"{label} changed: expected {expected} +/- {atol}, found {actual}"
    )


def _validate_notebook(path: Path, require_execution: bool) -> None:
    notebook = nbformat.read(path, as_version=4)
    code_cells = [cell for cell in notebook.cells if cell.cell_type == "code"]
    assert code_cells, f"No code cells found in {path.name}"
    for index, cell in enumerate(notebook.cells):
        if cell.cell_type != "code":
            continue
        errors = [output for output in cell.get("outputs", []) if output.get("output_type") == "error"]
        assert not errors, f"Error output in cell {index} of {path.name}: {errors}"
        if require_execution:
            assert cell.get("execution_count") is not None, f"Unexecuted cell {index} in {path.name}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help="Path to the extracted '532 and 1064' measurement folder",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    data = args.data_root.expanduser().resolve() if args.data_root else find_data_root(root)

    raw_files = sorted((data / "532nm_Spatial" / "Beam").glob("*.raw"))
    png_files = sorted(data.glob("532nm_*Focal*/*.png"))
    scope_files = sorted((data / "Pulse_Data").rglob("*.csv"))
    spectral_files = sorted((data / "Spectral_Data").rglob("*.txt"))

    assert len(raw_files) == 36, f"Expected 36 RAW files, found {len(raw_files)}"
    assert len(png_files) == 22, f"Expected 22 supplemental PNG files, found {len(png_files)}"
    assert len(scope_files) == 19, f"Expected 19 scope CSVs, found {len(scope_files)}"
    assert len(spectral_files) == 40, f"Expected 40 spectra, found {len(spectral_files)}"

    raw = load_bayer_raw(raw_files[0], CameraConfig())
    assert raw.shape == (2064, 3088)
    assert raw.dtype.kind == "u"
    green = green_superpixel(raw)
    assert green.shape == (1032, 1544)

    temporal, _ = analyze_scope_record(scope_files[0], TemporalConfig(), keep_waveform=False)
    assert temporal["sample_count"] > 100
    assert temporal["record_duration_s"] > 0

    spectrum, arrays = analyze_spectrum(spectral_files[0], SpectralConfig(), keep_arrays=True)
    assert spectrum["sample_count"] > 10
    assert spectrum["fwhm_nm"] > 0
    assert len(arrays["wavelength_nm"]) == spectrum["sample_count"]

    loaded: dict[str, pd.DataFrame] = {}
    for name, expected_rows in EXPECTED_OUTPUT_ROWS.items():
        path = root / "outputs" / name
        assert path.exists(), f"Missing output: {path}"
        frame = pd.read_csv(path)
        loaded[name] = frame
        assert len(frame) == expected_rows, (
            f"Unexpected row count in {name}: expected {expected_rows}, found {len(frame)}"
        )

    inventory = loaded["data_inventory.csv"].set_index("family")
    expected_inventory = {
        "native spatial RAW": 36,
        "supplemental spatial PNG": 22,
        "oscilloscope CSV": 19,
        "HORIBA spectrum": 40,
    }
    for family, count in expected_inventory.items():
        assert int(inventory.loc[family, "count"]) == count, f"Inventory mismatch for {family}"

    spatial_groups = loaded["spatial_group_metrics.csv"]
    assert int(spatial_groups["focal"].eq("yes").sum()) == 8
    assert int(spatial_groups["focal"].eq("no").sum()) == 4
    focused = spatial_groups[spatial_groups["focal"].eq("yes")]
    assert int(focused["primary_saturation_any"].astype(bool).sum()) == 5
    assert int((focused["n_strict_unsaturated"] >= 2).sum()) == 3

    caustic = loaded["spatial_caustic_fit_sensitivity.csv"]
    reportable = caustic["reportable_m2"].astype(str).str.lower().map({"true": True, "false": False})
    assert reportable.notna().all() and not reportable.any(), "No M2 fit should be reportable"
    principal = caustic[
        caustic["category"].eq("second_moment_principal_axis")
        & caustic["fit_weighting"].eq("repeatability_sd")
    ].sort_values("measurement")
    assert len(principal) == 2
    _assert_close(float(principal["m2"].min()), 1.8142842485, 1e-6, "provisional M2 minimum")
    _assert_close(float(principal["m2"].max()), 2.5386708269, 1e-6, "provisional M2 maximum")

    modes = loaded["temporal_mode_summary.csv"].set_index("mode")
    assert set(modes.index) == {"picosecond", "1_ns", "65_ns"}
    assert modes["rep_rate_measured_hz"].between(199999.0, 200001.0).all()
    _assert_close(float(modes.loc["picosecond", "fwhm_ns"]), 20.8895714452, 1e-4, "picosecond electrical FWHM")
    _assert_close(float(modes.loc["1_ns", "fwhm_ns"]), 32.3680398134, 1e-4, "1 ns electrical FWHM")
    _assert_close(float(modes.loc["65_ns", "fwhm_ns"]), 49.0249562247, 1e-4, "65 ns electrical FWHM")
    assert modes["timing_metric_scope"].notna().all()
    optical_width = modes["absolute_optical_pulse_width_reportable"].astype(str).str.lower()
    assert optical_width.eq("false").all()

    spectra = loaded["spectral_metrics.csv"]
    spectra_532 = spectra[spectra["wavelength_family_nm"].eq(532)]
    assert len(spectra_532) == 4
    assert (spectra_532["points_per_fwhm"] < 3.0).all()
    assert not spectra_532["intrinsic_narrow_linewidth_reportable"].astype(bool).any()
    ps = loaded["spectral_picosecond_transform_limits.csv"]
    assert ps["fwhm_nm"].between(8.0, 10.0).all()
    assert ps["transform_limit_gaussian_ps"].between(0.16, 0.20).all()

    integrated_path = root / "outputs" / "integrated_results.json"
    assert integrated_path.exists(), f"Missing output: {integrated_path}"
    integrated = json.loads(integrated_path.read_text(encoding="utf-8"))
    assert integrated["spatial"]["reportable_m2"] is False
    assert integrated["temporal"]["intrinsic_timing_jitter_reportable"] is False
    assert integrated["spectral"]["absolute_1064_wavelength_validated"] is False

    for stem in EXPECTED_NOTEBOOKS:
        clean = root / f"{stem}.ipynb"
        executed = root / f"{stem}.executed.ipynb"
        assert clean.exists(), f"Missing clean notebook: {clean.name}"
        assert executed.exists(), f"Missing executed notebook: {executed.name}"
        _validate_notebook(clean, require_execution=False)
        _validate_notebook(executed, require_execution=True)

    for name in EXPECTED_FIGURES:
        path = root / "figures" / name
        assert path.exists() and path.stat().st_size > 10_000, f"Missing or empty figure: {path}"

    print("Validation passed")
    print(f"Data root: {data}")
    print(
        "Files: "
        f"{len(raw_files)} RAW, {len(png_files)} supplemental PNG, "
        f"{len(scope_files)} scope CSV, {len(spectral_files)} spectra"
    )
    print(f"Outputs: {len(EXPECTED_OUTPUT_ROWS)} tabular products + integrated_results.json")
    print(f"Notebooks: {len(EXPECTED_NOTEBOOKS)} clean + {len(EXPECTED_NOTEBOOKS)} executed; no error outputs")
    print(f"Figures: {len(EXPECTED_FIGURES)} publication-resolution PNG files")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"Validation failed: {exc}", file=sys.stderr)
        raise
