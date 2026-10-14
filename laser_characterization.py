"""Enhanced laser-characterization analysis for the supplied 532/1064 nm dataset.

Version 2 keeps the validated numerical core from ``laser_characterization_core``
but adds experiment-specific provenance, conservative claim-status fields, legacy
PNG handling, trigger-channel metrology, the correct legacy PDA36A 0 dB settings,
and wavelength-calibration caveats appropriate to the actual experiment.

The guiding rule is that a numerical value and a physically supportable claim are
not the same thing.  Tables therefore carry explicit status fields such as
``provisional_for_caustic``, ``pulse_resolution_class``, and ``width_status``.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import csv
import math
import re

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage

import laser_characterization_core as _core
from laser_characterization_core import *  # re-export the validated v1 implementation


# -----------------------------------------------------------------------------
# Experiment metadata and claim vocabulary
# -----------------------------------------------------------------------------

CLAIM_LEVELS = {
    "supported": "Supported by this dataset and measurement chain.",
    "supported_with_caveat": "Numerically supported, with a stated instrumental or calibration caveat.",
    "provisional": "Useful for model comparison, but not a standards-grade reported result.",
    "not_supported": "The supplied measurement chain cannot support the physical claim.",
}


# -----------------------------------------------------------------------------
# Spatial extensions
# -----------------------------------------------------------------------------

CameraConfig = _core.CameraConfig


def _raw_saturation_map(raw: np.ndarray, threshold_dn: float) -> np.ndarray:
    """Map native Bayer saturation to one Boolean value per 2 x 2 Bayer cell."""
    return (
        (raw[0::2, 0::2] >= threshold_dn)
        | (raw[0::2, 1::2] >= threshold_dn)
        | (raw[1::2, 0::2] >= threshold_dn)
        | (raw[1::2, 1::2] >= threshold_dn)
    )


def analyze_spatial_frame(
    path: str | Path,
    config: CameraConfig = CameraConfig(),
    center_hint: tuple[float, float] | None = None,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Analyze a native RAW frame and add experiment-specific quality semantics.

    Secondary reflections are retained as an image-quality metric but are not, by
    themselves, allowed to invalidate the highest-intensity primary lobe specified
    by the experimenter.  Saturation *inside* the selected primary aperture remains
    a strict exclusion for an unsaturated second-moment caustic.
    """
    metrics, arrays = _core.analyze_spatial_frame(path, config=config, center_hint=center_hint)
    raw = arrays["raw"]
    primary_mask = arrays["primary_mask"]
    threshold = 0.995 * config.full_scale_dn
    sat_map = _raw_saturation_map(raw, threshold)
    sat_total = int(np.count_nonzero(sat_map))
    sat_primary = int(np.count_nonzero(sat_map & primary_mask))
    primary_pixels = int(np.count_nonzero(primary_mask))
    primary_sat_fraction = sat_primary / primary_pixels if primary_pixels else np.nan

    # The Gaussian fit is used to define the physical lobe selected by the user.
    # It is not used as evidence that the beam itself is Gaussian.
    primary_lobe_usable = bool(
        (not metrics["clipped"])
        and np.isfinite(metrics["d4sigma_major_um"])
        and np.isfinite(metrics["d4sigma_minor_um"])
        and bool(metrics.get("success", True))
        and metrics["gaussian_fit_nrmse"] <= 0.20
    )
    strict_unsaturated = bool(primary_lobe_usable and sat_primary == 0)
    provisional = bool(primary_lobe_usable)
    reasons: list[str] = []
    if metrics["clipped"]:
        reasons.append("selected primary aperture intersects the sensor boundary")
    if metrics["gaussian_fit_nrmse"] > 0.20 or not bool(metrics.get("success", True)):
        reasons.append("primary-lobe localization fit is poor")
    if sat_primary:
        reasons.append(f"{sat_primary} saturated Bayer cells fall inside the primary aperture")
    if metrics["secondary_structure"]:
        reasons.append("secondary reflections/structure are present outside the selected lobe")
    if not reasons:
        reasons.append("no automated exclusion flag")

    metrics.update(
        {
            "data_source": "native_12bit_bayer_raw",
            "saturated_bayer_cells_total": sat_total,
            "saturated_bayer_cells_primary": sat_primary,
            "primary_saturation_fraction": float(primary_sat_fraction),
            "primary_saturation_present": bool(sat_primary > 0),
            "primary_lobe_usable": primary_lobe_usable,
            "strict_unsaturated_candidate": strict_unsaturated,
            "provisional_for_caustic": provisional,
            # v2 strict validity no longer rejects a physically identified primary
            # lobe merely because a ghost exists elsewhere in the field.
            "valid_for_m2": strict_unsaturated,
            "m2_exclusion_reasons": "; ".join(reasons),
            "spatial_claim_level": (
                "supported_with_caveat" if strict_unsaturated else ("provisional" if provisional else "not_supported")
            ),
        }
    )
    arrays["saturation_map_green"] = sat_map
    return metrics, arrays


def analyze_spatial_directory(
    directory: str | Path,
    config: CameraConfig = CameraConfig(),
    center_hints: Mapping[str, tuple[float, float]] | None = None,
    return_diagnostics_for: Iterable[str] | None = None,
    group_consistent_lobe: bool = True,
) -> tuple[pd.DataFrame, dict[str, dict[str, np.ndarray]]]:
    """Analyze all native RAW frames using a common dominant-lobe hint per z group."""
    directory = Path(directory)
    files = sorted(directory.glob("*.raw"))
    if not files:
        raise FileNotFoundError(f"No .raw files found in {directory}")
    center_hints = center_hints or {}
    diagnostic_names = set(return_diagnostics_for or [])
    parsed = {path: _core.parse_spatial_filename(path) for path in files}
    groups: dict[tuple[str, float], list[Path]] = {}
    for path, meta in parsed.items():
        groups.setdefault((meta["focal"], meta["position_cm"]), []).append(path)

    group_hints: dict[tuple[str, float], tuple[float, float]] = {}
    if group_consistent_lobe:
        for key, group_files in groups.items():
            group_key = f"{key[0]}:{key[1]:g}"
            if group_key in center_hints:
                group_hints[key] = center_hints[group_key]
                continue
            explicit = [center_hints[p.name] for p in group_files if p.name in center_hints]
            if explicit:
                group_hints[key] = tuple(np.median(np.asarray(explicit, dtype=float), axis=0))
                continue
            frames = [_core.green_superpixel(_core.load_bayer_raw(p, config), config.bayer_pattern) for p in group_files]
            median_green = np.median(np.stack(frames, axis=0), axis=0)
            plane, noise, _ = _core.estimate_background_plane(median_green, border_px=config.border_px_green)
            corrected = np.clip(median_green - plane, 0, None)
            _, seed, _ = _core.locate_primary_roi(corrected, noise, config)
            group_hints[key] = seed

    rows: list[dict[str, Any]] = []
    diagnostics: dict[str, dict[str, np.ndarray]] = {}
    for path in files:
        meta = parsed[path]
        key = (meta["focal"], meta["position_cm"])
        group_key = f"{key[0]}:{key[1]:g}"
        hint = center_hints.get(path.name, center_hints.get(group_key, group_hints.get(key)))
        result, arrays = analyze_spatial_frame(path, config=config, center_hint=hint)
        rows.append(result)
        if path.name in diagnostic_names:
            diagnostics[path.name] = arrays
    frame = pd.DataFrame(rows).sort_values(["focal", "position_cm", "replicate"]).reset_index(drop=True)
    return frame, diagnostics


def aggregate_spatial_metrics(frame_metrics: pd.DataFrame) -> pd.DataFrame:
    """Aggregate independent frames; standard deviations are Type-A repeatability only."""
    base = _core.aggregate_spatial_metrics(frame_metrics)
    extras: list[dict[str, Any]] = []
    for keys, group in frame_metrics.groupby(["focal", "position_cm", "position_m"], dropna=False):
        extras.append(
            {
                "focal": keys[0],
                "position_cm": keys[1],
                "position_m": keys[2],
                "n_primary_lobe_usable": int(group["primary_lobe_usable"].sum()),
                "n_provisional_caustic": int(group["provisional_for_caustic"].sum()),
                "n_strict_unsaturated": int(group["strict_unsaturated_candidate"].sum()),
                "primary_saturation_any": bool(group["primary_saturation_present"].any()),
                "saturated_bayer_cells_primary_max": int(group["saturated_bayer_cells_primary"].max()),
                "primary_saturation_fraction_max": float(group["primary_saturation_fraction"].max()),
                "primary_ellipticity_mean": float(group["primary_ellipticity"].mean()),
                "primary_ellipticity_std": float(group["primary_ellipticity"].std(ddof=1)) if len(group) > 1 else np.nan,
                "primary_angle_deg_mean": float(group["primary_angle_deg"].mean()),
                "primary_angle_deg_std": float(group["primary_angle_deg"].std(ddof=1)) if len(group) > 1 else np.nan,
                "full_field_d4sigma_major_um_mean": float(group["full_field_d4sigma_major_um"].mean()),
                "full_field_d4sigma_minor_um_mean": float(group["full_field_d4sigma_minor_um"].mean()),
                "repeatability_note": "SD across independently acquired frames; not total measurement uncertainty",
            }
        )
    extra = pd.DataFrame(extras)
    merged = base.merge(extra, on=["focal", "position_cm", "position_m"], how="left")
    return merged.sort_values(["focal", "position_cm"]).reset_index(drop=True)


def parse_spatial_png_filename(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    match = re.search(r"Position_([-+]?\d+(?:\.\d+)?)mm", path.name, re.I)
    if not match:
        raise ValueError(f"Cannot parse position from {path.name}")
    z_mm = float(match.group(1))
    parent = path.parent.name.lower()
    focal = "yes" if "yesfocal" in parent else ("no" if "nofocal" in parent else "unknown")
    return {
        "file": path.name,
        "path": str(path),
        "focal": focal,
        "position_mm": z_mm,
        "position_cm": z_mm / 10.0,
        "position_m": z_mm / 1000.0,
        "replicate": 0,
        "data_source": "legacy_8bit_rgb_png",
    }


def _png_green_native(path: str | Path) -> np.ndarray:
    arr = np.asarray(Image.open(path))
    if arr.ndim == 2:
        return arr.astype(float)
    if arr.ndim == 3 and arr.shape[2] >= 2:
        return arr[..., 1].astype(float)
    raise ValueError(f"Unsupported PNG shape {arr.shape}")


def analyze_spatial_png_frame(
    path: str | Path,
    background_path: str | Path,
    config: CameraConfig = CameraConfig(),
    center_hint: tuple[float, float] | None = None,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Analyze one supplemental 8-bit PNG without mixing it with native RAW data."""
    path = Path(path)
    green_raw = _png_green_native(path)
    background = _png_green_native(background_path)
    if green_raw.shape != background.shape:
        raise ValueError("PNG frame and background-reference dimensions differ")
    bg_subtracted = np.clip(green_raw - background, 0, None)
    plane, noise_sigma, bg_coeff = _core.estimate_background_plane(
        bg_subtracted, border_px=config.border_px_green
    )
    corrected = np.clip(bg_subtracted - plane, 0, None)
    roi_slice, seed, roi_info = _core.locate_primary_roi(corrected, noise_sigma, config, center_hint=center_hint)
    roi = corrected[roi_slice]
    yoff = roi_slice[0].start or 0
    xoff = roi_slice[1].start or 0
    fit = _core.fit_rotated_gaussian(
        roi, noise_sigma=noise_sigma, center_hint_roi=(seed[0] - xoff, seed[1] - yoff)
    )
    x0 = xoff + fit["x0_roi_px"]
    y0 = yoff + fit["y0_roi_px"]
    wx, wy, theta = fit["wx_1e2_px"], fit["wy_1e2_px"], fit["theta_rad"]
    yy, xx = np.indices(corrected.shape, dtype=float)
    ct, st = np.cos(theta), np.sin(theta)
    dx, dy = xx - x0, yy - y0
    xp = ct * dx + st * dy
    yp = -st * dx + ct * dy
    rho2 = (xp / max(config.primary_aperture_scale * wx, 1e-9)) ** 2 + (
        yp / max(config.primary_aperture_scale * wy, 1e-9)
    ) ** 2
    primary_mask = rho2 <= 1.0
    floor = max(3.0 * noise_sigma, 0.0005 * corrected.max(initial=0.0))
    intensity = np.where(corrected > floor, corrected, 0.0)
    primary = _core.second_moments(intensity, mask=primary_mask)
    full = _core.second_moments(intensity)
    total = float(intensity.sum())
    primary_signal = float(intensity[primary_mask].sum())
    outside = ndimage.gaussian_filter(corrected, 1.5).copy()
    outside[primary_mask] = 0
    secondary_ratio = float(outside.max(initial=0.0) / max(corrected.max(initial=0.0), 1e-15))

    x_extent = config.primary_aperture_scale * math.sqrt((wx * ct) ** 2 + (wy * st) ** 2)
    y_extent = config.primary_aperture_scale * math.sqrt((wx * st) ** 2 + (wy * ct) ** 2)
    clipped = bool(
        x0 - x_extent <= 3
        or x0 + x_extent >= corrected.shape[1] - 4
        or y0 - y_extent <= 3
        or y0 + y_extent >= corrected.shape[0] - 4
    )
    sat_map = green_raw >= 255
    sat_primary = int(np.count_nonzero(sat_map & primary_mask))
    sat_total = int(np.count_nonzero(sat_map))
    pitch_um = config.sensor_pixel_pitch_um / config.magnification
    usable = bool(not clipped and fit["success"] and fit["fit_nrmse"] <= 0.20)
    meta = parse_spatial_png_filename(path)
    metrics: dict[str, Any] = {
        **meta,
        **{f"background_{k}": float(v) for k, v in bg_coeff.items()},
        "background_noise_dn": float(noise_sigma),
        "raw_min_dn": float(green_raw.min()),
        "raw_max_dn": float(green_raw.max()),
        "saturated_pixels_total": sat_total,
        "saturated_pixels_primary": sat_primary,
        "saturation_fraction": float(np.mean(sat_map)),
        "saturated": bool(sat_total > 0),
        "primary_saturation_present": bool(sat_primary > 0),
        "primary_center_x_px": float(x0),
        "primary_center_y_px": float(y0),
        "primary_center_x_mm": float(x0 * pitch_um / 1000),
        "primary_center_y_mm": float(y0 * pitch_um / 1000),
        "gaussian_w_major_um": float(wx * pitch_um),
        "gaussian_w_minor_um": float(wy * pitch_um),
        "gaussian_theta_deg": float(fit["theta_deg"]),
        "gaussian_fit_nrmse": float(fit["fit_nrmse"]),
        "d4sigma_x_um": float(primary["d4sigma_x"] * pitch_um),
        "d4sigma_y_um": float(primary["d4sigma_y"] * pitch_um),
        "d4sigma_major_um": float(primary["d4sigma_major"] * pitch_um),
        "d4sigma_minor_um": float(primary["d4sigma_minor"] * pitch_um),
        "primary_ellipticity": float(primary["ellipticity_major_minor"]),
        "primary_angle_deg": float(primary["major_axis_angle_deg"]),
        "full_field_d4sigma_major_um": float(full["d4sigma_major"] * pitch_um),
        "full_field_d4sigma_minor_um": float(full["d4sigma_minor"] * pitch_um),
        "primary_fraction": primary_signal / total if total > 0 else np.nan,
        "secondary_peak_ratio": secondary_ratio,
        "clipped": clipped,
        "primary_lobe_usable": usable,
        "provisional_for_caustic": bool(usable and sat_primary == 0),
        "valid_for_m2": False,
        "format_limitation": "8-bit RGB PNG is demosaiced/processed and is not equivalent to native BayerRG12",
        **roi_info,
    }
    arrays = {
        "green_raw": green_raw,
        "background": background,
        "corrected": corrected,
        "primary_mask": primary_mask,
        "saturation_map": sat_map,
    }
    return metrics, arrays


def analyze_spatial_png_directory(
    directory: str | Path,
    config: CameraConfig = CameraConfig(),
    return_diagnostics_for: Iterable[str] | None = None,
) -> tuple[pd.DataFrame, dict[str, dict[str, np.ndarray]]]:
    directory = Path(directory)
    background_path = directory / "Background_Reference.png"
    if not background_path.exists():
        raise FileNotFoundError(f"Missing {background_path}")
    files = sorted(
        [p for p in directory.glob("Position_*.png")],
        key=lambda p: parse_spatial_png_filename(p)["position_mm"],
    )
    keep = set(return_diagnostics_for or [])
    rows: list[dict[str, Any]] = []
    diagnostics: dict[str, dict[str, np.ndarray]] = {}
    for path in files:
        result, arrays = analyze_spatial_png_frame(path, background_path, config=config)
        rows.append(result)
        if path.name in keep:
            diagnostics[path.name] = arrays
    return pd.DataFrame(rows).sort_values("position_mm").reset_index(drop=True), diagnostics


def fit_beam_propagation(
    z_m: Sequence[float],
    radius_m: Sequence[float],
    wavelength_nm: float,
    radius_std_m: Sequence[float] | None = None,
    enforce_m2_ge_1: bool = True,
) -> dict[str, Any]:
    """Fit a caustic and attach a conservative, non-ISO coverage diagnostic."""
    result = _core.fit_beam_propagation(
        z_m,
        radius_m,
        wavelength_nm,
        radius_std_m=radius_std_m,
        enforce_m2_ge_1=enforce_m2_ge_1,
    )
    z = np.asarray(z_m, dtype=float)
    w = np.asarray(radius_m, dtype=float)
    valid = np.isfinite(z) & np.isfinite(w) & (w > 0)
    z = z[valid]
    z0 = result["z0_m"]
    zr = result["rayleigh_range_m"]
    left_span = (z0 - np.min(z)) / zr if zr > 0 else np.nan
    right_span = (np.max(z) - z0) / zr if zr > 0 else np.nan
    n_near = int(np.sum(np.abs(z - z0) <= zr)) if zr > 0 else 0
    # Deliberately conservative heuristic. It is a project QC rule, not a claim
    # that ISO 11146 compliance has been established.
    coverage_pass = bool(
        result["n_points"] >= 10
        and result["n_left_of_waist"] >= 3
        and result["n_right_of_waist"] >= 3
        and left_span >= 2.0
        and right_span >= 2.0
        and n_near >= 4
    )
    if coverage_pass:
        grade = "coverage_heuristic_passes; other systematic checks still required"
    elif result["well_constrained"]:
        grade = "two-sided numerical fit, but coverage is provisional"
    else:
        grade = "underconstrained caustic"
    result.update(
        {
            "left_span_in_rayleigh_ranges": float(left_span),
            "right_span_in_rayleigh_ranges": float(right_span),
            "n_points_within_one_rayleigh_range": n_near,
            "coverage_heuristic_pass": coverage_pass,
            "coverage_grade": grade,
            "fit_uncertainty_scope": "curve-fit covariance only; excludes z scale, pixel scale, background, saturation and model selection",
            "reportable_m2": False if not coverage_pass else True,
        }
    )
    return result


# -----------------------------------------------------------------------------
# Temporal extensions
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class TemporalConfig:
    expected_rep_rate_hz: float = 200_000.0
    detector_bandwidth_hz: float | None = 10e6
    scope_bandwidth_hz: float | None = 1e9
    align_cycles: bool = True
    minimum_single_cycle_snr_for_alignment: float = 5.0
    minimum_cycles_for_folding: int = 4
    local_alignment_fraction_of_period: float = 0.04
    positive_polarity: bool | None = None
    detector_responsivity_a_per_w: float | None = None
    # Legacy PDA36A, 0 dB, Rev E specifications. 50 ohm loading halves the
    # high-impedance voltage gain because of the detector's 50 ohm series output.
    transimpedance_v_per_a_50ohm: float | None = 0.75e3
    transimpedance_v_per_a_high_z: float | None = 1.51e3
    optical_path_transmission: float = 1.0
    user_reported_ch1_coupling: str = "DC"


def parse_scope_setup(path: str | Path) -> dict[str, Any]:
    """Extract normalized DSO6104A channel, trigger and acquisition settings."""
    path = Path(path)
    rows: list[tuple[str, str]] = []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        for i, row in enumerate(reader):
            if i > 100:
                break
            row = row + [""] * max(0, 2 - len(row))
            rows.append((row[0].strip().rstrip(":"), row[1].strip()))

    result: dict[str, Any] = {}
    # Instrument identity.
    for key, value in rows[:12]:
        if key in {"Model", "Serial Number", "Firmware Version", "Start Time", "Analog Sample Count"}:
            result[key.lower().replace(" ", "_")] = value

    headings: list[tuple[int, str]] = []
    for i, (key, value) in enumerate(rows):
        if key.startswith("Channel 1"):
            headings.append((i, "ch1"))
        elif key.startswith("Channel 3"):
            headings.append((i, "ch3"))
        elif key == "TRIGGER":
            headings.append((i, "trigger"))
        elif key == "HORIZONTAL":
            headings.append((i, "horizontal"))
        elif key == "ACQUISITION":
            headings.append((i, "acquisition"))
    headings.append((len(rows), "end"))
    for (start, label), (stop, _) in zip(headings[:-1], headings[1:]):
        block: dict[str, str] = {}
        for key, value in rows[start + 1 : stop]:
            if key and value:
                block[key.lower().replace(" ", "_")] = value
        for key, value in block.items():
            result[f"{label}_{key}"] = value
    return result


def trigger_edge_metrics(
    time_s: Sequence[float],
    trigger_v: Sequence[float],
    optical_v: Sequence[float] | None = None,
    expected_rep_rate_hz: float = 200_000.0,
) -> dict[str, Any]:
    """Measure CH3 edge timing and optional trigger-to-electrical-peak repeatability.

    The implementation is vectorized over waveform samples. A high-threshold
    crossing proposes each rising edge, a midpoint crossing supplies the
    interpolated timestamp, and a 0.70-period refractory interval rejects the
    multiple crossings caused by ringing on a single electrical edge.
    """
    t = np.asarray(time_s, dtype=float)
    trig = np.asarray(trigger_v, dtype=float)
    valid = np.isfinite(t) & np.isfinite(trig)
    t, trig = t[valid], trig[valid]
    if t.size < 10:
        return {"trigger_edge_count": 0}

    low, high = np.percentile(trig, [10, 90])
    span = high - low
    threshold = 0.5 * (low + high)
    low_threshold = low + 0.30 * span
    high_threshold = low + 0.70 * span
    smooth = ndimage.uniform_filter1d(trig, size=9, mode="nearest")

    # Candidate rising edges. The high threshold implements the Schmitt
    # trigger's arming margin; the refractory interval rejects edge ringing.
    high_crossings = np.flatnonzero(
        (smooth[:-1] < high_threshold) & (smooth[1:] >= high_threshold)
    ) + 1
    midpoint_crossings = np.flatnonzero(
        (smooth[:-1] < threshold) & (smooth[1:] >= threshold)
    )
    low_samples = np.flatnonzero(smooth <= low_threshold)

    minimum_separation = 0.70 / expected_rep_rate_hz
    crossing: list[float] = []
    last_accepted_index = -1
    for j in high_crossings:
        # Require a low state after the preceding accepted edge. This preserves
        # the Schmitt-state behavior without a Python loop over every sample.
        lo_pos = np.searchsorted(low_samples, j, side="left") - 1
        if lo_pos < 0:
            continue
        last_low = int(low_samples[lo_pos])
        if last_low <= last_accepted_index:
            continue

        mid_pos = np.searchsorted(midpoint_crossings, j, side="left") - 1
        if mid_pos < 0:
            continue
        k = int(midpoint_crossings[mid_pos])
        if k < last_low or k + 1 >= t.size:
            continue
        dy = smooth[k + 1] - smooth[k]
        tc = t[k] if dy == 0 else t[k] + (threshold - smooth[k]) * (t[k + 1] - t[k]) / dy
        if not crossing or tc - crossing[-1] >= minimum_separation:
            crossing.append(float(tc))
            last_accepted_index = j

    edges = np.asarray(crossing, dtype=float)
    result: dict[str, Any] = {
        "trigger_low_v": float(low),
        "trigger_high_v": float(high),
        "trigger_pp_v": float(high - low),
        "trigger_threshold_v": float(threshold),
        "trigger_edge_count": int(edges.size),
    }
    if edges.size >= 3:
        n = np.arange(edges.size, dtype=float)
        slope, intercept = np.polyfit(n, edges, 1)
        residual = edges - (intercept + slope * n)
        periods = np.diff(edges)
        result.update(
            {
                "trigger_rep_rate_hz": float(1.0 / slope),
                "trigger_period_mean_us": float(np.mean(periods) * 1e6),
                "trigger_period_std_ps": float(np.std(periods, ddof=1) * 1e12),
                "trigger_edge_residual_rms_ps": float(np.sqrt(np.mean(residual**2)) * 1e12),
            }
        )
    else:
        result.update(
            {
                "trigger_rep_rate_hz": np.nan,
                "trigger_period_mean_us": np.nan,
                "trigger_period_std_ps": np.nan,
                "trigger_edge_residual_rms_ps": np.nan,
            }
        )

    if optical_v is not None and edges.size >= 2:
        optical = np.asarray(optical_v, dtype=float)[valid]
        centered = optical - np.median(optical)
        polarity = 1.0 if np.percentile(centered, 99.9) >= -np.percentile(centered, 0.1) else -1.0
        signal = polarity * centered
        window_s = min(0.45 / expected_rep_rate_hz, 1.5e-6)
        delays: list[float] = []
        amplitudes: list[float] = []
        for edge in edges:
            lo_idx = int(np.searchsorted(t, edge, side="left"))
            hi_idx = int(np.searchsorted(t, edge + window_s, side="right"))
            if hi_idx - lo_idx < 5:
                continue
            j = lo_idx + int(np.argmax(signal[lo_idx:hi_idx]))
            delays.append(float(t[j] - edge))
            amplitudes.append(float(signal[j]))
        if delays:
            delay = np.asarray(delays)
            result.update(
                {
                    "trigger_to_electrical_peak_median_ns": float(np.median(delay) * 1e9),
                    "trigger_to_electrical_peak_std_ns": float(np.std(delay, ddof=1) * 1e9) if delay.size > 1 else np.nan,
                    "trigger_to_electrical_peak_n": int(delay.size),
                    "trigger_delay_interpretation": "repeatability only; includes laser, detector, cable and channel delay because no deskew calibration was performed",
                }
            )
    return result

def analyze_scope_record(
    path: str | Path,
    config: TemporalConfig = TemporalConfig(),
    keep_waveform: bool = False,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Analyze one scope record with the correct legacy PDA36A constants."""
    result, waveform = _core.analyze_scope_record(path, config=config, keep_waveform=True)
    setup = parse_scope_setup(path)
    if waveform is None:
        raise RuntimeError("Internal waveform retention failed")
    edge = trigger_edge_metrics(
        waveform["time_s"], waveform["ch3_v"], waveform["ch1_v"],
        expected_rep_rate_hz=(result["nominal_rep_rate_hz"] if np.isfinite(result["nominal_rep_rate_hz"]) else config.expected_rep_rate_hz),
    )

    # Remove the v1 Gaussian pulse-width heuristic. A -3 dB bandwidth can support
    # a step-rise-time scale, but it does not uniquely determine this detector's
    # impulse-response FWHM, especially near 1064 nm.
    result.pop("gaussian_impulse_fwhm_s_est", None)
    result.pop("gaussian_impulse_fwhm_ns_est", None)
    result.update(setup)
    result.update(edge)
    result["ch1_coupling_user_reported"] = config.user_reported_ch1_coupling
    result["photodetector_model"] = "Thorlabs PDA36A (legacy, not PDA36A2)"
    result["photodetector_gain_setting_db"] = 0
    result["detector_bandwidth_note"] = "10 MHz nominal at 0 dB; legacy manual warns NIR photodiode rise time may further reduce effective bandwidth"
    result["scope_bandwidth_note"] = "DSO6104A nominal 1 GHz; scope bandwidth is not the dominant stated limit"

    mode = result["mode"]
    if mode in {"picosecond", "1_ns"}:
        resolution_class = "unresolved_detector_limited"
        interpretation = "electrical detector-chain response; do not report as optical pulse duration"
        unresolved = True
    elif mode == "65_ns":
        resolution_class = "partially_resolved_detector_convolved"
        interpretation = "measured electrical FWHM is detector-convolved and is not a calibrated 65 ns optical width"
        unresolved = False
    else:
        resolution_class = "unknown"
        interpretation = "measurement-chain waveform"
        unresolved = True
    result["pulse_resolution_class"] = resolution_class
    result["measured_fwhm_interpretation"] = interpretation
    result["pulse_unresolved_by_detector"] = unresolved
    result["absolute_trigger_delay_reportable"] = False
    result["absolute_optical_pulse_width_reportable"] = False
    result["intrinsic_laser_timing_jitter_reportable"] = False
    delay_std = result.get("trigger_to_electrical_peak_std_ns", np.nan)
    rise_scale = result.get("instrument_rise_time_ns_est", np.nan)
    result["trigger_delay_repeatability_usable"] = bool(
        np.isfinite(delay_std) and np.isfinite(rise_scale) and delay_std < 0.25 * rise_scale
    )
    result["timing_metric_scope"] = "function generator + scope sampling + laser + detector-chain timing; not intrinsic GSL jitter"
    result["electrical_waveform_reportable"] = True
    result["chopper_role"] = "987 Hz mechanical chopper present in some filenames but user reports it did not chop the beam"

    meta_imp = str(setup.get("ch1_imp", ""))
    expected_imp = "50 Ohm" if result["termination"] == "50ohm" else "1M Ohm"
    result["scope_input_termination_consistent"] = bool(expected_imp.lower().replace(" ", "") in meta_imp.lower().replace(" ", ""))
    result["detector_output_saturation_concern"] = bool(
        result["amplitude_v"] >= (4.5 if result["termination"] == "50ohm" else 9.0)
    )
    # Prefer directly measured CH3 repetition rate when it is available.
    if np.isfinite(result.get("trigger_rep_rate_hz", np.nan)):
        result["rep_rate_measured_hz"] = result["trigger_rep_rate_hz"]

    if not keep_waveform:
        waveform = None
    return result, waveform


def analyze_scope_directory(
    directory: str | Path,
    config: TemporalConfig = TemporalConfig(),
    keep_waveforms_for: Iterable[str] | None = None,
) -> tuple[pd.DataFrame, dict[str, dict[str, Any]]]:
    directory = Path(directory)
    files = sorted(directory.rglob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No oscilloscope CSV files found in {directory}")
    keep = set(keep_waveforms_for or [])
    rows: list[dict[str, Any]] = []
    waveforms: dict[str, dict[str, Any]] = {}
    for path in files:
        result, waveform = analyze_scope_record(path, config=config, keep_waveform=path.name in keep)
        rows.append(result)
        if waveform is not None:
            waveforms[path.name] = waveform
    frame = pd.DataFrame(rows).sort_values(["termination", "mode", "record_duration_s"]).reset_index(drop=True)
    return frame, waveforms


def pulse_energy_from_average_power(average_power_w: float, repetition_rate_hz: float) -> float:
    """Return pulse energy from an independently measured average power."""
    if average_power_w < 0 or repetition_rate_hz <= 0:
        raise ValueError("average_power_w must be nonnegative and repetition_rate_hz positive")
    return average_power_w / repetition_rate_hz


# -----------------------------------------------------------------------------
# Spectral extensions
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class SpectralConfig:
    baseline_method: str = "edge_linear"
    baseline_edge_fraction: float = 0.12
    baseline_polynomial_order: int = 1
    baseline_lambda: float = 1e6
    baseline_asymmetry: float = 0.005
    baseline_iterations: int = 12
    peak_prominence_fraction: float = 0.04
    saturation_dn: float = 65_535.0
    near_saturation_fraction: float = 0.85
    minimum_points_for_voigt: int = 7


def _parse_horiba_datetime(text: Any) -> pd.Timestamp:
    try:
        return pd.to_datetime(str(text), format="%d.%m.%Y %H:%M:%S", errors="coerce")
    except Exception:
        return pd.NaT


def analyze_spectrum(
    path: str | Path,
    config: SpectralConfig = SpectralConfig(),
    keep_arrays: bool = False,
) -> tuple[dict[str, Any], dict[str, np.ndarray] | None]:
    result, arrays = _core.analyze_spectrum(path, config=config, keep_arrays=keep_arrays)
    raw_peak = float(result["peak_raw_counts"])
    near_saturation = bool(raw_peak >= config.near_saturation_fraction * config.saturation_dn)
    family = int(result["wavelength_family_nm"]) if np.isfinite(result["wavelength_family_nm"]) else None
    mode = result["mode"]
    points = float(result["points_per_fwhm"])
    narrow = bool(family == 532 or mode in {"cw", "65_ns"})
    instrument_limited = bool(narrow)
    broad_resolved = bool(mode == "picosecond" and np.isfinite(points) and points >= 10 and not near_saturation)

    if family == 532:
        calibration_status = "visible Hg calibration applies in the same spectral region, but residual/error data were not recorded"
        absolute_status = "moderate confidence; do not infer uncertainty from four displayed digits"
    elif family == 1064:
        calibration_status = "visible Hg calibration was extrapolated into the near-IR; no near-1064 reference line was measured"
        absolute_status = "absolute wavelength is unvalidated at 1064 nm; relative shifts within a fixed setting are more credible"
    else:
        calibration_status = "unknown"
        absolute_status = "unknown"

    peak_nm = float(result["peak_wavelength_nm"])
    centroid_nm = float(result["centroid_wavelength_nm"])
    in_gsl_peak = bool(1054.0 <= peak_nm <= 1074.0) if family == 1064 else False
    in_gsl_centroid = bool(1054.0 <= centroid_nm <= 1074.0) if family == 1064 and np.isfinite(centroid_nm) else False

    if near_saturation:
        width_status = "not_reportable_high_count_linearity_not_established"
    elif narrow and (not np.isfinite(points) or points < 3.0):
        width_status = "undersampled_and_instrument_response_unremoved"
    elif instrument_limited:
        width_status = "instrument_response_unremoved_not_intrinsic_linewidth"
    elif broad_resolved:
        width_status = "resolved_broadband_spectral_fwhm"
    else:
        width_status = "empirical_fwhm_with_configuration_caveat"

    result.update(
        {
            "near_saturation": near_saturation,
            "near_saturation_threshold_counts": config.near_saturation_fraction * config.saturation_dn,
            "detector_full_scale_fraction": raw_peak / config.saturation_dn,
            "detector_linearity_warning": "at or above the conservative 85% full-scale QC threshold; detector linearity was not independently established" if near_saturation else "below the conservative 85% full-scale QC threshold",
            "instrument_limited_likely": instrument_limited,
            "broadband_width_resolved": broad_resolved,
            "width_status": width_status,
            "intrinsic_narrow_linewidth_reportable": False if narrow else np.nan,
            "absolute_wavelength_calibration_status": calibration_status,
            "absolute_center_claim_status": absolute_status,
            "relative_wavelength_comparisons_preferred": bool(family == 1064),
            "gsl_nominal_band_low_nm": 1054.0 if family == 1064 else np.nan,
            "gsl_nominal_band_high_nm": 1074.0 if family == 1064 else np.nan,
            "peak_within_gsl_1064_plusminus10_nm": in_gsl_peak if family == 1064 else np.nan,
            "centroid_within_gsl_1064_plusminus10_nm": in_gsl_centroid if family == 1064 else np.nan,
            "spectrometer_focal_length_mm": 550.0,
            "spectrometer_entrance_aperture_f_number": 6.4,
            "slit_zero_interpretation": "minimum/fully closed setting; actual physical gap unknown" if result["slit_um_filename"] == 0 else "filename nominal slit setting",
            "spectral_response_corrected": False,
            "spectral_response_note": "broadband peak/centroid can be weighted by grating and detector response",
            "acquired_datetime": _parse_horiba_datetime(result.get("acquired_metadata", "")),
        }
    )
    return result, arrays


def analyze_spectral_directory(
    directory: str | Path,
    config: SpectralConfig = SpectralConfig(),
    keep_arrays_for: Iterable[str] | None = None,
) -> tuple[pd.DataFrame, dict[str, dict[str, np.ndarray]]]:
    directory = Path(directory)
    files = sorted(directory.rglob("*.txt"))
    if not files:
        raise FileNotFoundError(f"No spectrum text files found in {directory}")
    keep = set(keep_arrays_for or [])
    rows: list[dict[str, Any]] = []
    arrays: dict[str, dict[str, np.ndarray]] = {}
    for path in files:
        result, data = analyze_spectrum(path, config=config, keep_arrays=path.name in keep)
        rows.append(result)
        if data is not None:
            arrays[path.name] = data
    frame = pd.DataFrame(rows).sort_values(
        ["wavelength_family_nm", "mode", "grating_g_mm", "slit_um_filename", "acquired_datetime", "file"]
    ).reset_index(drop=True)
    return assign_spectral_sequences(frame), arrays


def assign_spectral_sequences(metrics: pd.DataFrame, gap_seconds: float = 120.0) -> pd.DataFrame:
    """Assign acquisition sequences without pretending all same-name files are replicates."""
    frame = metrics.copy()
    physical = [
        "wavelength_family_nm", "mode", "grating_g_mm", "slit_um_filename",
        "acquisition_s_filename", "spectrometer_center_nm_metadata",
    ]
    frame["sequence_id"] = ""
    for _, indices in frame.groupby(physical, dropna=False).groups.items():
        idx = list(indices)
        group = frame.loc[idx].sort_values(["acquired_datetime", "file"])
        sequence = 1
        previous_time: pd.Timestamp | None = None
        previous_sat: bool | None = None
        labels: dict[int, str] = {}
        for row_index, row in group.iterrows():
            current_time = row["acquired_datetime"]
            current_sat = bool(row["near_saturation"])
            if previous_time is not None and pd.notna(current_time) and pd.notna(previous_time):
                if (current_time - previous_time).total_seconds() > gap_seconds or current_sat != previous_sat:
                    sequence += 1
            labels[row_index] = f"seq{sequence}"
            previous_time = current_time if pd.notna(current_time) else previous_time
            previous_sat = current_sat
        for row_index, label in labels.items():
            frame.at[row_index, "sequence_id"] = label
    return frame


def summarize_replicate_spectra(
    metrics: pd.DataFrame,
    exclude_near_saturation: bool = True,
) -> pd.DataFrame:
    """Summarize physically matched, time-contiguous acquisition sequences."""
    frame = assign_spectral_sequences(metrics) if "sequence_id" not in metrics else metrics.copy()
    if exclude_near_saturation:
        frame = frame[~frame["near_saturation"]].copy()
    keys = [
        "wavelength_family_nm", "mode", "grating_g_mm", "slit_um_filename",
        "acquisition_s_filename", "sequence_id",
    ]
    rows: list[dict[str, Any]] = []
    for group_key, group in frame.groupby(keys, dropna=False):
        rows.append(
            {
                **dict(zip(keys, group_key)),
                "n": len(group),
                "acquired_first": group["acquired_datetime"].min(),
                "acquired_last": group["acquired_datetime"].max(),
                "peak_wavelength_mean_nm": float(group["peak_wavelength_nm"].mean()),
                "peak_wavelength_std_nm": float(group["peak_wavelength_nm"].std(ddof=1)) if len(group) > 1 else np.nan,
                "centroid_wavelength_mean_nm": float(group["centroid_wavelength_nm"].mean()),
                "fwhm_mean_nm": float(group["fwhm_nm"].mean()),
                "fwhm_std_nm": float(group["fwhm_nm"].std(ddof=1)) if len(group) > 1 else np.nan,
                "peak_counts_cv_percent": float(100 * group["peak_corrected_counts"].std(ddof=1) / group["peak_corrected_counts"].mean()) if len(group) > 1 and group["peak_corrected_counts"].mean() != 0 else np.nan,
                "repeatability_scope": "within this contiguous acquisition sequence only; excludes wavelength calibration uncertainty",
            }
        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(keys).reset_index(drop=True)


# Rebuild the public export list: core exports plus v2 additions/overrides.
__all__ = sorted(
    set(getattr(_core, "__all__", []))
    | {
        "CLAIM_LEVELS",
        "CameraConfig",
        "TemporalConfig",
        "SpectralConfig",
        "analyze_spatial_frame",
        "analyze_spatial_directory",
        "aggregate_spatial_metrics",
        "parse_spatial_png_filename",
        "analyze_spatial_png_frame",
        "analyze_spatial_png_directory",
        "fit_beam_propagation",
        "parse_scope_setup",
        "trigger_edge_metrics",
        "analyze_scope_record",
        "analyze_scope_directory",
        "pulse_energy_from_average_power",
        "analyze_spectrum",
        "analyze_spectral_directory",
        "assign_spectral_sequences",
        "summarize_replicate_spectra",
    }
)
