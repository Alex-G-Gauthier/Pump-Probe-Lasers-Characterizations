"""Reusable analysis tools for the 532/1064 nm laser-characterization notebooks.

The module is intentionally conservative: it carries quality-control flags forward
instead of silently reporting physically invalid M^2, pulse-width, or linewidth
values.  Functions are written for the native files in the supplied archive:

* Basler acA3088-57uc BayerRG12 full-frame RAW files (3088 x 2064, uint16)
* Keysight/Agilent DSO6104A CSV exports with metadata and analog channels
* HORIBA text spectra using decimal commas and tab-separated numeric columns

Author: OpenAI, prepared for the user's experimental dataset.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import csv
import math
import re
import warnings

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage, sparse
from scipy.optimize import curve_fit, least_squares
from scipy.signal import find_peaks
from scipy.sparse.linalg import spsolve
from scipy.special import voigt_profile

C_LIGHT = 299_792_458.0  # m/s, exact by SI definition


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

def robust_sigma(values: np.ndarray | Sequence[float]) -> float:
    """Return a median-absolute-deviation estimate of one standard deviation."""
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    med = np.median(x)
    sigma = 1.4826 * np.median(np.abs(x - med))
    if sigma == 0 and x.size > 1:
        # Quantized data can have MAD=0. Use a clipped standard deviation fallback.
        q1, q99 = np.percentile(x, [1, 99])
        core = x[(x >= q1) & (x <= q99)]
        sigma = float(np.std(core)) if core.size > 1 else 0.0
    return float(sigma)


def ensure_directory(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def find_data_root(start: str | Path | None = None) -> Path:
    """Locate the extracted ``532 and 1064`` data directory.

    The search is deliberately explicit so users can run the notebooks from the
    project directory, beside the data directory, or in the original sandbox.
    """
    starts: list[Path] = []
    if start is not None:
        starts.append(Path(start).expanduser())
    starts.extend(
        [
            Path.cwd(),
            Path.cwd().parent,
            Path("/mnt/data/laser_project"),
            Path("/mnt/data"),
        ]
    )
    checked: list[Path] = []
    for root in starts:
        for candidate in (
            root,
            root / "532 and 1064",
            root / "data" / "532 and 1064",
            root / "laser_project" / "532 and 1064",
        ):
            candidate = candidate.resolve()
            checked.append(candidate)
            if (
                candidate.is_dir()
                and (candidate / "Pulse_Data").is_dir()
                and (candidate / "Spectral_Data").is_dir()
                and (candidate / "532nm_Spatial").is_dir()
            ):
                return candidate
    checked_text = "\n".join(f"  - {p}" for p in checked)
    raise FileNotFoundError(
        "Could not locate the extracted '532 and 1064' folder. Edit DATA_ROOT in "
        f"the notebook. Checked:\n{checked_text}"
    )


def minmax_envelope(x: np.ndarray, y: np.ndarray, target_bins: int = 2500) -> pd.DataFrame:
    """Downsample a long waveform without hiding narrow extrema."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size != y.size:
        raise ValueError("x and y must have the same length")
    if x.size <= target_bins:
        return pd.DataFrame({"x": x, "y_min": y, "y_max": y, "y_median": y})
    block = max(1, x.size // target_bins)
    n = x.size // block
    xr = x[: n * block].reshape(n, block)
    yr = y[: n * block].reshape(n, block)
    return pd.DataFrame(
        {
            "x": xr[:, 0],
            "y_min": np.nanmin(yr, axis=1),
            "y_max": np.nanmax(yr, axis=1),
            "y_median": np.nanmedian(yr, axis=1),
        }
    )


# -----------------------------------------------------------------------------
# Spatial analysis (Basler BayerRG12)
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class CameraConfig:
    width_px: int = 3088
    height_px: int = 2064
    sensor_pixel_pitch_um: float = 2.4
    bit_depth: int = 12
    bayer_pattern: str = "RGGB"
    magnification: float = 1.0
    wavelength_nm: float = 532.0
    raw_dtype: str = "<u2"
    border_px_green: int = 40
    core_fraction: float = 0.20
    roi_expand: float = 4.0
    min_roi_halfwidth_green_px: int = 45
    max_roi_halfwidth_green_px: int = 420
    primary_aperture_scale: float = 2.0
    saturation_fraction_threshold: float = 1e-6
    secondary_peak_ratio_threshold: float = 0.10

    @property
    def full_scale_dn(self) -> int:
        return 2**self.bit_depth - 1

    @property
    def green_pixel_pitch_um(self) -> float:
        # Averaging the two green samples in each 2x2 Bayer cell creates one
        # green superpixel per cell.
        return 2.0 * self.sensor_pixel_pitch_um / self.magnification


_SPATIAL_POSITION_RE = re.compile(r"Position_(?P<position>[-+]?\d+(?:\.\d+)?)cm", re.I)
_SPATIAL_REP_RE = re.compile(r"_(?P<rep>\d+)\.raw$", re.I)


def parse_spatial_filename(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    name = path.name
    pm = _SPATIAL_POSITION_RE.search(name)
    rm = _SPATIAL_REP_RE.search(name)
    if not pm:
        raise ValueError(f"Cannot parse propagation position from {name}")
    focal = "yes" if "Focal_Yes" in name else "no"
    exposure_match = re.search(r"Exposure_(\d+(?:\.\d+)?)us", name, re.I)
    nd_match = re.search(r"ND_([\d_]+)", name, re.I)
    return {
        "file": name,
        "path": str(path),
        "position_cm": float(pm.group("position")),
        "position_m": float(pm.group("position")) / 100.0,
        "focal": focal,
        "replicate": int(rm.group("rep")) if rm else 0,
        "exposure_us": float(exposure_match.group(1)) if exposure_match else np.nan,
        "nd_label": nd_match.group(1).replace("_", ".") if nd_match else "",
    }


def load_bayer_raw(path: str | Path, config: CameraConfig = CameraConfig()) -> np.ndarray:
    """Load a full-frame unpacked BayerRG12 image stored as little-endian uint16."""
    path = Path(path)
    data = np.fromfile(path, dtype=np.dtype(config.raw_dtype))
    expected = config.width_px * config.height_px
    if data.size != expected:
        raise ValueError(
            f"{path.name}: expected {expected:,} pixels for "
            f"{config.width_px}x{config.height_px}, found {data.size:,}."
        )
    raw = data.reshape(config.height_px, config.width_px)
    # Some camera APIs store 12-bit values left-aligned in 16 bits. Detect that
    # case conservatively. The supplied files are right-aligned (0...4095).
    if raw.max(initial=0) > config.full_scale_dn:
        sampled = raw.ravel()[:: max(1, raw.size // 100_000)]
        if np.mean((sampled & 0xF) == 0) > 0.95:
            raw = raw >> 4
        else:
            raise ValueError(
                f"{path.name}: values exceed {config.full_scale_dn}; unknown 12-bit packing."
            )
    return raw


def green_superpixel(raw: np.ndarray, pattern: str = "RGGB") -> np.ndarray:
    """Return an interpolation-free green image from a Bayer mosaic.

    For BayerRG / RGGB, the two green sites in each 2x2 cell are averaged.  This
    avoids color-demosaicing artifacts in second-moment beam measurements.
    """
    if raw.ndim != 2:
        raise ValueError("raw Bayer data must be a 2D array")
    pattern = pattern.upper()
    if pattern != "RGGB":
        raise NotImplementedError("Only RGGB/BayerRG is implemented for this dataset")
    g1 = raw[0::2, 1::2].astype(np.float64)
    g2 = raw[1::2, 0::2].astype(np.float64)
    return 0.5 * (g1 + g2)


def load_png_green(path: str | Path, background_path: str | Path | None = None) -> np.ndarray:
    """Load the green channel of a saved RGB PNG and optionally subtract a background."""
    arr = np.asarray(Image.open(path)).astype(np.float64)
    if arr.ndim == 2:
        green = arr
    elif arr.ndim == 3 and arr.shape[2] >= 2:
        green = arr[..., 1]
    else:
        raise ValueError(f"Unsupported PNG array shape {arr.shape}")
    if background_path is not None:
        bg = np.asarray(Image.open(background_path)).astype(np.float64)
        bg_green = bg if bg.ndim == 2 else bg[..., 1]
        if bg_green.shape != green.shape:
            raise ValueError("background and image shapes do not match")
        green = green - bg_green
    return np.clip(green, 0, None)


def estimate_background_plane(
    image: np.ndarray,
    border_px: int = 40,
    stride: int = 4,
    iterations: int = 5,
) -> tuple[np.ndarray, float, dict[str, float]]:
    """Fit a robust planar background to perimeter pixels.

    Positive residuals are iteratively rejected so stray spots near an edge do
    not set the background. Returns the full plane, residual noise sigma, and
    plane coefficients in DN per normalized coordinate.
    """
    image = np.asarray(image, dtype=float)
    h, w = image.shape
    b = int(np.clip(border_px, 2, max(2, min(h, w) // 4)))
    border = np.zeros_like(image, dtype=bool)
    border[:b, :] = True
    border[-b:, :] = True
    border[:, :b] = True
    border[:, -b:] = True
    yy, xx = np.nonzero(border)
    select = (np.arange(xx.size) % max(1, stride)) == 0
    xx = xx[select]
    yy = yy[select]
    zz = image[yy, xx]
    xn = (xx - (w - 1) / 2) / max(w - 1, 1)
    yn = (yy - (h - 1) / 2) / max(h - 1, 1)
    design = np.column_stack([np.ones_like(xn), xn, yn])
    keep = np.isfinite(zz)
    beta = np.array([float(np.nanmedian(zz)), 0.0, 0.0])
    sigma = robust_sigma(zz)
    for _ in range(iterations):
        if keep.sum() < 10:
            break
        beta, *_ = np.linalg.lstsq(design[keep], zz[keep], rcond=None)
        resid = zz - design @ beta
        sigma = robust_sigma(resid[keep])
        if not np.isfinite(sigma) or sigma <= 0:
            break
        # Reject bright contamination more aggressively than negative outliers.
        keep = np.isfinite(resid) & (resid > -4.0 * sigma) & (resid < 2.5 * sigma)
    ygrid, xgrid = np.indices(image.shape, dtype=float)
    xgrid_n = (xgrid - (w - 1) / 2) / max(w - 1, 1)
    ygrid_n = (ygrid - (h - 1) / 2) / max(h - 1, 1)
    plane = beta[0] + beta[1] * xgrid_n + beta[2] * ygrid_n
    return plane, float(max(sigma, 0.0)), {"b0": beta[0], "bx": beta[1], "by": beta[2]}


def _component_containing(mask: np.ndarray, y: int, x: int) -> np.ndarray:
    labels, count = ndimage.label(mask)
    if count == 0:
        return np.zeros_like(mask, dtype=bool)
    label_id = labels[int(np.clip(y, 0, mask.shape[0] - 1)), int(np.clip(x, 0, mask.shape[1] - 1))]
    if label_id == 0:
        # Use the closest component centroid if the exact seed is just below threshold.
        objects = ndimage.find_objects(labels)
        candidates: list[tuple[float, int]] = []
        for idx, slc in enumerate(objects, start=1):
            if slc is None:
                continue
            cy = 0.5 * (slc[0].start + slc[0].stop - 1)
            cx = 0.5 * (slc[1].start + slc[1].stop - 1)
            candidates.append(((cy - y) ** 2 + (cx - x) ** 2, idx))
        if not candidates:
            return np.zeros_like(mask, dtype=bool)
        label_id = min(candidates)[1]
    return labels == label_id


def locate_primary_roi(
    corrected: np.ndarray,
    noise_sigma: float,
    config: CameraConfig,
    center_hint: tuple[float, float] | None = None,
    hint_search_radius_px: int = 100,
) -> tuple[tuple[slice, slice], tuple[float, float], dict[str, float]]:
    """Locate the dominant lobe and create a compact fit ROI.

    ``center_hint`` is (x, y) in green-superpixel coordinates and can be used to
    keep the same physical lobe selected when a secondary reflection becomes
    brighter at another z position.
    """
    smooth = ndimage.gaussian_filter(np.asarray(corrected, dtype=float), sigma=1.5)
    h, w = smooth.shape
    if center_hint is None:
        seed_y, seed_x = np.unravel_index(np.nanargmax(smooth), smooth.shape)
    else:
        hx, hy = center_hint
        hx = int(np.clip(round(hx), 0, w - 1))
        hy = int(np.clip(round(hy), 0, h - 1))
        r = int(max(5, hint_search_radius_px))
        y0, y1 = max(0, hy - r), min(h, hy + r + 1)
        x0, x1 = max(0, hx - r), min(w, hx + r + 1)
        local_y, local_x = np.unravel_index(np.nanargmax(smooth[y0:y1, x0:x1]), (y1 - y0, x1 - x0))
        seed_y, seed_x = y0 + local_y, x0 + local_x
    peak = float(smooth[seed_y, seed_x])
    threshold = max(config.core_fraction * peak, 8.0 * max(noise_sigma, 0.0))
    component = _component_containing(smooth >= threshold, seed_y, seed_x)
    if not component.any():
        threshold = max(0.08 * peak, 4.0 * max(noise_sigma, 0.0))
        component = _component_containing(smooth >= threshold, seed_y, seed_x)
    ys, xs = np.where(component)
    if xs.size == 0:
        core_w = core_h = 8
        core_bounds = (seed_x - 4, seed_x + 5, seed_y - 4, seed_y + 5)
    else:
        xlo, xhi = int(xs.min()), int(xs.max()) + 1
        ylo, yhi = int(ys.min()), int(ys.max()) + 1
        core_w, core_h = xhi - xlo, yhi - ylo
        core_bounds = (xlo, xhi, ylo, yhi)
    core_size = max(core_w, core_h, 8)
    half = int(np.clip(config.roi_expand * core_size, config.min_roi_halfwidth_green_px, config.max_roi_halfwidth_green_px))
    x0, x1 = max(0, seed_x - half), min(w, seed_x + half + 1)
    y0, y1 = max(0, seed_y - half), min(h, seed_y + half + 1)
    return (slice(y0, y1), slice(x0, x1)), (float(seed_x), float(seed_y)), {
        "seed_peak_dn": peak,
        "core_threshold_dn": float(threshold),
        "core_width_px": float(core_w),
        "core_height_px": float(core_h),
        "roi_x0": float(x0),
        "roi_x1": float(x1),
        "roi_y0": float(y0),
        "roi_y1": float(y1),
        "core_x0": float(core_bounds[0]),
        "core_x1": float(core_bounds[1]),
        "core_y0": float(core_bounds[2]),
        "core_y1": float(core_bounds[3]),
    }


def _gaussian_model(params: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    offset, grad_x, grad_y, amplitude, x0, y0, log_wx, log_wy, theta = params
    wx = np.exp(log_wx)
    wy = np.exp(log_wy)
    ct, st = np.cos(theta), np.sin(theta)
    dx, dy = x - x0, y - y0
    xp = ct * dx + st * dy
    yp = -st * dx + ct * dy
    return offset + grad_x * dx + grad_y * dy + amplitude * np.exp(-2.0 * ((xp / wx) ** 2 + (yp / wy) ** 2))


def fit_rotated_gaussian(
    roi: np.ndarray,
    noise_sigma: float = 0.0,
    max_fit_points: int = 25_000,
    center_hint_roi: tuple[float, float] | None = None,
) -> dict[str, float]:
    """Robustly fit an elliptical 2D Gaussian plus a planar local background."""
    data = np.asarray(roi, dtype=float)
    h, w = data.shape
    yy, xx = np.indices(data.shape, dtype=float)
    peak_y, peak_x = np.unravel_index(np.nanargmax(data), data.shape)
    if center_hint_roi is not None:
        hint_x = float(np.clip(center_hint_roi[0], 0, w - 1))
        hint_y = float(np.clip(center_hint_roi[1], 0, h - 1))
        search = max(5, int(0.18 * min(h, w)))
        ys0, ys1 = max(0, int(round(hint_y)) - search), min(h, int(round(hint_y)) + search + 1)
        xs0, xs1 = max(0, int(round(hint_x)) - search), min(w, int(round(hint_x)) + search + 1)
        py, px = np.unravel_index(np.nanargmax(data[ys0:ys1, xs0:xs1]), (ys1 - ys0, xs1 - xs0))
        peak_y, peak_x = ys0 + py, xs0 + px
    else:
        hint_x, hint_y = float(peak_x), float(peak_y)
    background0 = float(np.percentile(data, 10))
    positive = np.clip(data - background0, 0, None)
    threshold = max(0.10 * positive.max(initial=0.0), 4.0 * max(noise_sigma, 0.0))
    local_radius = max(12.0, 0.40 * min(h, w))
    local_region = (xx - hint_x) ** 2 + (yy - hint_y) ** 2 <= local_radius**2
    mask = (positive >= threshold) & local_region
    component = _component_containing(mask, peak_y, peak_x)
    if component.sum() >= 12:
        mask = component
    elif mask.sum() < 12:
        mask = (positive >= max(0.03 * positive.max(initial=0.0), max(noise_sigma, 0.0))) & local_region
    weights = positive * mask
    power = weights.sum()
    if power > 0:
        x0 = float((weights * xx).sum() / power)
        y0 = float((weights * yy).sum() / power)
        var_x = float((weights * (xx - x0) ** 2).sum() / power)
        var_y = float((weights * (yy - y0) ** 2).sum() / power)
        wx0 = max(2.0, 2.0 * math.sqrt(max(var_x, 1.0)))
        wy0 = max(2.0, 2.0 * math.sqrt(max(var_y, 1.0)))
    else:
        x0, y0 = float(peak_x), float(peak_y)
        wx0 = wy0 = max(3.0, min(h, w) / 10.0)
    amplitude0 = max(float(data[peak_y, peak_x] - background0), 1.0)
    p0 = np.array([background0, 0.0, 0.0, amplitude0, x0, y0, np.log(wx0), np.log(wy0), 0.0])
    max_width = float(max(h, w) * 1.5)
    gradient_bound = max(amplitude0, 1.0) / max(min(h, w), 1)
    center_tolerance = max(8.0, 0.22 * min(h, w))
    center_x_low = max(0.0, hint_x - center_tolerance)
    center_x_high = min(w - 1.0, hint_x + center_tolerance)
    center_y_low = max(0.0, hint_y - center_tolerance)
    center_y_high = min(h - 1.0, hint_y + center_tolerance)
    x0 = float(np.clip(x0, center_x_low, center_x_high))
    y0 = float(np.clip(y0, center_y_low, center_y_high))
    p0[4] = x0
    p0[5] = y0
    lower = np.array(
        [
            -0.25 * amplitude0,
            -gradient_bound,
            -gradient_bound,
            0.0,
            center_x_low,
            center_y_low,
            np.log(1.0),
            np.log(1.0),
            -np.pi / 2,
        ]
    )
    upper = np.array(
        [
            background0 + 0.5 * amplitude0,
            gradient_bound,
            gradient_bound,
            3.0 * amplitude0 + max(background0, 0.0),
            center_x_high,
            center_y_high,
            np.log(max_width),
            np.log(max_width),
            np.pi / 2,
        ]
    )
    signal_mask = (positive > max(0.008 * positive.max(initial=0.0), 2.0 * max(noise_sigma, 0.0))) & local_region
    background_sample = np.zeros_like(signal_mask)
    stride = max(1, int(math.sqrt(data.size / max(max_fit_points // 3, 1))))
    background_sample[::stride, ::stride] = True
    fit_mask = signal_mask | background_sample
    xfit, yfit, zfit = xx[fit_mask], yy[fit_mask], data[fit_mask]
    if zfit.size > max_fit_points:
        idx = np.linspace(0, zfit.size - 1, max_fit_points, dtype=int)
        xfit, yfit, zfit = xfit[idx], yfit[idx], zfit[idx]
    scale = max(robust_sigma(zfit[zfit <= np.percentile(zfit, 60)]), 0.005 * amplitude0, 1e-6)

    def residual(p: np.ndarray) -> np.ndarray:
        return (_gaussian_model(p, xfit, yfit) - zfit) / scale

    result = least_squares(
        residual,
        p0,
        bounds=(lower, upper),
        loss="soft_l1",
        f_scale=1.0,
        x_scale="jac",
        max_nfev=250,
    )
    p = result.x
    offset, grad_x, grad_y, amplitude, x0, y0, log_wx, log_wy, theta = p
    wx, wy = float(np.exp(log_wx)), float(np.exp(log_wy))
    # Canonicalize so wx >= wy and theta describes the major Gaussian axis.
    if wy > wx:
        wx, wy = wy, wx
        theta += np.pi / 2
    theta = (theta + np.pi / 2) % np.pi - np.pi / 2
    model = _gaussian_model(
        np.array([offset, grad_x, grad_y, amplitude, x0, y0, np.log(wx), np.log(wy), theta]),
        xx,
        yy,
    )
    ct, st = np.cos(theta), np.sin(theta)
    dx, dy = xx - x0, yy - y0
    xp = ct * dx + st * dy
    yp = -st * dx + ct * dy
    core = (xp / max(wx, 1e-9)) ** 2 + (yp / max(wy, 1e-9)) ** 2 <= 1.5**2
    rmse = float(np.sqrt(np.mean((data[core] - model[core]) ** 2))) if core.any() else np.nan
    nrmse = rmse / max(amplitude, 1e-12)
    return {
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "offset_dn": float(offset),
        "gradient_x_dn_per_px": float(grad_x),
        "gradient_y_dn_per_px": float(grad_y),
        "amplitude_dn": float(amplitude),
        "x0_roi_px": float(x0),
        "y0_roi_px": float(y0),
        "wx_1e2_px": wx,
        "wy_1e2_px": wy,
        "theta_rad": float(theta),
        "theta_deg": float(np.degrees(theta)),
        "fit_rmse_dn": rmse,
        "fit_nrmse": float(nrmse),
        "cost": float(result.cost),
        "nfev": int(result.nfev),
    }


def second_moments(
    intensity: np.ndarray,
    x_coords: np.ndarray | None = None,
    y_coords: np.ndarray | None = None,
    mask: np.ndarray | None = None,
) -> dict[str, float]:
    """Compute ISO-style centroid/covariance and D4-sigma beam diameters."""
    image = np.asarray(intensity, dtype=float)
    if mask is None:
        mask = np.ones_like(image, dtype=bool)
    weights = np.where(mask & np.isfinite(image) & (image > 0), image, 0.0)
    total = float(weights.sum())
    if total <= 0:
        raise ValueError("No positive intensity in the selected aperture")
    yy, xx = np.indices(image.shape, dtype=float)
    if x_coords is not None:
        xx = np.asarray(x_coords, dtype=float)
    if y_coords is not None:
        yy = np.asarray(y_coords, dtype=float)
    xbar = float((weights * xx).sum() / total)
    ybar = float((weights * yy).sum() / total)
    dx, dy = xx - xbar, yy - ybar
    cxx = float((weights * dx * dx).sum() / total)
    cyy = float((weights * dy * dy).sum() / total)
    cxy = float((weights * dx * dy).sum() / total)
    cov = np.array([[cxx, cxy], [cxy, cyy]])
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = np.maximum(eigvals[order], 0.0)
    eigvecs = eigvecs[:, order]
    major_vec = eigvecs[:, 0]
    angle = math.atan2(major_vec[1], major_vec[0])
    return {
        "power_arb": total,
        "centroid_x": xbar,
        "centroid_y": ybar,
        "sigma_x": math.sqrt(max(cxx, 0.0)),
        "sigma_y": math.sqrt(max(cyy, 0.0)),
        "sigma_xy": cxy,
        "d4sigma_x": 4.0 * math.sqrt(max(cxx, 0.0)),
        "d4sigma_y": 4.0 * math.sqrt(max(cyy, 0.0)),
        "d4sigma_major": 4.0 * math.sqrt(eigvals[0]),
        "d4sigma_minor": 4.0 * math.sqrt(eigvals[1]),
        "major_axis_angle_rad": angle,
        "major_axis_angle_deg": math.degrees(angle),
        "ellipticity_major_minor": math.sqrt(eigvals[0] / eigvals[1]) if eigvals[1] > 0 else np.inf,
    }


def analyze_spatial_frame(
    path: str | Path,
    config: CameraConfig = CameraConfig(),
    center_hint: tuple[float, float] | None = None,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Analyze one RAW frame and return metrics plus arrays for diagnostic plotting."""
    path = Path(path)
    raw = load_bayer_raw(path, config)
    green = green_superpixel(raw, config.bayer_pattern)
    plane, noise_sigma, bg_coeff = estimate_background_plane(
        green, border_px=config.border_px_green
    )
    corrected = np.clip(green - plane, 0, None)
    roi_slice, seed, roi_info = locate_primary_roi(
        corrected, noise_sigma, config, center_hint=center_hint
    )
    roi = corrected[roi_slice]
    y0g = roi_slice[0].start or 0
    x0g = roi_slice[1].start or 0
    fit = fit_rotated_gaussian(
        roi,
        noise_sigma=noise_sigma,
        center_hint_roi=(seed[0] - x0g, seed[1] - y0g),
    )
    x0 = x0g + fit["x0_roi_px"]
    y0 = y0g + fit["y0_roi_px"]
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
    noise_floor = max(3.0 * noise_sigma, 0.0005 * corrected.max(initial=0.0))
    intensity_for_moments = np.where(corrected > noise_floor, corrected, 0.0)
    primary_mom = second_moments(intensity_for_moments, mask=primary_mask)
    full_mom = second_moments(intensity_for_moments)
    total_signal = float(intensity_for_moments.sum())
    primary_signal = float(intensity_for_moments[primary_mask].sum())
    primary_fraction = primary_signal / total_signal if total_signal > 0 else np.nan
    outside = ndimage.gaussian_filter(corrected, 1.5).copy()
    outside[primary_mask] = 0.0
    secondary_peak = float(outside.max(initial=0.0))
    primary_peak = float(corrected.max(initial=0.0))
    secondary_ratio = secondary_peak / primary_peak if primary_peak > 0 else np.nan
    # Axis-aligned ellipse extents for clipping tests.
    scale = config.primary_aperture_scale
    x_extent = scale * math.sqrt((wx * ct) ** 2 + (wy * st) ** 2)
    y_extent = scale * math.sqrt((wx * st) ** 2 + (wy * ct) ** 2)
    margin = 3.0
    clipped = bool(
        x0 - x_extent <= margin
        or x0 + x_extent >= corrected.shape[1] - 1 - margin
        or y0 - y_extent <= margin
        or y0 + y_extent >= corrected.shape[0] - 1 - margin
    )
    saturation_dn = 0.995 * config.full_scale_dn
    saturation_fraction = float(np.mean(raw >= saturation_dn))
    saturated = bool(
        saturation_fraction > config.saturation_fraction_threshold
        or raw.max(initial=0) >= saturation_dn
    )
    poor_fit = bool(not fit["success"] or fit["fit_nrmse"] > 0.15)
    secondary_structure = bool(secondary_ratio > config.secondary_peak_ratio_threshold)
    valid_for_m2 = bool(
        (not clipped)
        and (not saturated)
        and (not poor_fit)
        and primary_fraction >= 0.80
        and (not secondary_structure)
    )
    pitch_um = config.green_pixel_pitch_um
    meta = parse_spatial_filename(path)
    metrics: dict[str, Any] = {
        **meta,
        **{f"background_{k}": float(v) for k, v in bg_coeff.items()},
        "background_noise_dn": float(noise_sigma),
        "raw_min_dn": int(raw.min(initial=0)),
        "raw_max_dn": int(raw.max(initial=0)),
        "saturation_fraction": saturation_fraction,
        "saturated": saturated,
        "primary_center_x_green_px": float(x0),
        "primary_center_y_green_px": float(y0),
        "primary_center_x_mm": float(x0 * pitch_um / 1000.0),
        "primary_center_y_mm": float(y0 * pitch_um / 1000.0),
        "gaussian_w_major_um": float(wx * pitch_um),
        "gaussian_w_minor_um": float(wy * pitch_um),
        "gaussian_theta_deg": fit["theta_deg"],
        "gaussian_fit_nrmse": fit["fit_nrmse"],
        "primary_fraction": float(primary_fraction),
        "secondary_peak_ratio": float(secondary_ratio),
        "clipped": clipped,
        "poor_gaussian_fit": poor_fit,
        "secondary_structure": secondary_structure,
        "valid_for_m2": valid_for_m2,
        "d4sigma_x_um": float(primary_mom["d4sigma_x"] * pitch_um),
        "d4sigma_y_um": float(primary_mom["d4sigma_y"] * pitch_um),
        "d4sigma_major_um": float(primary_mom["d4sigma_major"] * pitch_um),
        "d4sigma_minor_um": float(primary_mom["d4sigma_minor"] * pitch_um),
        "primary_ellipticity": float(primary_mom["ellipticity_major_minor"]),
        "primary_angle_deg": float(primary_mom["major_axis_angle_deg"]),
        "full_field_d4sigma_major_um": float(full_mom["d4sigma_major"] * pitch_um),
        "full_field_d4sigma_minor_um": float(full_mom["d4sigma_minor"] * pitch_um),
        **{k: v for k, v in fit.items() if k not in {"message"}},
        **roi_info,
        "seed_x_green_px": seed[0],
        "seed_y_green_px": seed[1],
        "camera_pitch_um": config.sensor_pixel_pitch_um,
        "green_pitch_um": pitch_um,
        "magnification": config.magnification,
    }
    arrays = {
        "raw": raw,
        "green": green,
        "background_plane": plane,
        "corrected": corrected,
        "primary_mask": primary_mask,
        "roi": roi,
    }
    return metrics, arrays


def analyze_spatial_directory(
    directory: str | Path,
    config: CameraConfig = CameraConfig(),
    center_hints: Mapping[str, tuple[float, float]] | None = None,
    return_diagnostics_for: Iterable[str] | None = None,
    group_consistent_lobe: bool = True,
) -> tuple[pd.DataFrame, dict[str, dict[str, np.ndarray]]]:
    """Analyze all RAW frames in a directory.

    By default, a median frame is formed for each (focal condition, z position)
    group and its dominant lobe is used as a common center hint for all
    replicates. This prevents replicate-to-replicate switching between a real
    beam spot and a similarly bright reflection. Manual hints can be supplied
    either by exact filename or by a group key such as ``"no:20"``.
    """
    directory = Path(directory)
    files = sorted(directory.glob("*.raw"))
    if not files:
        raise FileNotFoundError(f"No .raw files found in {directory}")
    center_hints = center_hints or {}
    diagnostic_names = set(return_diagnostics_for or [])
    parsed = {path: parse_spatial_filename(path) for path in files}
    groups: dict[tuple[str, float], list[Path]] = {}
    for path, meta in parsed.items():
        groups.setdefault((meta["focal"], meta["position_cm"]), []).append(path)
    automatic_group_hints: dict[tuple[str, float], tuple[float, float]] = {}
    if group_consistent_lobe:
        for key, group_files in groups.items():
            group_key = f"{key[0]}:{key[1]:g}"
            if group_key in center_hints:
                automatic_group_hints[key] = center_hints[group_key]
                continue
            # Exact filename hints take precedence; use their median as a group hint.
            explicit = [center_hints[p.name] for p in group_files if p.name in center_hints]
            if explicit:
                automatic_group_hints[key] = tuple(np.median(np.asarray(explicit, dtype=float), axis=0))
                continue
            green_frames = [green_superpixel(load_bayer_raw(p, config), config.bayer_pattern) for p in group_files]
            median_green = np.median(np.stack(green_frames, axis=0), axis=0)
            plane, noise, _ = estimate_background_plane(median_green, border_px=config.border_px_green)
            corrected = np.clip(median_green - plane, 0, None)
            _, seed, _ = locate_primary_roi(corrected, noise, config)
            automatic_group_hints[key] = seed
    rows: list[dict[str, Any]] = []
    diagnostics: dict[str, dict[str, np.ndarray]] = {}
    for path in files:
        meta = parsed[path]
        key = (meta["focal"], meta["position_cm"])
        group_key = f"{key[0]}:{key[1]:g}"
        hint = center_hints.get(path.name, center_hints.get(group_key, automatic_group_hints.get(key)))
        metrics, arrays = analyze_spatial_frame(path, config=config, center_hint=hint)
        rows.append(metrics)
        if path.name in diagnostic_names:
            diagnostics[path.name] = arrays
    df = pd.DataFrame(rows).sort_values(["focal", "position_cm", "replicate"]).reset_index(drop=True)
    return df, diagnostics


def aggregate_spatial_metrics(frame_metrics: pd.DataFrame) -> pd.DataFrame:
    """Aggregate frame-level metrics across replicate exposures."""
    numeric_cols = [
        "d4sigma_x_um",
        "d4sigma_y_um",
        "d4sigma_major_um",
        "d4sigma_minor_um",
        "gaussian_w_major_um",
        "gaussian_w_minor_um",
        "primary_center_x_mm",
        "primary_center_y_mm",
        "primary_fraction",
        "secondary_peak_ratio",
        "gaussian_fit_nrmse",
        "saturation_fraction",
    ]
    available = [c for c in numeric_cols if c in frame_metrics]
    grouped = frame_metrics.groupby(["focal", "position_cm", "position_m"], dropna=False)
    rows: list[dict[str, Any]] = []
    for keys, group in grouped:
        row: dict[str, Any] = {
            "focal": keys[0],
            "position_cm": keys[1],
            "position_m": keys[2],
            "n_frames": len(group),
            "n_valid_for_m2": int(group["valid_for_m2"].sum()),
            "any_saturated": bool(group["saturated"].any()),
            "any_clipped": bool(group["clipped"].any()),
            "any_secondary_structure": bool(group["secondary_structure"].any()),
        }
        for col in available:
            row[f"{col}_mean"] = float(group[col].mean())
            row[f"{col}_std"] = float(group[col].std(ddof=1)) if len(group) > 1 else np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["focal", "position_cm"]).reset_index(drop=True)


def beam_radius_model(z_m: np.ndarray, w0_m: float, z0_m: float, m2: float, wavelength_m: float) -> np.ndarray:
    """Second-moment beam-radius propagation model."""
    z_m = np.asarray(z_m, dtype=float)
    divergence_term = m2 * wavelength_m / (np.pi * w0_m)
    return np.sqrt(w0_m**2 + (divergence_term * (z_m - z0_m)) ** 2)


def fit_beam_propagation(
    z_m: Sequence[float],
    radius_m: Sequence[float],
    wavelength_nm: float,
    radius_std_m: Sequence[float] | None = None,
    enforce_m2_ge_1: bool = True,
) -> dict[str, Any]:
    """Fit w(z) and return waist, Rayleigh range, divergence, and M^2.

    The routine adds a ``well_constrained`` flag. A numerically converged fit is
    not considered physically well constrained unless the measured points span
    both sides of the fitted waist.
    """
    z = np.asarray(z_m, dtype=float)
    w = np.asarray(radius_m, dtype=float)
    valid = np.isfinite(z) & np.isfinite(w) & (w > 0)
    z, w = z[valid], w[valid]
    if z.size < 4:
        raise ValueError("At least four valid propagation points are required")
    sigma = None
    if radius_std_m is not None:
        s = np.asarray(radius_std_m, dtype=float)[valid]
        positive = np.isfinite(s) & (s > 0)
        if positive.any():
            fallback = np.nanmedian(s[positive])
            sigma = np.where(positive, s, fallback)
    wavelength_m = wavelength_nm * 1e-9
    i_min = int(np.argmin(w))
    w0_guess = max(float(w[i_min]), 1e-9)
    z0_guess = float(z[i_min])
    span = max(float(np.ptp(z)), 1e-6)
    p0 = [w0_guess, z0_guess, 1.2]
    lower_m2 = 1.0 if enforce_m2_ge_1 else 0.0
    bounds = ([1e-9, z.min() - 3 * span, lower_m2], [10 * w.max(), z.max() + 3 * span, 100.0])

    def model(z_values: np.ndarray, w0: float, z0: float, m2: float) -> np.ndarray:
        return beam_radius_model(z_values, w0, z0, m2, wavelength_m)

    popt, pcov = curve_fit(
        model,
        z,
        w,
        p0=p0,
        sigma=sigma,
        absolute_sigma=sigma is not None,
        bounds=bounds,
        maxfev=30_000,
    )
    perr = np.sqrt(np.diag(pcov))
    w0, z0, m2 = map(float, popt)
    w0_err, z0_err, m2_err = map(float, perr)
    fitted = model(z, *popt)
    residuals = w - fitted
    rmse = float(np.sqrt(np.mean(residuals**2)))
    r2 = float(1.0 - np.sum(residuals**2) / np.sum((w - w.mean()) ** 2)) if np.ptp(w) > 0 else np.nan
    z_rayleigh = np.pi * w0**2 / (m2 * wavelength_m)
    divergence = m2 * wavelength_m / (np.pi * w0)
    n_left = int(np.sum(z < z0))
    n_right = int(np.sum(z > z0))
    well_constrained = bool(n_left >= 2 and n_right >= 2 and z.min() <= z0 <= z.max())
    return {
        "w0_m": w0,
        "w0_um": w0 * 1e6,
        "w0_std_m": w0_err,
        "w0_std_um": w0_err * 1e6,
        "z0_m": z0,
        "z0_cm": z0 * 100,
        "z0_std_m": z0_err,
        "m2": m2,
        "m2_std": m2_err,
        "rayleigh_range_m": z_rayleigh,
        "rayleigh_range_cm": z_rayleigh * 100,
        "far_field_divergence_rad": divergence,
        "far_field_divergence_mrad": divergence * 1e3,
        "rmse_m": rmse,
        "rmse_um": rmse * 1e6,
        "r2": r2,
        "n_points": int(z.size),
        "n_left_of_waist": n_left,
        "n_right_of_waist": n_right,
        "well_constrained": well_constrained,
        "covariance": pcov,
    }


# -----------------------------------------------------------------------------
# Temporal analysis (Keysight DSO CSV)
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class TemporalConfig:
    expected_rep_rate_hz: float = 200_000.0
    detector_bandwidth_hz: float | None = 12e6
    scope_bandwidth_hz: float | None = 1e9
    align_cycles: bool = True
    minimum_single_cycle_snr_for_alignment: float = 5.0
    minimum_cycles_for_folding: int = 4
    local_alignment_fraction_of_period: float = 0.04
    positive_polarity: bool | None = None
    detector_responsivity_a_per_w: float | None = None
    transimpedance_v_per_a_50ohm: float | None = None
    transimpedance_v_per_a_high_z: float | None = None
    optical_path_transmission: float = 1.0


def parse_scope_metadata(path: str | Path, max_rows: int = 260) -> dict[str, str]:
    """Read the metadata stored in the first three CSV columns."""
    path = Path(path)
    metadata: dict[str, str] = {}
    section = ""
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        for i, row in enumerate(reader):
            if i >= max_rows:
                break
            row = row + [""] * max(0, 6 - len(row))
            key = row[0].strip().rstrip(":")
            value = row[1].strip()
            if not key:
                continue
            if key.isupper() and not value:
                section = key
                continue
            full_key = f"{section}/{key}" if section else key
            # Preserve duplicate keys with a numeric suffix.
            candidate = full_key
            suffix = 2
            while candidate in metadata:
                candidate = f"{full_key} [{suffix}]"
                suffix += 1
            metadata[candidate] = value
    return metadata


def read_scope_csv(path: str | Path) -> tuple[pd.DataFrame, dict[str, str]]:
    """Read time, CH1, and CH3 from a DSO6104A CSV export."""
    path = Path(path)
    # The first CSV row contains instrument metadata. The second row contains
    # the data-column headers in columns 4-6.
    data = pd.read_csv(
        path,
        usecols=[3, 4, 5],
        skiprows=1,
        encoding="utf-8-sig",
        engine="c",
    )
    data.columns = ["time_s", "ch1_v", "ch3_v"]
    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    data = data.dropna().reset_index(drop=True)
    if len(data) < 2:
        raise ValueError(f"No waveform samples found in {path}")
    metadata = parse_scope_metadata(path)
    return data, metadata


def parse_pulse_filename(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    name = path.name
    mode = "picosecond" if "Mode_ps" in name else ("65_ns" if "Mode_65ns" in name else ("1_ns" if "Mode_1ns" in name else "unknown"))
    time_div_match = re.search(r"TimeDiv_([\d.]+)(ns|us|ms)", name, re.I)
    if time_div_match:
        scale = {"ns": 1e-9, "us": 1e-6, "ms": 1e-3}[time_div_match.group(2).lower()]
        time_div_s = float(time_div_match.group(1)) * scale
    else:
        time_div_s = np.nan
    rep_match = re.search(r"ExtTrigger_(\d+(?:\.\d+)?)(kHz|MHz|Hz)", name, re.I)
    rep_rate = np.nan
    if rep_match:
        mult = {"hz": 1.0, "khz": 1e3, "mhz": 1e6}[rep_match.group(2).lower()]
        rep_rate = float(rep_match.group(1)) * mult
    nd_match = re.search(r"ND_(\d+(?:_\d+)?)", name, re.I)
    termination = path.parents[1].name if len(path.parents) >= 2 else ""
    return {
        "file": name,
        "path": str(path),
        "mode": mode,
        "termination": "50ohm" if "50ohm" in termination else ("1Mohm" if "1Mohm" in termination else termination),
        "time_div_s": time_div_s,
        "nominal_rep_rate_hz": rep_rate,
        "nd_label": nd_match.group(1).replace("_", ".") if nd_match else "",
    }


def _crossing_time(t: np.ndarray, y: np.ndarray, start: int, stop: int, level: float, rising: bool) -> float:
    """Linear-interpolated crossing in a directed interval."""
    if rising:
        indices = np.where((y[start:stop] <= level) & (y[start + 1 : stop + 1] > level))[0]
    else:
        indices = np.where((y[start:stop] >= level) & (y[start + 1 : stop + 1] < level))[0]
    if indices.size == 0:
        return np.nan
    j = start + (indices[-1] if rising else indices[0])
    dy = y[j + 1] - y[j]
    if dy == 0:
        return float(t[j])
    return float(t[j] + (level - y[j]) * (t[j + 1] - t[j]) / dy)


def pulse_shape_metrics(time_s: Sequence[float], signal_v: Sequence[float]) -> dict[str, float]:
    """Calculate nonparametric pulse-shape metrics around the dominant positive peak."""
    t = np.asarray(time_s, dtype=float)
    y = np.asarray(signal_v, dtype=float)
    if t.size != y.size or t.size < 5:
        raise ValueError("time and signal must be equal-length arrays with at least five samples")
    baseline = float(np.median(y))
    yc = y - baseline
    pos = np.percentile(yc, 99.9)
    neg = -np.percentile(yc, 0.1)
    polarity = 1.0 if pos >= neg else -1.0
    yp = polarity * yc
    peak_index = int(np.argmax(yp))
    amplitude = float(yp[peak_index])
    noise = robust_sigma(yp[yp < np.percentile(yp, 60)])
    if amplitude <= 0:
        raise ValueError("No positive pulse excursion found")
    half = 0.5 * amplitude
    t_left = _crossing_time(t, yp, 0, max(peak_index, 1), half, rising=True)
    t_right = _crossing_time(t, yp, peak_index, t.size - 1, half, rising=False)
    fwhm = t_right - t_left if np.isfinite(t_left) and np.isfinite(t_right) else np.nan
    t10_left = _crossing_time(t, yp, 0, max(peak_index, 1), 0.1 * amplitude, rising=True)
    t90_left = _crossing_time(t, yp, 0, max(peak_index, 1), 0.9 * amplitude, rising=True)
    t90_right = _crossing_time(t, yp, peak_index, t.size - 1, 0.9 * amplitude, rising=False)
    t10_right = _crossing_time(t, yp, peak_index, t.size - 1, 0.1 * amplitude, rising=False)
    rise = t90_left - t10_left if np.isfinite(t10_left) and np.isfinite(t90_left) else np.nan
    fall = t10_right - t90_right if np.isfinite(t90_right) and np.isfinite(t10_right) else np.nan
    positive = np.clip(yp, 0, None)
    area = float(np.trapezoid(positive, t))
    equivalent_width = area / amplitude if amplitude > 0 else np.nan
    after = yp[peak_index + 1 :]
    undershoot = float(-np.min(after, initial=0.0))
    second_positive = float(np.max(after, initial=0.0))
    return {
        "polarity": polarity,
        "baseline_v": baseline,
        "noise_rms_v": float(noise),
        "snr_peak": amplitude / max(noise, 1e-15),
        "amplitude_v": amplitude,
        "peak_time_s": float(t[peak_index]),
        "fwhm_s": float(fwhm),
        "fwhm_ns": float(fwhm * 1e9),
        "rise_10_90_s": float(rise),
        "rise_10_90_ns": float(rise * 1e9),
        "fall_90_10_s": float(fall),
        "fall_90_10_ns": float(fall * 1e9),
        "area_v_s": area,
        "equivalent_width_s": float(equivalent_width),
        "equivalent_width_ns": float(equivalent_width * 1e9),
        "undershoot_ratio": undershoot / amplitude,
        "postpulse_positive_ratio": second_positive / amplitude,
    }


def fold_periodic_waveform(
    time_s: Sequence[float],
    signal_v: Sequence[float],
    period_s: float,
    align_cycles: bool = True,
    alignment_fraction: float = 0.04,
    minimum_single_cycle_snr_for_alignment: float = 5.0,
) -> dict[str, Any]:
    """Coherently average a periodic pulse train and quantify cycle stability."""
    t = np.asarray(time_s, dtype=float)
    y = np.asarray(signal_v, dtype=float)
    dt = float(np.median(np.diff(t)))
    samples_per_period = int(round(period_s / dt))
    if samples_per_period < 8:
        raise ValueError("Too few samples per repetition period")
    n_cycles = y.size // samples_per_period
    if n_cycles < 2:
        raise ValueError("Record is too short for periodic folding")
    cycles = y[: n_cycles * samples_per_period].reshape(n_cycles, samples_per_period).astype(float)
    # Subtract a robust baseline per cycle. Low duty cycle makes the median stable.
    baselines = np.median(cycles, axis=1)
    cycles -= baselines[:, None]
    mean0 = cycles.mean(axis=0)
    pos = np.percentile(mean0, 99.9)
    neg = -np.percentile(mean0, 0.1)
    polarity = 1.0 if pos >= neg else -1.0
    cycles *= polarity
    mean0 *= polarity
    reference_peak = int(np.argmax(mean0))
    shifts = np.zeros(n_cycles, dtype=int)
    valid = np.ones(n_cycles, dtype=bool)
    baseline_noise = np.array([robust_sigma(c[c < np.percentile(c, 60)]) for c in cycles])
    radius = max(2, int(alignment_fraction * samples_per_period))
    lo0 = max(0, reference_peak - radius)
    hi0 = min(samples_per_period, reference_peak + radius + 1)
    local_amplitude = np.max(cycles[:, lo0:hi0], axis=1) if hi0 > lo0 else np.max(cycles, axis=1)
    single_cycle_snr = float(np.nanmedian(local_amplitude / np.maximum(baseline_noise, 1e-15)))
    alignment_used = bool(align_cycles and single_cycle_snr >= minimum_single_cycle_snr_for_alignment)
    if alignment_used:
        for i, cycle in enumerate(cycles):
            lo = max(0, reference_peak - radius)
            hi = min(samples_per_period, reference_peak + radius + 1)
            local = cycle[lo:hi]
            if local.size == 0:
                valid[i] = False
                continue
            local_peak = lo + int(np.argmax(local))
            shifts[i] = reference_peak - local_peak
        max_shift = max(2, radius)
        valid &= np.abs(shifts) <= max_shift
        aligned = np.stack([np.roll(cycles[i], shifts[i]) for i in range(n_cycles) if valid[i]], axis=0)
    else:
        aligned = cycles
    mean = aligned.mean(axis=0)
    std = aligned.std(axis=0, ddof=1) if aligned.shape[0] > 1 else np.zeros_like(mean)
    sem = std / math.sqrt(max(aligned.shape[0], 1))
    peak_index = int(np.argmax(mean))
    phase = (np.arange(samples_per_period) - peak_index) * dt
    mean = np.roll(mean, samples_per_period // 2 - peak_index)
    std = np.roll(std, samples_per_period // 2 - peak_index)
    sem = np.roll(sem, samples_per_period // 2 - peak_index)
    phase = (np.arange(samples_per_period) - samples_per_period // 2) * dt
    # Cycle-level amplitudes and areas in a window centered on the unrolled peak.
    window_half = max(5, int(0.10 * samples_per_period))
    indices = (np.arange(reference_peak - window_half, reference_peak + window_half + 1) % samples_per_period).astype(int)
    amplitudes = np.max(aligned[:, indices], axis=1)
    positive = np.clip(aligned[:, indices], 0, None)
    areas = np.trapezoid(positive, dx=dt, axis=1)
    arrival_residuals = -shifts[valid] * dt
    return {
        "phase_s": phase,
        "mean_v": mean,
        "std_v": std,
        "sem_v": sem,
        "n_cycles_total": int(n_cycles),
        "n_cycles_used": int(aligned.shape[0]),
        "samples_per_period": samples_per_period,
        "dt_s": dt,
        "polarity": polarity,
        "alignment_used": alignment_used,
        "single_cycle_snr_median": single_cycle_snr,
        "cycle_amplitude_v": amplitudes,
        "cycle_area_v_s": areas,
        "arrival_residual_s": arrival_residuals,
        "baseline_per_cycle_v": baselines,
    }


def effective_bandwidth_hz(detector_bandwidth_hz: float | None, scope_bandwidth_hz: float | None) -> float | None:
    """RSS-combine independent Gaussian-equivalent bandwidth limitations."""
    values = [v for v in (detector_bandwidth_hz, scope_bandwidth_hz) if v is not None and v > 0]
    if not values:
        return None
    inverse_square = sum(1.0 / float(v) ** 2 for v in values)
    return 1.0 / math.sqrt(inverse_square)


def analyze_scope_record(
    path: str | Path,
    config: TemporalConfig = TemporalConfig(),
    keep_waveform: bool = False,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Analyze one oscilloscope record, using coherent folding when possible."""
    path = Path(path)
    data, metadata = read_scope_csv(path)
    meta = parse_pulse_filename(path)
    t = data["time_s"].to_numpy(float)
    y = data["ch1_v"].to_numpy(float)
    trigger = data["ch3_v"].to_numpy(float)
    dt = float(np.median(np.diff(t)))
    duration = float(t[-1] - t[0])
    expected_rate = meta["nominal_rep_rate_hz"] if np.isfinite(meta["nominal_rep_rate_hz"]) else config.expected_rep_rate_hz
    period = 1.0 / expected_rate
    use_fold = duration >= config.minimum_cycles_for_folding * period
    folded: dict[str, Any] | None = None
    if use_fold:
        folded = fold_periodic_waveform(
            t,
            y,
            period,
            align_cycles=config.align_cycles,
            alignment_fraction=config.local_alignment_fraction_of_period,
            minimum_single_cycle_snr_for_alignment=config.minimum_single_cycle_snr_for_alignment,
        )
        shape = pulse_shape_metrics(folded["phase_s"], folded["mean_v"])
        amp = np.asarray(folded["cycle_amplitude_v"], dtype=float)
        area = np.asarray(folded["cycle_area_v_s"], dtype=float)
        arrivals = np.asarray(folded["arrival_residual_s"], dtype=float)
        rep_rate_measured = expected_rate
        if arrivals.size >= 4:
            n = np.arange(arrivals.size, dtype=float)
            # Arrival residual slope represents a small period correction.
            slope, intercept = np.polyfit(n, arrivals, 1)
            measured_period = period + slope
            rep_rate_measured = 1.0 / measured_period if measured_period > 0 else expected_rate
            detrended_arrivals = arrivals - (slope * n + intercept)
            timing_jitter = float(np.std(detrended_arrivals, ddof=1))
        else:
            timing_jitter = np.nan
        stability = {
            "n_cycles": int(folded["n_cycles_used"]),
            "cycle_alignment_used": bool(folded["alignment_used"]),
            "single_cycle_snr_median": float(folded["single_cycle_snr_median"]),
            "amplitude_mean_v": float(np.mean(amp)),
            "amplitude_std_v": float(np.std(amp, ddof=1)) if amp.size > 1 else np.nan,
            "amplitude_cv_percent": float(100 * np.std(amp, ddof=1) / np.mean(amp)) if amp.size > 1 and np.mean(amp) != 0 else np.nan,
            "area_mean_v_s": float(np.mean(area)),
            "area_std_v_s": float(np.std(area, ddof=1)) if area.size > 1 else np.nan,
            "area_cv_percent": float(100 * np.std(area, ddof=1) / np.mean(area)) if area.size > 1 and np.mean(area) != 0 else np.nan,
            "timing_jitter_rms_s": timing_jitter,
            "timing_jitter_rms_ps": timing_jitter * 1e12,
            "rep_rate_measured_hz": float(rep_rate_measured),
        }
    else:
        baseline = np.median(y)
        centered = y - baseline
        if config.positive_polarity is False:
            centered = -centered
        elif config.positive_polarity is None:
            if -np.percentile(centered, 0.1) > np.percentile(centered, 99.9):
                centered = -centered
        peak_index = int(np.argmax(centered))
        half_window = min(len(y) // 2, max(20, int(0.45 * len(y))))
        lo, hi = max(0, peak_index - half_window), min(len(y), peak_index + half_window + 1)
        shape = pulse_shape_metrics(t[lo:hi] - t[peak_index], y[lo:hi])
        stability = {
            "n_cycles": 1,
            "cycle_alignment_used": False,
            "single_cycle_snr_median": shape["snr_peak"],
            "amplitude_mean_v": shape["amplitude_v"],
            "amplitude_std_v": np.nan,
            "amplitude_cv_percent": np.nan,
            "area_mean_v_s": shape["area_v_s"],
            "area_std_v_s": np.nan,
            "area_cv_percent": np.nan,
            "timing_jitter_rms_s": np.nan,
            "timing_jitter_rms_ps": np.nan,
            "rep_rate_measured_hz": np.nan,
        }
    bw = effective_bandwidth_hz(config.detector_bandwidth_hz, config.scope_bandwidth_hz)
    rise_limit = 0.35 / bw if bw else np.nan
    impulse_fwhm = 0.44 / bw if bw else np.nan
    mode_expected_s = {"picosecond": 0.1e-9, "1_ns": 1e-9, "65_ns": 65e-9}.get(meta["mode"], np.nan)
    unresolved = bool(np.isfinite(mode_expected_s) and np.isfinite(rise_limit) and mode_expected_s < 0.5 * rise_limit)
    trigger_range = float(np.ptp(trigger))
    if meta["termination"] == "50ohm":
        transimpedance = config.transimpedance_v_per_a_50ohm
    else:
        transimpedance = config.transimpedance_v_per_a_high_z
    responsivity = config.detector_responsivity_a_per_w
    transmission = config.optical_path_transmission
    calibration_ready = bool(
        transimpedance is not None
        and transimpedance > 0
        and responsivity is not None
        and responsivity > 0
        and transmission > 0
    )
    if calibration_ready:
        calibration_factor = transimpedance * responsivity * transmission
        peak_power_w = shape["amplitude_v"] / calibration_factor
        pulse_energy_j = shape["area_v_s"] / calibration_factor
        average_power_w = pulse_energy_j * (stability["rep_rate_measured_hz"] if np.isfinite(stability["rep_rate_measured_hz"]) else expected_rate)
    else:
        peak_power_w = pulse_energy_j = average_power_w = np.nan
    result: dict[str, Any] = {
        **meta,
        "sample_count": int(len(data)),
        "sample_interval_s": dt,
        "sample_interval_ns": dt * 1e9,
        "sample_rate_gsa_s": 1.0 / dt / 1e9,
        "record_duration_s": duration,
        "record_duration_us": duration * 1e6,
        "trigger_range_v": trigger_range,
        "folded": use_fold,
        **shape,
        **stability,
        "detector_bandwidth_hz_assumed": config.detector_bandwidth_hz,
        "scope_bandwidth_hz_assumed": config.scope_bandwidth_hz,
        "combined_bandwidth_hz": bw,
        "instrument_rise_time_s_est": rise_limit,
        "instrument_rise_time_ns_est": rise_limit * 1e9 if np.isfinite(rise_limit) else np.nan,
        "gaussian_impulse_fwhm_s_est": impulse_fwhm,
        "gaussian_impulse_fwhm_ns_est": impulse_fwhm * 1e9 if np.isfinite(impulse_fwhm) else np.nan,
        "expected_mode_duration_s": mode_expected_s,
        "pulse_unresolved_by_detector": unresolved,
        "scope_model": metadata.get("Model", "DSO6104A"),
        "scope_serial": metadata.get("Serial Number", ""),
        "scope_start_time": metadata.get("Start Time", ""),
        "detector_responsivity_a_per_w": responsivity,
        "transimpedance_v_per_a_used": transimpedance,
        "optical_path_transmission": transmission,
        "absolute_power_calibrated": calibration_ready,
        "peak_power_w": peak_power_w,
        "pulse_energy_j": pulse_energy_j,
        "average_power_w": average_power_w,
    }
    waveform: dict[str, Any] | None = None
    if keep_waveform:
        waveform = {
            "time_s": t,
            "ch1_v": y,
            "ch3_v": trigger,
            "folded": folded,
            "metadata": metadata,
        }
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
    df = pd.DataFrame(rows).sort_values(["termination", "mode", "record_duration_s"]).reset_index(drop=True)
    return df, waveforms


# -----------------------------------------------------------------------------
# Spectral analysis (HORIBA text exports)
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
    minimum_points_for_voigt: int = 7


def read_horiba_spectrum(path: str | Path) -> tuple[pd.DataFrame, dict[str, str]]:
    """Parse a HORIBA spectrum with #metadata and decimal commas."""
    path = Path(path)
    metadata: dict[str, str] = {}
    wavelength: list[float] = []
    counts: list[float] = []
    with path.open("r", encoding="latin-1", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                text = line[1:]
                if "=" in text:
                    key, value = text.split("=", 1)
                    metadata[key.strip()] = value.strip().replace(",", ".")
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                wavelength.append(float(parts[0].replace(",", ".")))
                counts.append(float(parts[1].replace(",", ".")))
            except ValueError:
                continue
    if len(wavelength) < 3:
        raise ValueError(f"No numeric spectrum found in {path}")
    frame = pd.DataFrame({"wavelength_nm": wavelength, "counts": counts})
    return frame, metadata


def parse_spectrum_filename(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    name = path.name
    lower = name.lower()
    wavelength_family = 532 if "532nm" in lower else (1064 if "1064nm" in lower else np.nan)
    if "picosecond" in lower:
        mode = "picosecond"
    elif "65nano" in lower:
        mode = "65_ns"
    elif "cw" in lower:
        mode = "cw"
    else:
        mode = "unknown"
    grating_match = re.search(r"Grad_(\d+)", name, re.I)
    slit_match = re.search(r"Slit_(\d+)um", name, re.I)
    acq_match = re.search(r"Acq_(\d+)(?:_(\d+))?s", name, re.I)
    rtd_match = re.search(r"RTD_(\d+(?:_\d+)?)s", name, re.I)
    replicate_match = re.search(r"_\((\d+)\)\.txt$", name)
    if acq_match:
        acq_s = float(acq_match.group(1)) + (float("0." + acq_match.group(2)) if acq_match.group(2) else 0.0)
    elif rtd_match:
        acq_s = float(rtd_match.group(1).replace("_", "."))
    else:
        acq_s = np.nan
    base_name = re.sub(r"\s*\(\d+\)(?=\.txt$)", "", name)
    return {
        "file": name,
        "path": str(path),
        "wavelength_family_nm": wavelength_family,
        "mode": mode,
        "grating_g_mm": int(grating_match.group(1)) if grating_match else np.nan,
        "slit_um_filename": int(slit_match.group(1)) if slit_match else np.nan,
        "acquisition_s_filename": acq_s,
        "replicate": int(replicate_match.group(1)) if replicate_match else 0,
        "replicate_group": base_name,
    }


def asymmetric_least_squares(
    y: Sequence[float],
    lam: float = 1e6,
    p: float = 0.005,
    iterations: int = 12,
) -> np.ndarray:
    """Asymmetric least-squares baseline (Eilers/Boelens style)."""
    y = np.asarray(y, dtype=float)
    length = y.size
    if length < 3:
        return np.zeros_like(y)
    difference = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(length - 2, length), format="csc")
    penalty = lam * (difference.T @ difference)
    weights = np.ones(length)
    baseline = np.zeros_like(y)
    for _ in range(iterations):
        wmat = sparse.spdiags(weights, 0, length, length)
        baseline = spsolve(wmat + penalty, weights * y)
        weights = p * (y > baseline) + (1.0 - p) * (y <= baseline)
    return np.asarray(baseline)


def edge_polynomial_baseline(
    x: Sequence[float],
    y: Sequence[float],
    edge_fraction: float = 0.12,
    order: int = 1,
    iterations: int = 5,
) -> np.ndarray:
    """Fit a robust polynomial baseline to the two ends of a spectrum.

    This is preferable to a smoothness-penalty baseline for broadband laser
    spectra, where ALS can incorrectly classify the optical envelope itself as
    background.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n_edge = max(order + 3, int(round(edge_fraction * x.size)))
    n_edge = min(n_edge, max(order + 3, x.size // 2))
    idx = np.r_[np.arange(n_edge), np.arange(x.size - n_edge, x.size)]
    xe, ye = x[idx], y[idx]
    keep = np.isfinite(xe) & np.isfinite(ye)
    coeff = np.polyfit(xe[keep], ye[keep], order)
    for _ in range(iterations):
        fit = np.polyval(coeff, xe)
        resid = ye - fit
        sigma = robust_sigma(resid[keep])
        if not np.isfinite(sigma) or sigma <= 0:
            break
        # Reject bright spectral tails/outliers more strongly than negative noise.
        keep = np.isfinite(resid) & (resid > -4 * sigma) & (resid < 2.5 * sigma)
        if keep.sum() <= order + 1:
            break
        coeff = np.polyfit(xe[keep], ye[keep], order)
    return np.polyval(coeff, x)


def estimate_spectral_baseline(x: Sequence[float], y: Sequence[float], config: SpectralConfig) -> np.ndarray:
    method = config.baseline_method.lower()
    if method in {"edge", "edge_linear", "polynomial", "edge_polynomial"}:
        return edge_polynomial_baseline(
            x,
            y,
            edge_fraction=config.baseline_edge_fraction,
            order=config.baseline_polynomial_order,
        )
    if method in {"als", "asls", "asymmetric_least_squares"}:
        return asymmetric_least_squares(
            y,
            lam=config.baseline_lambda,
            p=config.baseline_asymmetry,
            iterations=config.baseline_iterations,
        )
    raise ValueError(f"Unknown spectral baseline method: {config.baseline_method}")


def _fwhm_crossings(x: np.ndarray, y: np.ndarray, peak_index: int | None = None) -> tuple[float, float, float]:
    """Return left crossing, right crossing, and FWHM around one peak."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if peak_index is None:
        peak_index = int(np.argmax(y))
    peak = float(y[peak_index])
    if peak <= 0:
        return np.nan, np.nan, np.nan
    half = peak / 2.0
    left_candidates = np.where(y[:peak_index] <= half)[0]
    right_candidates = np.where(y[peak_index:] <= half)[0]
    if left_candidates.size == 0 or right_candidates.size == 0:
        return np.nan, np.nan, np.nan
    il = int(left_candidates[-1])
    ir = int(peak_index + right_candidates[0])

    def interp(i0: int, i1: int) -> float:
        dy = y[i1] - y[i0]
        if dy == 0:
            return float(x[i0])
        return float(x[i0] + (half - y[i0]) * (x[i1] - x[i0]) / dy)

    xl = interp(il, il + 1)
    xr = interp(ir - 1, ir)
    return xl, xr, xr - xl


def _voigt_model(x: np.ndarray, offset: float, slope: float, area: float, center: float, sigma: float, gamma: float) -> np.ndarray:
    profile = voigt_profile(x - center, sigma, gamma)
    return offset + slope * (x - center) + area * profile


def fit_voigt_peak(x: np.ndarray, y: np.ndarray, peak_index: int, fwhm_guess: float) -> dict[str, float]:
    """Fit a Voigt peak in a local window. Intended only for adequately sampled peaks."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    center0 = float(x[peak_index])
    step = float(np.median(np.diff(x)))
    width = max(float(fwhm_guess), 5 * step)
    mask = np.abs(x - center0) <= 4.0 * width
    if mask.sum() < 7:
        mask = np.abs(x - center0) <= 8.0 * step
    xf, yf = x[mask], y[mask]
    baseline0 = float(np.percentile(yf, 10))
    amplitude0 = max(float(yf.max() - baseline0), 1e-12)
    sigma0 = max(width / 2.355, step / 2)
    gamma0 = max(width / 4.0, step / 4)
    area0 = amplitude0 / max(float(voigt_profile(0.0, sigma0, gamma0)), 1e-12)
    p0 = [baseline0, 0.0, area0, center0, sigma0, gamma0]
    bounds = (
        [-np.inf, -np.inf, 0.0, xf.min(), step / 20, 0.0],
        [np.inf, np.inf, np.inf, xf.max(), 5 * width, 5 * width],
    )
    popt, pcov = curve_fit(_voigt_model, xf, yf, p0=p0, bounds=bounds, maxfev=30_000)
    perr = np.sqrt(np.diag(pcov))
    offset, slope, area, center, sigma, gamma = map(float, popt)
    fwhm_l = 2.0 * gamma
    fwhm_g = 2.354820045 * sigma
    fwhm_v = 0.5346 * fwhm_l + math.sqrt(0.2166 * fwhm_l**2 + fwhm_g**2)
    fit = _voigt_model(xf, *popt)
    rmse = float(np.sqrt(np.mean((yf - fit) ** 2)))
    return {
        "voigt_center_nm": center,
        "voigt_center_std_nm": float(perr[3]),
        "voigt_sigma_nm": sigma,
        "voigt_gamma_nm": gamma,
        "voigt_fwhm_nm": fwhm_v,
        "voigt_area": area,
        "voigt_rmse": rmse,
    }


def analyze_spectrum(
    path: str | Path,
    config: SpectralConfig = SpectralConfig(),
    keep_arrays: bool = False,
) -> tuple[dict[str, Any], dict[str, np.ndarray] | None]:
    path = Path(path)
    frame, metadata = read_horiba_spectrum(path)
    meta = parse_spectrum_filename(path)
    x = frame["wavelength_nm"].to_numpy(float)
    y = frame["counts"].to_numpy(float)
    baseline = estimate_spectral_baseline(x, y, config)
    corrected = np.clip(y - baseline, 0, None)
    peak_index = int(np.argmax(corrected))
    peak = float(corrected[peak_index])
    step = float(np.median(np.diff(x)))
    xl, xr, fwhm = _fwhm_crossings(x, corrected, peak_index)
    noise = robust_sigma((y - baseline)[corrected < np.percentile(corrected, 50)])
    threshold = max(3.0 * max(noise, 0.0), 0.003 * peak)
    positive_mask = corrected > threshold
    # Keep the contiguous region around the dominant peak for line moments.
    labels, _ = ndimage.label(positive_mask)
    region_id = labels[peak_index]
    peak_region = labels == region_id if region_id > 0 else positive_mask
    if peak_region.sum() < 3:
        peak_region = corrected > 0.001 * peak
    weights = corrected * peak_region
    weight_sum = float(weights.sum())
    centroid = float(np.sum(x * weights) / weight_sum) if weight_sum > 0 else np.nan
    rms_width = float(np.sqrt(np.sum(weights * (x - centroid) ** 2) / weight_sum)) if weight_sum > 0 else np.nan
    equivalent_width = float(np.trapezoid(corrected, x) / peak) if peak > 0 else np.nan
    prominence = max(config.peak_prominence_fraction * peak, 5.0 * max(noise, 0.0))
    peak_trace = ndimage.gaussian_filter1d(corrected, sigma=0.8)
    peaks, properties = find_peaks(peak_trace, prominence=prominence, distance=max(1, int(round(0.5 * max(fwhm, step) / step))))
    points_per_fwhm = fwhm / step if np.isfinite(fwhm) else np.nan
    peak_nm = float(x[peak_index])
    delta_nu = C_LIGHT * (fwhm * 1e-9) / (peak_nm * 1e-9) ** 2 if np.isfinite(fwhm) else np.nan
    resolving_power = peak_nm / fwhm if np.isfinite(fwhm) and fwhm > 0 else np.nan
    saturated = bool(np.max(y) >= 0.98 * config.saturation_dn)
    undersampled = bool(np.isfinite(points_per_fwhm) and points_per_fwhm < 3.0)
    voigt: dict[str, float] = {
        "voigt_center_nm": np.nan,
        "voigt_center_std_nm": np.nan,
        "voigt_sigma_nm": np.nan,
        "voigt_gamma_nm": np.nan,
        "voigt_fwhm_nm": np.nan,
        "voigt_area": np.nan,
        "voigt_rmse": np.nan,
    }
    if (
        np.isfinite(fwhm)
        and points_per_fwhm >= config.minimum_points_for_voigt
        and meta["mode"] != "picosecond"
    ):
        try:
            voigt = fit_voigt_peak(x, y - baseline, peak_index, fwhm)
        except Exception as exc:  # Fit failure should not discard nonparametric metrics.
            warnings.warn(f"Voigt fit failed for {path.name}: {exc}")
    slit_metadata = metadata.get("Front entrance slit", "")
    acq_metadata = metadata.get("Acq. time (s)", "")
    grating_metadata = metadata.get("Grating", "")
    result: dict[str, Any] = {
        **meta,
        "sample_count": int(x.size),
        "wavelength_min_nm": float(x.min()),
        "wavelength_max_nm": float(x.max()),
        "sampling_step_nm": step,
        "peak_wavelength_nm": peak_nm,
        "centroid_wavelength_nm": centroid,
        "fwhm_left_nm": xl,
        "fwhm_right_nm": xr,
        "fwhm_nm": float(fwhm),
        "fwhm_frequency_hz": float(delta_nu),
        "fwhm_frequency_ghz": float(delta_nu / 1e9),
        "rms_width_nm": rms_width,
        "equivalent_width_nm": equivalent_width,
        "resolving_power_lambda_over_fwhm": resolving_power,
        "points_per_fwhm": float(points_per_fwhm),
        "peak_corrected_counts": peak,
        "peak_raw_counts": float(np.max(y)),
        "baseline_at_peak_counts": float(baseline[peak_index]),
        "noise_sigma_counts": float(noise),
        "snr_peak": peak / max(noise, 1e-15),
        "n_detected_peaks": int(len(peaks)),
        "saturated": saturated,
        "undersampled_peak": undersampled,
        "grating_metadata": float(grating_metadata) if str(grating_metadata).replace(".", "", 1).isdigit() else grating_metadata,
        "slit_um_metadata": float(slit_metadata) if str(slit_metadata).replace(".", "", 1).isdigit() else slit_metadata,
        "acquisition_s_metadata": float(acq_metadata) if str(acq_metadata).replace(".", "", 1).isdigit() else acq_metadata,
        "spectrometer_center_nm_metadata": metadata.get("Spectro (nm)", ""),
        "detector_temperature_C_metadata": metadata.get("Detector temperature (°C)", metadata.get("Detector temperature (Â°C)", "")),
        "acquired_metadata": metadata.get("Acquired", ""),
        **voigt,
    }
    arrays = None
    if keep_arrays:
        arrays = {
            "wavelength_nm": x,
            "raw_counts": y,
            "baseline_counts": baseline,
            "corrected_counts": corrected,
            "detected_peak_indices": peaks,
        }
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
        ["wavelength_family_nm", "mode", "grating_g_mm", "slit_um_filename", "file"]
    ).reset_index(drop=True)
    return frame, arrays


def transform_limited_duration(fwhm_frequency_hz: float, pulse_shape: str = "gaussian") -> float:
    """Return the transform-limited intensity FWHM for a chosen pulse shape."""
    if not np.isfinite(fwhm_frequency_hz) or fwhm_frequency_hz <= 0:
        return np.nan
    shape = pulse_shape.lower().replace("^", "")
    tbp = 0.441 if shape in {"gaussian", "gauss"} else (0.315 if shape in {"sech2", "sech²", "sech"} else None)
    if tbp is None:
        raise ValueError("pulse_shape must be 'gaussian' or 'sech2'")
    return tbp / fwhm_frequency_hz


def summarize_replicate_spectra(metrics: pd.DataFrame) -> pd.DataFrame:
    """Summarize repeated spectra grouped by filename stem without _(n)."""
    grouped = metrics.groupby("replicate_group", dropna=False)
    rows: list[dict[str, Any]] = []
    for group_name, group in grouped:
        rows.append(
            {
                "replicate_group": group_name,
                "n": len(group),
                "wavelength_family_nm": group["wavelength_family_nm"].iloc[0],
                "mode": group["mode"].iloc[0],
                "grating_g_mm": group["grating_g_mm"].iloc[0],
                "slit_um_filename": group["slit_um_filename"].iloc[0],
                "peak_wavelength_mean_nm": float(group["peak_wavelength_nm"].mean()),
                "peak_wavelength_std_nm": float(group["peak_wavelength_nm"].std(ddof=1)) if len(group) > 1 else np.nan,
                "fwhm_mean_nm": float(group["fwhm_nm"].mean()),
                "fwhm_std_nm": float(group["fwhm_nm"].std(ddof=1)) if len(group) > 1 else np.nan,
                "peak_counts_cv_percent": float(100 * group["peak_corrected_counts"].std(ddof=1) / group["peak_corrected_counts"].mean()) if len(group) > 1 and group["peak_corrected_counts"].mean() != 0 else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["wavelength_family_nm", "mode", "grating_g_mm"]).reset_index(drop=True)


__all__ = [
    "C_LIGHT",
    "CameraConfig",
    "TemporalConfig",
    "SpectralConfig",
    "find_data_root",
    "ensure_directory",
    "minmax_envelope",
    "parse_spatial_filename",
    "load_bayer_raw",
    "green_superpixel",
    "load_png_green",
    "estimate_background_plane",
    "fit_rotated_gaussian",
    "second_moments",
    "analyze_spatial_frame",
    "analyze_spatial_directory",
    "aggregate_spatial_metrics",
    "beam_radius_model",
    "fit_beam_propagation",
    "parse_scope_metadata",
    "read_scope_csv",
    "parse_pulse_filename",
    "pulse_shape_metrics",
    "fold_periodic_waveform",
    "effective_bandwidth_hz",
    "analyze_scope_record",
    "analyze_scope_directory",
    "read_horiba_spectrum",
    "parse_spectrum_filename",
    "asymmetric_least_squares",
    "edge_polynomial_baseline",
    "estimate_spectral_baseline",
    "fit_voigt_peak",
    "analyze_spectrum",
    "analyze_spectral_directory",
    "transform_limited_duration",
    "summarize_replicate_spectra",
]
