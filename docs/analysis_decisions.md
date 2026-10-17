# Analysis Decisions, Definitions, and Claim Boundaries

This document records why the code makes each major choice. It is intended to prevent a numerical output from being mistaken for a physically stronger claim.

## 1. General reporting policy

1. A calculated number is not automatically a reportable laser property.
2. Independent repeatability is reported when repeated acquisitions exist. It is not relabeled as total uncertainty.
3. Curve-fit covariance is conditional on the model, weights, selected data, and coordinate calibration. It is not a substitute for systematic uncertainty.
4. No confidence level, coverage factor, or calibration uncertainty is invented.
5. Extra digits remain in machine-readable CSVs so results can be regression-tested. Narrative rounding does not imply traceable accuracy.

## 2. Spatial definitions and choices

### 2.1 Bayer handling

The native files contain BayerRG12 values stored in 16-bit words. The analysis averages the two green photosites in each 2 x 2 cell. This creates an interpolation-free green superpixel and avoids demosaicing artifacts in centroids and second moments. The effective analysis pitch is 4.8 um because the physical camera pitch is 2.4 um and no optical relay was present.

### 2.2 Background treatment

No matched native dark frames were recorded. A robust planar background is estimated from perimeter pixels after excluding high-valued edge outliers. This is preferable to subtracting a single arbitrary constant, but it cannot remove fixed-pattern noise, exposure-dependent black-level offsets, or coherent stray light. The absence of a dark frame remains a systematic limitation.

### 2.3 Primary-lobe selection

The experimenter identified the highest-intensity spot as the real beam. The code therefore chooses the dominant connected lobe and carries a consistent center hint through replicates. Secondary spots, rings, and reflections do not automatically invalidate the selected lobe, but they are quantified through secondary-peak, primary-energy-fraction, and full-field metrics.

### 2.4 Gaussian fitting

A rotated two-dimensional Gaussian is used to localize the primary lobe and define a stable aperture. It is not used to claim that the field is a fundamental Gaussian mode. The primary quantitative widths are second-moment D4sigma diameters. Gaussian 1/e2 radii are exported as model-dependent diagnostics.

### 2.5 Saturation

A Bayer cell is flagged when any of its four native photosites reaches 99.5 percent of the 12-bit full scale. Saturation inside the selected primary aperture excludes that frame from the strict unsaturated caustic set. Saturation elsewhere remains a field-quality warning.

### 2.6 M-squared model

For a second-moment radius w(z), the fit uses

```text
w(z)^2 = w0^2 + [M^2 * lambda / (pi * w0)]^2 * (z - z0)^2
```

Principal-major and principal-minor axes are fit separately. Sensor x/y axes and Gaussian-fit radii are included as sensitivity variants. The project requires adequate points, both sides of the fitted waist, sufficient span in fitted Rayleigh ranges, unsaturated waist information, and known coordinate quality before a value is marked reportable.

The coverage test in the code is a conservative project heuristic, not a declaration of ISO 11146 compliance. The present scan fails it.

### 2.7 No-lens slopes and centroid motion

A linear radius-versus-z slope over the supplied no-lens range is called a local expansion slope. It is not called far-field divergence because the far-field condition was not demonstrated and z = 0 is arbitrary.

Centroid changes between different z positions describe a propagation trajectory and alignment. They are not pointing jitter because the camera was not held at one plane for a time-series acquisition.

### 2.8 Supplemental PNG data

The older PNG scans are 8-bit, demosaiced, and heavily saturated. They are analyzed in a separate table and never merged into native 12-bit results. Their purpose is qualitative coverage context only.

## 3. Temporal definitions and choices

### 3.1 Detector model

The detector is the legacy PDA36A. The analysis does not substitute PDA36A2 specifications. At 0 dB, the nominal values used are 10 MHz bandwidth, 0.750 kV/A into 50 ohm, and 1.51 kV/A into high impedance. The near-infrared photodiode rise-time warning is retained.

### 3.2 Proper termination

Only records whose exported metadata confirm CH1 at 50 ohm are used for quantitative waveform conclusions. The 1 Mohm data demonstrate loading and reflections but are not pooled with the 50 ohm data.

### 3.3 Electrical versus optical width

The analysis reports nonparametric FWHM and 10-90/90-10 transition times of the measured electrical waveform. It does not equate these with the optical pulse duration.

The common step-response relation tr approximately 0.35/B is used only as a bandwidth scale. A -3 dB bandwidth does not uniquely define an impulse-response FWHM. Therefore the project removes the earlier Gaussian-impulse approximation and performs no guessed deconvolution.

### 3.4 Coherent folding

For long records, the CH1 waveform is divided into nominal 5 us periods and coherently averaged. A per-cycle median removes the DC baseline. Cycle alignment is enabled only when the single-cycle signal-to-noise ratio is adequate. The exported amplitude and area CV values include detector and oscilloscope noise.

### 3.5 Trigger edge detection

CH3 edges contain sample-level noise and ringing. A naive midpoint-crossing detector can count multiple transitions. The final detector smooths over nine samples, uses low/high hysteretic thresholds, interpolates the midpoint crossing, and enforces a 0.70-period refractory interval. This yields the physically expected 200 kHz result.

### 3.6 Timing interpretation

CH3 edge residuals and CH3-to-CH1 peak delays include the function generator, trigger threshold, laser, photodetector, oscilloscope sampling, cable delays, and lack of deskew calibration. They are never labeled intrinsic GSL timing jitter.

### 3.7 Power and energy

Absolute detector-derived power remains disabled because responsivity and optical attenuation were not supplied. Pulse energy instead comes from the independently measured average power through Pavg/f. Peak power remains blocked until optical pulse duration and temporal shape are established.

## 4. Spectral definitions and choices

### 4.1 Baseline

A robust linear baseline is estimated from the wavelength-window edges. This is less likely than a high-order polynomial to create artificial narrow features or distort the broad picosecond spectra. Baseline-corrected negative values are retained for diagnostics, while positive-clipped values are used only where an integral requires nonnegative intensity.

### 4.2 Center metrics

Peak wavelength is the sampled maximum. Centroid is the first moment of the positive baseline-corrected spectrum. The peak is more sensitive to sampling; the centroid is more sensitive to baseline, window, and system-response correction. Both are reported.

### 4.3 FWHM and sampling adequacy

FWHM is found by linearly interpolating the half-maximum crossings around the dominant peak. The number of sampling intervals per FWHM is calculated explicitly. Fewer than three points per FWHM is flagged as undersampled for line-shape/linewidth interpretation.

### 4.4 Instrument line-spread function

No nominal catalog resolution is subtracted in quadrature. Catalog resolution is configuration specific and does not replace a line-spread-function measurement under the same grating, slit, wavelength, detector, and focus conditions. Narrow widths are therefore reported as measured system widths, not intrinsic laser linewidths.

### 4.5 Slit label and high-count warning

A filename slit value of 0 um is interpreted only as the software's minimum/fully closed setting. It is not treated as a known physical zero-width aperture.

The 85 percent full-scale threshold is a conservative project quality-control warning. It is not claimed to be the Synapse Plus saturation or linearity specification. It prevents high-count spectra from being silently pooled with ordinary repeats.

### 4.6 Repeatability grouping

Spectra are grouped as repeats only when wavelength family, mode, grating, slit label, and acquisition settings match and the files are time contiguous. Distinct high-count states are separated. Sequence standard deviation quantifies short-term acquisition repeatability only.

### 4.7 Wavelength calibration

The user calibrated two dominant visible mercury lines for all three gratings. This is regionally useful near 532 nm. Near 1064 nm, the result is an extrapolation because no near-infrared reference or saved residual curve is available. The GSL106A nominal 1064 +/- 10 nm band is used only as a plausibility test, not as a calibration standard.

### 4.8 Transform limits

The code converts wavelength-domain FWHM to optical-frequency FWHM and evaluates

```text
Gaussian: delta_t = 0.441 / delta_nu
sech^2:   delta_t = 0.315 / delta_nu
```

These are ideal transform-limited intensity-duration lower bounds. They do not determine the actual pulse without spectral phase information, autocorrelation with an assumed shape, or a phase-sensitive measurement such as FROG.

## 5. Claim-status vocabulary

- `supported`: the dataset and chain support the stated measurement-level claim.
- `supported_with_caveat`: the number is useful, but a named calibration or systematic caveat remains.
- `provisional`: the result is a model sensitivity or diagnostic, not a final specification.
- `not_supported`: the measurement chain cannot establish the physical quantity.

These labels appear in output tables and notebook narrative so downstream users do not have to infer validity from file names alone.


## Power and calibration reporting decisions

- Inequality symbols from the laboratory notes are preserved as bounds; lower-bound powers are never converted into equality estimates.
- The 65 ns power is represented as the observed approximate interval 330-340 uW, with 333 uW retained as the original approximate reading.
- Pulse energy is derived from average power and 200 kHz repetition rate. Peak power is blocked unless a trustworthy optical temporal profile is available.
- The PM100D/S121C manufacturer specification is not substituted for a measurement-specific uncertainty budget. No combined or expanded uncertainty is invented.
- The Hg calibration at 435.8 and 546.1 nm is treated as a visible-region offset calibration. Absolute accuracy near 1064 nm is explicitly unverified without a near-IR reference.
- The warm-up of more than two minutes is documented as a possible contributor to drift and session-to-session reproducibility.
- All supplied files are catalogued. Files that do not map one-to-one to a numerical table remain represented in the source inventory and domain contact sheets.
