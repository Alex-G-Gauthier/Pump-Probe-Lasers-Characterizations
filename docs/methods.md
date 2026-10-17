# Methods, Uncertainty and Reproducing the Analysis

## Analysis Principles

### Spatial Analysis

The Basler RAW files are treated as 3088 x 2064 arrays of right-aligned 12-bit values stored in 16-bit words. The two green photosites in each 2 x 2 Bayer cell are averaged to produce a 1544 x 1032 interpolation-free green image. With direct interception and no relay optics, one green superpixel corresponds to 4.8 um in each transverse direction.

A robust planar background is estimated from image-perimeter pixels because matched native dark frames were not recorded. The highest-intensity connected lobe is treated as the intended beam, following the experimenter's rule. A rotated Gaussian fit localizes that lobe and defines a stable aperture, but the reported beam diameters are second-moment quantities rather than Gaussian-fit widths. Reflections and satellite spots remain visible through full-field metrics and quality flags.

The focused scan is fit to the paraxial second-moment propagation model only as a sensitivity analysis. The project never promotes a numerically converged fit to a final M-squared result unless coverage, saturation, clipping, repeatability, and position metadata all support the claim.

### Temporal Analysis

The DSO6104A metadata embedded in each CSV are parsed and checked. The quantitative files confirm CH1 at 50 ohm, CH3 at 1 Mohm, bandwidth limit off, normal real-time acquisition, CH3 trigger, DC trigger coupling, and zero stored channel skew.

The legacy PDA36A at 0 dB is represented by the documented nominal 10 MHz bandwidth, 0.750 kV/A transimpedance into 50 ohm, and 1.51 kV/A into high impedance. The manufacturer also warns that the photodiode rise time can further limit effective bandwidth in the near infrared. A 10 MHz bandwidth provides a step-rise-time scale of about 35 ns; it does not uniquely specify an impulse-response FWHM. No arbitrary Gaussian deconvolution is performed.

Long records are coherently folded at the intended 200 kHz period. CH3 rising edges are detected with a hysteretic threshold and a refractory interval to prevent ringing from being counted as extra pulses. Timing residuals are explicitly labeled as function-generator plus oscilloscope plus laser plus detector-chain metrics, not intrinsic laser jitter.

### Spectral Analysis

HORIBA text files are parsed directly, including decimal commas and acquisition metadata. A robust edge-based linear baseline is used. Peak wavelength, centroid, nonparametric FWHM, RMS width, equivalent width, frequency bandwidth, resolving power, and points per FWHM are reported.

A narrow measured FWHM is not called an intrinsic linewidth without an instrument line-spread function measured at the same wavelength, grating, slit, detector, and acquisition settings. Fewer than three sampled points across FWHM is treated as undersampled for quantitative line-shape interpretation. The project's 85 percent full-scale threshold is a conservative data-quality warning, not a Synapse Plus linearity specification.

### Power and Pulse Energy

Average power was measured with a Thorlabs PM100D console and S121C sensor, recorded in `outputs/average_power_inputs.csv`. Pulse energy is

```text
pulse energy = average power / measured repetition rate
```

Peak power is intentionally not calculated because the optical pulse duration and shape are not established by the current detector chain.

## Uncertainty and Significant Figures

No confidence level or uncertainty budget is invented. The project uses:

- observed standard deviations across independent repeated frames or time-contiguous repeated spectra;
- conditional curve-fit covariance only as a fit diagnostic;
- explicit systematic limitations for missing dark frames, z uncertainty, wavelength-calibration residuals, spectral response, detector impulse response, and power-chain calibration.

CSV files retain extra digits to make numerical regeneration testable. Narrative results normally use at least three significant figures, but displayed digits are not presented as uncertainty.

## Current Limitations That Affect Interpretation

- A laser warm-up of more than two minutes may be insufficient for thermal stability.
- No matched native Bayer dark frames were recorded; camera black level and gamma were not recorded.
- Five of eight focused positions contain saturated pixels inside the selected primary aperture.
- Spatial z-position uncertainty and focusing-lens part number/orientation were not recorded.
- PDA36A responsivity at 1064 nm and ND-filter transmission at 1064 nm were not supplied.
- The picosecond and 1 ns average powers are lower bounds, so their pulse energies are lower bounds.
- Spectrometer calibration residuals, grating blaze/order-sorting details, actual minimum slit gap, detector variant, and near-1064 nm reference were not recorded.
- The spectrometer system response was not corrected.

These omissions do not make all data unusable. They define which statements are relative measurement-chain observations and which intrinsic laser claims remain unsupported.

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m ipykernel install --user --name laser-characterization
```

## Data Placement

The repository does not include the roughly 1.1 GB measurement archive. Extract the original data so that one of these layouts exists:

```text
Pump-Probe-Lasers-Characterizations/
    data/
        532 and 1064/
            532nm_Spatial/
            Pulse_Data/
            Spectral_Data/
```

or

```text
parent_directory/
    Pump-Probe-Lasers-Characterizations/
    532 and 1064/
```

`find_data_root()` searches both layouts. For a custom location, pass it explicitly:

```bash
python run_all_analyses.py --data-root "/absolute/path/to/532 and 1064"
```

The file checksums of the original archive are in `docs/source_data_manifest.sha256`.

## Reproduce the Outputs

From the project directory, with the raw data in place:

```bash
python run_all_analyses.py --force
cd notebooks
for nb in 0*.ipynb; do jupyter nbconvert --to notebook --execute --inplace "$nb"; done
cd ..
python validate_project.py --data-root "/absolute/path/to/532 and 1064"
```

The `--force` option reparses every raw file. Without it, the runner reuses existing per-file metric tables and regenerates derived summaries.
