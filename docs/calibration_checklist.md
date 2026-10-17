# Remaining Calibration and Remeasurement Checklist

This version contains only unresolved items. The setup facts already supplied by the experimenter are stored in `experiment_metadata.json` and are not asked again.

## 1. Inputs that can be added without repeating the experiment

### Pulsed average power

Average powers for the three GSL106A pulsed modes are recorded in `outputs/average_power_inputs.csv` (PM100D/S121C, after the collimator, before attenuation). The picosecond and 1 ns values are lower bounds. Still worth recording: the sensor's calibration date, the measurement range used, and the number and spread of repeated readings, which would turn those bounds into values with an uncertainty.

### Spatial hardware metadata

Record, when available:

- exact Thorlabs 100 mm plano-convex lens part number;
- lens orientation relative to the incoming beam;
- translation-stage model, scale resolution, backlash procedure, and estimated z uncertainty;
- camera black-level, gamma/digital-shift, and any saved camera configuration file;
- exact reflective ND filter part numbers and angles of incidence.

These do not rescue saturated frames, but they improve the physical scale and systematic-error discussion.

### Spectrometer hardware metadata

Record, when available:

- exact Synapse Plus detector model/serial number and detector temperature;
- grating blaze wavelength, diffraction order, turret position, and order-sorting filter for each acquisition;
- actual minimum entrance/exit slit gap corresponding to filename `0 um`;
- LabSpec version and any saved calibration/configuration files;
- mercury line wavelengths used, fitted calibration equation, residuals, and date/time.

## 2. Measurements required for a final spatial M-squared result

1. Warm the COMPASS laser until power and beam position are stable; document the stabilization criterion rather than using an arbitrary time alone.
2. Record a native BayerRG12 dark frame at the same exposure, gain, black level, and sensor temperature.
3. Adjust reflective attenuation so that no primary-lobe pixel exceeds roughly 80 to 90 percent of full scale at any z.
4. Keep the same beam path and identified primary lobe for all positions.
5. Measure at least ten well-distributed z positions, including several within about one Rayleigh range of the waist and several beyond about two Rayleigh ranges on each side after an initial estimate.
6. Record the lens plane, stage reading, direction, repeatability/backlash procedure, and z uncertainty.
7. Acquire at least three independent frames per position without changing exposure or attenuation. More repeats near the waist are useful only after saturation is eliminated.
8. Record an aperture/clipping check for the lens, ND mount, mirror, and sensor at every position.

## 3. Measurements required for optical pulse duration

### 65 ns mode

- Measure or obtain the complete PDA36A impulse response at approximately 1064 nm and the exact 0 dB/50 ohm configuration, or use a detector with substantially faster documented response.
- Confirm detector linearity and absence of output saturation.
- Calibrate CH3-to-CH1 delay/deskew only when absolute delay is required.
- A measured detector transfer function is required before deconvolution; a catalog bandwidth alone is insufficient.

### 1 ns mode

Use a detector and oscilloscope chain with a rise time comfortably below 1 ns, including adequate detector bandwidth and a controlled 50 ohm path. Record the system impulse response with a faster reference source when an absolute pulse width is required.

### Picosecond mode

Use an optical autocorrelator, cross-correlator, streak camera, or FROG-class measurement appropriate to the wavelength and expected duration. A silicon amplified detector with a 10 MHz nominal bandwidth cannot establish this duration.

## 4. Measurements required for intrinsic linewidth or absolute wavelength

1. Calibrate the wavelength scale with one or more known lines near 1064 nm and save fit residuals across the relevant detector window.
2. Recheck calibration after grating changes to quantify turret repeatability/backlash.
3. Measure a source substantially narrower than the expected laser line under each grating/slit/detector configuration to obtain the instrument line-spread function.
4. Record actual entrance and exit slit gaps rather than relying only on filename labels.
5. Use sufficient spectral sampling: target substantially more than three points per expected FWHM for line-shape fitting.
6. Prevent detector nonlinearity and document the usable count range with a linearity or attenuation series.
7. Record and apply the wavelength-dependent grating/detector response when comparing broadband centroids or shapes across gratings.
8. Record polarization or scramble/control it when grating-efficiency effects matter.

## 5. Inputs required for detector-derived absolute optical power

The current project intentionally does not reconstruct absolute optical power from CH1. That calculation would require:

- PDA36A responsivity at the actual wavelength for the specific detector;
- 0 dB transimpedance for the actual scope load;
- measured transmission of every ND filter at 1064 nm and its angle of incidence;
- confirmation that the full beam falls on the active area;
- detector linearity and peak/average power limits;
- cable/load transfer function if waveform area is used;
- uncertainty and calibration date for each term.

An independent calibrated average-power measurement is the preferred route for pulse energy in this dataset.

## 6. Minimum additions before final thesis tables

For each reported result, add a column labeled `claim scope` with one of:

```text
intrinsic laser property
sensor-plane relative measurement
measurement-chain electrical result
instrument-limited spectral result
model sensitivity only
```

Keep Type A repeatability separate from systematic/calibration terms. Do not assign an expanded uncertainty or confidence interval unless the required coverage convention and component uncertainties are actually available.
