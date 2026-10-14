# Results Summary and Confidence Assessment

## Executive conclusion

The dataset supports a useful, careful characterization, but it does not support every common headline specification. The strongest conclusions are:

1. The 532 nm selected primary lobe changes size reproducibly with axial position and clearly passes through a focused minimum.
2. A numerical caustic fit converges, but a final M-squared value is not defensible from this scan because the waist region is saturated and the scan does not extend far enough in fitted Rayleigh-range units.
3. The GSL106A trigger rate is confirmed at approximately 200.000 kHz.
4. The recorded temporal traces are electrical responses of the PDA36A-cable-oscilloscope chain. They do not directly measure the picosecond or nominal 1 ns optical pulse widths, and they do not calibrate the 65 ns optical width.
5. The narrow spectral peaks are configuration- and sampling-limited. The broad picosecond-mode spectra are well sampled and support a spectral bandwidth and ideal transform-limit lower bound.
6. A visible mercury calibration makes the 532 nm center values plausible, but the absolute near-1064 nm scale remains unvalidated without a near-infrared reference.

## 1. Spatial characterization of the 532 nm COMPASS laser

### Measurement basis

- 36 native BayerRG12 RAW frames.
- Four no-focusing-lens positions: 0, 10, 20, and 27 cm from an arbitrary reference.
- Eight focused positions: 5, 8, 9, 10, 11, 12, 15, and 20 cm from the 100 mm focusing-lens plane.
- Three independently acquired frames per position.
- Direct sensor interception; no camera lens or relay; nominal magnification 1.
- Basler physical pixel pitch 2.4 um; green-superpixel pitch 4.8 um after interpolation-free 2 x 2 Bayer reduction.
- The highest-intensity lobe is the intended beam. Other spots and rings are retained as secondary-structure diagnostics.

### Focused selected-lobe dimensions

Second-moment diameters are reported for the selected primary lobe. The values below are means of three independent frames.

| z from focusing lens (cm) | Major D4sigma (um) | Minor D4sigma (um) | Ellipticity | Primary saturation |
|---:|---:|---:|---:|---|
| 5 | 477.533 | 399.220 | 1.196 | No |
| 8 | 454.544 | 368.813 | 1.232 | Yes |
| 9 | 419.788 | 354.743 | 1.183 | Yes |
| 10 | 418.953 | 347.324 | 1.206 | Yes |
| 11 | 423.164 | 357.965 | 1.182 | Yes |
| 12 | 425.607 | 332.457 | 1.280 | Yes |
| 15 | 511.642 | 382.211 | 1.339 | No |
| 20 | 569.434 | 470.321 | 1.211 | No |

The smallest measured minor-axis diameter is 332.457 um at 12 cm, but that position contains saturated primary-lobe pixels. The smallest unsaturated focused minor-axis diameter among the supplied positions is 382.211 um at 15 cm. Neither value should automatically be called the true waist diameter because the discrete z sampling and saturation can move the observed minimum.

### Provisional caustic fits

Weighted fits to all eight focused positions produce the following sensitivity values:

| Axis | Fitted waist radius w0 (um) | Fitted waist z0 (cm) | Provisional M-squared | Fitted zR (cm) | Fitted divergence (mrad) | R-squared |
|---|---:|---:|---:|---:|---:|---:|
| Principal major D4sigma/2 | 209.386 | 10.3428 | 2.53867 | 10.1983 | 2.05315 | 0.838990 |
| Principal minor D4sigma/2 | 171.871 | 11.2820 | 1.81428 | 9.61473 | 1.78758 | 0.943313 |

These numbers are **not reportable final M-squared values**. The reasons are physical, not merely statistical:

- only eight focused z positions are available;
- the 8, 9, 10, 11, and 12 cm positions are saturated inside the selected primary aperture;
- only the 5, 15, and 20 cm positions are strict unsaturated candidates;
- the measured scan spans only about 0.52 to 0.65 fitted zR on one side and 0.91 to 0.95 fitted zR on the other;
- z-position uncertainty, stage backlash, exact lens part number, and lens orientation were not recorded;
- the major-axis model fit is visibly less adequate than the minor-axis fit.

The correct report language is: *the focused scan is consistent with a non-ideal, elliptical beam and gives a provisional principal-axis M-squared sensitivity range of about 1.81 to 2.54, but the present acquisition does not meet the project's validity conditions for a final M-squared claim.*

### No-lens expansion

The no-lens second-moment radii increase almost linearly over 0 to 27 cm. Local slopes are 1.08531 mrad for the principal major radius and 1.05455 mrad for the principal minor radius. These are high-quality measured-range slopes, not automatically the far-field divergences: the arbitrary z origin and absence of a demonstrated far-field regime prevent that stronger interpretation.

### Spatial confidence

| Quantity | Confidence | Reason |
|---|---|---|
| Physical transverse scale | High with manufacturer-spec caveat | Direct interception and documented 2.4 um camera pixels; no relay magnification |
| Relative selected-lobe size versus z | Moderate to high | Three independent frames per point and clear propagation trend |
| Absolute size near the focused minimum | Moderate to low | Saturation, no native dark frame, unknown black level/gamma |
| Ellipticity of selected lobe | Moderate | Stable primary selection, but coherent reflections and aperture choice matter |
| Time-domain pointing jitter | Unsupported | The images are different z positions, not a fixed-plane time series |
| Final M-squared | Unsupported | Inadequate unsaturated caustic coverage and missing z uncertainty |

The supplemental 8-bit PNG scans are analyzed separately and never merged with the RAW data. They contain more z positions but are demosaiced and heavily saturated, so they are useful only for qualitative propagation context.

## 2. Temporal characterization of the 1064 nm GSL106A

### Verified measurement chain

- Legacy Thorlabs PDA36A, not PDA36A2.
- Gain setting 0 dB.
- Nominal 10 MHz bandwidth at 0 dB.
- Nominal transimpedance 0.750 kV/A into 50 ohm and 1.51 kV/A into high impedance.
- Quantitative CH1 records use the oscilloscope's internal 50 ohm input.
- CH3 monitors the Agilent 33120A output split by a BNC T to the GSL106A trigger input.
- Scope exports report bandwidth limit off, normal real-time acquisition, CH3 trigger, DC trigger coupling, and no stored skew correction.
- Long records have 0.25 ns sample spacing, corresponding to 4 GSa/s.

### Long-record electrical metrics

| Mode label | Electrical FWHM (ns) | 10-90 rise (ns) | 90-10 fall (ns) | Coherent amplitude (mV) | Amplitude CV (%) | Area CV (%) | Measured repetition rate (Hz) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Picosecond | 20.8896 | 21.6170 | 25.5227 | 2.210 | 13.5448 | 7.90290 | 200000.383722 |
| 1 ns | 32.3680 | 18.7456 | 20.6640 | 17.062 | 14.4601 | 12.0378 | 200000.295107 |
| 65 ns | 49.0250 | 32.6072 | 23.7330 | 92.469 | 2.54305 | 5.10384 | 200000.370400 |

The repetition-rate result is robust. A simple midpoint-crossing algorithm initially over-counts noisy/ringing CH3 transitions; the final analysis uses hysteresis and a 0.70-period refractory interval, yielding 40 or 100 edges as expected in the long records and agreeing with the 200 kHz spectral fundamental.

The pulse-width interpretation is intentionally conservative:

- **Picosecond mode:** unresolved. A 20.9 ns electrical FWHM cannot be interpreted as a picosecond optical duration.
- **Nominal 1 ns mode:** unresolved. A 32.4 ns electrical FWHM is dominated by the detector-chain response.
- **Nominal 65 ns mode:** partially resolved but detector-convolved. The 49.0 ns electrical FWHM is not a calibrated optical width and should not be deconvolved from a guessed Gaussian response.

The strong alternating lobes after the main peak show ringing in the measurement chain. The 1 Mohm records alter the amplitude and transfer function and are retained only as a termination demonstration, not as quantitative pulse-shape measurements.

### Timing and stability

The hysteretic CH3 edge-fit residual is 343 to 464 ps RMS in the long records. CH3-to-electrical-peak repeatability is about 2.03 ns for the 1 ns record and 2.92 ns for the 65 ns record; the low-SNR picosecond record gives 45.7 ns scatter and is not usable for delay repeatability.

None of these values is intrinsic GSL106A jitter. They combine the function generator, trigger threshold, scope sampling, laser response, detector response, cable delays, and lack of channel deskew. The 0.25 ns sampling interval alone is much larger than a picosecond-scale source-jitter specification.

Amplitude and area CV values are also system-level stability metrics. The high CV in the fast modes is strongly affected by their small electrical amplitude and detector/scope noise.

### Pulse energy and peak power

Pulse energy can be calculated from an independent average-power measurement:

```text
E_pulse = P_average / f_rep
```

Numerical average powers were not supplied, so the project leaves pulse energies blank in `outputs/average_power_inputs.csv`. Peak power is not supported because a trustworthy optical pulse duration and shape are missing.

## 3. Spectral characterization

### 532 nm CW spectra

| Grating (groove/mm) | Acquisition (s) | Peak (nm) | Centroid (nm) | Measured FWHM (nm) | Sampling step (nm) | Points/FWHM |
|---:|---:|---:|---:|---:|---:|---:|
| 300 | 0.5 | 531.952 | 531.985766 | 0.234848 | 0.149 | 1.576 |
| 300 | 1.0 | 531.988 | 532.097468 | 0.319945 | 0.149 | 2.147 |
| 1200 | 0.5 | 532.102 | 532.095846 | 0.060464 | 0.034 | 1.778 |
| 1800 | 0.5 | 532.228 | 532.218954 | 0.036266 | 0.021 | 1.727 |

Every peak has fewer than 2.2 points across FWHM. The numerical interpolation is reproducible, but the line shape and intrinsic linewidth are not resolved. The change in measured width with grating is direct evidence that the result is instrument/sampling dependent.

The centers are broadly consistent with a nominal 532 nm source. Because the Hg calibration residuals and uncertainty were not saved, differences at the 0.1 nm scale should not be over-interpreted as laser drift or a true source offset.

### 1064 nm CW slit series

For the 1200 groove/mm grating:

| Nominal slit label (um) | Peak (nm) | Measured FWHM (nm) | Points/FWHM |
|---:|---:|---:|---:|
| 10 | 1060.58 | 0.055884 | 1.863 |
| 50 | 1060.58 | 0.061156 | 2.039 |
| 100 | 1060.61 | 0.111227 | 3.708 |
| 200 | 1060.69 | 0.224608 | 7.487 |

The approximately fourfold width increase between the 10 and 200 um settings is strong empirical evidence of slit-limited instrumental broadening. It prevents an intrinsic narrow-line linewidth claim without a measured line-spread function.

The filename label `0 um` is treated as minimum/fully closed separation, not a physically calibrated zero-width slit.

### 65 ns spectral acquisition states

The 300 groove/mm, 10 um series contains two distinct acquisition states:

- an earlier lower-count sequence with FWHM around 0.517 to 0.580 nm;
- a later high-count sequence with FWHM around 2.28 to 2.38 nm and four files above the project's conservative 85 percent full-scale warning threshold.

These files are not pooled as ordinary repeats. The threshold is not a manufacturer linearity limit; it signals that detector linearity and baseline behavior near the exported full scale were not established.

### Picosecond-mode spectra

Across the picosecond files:

- peak wavelength range: 1052.89 to 1055.20 nm;
- centroid range: 1052.943 to 1053.672 nm;
- spectral FWHM range: 8.49208 to 9.60616 nm;
- frequency FWHM range: 2296.51 to 2591.49 GHz;
- Gaussian transform-limit lower bound: 0.170173 to 0.192030 ps;
- sech-squared transform-limit lower bound: 0.121552 to 0.137165 ps.

These broad spectra have 61 to 284 sampled points per FWHM and are therefore much better resolved than the narrow lines. The bandwidth claim is supported with the caveat that the spectrometer/grating/CCD response and polarization dependence were not corrected.

The transform-limit values are not measured pulse durations. They are ideal lower bounds derived from the intensity spectrum under assumed pulse shapes and zero spectral phase. The actual pulse can be longer because of chirp or arbitrary spectral phase.

### Absolute wavelength confidence near 1064 nm

The narrow CW and 65 ns peaks near 1059 to 1061 nm lie within the GSL106A manufacturer's nominal 1064 +/- 10 nm band and are physically plausible. The broad picosecond peaks near 1053 to 1055 nm include a 1200 groove/mm value slightly below that nominal band. Neither observation proves an incorrect laser wavelength or an incorrect spectrometer. The calibration used two visible Hg lines and no near-1064 nm reference, so wavelength-model extrapolation, turret/grating effects, order sorting, and detector-response effects remain unquantified.

## 4. Direct answer: which values can be trusted?

| Result | Can it be used? | Proper interpretation |
|---|---|---|
| 532 nm selected-lobe dimensions in each RAW frame | Yes, with caveats | Direct sensor-plane second moments; saturation flags must accompany waist-region values |
| 532 nm relative propagation trend | Yes | Reproducible selected-lobe trend over the measured z coordinates |
| Provisional M-squared values 1.81 to 2.54 | Only as sensitivity analysis | Demonstrates model scale and axis dependence; not a final source specification |
| 200.000 kHz GSL repetition rate | Yes | Electrical trigger-chain repetition rate |
| 20.9, 32.4, and 49.0 ns FWHM | Yes as electrical widths only | Do not label them picosecond, 1 ns, or 65 ns optical pulse durations |
| CH3 edge residual or trigger-to-peak scatter | Yes as chain repeatability | Not intrinsic laser jitter |
| 532 nm narrow FWHM | No as intrinsic linewidth | Undersampled and instrument limited |
| 1064 nm narrow FWHM | No as intrinsic linewidth | Slit dependent and instrument line-spread function unknown |
| Picosecond spectral FWHM 8.49 to 9.61 nm | Yes, with response caveat | Broadband spectral intensity width |
| 0.170 to 0.192 ps Gaussian transform limit | Yes as a lower bound | Not an actual measured temporal pulse width |
| Absolute near-1064 nm peak position | Not yet | Requires a near-IR reference or calibration residuals spanning the region |
| Pulse energy | After average powers are entered | Use average power divided by measured repetition rate |
| Peak power | No | Optical duration/shape and full calibration are unavailable |

## 5. Most valuable follow-up measurements

1. Repeat the focused spatial scan with attenuation chosen to eliminate primary-lobe saturation at every z. Record a matched native dark frame, stage resolution/backlash, lens part number/orientation, and at least ten positions extending beyond two estimated Rayleigh ranges on both sides of the waist.
2. Use a detector with a measured response appropriate to the 65 ns mode. Use a fast photodiode/sampling system for the 1 ns mode and an autocorrelator, cross-correlator, or FROG-class instrument for the picosecond mode.
3. Enter the independently measured average powers for each pulsed mode to obtain pulse energy. Record power-meter model/calibration and measurement repeatability.
4. Calibrate the spectrometer with one or more traceable lines near 1064 nm and save calibration residuals. Measure an instrument line-spread function at each grating/slit setting used for linewidth claims.
5. For broadband spectral shape, record grating blaze/order-sorting configuration and apply a characterized wavelength-dependent system-response correction.

## 6. Reporting language recommended for a thesis or laboratory report

Use statements such as:

- "The selected 532 nm beam lobe exhibited a focused second-moment diameter minimum in the 10 to 12 cm region; however, saturation in that region prevented a standards-grade M-squared determination."
- "The GSL106A was externally triggered at 200.000 kHz. The measured waveforms are detector-chain electrical responses; the PDA36A did not resolve the picosecond or 1 ns optical pulses."
- "The narrow spectral features were instrument limited and undersampled. The picosecond mode exhibited an 8.49 to 9.61 nm spectral FWHM, corresponding to an ideal Gaussian transform-limit lower bound of 0.170 to 0.192 ps."
- "Absolute wavelength accuracy near 1064 nm was not independently validated because calibration used visible mercury lines only."


# Final power, calibration, and completeness update (Version 3)

## Power measurements and pulse energies

Power was measured with a Thorlabs PM100D console and S121C sensor. For the GSL106A, the sensor was placed in free space immediately after the F280APC-1064 fiber collimator and before any reflective ND attenuation. The notation is preserved exactly: `>` is a strict lower bound and `~` denotes an approximate value.

| Source/mode | Reported average power | Derived quantity |
|---|---:|---:|
| COMPASS 532 nm CW | >130 mW | Photon flux >3.482e+17 s^-1 at 532 nm |
| GSL106A 1064 nm CW | >16 mW | Photon flux >8.570e+16 s^-1 at nominal 1064 nm |
| GSL106A ps, 200 kHz | >1.35 uW | Pulse energy >6.75 pJ; >3.616e+07 photons/pulse |
| GSL106A 1 ns, 200 kHz | >9.4 uW | Pulse energy >47.0 pJ; >2.517e+08 photons/pulse |
| GSL106A 65 ns, 200 kHz | ~333 uW; observed 330-340 uW | Pulse energy 1.65-1.70 nJ; 8.838e+09-9.106e+09 photons/pulse |

These pulse energies follow directly from E_pulse = P_avg / f_rep. Peak powers are intentionally not reported as measurements because the optical pulse durations and pulse shapes were not resolved by the PDA36A/oscilloscope chain.

## Wavelength calibration consequence

Each grating was checked with Hg lines at 435.8 nm and 546.1 nm by applying wavelength offsets until the peaks matched to four significant figures. This makes the 532 nm region much better anchored than the 1064 nm region. With no traceable line near 1064 nm, the absolute center wavelengths around 1059-1061 nm remain plausible but unverified. The broadband picosecond spectral width is more robust than its absolute center because a nearly constant wavelength offset has little effect on FWHM.

## Completeness

The original archive SHA-256 is `c58900255d515fab38bb8883aa2cf632de8d63608b10f0956d7ba76f622d03f2`. The archive contains 119 non-directory entries. The complete file-level catalog, CRC values, available SHA-256 values, domain assignments, and traceability notes are provided in `outputs/source_data_inventory.csv` and `outputs/file_to_result_traceability.csv`. No supplied acquisition file was intentionally discarded.
