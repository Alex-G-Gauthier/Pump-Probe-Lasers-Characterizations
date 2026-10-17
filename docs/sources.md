# Instrument and Method Sources

Accessed 2026-08-30 unless otherwise noted. Manufacturer facts are separated from user-reported experimental metadata.

## Manufacturer and archived instrument documentation

1. Basler, **acA3088-57uc product page and product documentation**. Used for the Sony IMX178 sensor, 3088 x 2064 resolution, 2.4 um x 2.4 um physical pixel pitch, color sensor, and 10/12-bit pixel-depth capability.
   - https://www.baslerweb.com/en-us/shop/aca3088-57uc/
   - https://docs.baslerweb.com/aca3088-57uc

2. HORIBA Scientific, **iHR550 Imaging Spectrometer** brochure and iHR Series product page. Used for 550 mm focal length, f/6.4 aperture, three-grating turret, and the configuration dependence of nominal resolution specifications. The brochure states that its 0.025 nm resolution, +/-0.20 nm wavelength accuracy, and related values are specified for a 1200 groove/mm grating at 435 nm; those values are not treated as the experiment's uncertainty at 532 or 1064 nm.
   - https://static.horiba.com/fileadmin/Horiba/Products/Scientific/Optical_Components_and_OEM/iHR_Series/iHR550.pdf
   - https://www.horiba.com/int/scientific/products/detail/action/show/Product/ihr-series-1590/

3. Thorlabs, **Picosecond and Nanosecond Gain-Switched Lasers**. Used for GSL106A nominal center wavelength 1064 +/- 10 nm, HI1060 output fiber, external-trigger capability, and trigger-input context.
   - https://www.thorlabs.com/picosecond-and-nanosecond-gain-switched-lasers

4. Thorlabs, **PDA36A legacy Si Switchable Gain Detector User Guide**, Rev. E, 2017, and archived catalog page. Used for the legacy model's 350 to 1100 nm range, 0 dB nominal 10.0 MHz bandwidth, 1.51 x 10^3 V/A high-impedance gain, 0.75 x 10^3 V/A gain into 50 ohm, internal 50 ohm series output resistor, and warning that near-infrared photodiode rise time may reduce effective bandwidth.
   - Official discontinued-product page: https://www.thorlabs.com/zh/item/pda36a-ec
   - Official archived catalog: https://www.thorlabs.com/catalogpages/V21/1573.PDF
   - Accessible manual mirror used to cross-check the specification table: https://wiki.physics.wisc.edu/yavuz/images/5/59/PDA36A-Manual.pdf

5. Keysight, **DSO6104A Oscilloscope product page**. Used for the nominal 1 GHz bandwidth and 4 GSa/s maximum sample rate. The actual exported records, not only the catalog, determine the 0.25 ns sample interval and stored acquisition settings.
   - https://www.keysight.com/us/en/product/DSO6104A/oscilloscope-1-ghz-4-analog-channels.html

6. Keysight, **33120A Function Generator User's Guide**. Used only as instrument provenance. Trigger amplitude, offset, waveform, and output-termination setting were not recoverable from the supplied data and are not invented.
   - https://www.keysight.com/us/en/assets/9018-04436/user-manuals/9018-04436.pdf

7. Coherent, **COMPASS 315M-100**. A current exact manufacturer datasheet was not located in the accessible Coherent resource index during this revision. The model, 532 nm CW operation, and 130 mW output are therefore stored as user-reported experimental metadata. The 130 mW value is not used as a verified catalog specification or as a measured power in calculations.

## Method references and conventions

- Second-moment beam diameters and paraxial caustic notation follow the conventions commonly associated with ISO 11146 beam-width and propagation analysis. This project calls the implementation `ISO-style` and does not claim formal compliance because acquisition geometry, saturation, sampling, and uncertainty requirements are not all met.
- The Gaussian and sech-squared time-bandwidth products used for transform-limit lower bounds are 0.441 and 0.315, respectively, for intensity FWHM quantities.
- The exact SI value of the speed of light, 299792458 m/s, is used for wavelength-to-frequency conversion.

## User-reported metadata used directly

The following are not manufacturer specifications; they are experimenter-provided facts recorded in `experiment_metadata.json`:

- COMPASS 315M-100 operated as a separate 532 nm CW source, user-reported 130 mW.
- GSL106A fiber output operated in picosecond, 1 ns, 65 ns, and CW spectral modes at intended 200 kHz for pulsed acquisitions.
- Warm-up of more than two minutes.
- Direct camera interception; no camera/relay lens; BayerRG12; automatic gain/exposure and color correction disabled.
- Highest-intensity spatial lobe designated as the real beam.
- 100 mm Thorlabs plano-convex focusing lens; positive z farther from the lens.
- Legacy PDA36A at 0 dB; DC-coupled scope input; BNC cables shorter than 1 m; no external termination; CH3 function-generator output; no beam splitter.
- Reflective ND filters used for the 1 ns and 65 ns temporal measurements.
- iHR550 front entrance, free-space coupling, minimum/fully closed slit interpretation for `0 um`, Hg calibration of two dominant visible lines, LabSpec signal-dark correction, binning 1, and one accumulation.

All unknown settings remain explicitly null or `not recorded` rather than being filled from a similar instrument model.


## Additional manufacturer documentation (power and calibration)

- Thorlabs, PM100D Compact Power and Energy Meter Console, official product page and manual.
- Thorlabs, S121C Standard Photodiode Power Sensor, official product page and calibration documentation.
- Thorlabs, F280APC-1064 Fiber Collimation Package, official product page.
- Coherent, COMPASS 315M series legacy product documentation where available.

Manufacturer specifications are used to describe nominal instrument capability only. They do not replace experiment-specific calibration, uncertainty, or traceability.
