# Stellar irradiation

The flux that enters the top interface of the atmosphere comes from the stellar spectrum in `sflux_file`, given at the stellar surface, diluted by the geometric factor (R_star / a)^2 and interpolated onto the wavelength bins of the model (`build_atm.read_sflux`). The [overview](index.md) lists every case of the tab.

## Top-of-atmosphere flux

This case is an input-to-output identity of the top boundary on a real run: the flux in the top interface must equal the input spectrum times the dilution. It is not a radiative-transfer benchmark. The test runs the default model on 20 layers with photochemistry (the HD 189733b setup: `sflux-HD189_Moses11.txt`, R_star = 0.805 R_sun, a = 0.03142 au, a dilution of 0.0142). It reads the stellar file with numpy and interpolates it with `np.interp`, independently of the `scipy.interpolate.interp1d` call in VULCAN; the geometric factor is computed in the test from the configuration and the constants of `vulcan.phy_const`.

![Top-of-atmosphere flux against wavelength and its relative difference from the diluted input spectrum](../assets/validation/toa_flux_diluted_stellar_spectrum.svg#only-light){ width="80%" }
![Top-of-atmosphere flux against wavelength and its relative difference from the diluted input spectrum](../assets/validation/toa_flux_diluted_stellar_spectrum_dark.svg#only-dark){ width="80%" }

Top: the diluted input spectrum (line, every 25th of the 2588 bins from 2 to 654 nm) and the flux in the top interface of the run (dots, every 100th bin). Bottom: their relative difference at every 50th bin and the tolerance (dashed); bins that agree exactly are drawn at 1e-18. Over all bins the largest relative difference is 2.4e-16, and 1803 of the 2588 bins agree exactly.

- Reference: the stellar input file times (R_star / a)^2, interpolated with `np.interp`
- Tolerance: `rtol=1e-12` on every bin
- Test id: `tests/test_integration_quick.py::test_top_of_atmosphere_flux_is_diluted_stellar_spectrum`
- Command: `make` in `fastchem_vulcan/`, then `pytest tests/test_integration_quick.py::test_top_of_atmosphere_flux_is_diluted_stellar_spectrum` (add `--validation-plots` to write the plot pair)
