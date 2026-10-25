# Stellar irradiation

## What this page checks

This page checks the starlight that enters the top of the model atmosphere. VULCAN reads the stellar spectrum from the file `sflux_file`, which gives the flux at the surface of the star. It dilutes this flux by the geometric factor $(R_\star / a)^2$, with $R_\star$ the stellar radius and $a$ the orbital distance, and interpolates it onto the wavelength bins of the model (`build_atm.read_sflux`). Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. If this part of the code is wrong, every photolysis rate is wrong, and so is the chemistry of the upper atmosphere. Each case below says why the check exists, what it tests, how we know the right answer and why it matters, and ends with a box "Test details". The [overview](index.md) lists every case of the tab.

## Top-of-atmosphere flux

**Why this check exists.** The flux at the top of the atmosphere passes through a file reader, a dilution factor and an interpolation. A missing square in the dilution, a wrong constant or a shifted wavelength grid gives the wrong stellar flux without an error message.

**What it tests.** The test runs the default model on 20 layers with photochemistry, in the HD 189733b setup: the spectrum `sflux-HD189_Moses11.txt`, $R_\star$ = 0.805 $R_\odot$ and $a$ = 0.03142 au, which give a dilution of 0.0142. The flux in the top interface of the run must equal the input spectrum times the dilution in every wavelength bin. This case is an input-to-output identity of the top boundary on a real run. It is not a test of the radiative transfer below the top.

**How we know the right answer.** The test reads the stellar file again with numpy and interpolates it with `np.interp`, independently of the `scipy.interpolate.interp1d` call in VULCAN. It computes the geometric factor from the configuration and the constants of `vulcan.phy_const`. Without the square, the factor $R_\star / a$ is 8.4 times larger, so a missing square fails the test in every bin.

**Why it matters.** The top-of-atmosphere flux drives every photolysis rate in the model. A wrong flux changes the photochemical products of the upper atmosphere.

![Top-of-atmosphere flux against wavelength and its relative difference from the diluted input spectrum](../assets/validation/toa_flux_diluted_stellar_spectrum.svg#only-light){ width="80%" }
![Top-of-atmosphere flux against wavelength and its relative difference from the diluted input spectrum](../assets/validation/toa_flux_diluted_stellar_spectrum_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The top panel shows the flux against wavelength on a log scale: the line is the diluted input spectrum (every 25th of the 2588 bins from 2 to 654 nm), and the dots are the flux in the top interface of the run (every 100th bin). The bottom panel shows their relative difference at every 50th bin, with the tolerance as a dashed line; bins that agree exactly are drawn at 1e-18. A pass means that every point lies below the dashed line. All bins agree to machine precision, well inside the tolerance.

??? info "Test details"
    - Reference: the stellar input file times (R_star / a)^2, interpolated with `np.interp`
    - Tolerance: `rtol=1e-12` on every bin
    - Test id: `tests/test_integration_quick.py::test_top_of_atmosphere_flux_is_diluted_stellar_spectrum`
    - Command: `make` in `fastchem_vulcan/`, then `pytest tests/test_integration_quick.py::test_top_of_atmosphere_flux_is_diluted_stellar_spectrum` (add `--validation-plots` to write the plot pair)
