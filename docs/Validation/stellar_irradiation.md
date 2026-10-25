# Stellar irradiation

## What this page checks

This page checks the starlight that enters the top interface, the upper boundary of the model atmosphere. VULCAN reads the stellar spectrum from the file `sflux_file`, which gives the flux at the surface of the star in erg s$^{-1}$ cm$^{-2}$ nm$^{-1}$. It dilutes this flux by the geometric factor $(R_\star / a)^2$, with $R_\star$ the stellar radius and $a$ the orbital distance. It then interpolates the result onto the wavelength bins of the model (`build_atm.read_sflux`). Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. If this part of the code is wrong, every photolysis rate would be wrong, and so would the chemistry of the upper atmosphere; the case below has the 4 parts of each case and a box "Test details" ([overview](index.md)).

## Top-of-atmosphere flux

**Why this check exists.** A missing square in the dilution, or a unit slip in the stellar radius or the orbital distance, gives the wrong stellar flux without an error message.

**What it tests.** The test runs the default model on 20 layers with photochemistry. It uses the HD 189733b setup: the spectrum `sflux-HD189_Moses11.txt`, $R_\star$ = 0.805 $R_\odot$ and $a$ = 0.03142 au (astronomical units), which give a dilution of 0.0142. The flux in the top interface of the run must equal the input spectrum times the dilution, with a relative difference below 1e-12 in every wavelength bin. This case is an input-to-output identity of the top boundary on a real run. It is not a test of the radiative transfer below the top.

**How we know the right answer.** The test reads the stellar file again with numpy and interpolates it with `np.interp`, independently of the `scipy.interpolate.interp1d` call in VULCAN. It computes the geometric factor from the configuration and the constants of `vulcan.phy_const`. It also pins the dilution to 0.0141962 (relative 1e-5), the value from the IAU 2015 solar radius and the IAU 2012 au, so a wrong constant in `vulcan.phy_const` fails the test. Without the square, the factor $R_\star / a$ is 8.4 times larger. The test asserts that this wrong flux differs from the right one by more than the tolerance in every bin, so the tolerance can tell the two apart.

**Why it matters.** VULCAN computes the flux in the atmosphere from the flux in the top interface, and the photolysis rates from that flux. A wrong top flux would change the photochemical products of the upper atmosphere.

![Top-of-atmosphere flux against wavelength and its relative difference from the diluted input spectrum](../assets/validation/toa_flux_diluted_stellar_spectrum.svg#only-light){ width="80%" }
![Top-of-atmosphere flux against wavelength and its relative difference from the diluted input spectrum](../assets/validation/toa_flux_diluted_stellar_spectrum_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The top panel shows the flux against wavelength on a log scale: the line is the diluted input spectrum and the dots are the flux in the top interface of the run. The bottom panel shows their relative difference, and a pass means that every point lies below the dashed tolerance line. The line shows every 25th of the 2588 bins from 2 to 654 nm, the dots every 100th bin and the bottom panel every 50th bin; bins that agree exactly are drawn at 1e-18. All bins agree to machine precision, well inside the tolerance.

??? info "Test details"
    - Reference: the stellar input file times (R_star / a)^2, interpolated with `np.interp`
    - Tolerance: `rtol=1e-12` on every bin
    - Test id: `tests/test_integration_quick.py::test_top_of_atmosphere_flux_is_diluted_stellar_spectrum`
    - Command: `make` in `fastchem_vulcan/`, then `pytest tests/test_integration_quick.py::test_top_of_atmosphere_flux_is_diluted_stellar_spectrum` (add `--validation-plots` to write the plot pair)
