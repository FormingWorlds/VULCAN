# Computing photochemistry

Stellar radiation drives disequilibrium chemistry through photodissociation. VULCAN computes
the wavelength-dependent radiation field, then integrates it against absorption cross sections
to obtain photolysis rates (Section 2.4 of Tsai et al. 2021 [^tsai2021]).

## Photodissociation

A photolysis reaction is written schematically as the unimolecular process

$$A + h\nu \rightarrow B + C$$

producing reactive radicals that initiate chains essential to atmospheric chemistry (the
ozone cycle on Earth, organic-haze formation on Titan).

## Actinic flux

The radiation field is described by the **actinic flux** $J(z,\lambda)$: photons per unit
time, area, and wavelength from all directions. It has a direct-beam and a diffuse component:

$$J(z,\lambda) = J_\infty(\lambda)\,e^{-\tau(z,\lambda)/\mu} + J_\mathrm{diff}(z,\lambda) \tag{8}$$

with $\mu = \cos\theta$ for stellar zenith angle $\theta$. The direct term has no cosine
prefactor (unlike radiative heating) because the number of intercepted molecules is
independent of beam direction.

The optical depth includes absorption and scattering,

$$\tau(z,\lambda) = \int \sum_i \big[\sigma_{a,i}(\lambda) + \sigma_{s,i}(\lambda)\big]\,n_i\,\mathrm{d}z \tag{9}$$

where the absorption cross section $\sigma_{a,i}$ can differ from the photodissociation cross
section (absorption is not always followed by dissociation).

The diffuse flux is obtained with the two-stream approximation of Malik et al. (2019) [^malik2019]
and converted to total intensity using the first Eddington coefficient $\bar\epsilon$
(Heng et al. 2018 [^heng2018]):

$$J_\mathrm{diff}(z,\lambda) = \frac{F_\mathrm{diff}(z,\lambda)}{\bar\epsilon}, \qquad \bar\epsilon = 0.5 \tag{10}$$

Multiple scattering is handled iteratively; the converged state is typically reached within
~200 iterations for a strongly irradiated hot Jupiter.

## Photolysis rates

The photolysis rate coefficient integrates the actinic flux and absorption cross section over
wavelength,

$$k = \int \sigma_a(\lambda)\,q(\lambda)\,J(z,\lambda)\,\mathrm{d}\lambda \tag{11}$$

with $q(\lambda)$ the quantum yield (probability of a given branch per absorbed photon), and
the photolysis rate of the reaction is $\mathrm{d}n_A/\mathrm{d}t = -k\,n_A$. Cross sections
are taken from the Leiden Observatory database (Heays et al. 2017 [^heays2017]) where available.

## Implementation in VULCAN

The radiative-transfer and photolysis routines live in [`op.ODESolver`](../Reference/api/op.md#vulcan.op.ODESolver):

| Step | Routine | Notes |
|---|---|---|
| Optical depth, Eq. (9) | `compute_tau` | accumulates `y·dz·cross` (and scattering `cross_scat`) downward; uses `cross_T` for T-dependent species |
| Radiation field, Eqs. (8),(10) | `compute_flux` | two-stream with $\zeta_\pm$, transmission `tran`, Eddington `edd`; direct beam by Beer's law `sflux_top·exp(-tau/cos(sl_angle))` |
| Photolysis rates, Eq. (11) | `compute_J` | trapezoidal sum over the two-resolution bin grid; writes `k[pho_rate_index]` |
| Photoionization rates | `compute_Jion` | only when `use_ion = True` |

The stellar spectrum is read and scaled to the planet in [`build_atm.Atm.read_sflux`](../Reference/api/build_atm.md#vulcan.build_atm.Atm.read_sflux), by the
factor $(R_\star R_\odot / a\,\mathrm{au})^2$, and interpolated onto the bin grid with an
energy-conservation check logged as a percentage. The wavelength grid is built in
[`ReadRate.make_bins_read_cross`](../Reference/api/op.md#vulcan.op.ReadRate.make_bins_read_cross) using two uniform resolutions, `dbin1` below `dbin_12trans`
and `dbin2` above, so that the fine structure in the (X)UV is resolved without oversampling
the longer wavelengths (resolution errors are quantified in Appendix B of the paper).

Each computed branch rate is multiplied by `f_diurnal`, the diurnal-averaging factor
(`1.0` for a tidally locked planet, `0.5` for a rotating planet such as Earth). The actinic
flux is recomputed periodically rather than every step: every `ini_update_photo_frq` steps
initially, switching to `final_update_photo_frq` once the solution is near convergence.

### Relevant parameters

| Parameter | Meaning |
|---|---|
| `sl_angle` | stellar zenith angle $\theta$ (radians); dayside/terminator-average values discussed in Appendix C of the paper ($\approx 58°$ and $67°$) |
| `edd` | first Eddington coefficient $\bar\epsilon$ (0.5) |
| `f_diurnal` | diurnal-averaging factor |
| `dbin1`, `dbin2`, `dbin_12trans` | bin widths (nm) and the switch wavelength |
| `scat_sp` | species included in Rayleigh scattering (e.g. H$_2$, He) |
| `r_star`, `orbit_radius` | stellar radius and orbital distance for flux scaling |

[^tsai2021]: Tsai, S.-M., Malik, M., Kitzmann, D., et al. (2021). A comparative study of atmospheric chemistry with VULCAN. *The Astrophysical Journal, 923*(2), 264. https://doi.org/10.3847/1538-4357/ac29bc

[^malik2019]: Malik, M., Kitzmann, D., Mendonça, J. M., et al. (2019). Self-luminous and irradiated exoplanetary atmospheres explored with HELIOS. *The Astronomical Journal, 157*(5), 170. https://doi.org/10.3847/1538-3881/ab1084

[^heng2018]: Heng, K., Malik, M., & Kitzmann, D. (2018). Analytical models of exoplanetary atmospheres. VI. VI. Full Solutions for Improved Two-stream Radiative Transfer, Including Direct Stellar Beam. *The Astrophysical Journal Supplement Series, 237*(2), 29. https://doi.org/10.3847/1538-4365/aad199

[^heays2017]: Heays, A. N., Bosman, A. D., & van Dishoeck, E. F. (2017). Photodissociation and photoionisation of atoms and molecules of astrophysical interest. *Astronomy & Astrophysics, 602*, A105. https://doi.org/10.1051/0004-6361/201628742