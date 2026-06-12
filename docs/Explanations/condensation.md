# Condensation and rainout

VULCAN treats condensation and evaporation through a particle growth rate, and lets the
resulting particles fall under gravity (Section 2.6 of Tsai et al. 2021 [^tsai2021]).

## Growth rate

For the schematic condensation/evaporation reaction
$A_\mathrm{(gas)} \leftrightarrow A_\mathrm{(particle)}$, the growth rate in the
**continuum regime** (particles larger than the mean free path, Knudsen number $\mathrm{Kn} < 1$)
follows the mass-balance expression (Seinfeld & Pandis 2016 [^sp2016]):

$$\frac{\mathrm{d}n_A}{\mathrm{d}t} = -\frac{D_A\,m_A}{\rho_\text{p}\,r_\text{p}^{2}}\left(n_A - n_A^\mathrm{sat}\right) n_A \tag{14}$$

where $D_A$ and $m_A$ are the molecular diffusion coefficient and mass of gas $A$, $\rho_\text{p}$ and
$r_\text{p}$ the particle density and radius, and $n_A^\mathrm{sat}$ the saturation number density.
A negative value (when $n_A > n_A^\mathrm{sat}$) is condensation; a positive value
(when $n_A < n_A^\mathrm{sat}$) is evaporation.

!!! note "Difference from earlier models"
    Hu et al. (2012) [^hu2012] and Rimmer & Helling (2016) [^rh2016] use the *kinetic*-regime
    growth rate. VULCAN uses the **continuum** form because, for most applications,
    condensation occurs in the lower atmosphere with micron-size or larger particles. The
    code does contain a commented-out kinetic-regime expression (`rate_c`) for reference.

## Particle settling

Once gas condenses to particles, they fall at the terminal velocity from Stokes' law:

$$v_\text{s} = \frac{2\,r_\text{p}^{2}\,\rho_\text{p}\,g}{9\,\mu} \tag{15}$$

with $\mu$ the atmospheric dynamic viscosity. The slip-correction factor is taken as unity
(large-particle limit).

## Saturation vapor pressures

[`build_atm.Atm.sp_sat`](../Reference/api/build_atm.md#vulcan.build_atm.Atm.sp_sat) precomputes the saturation vapor pressure (dyne cm$^{-2}$) for each
condensable species:

| Species | Parameterization |
|---|---|
| H$_2$O | Ackerman & Marley (2001), separate ice ($T<0$ °C) and liquid branches |
| NH$_3$ | Weast (1971) |
| H$_2$SO$_4$ | Ayers form |
| S$_2$, S$_8$ | Zahnle et al. (2016), refit from Lyons (2008) [^lyons2008]; piecewise at 413 K |
| S$_4$ | Lyons (2008) |
| C | NIST fit |
| H$_2$S | Giauque & Blue (1936), separate ice/liquid branches |

The dynamic viscosity used in Eq. (15) is taken per background gas (`atm_base`) from the
Cloutman (2000) combustion-coefficient table in [`build_atm.Atm.f_mu_dz`](../Reference/api/build_atm.md#vulcan.build_atm.Atm.f_mu_dz).

## Implementation in VULCAN

Condensation is driven from [`op.Integration`](../Reference/api/op.md#vulcan.op.Integration) and applied within the transport operator:

- `Integration.conden` updates the condensation/evaporation rate coefficients each step
  using the continuum form of Eq. (14), one branch per species
  (`H2O → H2O_l_s`, `NH3 → NH3_l`, `H2SO4 → H2SO4_l`, `S2/S4/S8 → *_l_s`, `C → C_s`).
- Settling enters the transport operator and Jacobian through the
  `diffdf_settling` / `lhs_jac_settling` variants when `use_settling = True`; the settling
  velocity is stored in `atm.vs` (computed in `f_mu_dz`).
- An alternative **relaxation** scheme (`h2o_conden_evap_relax`, `nh3_conden_evap_relax`) is
  used when `use_relax` is set, removing super-saturated vapor with an implicit-Euler step.
- Because condensation operates on a short timescale, the code can fix the condensing species
  and condensates after dynamic equilibrium is reached (`fix_species`,
  `fix_species_from_coldtrap_lev`) — a quasi-steady-state approach that decouples fast and
  slow processes. The cold-trap level is found and species are held fixed above it.

Non-gaseous (condensate) species are excluded when forming mixing ratios from number
densities (`atm.gas_indx`), so hydrostatic balance is maintained over the gas phase only.

### Relevant parameters

| Parameter | Meaning |
|---|---|
| `use_condense` | enable condensation chemistry |
| `use_settling` | enable gravitational settling |
| `condense_sp` | species allowed to condense |
| `non_gas_sp` | condensate (non-gaseous) species |
| `r_p`, `rho_p` | per-condensate particle radius (cm) and density (g cm$^{-3}$) |
| `humidity` | relative-humidity multiplier applied to the H$_2$O saturation |
| `fix_species`, `start_conden_time`, `stop_conden_time` | quasi-steady-state fixing controls |

[^tsai2021]: Tsai, S.-M., Malik, M., Kitzmann, D., et al. (2021). A comparative study of atmospheric chemistry with VULCAN. *The Astrophysical Journal, 923*(2), 264. https://doi.org/10.3847/1538-4357/ac29bc

[^sp2016]: Seinfeld, J. H., & Pandis, S. N. (2016). *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 3rd ed. Wiley.

[^hu2012]: Hu, R., Seager, S., & Bains, W. (2012). Photochemistry in terrestrial exoplanet atmospheres. I. *The Astrophysical Journal, 761*(2), 166. https://doi.org/10.1088/0004-637X/761/2/166

[^rh2016]: Rimmer, P. B., & Helling, C. (2016). A chemical kinetics network for lightning and life in planetary atmospheres. *The Astrophysical Journal Supplement Series, 224*(1), 9. https://doi.org/10.3847/0067-0049/224/1/9

[^lyons2008]: Lyons, J. R. (2008). An estimate of the equilibrium speciation of sulfur vapor over solid sulfur and implications for planetary atmospheres. *Journal of Sulfur Chemistry, 29*(3–4), 269–279. https://doi.org/10.1080/17415990802208743