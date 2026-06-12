# Boundary conditions

The solution of the ODE system (Eq. 1 in [basic equations](basic_equations.md)) requires
boundary conditions at the top and bottom of the column. VULCAN supports the three standard
quantities used in photochemical models (Hu et al. 2012 [^hu2012]): **flux**, **velocity**,
and **mixing ratio**, plus a diffusion-limited escape option (Section 2.2 of
Tsai et al. 2021 [^tsai2021]).

## Flux boundary conditions

Because the discretized flux at an interface depends on the layers above and below, the
fluxes at the very top and bottom are otherwise unspecified and must be prescribed.
Constant fluxes represent, for example, surface emission for rocky planets at the bottom,
or meteoritic infall / escape outflow at the top.

Without additional constraints the flux is set to **zero**, i.e. no net material exchange.
Zero flux is the default lower boundary condition for gas giants, where the bottom is placed
deep enough to lie in the thermochemical-equilibrium region.

## Diffusion-limited escape

For light escaping species, the top flux can instead be set to the diffusion-limited value

$$\Phi_{\text{i},\mathrm{top}} = -D_\text{i}\left(\frac{1}{H_\text{i}} - \frac{1}{H_0}\right) n_\text{i} \tag{6}$$

which assumes escape is throttled by diffusion into the exosphere.

## Velocity boundary conditions

Sources and sinks that scale with abundance are represented by a velocity. A (dry/wet)
**deposition velocity** $v_\mathrm{dep}$ parameterizes surface uptake at the bottom; an
upward velocity at the top can represent escape or other outflow.

## Fixed mixing ratio

When the detailed exchange is complex but the abundance is known, a constant mixing ratio
can be prescribed (e.g. surface water set by relative humidity on an ocean world). Because a
fixed mixing ratio does not allow the composition at the boundary to change, it cannot be
combined with flux or velocity conditions for the same species.

## Implementation in VULCAN

Boundary conditions are read in [`build_atm.Atm.BC_flux`](../Reference/api/build_atm.md#vulcan.build_atm.Atm.BC_flux) and applied inside the transport
operator and Jacobian:

| Condition | Config flag | Input file column(s) | Applied in |
|---|---|---|---|
| Constant top flux | `use_topflux` | `top_BC_flux_file`: species, flux | `diff[-1] += top_flux/dzi[-1]` |
| Constant bottom flux | `use_botflux` | `bot_BC_flux_file`: species, flux, $v_\mathrm{dep}$ | `diff[0] += (bot_flux - y[0]·bot_vdep)/dzi[0]` |
| Deposition velocity | `use_botflux` | `bot_BC_flux_file` (3rd column) | `lhs_jac_*` diagonal term `-bot_vdep/dzi[0]` |
| Diffusion-limited escape | `diff_esc` (e.g. `['H']`) | — | `op.Integration.update_phi_esc` |
| Fixed bottom mixing ratio | `use_fix_sp_bot` | `bot_BC_flux_file` (4th column) | `op.Ros2.solver` clamps `sol[0]`; `solver_fix_all_bot` |

The diffusion-limited flux in `update_phi_esc` is computed as
`-Dzz[-1]·y[-1]·(1/Hp[-1] - ms·g[-1]/(Navo·kb·T[-1]))`, i.e. Eq. (6), and its magnitude is
capped at `max_flux` to keep the escape rate bounded.

A constant-flux term contributes nothing to the Jacobian, but the deposition term **does**
and is therefore included on the diagonal of the left-hand-side matrix in every `lhs_jac_*`
variant. When all bottom-layer species are held fixed (`solver_fix_all_bot`,
`lhs_jac_*_fix_all_bot`), the entire bottom column of the Jacobian is zeroed and the layer
is reset to its initial (equilibrium) composition each step.


[^tsai2021]: Tsai, S.-M., Malik, M., Kitzmann, D., et al. (2021). A comparative study of atmospheric chemistry with VULCAN. *The Astrophysical Journal, 923*(2), 264. https://doi.org/10.3847/1538-4357/ac29bc

[^hu2012]: Hu, R., Seager, S., & Bains, W. (2012). Photochemistry in terrestrial exoplanet atmospheres. I. Photochemistry model and benchmark cases. *The Astrophysical Journal, 761*(2), 166. https://doi.org/10.1088/0004-637X/761/2/166