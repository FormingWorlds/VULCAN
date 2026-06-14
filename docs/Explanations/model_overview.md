# VULCAN model overview

VULCAN is an open-source, one-dimensional chemical kinetics code for planetary and
exoplanetary atmospheres. It solves the coupled set of vertical continuity equations
for every chemical species, accounting for thermochemical reactions, photochemistry,
vertical transport (advection, eddy diffusion, molecular and thermal diffusion),
condensation and particle settling. The code was introduced in Tsai et al. (2017) [^tsai2017]
with a reduced C–H–O thermochemical network and extended in Tsai et al. (2021) [^tsai2021]
to C–H–N–O–S networks with photochemistry, condensation, advection, flexible boundary
conditions and temperature-dependent UV cross sections.

A reference for configuration parameters can be found [here](../Reference/config.md).

## Governing equation

VULCAN integrates the Eulerian continuity equation for the number density $n_i$
(cm$^{-3}$) of each species $i$:

$$\frac{\partial n_\text{i}}{\partial t} = P_\text{i} - L_\text{i} - \frac{\partial \Phi_\text{i}}{\partial z} \tag{1}$$

where $P_\text{i}$ and $L_\text{i}$ are the chemical production and loss rates (cm$^{-3}$ s$^{-1}$)
from both thermochemical and photochemical reactions, and $\Phi_i$ is the vertical
transport flux. The transport flux combines advection, eddy diffusion, and molecular
plus thermal diffusion; its form and discretization are described in
[basic equations and numerics](basic_equations.md).

The system is integrated forward in time with a second-order Rosenbrock solver
until a steady state is reached, rather than being solved as a boundary-value problem.
This makes the same machinery usable from hot, thermochemistry-dominated deep
atmospheres down to cold, photochemistry-dominated stratospheres.

## Scope and validation

The released model has been validated against three reference cases spanning a wide
range of temperature and oxidation state [^tsai2021]: the hot Jupiter HD 189733b
(intercompared with Moses et al. 2011 and Venot et al. 2012), Jupiter (low-temperature
hydrocarbon and condensation chemistry), and present-day Earth (an oxidizing atmosphere
with surface sources and sinks). The default configuration shipped in `config.py`
reproduces the HD 189733b setup ($R_\star = 0.805\,R_\odot$, $R_p = 1.138\,R_\mathrm{J}$,
$a = 0.03142$ au, $g = 2140$ cm s$^{-2}$).

## Code structure

The physics described in the Methods pages maps onto the source code, located in [src/vulcan](https://github.com/FormingWorlds/VULCAN/tree/main/src/vulcan), as follows.

| Module | Role |
|---|---|
| `vulcan.py` | Top-level driver (`run_vulcan`): builds the atmosphere, reads rates, runs the integration loop, saves output |
| `config.py` | The `Config` class holding every run parameter (replaces the legacy `vulcan_cfg.py`) |
| `store.py` | `Variables`, `AtmData`, `Parameters` containers for state, atmospheric structure, and numerical counters |
| `build_atm.py` | `InitialAbun` (initial mixing ratios, FastChem equilibrium, cold-trap setup) and `Atm` (T–P, $K_{zz}$, scale heights, molecular diffusion $D_{zz}$, settling, saturation pressures) |
| `make_chem_funs.py` | Parses the network file and **generates** `chem_funs.py` (the `chemdf` source term, the `Gibbs` reversal function, and the analytic Jacobians) |
| `op.py` | `ReadRate` (rate constants and cross sections), `Integration` (time stepping, condensation), `ODESolver`/`Ros2` (transport discretization, radiative transfer, photolysis, banded solve), `Output` (plots and saving) |
| `agni.py` | Optional coupling to the AGNI radiative–convective solver for self-consistent $T(p)$ (post-2021 addition; off by default) |
| `phy_const.py`, `paths.py`, `logs.py` | Physical constants (cgs), path resolution, logging |

## Model pages

1. [Basic equations and numerics](basic_equations.md)
2. [Boundary conditions](boundary_conditions.md)
3. [Chemical networks](chemical_networks.md)
4. [Computing photochemistry](photochemistry.md)
5. [Temperature-dependent UV cross sections](uv_cross_sections.md)
6. [Condensation and rainout](condensation.md)
7. [Chemistry of Ti and V compounds](titanium_vanadium.md)
8. [Photochemical haze precursors](haze_precursors.md)

[^tsai2017]: Tsai, S.-M., Lyons, J. R., Grosheintz, L., et al. (2017). VULCAN: An open-source, validated chemical kinetics Python code for exoplanetary atmospheres. *The Astrophysical Journal Supplement Series, 228*(2), 20. https://doi.org/10.3847/1538-4365/228/2/20

[^tsai2021]: Tsai, S.-M., Malik, M., Kitzmann, D., et al. (2021). A comparative study of atmospheric chemistry with VULCAN. *The Astrophysical Journal, 923*(2), 264. https://doi.org/10.3847/1538-4357/ac29bc