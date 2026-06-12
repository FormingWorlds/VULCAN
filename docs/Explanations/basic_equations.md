# Basic equations and numerics

VULCAN solves the one-dimensional Eulerian continuity equation (Section 2.1 of
Tsai et al. 2021 [^tsai2021]) for the number density of each species. This page covers
the transport flux, its spatial discretization, and the time-integration scheme.

## Continuity equation

For each species $\text{i}$,

$$\frac{\partial n_{\text{i}}}{\partial t} = P_{\text{i}} - L_{\text{i}} - \frac{\partial \Phi_{\text{i}}}{\partial z} \tag{1}$$

with $n_{\text{i}}$ the number density (cm$^{-3}$), $P_{\text{i}}$ and $L_{\text{i}}$ the chemical production and
loss rates (cm$^{-3}$ s$^{-1}$), and $\Phi_{\text{i}}$ the vertical transport flux. $P_{\text{i}} - L_{\text{i}}$
is evaluated by `chemdf` in the generated module `chem_funs.py`; see
[chemical networks](chemical_networks.md).

## Transport flux

Assuming hydrostatic balance, the flux combines advection, eddy diffusion, and molecular
plus thermal diffusion (following Chamberlain & Hunten 1987 [^ch1987]):

$$\Phi_{\text{i}} = n_{\text{i}} v
\;-\; K_{zz}\,n_\mathrm{tot}\,\frac{\partial X_{\text{i}}}{\partial z}
\;-\; D_{\text{i}}\left[\frac{\partial n_{\text{i}}}{\partial z}
+ n_{\text{i}}\left(\frac{1}{H_{\text{i}}} + \frac{1+\alpha_T}{T}\frac{\mathrm{d}T}{\mathrm{d}z}\right)\right] \tag{2}$$

where $v$ is the vertical wind, $K_{zz}$ the eddy diffusion coefficient, $D_{\text{i}}$ the molecular
diffusion coefficient, $X_{\text{i}} = n_{\text{i}}/n_\mathrm{tot}$ the mixing ratio, $H_{\text{i}} = k_B T/(m_{\text{i}} g)$ the
species scale height, and $\alpha_T$ the thermal diffusion factor. The four contributions are,
in order: advection along the wind, eddy diffusion smearing out compositional gradients,
molecular diffusion driving each species toward its own diffusive equilibrium, and thermal
diffusion (a secondary effect except for light species in the thermosphere).

!!! info "Two forms of Equation 2"
    Eq. 2 given above is the version as published in Tsai et al. (2021) [^tsai2021]. In the code, the term in round brackets is written as $-1/H_0 + 1/H_{\text{i}} + (\alpha_T/T)\,\mathrm{d}T/\mathrm{d}z$, where $H_0 = k_B T/(\bar m g)$ is the
    atmospheric (mean) scale height. The two forms are related by the hydrostatic identity
    $\,n_\mathrm{tot}^{-1}\,\partial n_\mathrm{tot}/\partial z = -1/H_0 - T^{-1}\,\mathrm{d}T/\mathrm{d}z\,$.

## Spatial discretization

The diffusive terms use a **second-order central difference** and the advective term a
**first-order upwind** scheme (Brasseur & Jacob 2017 [^bj2017]). On the staggered grid the
flux divergence in layer $\text{j}$ is

$$\frac{\partial \Phi_{\text{i}}}{\partial z}\bigg|_{\text{j}} \approx \frac{\Phi_{\text{i},\,\text{j}+1/2} - \Phi_{\text{i},\,\text{j}-1/2}}{\Delta z_{\text{j}}}$$

Substituting Eq. (2) reduces Eq. (1) to a system of ODEs of the form
$A_{\text{j}} y_{\text{j}} + B_{\text{j}+1} y_{\text{j}+1} + C_{\text{j}-1} y_{\text{j}-1}$ per layer.

### Implementation in VULCAN

The transport operator is assembled in [`op.ODESolver`](../Reference/api/op.md#vulcan.op.ODESolver), with four variants selected at
runtime in [`Ros2.solver`](../Reference/api/op.md#vulcan.op.Ros2.solver) according to the [config flags](../Reference/config.md):

| Method | `use_moldiff` | `use_settling` | `use_vm_mol` |
|---|---|---|---|
| `diffdf_no_mol` | `False` | `False` | `False` |
| `diffdf` | `True` | `False` | `False` |
| `diffdf_settling` | `True` | `True` | `False` |
| `diffdf_settling_vm` | `True` | `True` | `True` |

The upwind advection appears as the `(vz>0)`/`(vz<0)` masks. Interface quantities ($T_{\text{i}}$, $H_{p,\text{i}}$, $\Delta z_{\text{i}}$) are
precomputed in [`build_atm.Atm.f_mu_dz`](../Reference/api/build_atm.md), and the molecular
diffusion coefficients $D_{zz}$ and thermal diffusion factors $\alpha_T$ in
[`build_atm.Atm.mol_diff`](../Reference/api/build_atm.md). The latter scales a reference binary
diffusion coefficient (per background gas `atm_base`: H$_2$, N$_2$, O$_2$, or CO$_2$, after
Marrero & Mason 1972 [^mm1972] / Banks & Kockarts 1973 [^bk1973]) by molecular mass,
corresponding to Appendix A of the paper [^tsai2021].

## Time integration

The stiff ODE system is integrated with a **second-order Rosenbrock** method [^verwer1999] in
[`op.Ros2.solver`](../Reference/api/op.md#vulcan.op.Ros2.solver). For $\mathrm{d}n/\mathrm{d}t = f(n)$,
explicit (forward-Euler) differencing is cheap but unstable at the large steps stiff chemistry
needs, while fully implicit (backward-Euler) differencing is stable but requires an expensive
Newton–Raphson iteration for $n_{\text{k}+1}$. The Rosenbrock method is *linearly implicit*: it uses
the Jacobian $J \equiv \partial f/\partial n$ directly, avoiding iteration. The two-stage,
second-order scheme is

$$n_{\text{k}+1} = n_{\text{k}} + \tfrac{3}{2}\Delta t\, g_1 + \tfrac{1}{2}\Delta t\, g_2$$

$$(I - \gamma \Delta t J)\,g_1 = f(n_{\text{k}}), \qquad
(I - \gamma \Delta t J)\,g_2 = f(n_{\text{k}} + \Delta t\, g_1) - 2 g_1$$

with $\gamma = 1 + 1/\sqrt{2}$. Verwer et al. (1999) [^verwer1999] recommend this method for
its stability over large stepsizes, which suits the needs of chemical kinetics. Both stages
share the same left-hand-side matrix

$$\mathrm{LHS} = \frac{1}{\gamma\,\Delta t}\,I - J$$

so it is factored only once per step. The chemical part of $J$ is **analytic** (generated
symbolically by [`make_chem_funs.make_jac`/`make_neg_jac`](../Reference/api/make_chem_funs.md#vulcan.make_chem_funs.make_jac)
into `chem_funs.symjac`/`neg_symjac`); the transport part is added in the `lhs_jac_*` methods.
The combined block-tridiagonal matrix is stored in banded form (`store_bandM`) and solved with
`scipy.linalg.solve_banded`.

## Step-size control and convergence
 
The step size adapts to the truncation error $\mathcal{E} = |n_{\text{k}+1} - n^*_{\text{k}+1}|$, where
$n^*_{\text{k}+1} = n_{\text{k}} + \Delta t\, g_1$ is the embedded first-order estimate. In
[`Ros2.step_size`](../Reference/api/op.md#vulcan.op.Ros2.step_size),
 
$$\Delta t_{\text{k}+1} = \Delta t_{\text{k}} \cdot \mathrm{clamp}\!\left(0.9\,(\mathrm{rtol}/\mathcal{E})^{1/2},\ r_{\min},\ r_{\max}\right)$$
 
where `rtol` is the desired relative tolerance and 0.9 is a safety factor. The step-change
factor is bounded to $[r_{\min}, r_{\max}] = [0.5,\ 2]$ — the config parameters `dt_var_min`
and `dt_var_max` — and the step itself to $[10^{-8},\ 10^{18}]$ s (`dt_min`, `dt_max`). A step
with $\mathcal{E} > \mathrm{rtol}$ is rejected and retried with a smaller stepsize.
 
Steady state is declared in [`op.Integration`](../Reference/api/op.md#vulcan.op.Integration.conv)
from the long-term relative variation, following the criteria of Tsai et al. (2017) [^tsai2017],
evaluated over a look-back window from $f\tau$ to the current integration time $\tau$ (with $f$
set by `st_factor`):

$$\Delta\hat{n} \equiv \max_{\text{i},\text{j}} \frac{|n_{\text{i},\text{j},\text{k}} - n_{\text{i},\text{j},\text{k}'}|}{n_{\text{i},\text{j},\text{k}}},
\qquad \Delta t \equiv t_{\text{k}} - t_{\text{k}'}$$
 
where $\text{k}'$ is the timestep at $f\tau$. Species with mixing ratios below `mtol_conv` or number
densities below `atol` are excluded from the maximum so that negligible trace species do not
trigger it spuriously. Convergence requires both $\Delta\hat{n}$ and its rate
$\Delta\hat{n}/\Delta t$ to fall below either the primary thresholds (`yconv_cri`, `slope_cri`)
or the relaxed thresholds (`yconv_min`, `slope_min`), where the relaxed slope adapts to the
local vertical-diffusion timescale, $\sim \min_{\text{j}} K_{zz,\text{j}}/(0.1\,H_{p,\text{j}})^2$. With photochemistry
enabled, both conditions additionally require the actinic-flux change to fall below `flux_cri`.


[^tsai2021]: Tsai, S.-M., Malik, M., Kitzmann, D., et al. (2021). A comparative study of atmospheric chemistry with VULCAN. *The Astrophysical Journal, 923*(2), 264. https://doi.org/10.3847/1538-4357/ac29bc

[^tsai2017]: Tsai, S.-M., Lyons, J. R., Grosheintz, L., et al. (2017). VULCAN: An open-source, validated chemical kinetics Python code for exoplanetary atmospheres. *The Astrophysical Journal Supplement Series, 228*(2), 20. https://doi.org/10.3847/1538-4365/228/2/20

[^ch1987]: Chamberlain, J. W., & Hunten, D. M. (1987). *Theory of Planetary Atmospheres*, 2nd ed. Academic Press.

[^bj2017]: Brasseur, G. P., & Jacob, D. J. (2017). *Modeling of Atmospheric Chemistry*. Cambridge University Press.

[^verwer1999]: Verwer, J. G., Spee, E. J., Blom, J. G., & Hundsdorfer, W. (1999). A second-order Rosenbrock method applied to photochemical dispersion problems. *SIAM Journal on Scientific Computing, 20*(4), 1456–1480. https://doi.org/10.1137/S1064827597326651

[^mm1972]: Marrero, T. R., & Mason, E. A. (1972). Gaseous diffusion coefficients. *Journal of Physical and Chemical Reference Data, 1*(1), 3–118. https://doi.org/10.1063/1.3253094

[^bk1973]: Banks, P. M., & Kockarts, G. (1973). *Aeronomy*. Academic Press.