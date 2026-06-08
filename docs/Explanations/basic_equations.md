# Basic equations and numerics

VULCAN solves the one-dimensional Eulerian continuity equation (Section 2.1 of
Tsai et al. 2021 [^tsai2021]) for the number density of each species. This page covers
the transport flux, its spatial discretization, and the time-integration scheme.

## Continuity equation

For each species $i$,

$$\frac{\partial n_i}{\partial t} = P_i - L_i - \frac{\partial \Phi_i}{\partial z} \tag{1}$$

with $n_i$ the number density (cm$^{-3}$), $P_i$ and $L_i$ the chemical production and
loss rates (cm$^{-3}$ s$^{-1}$), and $\Phi_i$ the vertical transport flux. $P_i - L_i$
is evaluated by `chemdf` in the generated module `chem_funs.py`; see
[chemical networks](chemical_networks.md).

## Transport flux

Assuming hydrostatic balance, the flux combines advection, eddy diffusion, and molecular
plus thermal diffusion (following Chamberlain & Hunten 1987 [^ch1987]):

$$\Phi_i = n_i v
\;-\; K_{zz}\,n_\mathrm{tot}\,\frac{\partial X_i}{\partial z}
\;-\; D_i\,\left[\frac{\partial n_i}{\partial z}
+ n_i\left(\frac{1}{H_i} - \frac{1}{H_0} + \frac{\alpha_T}{T}\frac{\mathrm{d}T}{\mathrm{d}z}\right)\right] \tag{2}$$

where $v$ is the vertical wind, $K_{zz}$ the eddy diffusion coefficient, $D_i$ the molecular
diffusion coefficient, $X_i = n_i/n_\mathrm{tot}$ the mixing ratio, $H_i = k_B T/(m_i g)$ the
species scale height, $H_0 = k_B T/(\bar{m} g)$ the atmospheric (mean) scale height, and
$\alpha_T$ the thermal diffusion factor. The four contributions are, in order: advection
along the wind, eddy diffusion smearing out compositional gradients, molecular diffusion
driving each species toward its own diffusive equilibrium, and thermal diffusion (a
secondary effect except for light species in the thermosphere).

## Spatial discretization

The diffusive terms use a **second-order central difference** and the advective term a
**first-order upwind** scheme (Brasseur & Jacob 2017 [^bj2017]). On the staggered grid the
flux divergence in layer $j$ is

$$\frac{\partial \Phi_i}{\partial z}\bigg|_j \approx \frac{\Phi_{i,\,j+1/2} - \Phi_{i,\,j-1/2}}{\Delta z_j}$$

Substituting Eq. (2) reduces Eq. (1) to a system of ODEs of the form
$A_j y_j + B_{j+1} y_{j+1} + C_{j-1} y_{j-1}$ per layer.

### Implementation in VULCAN

The transport operator is assembled in [`op.ODESolver`](../Reference/api/op.md), with four variants selected at
runtime in `Ros2.solver` according to the [config flags](../Reference/config.md):

| Method | `use_moldiff` | `use_settling` | `use_vm_mol` |
|---|---|---|---|
| `diffdf_no_mol` | `False` | `False` | `False` |
| `diffdf` | `True` | `False` | `False` |
| `diffdf_settling` | `True` | `True` | `False` |
| `diffdf_settling_vm` | `True` | `True` | `True` |

In each, the molecular-diffusion advective velocity carries the term $-1/H_0 + 1/H_i + (\alpha_T/T)\,\mathrm{d}T/\mathrm{d}z$ from Eq. (2). The upwind advection
appears as the `(vz>0)`/`(vz<0)` masks. Interface quantities ($T_i$, $H_{p,i}$, $\Delta z_i$)
are precomputed in [`build_atm.Atm.f_mu_dz`](../Reference/api/build_atm.md), and the molecular diffusion coefficients $D_{zz}$
and thermal diffusion factors $\alpha_T$ in [`build_atm.Atm.mol_diff`](../Reference/api/build_atm.md). The latter scales a
reference binary diffusion coefficient (per background gas `atm_base`: H$_2$, N$_2$, O$_2$,
or CO$_2$, after Marrero & Mason 1972 [^mm1972] / Banks & Kockarts 1973 [^bk1973]) by molecular
mass, corresponding to Appendix A of the paper [^tsai2021].

## Time integration

The stiff ODE system is integrated with a **second-order Rosenbrock**
method [^verwer1997] implemented in
[`op.Ros2.solver`](../Reference/api/op.md). Consider the differential equation

$$\frac{dn}{dt} = f(n)$$

where $n$ is the dependent variable, $t$ is the independent variable, and $f(n)$ is an
arbitrary function of $n$. If we discretise this equation with stepsize $\Delta t$, then
explicit differencing (i.e., the forward Euler method) yields

$$n_{k+1} = n_k + \Delta t\, f(n_k)$$

where $k$ denotes the $k$th index of the discretised independent variable. Fully implicit
differencing (i.e., the backward Euler method) gives

$$n_{k+1} = n_k + \Delta t\, f(n_{k+1})$$

Solving these nonlinear equations for $n_{k+1}$ generally involves using a
Newton–Raphson-like iteration method, which is the most computationally expensive part in
chemical kinetics. A more efficient approach is to perform a Taylor expansion of $f$ for the
linear term (i.e., the semi-implicit Euler method [^press2007]),

$$n_{k+1} = n_k + \Delta t \left[ f(n_k) + \frac{\partial f}{\partial n}\bigg|_{n_k} (n_{k+1} - n_k) \right]$$

which may be expressed more compactly as

$$n_{k+1} = n_k + \Delta t \left(I - \Delta t J\right)^{-1} f(n_k)$$

where $J \equiv \partial f / \partial n \big|_{n_k}$ is the Jacobian matrix evaluated at
$n_k$. The $s$-stage Rosenbrock method generalises this. The second-order ($s = 2$) scheme takes the form

$$n_{k+1} = n_k + \tfrac{3}{2}\Delta t\, g_1 + \tfrac{1}{2}\Delta t\, g_2$$

$$\left(I - \gamma \Delta t J\right) g_1 = f(n_k)$$

$$\left(I - \gamma \Delta t J\right) g_2 = f(n_k + \Delta t\, k_1) - 2k_1$$

Verwer et al. (1998) [^verwer1997] recommended the second-order Rosenbrock method for being stable over
large stepsizes, which suits the needs of chemical kinetics. Each stage solves a linear
system with the same left-hand-side matrix

$$\mathrm{LHS} = \frac{1}{\gamma\,\Delta t}\,I - J$$

where $\gamma = 1 + 1/\sqrt{2}$.

The chemical part of the Jacobian is **analytic** (generated symbolically by
[`make_chem_funs.make_jac`/`make_neg_jac`](../Reference/api/make_chem_funs.md#vulcan.make_chem_funs.make_jac) into `chem_funs.symjac`/`neg_symjac`); the
transport part is added in the `lhs_jac_*` methods. The Jacobian has a block-tridiagonal
structure:

$$J_{\alpha\beta} = \frac{\partial f_\alpha}{\partial n_\beta}$$

where $\alpha$ and $\beta$ refer to the location of a block or submatrix. The
resulting block-tridiagonal matrix is stored in banded form (`store_bandM`) and solved with
`scipy.linalg.solve_banded`.

## Step-size control and convergence

The truncation error $\mathcal{E} = |n_{k+1} - n^*_{k+1}|$, where $n^*_{k+1} = n_k +
\Delta t\, k_1$ is the first-order estimate, is used to adjust the stepsize in
[`Ros2.step_size`](../Reference/api/op.md#vulcan.op.Ros2.step_size) according to

$$\Delta t_{k+1} = \Delta t_k \times \text{clamp}\!\left(0.9\left(\mathcal{T}/\mathcal{E}
\right)^{0.5},\ h_{\min},\ h_{\max}\right)$$

where $\mathcal{T}$ is the desired relative error tolerance (`rtol`) and 0.9 is a safety
factor. The scaling factor is clamped to the range $[0.5,\ 2]$ to prevent excessively
large or small step changes, and the resulting $\Delta t$ is further bounded to
$[10^{-10},\ 10^{18}]$ s. After each timestep, if $\mathcal{E}$ is greater than
$\mathcal{T}$ the solution is rejected and the stepsize is reduced; after a successful
integration, the previous scheme is applied to adjust the stepsize.

Steady state is declared in
[`op.Integration.conv`](../Reference/api/op.md#vulcan.op.Integration.conv) when the
long-term relative variation and its rate both fall below threshold values. After the $k$th
timestep and total integration time $\tau$, we compute the relative variation

$$\Delta\hat{n} \equiv \max_{i,j} \frac{|n_{i,j,k} - n_{i,j,k'}|}{X_{i,j,k}},
\quad \Delta t \equiv t_k - t_{k'}$$

where $X_{i,j,k}$ is the mixing ratio of species $i$ at level $j$, and $k'$ refers to the
timestep at $f\tau$ (i.e., the variation of the solution from $f\tau$ to $\tau$ is
examined). Species with mixing ratios below `mtol_conv` or number densities below `atol`
are excluded from the maximum to avoid spurious triggering on negligible trace species.

We declare steady state when either the primary or relaxed convergence condition is
satisfied:

$$\left(\Delta\hat{n} < \delta \;\text{ and }\; \frac{\Delta\hat{n}}{\Delta t} < \epsilon
\right) \quad\text{or}\quad \left(\Delta\hat{n} < \delta_{\min} \;\text{ and }\;
\frac{\Delta\hat{n}}{\Delta t} < \epsilon_{\min}\right)$$

where $\delta$ and $\epsilon$ are set by `yconv_cri` and `slope_cri`, and $\delta_{\min}$
and $\epsilon_{\min}$ are set by `yconv_min` and `slope_min`. The relaxed slope threshold
$\epsilon_{\min}$ is derived from the local diffusion timescale,

$$\epsilon_{\min} = \min\!\left(\min_{j} \frac{K_{zz,j}}{(0.1\, H_{p,j})^2},\
10^{-8}\right)$$

bounded from below at $10^{-10}\ \mathrm{s}^{-1}$, so that the convergence criterion
naturally adapts to the diffusive timescale of each atmosphere. With photochemistry
enabled, both conditions additionally require that the actinic-flux change falls below
`flux_cri`.

[^tsai2021]: Tsai, S.-M., Malik, M., Kitzmann, D., et al. (2021). A comparative study of atmospheric chemistry with VULCAN. *The Astrophysical Journal, 923*(2), 264. https://doi.org/10.3847/1538-4357/ac29bc

[^tsai2017]: Tsai, S.-M., Lyons, J. R., Grosheintz, L., et al. (2017). VULCAN: An open-source, validated chemical kinetics Python code for exoplanetary atmospheres. *The Astrophysical Journal Supplement Series, 228*(2), 20. https://doi.org/10.3847/1538-4365/228/2/20

[^ch1987]: Chamberlain, J. W., & Hunten, D. M. (1987). *Theory of Planetary Atmospheres*, 2nd ed. Academic Press.

[^bj2017]: Brasseur, G. P., & Jacob, D. J. (2017). *Modeling of Atmospheric Chemistry*. Cambridge University Press.

[^verwer1997]: Verwer, J. G., Spee, E. J., Blom, J. G., & Hundsdorfer, W. (1997). A second-order Rosenbrock method applied to photochemical dispersion problems. *SIAM Journal on Scientific Computing, 20*(4), 1456–1480. https://doi.org/10.1137/S1064827597326651

[^press2007]: Press, W. H. (2007). Numerical recipes 3rd edition: The art of scientific computing. Cambridge university press. https://doi.org//10.5555/1403886

[^mm1972]: Marrero, T. R., & Mason, E. A. (1972). Gaseous diffusion coefficients. *Journal of Physical and Chemical Reference Data, 1*(1), 3–118. https://doi.org/10.1063/1.3253094

[^bk1973]: Banks, P. M., & Kockarts, G. (1973). *Aeronomy*. Academic Press.