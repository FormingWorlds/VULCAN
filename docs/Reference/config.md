# Configuration parameters

Every VULCAN run is controlled by the `Config` class in `src/vulcan/config.py`. Each setting
is an attribute (`self.<name>`); edit the file, then run the model. The defaults shown below
are the shipped configuration, which describes the hot Jupiter **HD 189733b**.

!!! info "Conventions"
    - **Units are cgs** unless stated otherwise: pressure in dyne cm⁻² (1 bar = 10⁶ dyne cm⁻²),
      gravity in cm s⁻², $K_{zz}$ and diffusion in cm² s⁻¹, fluxes in molecules cm⁻² s⁻¹,
      number density in cm⁻³, time in seconds.
    - Stellar/orbital quantities are the exception: `r_star` in $R_\odot$, `orbit_radius` in au,
      `Rp` in cm, `sl_angle` in radians.
    - Calling `vulcan_cfg.write_file()` dumps the **active** configuration to
      `output/vulcan_cfg.txt`, which is the surest record of what a run actually used.

---

## Chemistry & network

| Parameter | Default | Description |
|---|---|---|
| `atom_list` | `['H', 'O', 'C']` | Elements tracked by the run. **Must match the network and composition file.** |
| `network` | `thermo/CHO_photo_network.txt` | Path to the reaction network file (prefixed with `VULCAN_DIR`). |
| `remove_list` | `[]` | Reaction indices to switch off, given in forward/reverse pairs, e.g. `[1, 2]`. |

!!! warning "Keep `atom_list` and `network` consistent"
    Initialization, the FastChem setup, and the element-conservation check all iterate over
    `atom_list`. Switching to a network with more elements (e.g. adding N or S) requires
    updating `atom_list` **and** the relevant elemental abundances.

---

## Atmospheric grid & background

| Parameter | Default | Description |
|---|---|---|
| `nz` | `50` | Number of vertical layers. Higher = finer resolution but slower. |
| `P_b` | `1e8` | Pressure at the bottom boundary (dyne cm⁻²); default ≈ 100 bar. |
| `P_t` | `1e-2` | Pressure at the top boundary (dyne cm⁻²); default ≈ 10⁻⁸ bar. |
| `atm_base` | `'H2'` | Background gas (`'H2'`, `'N2'`, `'O2'`, `'CO2'`); sets the molecular-diffusion reference coefficients and viscosity. |
| `rocky` | `False` | If `False`, gravity is referenced at 1 bar (gas giants); if `True`, `gs` is the surface gravity. |
| `gs` | `2140.0` | Gravitational acceleration (cm s⁻²) at the surface or at 1 bar. |
| `Rp` | `1.138 × R_J` | Planetary radius (cm). |

---

## Temperature–pressure profile

| Parameter | Default | Description |
|---|---|---|
| `atm_type` | `'file'` | T–P source: `'file'`, `'analytical'`, `'isothermal'`, `'vulcan_ini'`, `'table'`. |
| `atm_file` | `atm/atm_HD189_Kzz.txt` | Atmosphere file (P, T, and optionally $K_{zz}$, $v_z$) when `atm_type='file'`. |
| `Tiso` | `1000.0` | Isothermal temperature (K) when `atm_type='isothermal'`. |
| `para_anaTP` | `[120, 1500, 0.1, 0.02, 1.0, 1.0]` | Analytic T(P) parameters `[T_int, T_irr, ka_0, ka_s, beta_s, beta_l]` (Heng et al. 2014) when `atm_type='analytical'`. |

---

## Input files

| Parameter | Default | Description |
|---|---|---|
| `sflux_file` | `atm/stellar_flux/sflux-HD189_Moses11.txt` | Stellar flux at the stellar surface (scaled to the planet internally). |
| `top_BC_flux_file` | `atm/BC_top.txt` | Top boundary-condition file (read only if `use_topflux`). |
| `bot_BC_flux_file` | `atm/BC_bot.txt` | Bottom boundary-condition file (read if `use_botflux` or `use_fix_sp_bot`). |

---

## Initialization & elemental abundances

| Parameter | Default | Description |
|---|---|---|
| `ini_mix` | `'eq'` | Initial composition: `'eq'` (FastChem equilibrium), `'const_mix'`, `'const_lowt'`, `'vulcan_ini'`, `'table'`. |
| `const_mix` | `{...}` | Well-mixed initial mixing ratios, used when `ini_mix='const_mix'`. |
| `vul_ini` | `'_unset_'` | Path to a previous run (pickle/table) for `ini_mix='vulcan_ini'`/`'table'` and `atm_type='vulcan_ini'`/`'table'`. |
| `use_solar` | `True` | If `True`, FastChem uses built-in solar abundances; the custom `*_H` values below are ignored. |
| `fastchem_met_scale` | `1.0` | Metallicity scaling for elements not tracked by VULCAN (used when `use_solar=False`). |
| `O_H`, `C_H`, `N_H`, `S_H` | `5.37e-4`, `2.95e-4`, `7.08e-5`, `1.41e-5` | Custom elemental abundances relative to H (used only when `ini_mix='eq'` and `use_solar=False`). |
| `He_H` | `0.0838` | Helium-to-hydrogen ratio. |

!!! warning "Custom abundances require `use_solar = False`"
    With the default `use_solar = True`, `C_H`, `O_H`, `N_H`, `S_H` are not used — FastChem
    initializes from solar. The C/O ratio is simply `C_H / O_H`.

---

## Photochemistry & radiation

| Parameter | Default | Description |
|---|---|---|
| `use_photo` | `True` | Enable photochemistry (optical depth, actinic flux, photolysis). |
| `use_ion` | `False` | Enable photoionization and charge balance. |
| `r_star` | `0.805` | Stellar radius ($R_\odot$), for flux scaling. |
| `orbit_radius` | `0.03142` | Orbital distance (au), for flux scaling. |
| `sl_angle` | `48° (rad)` | Stellar zenith angle θ (radians). |
| `f_diurnal` | `1.0` | Diurnal-averaging factor (`1.0` tidally locked, `0.5` rotating, e.g. Earth). |
| `scat_sp` | `['H2', 'He']` | Species included in Rayleigh scattering. |
| `T_cross_sp` | `[]` | Species using temperature-dependent UV cross sections (need per-temperature files). |
| `edd` | `0.5` | First Eddington coefficient $\bar\epsilon$ for the two-stream radiative transfer. |
| `dbin1` | `0.1` | Wavelength bin width (nm) below `dbin_12trans`. |
| `dbin2` | `2.0` | Wavelength bin width (nm) above `dbin_12trans`. |
| `dbin_12trans` | `240.0` | Wavelength (nm) where the bin width switches from `dbin1` to `dbin2`. |
| `ini_update_photo_frq` | `100` | Initial step interval for recomputing the actinic flux/optical depth. |
| `final_update_photo_frq` | `5` | Step interval near convergence. |

---

## Vertical mixing & transport

| Parameter | Default | Description |
|---|---|---|
| `use_Kzz` | `True` | Enable eddy diffusion. |
| `Kzz_prof` | `'file'` | $K_{zz}$ source: `'file'`, `'const'`, `'Pfunc'` (Tsai 2021 form), `'JM16'` (Moses 2016 form). |
| `const_Kzz` | `1e10` | Constant $K_{zz}$ (cm² s⁻¹), used **only** when `Kzz_prof='const'`. |
| `K_max`, `K_p_lev` | `1e5`, `0.1` | Parameters for `Kzz_prof='Pfunc'`. |
| `use_vz` | `True` | Enable vertical advection. |
| `vz_prof` | `'const'` | Vertical-wind source: `'const'` or `'file'`. |
| `const_vz` | `0.0` | Constant vertical velocity (cm s⁻¹) when `vz_prof='const'`. |
| `use_moldiff` | `True` | Include molecular + thermal diffusion. |
| `use_vm_mol` | `False` | Upwind molecular-diffusion scheme — **under development, must stay `False`** (raises an error otherwise). |
| `update_frq` | `100` | Step interval for updating mean molecular weight and layer thicknesses. |

!!! warning "`const_Kzz` only applies with `Kzz_prof = 'const'`"
    Under the default `Kzz_prof = 'file'`, the value comes from `atm_file` and `const_Kzz` is
    ignored. (The `'JM16'` profile additionally expects a `K_deep` value to be set.)

---

## Boundary conditions

| Parameter | Default | Description |
|---|---|---|
| `use_topflux` | `False` | Apply constant top-boundary fluxes from `top_BC_flux_file`. |
| `use_botflux` | `False` | Apply bottom fluxes and deposition velocities from `bot_BC_flux_file`. |
| `use_fix_sp_bot` | `{}` | Fixed mixing ratios at the bottom for selected species. |
| `diff_esc` | `['H']` | Species given diffusion-limited escape at the top boundary. |
| `max_flux` | `1e13` | Cap on the diffusion-limited escape flux (molecules cm⁻² s⁻¹). |

---

## Condensation & settling

| Parameter | Default | Description |
|---|---|---|
| `use_condense` | `False` | Enable condensation/evaporation. |
| `use_settling` | `False` | Enable gravitational particle settling. |
| `condense_sp` | `[]` | Species allowed to condense (e.g. `'H2O'`, `'NH3'`, `'S2'`, `'S8'`). |
| `non_gas_sp` | `[]` | Condensate (non-gaseous) species, excluded from gas-phase mixing ratios. |
| `r_p` | `{'H2O_l_s': 5e-3}` | Particle radius (cm) per condensate. |
| `rho_p` | `{'H2O_l_s': 1}` | Particle density (g cm⁻³) per condensate. |
| `humidity` | `1.0` | Relative-humidity multiplier on the H₂O saturation. |
| `start_conden_time`, `stop_conden_time` | `0.0`, `1e5` | Simulation times (s) bracketing the condensation/fixing phase. |
| `fix_species`, `fix_species_time` | `[]`, `0` | Species to hold fixed after condensation equilibrium, and when. |

---

## AGNI climate coupling (optional)

| Parameter | Default | Description |
|---|---|---|
| `agni_call_frq` | `0` | Step interval for calling AGNI; `0` disables the coupling. |
| `solve_rce` | `False` | If `True`, solve radiative–convective equilibrium; else apply a prescribed profile and fill $K_{zz}$. |
| `spectral_file` | `'greygas'` | AGNI spectral file (or grey-gas solution). |
| `use_rayleigh` | `False` | Include Rayleigh scattering in AGNI. |
| `surf_albedo` | `0.0` | Surface albedo. |
| `Tsurf_guess` | `2000.0` | Initial surface-temperature guess (K). |
| `agni_atol`, `agni_rtol` | `1e-3`, `1e-3` | AGNI convergence tolerances. |
| `agni_nlev` | `60` | Number of AGNI atmospheric levels. |

---

## Steady-state & convergence

| Parameter | Default | Description |
|---|---|---|
| `yconv_cri` | `0.01` | Threshold on the long-term relative change of mixing ratios. |
| `slope_cri` | `1e-4` | Threshold on the rate of that change. |
| `yconv_min` | `0.1` | Relaxed convergence threshold (paired with an internal minimum slope). |
| `flux_cri` | `0.1` | Threshold on the actinic-flux change for convergence (photochemistry runs). |
| `flux_atol` | `1.0` | Actinic-flux floor below which changes are ignored (photons cm⁻² s⁻¹ nm⁻¹). |
| `st_factor` | `0.5` | Fraction of the runtime used as the convergence look-back window. |
| `conv_step` | `500` | Maximum number of steps to look back for convergence. |
| `conver_ignore` | `[]` | Species excluded from the convergence test (e.g. sinkless species like HC₃N). |

---

## ODE solver & time stepping

| Parameter | Default | Description |
|---|---|---|
| `ode_solver` | `'Ros2'` | Integration scheme (second-order Rosenbrock). |
| `dttry` | `1e-6` | Initial time step (s). |
| `dt_min`, `dt_max` | `1e-8`, `1e18` | Bounds on the time step (s). |
| `dt_var_min`, `dt_var_max` | `0.5`, `2.0` | Per-step shrink/grow factors. |
| `trun_min` | `1e2` | Minimum simulated time (s) before convergence is allowed. |
| `runtime` | `1e22` | Maximum simulated time (s). |
| `count_min`, `count_max` | `120`, `30000` | Minimum/maximum number of steps. |
| `rtol` | `0.25` | Relative tolerance controlling Rosenbrock step size. |
| `use_adapt_rtol` | `True` | Adapt `rtol` from element loss during the run. |
| `rtol_min`, `rtol_max` | `0.02`, `2.5` | Bounds for the adaptive `rtol`. |
| `post_conden_rtol` | `0.1` | `rtol` used after condensable species are fixed. |
| `atol` | `1e-1` | Absolute number-density floor for error control. |
| `mtol` | `1e-22` | Minimum mixing ratio considered. |
| `mtol_conv` | `1e-20` | Minimum mixing ratio considered for convergence. |
| `pos_cut`, `nega_cut` | `0`, `-1.0` | Clipping thresholds for small/negative number densities. |
| `loss_eps` | `1e-1` | Per-step element-conservation tolerance (steps exceeding it are rejected). |

---

## Output, logging & plotting

| Parameter | Default | Description |
|---|---|---|
| `output_dir` | `output/` | Output directory (created at runtime). |
| `plot_dir` | `output/plot/` | Plot directory. |
| `out_name` | `'example.pkl'` | Output filename. |
| `clean_output` | `True` | Wipe the output directory at startup. |
| `output_humanread` | `False` | Write human-readable text output instead of a pickle. |
| `save_evolution` | `False` | Save the time evolution (`y_time`, `t_time`). |
| `save_evo_frq` | `10` | Step interval for storing evolution data. |
| `log_level` | `'INFO'` | Logging verbosity (`'DEBUG'`, `'INFO'`, `'WARNING'`, `'ERROR'`). |
| `use_print_prog`, `print_prog_num` | `True`, `20` | Print progress, and how often (steps). |
| `use_print_delta` | `False` | Print the largest truncation error each progress step. |
| `plot_TP` | `True` | Plot the T–P/$K_{zz}$ profile at startup. |
| `use_live_plot`, `live_plot_frq` | `True`, `50` | Live mixing-ratio plotting and its step interval. |
| `use_live_flux` | `False` | Live flux plotting. |
| `use_plot_end`, `use_plot_evo` | `True`, `True` | Plot final profiles / time evolution at the end. |
| `use_flux_movie` | `False` | Save flux frames for a movie. |
| `plot_height` | `False` | Plot against height instead of pressure. |
| `plot_spec` | `['H2','H','H2O','CH4','CO','CO2','C2H2']` | Species to plot. |
| `plot_dpi` | `130` | Figure resolution. |
| `y_time_freq` | `1` | Step interval for storing time-series data. |

!!! warning "`clean_output` on networked filesystems"
    On NFS home directories (e.g. HPC clusters), `clean_output = True` can fail with
    `OSError: [Errno 39] Directory not empty` because the open log file is silly-renamed during
    the directory wipe. Set `clean_output = False`, or point `output_dir`/`plot_dir` at local
    scratch.