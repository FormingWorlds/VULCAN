# Output files

VULCAN writes its final state to `output_dir` when the integration ends. The format is chosen
with `save_type` in `config.py`, and the file name with `out_name`:

| `save_type` | File(s) | Contents |
|---|---|---|
| `'netcdf'` (default) | `<out_name>.nc` | Self-describing NetCDF4 file, described below. |
| `'pickle'` | `<out_name>.pkl` | Python pickle of `{'variable', 'atm', 'parameter'}` dicts. Required by `ini_mix = 'vulcan_ini'`, `atm_type = 'vulcan_ini'` and the scripts in `tools/`. |
| `'plain'` | `<out_name>.var`, `.atm`, `.par` | Human-readable text dumps (`str()` of each field). Not intended for parsing. |

`save_type` may also be a list, e.g. `['netcdf', 'pickle']`, to write several formats from
one run. The configuration is always dumped to `vulcan_cfg.txt` as well.

## NetCDF file

The file can be read with any NetCDF library, for example:

```python
import xarray as xr

ds = xr.open_dataset('output/example.nc')
ds.mixing_ratio.sel(species='CH4')
```

or with `netCDF4.Dataset('output/example.nc')`.

All quantities are in the **CGS units** used inside VULCAN, and every variable carries
`units` and `long_name` attributes. Along the vertical dimensions, **index 0 is the bottom of
the atmosphere** (highest pressure).

### Dimensions

| Dimension | Size | Present |
|---|---|---|
| `layer` | `nz` (layer centres) | always |
| `level` | `nz + 1` (layer interfaces) | always |
| `species` | number of species in the network | always |
| `reaction` | number of reactions `nr`, forward and reverse | always |
| `atom` | `len(atom_list)` | always |
| `wavelength` | number of spectral bins | `use_photo = True` |
| `time` | number of stored steps | `save_evolution = True` |

### Atmospheric structure

| Variable | Dimensions | Units | Description |
|---|---|---|---|
| `pressure` | layer | dyn cm-2 | Pressure at layer centres. |
| `pressure_level` | level | dyn cm-2 | Pressure at interfaces. |
| `temperature` | layer | K | Temperature. |
| `altitude` | layer | cm | Altitude at layer centres. |
| `altitude_level` | level | cm | Altitude at interfaces. Zero at the 1 bar level for gas giants, at the surface if `rocky`. |
| `Kzz` | level | cm2 s-1 | Eddy diffusion coefficient. Defined on interior interfaces only; the two boundary levels are fill values. |
| `vz` | level | cm s-1 | Vertical advection velocity, as `Kzz`. |
| `mean_molecular_weight` | layer | g mol-1 | Mean molecular weight. |
| `gravity` | layer | cm s-2 | Gravitational acceleration. |
| `scale_height` | layer | cm | Pressure scale height. |
| `number_density_total` | layer | cm-3 | Total number density. |
| `third_body_density` | layer | cm-3 | Third-body density M used in the rate equations. |

### Species

| Variable | Dimensions | Units | Description |
|---|---|---|---|
| `species` | species | | Species names (strings). |
| `molecular_mass` | species | g mol-1 | Molar mass of each species. |
| `number_density` | layer, species | cm-3 | Final number densities. |
| `mixing_ratio` | layer, species | 1 | Final volume mixing ratios. |
| `number_density_initial` | layer, species | cm-3 | Initial number densities. |
| `time` | time | s | Model time of stored steps (`save_evolution` only). |
| `number_density_evolution` | time, layer, species | cm-3 | Number densities every `save_evo_frq` steps (`save_evolution` only). |

### Reactions

Reactions are numbered 1 to `nr` as in the (renumbered) network file. Odd IDs are forward
reactions, and each even ID is the reverse of the reaction before it.

| Variable | Dimensions | Units | Description |
|---|---|---|---|
| `reaction_id` | reaction | | Reaction index. |
| `reaction_equation` | reaction | | Equation string, e.g. `H + H + M -> H2 + M`, written explicitly for reverse reactions too. |
| `reaction_type` | reaction | | `thermal`, `photolysis`, `photoionisation` or `condensation`. Reverse reactions inherit the type of their forward reaction. |
| `reaction_is_reverse` | reaction | | 1 for reverse reactions, else 0. |
| `reaction_order` | reaction | | Number of reactants including M. |
| `rate_coefficient` | layer, reaction | s-1, cm3 s-1 or cm6 s-1 | Rate coefficient k, with units set by `reaction_order` (1, 2 or 3). Photolysis coefficients include `f_diurnal`. Removed and unused reverse reactions are zero. |
| `reaction_rate` | layer, reaction | cm-3 s-1 | Rate k × Π(reactant number densities), with M taken from `third_body_density`. |

The net chemical tendency of each species is recovered by summing `reaction_rate` weighted by
its stoichiometric coefficients (products minus reactants). This reproduces the chemical
source term that VULCAN integrates.

### Radiation (`use_photo = True`)

| Variable | Dimensions | Units | Description |
|---|---|---|---|
| `wavelength` | wavelength | nm | Spectral grid (non-uniform: `dbin1` below `dbin_12trans`, `dbin2` above). |
| `actinic_flux` | layer, wavelength | photons cm-2 s-1 nm-1 | Actinic flux at layer centres. |
| `stellar_flux_toa` | wavelength | erg cm-2 s-1 nm-1 | Stellar flux at the top of the atmosphere. |
| `direct_flux` | level, wavelength | erg cm-2 s-1 nm-1 | Attenuated direct beam, not multiplied by cos(zenith angle). |
| `diffuse_flux_up` | level, wavelength | erg cm-2 s-1 nm-1 | Upward diffuse flux. |
| `diffuse_flux_down` | level, wavelength | erg cm-2 s-1 nm-1 | Downward diffuse flux. |
| `optical_depth` | level, wavelength | 1 | Vertical optical depth from the top of the atmosphere. |

### Element conservation

| Variable | Dimensions | Units | Description |
|---|---|---|---|
| `atom` | atom | | Element names. |
| `atom_initial` | atom | cm-3 | Initial element number density, summed over layers. |
| `atom_sum` | atom | cm-3 | Final element number density, summed over layers. |
| `atom_loss` | atom | 1 | Fractional change, (`atom_sum` − `atom_initial`) / `atom_initial`. |
| `atom_conden` | atom | cm-3 | Element number density in condensates. |

### Global attributes

| Attribute | Description |
|---|---|
| `vulcan_version`, `date_created` | Code version and UTC creation time. |
| `network`, `atom_list` | Chemical network file and elements. |
| `end_case`, `termination` | How the run ended: 1 converged, 2 `runtime` limit, 3 `count_max` limit. |
| `steps`, `model_time_s`, `dt_s` | Number of steps, final model time and step size. |
| `longdy`, `longdydt` | Final values of the convergence measures. |
| `wall_time_s` | Wall-clock run time. |
| `Rp`, `gs`, `r_star`, `orbit_radius`, `sl_angle`, `f_diurnal`, `rocky`, `use_photo`, ... | Selected configuration values (booleans as 0/1). |
| `config` | The full configuration at the end of the run, as a JSON string. |

!!! note "Consistency of rate coefficients"
    Photolysis rate coefficients are only updated every `final_update_photo_frq` steps, so
    they can lag the final abundances by a few steps. When coupled to AGNI
    (`agni_call_frq > 0`), rate coefficients are not re-evaluated after the last temperature
    update.
