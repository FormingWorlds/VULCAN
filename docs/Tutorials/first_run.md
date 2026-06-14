# Tutorial: first VULCAN run

This tutorial walks through your first VULCAN computation. The default configuration is the
hot Jupiter **HD 189733b**, a hydrogen-dominated atmosphere with photochemistry, so by the
end you will have run a full photochemical-kinetics model and produced vertical mixing-ratio
profiles.

If you haven't installed VULCAN yet, follow the [installation guide](../How-to/installation.md) first, and make sure FastChem is compiled (the default initialization uses it).

!!! info "What you'll do"
    - Run the default HD 189733b model with `run_vulcan.py`
    - Find and inspect the output
    - Change a few configuration parameters (mixing, C/O ratio, plotting)
    - Plot the resulting profiles

---

## 1. Run the default model

From the top-level `VULCAN/` directory:

```sh
python run_vulcan.py
```

This will:

- generate `chem_funs.py` from the default reaction network,
- initialize the composition at chemical equilibrium with FastChem (`ini_mix = 'eq'`),
- integrate the atmosphere to steady state, with live plotting of the mixing ratios on by
  default.

It runs until the convergence criteria are met, typically a few minutes. When it finishes, check the output directory:

```sh
ls output/
```

You should find the pickled result `example.pkl` (the name is set by `out_name`), the
diagnostic plots under `output/plot/`, a copy of the configuration `vulcan_cfg.txt`, and the
run log `vulcan.log`.

!!! note "The default model is a hot Jupiter, not Earth"
    The defaults in `config.py` describe HD 189733b: the T–P and $K_{zz}$ profile come from
    `atm/atm_HD189_Kzz.txt`, the stellar spectrum from `sflux-HD189_Moses11.txt`, with
    $g = 2140$ cm s$^{-2}$, $R_p = 1.138\,R_\mathrm{J}$, and an orbital distance of 0.03142 au.
    To model a different planet you change the atmosphere file, stellar flux, and planetary
    parameters.

---

## 2. Modify the configuration

Simulation parameters live in the `Config` class in `config.py` (each as `self.<name>`).
Edit a value, then re-run `python run_vulcan.py`.

### Weaker vertical mixing

The eddy diffusion profile source is chosen by `Kzz_prof`. The default reads $K_{zz}$ from the
atmosphere file (`Kzz_prof = 'file'`), in which case `const_Kzz` is **ignored**. To impose a
single constant value you must set both:

```python
self.Kzz_prof  = 'const'
self.const_Kzz = 1.E7      # cm^2 s^-1
```

!!! warning "`const_Kzz` only applies when `Kzz_prof = 'const'`"
    Setting `const_Kzz` alone has no effect under the default `Kzz_prof = 'file'`. Change `Kzz_prof` as well.

### Carbon-rich atmosphere (C/O ratio)

The C/O ratio is just `C_H / O_H`. But the custom elemental abundances are only used when
solar abundances are switched **off**:

```python
self.use_solar = False     # otherwise FastChem uses built-in solar abundances
self.O_H       = 5.37e-4
self.C_H       = 5.37e-4   # C/O = 1.0  (raise C_H above O_H for C/O > 1)
```

Set `use_solar = False` for your `*_H` values to take effect.

### Disable live plotting (faster)

```python
self.use_live_plot = False
```

---

## 3. Plot the results

A plotting helper lives in `tools/`:

```sh
cd tools/
python plot_vulcan.py <vulcan_output.pkl> <species> <plot_name> [-h for height axis]
```

For example, to plot the default species from the HD 189733b run:

```sh
python plot_vulcan.py ../output/example.pkl H2,H,H2O,CH4,CO,CO2,C2H2 hd189b
```

This reads the pickled output and writes `hd189b.png` to `plot_dir`, showing mixing ratio
versus pressure (or height with `-h`), with the initial (equilibrium) abundances shown as
dotted lines.

---

## 4. Speed tips

- **Skip chemistry regeneration.** If you only changed `config.py` (not the reaction network),
  reuse the compiled `chem_funs.py`:
  ```sh
  python run_vulcan.py -n
  ```
  Only use `-n` if you have **not** edited the network, see
  [editing the chemical network](../How-to/chem_network.md).
- **Turn off live plotting:** `self.use_live_plot = False`.
- **Lower the vertical resolution.** The number of layers is `nz` (default 50); fewer layers
  run faster at the cost of vertical detail:
  ```python
  self.nz = 40
  ```

!!! tip "Layer count `self.nz`"
    The layer count is `self.nz`. Increasing it (e.g. to 100) makes the run *slower* and more finely resolved; to speed things up, decrease it.

---

## Troubleshooting

!!! failure "FastChem not found"
    Compile it once, then re-run:
    ```sh
    cd fastchem_vulcan && make && cd ..
    ```
    VULCAN looks for the `fastchem` binary inside `fastchem_vulcan/`.

!!! failure "Missing Python modules"
    Reinstall the package (and its dependencies) in editable mode from `VULCAN/`:
    ```sh
    pip install -e .
    ```

!!! failure "`OSError: [Errno 39] Directory not empty` on a cluster"
    On NFS home directories, the automatic wipe of `output/` collides with the open log file.
    Either set `self.clean_output = False` in `config.py`, or point `self.output_dir` (and
    `self.plot_dir`) at local scratch.

!!! failure "Simulation too slow"
    Disable live plotting (`self.use_live_plot = False`) and/or reduce `self.nz`.

---

## Things to try from here

- **Switch planet:** change `atm_file`, `sflux_file`, and the planetary parameters in
  `config.py` to model a different target.
- **Compare mixing strengths:** run with `Kzz_prof = 'const'` at a few values of `const_Kzz`
  and compare the quench levels of CH$_4$ and NH$_3$.
- **Sweep the C/O ratio:** with `use_solar = False`, step `C_H` across solar-to-carbon-rich
  values and watch the hydrocarbon and oxide abundances respond.
- **Toggle photochemistry:** set `use_photo = False` for a thermochemistry-only run and
  compare the upper-atmosphere profiles to the full photochemical result.