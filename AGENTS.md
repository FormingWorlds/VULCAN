# VULCAN agent instructions

VULCAN computes the chemical kinetics of planetary atmospheres: photochemistry, vertical mixing, condensation and the coupling to AGNI. PROTEUS selects it as `atmos_chem.module = vulcan`. Before a first edit:

- Run every command from the repository root. Tests for `src/vulcan/<file>.py` go in `tests/test_<file>.py`; the test rules are in `tests/AGENTS.md`.
- These commands decide whether a change is ready (CI runs the tests with coverage and pre-commit on the changed files):

```bash
python -m pip install -e '.[develop]'    # the package alone: pip install -e .
(cd fastchem_vulcan && make)        # FastChem, needed by ini_mix = 'eq'; not make -j
pytest tests/
ruff check <changed .py files>
pre-commit run --files <changed files>
```

Run the model with `python run_vulcan.py` (the default demo); `python run_vulcan.py -n` skips the regeneration of the chemistry functions, which is safe only when the network is unchanged. One test: `pytest tests/test_imports.py::test_package_imports -q`; a group by name: `pytest -k import`. `ruff check .` lints the whole tree. Coverage as CI measures it: `pytest --cov-report term --cov-report json:coverage.json --cov=vulcan tests/`. Docs: `pip install -e '.[docs]'`, then `zensical build --clean`.

<!-- fwl-core:begin sha256=1adfa510a4f117ec -->
## PROTEUS ecosystem rules

PROTEUS couples separate module repositories into one model of the evolution of rocky planets and their atmospheres. Its config selects one module per process; the config values are the lowercase names:

- `interior_energetics.module`: `aragog`, `spider` (or the built-in `boundary`); `interior_struct.module`: `zalmoxis`, `spider`
- `atmos_clim.module`: `agni`, `janus`, both with SOCRATES spectral radiative transfer
- `outgas.module`: `calliope`, `atmodeller`; `atmos_chem.module`: `vulcan`
- `star.module`: `mors`; `escape.module`: `zephyrus`, `boreas`; `orbit.module` (tides): `obliqua`, `lovepy`
- fwl-io downloads the reference data the modules use

A module change to anything PROTEUS calls or reads (a function signature, a config key, an output column, a unit) can break the coupled model while the module's own tests pass. Search `src/proteus/` in PROTEUS for the name before you change it, and name the affected PROTEUS call sites in the pull request.

### Physics and numerics

- Units differ between and within modules: MORS and VULCAN use cgs; CALLIOPE and atmodeller work with pressures in bar; AGNI takes surface pressure in bar at setup and holds Pa inside; PROTEUS states the unit of every `hf_row` key in `GetHelpfileKeys` (`src/proteus/utils/coupler.py`: atmospheric pressures in bar, interior pressures in Pa, model time in years, orbital periods in s), and its config reference pages state each key's unit (for example stellar mass in M_sun, initial partial surface pressures in bar, stellar age in Gyr). State the unit of every physical quantity in its docstring and convert explicitly at the boundary: a unit mismatch between two modules passes the tests of both and shows only in the coupled run.
- A physics test must fail for the most plausible wrong formula: check a conservation law, a bound, a monotonicity, or a published or analytical value, and assert that the wrong result falls outside the tolerance.
- Take each physical constant from one source per repository (in Python `scipy.constants` or the module's constants file). Two retyped values of one constant differ at round-off and hide real differences between code paths.
- Do not loosen a solver tolerance, a conservation check or a clamp to make a run or a test pass. Find the cause; a loosened check also hides the next defect.

### Validation

A new or changed physics routine gets the validation cases named in its confirmed plan: an analytical law, a published figure or table, or a cross-code check. Each case is a `reference_pinned` test that runs the real code on real data, with no mocks, asserts the tolerance, and writes its light and dark plot through the `validation_plot` fixture of proteus-mpl. `pytest -m reference_pinned --validation-plots` makes all plots again. Line plots are SVG, maps are AVIF; each file is at most 60 kB, under `docs/assets/validation/`. Each case has a subsection on its page in the Validation tab, with the reference, the tolerance, the test id and the command. Data read from a published figure comes from 2 independent figlab reads that agree within the stated error, as a CSV with its provenance under `tests/validation/data/`.

### Code

- Keep an inline comment block to 2 lines or fewer, and never more than 4, because a long block drifts from the code it describes. An explanation that needs more goes in the docstring, the commit message or the pull-request description.

### Branches and pull requests

- Work on a feature branch `<initials>/<short-description>`; `main` changes only through a reviewed pull request.
- Fill every section of the repository's pull-request template, where it has one.
- Before you push, run the checks listed at the top of this `AGENTS.md` and in `tests/AGENTS.md`.

### Where knowledge goes

Rules for every contributor go in the `AGENTS.md` files, the reason for a change in its commit message and pull-request description, and the scientific validation of a module in its `docs/Validation/` pages, where the repository has them. Do not add memory, notes or decision-log files to the repository: nobody maintains them, and they go stale.
<!-- fwl-core:end -->

<!-- fwl-voice:begin sha256=b4fe24baeef8b734 -->
### Commit messages and public text

Commit messages, pull-request text, code comments, docstrings, test names, test skip reasons, parametrize ids, log strings that ship with the code, and CI job and step names describe the change and the current state of the code. They name no tool used to write the change, carry no tool-attribution trailer, use no internal plan, phase or work-group labels, and use no em or en dashes (a page range in a citation is the exception).
<!-- fwl-voice:end -->

## Architecture

`run_vulcan.py` at the repository root is the command-line entry point: it puts `src/` first on `sys.path` and calls `vulcan.vulcan.run_cli()`. `src/vulcan/vulcan.py` runs a model in this order:

1. Build the `Config` (`config.py`) and regenerate the chemistry functions (`make_chem_funs.make_all`) unless `-n` is given.
2. Build the state containers `store.Variables`, `store.AtmData` and `store.Parameters`.
3. Build the atmosphere (`build_atm.Atm`): pressure levels, temperature, Kzz, molecular diffusion, and the boundary fluxes.
4. Read the reactions (`op.ReadRate`) and set the initial abundances (`build_atm.InitialAbun`).
5. Integrate with `op.Integration` and the `op.Ros2` solver, with photochemistry when `use_photo` is set (the cross sections and the stellar flux at the top, `build_atm.Atm.read_sflux`, are read first) and AGNI coupling (`agni.py`) when `agni_call_frq > 0`.
6. Write the output (`op.Output.save_out`) and the plots (`plot_dir`, default `output/plot/`).

Module roles: `config.py` holds every runtime setting; `paths.py` the absolute paths to the repository data (thermo, cross sections, FastChem, AGNI, the generated chemistry file); `make_chem_funs.py` writes `chem_funs.py` from the network and the thermodynamic data; `op.py` reads the rates, integrates, solves and writes output; `build_atm.py` builds the atmosphere and the initial abundances (with FastChem equilibrium); `store.py` holds the state passed through the solver.

## Conventions

- `chem_funs.py` is generated: edit the network or the thermodynamic input, then run without `-n` so that `make_chem_funs.py` writes it again.
- `make_chem_funs.read_network` renumbers the reactions; reaction numbers in a network file are not stable identifiers.
- `Config` is an object passed explicitly (`vulcan_cfg`) through nearly every part; a new runtime setting is added there first. Each run writes the configuration to `output_dir/vulcan_cfg.txt` (`Config.write_file`).
- The default initial composition is `ini_mix = 'eq'`, which needs the `fastchem_vulcan/fastchem` executable.
- Output is controlled by `save_type` (default `['pickle', 'netcdf', 'plain']`), `out_name` (default `example`), `output_dir` (default `output/`), `clean_output` and `save_evolution`. `output_humanread` is deprecated; use `save_type = 'plain'`.
- AGNI coupling goes through the file system: AGNI is expected in `VULCAN/AGNI/` with `agni.jl`, called from Python through `juliacall`.

## Agent files

This repository keeps `AGENTS.md` and `tests/AGENTS.md`. `CLAUDE.md` stays local and untracked (`.gitignore`), and `.github/copilot-instructions.md` is a link to it, as the maintainer set it up.
