# Copilot instructions for VULCAN

## Build, test, and lint commands

Use the repo root (`VULCAN/`).

```bash
# install package
python -m pip install -e .

# install developer tooling used in CI/docs
python -m pip install -e '.[develop]'
python -m pip install -e '.[docs]'
```

```bash
# run model (default demo)
python vulcan.py

# skip chemistry regeneration (only safe when network is unchanged)
python vulcan.py -n

# compile FastChem executable used by eq initialisation
cd fastchem_vulcan && make
```

```bash
# tests
pytest

# single test
pytest tests/test_imports.py::test_package_imports -q

# keyword-filtered tests
pytest -k import

# coverage style used in CI
pytest --cov-report term --cov-report json:coverage.json --cov=vulcan tests/
```

```bash
# lint/format checks used by contributors
ruff check .
pre-commit run --files <changed_file_1> <changed_file_2>
```

```bash
# docs build
zensical build --clean
```

## High-level architecture

`vulcan.py` at the repository root is the CLI entrypoint; it prepends `src/` to `sys.path` and calls `vulcan.vulcan.run_cli()`.

`src/vulcan/vulcan.py` orchestrates a run in this order:
1. Create `Config` (`src/vulcan/config.py`), optionally regenerate chemistry (`make_chem_funs.make_all`) unless `-n`.
2. Build runtime containers: `store.Variables`, `store.AtmData`, `store.Parameters`.
3. Build atmospheric structure (`build_atm.Atm`): pressure interfaces, T-P-Kzz, molecular diffusion, boundary flux setup.
4. Read and construct chemistry (`op.ReadRate`), initialise abundances (`build_atm.InitialAbun`), then set derived atmospheric quantities.
5. Integrate with `op.Integration` + `op.Ros2` solver loop; optional photochemistry updates and optional AGNI coupling (`agni.py`) when `agni_call_frq > 0`.
6. Save outputs through `op.Output.save_out` (pickle by default) and plots in `output/plot/`.

Core module roles:
- `config.py`: all model/runtime settings (single source of configuration truth).
- `paths.py`: canonical absolute paths to repo data (thermo, cross-sections, FastChem, AGNI, generated chemistry file).
- `make_chem_funs.py`: generates `chem_funs.py` from the selected network and thermodynamic data.
- `op.py`: reaction-rate parsing, integration loop, ODE solver, IO/plot output.
- `build_atm.py`: atmosphere profile construction + initial abundance setup (including FastChem equilibrium option).
- `store.py`: in-memory state containers passed through the solver pipeline.

## Key conventions in this codebase

- Chemistry is generated code: edit network/thermo inputs, then rerun without `-n` so `make_chem_funs.py` regenerates `src/vulcan/chem_funs.py`.
- Network reaction IDs are rewritten automatically by `make_chem_funs.read_network`; do not treat numeric IDs in network files as stable identifiers.
- `Config` is object-based and passed explicitly (`vulcan_cfg`) through nearly every subsystem; add new runtime parameters there first.
- Runtime config is persisted each run (`output/vulcan_cfg.txt`) via `Config.write_file()`.
- Default initialisation mode is `ini_mix='eq'`, which requires the `fastchem_vulcan/fastchem` executable to exist.
- Output behavior is controlled by config flags (`clean_output`, `output_humanread`, `save_evolution`); default output file is `output/example.pkl`.
- AGNI coupling is filesystem-coupled: AGNI is expected at `VULCAN/AGNI/agni.jl` and activated through Julia via `juliacall`.
