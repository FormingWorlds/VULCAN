# API overview

This is an overview of VULCAN's API for reference while developing or extending the code.
If you want to understand the underlying model rather than the code structure, see the
[model overview](../../Explanations/model_overview.md).
## Source tree

```text
src/vulcan
├── __init__.py        # package entry point: exposes Config, run_vulcan, make_chem_funs, paths, __version__
├── vulcan.py           # CLI driver: run_vulcan / run_cli (the `-n` flag skips network regeneration)
├── config.py           # Config class: every run parameter (see Configuration parameters)
├── store.py            # Variables, AtmData, Parameters: the runtime state containers
├── build_atm.py        # InitialAbun, Atm: atmosphere setup, composition, diffusion coefficients
├── make_chem_funs.py   # parses the network file and generates chem_funs.py
├── chem_funs.py         # generated: chemdf, Gibbs, symjac/neg_symjac, spec_list, ni, nr
├── op.py                # ReadRate, Integration, ODESolver, Ros2, Output: the numerical core
├── agni.py              # optional coupling to the AGNI radiative-convective model
├── phy_const.py         # physical constants (cgs)
├── paths.py             # VULCAN_DIR, THERMO_DIR, CROSS_DIR and other path resolution
└── logs.py              # logging setup
```

!!! warning "`chem_funs.py` is generated, not version-controlled source"
    `chem_funs.py` does not contain hand-written code. It is (re)written by
    `make_chem_funs.make_all` from the active reaction network every time the network
    changes (i.e. whenever you run without `-n`). If you are reading the API reference for
    this module, keep in mind its contents depend on `config.network`.

## Modules by role

### Running a model

| Module | Key entry points | Reference |
|---|---|---|
| `vulcan` | `run_vulcan`, `run_cli` | [API](../Reference/api/vulcan.md) |
| `config` | `Config` | [API](../Reference/api/config.md) · [Parameter reference](../Reference/parameters.md) |
| `make_chem_funs` | `make_all` (and the lower-level `make_chemdf`, `make_jac`, `make_Gibbs`, `check_conserv`, `check_duplicate`) | [API](../Reference/api/make_chem_funs.md) |

### State and atmosphere

| Module | Key classes | Reference |
|---|---|---|
| `store` | `Variables`, `AtmData`, `Parameters` | [API](../Reference/api/store.md) |
| `build_atm` | `InitialAbun`, `Atm` | [API](../Reference/api/build_atm.md) |
| `chem_funs` *(generated)* | `chemdf`, `Gibbs`, `symjac`, `neg_symjac`, `spec_list` | [API](../Reference/api/chem_funs.md) |

### Numerical core

| Module | Key classes | Reference |
|---|---|---|
| `op` | `ReadRate`, `Integration`, `ODESolver`, `Ros2`, `Output` | [API](../Reference/api/op.md) |
| `agni` | `run_agni` (optional, requires `agni_call_frq > 0`) | [API](../Reference/api/agni.md) |

### Utilities

| Module | Contents | Reference |
|---|---|---|
| `phy_const` | Physical constants in cgs units (`kb`, `Navo`, `hc`, …) | [API](../Reference/api/phy_const.md) |
| `paths` | `VULCAN_DIR`, `THERMO_DIR`, `CROSS_DIR` and related path constants | [API](../Reference/api/paths.md) |
| `logs` | Logger configuration (`fwl.*` loggers) | [API](../Reference/api/logs.md) |
