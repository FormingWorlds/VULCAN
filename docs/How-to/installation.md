# Installation

This guide covers the installation of VULCAN, and running a quick verification.

!!! note "PROTEUS users"
    VULCAN is integrated into the [PROTEUS](https://proteus-framework.org) ecosystem. The standalone instructions below set VULCAN up on its own; if you are running it as part of PROTEUS, follow the [PROTEUS installation guide](https://proteus-framework.org/PROTEUS/How-to/installation.html) instead, which configures VULCAN alongside the other modules.

!!! info "Prerequisites"
    - **Python** ≥ 3.11 (3.12 and 3.13 also supported)
    - **C/C++ compiler** and **make**: required to build FastChem
    - **Git**: to clone the repository
    - **Internet access**: once, to clone the repo and pull Python dependencies
    - *(Optional)* **Conda** or **venv**: recommended to isolate the installation
    - *(Optional)* **Julia** ≥ 1.9: only needed for self-consistent climate (RCE) calculations with AGNI

---

## 1. Set up an environment (recommended)

Isolating the install keeps VULCAN's dependencies separate from your system Python. We
recommend [conda-forge](https://conda-forge.org/).

=== "Conda"
    ```sh
    conda create -n vulcan python=3.12 -y
    conda activate vulcan
    ```
=== "venv"
    ```sh
    python -m venv .venv
    source .venv/bin/activate   # Windows: .venv\Scripts\activate
    ```

## 2. Clone the repository

```sh
git clone https://github.com/FormingWorlds/VULCAN.git
cd VULCAN
```

## 3. Install VULCAN in editable mode

The Python dependencies are declared in
[`pyproject.toml`](https://github.com/FormingWorlds/VULCAN/blob/main/pyproject.toml) and are
installed automatically. From the `VULCAN/` directory:

```sh
pip install -U -e .
```

Editable mode keeps the source directory live, so any edits you make are picked up without
reinstalling.

## 4. Compile FastChem

FastChem is the equilibrium-chemistry code VULCAN uses to initialize compositions when
`ini_mix = 'eq'` (the default). Build it before the first run:

```sh
cd fastchem_vulcan
make
cd ..
```

On success you will see:

```
compiling done.
linking done.
everything is done and fine. enjoy your day!
```

This produces the `fastchem` executable **inside the `fastchem_vulcan/` directory** (object
files are placed in `obj/`). VULCAN looks for the binary at `fastchem_vulcan/fastchem`; if it
is missing, `build_atm.py` raises a clear error telling you to run `make` there.

!!! tip "When is FastChem needed?"
    FastChem is only invoked for equilibrium initialization (`ini_mix = 'eq'`). If you
    initialize from constant abundances (`const_mix`, `const_lowt`), a previous run
    (`vulcan_ini`), or an ASCII table (`table`), FastChem is not required and this step can be
    skipped.

## 5. (Optional) Install AGNI for RCE modelling

To couple VULCAN to the Julia-based 1D radiative–convective model
[AGNI](https://www.h-nicholls.space/AGNI), which lets VULCAN update $T(p)$ self-consistently rather than using a fixed profile, install Julia:

```sh
curl -fsSL https://install.julialang.org | sh
```

Then install AGNI following [the online instructions](https://www.h-nicholls.space/AGNI), and
create a symbolic link from your AGNI install into the `VULCAN/` folder. Verify the link from
inside `VULCAN/`:

```sh
file AGNI/agni.jl
```

VULCAN talks to AGNI through `juliacall` and expects to find `AGNI/agni.jl` at that path
(`agni.py` checks for it). The coupling is **off by default** — enable it by setting
`agni_call_frq > 0` (and `solve_rce` as desired) in your configuration.

---

## Verifying the installation

Run the bundled demo from the `VULCAN/` directory:

```sh
python run_vulcan.py
```

This will:

- generate `chem_funs.py` from the default chemical network, and
- run the default model (the HD 189733b setup defined in `config.py`).

If it succeeds, results are written to the `output/` folder (the pickled output
`output/example.pkl` and diagnostic plots under `output/plot/`).

!!! tip "Skip network regeneration"
    Building `chem_funs.py` is only necessary when the network changes. Pass `-n` to reuse the existing file and start faster:
    ```sh
    python run_vulcan.py -n
    ```

!!! tip "Threading"
    VULCAN reads `OMP_NUM_THREADS` (default 2) and propagates it to the BLAS/threading backends. Set it before launching to control the linear-algebra thread count, e.g.:
    ```sh
    OMP_NUM_THREADS=4 python run_vulcan.py
    ```

---

## Testing

Install the developer extras and run the test suite:

```sh
pip install -e ".[develop]"
pytest
```

If you run into any issues, please open an issue on
[GitHub](https://github.com/FormingWorlds/VULCAN/issues).