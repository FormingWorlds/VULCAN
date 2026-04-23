# Installation

This guide covers the installation and setup of VULCAN, a photochemical kinetics code for exoplanetary atmospheres.

## Prerequisites

Before installing VULCAN, ensure you have the following:

- **Python 3.11 or higher** (Python 3.12 and 3.13 are also supported)
- **C++ compiler** (for compiling FastChem)
- **Git** (for cloning the repository)
- **Basic build tools** (make, cmake, etc.)
- **Julia** (optional, for running self-consistent climate calculations)

## Installation steps

### 0. Install Python

You must install Python. We recommend using [conda forge](https://conda-forge.org/).


### 1. Clone the VULCAN repository

```bash
git clone https://github.com/FormingWorlds/VULCAN.git
cd VULCAN
```

### 2. Install VULCAN and its dependencies

VULCAN requires Python packages specified in `pyproject.toml`. These will be installed automatically.

To install the VULCAN package, run the following command in your `VULCAN/` folder.

```bash
pip install -U -e .
```

### 3. Compile FastChem

FastChem is an equilibrium chemistry code used to initialize compositions. It must be compiled before running VULCAN.

```bash
cd fastchem_vulcan
make
cd ..
```

The compilation process:
- Compiles the FastChem source code in `fastchem_src/`
- Builds the model main executable
- Links everything together
- Creates the necessary executables in the `obj/` folder

If compilation succeeds, you'll see:
```
compiling done.
linking done.
everything is done and fine. enjoy your day!
```

### 4. (Optional for RCE modelling) Install AGNI

If you plan on the coupling to the julia-based 1D radiative-convective model AGNI, julia needs to be installed.

To do this, run:

```bash
curl -fsSL https://install.julialang.org | sh
```

!!! warning "Pin Julia to version 1.11"
    Julia 1.12+ is **not yet supported** due to OpenSSL library incompatibilities with Python. After installing Julia, pin it to version 1.11:

    ```console
    juliaup add 1.11
    juliaup default 1.11
    ```

Then, install AGNI following [the instructions online](https://www.h-nicholls.space/AGNI).

Finally, create a symbolic link from your AGNI install folder to your VULCAN folder. Test this by running `file AGNI/agni.jl` inside your local `VULCAN/` folder.

### 5. Verify installation

Test that everything is working by running the demo:

```bash
python vulcan.py
```

This will:
- Generate `chem_funs.py` based on the default chemical network
- Run the default model

If successful, you'll see output files in the `output/` folder.

## Testing

Install test dependencies and run the pytest suite:

```bash
pip install -e ".[develop]"
pytest
```
