# Installation

This guide covers the installation and setup of VULCAN, a photochemical kinetics code for exoplanetary atmospheres.

## Prerequisites

Before installing VULCAN, ensure you have the following:

- **Python 3.11 or higher** (Python 3.12 and 3.13 are also supported)
- **C++ compiler** (for compiling FastChem)
- **Git** (for cloning the repository)
- **Basic build tools** (make, cmake, etc.)

## Installation steps

### 0. Install julia

For the coupling to the julia-based 1D radiative-convective model AGNI, julia needs to be installed. Do not install julia via your computer's package manager. Instead:

```bash
curl -fsSL https://install.julialang.org | sh
```

!!! warning "Pin Julia to version 1.11"
    Julia 1.12+ is **not yet supported** due to OpenSSL library incompatibilities with Python. After installing Julia, pin it to version 1.11:

    ```console
    juliaup add 1.11
    juliaup default 1.11
    ```


### 1. Clone the VULCAN repository

```bash
git clone https://github.com/FormingWorlds/VULCAN.git
cd VULCAN
```

### 2. Install python dependencies

VULCAN requires the following Python packages (specified in `pyproject.toml`):

- numpy (≥2.0.0)
- scipy
- matplotlib
- pandas
- sympy
- astropy
- juliacall

#### Option A: Install with pip (Recommended)

```bash
pip install -e .
```

This installs VULCAN in development mode along with all dependencies.

#### Option B: Install dependencies manually

```bash
pip install numpy scipy matplotlib pandas sympy astropy juliacall
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

### 4. Verify installation

Test that everything is working by running the demo:

```bash
python vulcan.py
```

This will:
- Generate `chem_funs.py` based on the default chemical network
- Run the default model (Earth-like atmosphere)
- Display real-time plotting (if `use_live_plot = True` in `config.py`)
- Take approximately 10-15 minutes to complete

If successful, you'll see output files in the `/output` folder.