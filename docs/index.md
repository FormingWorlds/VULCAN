# VULCAN
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Photochemical kinetics for exoplanetary atmospheres, a fast and easy-to-use python code. 

VULCAN can be used on its own, but is also integrated into the [PROTEUS framework](https://proteus-framework.org/PROTEUS/), a coupled simulation tool for the long-term evolution of atmospheres and interiors of rocky planets.

This distribution of VULCAN contains a number of performance and usability improvements.

* Demo with realtime plotting:

![Running with realtime plotting](assets/demo.gif)

## Project Structure

```
VULCAN/
├── vulcan.py                 # Main entry point
├── config.py                 # Configuration file (edit this!)
├── make_chem_funs.py         # Generates chemical functions
├── build_atm.py              # Atmospheric structure setup
├── op.py                     # Numerical operations & ODE solvers
├── store.py                  # Data storage & class definitions
├── phy_const.py              # Physical constants
├── agni.py                   # Additional utilities
├── paths.py                  # Path definitions
├── atm/                      # Atmospheric input files & stellar fluxes
├── fastchem_vulcan/          # FastChem equilibrium chemistry code
├── thermo/                   # Chemical networks & thermodynamic data
│   ├── NASA9/               # NASA-9 polynomial coefficients
│   ├── photo_cross/         # Photochemical cross sections
│   └── NCHO_photo_network.txt  # Default reaction network
├── cfg_examples/             # Example configuration files
├── tools/                    # Analysis and plotting utilities
├── output/                   # Output files (created at runtime)
└── plot/                     # Output plots (created at runtime)
```

Typically ```config.py``` is the only file you need to edit for each specific run. If you want to look inside or modify the code, `store.py` is where almost all classes and variables are declared.

## Publications

The theory papers of VULCAN can be found here: 

- [Tsai et al. 2021](https://arxiv.org/abs/2108.01790) (with photochemistry)
- [Tsai et al. 2017](https://arxiv.org/abs/1607.00409) (without photochemistry)

