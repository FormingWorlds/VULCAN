# VULCAN

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/VULCAN/tests.yaml?branch=main&label=Tests)](https://github.com/FormingWorlds/VULCAN/actions/workflows/tests.yaml)
[![Coverage](https://gist.githubusercontent.com/nichollsh/59f094e7d22cd6af9a1cb3ea665b4260/raw/covbadge.svg)](https://github.com/FormingWorlds/VULCAN/actions/workflows/tests.yaml)
[![Docs](https://img.shields.io/website?url=https%3A%2F%2Fproteus-framework.org%2FVULCAN%2F&label=Docs&up_message=online&down_message=offline)](https://proteus-framework.org/VULCAN/)

Photochemical and thermochemical kinetics for (exo)planetary atmospheres; a fast,
open-source 1D chemical kinetics code. This distribution of VULCAN contains a number of
performance and usability improvements, and is integrated into the
[PROTEUS framework](https://proteus-framework.org/PROTEUS/).

![Running with realtime plotting](docs/assets/demo.gif)

Full documentation: [proteus-framework.org/VULCAN](https://proteus-framework.org/VULCAN/). When contributing to this repository, please consult the [contributing guide](https://proteus-framework.org/VULCAN/Community/CONTRIBUTING.html).

Theory papers:
* [Tsai et al. 2021](https://arxiv.org/abs/2108.01790) (with photochemistry)
* [Tsai et al. 2017](https://arxiv.org/abs/1607.00409) (without photochemistry)

## Quick start

```sh
pip install -U -e .
cd fastchem_vulcan && make && cd ..
python run_vulcan.py
```

This installs VULCAN, compiles FastChem (used for equilibrium initialization), and runs
the default model (the hot Jupiter HD 189733b) with live plotting. It takes a few
minutes; results are written to `output/`.

Settings live in `src/vulcan/config.py`. For example, weaker vertical mixing:

```python
self.Kzz_prof  = 'const'
self.const_Kzz = 1.e7   # cm^2/s
```

or a carbon-rich (C/O = 1) atmosphere:

```python
self.use_solar = False
self.O_H = 5.37e-4
self.C_H = 5.37e-4
```

Set `self.use_live_plot = False` to skip real-time plotting. Pass `-n` to skip
regenerating the chemistry from the network (only valid if you haven't edited it).

## Plotting results

```sh
cd tools/
python plot_vulcan.py <output.pkl> <species> <plot_name> [-h]
```

Reads a pickled output (e.g. `output/example.pkl`) and plots mixing-ratio profiles for the
given comma-separated species (no spaces). Plots against pressure by default, or height
with `-h`.
