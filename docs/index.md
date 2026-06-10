---
title: VULCAN
---

<h1 align="center" style="margin-bottom:0">VULCAN</h1>
<h3 align="center" style="margin-top:.5em ; color:#6c6c6c">Photochemical &amp; thermochemical kinetics for (exo)planetary atmospheres</h3>

<p align="center">
  <a href="https://github.com/FormingWorlds/VULCAN/actions/workflows/tests.yaml">
    <img src="https://img.shields.io/github/actions/workflow/status/FormingWorlds/VULCAN/tests.yaml?branch=main&label=Tests">
  </a>
  <a href="https://github.com/FormingWorlds/VULCAN/actions/workflows/tests.yaml">
    <img src="https://gist.githubusercontent.com/nichollsh/59f094e7d22cd6af9a1cb3ea665b4260/raw/covbadge.svg">
  </a>
  <a href="https://proteus-framework.org/VULCAN/">
    <img src="https://img.shields.io/website?url=https%3A%2F%2Fproteus-framework.org%2FVULCAN%2F&label=Docs&up_message=online&down_message=offline">
  </a>
  <a href="https://www.gnu.org/licenses/gpl-3.0">
    <img src="https://img.shields.io/badge/License-GPLv3-blue.svg">
  </a>
</p>

**VULCAN** is an open-source 1D chemical kinetics code for planetary and exoplanetary
atmospheres. It solves the coupled vertical continuity equations for every species, combining
thermochemistry, photochemistry, and vertical transport to predict atmospheric composition far
from chemical equilibrium. VULCAN can be used on its own, and is also integrated into the
[PROTEUS framework](https://proteus-framework.org/PROTEUS/), a coupled simulation tool for the
long-term evolution of rocky-planet atmospheres and interiors.

<p align="center">
  <img src="assets/demo.gif" alt="VULCAN running with real-time plotting" width="50%"/>
</p>

<p style="text-align: center;"><strong>A VULCAN run converging to steady state.</strong></p>

## Key features

- **Editable  networks**: chemistry is read from plain-text reaction files
  (hierarchical C–H–O, C–H–N–O, C–H–N–O–S); VULCAN
  auto-generates the source term and analytic Jacobian from the network at runtime.
- **Thermochemistry *and* photochemistry**: forward reactions are reversed from NASA-9
  thermodynamic data to recover equilibrium, while a two-stream radiative-transfer scheme
  drives photodissociation with temperature-dependent UV cross sections.
- **Full vertical transport**: advection, eddy diffusion, and molecular plus thermal
  diffusion, with flexible boundary conditions.
- **Condensation and settling**: saturation-based condensation/evaporation with gravitational particle settling for species such as H$_2$O, NH$_3$, and sulfur.
- **Equilibrium initialization**: compositions are initialized at chemical equilibrium with
  the embedded [FastChem](https://github.com/NewStrangeWorlds/FastChem) code.
- **Validated across regimes**: benchmarked against hot Jupiters, Jupiter, and modern Earth,
  from reducing to oxidizing atmospheres.
- **Optional self-consistent climate**: couples to the radiative–convective model
  [AGNI](https://www.h-nicholls.space/AGNI) to update the temperature structure during a run


## Get started

<div class="grid cards" markdown>

-   :material-download: **Install VULCAN**

    Set up VULCAN and compile FastChem, with optional Julia/AGNI for climate coupling.

    [Installation guide](How-to/installation.md){ .md-button .md-button--primary }

-   :material-rocket-launch: **Run your first model**

    Run the default model, edit the configuration, and plot the resulting profiles.

    [First run](Tutorials/first_run.md){ .md-button .md-button--primary }

</div>

## Citation and license
 
If you use VULCAN through PROTEUS, please cite the VULCAN theory papers listed in the [publications](Reference/bibliography.md) and state the code version used. Citation metadata is provided in [`CITATION.cff`](https://github.com/FormingWorlds/VULCAN/blob/main/CITATION.cff).
 
VULCAN is released under the [GPL-3.0 license](https://github.com/FormingWorlds/VULCAN/blob/main/LICENSE.txt).
 
### Embedded dependency: FastChem
 
VULCAN embeds [FastChem](https://github.com/NewStrangeWorlds/FastChem) (Daniel Kitzmann & Joachim Stock) for equilibrium-chemistry initialisation of mixing ratios. FastChem carries its own license and should be cited independently. The relevant papers are listed in the [publications](Reference/bibliography.md).
 
!!! info "Licenses across the PROTEUS framework"
    Different components within the PROTEUS framework carry different licenses. Please find information about the use of licenses within the PROTEUS framework on the PROTEUS [licence page](https://proteus-framework.org/license).
