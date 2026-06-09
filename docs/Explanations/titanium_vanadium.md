# Chemistry of Ti and V compounds

Titanium and vanadium oxides (TiO, VO) are strong optical absorbers thought to cause
thermal inversions in highly irradiated atmospheres. VULCAN includes a first kinetic
treatment of titanium and vanadium species, used in the WASP-33b case study
(Section 2.7 of Tsai et al. 2021 [^tsai2021]).

## Species and thermodynamic data

The species list is expanded with **Ti, TiO, TiO$_2$, TiH, TiC, TiN, V, and VO**. Of these,
only Ti, TiO, and TiO$_2$ have NASA-polynomial thermodynamic data; the remainder are sourced
as follows:

| Species | Thermodynamic source |
|---|---|
| Ti, TiO, TiO$_2$ | NASA polynomials |
| TiH | Burrows et al. (2005) [^burrows2005], from ab-initio Gibbs free energy |
| TiC | Woitke et al. (2018) [^woitke2018] |
| TiN, V, VO | Tsuji (1973) [^tsuji1973] |

## Rate-coefficient estimates

High-temperature kinetics data for Ti/V species are essentially nonexistent, so the rates are
**estimated**:

1. Where an analogous transition-metal (e.g. Fe) reaction is measured at high temperature, the
   same rate coefficient is adopted.
2. Otherwise, the temperature dependence is estimated from transition-state theory: for an
   endothermic reaction the activation energy (the exponential term in the Arrhenius form) is
   approximated by the **enthalpy difference** between products and reactants, assuming the
   transition-state energy increase is small relative to the enthalpy difference for radical
   reactions.
3. The pre-exponential factor is then adjusted to match the reference value at low temperature.

The adopted Ti/V reactions and rates are tabulated in Table B1 of the paper [^tsai2021]. For photolysis,
TiO, TiO$_2$, TiH, TiC, and VO are included, with UV cross sections estimated from FeO
(Chestakov et al. 2005 [^chestakov2005]) at 252.39 nm and the threshold scaled by each species' bond
dissociation energy.

## Implementation in VULCAN

Unlike the transport, photochemistry, and condensation routines, the Ti/V chemistry is **not
hard-coded in the Python modules**. It enters entirely through the data layer that the general
machinery already supports:

- the Ti/V reactions (Table B1) live in the network file and are parsed by
  [`op.ReadRate.read_rate`](../Reference/api/op.md#vulcan.op.ReadRate.read_rate) like any other reaction, including the modified-Arrhenius and
  three-body forms;
- the thermodynamic data sit in `thermo/NASA9/` (with the externally sourced fits for TiH,
  TiC, etc.), and are reversed to enforce equilibrium through the same
  [`make_chem_funs.make_Gibbs`](../Reference/api/make_chem_funs.md#vulcan.make_chem_funs.make_Gibbs) / [`op.ReadRate.rev_rate`]((../Reference/api/op.md#vulcan.op.ReadRate.rev_rate)) path described in
  [Chemical networks](chemical_networks.md);
- the Ti/V photodissociation branches are handled by `compute_J` using their tabulated cross
  sections, exactly as for the C–H–N–O–S photolysis reactions.

In other words, modeling Ti/V is a matter of selecting a network and thermo dataset that
include these species; no dedicated code path is required.

[^tsai2021]: Tsai, S.-M., Malik, M., Kitzmann, D., et al. (2021). A comparative study of atmospheric chemistry with VULCAN. *The Astrophysical Journal, 923*(2), 264. https://doi.org/10.3847/1538-4357/ac29bc

[^burrows2005]: Burrows, A., Dulick, M., Bauschlicher, C. W., et al. (2005). Spectroscopic constants, abundances, and opacities of the TiH molecule. *The Astrophysical Journal, 624*(2), 988. https://doi.org/10.1086/429366

[^woitke2018]: Woitke, P., Helling, C., Hunter, G. H., et al. (2018). Equilibrium chemistry down to 100 K. *Astronomy & Astrophysics, 614*, A1. https://doi.org/10.1051/0004-6361/201732193

[^tsuji1973]: Tsuji, T. (1973). Molecular abundances in stellar atmospheres. II. *Astronomy & Astrophysics, 23*, 411.

[^chestakov2005]: Chestakov, D. A., Parker, D. H., & Baklanov, A. V. (2005). Iron monoxide photodissociation. *The Journal of chemical physics, 122*(8). https://doi.org/10.1063/1.1844271