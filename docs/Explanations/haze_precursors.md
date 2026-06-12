# Photochemical haze precursors

Clouds and photochemical hazes are ubiquitous across planetary atmospheres, but their
microphysical formation is complex and uncertain. Rather than model particle growth directly,
VULCAN tracks a set of gas-phase **precursor species** as proxies for haze formation
(Section 2.8 of Tsai et al. 2021 [^tsai2021]).

## Precursor selection

The model preferentially considers precursors related to polycyclic aromatic hydrocarbon
(PAH) or nitrile formation. The full precursor set is:

$$\mathrm{C_2H_2,\ C_2H_6,\ C_4H_2,\ C_6H_6,\ HCN,\ HC_3N,\ CH_2NH,\ CH_3CN,\ CS_2}$$

- **Benzene (C$_6$H$_6$)** is treated as the key PAH proxy: once the first aromatic ring forms,
  the thermodynamics of attaching further rings changes little, and benzene formation is
  argued to be the rate-limiting step in building complex hydrocarbons.
- **Nitriles** are represented by HCN together with the less abundant H$_2$CN, CH$_2$NH,
  CH$_3$CN, and HC$_3$N. HCN itself is rarely the limiting factor, so the rarer nitriles are
  more indicative of complex-nitrile formation.
- **CS$_2$** represents sulfur-bearing haze precursors, following laboratory work
  (He et al. 2020).

## Benzene mechanism

The simplified benzene-forming pathway uses propargyl recombination,

$$\mathrm{C_3H_3 + C_3H_3 \xrightarrow{M} C_6H_6}$$

with $\mathrm{C_3H_3}$ produced from $\mathrm{CH_3 + C_2H \rightarrow C_3H_3 + H}$
(Frenklach 2002 [^frenklach2002]). The intention is to capture the main formation pathway at
minimum network cost; benzene photodissociation branches are poorly constrained, so the
dominant channel is assumed to yield phenyl (C$_6$H$_5$) with a small (~15%) fraction to
C$_3$H$_3$ (Kislov et al. 2004 [^kislov2004]). See [chemical networks](chemical_networks.md)
for the supporting species.

## Diagnostic: column density above 1 mbar

The model does not grow particles; instead the **column number density of each precursor
above 1 mbar** is used as the diagnostic of haze-forming potential (Figure 37 of the paper [^tsai2021]).
Across the irradiated H$_2$-dominated atmospheres studied, HCN is consistently the most
abundant precursor, but this does not by itself imply complex-nitrile formation. HC$_3$N,
C$_4$H$_2$, and C$_6$H$_6$ tend to increase with decreasing temperature, while CH$_3$CN shows
the opposite trend; CS$_2$ (containing no hydrogen) is most favored in hot Jupiter conditions.

## Implementation in VULCAN

As with the [Ti/V chemistry](titanium_vanadium.md), haze precursors require no dedicated code
path: the precursor species and the benzene mechanism are part of the network file and are
integrated like any other species through `chemdf` and the photolysis routines. The benzene
photodissociation branches and cross sections (Boechat-Roberty et al. 2004; Capalbo et al.
2016) are read by [`op.ReadRate.make_bins_read_cross`](../Reference/api/op.md#vulcan.op.ReadRate.make_bins_read_cross) and applied in `compute_J`.


!!! note "Caveat"
    Benzene's photodissociation branches are poorly constrained across their various products, and the predicted abundances of  C$_4$H$_2$ and
    C$_6$H$_6$ are not considered accurate. They serve to assess relative haze-precursor
    trends rather than to predict haze mass.

[^tsai2021]: Tsai, S.-M., Malik, M., Kitzmann, D., et al. (2021). A comparative study of atmospheric chemistry with VULCAN. *The Astrophysical Journal, 923*(2), 264. https://doi.org/10.3847/1538-4357/ac29bc

[^frenklach2002]: Frenklach, M. (2002). Reaction mechanism of soot formation in flames. *Physical Chemistry Chemical Physics, 4*(11), 2028–2037. https://doi.org/10.1039/B110045A

[^kislov2004]: Kislov, V. V., Nguyen, T. L., Mebel, A. M., Lin, S. H., & Smith, S. C. (2004). Photodissociation of benzene. *The Journal of Chemical Physics, 120*(15), 7008. https://doi.org/10.1063/1.1676275