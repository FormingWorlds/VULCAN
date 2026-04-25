This folder contains photochemical network and thermodynamic files.
photo_cross has the UV cross sections for photochemistry.
NASA9 has the Gibbs free energy in the form of NASA-9 polynomials (http://garfield.chem.elte.hu/Burcat/burcat.html).
all_compose.txt containes the basic physical data for all species and stoichiometric numbers.

## CRAHCN Networks (thermo/CRAHCN/)

The `CRAHCN/` subdirectory contains the Consistent Reduced Atmospheric Hybrid
Chemical Network and its oxygen extension, developed by Dr. Ben K. D. Pearce.
These are reduced networks for modelling prebiotic chemistry precursors in
planetary atmospheres at low temperatures (50–400 K).

Files:
- `CRAHCN.reac`       — 104-reaction HCN network for reducing atmospheres (H2/N2/CH4)
- `CRAHCN.spec`       — species definitions for CRAHCN
- `CRAHCN-O_v3.reac`  — extended network with oxygen species (H2/CO2/N2/CO/H2O/CH4)
- `CRAHCN-O.spec`     — species definitions for CRAHCN-O

**Format note:** CRAHCN files use the kinetic form `k = α (T/300)^β exp(−γ/T)`,
which differs from VULCAN's standard `k = A T^B exp(−C/T)` Arrhenius form.
The algebraic conversion is: `A = α / 300^β`, `B = β`, `C = γ`.
Three-body reactions in CRAHCN-O also include per-gas enhancement factors
for N2, CO2, and H2 colliders. CRAHCN(-O) also uses species names that differ
from VULCAN's existing network conventions for some excited states and radicals
(e.g. `N2D`, `O1D`, `O3P`, `N4S` here vs. the underscore-delimited forms in
`all_compose.txt` and shipped networks); species name mapping is a required
additional step when integrating these files. Full integration with VULCAN's
network parser requires format conversion as a separate follow-up.

**License:** MIT — originally from https://github.com/bennski/CRAHCN

**References:**
- Pearce et al. (2020), ApJ — CRAHCN applied to Titan
- Pearce et al. (2022), ApJ — CRAHCN-O applied to early Earth
- Pearce et al. (2022), ACS Earth Space Chem — CRAHCN-O plasma experiment
