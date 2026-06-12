# Chemical networks

VULCAN's chemistry is defined by editable network files listing forward reactions and their
rate coefficients. The code reverses every thermochemical reaction internally so that
chemical equilibrium is recovered kinetically, and it generates the source term and Jacobian
from the network at runtime (Section 2.3 of Tsai et al. 2021 [^tsai2021]).

## Network hierarchy

Networks are provided hierarchically (C–H–O, C–H–N–O, C–H–N–O–S), each in a **reduced** and
a **full** version. "Reduced" refers both to oxidation state and to network size: the reduced
version strips out species and mechanisms (e.g. the ozone cycle) that matter only in
oxidizing conditions, and is more efficient for hydrogen-dominated atmospheres. The full
version spans reducing to oxidizing conditions.

The full C–H–N–O–S network contains **96 species**, about **570 forward thermochemical
reactions**, and **69 photodissociation branches**. Hydrocarbons are truncated at two
carbons, except for a few higher-order species retained as sinks or haze precursors.

The active network is selected with [`config.network`](../Reference/config.md#chemistry-network) (default
`thermo/CHO_photo_network.txt`) and `config.atom_list`. The reaction files live under
`atm/` in the repository.

## Rate coefficients

Forward rate coefficients follow the generalised Arrhenius equation:

$$
k = a\,T^{n}\exp(-E/T), \tag{7}
$$ 

where $k$ is the rate coefficient in cm$^3$s$^{-1}$ for bimolecular reactions and  cm$^6$s$^{-1}$ for termolecular reactions. Three-body reactions additionally carry high-pressure-limit
parameters $(a_\infty, n_\infty, E_\infty)$ and are combined through a Lindemann-type
falloff. Rate data are drawn from the NIST and KIDA databases and the combustion/planetary
literature; some coefficients are only measured over limited temperature ranges, which is a
known source of uncertainty (especially for sulfur).

### Implementation in VULCAN

[`op.ReadRate.read_rate`](../Reference/api/op.mdl#vulcan.op.ReadRate.read_rate) parses the network file block by block, recognizing the section
headers `# 3-body`, `# 3-body reactions without high-pressure rates`, `# special`,
`# condensation`, `# radiative`, `# photo`, and `# ionisation`. For each reaction it stores
$a, n, E$ (and $a_\infty, n_\infty, E_\infty$ for three-body reactions) and builds a callable
`k_fun[i]`. The three-body rate is assembled as
`k = k0 / (1 + k0·M / k_inf)` with $M$ the third-body number density ($M = p/k_B T$).

A few reactions use **hard-coded special forms**: the methanol-forming reaction
`OH + CH3 + M -> CH3OH + M` uses the Jasper (2007) [^jasper2007] rate coefficients because of its thorough treatment of the reaction. For the reaction `CH3OH + H -> CH3 + H2O` the rate coefficient by Moses et al. (2011) [^moses2011] is used. 

## Reversing reactions

All thermochemical reactions are reversed using the equilibrium constant derived from the
NASA polynomials. The forward/reverse reaction rate coefficients follow the following relation: 


$$\frac{k_\text{f}}{k_\text{r}} = K_\mathrm{eq} \left(\frac{k_\text{B} T}{P_0}\right)^{\Delta\mu} \tag{8}$$

where $K_\mathrm{eq}$ is the dimensional equilibrium constant, expressed via the standard Gibbs free energy as

$$K_\mathrm{eq} = \exp\!\left(-\frac{\Delta G^0}{\mathcal{R}T}\right) = \exp\!\left(-\frac{\Delta H^0 - T\Delta s^0}{\mathcal{R}T}\right) \tag{9}$$

with $\Delta G^0 = \Delta G[\text{products}] - \Delta G[\text{reactants}]$, so the reverse rate coefficient is

$$k_\text{r} = \frac{k_\text{f}}{\exp\!\left[-(\Delta H^0 - T\Delta s^0)/\mathcal{R}T\right]} \left(\frac{k_\text{B} T}{P_0}\right)^{-\Delta\mu} \tag{10}$$

The standard enthalpy and entropy from the NASA polynomials are

$$\frac{H^0}{\mathcal{R}T} = -a_1 T^{-2} + a_2 \frac{\ln T}{T} + a_3 + \frac{a_4}{2}T + \frac{a_5}{3}T^2 + \frac{a_6}{4}T^3 + \frac{a_7}{5}T^4 + \frac{a_8}{T} \tag{11}$$

### Implementation in VULCAN

[`make_chem_funs.make_Gibbs`](../Reference/api/make_chem_funs.md#vulcan.make_chem_funs.make_Gibbs) writes a `Gibbs(i, T)` function into `chem_funs.py` that
evaluates $K_\mathrm{eq}$ from the per-species Gibbs free energy. NASA-9 polynomials are loaded
from `thermo/NASA9/`. [`op.ReadRate.rev_rate`](../Reference/api/op.md#l#vulcan.op.ReadRate.rev_rate) then sets `k[i] = k[i-1] / Gibbs(i-1, Tco)` for
every reversible reaction up to the `# reverse stops` marker. Reactions in
`config.remove_list` are zeroed by `remove_rate`.

!!! warning "CH$_2$NH discrepancy"
    There is a significant discrepancy in the newer NASA-9 polynomial for CH$_2$NH
    relative to the NASA-7 fit, which can cause errors of several orders of magnitude; the
    NASA-7 fit is used for that species instead.

## Code generation

Rather than evaluating the network interpretively, VULCAN **writes out Python source** for
the chemistry. [`make_chem_funs.make_chemdf`](../Reference/api/make_chem_funs.md#vulcan.make_chem_funs.make_chemdf) parses the reaction table and emits
`chem_funs.py` containing:

- `chemdf(y, M, k)`: the vectorized production-minus-loss term $P_\text{i} - L_\text{i}$;
- `df(y, M, k)`: the explicit per-reaction expressions used to build the Jacobian;
- `re_dict` / `re_wM_dict`: reactant/product maps (with and without $M$);
- `spec_list`, `ni` (number of species), `nr` (number of reactions).

`make_jac` and `make_neg_jac` use SymPy to differentiate `df` symbolically and write the
analytic Jacobians `symjac` and `neg_symjac`. The driver `make_all` (in `make_chem_funs.py`)
runs this pipeline and then verifies the network with `check_conserv` (element conservation
per reaction) and `check_duplicate` (no repeated reactions).

## Benzene and modular networks

A simplified benzene mechanism is embedded for haze-precursor studies (see
[haze precursors](haze_precursors.md)): benzene forms through propargyl recombination
$\mathrm{C_3H_3 + C_3H_3 \xrightarrow{M} C_6H_6}$ (Frenklach 2002 [^frenklach2002]), with
$\mathrm{C_3H_3}$ produced from $\mathrm{CH_3 + C_2H}$ and supporting species
(C$_3$H$_2$, C$_3$H$_4$, C$_6$H$_5$) added for the hydrogen-abstraction steps.

A user can also build a **modular** network by picking a subset of species. Only reactions
involving the selected species are retained. Unlike Gibbs-minimization equilibrium codes,
care is needed to keep important trace intermediates in the subset.

[^tsai2021]: Tsai, S.-M., Malik, M., Kitzmann, D., et al. (2021). A comparative study of atmospheric chemistry with VULCAN. *The Astrophysical Journal, 923*(2), 264. https://doi.org/10.3847/1538-4357/ac29bc

[^frenklach2002]: Frenklach, M. (2002). Reaction mechanism of soot formation in flames. *Physical Chemistry Chemical Physics, 4*(11), 2028–2037. https://doi.org/10.1039/B110045A

[^jasper2007]: Jasper, A. W., Klippenstein, S. J., Harding, L. B., & Ruscic, B. (2007). Kinetics of the reaction of methyl radical with hydroxyl radical and methanol decomposition. *The Journal of Physical Chemistry A, 111*(19), 3932-3950. https://doi.org/10.1021/jp067585p

[^moses2011]: Moses, J. I., Visscher, C., Fortney, J. J., Showman, A. P., Lewis, N. K., Griffith, C. A., ... & Freedman, R. S. (2011). Disequilibrium carbon, oxygen, and nitrogen chemistry in the atmospheres of HD 189733b and HD 209458b. *The Astrophysical Journal, 737*(1), 15. https://doi.org/10.1088/0004-637X/737/1/15