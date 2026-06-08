# Editing or using a different chemical network

VULCAN's chemistry is **not hard-coded**: the network is read from a plain-text file, and the source term and Jacobian are generated from it at runtime (`make_chem_funs.py` writes
`chem_funs.py`). You can switch to a different shipped network or edit one yourself.

!!! note "How the network becomes code"
    On each run (unless you pass `-n`), VULCAN parses the network file and regenerates `chem_funs.py`: the production/loss term `chemdf`, the `Gibbs` reversal function, and the analytic Jacobians. It then checks the network for element conservation and duplicate reactions before integrating.

---

## Switching networks

The active network and the elements it spans are both set in the `Config` class
(`config.py`). The two **must be consistent**:

```python
self.network   = VULCAN_DIR + 'thermo/CHO_photo_network.txt'   # C, H, O only
self.atom_list = ['H', 'O', 'C']
```

For a richer network, update both, and the corresponding elemental abundances:

```python
self.network   = VULCAN_DIR + 'thermo/NCHO_photo_network.txt'  # adds nitrogen
self.atom_list = ['H', 'O', 'C', 'N']
# make sure N_H (and S_H for an S network, etc.) are set
```

Then run **without** `-n`:

```sh
python run_vulcan.py
```

!!! warning "Keep `atom_list` in sync"
    Initialization, element-conservation checking, and the FastChem equilibrium setup all loop over `atom_list`. If it does not match the elements present in the network (and in the composition file), the run will fail the conservation check or misinitialize abundances.

---

## Network file format

Each two-body thermochemical reaction is one line:

```text
[ Reactant1 + Reactant2 -> Product1 + Product2 ]  A  B  C
```

where `A B C` are the modified-Arrhenius coefficients (called `a n E` in the code), with

$$k = A\,T^{B}\exp(-C/T)$$

and `C` (the activation term) expressed in kelvin. Example:

```text
[ OH + H2 -> H2O + H ]  3.57E-16   1.520   1740.0
[ O  + H2 -> OH  + H ]  8.52E-20   2.670   3160.0
```

!!! note "Reactions are auto-numbered"
    VULCAN **rewrites the network file in place** on each run, prefixing every reaction with a managed integer index (`read_network` in `make_chem_funs.py`). Any leading number you type is overwritten, so just write the `[ ... ] A B C` line and let the code handle indexing.

### Reaction sections

Reaction types are grouped under comment headers, which the parser uses to switch modes. The
**column format differs by section**:

| Section header | Contents | Numeric columns |
|---|---|---|
| (top, default) | two-body thermochemical reactions | `A B C` (3) |
| `# 3-body` | three-body reactions with a high-pressure limit | `A B C  A∞ B∞ C∞` (6) |
| `# 3-body reactions without high-pressure rates` | three-body, low-pressure only | `A B C` (3) |
| `# special` | hard-coded rate forms (e.g. `OH + CH3 + M -> CH3OH + M`) | — |
| `# condensation` | condensation/evaporation pairs | (set internally) |
| `# radiative` | radiative recombination | `A B C` |
| `# photo` | photodissociation branches | branch index (not Arrhenius) |
| `# ionisation` | photoionization branches | branch index |

---

## Editing the default network

1. Open the network file:
    ```sh
    nano thermo/NCHO_photo_network.txt
    ```
2. **Add a reaction** — insert a line in the appropriate section (two-body vs. three-body vs.
   photo), with the right number of coefficient columns for that section:
    ```text
    [ Species1 + Species2 -> Product1 + Product2 ]  A  B  C
    ```
3. **Remove a reaction** — delete its line.
4. **Change rate coefficients** — edit the `A B C` (and `A∞ B∞ C∞` for three-body) values.
5. **Save and run** (without `-n`):
    ```sh
    python run_vulcan.py
    ```
    VULCAN regenerates `chem_funs.py` with your changes.

!!! info "Forward reactions only"
    List only **forward** reactions up to the `# reverse stops` marker. VULCAN computes each reverse rate from the equilibrium constant (`k_reverse = k_forward / K_eq(T)`), so adding a reverse reaction by hand would double-count it. Photodissociation and condensation reactions are not reversed this way.

---

## Adding new species

A new species needs thermodynamic data (for the equilibrium/reversal machinery) and an entry
in the composition file.

### 1. NASA-9 thermodynamic data

Create `thermo/NASA9/<SPECIES>.txt` (the filename must match the species name used in the
network). VULCAN loads it with `np.loadtxt`, flattens it, and reads coefficients **`[0:10]` as the low-temperature range and `[10:20]` as the high-temperature range** (`make_Gibbs`), so the file must contain **both** ranges; 20 values total. Each range is laid out as:

```text
a1 a2 a3 a4 a5
a6 a7 0. a8 a9
```

!!! note
    Mind the `0.` placeholder between `a7` and `a8` in each range. You can copy the layout from an existing file or from the bulk tables `thermo/NASA9/nasa9_2002_E.txt` / `new_nasa9.txt`.

### 2. Composition (element) entry

Add the species to the composition file (`COM_FILE`, by default `thermo/all_compose.txt`),
giving its per-element atom counts, charge (`e`), and molar mass. Match the existing column
order exactly:

```text
species         H   O   C   He  N   S   P   Na  K   Si  Fe  Ar  Ti  V   Mg  Ca  e   mass
H2O             2   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   18.016
```

This file drives the element-conservation check, so an incorrect atom count will be caught at startup.

### 3. (Photodissociating species only) cross sections

A species that appears in the `# photo` section also needs photoabsorption / dissociation
cross sections and branching ratios under `thermo/photo_cross/<SPECIES>/` (read by
`op.ReadRate.make_bins_read_cross`). Without them the run will stop with a "missing cross
section" error.

### 4. Run

```sh
python run_vulcan.py
```

---

!!! warning "Editing networks requires understanding the chemistry"
    Reverse rates and element balance depend on the network being physically consistent. Start from an existing network, make small changes, and confirm the element-conservation and duplicate checks pass before trusting the results.