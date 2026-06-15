# Editing or using a different chemical network

VULCAN's chemistry is **not hard-coded**: the network is read from a plain-text file, and the source term and Jacobian are generated from it at runtime (`make_chem_funs.py` writes
`chem_funs.py`). You can switch to a different shipped network or edit one yourself.

!!! note "How the network becomes code"
    On each run (unless you pass `-n`), VULCAN parses the network file and regenerates `chem_funs.py`: the production/loss term `chemdf`, the `Gibbs` reversal function, and the analytic Jacobians. It then checks the network for element conservation and duplicate reactions before integrating.

---

## Available networks

The shipped networks live in **`thermo/`**. The default is
`thermo/CHO_photo_network.txt` (set in `config.py`). The file name encodes what a network
covers and how it is used:

- **Element prefix**: the elements the network spans, e.g. `CHO` (C, H, O), `NCHO` (adds N),
  `SNCHO` (adds S), `PHO` (P, H, O), `TiSNCHO` (N-C-H-O plus TiO/VO kinetics). The file's first comment line states this explicitly (e.g. `# VULCAN S-N-C-H-O network…`).
- **`_photo_` vs `_thermo_`**: `photo` networks include a photochemistry section;
  `thermo` networks are thermochemistry-only (no photodissociation).
- **`_full`**: a more complete set than the standard variant (e.g.
  `NCHO_full_photo_network.txt` adds more oxidizing species).
- **`_lowT`**: extended with low-temperature rate fits (pair these with
  `use_lowT_limit_rates` in `config.py`).
- **Target/variant suffixes**: `_Jupiter`, `_earth`, `_oxidising`, `_DMS`, `_C3`, `_2024`,
  `_2025`, etc. denote application- or version-specific variants.

The main networks, grouped by element coverage:

| Elements | Networks (`thermo/…`) |
|---|---|
| C–H–O | `CHO_photo_network.txt` (default), `CHO_photo_network_lowT.txt`, `CHO_thermo_network.txt` |
| N–C–H–O | `NCHO_photo_network.txt`, `NCHO_full_photo_network.txt`, `NCHO_thermo_network.txt`, `NCHO_photo_network_lowT.txt`, `NCHO_photo_network_lowT_Jupiter.txt`, `NCHO_earth_photo_network.txt`, `NCHO_photo_oxidising_network.txt` |
| S–N–C–H–O | `SNCHO_photo_network.txt`, `SNCHO_full_photo_network.txt`, `SNCHO_photo_network_2024.txt`, `SNCHO_photo_network_2025.txt`, `SNCHO_photo_network_C3.txt`, `SNCHO_DMS_photo_network_Tsai2024.txt` |
| P–H–O | `PHO_full_photo_network.txt` |
| N–C–H–O + Ti/V | `TiSNCHO_photo_network.txt` |
| older C/H/N/O sets | `CRAHCNO_network_old.txt`, `CRAHCNO_network_v3.txt` |
| specialized mechanism | `SO3-H2SO4_mechanism.txt` (S oxidation sub-mechanism, not a standalone network) |

Additional variants live in `thermo/Test_networks/`. The other files in `thermo/` are
**support data, not networks**: `all_compose.txt` (species composition table), `NASA9/`
(thermodynamic polynomials), and `photo_cross/` (UV cross sections).

!!! tip "Confirm a network's elements before using it"
    Always read the chosen network's first comment line (and `thermo/README.md`) to confirm
    which elements it actually spans, then set `atom_list` to match (see below). The element
    prefix is a reliable guide, but the header is authoritative.

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
    VULCAN **rewrites the network file in place** on each run, prefixing every reaction with a managed integer index (`read_network` in `make_chem_funs.py`). Any leading number you type is overwritten, so just write the `[ ... ] A B C` line and let the code handle indexing. Forward reactions receive odd numbers (1, 3, 5, …).

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
`op.ReadRate.make_bins_read_cross`). Without them the run stops with a
`Missing the cross section from <SPECIES>` error.

### 4. Run

```sh
python run_vulcan.py
```

---

!!! warning "Editing networks requires understanding the chemistry"
    Reverse rates and element balance depend on the network being physically consistent. Start from an existing network, make small changes, and confirm the element-conservation and duplicate checks pass before trusting the results.