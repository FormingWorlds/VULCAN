# Editing or using a different chemical network

VULCAN is developed in a flexible way, so that the chemical network is _not_ hard coded. You can switch networks or edit the default one.

## Switching networks

Edit `config.py`:

```python
network = 'thermo/CHO_photo_network.txt'  # Use C, H, O only
```

Then run:
```bash
python vulcan.py
```

Do **NOT** use `-n` flag when changing networks!

## Editing the default network

### Network file format

Each reaction is written as:
```
[ Reactant1 + Reactant2 -> Product1 + Product2 ]  A  B  C
```

Where `A`, `B`, `C` are Arrhenius coefficients: $k = A T^B \exp(-C/T)$

Example from `NCHO_photo_network.txt`:
```
[ OH + H2 -> H2O + H ]  3.57E-16    1.520   1740.0
[ O + H2 -> OH + H   ]  8.52E-20    2.670   3160.0
```

### Steps to edit

1. Open the network file:
```bash
nano thermo/NCHO_photo_network.txt
```

2. **Add a reaction:** Insert at appropriate location:
```
[ New_species1 + New_species2 -> Product1 + Product2 ]  A   B   C
```
where you replace `A`, `B`, `C` with the appropriate Arrhenius coefficients. 

3. **Remove a reaction:** Delete the line

4. **Change rate coefficients:** Edit the A, B, C values

5. **Save and run:**
```bash
python vulcan.py
```

VULCAN will regenerate `chem_funs.py` with your changes.

## Adding new species

If you add new species, you must add NASA-9 thermodynamic data:

1. Find the coefficients in `thermo/NASA9/nasa9_2002_E.txt` or `new_nasa9.txt`

2. Create a file in `thermo/NASA9/` with the species name, e.g., `MY_SPECIES.txt`. The format of the NASA 9 polynomials is as follows:
```
a1 a2 a3 a4 a5
a6 a7 0. a8 a9
```

**Note the "0." separator between a7 and a8!**

3. Check `thermo/all_compose.txt` for your new species. If you cannot find your species, add it in the same format, e.g.:

```
species			H	O	C	He	N	S	P	Na	K	Si	Fe  Ar 	Ti	V	Mg	Ca	e 	mass    
H2O				2	1	0	0	0	0	0	0	0	0	0	0	0   0	0	0   0   18.016
```

4. Run VULCAN:
```bash
python vulcan.py
```


## Important notes

- **Forward reactions only** - List only forward reactions; VULCAN automatically computes reverse reactions from thermodynamics
- **Three-body reactions and photodissociation reactions** - Listed separately after a comment line in default networks
- **Don't skip chemical regeneration** - Always run WITHOUT `-n` after editing networks!

---

**Remember:** Editing reaction networks requires understanding the chemistry. Start with existing networks and make small changes!
