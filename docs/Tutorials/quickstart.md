# VULCAN Quick Start 

Run your first photochemical simulation!

## 1. Run the default demo

From the main VULCAN directory:

```bash
python vulcan.py
```

This will:

- Generate chemical functions from the reaction network
- Initialize an Earth-like atmosphere
- Run the simulation with real-time plotting
- Complete in ~10-15 minutes

Check the results in `output/`:
```bash
ls output/
```

---

## 2. Modify the configuration

Edit `config.py` to change simulation parameters:

**Weaker vertical mixing:**
```python
const_Kzz = 1.E7
```

**Carbon-rich atmosphere (C/O = 1):**
```python
C_H = 6.0618E-4
```

**Disable plotting (faster):**
```python
use_live_plot = False
```

Then run again:
```bash
python vulcan.py
```

---

## 3. Try planet examples

VULCAN includes examples for different planets, listed in `cgf_examples`. 

Available examples:

- `cfg_Earth.txt` - Earth
- `cfg_HD189.txt` - HD 189733b (hot Jupiter)
- `cfg_HD209.txt` - HD 209458b (hot Jupiter)
- `cfg_Jupiter.txt` - Jupiter

To try one of these examples, update `config.py` with the parameters as listed in the example `.txt` files. 

---

## 4. Plot results

Run ```plot_vulcan.py``` within ```tools```:

```
python plot_vulcan.py [vulcan output] [species] [plot name] [-h (for plotting height)]
```

---

## 5. Speed tips

Skip regenerating chemistry if you only changed `config.py`:
```bash
python vulcan.py -n
```

Only use `-n` if you haven't modified the reaction network!

---

## Troubleshooting

**FastChem not found:**
```bash
cd fastchem_vulcan && make && cd ..
```

**Missing Python modules:**
```bash
pip install -e .
```

**Simulation too slow:**
```python
# In config.py:
use_live_plot = False
n_layer = 100  # Reduce resolution
```