from __future__ import annotations

import os

# Path to VULCAN repository folder
_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
VULCAN_DIR = os.path.abspath(os.path.join(_PACKAGE_DIR, '..', '..')) + '/'

# Path to FastChem folder (custom version packaged with VULCAN)
FASTCHEM_DIR = os.path.join(VULCAN_DIR, 'fastchem_vulcan') + '/'

# Thermodynamic and chemistry data
THERMO_DIR = os.path.join(
    VULCAN_DIR,
    'thermo',
)

# Path to AGNI folder
AGNI_DIR = os.path.join(VULCAN_DIR, 'AGNI') + '/'

# Composition and gas properties
COM_FILE = os.path.join(VULCAN_DIR, 'thermo', 'all_compose.txt')

# (all the nasa9 files must be placed in the folder: thermo/NASA9/)
GIBBS_FILE = os.path.join(VULCAN_DIR, 'thermo', 'gibbs_text.txt')

# Photochemistry cross-sections
CROSS_DIR = os.path.join(VULCAN_DIR, 'thermo', 'photo_cross') + '/'

# Symbolic chemical functions
CHEM_FUNS_FILE = os.path.join(VULCAN_DIR, 'src', 'vulcan', 'chem_funs.py')
