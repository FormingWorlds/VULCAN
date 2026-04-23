#!/usr/bin/env python3
from __future__ import annotations

import os
import sys

_SRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src')
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

from vulcan.vulcan import run_cli, __version__

if __name__ == '__main__':
    run_cli()
