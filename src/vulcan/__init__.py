from __future__ import annotations

__version__ = '26.10.16'

from .config import Config
from .vulcan import run_vulcan

__all__ = ['__version__', 'run_vulcan', 'Config', 'plot']
