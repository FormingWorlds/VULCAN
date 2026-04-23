from __future__ import annotations

import subprocess as sp
from pathlib import Path

import pytest

from vulcan.paths import FASTCHEM_DIR


@pytest.mark.unit
def test_fastchem_exists():

    # Fastchem dir
    assert Path(FASTCHEM_DIR).exists(), 'The FASTCHEM_DIR does not exist.'
    assert Path(FASTCHEM_DIR).is_dir(), 'The FASTCHEM_DIR is not a directory.'

    # Fastchem executable
    FASTCHEM_EXE = Path(FASTCHEM_DIR) / 'fastchem'
    assert FASTCHEM_EXE.exists(), 'The FC executable does not exist.'
    assert FASTCHEM_EXE.is_file(), 'The FC executable is not a file.'

    # Try running it to check if it is executable
    FASTCHEM_CFG = Path(FASTCHEM_DIR) / 'input' / 'config.input'
    try:
        result = sp.run(
            [str(FASTCHEM_EXE), str(FASTCHEM_CFG)],
            capture_output=True,
            text=True,
            cwd=FASTCHEM_DIR,
            check=True,
        )
    except sp.CalledProcessError as e:
        pytest.fail(f'Running the FC executable failed: {e}\nstderr: {e.stderr}')

    # Check log result
    assert 'Equilibrium calculation (Fastchem) finished.' in result.stdout
