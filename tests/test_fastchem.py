from __future__ import annotations

import os
import sys
import types
from pathlib import Path

import numpy as np
import pytest

try:
    from vulcan import build_atm
except ModuleNotFoundError as exc:
    if exc.name != 'vulcan.chem_funs':
        raise
    # chem_funs.py is generated at runtime, so provide a minimal test stub.
    sys.modules['vulcan.chem_funs'] = types.SimpleNamespace(ni=0, spec_list=())
    from vulcan import build_atm

from vulcan.build_atm import InitialAbun
from vulcan.config import Config


@pytest.mark.unit
def test_ini_fc_raises_when_fastchem_binary_missing(monkeypatch):
    cfg = Config()
    init = InitialAbun(cfg)
    atm = types.SimpleNamespace(pco=np.array([1.0e6]), Tco=np.array([1000.0]))

    monkeypatch.setattr(build_atm.os.path, 'isfile', lambda _: False)

    with pytest.raises(RuntimeError, match='FastChem cannot be found'):
        init.ini_fc(data_var=types.SimpleNamespace(), data_atm=atm)


@pytest.mark.smoke
def test_ini_fc_invokes_fastchem_with_explicit_args(monkeypatch, tmp_path):
    fastchem_dir = tmp_path / 'fastchem'
    input_dir = fastchem_dir / 'input'
    output_dir = fastchem_dir / 'output'
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)

    (fastchem_dir / 'fastchem').write_text('#!/bin/sh\n')
    (input_dir / 'config.input').write_text('config\n')
    (input_dir / 'parameters_wo_ion.dat').write_text('wo ion\n')
    (input_dir / 'parameters_ion.dat').write_text('ion\n')
    (input_dir / 'solar_element_abundances.dat').write_text('H 12.00\nC 8.43\nO 8.69\n')

    cfg = Config()
    cfg.use_solar = True
    cfg.use_ion = False
    init = InitialAbun(cfg)
    atm = types.SimpleNamespace(pco=np.array([1.0e6]), Tco=np.array([1000.0]))

    monkeypatch.setattr(build_atm, 'FASTCHEM_DIR', f'{fastchem_dir}{os.sep}')

    def fake_copyfile(src, dst):
        Path(dst).write_text(Path(src).read_text())
        return dst

    called: dict[str, object] = {}

    class _FakeStdout:
        @staticmethod
        def readline():
            return b''

    class _FakeProcess:
        stdout = _FakeStdout()

        @staticmethod
        def poll():
            return 0

    def fake_popen(cmd, stdout, stderr, cwd):
        called['cmd'] = cmd
        called['cwd'] = cwd
        called['stdout'] = stdout
        called['stderr'] = stderr
        return _FakeProcess()

    monkeypatch.setattr(build_atm, 'copyfile', fake_copyfile)
    monkeypatch.setattr(build_atm.subprocess, 'Popen', fake_popen)

    init.ini_fc(data_var=types.SimpleNamespace(), data_atm=atm)

    assert called['cmd'] == [
        os.path.join(f'{fastchem_dir}{os.sep}', 'fastchem'),
        os.path.join(f'{fastchem_dir}{os.sep}', 'input', 'config.input'),
    ]
    assert called['cwd'] == f'{fastchem_dir}{os.sep}'
    assert (input_dir / 'element_abundances_vulcan.dat').exists()
    assert (input_dir / 'vulcan_TP' / 'vulcan_TP.dat').exists()
