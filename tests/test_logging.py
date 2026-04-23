from __future__ import annotations

import pytest

from vulcan.logs import setup_logger


def test_setup_logger_rejects_invalid_level(tmp_path):
    with pytest.raises(ValueError, match='Invalid log level'):
        setup_logger(logpath=str(tmp_path / 'vulcan.log'), level='verbose', logterm=False)


@pytest.mark.unit
def test_setup_logger_writes_to_file(tmp_path):
    logpath = tmp_path / 'vulcan.log'
    logger = setup_logger(logpath=str(logpath), level='INFO', logterm=False)

    logger.info('hello from test')

    assert logpath.exists()
    assert 'hello from test' in logpath.read_text()
