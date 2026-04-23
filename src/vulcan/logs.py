# Setup logger for VULCAN
from __future__ import annotations

import logging
import os
import sys


class CustomFormatter(logging.Formatter):
    part1 = '[\033['

    info = '32'
    warn = '93'
    debug = '96'
    error = '91'

    part2 = 'm\033[1m %(levelname)-5s \033[21m\033[0m] %(message)s'

    FORMATS = {
        logging.DEBUG: part1 + debug + part2,
        logging.INFO: part1 + info + part2,
        logging.WARNING: part1 + warn + part2,
        logging.ERROR: part1 + error + part2,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


# Custom logger instance
def setup_logger(logpath: str = 'new.log', level: str = 'INFO', logterm: bool = True):
    logger_name = 'fwl'

    # https://stackoverflow.com/a/61457119
    custom_logger = logging.getLogger(logger_name)

    # Clear existing handlers to prevent accumulation when setup_logger is called multiple times
    # (e.g., when CLI command groups are invoked multiple times in tests)
    custom_logger.handlers.clear()

    if os.path.exists(logpath):
        os.remove(logpath)

    level = str(level).strip().upper()
    if level not in ['INFO', 'DEBUG', 'ERROR', 'WARNING']:
        raise ValueError(f'Invalid log level: {level}')
    level_code = logging.getLevelNamesMapping()[level]

    # Add terminal output to logger
    if logterm:
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(CustomFormatter())
        sh.setLevel(level_code)
        custom_logger.addHandler(sh)

    # Add file output to logger
    fh = logging.FileHandler(logpath)
    fh.setFormatter(logging.Formatter('[ %(levelname)-5s ] %(message)s'))

    fh.setLevel(level)
    custom_logger.addHandler(fh)
    custom_logger.setLevel(level_code)

    # Capture unhandled exceptions
    # https://stackoverflow.com/a/16993115
    def handle_exception(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            custom_logger.error('KeyboardInterrupt')
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        custom_logger.critical(
            'Uncaught exception', exc_info=(exc_type, exc_value, exc_traceback)
        )

    sys.excepthook = handle_exception

    return custom_logger
