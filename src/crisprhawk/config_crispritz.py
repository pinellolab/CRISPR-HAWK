""" """

from .config import Config
from .config_utils import CONFIG, set_command
from .crisprhawk_error import CrisprHawkCrispritzConfigError
from .exception_handlers import exception_handler
from .utils import warning, suppress_stdout, suppress_stderr

from typing import Optional, Dict

import subprocess
import json
import sys
import os


# define crispritz config identifier
CRISPRITZ = "crispritz"

# define crispritz bin
CRISPRITZ = "crispritz.py"


class CrispritzConfig:
    """Manages configuration for CRISPRitz environment and output directory.

    This class loads, saves, and validates configuration settings for the CRISPRitz
    tool, including environment name and output directory. It provides methods to
    update, reset, and display configuration, as well as to set the conda/mamba
    command for environment management.

    Attributes:
        _config_file (str): Path to the configuration JSON file.
        _conda (str): The conda or mamba command used for environment management.
        _config (dict): The loaded configuration dictionary.
    """

    def __init__(self) -> None:
        """Initializes a CrispritzConfig object for managing CRISPRitz configuration.

        Loads the configuration file, sets the default conda command, and prepares
        the configuration dictionary.
        """
        self._config_file = CONFIG
        self._conda = set_command()  # decide whether to use mamba or conda
        self._config = Config(CRISPRITZ)
        self._config.load(CONFIG)  # load config file

    def __repr__(self) -> str:
        """Returns a detailed string representation of the CrispritzConfig object.

        The representation includes the config file path, environment name, and
        output directory.

        Returns:
            str: A detailed string representation of the CrispritzConfig object.
        """
        return (
            f"<{self.__class__.__name__} object; "
            f"config_file='{self._config_file}', env_name='{self.env_name}', "
            f"outdir='{self.outdir}'>"
        )

    @property
    def env_name(self) -> str:
        return self._config.env_name

    @property
    def outdir(self) -> str:
        return self._config.outdir

    @property
    def conda(self) -> str:
        return self._conda

    def validate(self) -> None:
        # ensure environment and outdir are configured
        if not self._config.validate():
            exception_handler(
                CrisprHawkCrispritzConfigError,
                "CRISPRon environment configuration missing arguments",
                os.EX_DATAERR,
                True,
            )


def check_crispritz_env(env_name: str, conda: str) -> bool:
    """Checks if CRISPRitz is installed in the specified conda or mamba environment.

    Attempts to run CRISPRitz in the given environment and returns True if successful,
    otherwise issues a warning and returns False.

    Args:
        env_name (str): The name of the conda or mamba environment to check.
        conda (str): The conda or mamba command to use.

    Returns:
        bool: True if CRISPRitz is installed and runnable, False otherwise.
    """
    # check whether crispritz is installed within the environment
    # if not, skip off-targets estimation
    try:
        with suppress_stdout(), suppress_stderr():
            subprocess.check_call(
                f"{conda} run -n {env_name} {CRISPRITZ}",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
    except (FileNotFoundError, subprocess.CalledProcessError):
        warning(
            f"CRISPRitz not installed in environment {env_name}, skipping off-targets estimation",
            1,
        )
        return False
    return True
