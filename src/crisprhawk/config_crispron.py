""" """

from .config import Config
from .config_utils import CONFIG, set_command, create_mamba_env
from .crisprhawk_error import CrisprHawkCrisprOnConfigError
from .exception_handlers import exception_handler
from .utils import OSSYSTEMS, warning, suppress_stdout, suppress_stderr

from typing import Optional

import platform
import subprocess
import os


# define crispron config identifier
CRISPRON = "crispron"

# define crispron script
CRISPRONSH = "crispron.sh"

# packages required to create crispron environment
CRISPRON_PACKAGES = [
    "bedtools=2.31.1",
    "biopython=1.83",
    "blat=35",
    "bowtie=1.3.1",
    "gff3sort=0.1.*",
    "gffread=0.12.7",
    "matplotlib=3.8.4",
    "matplotlib-base=3.8.4",
    "mysql-connector-python=8.3.0",
    "pandas=2.2.2",
    "pyfaidx=0.8.1.1",
    "samtools=1.20",
    "scikit-learn=1.4.2",
    "tensorboard=2.14",
    "tensorflow=2.14",
    "trf=4.09.1",
    "viennarna=2.6.4",
    "gffutils",
    "pysqlite3",
]


class CrisprOnConfig:
    """Manages configuration for crispron environment and output directory.

    This class loads, saves, and validates configuration settings for the crispron
    tool, including environment name and output directory. It provides methods to
    update, reset, and display configuration, as well as to set the conda/mamba
    command for environment management.

    Attributes:
        _config_file (str): Path to the configuration JSON file.
        _conda (str): The conda or mamba command used for environment management.
        _config (dict): The loaded configuration dictionary.
    """

    def __init__(self) -> None:
        """Initializes a CRISPRonConfig object for managing crispron configuration.

        Loads the configuration file, sets the default conda command, and prepares
        the configuration dictionary.
        """
        self._config_file = CONFIG
        self._conda = set_command()  # decide whether using mamba or conda
        self._config = Config(CRISPRON)
        self._config.load(CONFIG)  # load config file

    def __repr__(self) -> str:
        """Returns a detailed string representation of the CRISPRonConfig object.

        The representation includes the config file path, environment name, and
        output directory.

        Returns:
            str: A detailed string representation of the CRISPRonConfig object.
        """
        return (
            f"<{self.__class__.__name__} object; "
            f"config_file='{self._config_file}', env_name='{self.env_name}', "
            f"outdir='{self.outdir}'>"
        )

    def validate(self) -> None:
        # ensure environment and outdir are configured
        if not self._config.validate():
            exception_handler(
                CrisprHawkCrisprOnConfigError,
                "CRISPRon environment configuration missing arguments",
                os.EX_DATAERR,
                True,
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


def check_crispron_env(env_name: str, conda: str) -> bool:
    """Checks if the CRISPRon environment is available and the executable exists.

    Verifies that the CRISPRon shell script exists on disk and that the specified
    conda or mamba environment is runnable. Issues a warning and returns False if
    either check fails.

    Args:
        env_name (str): The name of the conda or mamba environment to check.
        conda (str): The conda or mamba command to use.

    Returns:
        bool: True if the CRISPRon environment is ready, False otherwise.
    """
    crispron_dir = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), "scores", "crispron"
    )
    crispron_bin = os.path.join(crispron_dir, "bin", "CRISPRon.sh")
    if not os.path.isfile(crispron_bin):
        exception_handler(
            CrisprHawkCrisprOnConfigError,
            f"Unable to locate CRISPRon executable: {crispron_bin}",
            os.EX_DATAERR,
            True,
        )
    try:
        with suppress_stdout(), suppress_stderr():
            subprocess.check_call(
                [conda, "run", "-n", env_name, "python", "-c", "pass"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
    except (FileNotFoundError, subprocess.CalledProcessError):
        warning(
            f"CRISPRon environment '{env_name}' not found, attempting to create it",
            1,
        )
        return False
    return True


def prepare_crispron_env() -> Optional[CrisprOnConfig]:
    if platform.system() != OSSYSTEMS[0]:  # if system is not Linux
        warning(
            f"CRISPRon scoring is only supported on {OSSYSTEMS[0]} "
            "systems. Off-target estimation automatically disabled",
            1,
        )  # always disply this warning
    config = CrisprOnConfig()  # loads config.json
    # look for crispron environment, if not available create it
    if not check_crispron_env(config.env_name, config.conda):
        warning("CRISPRon environment not available, creating environment...", 1)
        if not create_mamba_env(
            config.conda, config.env_name, CRISPRON_PACKAGES, python_version="3.10"
        ):
            warning(
                "CRISPRon environment creation failed, skipping CRISPRon scoring", 1
            )
            return None
    return config
