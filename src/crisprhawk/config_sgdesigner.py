""" """

from .config import Config
from .config_utils import CONFIG, set_command, create_mamba_env
from .crisprhawk_error import CrisprHawkSgDesignerError
from .exception_handlers import exception_handler
from .utils import OSSYSTEMS, warning, suppress_stderr, suppress_stdout

from typing import Optional

import platform
import subprocess
import os


# define sgdesigner config identifier
SGDESIGNER = "sgdesigner"

# packages required to create sgdesigner environment
SGDESIGNER_PACKAGES = [
    "numpy=1.15.2", 
    "scipy=1.1.0", 
    "scikit-learn=0.20.0",
    "xgboost=0.80",
    "joblib=0.13.2",
]


class sgDesignerConfig:
    """Manages configuration for sgDesigner environment and output directory.

    This class loads, saves, and validates configuration settings for the sgdesigner
    tool, including environment name and output directory. It provides methods to
    update, reset, and display configuration, as well as to set the conda/mamba
    command for environment management.

    Attributes:
        _config_file (str): Path to the configuration JSON file.
        _conda (str): The conda or mamba command used for environment management.
        _config (dict): The loaded configuration dictionary.
    """

    def __init__(self) -> None:
        """Initializes a sgDesignerConfig object for managing sgDesigner configuration.

        Loads the configuration file, sets the default conda command, and prepares
        the configuration dictionary.
        """
        self._config_file = CONFIG
        self._conda = set_command()  # decide whether using mamba or conda
        self._config = Config(SGDESIGNER)
        self._config.load(CONFIG)  # load config file

    def __repr__(self) -> str:
        """Returns a detailed string representation of the sgdesignerConfig object.

        The representation includes the config file path, environment name, and
        output directory.

        Returns:
            str: A detailed string representation of the sgdesignerConfig object.
        """
        return (
            f"<{self.__class__.__name__} object; "
            f"config_file='{self._config_file}', env_name='{self.env_name}', "
            f"outdir='{self.outdir}'>"
        )

    def validate(self) -> None:
        # ensure environment and outdir are configured
        if not self._config.validate():
            exception_handler(CrisprHawkSgDesignerError, "SgDesigner environment configuration missing arguments", os.EX_DATAERR, True)

    @property
    def env_name(self) -> str:
        return self._config.env_name
    
    @property
    def outdir(self) -> str:
        return self._config.outdir
    
    @property
    def conda(self) -> str:
        return self._conda


def check_sgdesigner_env(env_name: str, conda: str) -> bool:
    sgdesigner_dir = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), "scores", "sgdesigner"
    )
    sgdesigner_bin = os.path.join(sgdesigner_dir, "sgDesigner.pl")
    if not os.path.isfile(sgdesigner_bin):
        exception_handler(
            CrisprHawkSgDesignerError,
            f"Unable to locate sgDesigner executable: {sgdesigner_bin}",
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
            f"sgDesigner environment '{env_name}' not found, attempting to create it",
            1,
        )
        return False
    return True

def prepare_sgdesigner_env() -> Optional[sgDesignerConfig]:
    if platform.system() != OSSYSTEMS[0]:  # if system is not Linux
        warning(
            f"sgDesigner scoring is only supported on {OSSYSTEMS[0]} "
            "systems. Off-target estimation automatically disabled",
            1,
        )  # always disply this warning
    config = sgDesignerConfig()  # loads config.json
    # look for crispron environment, if not available create it
    if not check_sgdesigner_env(config.env_name, config.conda):
        warning(
            "sgDesigner environment not available, creating environment...", 1
        )
        if not create_mamba_env(
            config.conda, config.env_name, SGDESIGNER_PACKAGES, python_version="3.7"
        ):
            warning("sgDesigner environment creation failed, skipping CRISPRon scoring", 1)
            return None
    return config

