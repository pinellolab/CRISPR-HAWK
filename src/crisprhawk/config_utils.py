""" """

from .config_crispron import CRISPRON_PACKAGES, CrisprOnConfig, check_crispron_env
from .exception_handlers import exception_handler
from .utils import OSSYSTEMS, warning

from typing import List, Optional

import subprocess
import platform
import shutil
import os

# ------------------------------------------------------------------------------
# 
# Define constant variables (config files)
#
# ------------------------------------------------------------------------------

# config file location
CONFIG = os.path.abspath(os.path.join(os.path.dirname(__file__), "config/config.json"))

# mamba/conda commands
MAMBA = "mamba"
CONDA = "conda"

# ------------------------------------------------------------------------------
#
# Define utilities classes (config files)
#
# ------------------------------------------------------------------------------

class ScoringEnvs:

    def __init__(self) -> None:
        self._crispron = None
        self._sgdesigner = None

    @property
    def crispron_env(self) -> Optional[CrisprOnConfig]:
        return self._crispron
    
    @crispron_env.setter
    def crispron_env(self, value: CrisprOnConfig) -> None:
        if isinstance(value, CrisprOnConfig):
            self._crispron = value

    @property
    def sgdesigner_env(self) -> Optional[CrisprOnConfig]:
        return self._sgdesigner
    
    @sgdesigner_env.setter
    def sgdesigner_env(self, value: CrisprOnConfig) -> None:
        if isinstance(value, CrisprOnConfig):
            self._sgdesigner = value


# ------------------------------------------------------------------------------
# 
# Define utilities functions (config files)
#
# ------------------------------------------------------------------------------

def command_exists(command: str) -> bool:
    """Check if a command exists in the system's PATH.

    Returns True if the specified command is found in the system's executable
    search path, otherwise False.

    Args:
        command (str): The command to check for existence.

    Returns:
        bool: True if the command exists, False otherwise.
    """
    return bool(shutil.which(command))

def set_command() -> str:
    if command_exists(MAMBA):
        return MAMBA
    return CONDA if command_exists(CONDA) else ""

def create_mamba_env(conda: str, env_name: str, packages: Optional[List[str]] = None, python_version: str = "3.8") -> bool:
    cmd = [conda, "create", "-y", "-n", env_name, f"python={python_version}"]
    if packages:
        cmd.extend(packages)
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0
        

def prepare_crispron_env() -> CrisprOnConfig:
    if platform.system() != OSSYSTEMS[0]:  # if system is not Linux
        warning(
            f"CRISPRon scoring is only supported on {OSSYSTEMS[0]} "
            "systems. Off-target estimation automatically disabled",
            1,
        )  # always disply this warning 
    config = CrisprOnConfig()  # loads config.json
    # look for crispron environment, if not available create it
    if not check_crispron_env(config.env_name, config.conda):
        warning("Impossible to create CRISPRon environment, skipping CRISPRon scoring", 1)
    if not create_mamba_env(config.conda, config.env_name, CRISPRON_PACKAGES, python_version="3.10"):
        warning("CRISPRon environment creation failed, skipping CRISPRon scoring", 1)
    return config


def prepare_scoring_envs() -> ScoringEnvs:
    # prepare scoring environments
    scoring_envs = ScoringEnvs()
    scoring_envs.crispron_env = prepare_crispron_env()  # crispron
    return scoring_envs
    
