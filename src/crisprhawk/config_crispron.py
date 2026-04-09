""" """

from .crisprhawk_error import CrisprHawkCRISPRonConfigError
from .exception_handlers import exception_handler
from .utils import warning, suppress_stdout, suppress_stderr, command_exists

from typing import Optional, Dict

import subprocess
import json
import sys
import os

# config file location
CONFIG = os.path.abspath(os.path.join(os.path.dirname(__file__), "config/config.json"))
# mamba/conda commands
MAMBA = "mamba"
CONDA = "conda"
# crispron
crispron = "crispron.py"


class CRISPRonConfig:
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
        self._conda = CONDA  # default to conda
        self._config = self._load_config()

    def __str__(self) -> str:
        """Returns a human-readable string representation of the CRISPRonConfig
        object.

        The string includes the current environment name and output directory.

        Returns:
            str: A string representation of the CRISPRonConfig object.
        """
        return f"CRISPRonConfig(env_name='{self.env_name}', outdir='{self.outdir}')"

    def __repr__(self) -> str:
        """Returns a detailed string representation of the CRISPRonConfig object.

        The representation includes the config file path, environment name, and
        output directory.

        Returns:
            str: A detailed string representation of the CRISPRonConfig object.
        """
        return f"<CRISPRonConfig object; config_file='{self._config_file}', env_name='{self.env_name}', outdir='{self.outdir}'>"

    def _load_config(self) -> dict:
        """Loads the CRISPRon configuration from the config file.

        Reads the configuration JSON file, creates a default config if missing,
        and ensures required fields are present.

        Returns:
            dict: The loaded configuration dictionary.

        Raises:
            CrisprHawkCRISPRonConfigError: If the config file cannot be loaded
                or parsed.
        """
        if not os.path.exists(self._config_file):
            default_config = {
                "crispritz": {
                    "env_name": "crispritz-crisprhawk",
                    "outdir": ".crispritz_targets",
                },
                "crispron": {
                    "env_name": "crispron-crisprhawk",
                    "outdir": ".crispron_outdir",
                },
            }
            self._save_config(default_config)
            return default_config
        try:
            with open(self._config_file, mode="r") as infile:
                config = json.load(infile)
            if "crispron" not in config:  # ensure 'crispron' section exists
                config["crispron"] = {
                    "env_name": "crispron-crisprhawk",
                    "outdir": ".crispron_outdir",
                }
                self._save_config(config)
            return config
        except (json.JSONDecodeError, IOError) as e:  
            exception_handler(
                CrisprHawkCRISPRonConfigError,
                f"Error loading config file {self._config_file}",
                os.EX_IOERR,
                True,
                e,
            )

    def _save_config(self, config: Optional[Dict[str, Dict[str, str]]] = None) -> None:
        """Saves the crispron configuration to the config file.

        Writes the provided configuration dictionary (or the current config) to
        the JSON config file, creating the directory if necessary.

        Args:
            config (Optional[Dict[str, Dict[str, str]]]): The configuration dictionary
                to save. If None, saves the current config.

        Raises:
            CrisprHawkCRISPRonConfigError: If the config file cannot be saved.
        """
        config_to_save = config if config is not None else self._config
        try:  # ensure config directory exists
            if not os.path.isdir(os.path.dirname(self._config_file)):
                os.makedirs(os.path.dirname(self._config_file), exist_ok=True)
            with open(self._config_file, mode="w") as outfile:  # write config json
                json.dump(config_to_save, outfile, indent=2, sort_keys=True)
        except (IOError, Exception) as e:  # always traced
            exception_handler(
                CrisprHawkCRISPRonConfigError,
                f"Error saving config file {self._config}",
                os.EX_IOERR,
                True,
                e,
            )

    @property
    def env_name(self) -> str:
        return self._config["crispron"].get("env_name")

    @env_name.setter
    def env_name(self, value: str) -> None:
        """Sets the environment name for the crispron configuration.

        Validates that the provided value is a non-empty string and updates the
        configuration.

        Args:
            value (str): The new environment name to set.

        Raises:
            CrisprHawkCRISPRonConfigError: If the value is not a non-empty string.
        """
        if not isinstance(value, str) or not value.strip():
            exception_handler(
                CrisprHawkCRISPRonConfigError,
                f"Environment name must be a non-empty string, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._config["crispron"]["env_name"] = value.strip()
        self._save_config()

    @property
    def outdir(self) -> str:
        return self._config["crispron"].get("outdir")

    @outdir.setter
    def outdir(self, value: str) -> None:
        """Sets the output directory for the crispron configuration.

        Validates that the provided value is a string and updates the configuration.

        Args:
            value (str): The new output directory to set.

        Raises:
            CrisprHawkCRISPRonConfigError: If the value is not a string.
        """
        if not isinstance(value, str):
            exception_handler(
                CrisprHawkCRISPRonConfigError,
                f"Output directory must be a string, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._config["crispron"]["outdir"] = value
        self._save_config()

    @property
    def conda(self) -> str:
        return self._conda

    @conda.setter
    def conda(self, command: str) -> None:
        """Sets the conda or mamba command for environment management.

        Updates the internal command used to manage CRISPRon environments.

        Args:
            command (str): The conda or mamba command to set.
        """
        self._conda = command

    def set_command(self) -> bool:
        """Sets the conda or mamba command for CRISPRon environment management.

        This method checks for the availability of mamba or conda and sets the
        internal command accordingly. If neither is found, it issues a warning and
        returns False.

        Returns:
            bool: True if a valid command was set, False otherwise.
        """
        if command_exists(MAMBA):
            self.conda = MAMBA
            return True
        elif command_exists(CONDA):
            self.conda = CONDA
            return True
        warning(f"Neither {CONDA} or {MAMBA} found, skipping off-targets estimation", 1)
        return False

    def update_config(
        self, env_name: Optional[str] = None, outdir: Optional[str] = None
    ) -> None:
        """Updates the crispron configuration with new environment name and/or
        output directory.

        This method sets the environment name and output directory if provided,
        updating the configuration accordingly.

        Args:
            env_name (Optional[str]): The new environment name to set. If None,
                the environment name is not changed.
            outdir (Optional[str]): The new output directory to set. If None, the
                output directory is not changed.
        """
        if env_name is not None:
            self.env_name = env_name
        if outdir is not None:
            self.outdir = outdir

    def reset_to_defaults(self) -> None:
        """Resets the crispron configuration to default values.

        This method restores the environment name and output directory to their
        default settings and saves the updated configuration.
        """
        self._config["crispron"] = {
            "env_name": "crispron-crisprhawk",
            "outdir": ".crispron_outdir",
        }
        self._save_config()

    def show_config(self) -> str:
        """Returns the current crispron configuration as a formatted JSON string.

        This method provides a human-readable representation of the loaded
        configuration.

        Returns:
            str: The configuration as a formatted JSON string.
        """
        return json.dumps(self._config, indent=2, sort_keys=True)

    def validate_config(self) -> None:
        """Validates the current crispron configuration for required fields and
        types.

        This method checks that the configuration contains the necessary fields
        and that they are properly formatted.

        Raises:
            CrisprHawkCRISPRonConfigError: If any required field is missing or
                has an invalid type or value.
        """
        if "crispron" not in self._config:  # check if crispron section exists
            exception_handler(
                CrisprHawkCRISPRonConfigError,
                "Field 'crispron' not found in config file",
                os.EX_DATAERR,
                False,
            )
        crispritz_config = self._config["crispron"]  # check required fields
        required_fields = ["env_name", "outdir"]
        for field in required_fields:
            if field not in crispritz_config:
                exception_handler(
                    CrisprHawkCRISPRonConfigError,
                    f"Missing field '{field}' from config file",
                    os.EX_DATAERR,
                    False,
                )
            if not isinstance(crispritz_config[field], str):
                exception_handler(
                    CrisprHawkCRISPRonConfigError,
                    f"Wrong type format in field '{field}'",
                    os.EX_DATAERR,
                    False,
                )
        if not crispritz_config["env_name"].strip():  # check env_name is not empty
            exception_handler(
                CrisprHawkCRISPRonConfigError,
                "Missing 'env_name' in config file",
                os.EX_DATAERR,
                False,
            )
        if not crispritz_config["outdir"].strip():  # check env_name is not empty
            exception_handler(
                CrisprHawkCRISPRonConfigError,
                "Missing 'outdir' in config file",
                os.EX_DATAERR,
                False,
            )


def check_crispron_env(env_name: str, conda: str) -> bool:
    crispron_dir = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), "scores", "crispron"
    )
    crispron_bin = os.path.join(crispron_dir, "bin", "CRISPRon.sh")

    if not os.path.isfile(crispron_bin):
        warning(f"Cannot find CRISPRon executable: {crispron_bin}", 1)
        return False

    try:
        with suppress_stdout(), suppress_stderr():
            subprocess.check_call(
                [conda, "run", "-n", env_name, "python", "-c", "pass"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                cwd=crispron_dir,
            )
    except (FileNotFoundError, subprocess.CalledProcessError):
        warning(
            f"CRISPRon environment '{env_name}' not found or not runnable, skipping CRISPRon scoring",
            1,
        )
        return False

    return True


def config_crispron(config: CRISPRonConfig, env_name: str, targets_dir: str) -> None:
    """Configures the CRISPRon environment and output directory.

    This function updates the CRISPRon configuration with a new environment name
    and/or output directory, reloads the configuration, and validates it.

    Args:
        config (CRISPRonConfig): The CRISPRon configuration object to update.
        env_name (str): The new environment name to set.
        targets_dir (str): The new output directory to set.
    """
    if env_name:  # set new enviornment name
        sys.stdout.write(f"Setting 'env_name' = {env_name}\n")
        config.env_name = env_name
    if targets_dir:  # set new targets output folder
        sys.stdout.write(f"Setting 'outdir' = {targets_dir}\n")
        config.outdir = targets_dir
    config = CRISPRonConfig()  # reload new config file
    config.validate_config()  # validate config
