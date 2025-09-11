""" """

from .crisprhawk_error import CrisprHawkCrispritzConfigError
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
# CRISPRitz
CRISPRITZ = "crispritz.py"


class CrispritzConfig:

    def __init__(self):
        self._config_file = CONFIG
        self._conda = CONDA  # default to conda
        self._config = self._load_config()

    def __str__(self) -> str:
        return f"CrispritzConfig(env_name='{self.env_name}', outdir='{self.outdir}')"

    def __repr__(self) -> str:
        return f"<CrispritzConfig object; config_file='{self._config_file}', env_name='{self.env_name}', outdir='{self.outdir}'>"

    def _load_config(self) -> dict:
        if not os.path.exists(self._config_file) and os.path.isfile(self._config_file):
            # create default config if file doesn't exist
            default_config = {
                "crispritz": {"env_name": "crispritz", "outdir": ".crispritz_targets"}
            }
            self._save_config(default_config)
            return default_config
        try:
            with open(self._config_file, mode="r") as infile:
                config = json.load(infile)
            if "crispritz" not in config:  # ensure 'crispritz' section exists
                config["crispritz"] = {
                    "env_name": "crispritz",
                    "outdir": ".crispritz_targets",
                }
                self._save_config(config)
            return config
        except (json.JSONDecodeError, IOError) as e:  # always traced
            exception_handler(
                CrisprHawkCrispritzConfigError,
                f"Error loading config file {self._config_file}",
                os.EX_IOERR,
                True,
                e,
            )

    def _save_config(self, config: Optional[Dict[str, Dict[str, str]]] = None) -> None:
        config_to_save = config if config is not None else self._config
        try:  # ensure config directory exists
            if not os.path.isdir(os.path.dirname(self._config_file)):
                os.makedirs(os.path.dirname(self._config_file), exist_ok=True)
            with open(self._config_file, mode="w") as outfile:  # write config json
                json.dump(config_to_save, outfile, indent=2, sort_keys=True)
        except (FileNotFoundError, IOError, Exception) as e:  # always traced
            exception_handler(
                CrisprHawkCrispritzConfigError,
                f"Error saving config file {self._config}",
                os.EX_IOERR,
                True,
                e,
            )

    @property
    def env_name(self) -> str:
        return self._config["crispritz"].get("env_name")

    @env_name.setter
    def env_name(self, value: str) -> None:
        if not isinstance(value, str) or not value.strip():
            exception_handler(
                CrisprHawkCrispritzConfigError,
                f"Environment name must be a non-empty string, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._config["crispritz"]["env_name"] = value.strip()
        self._save_config()

    @property
    def outdir(self) -> str:
        return self._config["crispritz"].get("outdir")

    @outdir.setter
    def outdir(self, value: str) -> None:
        if not isinstance(value, str):
            exception_handler(
                CrisprHawkCrispritzConfigError,
                f"Output directory must be a string, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._config["crispritz"]["outdir"] = value
        self._save_config()

    @property
    def conda(self) -> str:
        return self._conda

    @conda.setter
    def conda(self, command: str) -> None:
        self._conda = command

    def set_env_name(self, env_name: str) -> None:
        self.env_name = env_name

    def set_outdir(self, outdir: str) -> None:
        self.outdir = outdir

    def set_command(self) -> bool:
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
        if env_name is not None:
            self.env_name = env_name
        if outdir is not None:
            self.outdir = outdir

    def reset_to_defaults(self) -> None:
        self._config["crispritz"] = {
            "env_name": "crispritz",
            "outdir": ".crispritz_targets",
        }
        self._save_config()

    def show_config(self) -> str:
        return json.dumps(self._config, indent=2, sort_keys=True)

    def validate_config(self) -> None:
        if "crispritz" not in self._config:  # check if crispritz section exists
            exception_handler(
                CrisprHawkCrispritzConfigError,
                "Field 'crispritz' not found in config file",
                os.EX_DATAERR,
                False,
            )
        crispritz_config = self._config["crispritz"]  # check required fields
        required_fields = ["env_name", "outdir"]
        for field in required_fields:
            if field not in crispritz_config:
                exception_handler(
                    CrisprHawkCrispritzConfigError,
                    f"Missing field '{field}' from config file",
                    os.EX_DATAERR,
                    False,
                )
            if not isinstance(crispritz_config[field], str):
                exception_handler(
                    CrisprHawkCrispritzConfigError,
                    f"Wrong type format in field '{field}'",
                    os.EX_DATAERR,
                    False,
                )
        if not crispritz_config["env_name"].strip():  # check env_name is not empty
            exception_handler(
                CrisprHawkCrispritzConfigError,
                "Missing 'env_name' in config file",
                os.EX_DATAERR,
                False,
            )
        if not crispritz_config["outdir"].strip():  # check env_name is not empty
            exception_handler(
                CrisprHawkCrispritzConfigError,
                "Missing 'outdir' in config file",
                os.EX_DATAERR,
                False,
            )


def check_crispritz_env(env_name: str, conda: str) -> bool:
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


def config_crispritz(config: CrispritzConfig, env_name: str, targets_dir: str) -> None:
    if env_name:  # set new enviornment name
        sys.stdout.write(f"Setting 'env_name' = {env_name}\n")
        config.env_name = env_name
    if targets_dir:  # set new targets output folder
        sys.stdout.write(f"Setting 'outdir' = {targets_dir}\n")
        config.outdir = targets_dir
    config = CrispritzConfig()  # reload new config file
    config.validate_config()  # validate config
