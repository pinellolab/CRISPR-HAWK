"""Configuration handling for CRISPR-HAWK tools.

This module defines the Config class, which loads tool-specific settings from
JSON files and validates that required options such as conda environment name
and output directory are provided.
"""

from .exception_handlers import exception_handler
from .crisprhawk_error import CrisprHawkConfigError

import json
import os


class Config:
    """Represent tool-specific configuration loaded from a JSON file.

    This class stores the conda environment name and output directory for a
    given tool identifier and provides helpers to load and validate them.
    """

    def __init__(self, identifier: str, env_name: str = "", outdir: str = "") -> None:
        self._identifier = identifier
        self._env_name = env_name
        self._outdir = outdir

    def load(self, config_file: str) -> None:
        """Load configuration settings for this identifier from a JSON file.

        This method updates the stored conda environment name and output
        directory if the identifier is present in the configuration file.

        Args:
            config_file: Path to the JSON configuration file to load.
        """
        try:
            with open(config_file, mode="r") as fin:
                config = json.load(fin)
            if self._identifier in config:
                self._env_name = config[self._identifier]["env_name"]
                self._outdir = config[self._identifier]["outdir"]
        except (json.JSONDecodeError, IOError) as e:
            exception_handler(
                CrisprHawkConfigError,
                f"Error loading config file {config_file}",
                os.EX_IOERR,
                True,
                e,
            )

    def validate(self) -> bool:
        """Validate that the configuration is complete and usable.

        This method checks that both the conda environment name and output
        directory have been set for this configuration.

        Returns:
            bool: True if both env_name and outdir are non-empty, otherwise False.
        """
        return bool(self._env_name and self._outdir)

    @property
    def env_name(self) -> str:
        return self._env_name

    @property
    def outdir(self) -> str:
        return self._outdir
