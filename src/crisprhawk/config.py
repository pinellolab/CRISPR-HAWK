""" """

from .exception_handlers import exception_handler
from .crisprhawk_error import CrisprHawkConfigError

import json
import os

class Config:

    def __init__(self, identifier: str, env_name: str = "", outdir: str = "") -> None:
        self._identifier = identifier
        self._env_name = env_name
        self._outdir = outdir

    def load(self, config_file: str) -> None:
        try:
            with open(config_file, mode="r") as fin:
                config = json.load(fin)
            if self._identifier in config:
                self._env_name = config[self._identifier]["env_name"]
                self._outdir = config[self._identifier]["outdir"]
        except (json.JSONDecodeError, IOError) as e:
            exception_handler(CrisprHawkConfigError, f"Error loading config file {config_file}", os.EX_IOERR, True, e)

    def validate(self) -> bool:
        return bool(self._env_name and self._outdir)
        

    @property
    def env_name(self) -> str:
        return self._env_name
    
    @property
    def outdir(self) -> str:
        return self._outdir
