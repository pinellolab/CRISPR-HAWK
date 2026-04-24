""" """

from .config_crispron import CrisprOnConfig, prepare_crispron_env
from .config_sgdesigner import sgDesignerConfig, prepare_sgdesigner_env

from typing import Optional


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
    def sgdesigner_env(self) -> Optional[sgDesignerConfig]:
        return self._sgdesigner

    @sgdesigner_env.setter
    def sgdesigner_env(self, value: sgDesignerConfig) -> None:
        if isinstance(value, sgDesignerConfig):
            self._sgdesigner = value


def prepare_scoring_envs() -> ScoringEnvs:
    # prepare scoring environments
    scoring_envs = ScoringEnvs()
    if crispron_env := prepare_crispron_env():  # crispron
        scoring_envs.crispron_env = crispron_env
    if sgdesigner_env := prepare_sgdesigner_env():  # sgdesigner
        scoring_envs.sgdesigner_env = sgdesigner_env
    return scoring_envs
