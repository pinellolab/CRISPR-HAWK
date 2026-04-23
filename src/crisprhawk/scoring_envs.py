""" """

from .config_crispron import CrisprOnConfig, prepare_crispron_env

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

    # @property
    # def sgdesigner_env(self) -> Optional[CrisprOnConfig]:
    #     return self._sgdesigner

    # @sgdesigner_env.setter
    # def sgdesigner_env(self, value: CrisprOnConfig) -> None:
    #     if isinstance(value, CrisprOnConfig):
    #         self._sgdesigner = value


def prepare_scoring_envs() -> ScoringEnvs:
    # prepare scoring environments
    scoring_envs = ScoringEnvs()
    scoring_envs.crispron_env = prepare_crispron_env()  # crispron
    return scoring_envs
