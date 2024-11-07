"""
"""

from utils import COMMAND, IUPAC, VERBOSITYLVL
from crisprhawk_version import __version__

from argparse import (
    SUPPRESS,
    ArgumentParser,
    HelpFormatter,
    Action,
    _MutuallyExclusiveGroup,
    Namespace,
)
from typing import Iterable, Optional, TypeVar, Tuple, Dict, NoReturn
from colorama import Fore
from glob import glob


import sys
import os

# define abstract generic types for typing
_D = TypeVar("_D")
_V = TypeVar("_V")


class CrisprHawkArgumentParser(ArgumentParser):

    class CrisprHawkHelpFormatter(HelpFormatter):

        def add_usage(
            self,
            usage: str,
            actions: Iterable[Action],
            groups: Iterable[_MutuallyExclusiveGroup],
            prefix: Optional[str] = None,
        ) -> None:
            # add usage description for help only if the set action is not to
            # suppress the display of the help formatter
            if usage != SUPPRESS:
                args = (usage, actions, groups, "")
                self._add_item(self._format_usage, args)  # initialize the formatter

    def __init__(self, *args: Tuple[_D], **kwargs: Dict[_D, _V]) -> None:
        # set custom help formatter defined as
        kwargs["formatter_class"] = self.CrisprHawkHelpFormatter
        # replace the default version display in usage help with a custom
        # version display formatter
        kwargs["usage"] = kwargs["usage"].replace("{version}", __version__)
        # initialize argument parser object with input parameters for
        # usage display
        super().__init__(*args, **kwargs)

    def error(self, error: str) -> NoReturn:
        # display error messages raised by argparse in red
        errormsg = (
            f"{Fore.RED}\nERROR: {error}.{Fore.RESET}"
            + f"\n\nRun {COMMAND} -h for usage\n\n"
        )
        sys.stderr.write(errormsg)  # write error to stderr
        sys.exit(os.EX_USAGE)  # exit execution -> usage error

    def error_noargs(self) -> None:
        self.print_help()  # if no input argument, print help
        sys.exit(os.EX_NOINPUT)  # exit with no input code


class CrisprHawkInputArgs:
    def __init__(self, args: Namespace, parser: CrisprHawkArgumentParser) -> None:
        self._args = args
        self._parser = parser
        self._check_consistency()  # check input args consistency

    def _check_consistency(self):
        # fasta file
        if not os.path.exists(self._args.fasta) or not os.path.isfile(self._args.fasta):
            self._parser.error(f"Cannot find input FASTA {self._args.fasta}")
        if os.stat(self._args.fasta).st_size <= 0:
            self._parser.error(f"{self._args.fasta} is empty")
        # fasta index
        if self._args.fasta_idx and (
            not os.path.exists(self._args.fasta) or not os.path.isfile(self._args.fasta)
        ):
            self._parser.error(f"Cannot find input FASTA index {self._args.fasta_idx}")
        if self._args.fasta_idx and os.stat(self._args.fasta_idx).st_size <= 0:
            self._parser.error(f"{self._args.fasta_idx} is empty")
        # bed file
        if not os.path.exists(self._args.bedfile) or not os.path.isfile(
            self._args.bedfile
        ):
            self._parser.error(f"Cannot find input BED {self._args.bedfile}")
        if os.stat(self._args.bedfile).st_size <= 0:
            self._parser.error(f"{self._args.bdefile} is empty")
        # vcf folder
        if self._args.vcf and (not os.path.isdir(self._args.vcf)):
            self._parser.error(f"Cannot find VCF folder {self._args.vcf}")
        self._vcfs = glob(os.path.join(self._args.vcf, "*.vcf.gz"))
        if self._args.vcf and not self._vcfs:
            self._parser.error(f"No VCF file found in {self._args.vcf}")
        # pam
        if any(nt not in IUPAC for nt in self._args.pam):
            self._parser.error(f"PAM {self._args.pam} contains non IUPAC characters")
        # guide length
        if self._args.guidelen < 1:
            self._parser.error(f"Forbidden guide length ({self._args.guidelen})")
        # output folder
        if not os.path.exists(self._args.outdir) or not os.path.isdir(
            self._args.outdir
        ):
            self._parser.error(f"Cannot find output folder {self._args.outdir}")
        # verbosity
        if self._args.verbosity not in VERBOSITYLVL:
            self._parser.error(
                f"Forbidden verbosity level selected ({self._args.verbosity})"
            )

    @property
    def fasta(self) -> str:
        return self._args.fasta

    @property
    def fasta_idx(self) -> str:
        return self._args.fasta_idx

    @property
    def bedfile(self) -> str:
        return self._args.bedfile

    @property
    def vcfs(self) -> str:
        return self._vcfs

    @property
    def pam(self) -> str:
        return self._args.pam

    @property
    def guidelen(self) -> int:
        return self._args.guidelen

    @property
    def right(self) -> bool:
        return self._args.right

    @property
    def outdir(self) -> str:
        return self._args.outdir

    @property
    def no_filter(self) -> bool:
        return self._args.no_filter

    @property
    def verbosity(self) -> int:
        return self._args.verbosity

    @property
    def debug(self) -> bool:
        return self._args.debug
