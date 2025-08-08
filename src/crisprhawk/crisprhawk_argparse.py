""" """

from .utils import COMMAND, IUPAC, VERBOSITYLVL
from .crisprhawk_version import __version__

from argparse import (
    SUPPRESS,
    ArgumentParser,
    HelpFormatter,
    Action,
    _MutuallyExclusiveGroup,
    Namespace,
)
from typing import Iterable, Optional, TypeVar, Tuple, Dict, NoReturn, List
from colorama import Fore
from glob import glob


import sys
import os

# define abstract generic types for typing
_D = TypeVar("_D")
_V = TypeVar("_V")


class CrisprHawkArgumentParser(ArgumentParser):
    """Custom argument parser for CRISPR-HAWK command-line interface.

    This class extends argparse.ArgumentParser to provide custom help formatting,
    error handling, and version display for the CRISPR-HAWK tool.

    Attributes:
        usage (str): The usage string for the parser, with version information.
        formatter_class (type): The custom help formatter class.
    """

    class CrisprHawkHelpFormatter(HelpFormatter):
        """Custom help formatter for CRISPR-HAWK argument parser.

        This formatter customizes the usage message display for the help output.

        Attributes:
            None
        """

        def add_usage(  # type: ignore
            self,
            usage: str,
            actions: Iterable[Action],
            groups: Iterable[_MutuallyExclusiveGroup],
            prefix: Optional[str] = None,
        ) -> None:
            """Add a usage message to the help output.

            Displays the usage description unless suppressed.

            Args:
                usage (str): The usage string to display.
                actions (Iterable[Action]): The actions associated with the parser.
                groups (Iterable[_MutuallyExclusiveGroup]): Mutually exclusive
                    groups.
                prefix (Optional[str]): Optional prefix for the usage message.
            """
            # add usage description for help only if the set action is not to
            # suppress the display of the help formatter
            if usage != SUPPRESS:
                args = (usage, actions, groups, "")
                self._add_item(self._format_usage, args)  # initialize the formatter

    def __init__(self, *args: Tuple[_D], **kwargs: Dict[_D, _V]) -> None:
        """Initialize the CRISPR-HAWK argument parser.

        Sets up the parser with a custom help formatter and version display.

        Args:
            *args: Positional arguments for ArgumentParser.
            **kwargs: Keyword arguments for ArgumentParser.
        """
        # set custom help formatter defined as
        kwargs["formatter_class"] = self.CrisprHawkHelpFormatter  # type: ignore
        # replace the default version display in usage help with a custom
        # version display formatter
        kwargs["usage"] = kwargs["usage"].replace("{version}", __version__)  # type: ignore
        # initialize argument parser object with input parameters for
        # usage display
        super().__init__(*args, **kwargs)  # type: ignore

    def error(self, error: str) -> NoReturn:  # type: ignore
        """Display an error message and exit.

        Shows the error in red and suggests running the help command.

        Args:
            error (str): The error message to display.

        Raises:
            SystemExit: Exits the program with a usage error code.
        """
        # display error messages raised by argparse in red
        errormsg = (
            f"{Fore.RED}\nERROR: {error}.{Fore.RESET}"
            + f"\n\nRun {COMMAND} -h for usage\n\n"
        )
        sys.stderr.write(errormsg)  # write error to stderr
        sys.exit(os.EX_USAGE)  # exit execution -> usage error

    def error_noargs(self) -> None:
        """Display help and exit when no arguments are provided.

        Prints the help message and exits with a no input code.

        Raises:
            SystemExit: Exits the program with a no input error code.
        """
        self.print_help()  # if no input argument, print help
        sys.exit(os.EX_NOINPUT)  # exit with no input code


class CrisprHawkInputArgs:
    """Handles and validates parsed command-line arguments for CRISPR-HAWK.

    This class checks the consistency of input arguments and provides convenient
    access to validated argument values as properties.

    Attributes:
        _args (Namespace): The parsed arguments namespace.
        _parser (CrisprHawkArgumentParser): The argument parser instance.
    """

    def __init__(self, args: Namespace, parser: CrisprHawkArgumentParser) -> None:
        """Initialize CrisprHawkInputArgs with parsed arguments and parser.

        Stores the parsed arguments and parser, then checks argument consistency.

        Args:
            args (Namespace): The parsed arguments namespace.
            parser (CrisprHawkArgumentParser): The argument parser instance.
        """
        self._args = args
        self._parser = parser
        self._check_consistency()  # check input args consistency

    def _check_consistency(self):
        """Check the consistency and validity of parsed input arguments.

        Validates the existence, type, and content of input files and directories,
        and sets the list of VCF files found in the specified directory.

        Returns:
            None
        """
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
        # functional annotation bed
        if self._args.functional_annotation and (not os.path.exists(self._args.functional_annotation) or not os.path.isfile(self._args.functional_annotation)):
            self._parser.error(f"Cannot find functional annotation BED {self._args.functional_annotation}")
        if self._args.functional_annotation and os.stat(self._args.functional_annotation).st_size <= 0:
            self._parser.error(f"{self._args.functional_annotation} is empty")
        # gene annotation bed
        if self._args.gene_annotation and (not os.path.exists(self._args.gene_annotation) or not os.path.isfile(self._args.gene_annotation)):
            self._parser.error(f"Cannot find functional annotation BED {self._args.functional_annotation}")
        if self._args.gene_annotation and os.stat(self._args.gene_annotation).st_size <= 0:
            self._parser.error(f"{self._args.gene_annotation} is empty")
        # off-targets estimation
        if self._args.write_offtargets_report and not self._args.estimate_offtargets:
            self._parser.error("Cannot write off-targets report if off-target estimation not enabled")
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
    def vcfs(self) -> List[str]:
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
    def functional_annotation(self) -> str:
        return self._args.functional_annotation

    @property
    def gene_annotation(self) -> str:
        return self._args.gene_annotation

    @property
    def haplotype_table(self) -> bool:
        return self._args.haplotype_table
    
    @property
    def estimate_offtargets(self) -> bool:
        return self._args.estimate_offtargets
    
    @property
    def write_offtargets_report(self) -> bool:
        return self._args.write_offtargets_report

    @property
    def verbosity(self) -> int:
        return self._args.verbosity

    @property
    def debug(self) -> bool:
        return self._args.debug
