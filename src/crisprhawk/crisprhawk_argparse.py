""" """

from .utils import (
    warning,
    command_exists,
    COMMAND,
    IUPAC,
    VERBOSITYLVL,
    TOOLNAME,
    OSSYSTEMS,
)
from .crisprhawk_version import __version__
from .config_crispritz import CrispritzConfig, check_crispritz_env, MAMBA

from argparse import (
    SUPPRESS,
    ArgumentParser,
    HelpFormatter,
    Action,
    _MutuallyExclusiveGroup,
    Namespace,
)
from typing import Iterable, Optional, TypeVar, Tuple, Dict, NoReturn, List, Union
from colorama import Fore
from glob import glob

import multiprocessing
import platform
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
        if "usage" in kwargs:
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


class CrisprHawkSearchInputArgs:
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
        if not os.path.exists(self._args.fasta) or not os.path.isdir(self._args.fasta):
            self._parser.error(f"Cannot find input FASTA folder {self._args.fasta}")
        self._fastas = glob(os.path.join(self._args.fasta, "*.fa")) + glob(
            os.path.join(self._args.fasta, "*.fasta")
        )
        if not self._fastas:
            self._parser.error(f"No FASTA file found in {self._args.fasta}")
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
        if self._args.annotations and (
            any(not os.path.isfile(f) for f in self._args.annotations)
        ):
            annfiles = ", ".join(self._args.annotations)
            self._parser.error(
                f"Cannot find the specified annotation BED files {annfiles}"
            )
        if self._args.annotations and any(
            os.stat(f).st_size <= 0 for f in self._args.annotations
        ):
            annfiles = ", ".join(self._args.annotations)
            self._parser.error(f"{annfiles} look empty")
        # functional annotation colnames
        if self._args.annotation_colnames and not self._args.annotations:
            self._parser.error(
                "Annotation column names provided, but no input annotation file"
            )
        if self._args.annotation_colnames and (
            len(self._args.annotation_colnames) != len(self._args.annotations)
        ):
            self._parser.error(
                f"Mismatching number of annotation files and annotation column names"
            )
        # gene annotation bed
        if self._args.gene_annotations and (
            any(not os.path.isfile(f) for f in self._args.gene_annotations)
        ):
            annfiles = ", ".join(self._args.gene_annotations)
            self._parser.error(f"Cannot find gene annotation BED files {annfiles}")
        if self._args.gene_annotations and any(
            os.stat(f).st_size <= 0 for f in self._args.gene_annotations
        ):
            annfiles = ", ".join(self._args.gene_annotations)
            self._parser.error(f"{annfiles} look empty")
        # gene annotation colnames
        if self._args.gene_annotation_colnames and not self._args.gene_annotations:
            self._parser.error(
                "Gene annotation column names provided, but no input gene annotation file"
            )
        if self._args.gene_annotation_colnames and (
            len(self._args.gene_annotation_colnames) != len(self._args.gene_annotations)
        ):
            self._parser.error(
                f"Mismatching number of gene annotation files and gene annotation column names"
            )
        # elevation score
        if self._args.compute_elevation and (
            self._args.guidelen + len(self._args.pam) != 23 or self._args.right
        ):
            self._parser.error(
                "Elevation score requires that the combined length of the guide "
                "and PAM is exactly 23 bp, and that the guide sequence is located "
                "downstream of the PAM"
            )
        # off-targets estimation
        if self._args.estimate_offtargets and platform.system() != OSSYSTEMS[0]:
            warning(
                f"Off-target estimation is only supported on {OSSYSTEMS[0]} "
                "systems. Off-target estimation automatically disabled",
                1,
            )  # always disply this warning
            self._estimate_offtargets = False
            self._crispritz_config = None
        if self._args.estimate_offtargets:
            self._estimate_offtargets = self._args.estimate_offtargets
            self._crispritz_config = CrispritzConfig()  # read crispritz config
            if not self._crispritz_config.set_command() or not check_crispritz_env(
                self._crispritz_config.env_name, self._crispritz_config.conda
            ):  # check if mamba/conda and crispritz environment are available
                self._estimate_offtargets = False
                self._crispritz_config = None
        else:
            self._estimate_offtargets = False
            self._crispritz_config = None
        # crispritz genome index
        if self._args.estimate_offtargets and not self._args.crispritz_index:
            self._parser.error("Genome index required for off-targets estimation")
        # check mm, bdna and brna arguments
        if self._args.estimate_offtargets and self._args.mm < 0:
            self._parser.error(f"Forbidden number of mismatches given: {self._args.mm}")
        if self._args.estimate_offtargets and self._args.bdna < 0:
            self._parser.error(
                f"Forbidden number of DNA bulges given: {self._args.bdna}"
            )
        if self._args.estimate_offtargets and self._args.brna < 0:
            self._parser.error(
                f"Forbidden number of RNA bulges given: {self._args.brna}"
            )
        # threads number
        if self._args.threads < 0 or self._args.threads > multiprocessing.cpu_count():
            self._parser.error(
                f"Forbidden number of threads provided ({self._args.threads}). Max number of available cores: {multiprocessing.cpu_count()}"
            )
        if self._args.threads == 0:  # use all cores
            self._args.threads = multiprocessing.cpu_count()
        # verbosity
        if self._args.verbosity not in VERBOSITYLVL:
            self._parser.error(
                f"Forbidden verbosity level selected ({self._args.verbosity})"
            )

    @property
    def fastas(self) -> List[str]:
        return self._fastas

    @property
    def fastadir(self) -> str:
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
    def annotations(self) -> List[str]:
        return self._args.annotations

    @property
    def annotation_colnames(self) -> List[str]:
        return self._args.annotation_colnames

    @property
    def gene_annotations(self) -> List[str]:
        return self._args.gene_annotations

    @property
    def gene_annotation_colnames(self) -> List[str]:
        return self._args.gene_annotation_colnames

    @property
    def haplotype_table(self) -> bool:
        return self._args.haplotype_table

    @property
    def compute_elevation(self) -> bool:
        return self._args.compute_elevation

    @property
    def estimate_offtargets(self) -> bool:
        return self._estimate_offtargets

    @property
    def crispritz_config(self) -> Union[None, CrispritzConfig]:
        return self._crispritz_config

    @property
    def crispritz_index(self) -> str:
        return self._args.crispritz_index

    @property
    def mm(self) -> int:
        return self._args.mm

    @property
    def bdna(self) -> int:
        return self._args.bdna

    @property
    def brna(self) -> int:
        return self._args.brna

    @property
    def graphical_reports(self) -> bool:
        return self._args.graphical_reports

    @property
    def threads(self) -> int:
        return self._args.threads

    @property
    def verbosity(self) -> int:
        return self._args.verbosity

    @property
    def debug(self) -> bool:
        return self._args.debug


class CrisprHawkConverterInputArgs:

    def __init__(self, args: Namespace, parser: CrisprHawkArgumentParser) -> None:
        self._args = args
        self._parser = parser
        self._check_consistency()  # check input args consistency

    def _check_consistency(self):
        # vcf folder
        if self._args.gnomad_vcf_dir and (not os.path.isdir(self._args.gnomad_vcf_dir)):
            self._parser.error(f"Cannot find VCF folder {self._args.gnomad_vcf_dir}")
        self._gnomad_vcfs = glob(
            os.path.join(self._args.gnomad_vcf_dir, "*.vcf.bgz")
        ) + glob(os.path.join(self._args.gnomad_vcf_dir, "*.vcf.gz"))
        if self._args.gnomad_vcf_dir and not self._gnomad_vcfs:
            self._parser.error(
                f"No gnomAD VCF file found in {self._args.gnomad_vcf_dir}"
            )
        # output folder
        if not os.path.exists(self._args.outdir) or not os.path.isdir(
            self._args.outdir
        ):
            self._parser.error(f"Cannot find output folder {self._args.outdir}")
        # threads number
        if self._args.threads < 0 or self._args.threads > multiprocessing.cpu_count():
            self._parser.error(
                f"Forbidden number of threads provided ({self._args.threads}). Max number of available cores: {multiprocessing.cpu_count()}"
            )
        if self._args.threads == 0:  # use all cores
            self._args.threads = multiprocessing.cpu_count()
        # verbosity
        if self._args.verbosity not in VERBOSITYLVL:
            self._parser.error(
                f"Forbidden verbosity level selected ({self._args.verbosity})"
            )

    @property
    def gnomad_vcfs(self) -> List[str]:
        return self._gnomad_vcfs

    @property
    def outdir(self) -> str:
        return self._args.outdir

    @property
    def joint(self) -> bool:
        return self._args.joint

    @property
    def keep(self) -> bool:
        return self._args.keep

    @property
    def suffix(self) -> str:
        return self._args.suffix

    @property
    def threads(self) -> int:
        return self._args.threads

    @property
    def verbosity(self) -> int:
        return self._args.verbosity

    @property
    def debug(self) -> bool:
        return self._args.debug


class CrisprHawkPrepareDataInputArgs:

    def __init__(self, args: Namespace, parser: CrisprHawkArgumentParser) -> None:
        self._args = args
        self._parser = parser
        self._check_consistency()  # check input args consistency

    def _check_consistency(self):
        # crisprhawk report
        if self._args.report and (not os.path.isfile(self._args.report)):
            self._parser.error(f"Cannot find {TOOLNAME} report {self._args.report}")
        # output folder
        if not os.path.exists(self._args.outdir) or not os.path.isdir(
            self._args.outdir
        ):
            self._parser.error(f"Cannot find output folder {self._args.outdir}")

    @property
    def report(self) -> str:
        return self._args.report

    @property
    def create_pam(self) -> bool:
        return self._args.create_pam

    @property
    def outdir(self) -> str:
        return self._args.outdir

    @property
    def debug(self) -> bool:
        return self._args.debug


class CrisprHawkCrispritzConfigInputArgs:

    def __init__(self, args: Namespace, parser: CrisprHawkArgumentParser) -> None:
        self._args = args
        self._parser = parser
        self._check_consistency()  # check input args consistency

    def _check_consistency(self):
        # crispritz config file
        if self._args.targets_dir:
            if not os.path.exists(self._args.targets_dir) and not os.path.isdir(
                self._args.targets_dir
            ):
                self._parser.error(
                    f"Cannot find targets directory {self._args.targets_dir}"
                )
        # show option
        if (
            self._args.env_name
            or self._args.targets_dir
            or self._args.reset
            or self._args.validate
        ) and self._args.show:
            self._parser.error(
                f"--show options cannot be used with other input arguments"
            )
        # reset option
        if (
            self._args.env_name
            or self._args.targets_dir
            or self._args.show
            or self._args.validate
        ) and self._args.reset:
            self._parser.error(
                f"--reset options cannot be used with other input arguments"
            )
        # validate option
        if (
            self._args.env_name
            or self._args.targets_dir
            or self._args.reset
            or self._args.show
        ) and self._args.validate:
            self._parser.error(
                f"--validate options cannot be used with other input arguments"
            )

    @property
    def env_name(self) -> str:
        return self._args.env_name

    @property
    def targets_dir(self) -> str:
        return self._args.targets_dir

    @property
    def show(self) -> bool:
        return self._args.show

    @property
    def reset(self) -> bool:
        return self._args.reset

    @property
    def validate(self) -> bool:
        return self._args.validate
