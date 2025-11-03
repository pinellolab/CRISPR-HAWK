"""Argument parsing and validation for the CRISPR-HAWK command-line interface.

This module defines custom argument parsers and input argument handler classes for
the CRISPR-HAWK tool, supporting search, VCF conversion, data preparation, and
CRISPRitz configuration workflows. It ensures input consistency, provides helpful
error messages, and exposes validated arguments as convenient properties.
"""

from .utils import (
    warning,
    COMMAND,
    IUPAC,
    VERBOSITYLVL,
    TOOLNAME,
    OSSYSTEMS,
)
from .crisprhawk_version import __version__
from .config_crispritz import CrispritzConfig, check_crispritz_env

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

    def _validate_fasta_files(self) -> None:
        """Validates the existence and content of the input FASTA directory.

        This function checks that the specified FASTA directory exists and contains
        at least one FASTA file.

        Returns:
            None
        """
        if not os.path.exists(self._args.fasta) or not os.path.isdir(self._args.fasta):
            self._parser.error(f"Cannot find input FASTA folder {self._args.fasta}")
        self._fastas = glob(os.path.join(self._args.fasta, "*.fa")) + glob(
            os.path.join(self._args.fasta, "*.fasta")
        )
        if not self._fastas:
            self._parser.error(f"No FASTA file found in {self._args.fasta}")

    def _validate_bed(self) -> None:
        """Validates the existence and content of the input BED file.

        This function checks that the specified BED file exists and is not empty.

        Returns:
            None
        """
        if not os.path.exists(self._args.bedfile) or not os.path.isfile(
            self._args.bedfile
        ):
            self._parser.error(f"Cannot find input BED {self._args.bedfile}")
        if os.stat(self._args.bedfile).st_size <= 0:
            self._parser.error(f"{self._args.bdefile} is empty")

    def _validate_vcf_folder(self) -> None:
        """Validates the existence and content of the input VCF folder.

        This function checks that the specified VCF directory exists and contains
        at least one VCF file.

        Returns:
            None
        """
        if not self._args.vcf:
            self._vcfs = []
            return
        if not os.path.isdir(self._args.vcf):
            self._parser.error(f"Cannot find VCF folder {self._args.vcf}")
        self._vcfs = glob(os.path.join(self._args.vcf, "*.vcf.gz"))
        if not self._vcfs:
            self._parser.error(f"No VCF file found in {self._args.vcf}")

    def _validate_pam(self) -> None:
        """Validates the PAM sequence for allowed IUPAC nucleotide characters.

        This function checks that the PAM sequence contains only valid IUPAC characters.

        Returns:
            None
        """
        if any(nt not in IUPAC for nt in self._args.pam):
            self._parser.error(f"PAM {self._args.pam} contains non IUPAC characters")

    def _validate_guide_length(self) -> None:
        """Validates the guide length argument for minimum allowed value.

        This function checks that the guide length is at least 1.

        Returns:
            None
        """
        if self._args.guidelen < 1:
            self._parser.error(f"Forbidden guide length ({self._args.guidelen})")

    def _validate_output_folder(self) -> None:
        """Validates the existence of the output folder and creates it if necessary.

        This function checks that the specified output directory exists, creates
        it if missing, and sets its absolute path.

        Returns:
            None
        """
        if not os.path.exists(self._args.outdir) or not os.path.isdir(
            self._args.outdir
        ):
            if not os.path.isdir(
                os.path.dirname(self._args.outdir)
            ):  # parent doesn't exist
                self._parser.error(f"Cannot find output folder {self._args.outdir}")
            os.makedirs(self._args.outdir)  # create output folder
        self._outdir = os.path.abspath(self._args.outdir)
        assert os.path.isdir(self._outdir)

    def _validate_annotation_colnames(
        self, colnames: List[str], annotation_files: List[str], annotation_type: str
    ) -> None:
        """Validates the consistency between annotation column names and annotation files.

        This function checks that annotation column names are provided only when
        annotation files exist, and that the number of column names matches the
        number of annotation files.

        Args:
            colnames (List[str]): List of annotation column names.
            annotation_files (List[str]): List of annotation file paths.
            annotation_type (str): Type of annotation (e.g., "Annotation", "Gene
                Annotation").

        Returns:
            None
        """
        if colnames and not annotation_files:
            self._parser.error(
                f"{annotation_type} column names provided, but no input {annotation_type.lower()} file"
            )

        if colnames and len(colnames) != len(annotation_files):
            self._parser.error(
                f"Mismatching number of {annotation_type.lower()} files and {annotation_type.lower()} column names"
            )

    def _validate_annotations(self) -> None:
        """Validates the existence and content of input annotation BED files.

        This function checks that all specified annotation files exist, are not empty,
        and validates their column names.

        Returns:
            None
        """
        if not self._args.annotations:
            return  # no input annotation file
        if any(not os.path.isfile(f) for f in self._args.annotations):
            annfiles = ", ".join(self._args.annotations)
            self._parser.error(
                f"Cannot find the specified annotation BED files {annfiles}"
            )
        if any(os.stat(f).st_size <= 0 for f in self._args.annotations):
            annfiles = ", ".join(self._args.annotations)
            self._parser.error(f"{annfiles} look empty")
        # validate annotation colnames
        self._validate_annotation_colnames(
            self._args.annotation_colnames, self._args.annotations, "Annotation"
        )

    def _validate_gene_annotations(self) -> None:
        """Validates the existence and content of input gene annotation BED files.

        This function checks that all specified gene annotation files exist, are
        not empty, and validates their column names.

        Returns:
            None
        """
        if not self._args.gene_annotations:
            return  # no input gene annotation file
        if any(not os.path.isfile(f) for f in self._args.gene_annotations):
            annfiles = ", ".join(self._args.gene_annotations)
            self._parser.error(f"Cannot find gene annotation BED files {annfiles}")
        if any(os.stat(f).st_size <= 0 for f in self._args.gene_annotations):
            annfiles = ", ".join(self._args.gene_annotations)
            self._parser.error(f"{annfiles} look empty")
        self._validate_annotation_colnames(
            self._args.gene_annotation_colnames,
            self._args.gene_annotations,
            "Gene annotation",
        )

    def _validate_elevation_score(self) -> None:
        """Validates the input arguments for Elevation score computation.

        This function checks that the guide and PAM lengths sum to 23 and that
        the guide is downstream of the PAM.

        Returns:
            None
        """
        if not self._args.compute_elevation:
            return  # no elevation score request
        if self._args.guidelen + len(self._args.pam) != 23 or self._args.right:
            self._parser.error(
                "Elevation score requires that the combined length of the guide "
                "and PAM is exactly 23 bp, and that the guide sequence is located "
                "downstream of the PAM"
            )

    def _validate_offtargets_parameters(self) -> None:
        """Validates the input parameters for off-targets estimation.

        This function checks that the CRISPRitz genome index is provided and that
        the mismatch, DNA bulge, and RNA bulge arguments are non-negative.

        Returns:
            None
        """
        if not self._args.crispritz_index:  # check crispritz genome index
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

    def _validate_offtargets_annotations(self) -> None:
        """Validates the consistency of off-targets annotation arguments.

        This function checks that off-targets annotation is only requested when
        off-targets estimation is enabled, and validates the annotation column
        names and files.

        Returns:
            None
        """
        if not self._args.estimate_offtargets and (
            self._args.offtargets_annotation_colnames
            or self._args.offtargets_annotations
        ):
            self._parser.error(
                "Off-targets annotation requested, but missing off-targets estimation argument"
            )
        if self._args.estimate_offtargets:
            self._validate_annotation_colnames(
                self._args.offtargets_annotation_colnames,
                self._args.offtargets_annotations,
                "Off-targets annotation",
            )

    def _validate_offtargets_estimation(self) -> None:
        """Validates the input arguments and environment for off-targets estimation.

        This function checks that off-targets estimation is supported on the current
        system, initializes the CRISPRitz configuration, and validates annotation
        and parameter arguments.

        Returns:
            None
        """
        if self._args.estimate_offtargets and platform.system() != OSSYSTEMS[0]:
            warning(
                f"Off-target estimation is only supported on {OSSYSTEMS[0]} "
                "systems. Off-target estimation automatically disabled",
                1,
            )  # always disply this warning
            self._estimate_offtargets = False
            self._crispritz_config = None
            return  # skip off-targets estimation
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
        self._validate_offtargets_annotations()
        if self._args.estimate_offtargets:
            self._validate_offtargets_parameters()

    def _validate_candidate_guides(self) -> None:
        if any(len(g.split(":")) != 2 for g in self._args.candidate_guides):
            self._parser.error(
                "Candidate guides appear to not follow <chr>:<position> format"
            )
        if any(int(g.split(":")[1]) < 1 for g in self._args.candidate_guides):
            self._parser.error("Do candidate guides have negative position?")

    def _validate_threads(self) -> None:
        """Validates the thread count argument for allowed range.

        This function checks that the number of threads is non-negative and does
        not exceed the number of available CPU cores.

        Returns:
            None
        """
        max_threads = multiprocessing.cpu_count()
        if self._args.threads < 0 or self._args.threads > max_threads:
            self._parser.error(
                f"Forbidden number of threads provided ({self._args.threads}). "
                f"Max number of available cores: {max_threads}"
            )
        self._threads = max_threads if self._args.threads == 0 else self._args.threads

    def _validate_verbosity(self) -> None:
        """Validates the verbosity level argument for allowed values.

        This function checks that the verbosity level is one of the accepted values.

        Returns:
            None
        """
        if self._args.verbosity not in VERBOSITYLVL:
            self._parser.error(
                f"Forbidden verbosity level selected ({self._args.verbosity})"
            )

    def _check_consistency(self):
        """Checks the consistency and validity of all parsed input arguments.

        This function runs all validation routines for input files, parameters,
        and options, ensuring that the command-line arguments are correct and
        compatible for CRISPR-HAWK analysis.

        Returns:
            None
        """
        self._validate_fasta_files()  # check fasta file
        self._validate_bed()  # check bed file
        self._validate_vcf_folder()  # check vcf folder
        self._validate_pam()  # check pam
        self._validate_guide_length()  # check guide length
        self._validate_output_folder()  # check output folder
        self._validate_annotations()  # check functional annotation bed
        self._validate_gene_annotations()  # check gene annotation bed
        self._validate_elevation_score()  # check elevation score
        self._validate_offtargets_estimation()  # check off-targets estimation
        self._validate_candidate_guides()  # check candidate guides
        self._validate_threads()  # check threads number
        self._validate_verbosity()  # check verbosity

    @property
    def fastas(self) -> List[str]:
        return self._fastas

    @property
    def fastadir(self) -> str:
        return self._args.fasta

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
        return self._outdir

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
    def crispritz_config(self) -> Optional[CrispritzConfig]:
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
    def offtargets_annotations(self) -> List[str]:
        return self._args.offtargets_annotations

    @property
    def offtargets_annotation_colnames(self) -> List[str]:
        return self._args.offtargets_annotation_colnames

    @property
    def candidate_guides(self) -> List[str]:
        return self._args.candidate_guides

    @property
    def graphical_reports(self) -> bool:
        return self._args.graphical_reports

    @property
    def threads(self) -> int:
        return self._threads

    @property
    def verbosity(self) -> int:
        return self._args.verbosity

    @property
    def debug(self) -> bool:
        return self._args.debug


class CrisprHawkConverterInputArgs:
    """Handles and validates parsed command-line arguments for CRISPR-HAWK VCF
    conversion.

    This class checks the consistency of input arguments for VCF conversion and
    provides convenient access to validated argument values as properties.

    Attributes:
        _args (Namespace): The parsed arguments namespace.
        _parser (CrisprHawkArgumentParser): The argument parser instance.
    """

    def __init__(self, args: Namespace, parser: CrisprHawkArgumentParser) -> None:
        """Initializes the CrisprHawkConverterInputArgs with parsed arguments and
        parser.

        Stores the parsed arguments and parser, then checks argument consistency.
        """
        self._args = args
        self._parser = parser
        self._check_consistency()  # check input args consistency

    def _validate_gnomad_vcf_folder(self) -> None:
        """Validates the existence and content of the input gnomAD VCF folder.

        This function checks that the specified gnomAD VCF directory exists and
        contains at least one VCF file.

        Returns:
            None
        """
        if not os.path.isdir(self._args.gnomad_vcf_dir):
            self._parser.error(f"Cannot find VCF folder {self._args.gnomad_vcf_dir}")
        self._gnomad_vcfs = glob(
            os.path.join(self._args.gnomad_vcf_dir, "*.vcf.bgz")
        ) + glob(os.path.join(self._args.gnomad_vcf_dir, "*.vcf.gz"))
        if self._args.gnomad_vcf_dir and not self._gnomad_vcfs:
            self._parser.error(
                f"No gnomAD VCF file found in {self._args.gnomad_vcf_dir}"
            )

    def _validate_output_folder(self) -> None:
        """Validates the existence of the output folder and creates it if necessary.

        This function checks that the specified output directory exists, creates
        it if missing, and sets its absolute path.

        Returns:
            None
        """
        if not os.path.exists(self._args.outdir) or not os.path.isdir(
            self._args.outdir
        ):
            if not os.path.isdir(
                os.path.dirname(self._args.outdir)
            ):  # parent doesn't exist
                self._parser.error(f"Cannot find output folder {self._args.outdir}")
            os.makedirs(self._args.outdir)  # create output folder
        self._outdir = os.path.abspath(self._args.outdir)
        assert os.path.isdir(self._outdir)

    def _validate_threads(self) -> None:
        """Validates the thread count argument for allowed range.

        This function checks that the number of threads is non-negative and does
        not exceed the number of available CPU cores.

        Returns:
            None
        """
        max_threads = multiprocessing.cpu_count()
        if self._args.threads < 0 or self._args.threads > max_threads:
            self._parser.error(
                f"Forbidden number of threads provided ({self._args.threads}). "
                f"Max number of available cores: {max_threads}"
            )
        self._threads = max_threads if self._args.threads == 0 else self._args.threads

    def _validate_verbosity(self) -> None:
        """Validates the verbosity level argument for allowed values.

        This function checks that the verbosity level is one of the accepted values.

        Returns:
            None
        """
        if self._args.verbosity not in VERBOSITYLVL:
            self._parser.error(
                f"Forbidden verbosity level selected ({self._args.verbosity})"
            )

    def _check_consistency(self) -> None:
        """Checks the consistency and validity of all parsed input arguments for
        gnomAD VCF conversion.

        This function runs all validation routines for gnomAD VCF folder, output
        folder, thread count, and verbosity, ensuring that the command-line arguments
        are correct and compatible for CRISPR-HAWK VCF conversion analysis.

        Returns:
            None
        """
        self._validate_gnomad_vcf_folder()  # check gnomad vcf folder
        self._validate_output_folder()  # check output folder
        self._validate_threads()  # check threads number
        self._validate_verbosity()  # check verbosity

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
    """Handles and validates parsed command-line arguments for CRISPR-HAWK
    CRISPRme's input data preparation.

    This class checks the consistency of input arguments for data preparation and
    provides convenient access to validated argument values as properties.

    Attributes:
        _args (Namespace): The parsed arguments namespace.
        _parser (CrisprHawkArgumentParser): The argument parser instance.
    """

    def __init__(self, args: Namespace, parser: CrisprHawkArgumentParser) -> None:
        """Initialize CrisprHawkPrepareDataInputArgs with parsed arguments and
        parser.

        Stores the parsed arguments and parser, then checks argument consistency.

        Args:
            args (Namespace): The parsed arguments namespace.
            parser (CrisprHawkArgumentParser): The argument parser instance.
        """
        self._args = args
        self._parser = parser
        self._check_consistency()  # check input args consistency

    def _validate_report(self) -> None:
        """Validates the existence of the CRISPR-HAWK report file.

        This function checks that the specified report file exists before proceeding
        with data preparation.

        Returns:
            None
        """
        if self._args.report and (not os.path.isfile(self._args.report)):
            self._parser.error(f"Cannot find {TOOLNAME} report {self._args.report}")

    def _validate_output_folder(self) -> None:
        """Validates the existence of the output folder and creates it if necessary.

        This function checks that the specified output directory exists, creates
        it if missing, and sets its absolute path.

        Returns:
            None
        """
        if not os.path.exists(self._args.outdir) or not os.path.isdir(
            self._args.outdir
        ):
            if not os.path.isdir(
                os.path.dirname(self._args.outdir)
            ):  # parent doesn't exist
                self._parser.error(f"Cannot find output folder {self._args.outdir}")
            os.makedirs(self._args.outdir)  # create output folder
        self._outdir = os.path.abspath(self._args.outdir)
        assert os.path.isdir(self._outdir)

    def _check_consistency(self) -> None:
        """Checks the consistency and validity of all parsed input arguments for
        CRISPRme data preparation.

        This function runs all validation routines for the CRISPR-HAWK report and
        output folder, ensuring that the command-line arguments are correct and
        compatible for CRISPRme data preparation.

        Returns:
            None
        """
        self._validate_report()  # check crisprhawk report
        self._validate_output_folder()  # check output folder

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
    """Handles and validates parsed command-line arguments for CRISPR-HAWK Crispritz
    configuration.

    This class checks the consistency of input arguments for configuring Crispritz
    and provides convenient access to validated argument values as properties.

    Attributes:
        _args (Namespace): The parsed arguments namespace.
        _parser (CrisprHawkArgumentParser): The argument parser instance.
    """

    def __init__(self, args: Namespace, parser: CrisprHawkArgumentParser) -> None:
        """Initialize CrisprHawkCrispritzConfigInputArgs with parsed arguments and
        parser.

        Stores the parsed arguments and parser, then checks argument consistency.

        Args:
            args (Namespace): The parsed arguments namespace.
            parser (CrisprHawkArgumentParser): The argument parser instance.
        """
        self._args = args
        self._parser = parser
        self._check_consistency()  # check input args consistency

    def _validate_targets_dir(self) -> None:
        """Validates the existence of the CRISPRitz targets directory.

        This function checks that the specified targets directory exists before
        proceeding with configuration.

        Returns:
            None
        """
        if not self._args.targets_dir:
            return  # no targets folder specified, use default
        if not os.path.exists(self._args.targets_dir) and not os.path.isdir(
            self._args.targets_dir
        ):
            self._parser.error(
                f"Cannot find targets directory {self._args.targets_dir}"
            )

    def _validate_show_option(self) -> None:
        """Validates the --show option for CRISPRitz configuration argument parsing.

        This function checks that the --show option is not used in combination with
        other input arguments.

        Returns:
            None
        """
        if (
            self._args.env_name
            or self._args.targets_dir
            or self._args.reset
            or self._args.validate
        ) and self._args.show:
            self._parser.error(
                "--show options cannot be used with other input arguments"
            )

    def _validate_validate_option(self) -> None:
        """Validates the --validate option for CRISPRitz configuration argument
        parsing.

        This function checks that the --validate option is not used in combination
        with other input arguments.

        Returns:
            None
        """
        if (
            self._args.env_name
            or self._args.targets_dir
            or self._args.reset
            or self._args.show
        ) and self._args.validate:
            self._parser.error(
                "--validate options cannot be used with other input arguments"
            )

    def _validate_reset_option(self) -> None:
        """Validates the --reset option for CRISPRitz configuration argument parsing.

        This function checks that the --reset option is not used in combination with
        other input arguments.

        Returns:
            None
        """
        if (
            self._args.env_name
            or self._args.targets_dir
            or self._args.show
            or self._args.validate
        ) and self._args.reset:
            self._parser.error(
                "--reset options cannot be used with other input arguments"
            )

    def _check_consistency(self) -> None:
        """Checks the consistency and validity of all parsed input arguments for
        CRISPRitz configuration.

        This function runs all validation routines for the CRISPRitz targets directory
        and configuration options, ensuring that the command-line arguments are correct
        and compatible for CRISPRitz configuration management.

        Returns:
            None
        """
        self._validate_targets_dir()  # check crispritz config file
        self._validate_show_option()  # check show option
        self._validate_reset_option()  # check reset option
        self._validate_validate_option()  # check validate option

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
