""" """

from .crisprhawk_error import CrisprHawkConverterError
from .exception_handlers import exception_handler
from .utils import print_verbosity, remove_file, VERBOSITYLVL
from .variant import find_tbi, TBI

from pysam import VariantFile, VariantHeader, VariantRecord
from typing import List, Tuple
from functools import partial
from time import time

import multiprocessing
import pysam
import os


# define gnomAD populations (9 superpopulations and 1 collecting minority samples)
GNOMADPOPS = [
    "afr",
    "ami",
    "amr",
    "asj",
    "eas",
    "fin",
    "nfe",
    "mid",
    "sas",
    "remaining",
]

# define genotype header line
GTLINE = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Sample Collapsed Genotype">'

# genotype values (unphased)
GT = ["0/0", "0/1"]  # 0/0 -> absence, 0/1 -> occurrence


def tabix_index(vcf_fname: str, verbosity: int, debug: bool) -> str:
    """Indexes a VCF file using Tabix and returns the path to the index file.

    This function creates a Tabix index for the specified VCF file, checks for
    errors, and returns the index file path.

    Args:
        vcf_fname (str): Path to the VCF file to index.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        str: The path to the created Tabix index file.

    Raises:
        CrisprHawkConverterError: If indexing fails or the index file is not created.
    """
    print_verbosity(
        f"Indexing VCF {os.path.basename(vcf_fname)}", verbosity, VERBOSITYLVL[2]
    )
    start = time()
    pysam.tabix_index(vcf_fname, preset="vcf")  # index input vcf
    if not find_tbi(vcf_fname):  # error occurred while indexing file
        exception_handler(
            CrisprHawkConverterError,
            f"Failed indexing VCF {vcf_fname}",
            os.EX_DATAERR,
            debug,
        )
    print_verbosity(
        f"VCF {os.path.basename(vcf_fname)} indexed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return f"{vcf_fname}.{TBI}"


def load_vcf(vcf_fname: str, verbosity: int, debug: bool) -> VariantFile:
    """Loads a VCF file with an associated Tabix index.

    This function checks for an existing Tabix index for the VCF file, creates
    one if necessary, and returns a VariantFile object for reading the VCF.

    Args:
        vcf_fname (str): Path to the VCF file to load.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        VariantFile: A pysam VariantFile object for the indexed VCF file.
    """
    # search for tabix index, if not found index the current vcf
    tbi_index = (
        f"{vcf_fname}.{TBI}"
        if find_tbi(vcf_fname)
        else tabix_index(vcf_fname, verbosity, debug)
    )
    return VariantFile(vcf_fname, index_filename=tbi_index)


def format_ac(joint: bool) -> List[str]:
    """Constructs the allele count (AC) field names for gnomAD VCFs by population.

    Returns a list of formatted AC field names for each population, using either
    the joint or standard prefix.

    Args:
        joint (bool): Whether to use the joint AC field format.

    Returns:
        List[str]: A list of formatted AC field names for all gnomAD populations.
    """
    # construct the AC field format in gnomad vcf for each population
    return [f"AC_joint_{p}" if joint else f"AC_{p}" for p in GNOMADPOPS]


def variant_observed(allele_count: Tuple[int]) -> bool:
    """Determines if any allele in the count tuple is observed.

    Returns True if at least one allele count is greater than zero, otherwise False.

    Args:
        allele_count (Tuple[int]): A tuple of allele counts.

    Returns:
        bool: True if any count is greater than zero, False otherwise.
    """
    return any(ac > 0 for ac in allele_count)  # at least one count for allele


def _update_header(header: VariantHeader, joint: bool) -> str:
    """Updates the VCF header with genotype and population sample information.

    Adds the GT FORMAT field and population samples to the header, and adjusts
    the AF field if joint is True.

    Args:
        header (VariantHeader): The VCF header to update.
        joint (bool): Whether to adjust the AF field for joint VCFs.

    Returns:
        str: The updated VCF header as a string.
    """
    header.add_line(GTLINE)  # add FORMAT metadata field (genotype)
    header.add_samples(GNOMADPOPS)  # add populations as samples in header
    header = str(header).replace("<ID=AF_joint,", "<ID=AF,") if joint else str(header)  # type: ignore
    return header  # type: ignore


def _asses_genotype(
    variant: VariantRecord, ac_formatted: List[str], debug: bool
) -> str:
    """Assesses the genotype for each population based on allele counts in a variant.

    Returns a tab-separated string of genotype calls for each population, using
    GT[1] if the allele is observed and GT[0] otherwise.

    Args:
        variant (VariantRecord): The variant record containing allele count information.
        ac_formatted (List[str]): List of allele count field names for each population.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        str: Tab-separated genotype calls for each population.

    Raises:
        CrisprHawkConverterError: If genotype assessment fails due to missing or
            invalid data.
    """
    try:
        return "\t".join(
            [
                GT[1] if variant_observed(variant.info[acf]) else GT[0]
                for acf in ac_formatted
            ]
        )
    except Exception as e:
        exception_handler(
            CrisprHawkConverterError,
            f"Failed genotyoe assessment on variant {variant}",
            os.EX_DATAERR,
            debug,
            e,
        )


def _format_vrecord(variant: VariantRecord, genotypes: str) -> str:
    """Formats a VariantRecord and genotype string into a VCF record line.

    Constructs a tab-separated VCF record string from the variant's fields and
    provided genotype information, handling missing allele frequencies as needed.

    Args:
        variant (VariantRecord): The variant record to format.
        genotypes (str): The genotype string to include in the record.

    Returns:
        str: The formatted VCF record as a tab-separated string.
    """
    try:
        af = ",".join(list(map(str, variant.info["AF"])))
    except KeyError:  # catch potential AF missing in variant INFO
        af = ",".join(["0.0" for _ in variant.alts])  # type: ignore
    variant_format = [
        variant.chrom,
        variant.pos,
        variant.id,
        variant.ref,
        ",".join(variant.alts),  # handle multiallelic sites # type: ignore
        variant.qual,
        ";".join(variant.filter.keys()),  # type: ignore
        f"AF={af}",  # keep allele frequencies
        "GT",  # add genotype to format
        genotypes,
    ]
    return "\t".join([f"{e}" if e is not None else "." for e in variant_format])


def _convert(
    vcf: VariantFile,
    vcf_outfname: str,
    ac_formatted: List[str],
    joint: bool,
    keep: bool,
    debug: bool,
) -> None:
    """Converts a VCF file to a new format with updated header and genotype information.

    Writes a new VCF file with an updated header, filters variants based on the
    'keep' flag, and includes formatted genotype information for each variant.

    Args:
        vcf (VariantFile): The input VCF file to convert.
        vcf_outfname (str): The output filename for the converted VCF.
        ac_formatted (List[str]): List of allele count field names for each population.
        joint (bool): Whether the VCF is a joint file (affects header formatting).
        keep (bool): Whether to retain all variants or only those with 'PASS' in
            the filter.
        debug (bool): Flag to enable debug mode for error handling.

    Raises:
        CrisprHawkConverterError: If an I/O error occurs or the output VCF is empty.
    """
    try:
        with open(vcf_outfname, mode="w") as outfile:
            # write the updated header
            outfile.write(_update_header(vcf.header.copy(), joint))
            for variant in vcf:  # iterate over vcf variants
                if not keep and "PASS" not in variant.filter.keys():
                    continue
                genotypes = _asses_genotype(variant, ac_formatted, debug)
                outfile.write(f"{_format_vrecord(variant, genotypes)}\n")
    except IOError as e:
        exception_handler(
            CrisprHawkConverterError,
            f"An error occurred while converting {vcf_outfname}",
            os.EX_DATAERR,
            debug,
            e,
        )
    if os.stat(vcf_outfname).st_size <= 0:
        exception_handler(
            CrisprHawkConverterError,
            f"Empty converted VCF {vcf_outfname}",
            os.EX_DATAERR,
            debug,
        )


def _compress(vcf_fname_tmp: str, vcf_fname: str, verbosity: int, debug: bool) -> None:
    """Compresses a VCF file using bgzip and removes the temporary file.

    This function compresses the specified VCF file with bgzip, handles errors,
    removes the temporary file, and logs the operation.

    Args:
        vcf_fname_tmp (str): Path to the temporary VCF file to compress.
        vcf_fname (str): Path to the output compressed VCF file.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Raises:
        CrisprHawkConverterError: If bgzip compression fails.
    """
    print_verbosity(
        f"Compressing {os.path.basename(vcf_fname)}", verbosity, VERBOSITYLVL[2]
    )
    start = time()
    try:  # compress converted vcf using bgzip (compliant with tabix rules)
        pysam.tabix_compress(vcf_fname_tmp, vcf_fname, force=True)
    except OSError as e:
        exception_handler(
            CrisprHawkConverterError,
            f"BGZIP compression failed for {os.path.basename(vcf_fname)}",
            os.EX_DATAERR,
            debug,
            e,
        )
    remove_file(vcf_fname_tmp)  # remove temporary converted vcf
    print_verbosity(
        f"VCF {os.path.basename(vcf_fname)} compressed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )


def convert_vcf(
    vcf_fname: str,
    joint: bool,
    keep: bool,
    suffix: str,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> None:
    """Converts a VCF file to a CRISPR-HAWK compatible format and compresses it.

    This function loads a VCF file, reformats it with updated genotype and population
    information, compresses the result, and logs the process.

    Args:
        vcf_fname (str): Path to the input VCF file.
        joint (bool): Whether the VCF is a joint file (affects AC field formatting).
        keep (bool): Whether to retain all variants or only those with 'PASS' in
            the filter.
        suffix (str): Suffix to append to the output file name.
        outdir (str): Output directory for the converted VCF.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.
    """
    print_verbosity(
        f"Converting VCF {os.path.basename(vcf_fname)}", verbosity, VERBOSITYLVL[2]
    )
    start = time()
    try:  # optimized vcf loading through pysam
        vcf = load_vcf(vcf_fname, verbosity, debug)
    except OSError as e:
        exception_handler(
            CrisprHawkConverterError,
            f"Failed loading VCF {vcf_fname}",
            os.EX_DATAERR,
            debug,
            e,
        )
    ac_formatted = format_ac(joint)  # guess AC format in vcfs
    vcf_outfname_prefix = os.path.join(
        outdir,
        f"{os.path.splitext(os.path.splitext(os.path.basename(vcf_fname))[0])[0]}.{suffix}",
    )
    vcf_outfname_tmp = f"{vcf_outfname_prefix}.tmp.vcf"
    _convert(
        vcf, vcf_outfname_tmp, ac_formatted, joint, keep, debug
    )  # run vcf conversion
    vcf_outfname = f"{vcf_outfname_prefix}.vcf.gz"
    _compress(
        vcf_outfname_tmp, vcf_outfname, verbosity, debug
    )  # compress converted vcf
    print_verbosity(
        f"{os.path.basename(vcf_fname)} converted in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )


def convert_gnomad_vcf(
    gnomad_vcfs: List[str],
    joint: bool,
    keep: bool,
    suffix: str,
    outdir: str,
    threads: int,
    verbosity: int,
    debug: bool,
) -> None:
    """Converts multiple gnomAD VCF files to CRISPR-HAWK compatible format in parallel.

    This function uses multiprocessing to convert a list of gnomAD VCF files, applying
    the specified options to each file and handling errors as needed.

    Args:
        gnomad_vcfs (List[str]): List of input gnomAD VCF file paths.
        joint (bool): Whether the VCFs are joint files (affects AC field formatting).
        keep (bool): Whether to retain all variants or only those with 'PASS' in the
            filter.
        suffix (str): Suffix to append to each output file name.
        outdir (str): Output directory for the converted VCFs.
        threads (int): Number of parallel processes to use.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Raises:
        CrisprHawkConverterError: If an error occurs during the conversion process.
    """
    pool = multiprocessing.Pool(processes=threads)  # create processes pool
    try:
        convert_vcf_ = partial(
            convert_vcf,
            joint=joint,
            keep=keep,
            suffix=suffix,
            outdir=outdir,
            verbosity=verbosity,
            debug=debug,
        )
        pool.map(convert_vcf_, gnomad_vcfs)  # map vcfs to processes
        pool.close()
        pool.join()
    except OSError as e:
        exception_handler(
            CrisprHawkConverterError,
            "GnomAD VCF conversion failed",
            os.EX_DATAERR,
            debug,
            e,
        )
