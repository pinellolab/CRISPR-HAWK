"""Utility functions for configuring CRISPR-HAWK environments and models.

This module provides helpers to locate configuration files, manage conda/mamba
environments, and download and prepare scoring model data required by
CRISPR-HAWK.
"""

from .exception_handlers import exception_handler
from .utils import OSSYSTEMS, TOOLNAME, remove_file_silent, warning, rename_file

from typing import List, Optional

import shutil
import subprocess
import sys
import os
import urllib.request
import urllib.error
import zipfile


# ==============================================================================
# Define constant variables (config files)
# ==============================================================================

# config file location
CONFIG = os.path.abspath(os.path.join(os.path.dirname(__file__), "config/config.json"))

# mamba/conda commands
MAMBA = "mamba"
CONDA = "conda"

# define crispr-hawk scoring models catalogue
_MODEL_CATALOGUE = [
    (
        "https://zenodo.org/records/19680449/files/azimuth.zip",
        "azimuth",
        "saved_models",
    ),
    (
        "https://zenodo.org/records/19680449/files/cfdscore.zip",
        "cfdscore",
        "models",
    ),
    (
        "https://zenodo.org/records/19680449/files/crispron.zip",
        "crispron",
        "data",
    ),
    (
        "https://zenodo.org/records/19680449/files/deepCpf1.zip",
        "deepCpf1",
        "weights",
    ),
    (
        "https://zenodo.org/records/19680449/files/elevation.zip",
        "elevation",
        "models",
    ),
    (
        "https://zenodo.org/records/19680449/files/plm_crispr.zip",
        "plm_crispr",
        "models",
    ),
    (
        "https://zenodo.org/records/19680449/files/sgdesigner.zip",
        "sgdesigner",
        "models",
    ),
]

# define download tuning constants
_CONNECT_TIMEOUT = 30  # seconds to wait for the TCP connection
_READ_TIMEOUT = 120  # seconds to wait between data chunks
_MAX_RETRIES = 3  # total attempts per file
_RETRY_BACKOFF = 5  # seconds to wait before first retry (doubles each time)
_CHUNK_SIZE = 1024 * 256  # 256 KiB read chunks


# ==============================================================================
# Define utilities functions (config files)
# ==============================================================================


def command_exists(command: str) -> bool:
    """Check if a command exists in the system's PATH.

    Returns True if the specified command is found in the system's executable
    search path, otherwise False.

    Args:
        command (str): The command to check for existence.

    Returns:
        bool: True if the command exists, False otherwise.
    """
    return bool(shutil.which(command))


def set_command() -> str:
    """Select the preferred package manager command for CRISPR-HAWK.

    This function checks for the availability of mamba and conda in the system
    PATH and returns the most suitable executable name to use, or an empty
    string if neither is found.

    Returns:
        str: The command name to use ('mamba', 'conda', or an empty string
            if neither are available).
    """
    if command_exists(MAMBA):
        return MAMBA
    return CONDA if command_exists(CONDA) else ""


def create_mamba_env(
    conda: str,
    env_name: str,
    packages: Optional[List[str]] = None,
    python_version: str = "3.8",
) -> bool:
    """Create a new conda or mamba environment for CRISPR-HAWK tools.

    This function invokes the specified package manager to create an environment
    with a given Python version and an optional list of additional packages.

    Args:
        conda: The command to invoke (typically 'mamba' or 'conda').
        env_name: Name of the environment to create.
        packages: Optional list of additional packages to install in the environment.
        python_version: Python version string to use for the environment.

    Returns:
        bool: True if the environment creation command succeeded, otherwise False.
    """
    cmd = [conda, "create", "-y", "-n", env_name, f"python={python_version}"]
    if packages:
        cmd.extend(packages)
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def _scores_dir() -> str:
    """Return the absolute path to the CRISPR-HAWK scores directory.

    This helper resolves the scores subdirectory relative to the current
    module location for use when locating or installing model data.

    Returns:
        str: Absolute path to the 'scores' directory.
    """
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), "scores")


def _sentinel_present(destdir: str, sentinel: str) -> bool:
    """Check whether a sentinel file or directory exists in a destination folder.

    This helper is used to determine if model data has already been installed
    by looking for a specific marker path inside the target directory.

    Args:
        destdir: Base directory in which to search for the sentinel.
        sentinel: File or directory name that signals the presence of model data.

    Returns:
        bool: True if the sentinel path exists as a file or directory, otherwise False.
    """
    full = os.path.join(destdir, sentinel)
    return os.path.isdir(full) or os.path.isfile(full)


def _download_with_progress(url: str, dest_path: str, label: str) -> None:
    """Download a file from a URL while reporting progress to stderr.

    This helper streams data in chunks, writes them to disk, and prints
    approximate download progress based on the reported content length.

    Args:
        url: The URL from which to download the file.
        dest_path: Local filesystem path where the downloaded file will be saved.
        label: Human-readable label used in progress messages.
    """
    req = urllib.request.Request(url, headers={"User-Agent": "crisprhawk/1.0"})
    with urllib.request.urlopen(req, timeout=_CONNECT_TIMEOUT) as response:
        total = response.headers.get("Content-Length")
        total_mb = f"{int(total) / 1048576:.1f} MB" if total else "unknown size"
        sys.stderr.write(f"Downloading {label} ({total_mb})...\n")
        downloaded = 0
        with open(dest_path, mode="wb") as fh:
            while chunk := response.read(_CHUNK_SIZE):
                fh.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded * 100 // int(total)
                    sys.stderr.write(f"\rProgress: {pct:3d}%")
                    sys.stderr.flush()
        if total:
            sys.stderr.write("\n")


def _download_with_retry(url: str, dest_path: str, label: str) -> None:
    """Download a file with automatic retries and error handling.

    This helper repeatedly attempts to download a file, applying a backoff
    delay between failures and raising a fatal error if all attempts fail.

    Args:
        url: The URL from which to download the file.
        dest_path: Local filesystem path where the downloaded file will be saved.
        label: Human-readable label used in log and warning messages.
    """
    delay = _RETRY_BACKOFF
    exc: Exception = RuntimeError("No attempts made")
    for attempt in range(1, _MAX_RETRIES + 1):
        try:
            _download_with_progress(url, dest_path, label)
            return  # success
        except (
            urllib.error.URLError,
            urllib.error.HTTPError,
            OSError,
            TimeoutError,
            ConnectionError,
        ) as e:
            exc = e
            if attempt < _MAX_RETRIES:
                warning(
                    f"Download attempt {attempt}/{_MAX_RETRIES} for {label} failed"
                    f"({exc}). Retrying in {delay}s",
                    1,
                )
            else:
                warning(
                    f"Download attempt {attempt}/{_MAX_RETRIES} for {label} failed "
                    f"({exc}). No more retries",
                    1,
                )
        # remove partially-written file so a future run starts clean
        if os.path.isfile(dest_path):
            remove_file_silent(dest_path)
        exception_handler(
            RuntimeError,
            f"Failed to download {label} after {_MAX_RETRIES} attempts. Last error: {exc}",
            os.EX_DATAERR,
            True,
        )


def _extract_zip(zip_path: str, destdir: str, label: str) -> None:
    """Extract a ZIP archive containing model data into the target directory.

    This helper unpacks the archive, logs progress, and cleans up the ZIP file
    even if extraction fails so that future runs start from a clean state.

    Args:
        zip_path: Path to the ZIP file to extract.
        destdir: Destination directory where the archive contents will be unpacked.
        label: Human-readable label used in log and error messages.
    """
    sys.stderr.write(f"Extracting {label}...\n")
    try:
        with zipfile.ZipFile(zip_path, mode="r") as zf:
            zf.extractall(path=destdir)
    except (zipfile.BadZipFile, RuntimeError) as e:
        exception_handler(
            zipfile.BadZipFile,
            f"Extraction of {zip_path} failed",
            os.EX_DATAERR,
            True,
            e,
        )
    finally:
        # always remove the zip, even when extraction fails, so that a re-run
        # does not try to extract an incomplete archive
        if os.path.isfile(zip_path):
            remove_file_silent(zip_path)


def _ensure_model(scoresdir: str, url: str, subdir: str, sentinel: str) -> None:
    """Ensure that scoring model data for a given submodule are available.

    This helper checks for existing model files, downloads and extracts them
    if necessary, and validates that the expected sentinel resource is present.

    Args:
        scoresdir: Base directory where all scoring model subdirectories reside.
        url: Remote URL from which to download the model archive if needed.
        subdir: Name of the scoring model subdirectory to inspect and populate.
        sentinel: File or directory name used to confirm that the model data
            have been correctly installed.
    """
    destdir = os.path.join(scoresdir, subdir)
    # the sub-folder itself must exist it contains __init__.py etc.
    if not os.path.isdir(destdir):
        exception_handler(
            FileNotFoundError,
            f"Cannot find score module directory: {destdir}. Please re-install {TOOLNAME}",
            os.EX_DATAERR,
            True,
        )
    # nothing to do if the model data are already in place
    if _sentinel_present(destdir, sentinel):
        return
    zip_path = os.path.join(destdir, f"{subdir}.zip")
    warning(
        f"Model data for '{subdir}' not found. Downloading model. "
        "This may take a while",
        1,
    )
    _download_with_retry(url, zip_path, subdir)
    _extract_zip(zip_path, destdir, subdir)
    orig = os.path.join(destdir, subdir, sentinel)
    dest = os.path.join(destdir, sentinel)
    rename_file(orig, dest)
    # final sanity-check
    if not _sentinel_present(destdir, sentinel):
        exception_handler(
            RuntimeError,
            f"Model dat for '{destdir}' still missing after extraction. Expected '{sentinel}' inside {destdir}",
            os.EX_DATAERR,
            True,
        )
    sys.stderr.write(f"'{destdir}' model data ready\n")


def prepare_scoring_models() -> None:
    """Ensure that all CRISPR-HAWK scoring models are downloaded and ready to use.

    This function verifies the presence of model data for each supported scoring
    module and automatically downloads and installs any missing resources.

    """
    scores_dir = _scores_dir()  # retrieve scores location
    if not os.path.isdir(scores_dir):
        exception_handler(
            FileNotFoundError,
            f"Scores directory not found: {scores_dir}. Please re-install {TOOLNAME}",
            os.EX_DATAERR,
            True,
        )
    for url, subdir, sentinel in _MODEL_CATALOGUE:
        _ensure_model(scores_dir, url, subdir, sentinel)
