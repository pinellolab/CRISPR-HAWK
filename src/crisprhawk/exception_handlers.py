"""Exception handling utilities for CRISPR-HAWK.

This module provides functions to handle system interrupts and exceptions in a 
consistent manner.

It includes handlers for SIGINT signals and for printing error messages and 
exiting the program gracefully.
"""

from typing import NoReturn, Optional
from colorama import init, Fore

import sys
import os


def sigint_handler() -> None:
    """Handle SIGINT (interrupt signal) to exit the program gracefully.

    Prints a message to standard error and exits the program with an OS error 
    code when SIGINT is received.
    """
    # print message when SIGINT is caught to exit gracefully from the execution
    sys.stderr.write(f"\nCaught SIGINT. Exit CRISPR-HAWK")
    sys.exit(os.EX_OSERR)  # mark as os error code


def exception_handler(
    exception_type: type,
    exception: str,
    code: int,
    debug: bool,
    e: Optional[Exception] = None,
) -> NoReturn:
    """Handle exceptions by printing an error message and exiting the program.

    Raises the specified exception with a formatted message in debug mode, or 
    prints an error and exits with the given code otherwise.

    Args:
        exception_type: The type of exception to raise.
        exception: The error message to display.
        code: The exit code to use when terminating the program.
        debug: Flag to enable debug mode for full stack trace.
        e: An optional previous exception to chain.

    Returns:
        This function does not return; it exits the program or raises an 
            exception.
    """
    init()  # initialize colorama render
    if debug:  # debug mode -> always trace back the full error stack
        if e:  # inherits from previous error
            raise exception_type(f"\n\n{exception}") from e
        raise exception_type(f"\n\n{exception}")  # divide exception message from stack
    # gracefully trigger error and exit execution
    sys.stderr.write(f"{Fore.RED}\n\nERROR: {exception}\n{Fore.RESET}")
    sys.exit(code)  # exit execution returning appropriate error code
