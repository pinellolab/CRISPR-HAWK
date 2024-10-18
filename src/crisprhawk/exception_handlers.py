"""
"""

from utils import TOOLNAME

from typing import NoReturn, Optional
from colorama import init, Fore

import sys
import os


def sigint_handler() -> None:
    # print message when SIGINT is caught to exit gracefully from the execution
    sys.stderr.write(f"\nCaught SIGINT. Exit {TOOLNAME}")
    sys.exit(os.EX_OSERR)  # mark as os error code


def exception_handler(
    exception_type: Exception,
    exception: str,
    code: int,
    debug: bool,
    e: Optional[Exception] = None,
) -> NoReturn:
    init()  # initialize colorama render
    if debug:  # debug mode -> always trace back the full error stack
        if e:  # inherits from previous error
            raise exception_type(f"\n\n{exception}") from e
        raise exception_type(f"\n\n{exception}")  # divide exception message from stack
    # gracefully trigger error and exit execution
    sys.stderr.write(f"{Fore.RED}\n\nERROR: {exception}\n{Fore.RESET}")
    sys.exit(code)  # exit execution returning appropriate error code
