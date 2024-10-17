"""
"""

from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkError

from argparse import Namespace

import os

def crisprhawk(args: Namespace) -> None:
    print("congrats")

    exception_handler(CrisprHawkError, "error raised", os.EX_USAGE, False)