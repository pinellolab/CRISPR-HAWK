"""
"""

from colorama import Fore

import sys

# define verbosity levels
VERBOSITYLVL = [0, 1, 2, 3]

def warning(message: str, verbosity: int) -> None:
    if verbosity >= VERBOSITYLVL[1]:
        sys.stderr.write(f"{Fore.YELLOW}WARNING: {message}.{Fore.RESET}\n")
    return