"""
"""

from utils import COMMAND
from crisprhawk_version import __version__  

from argparse import SUPPRESS, ArgumentParser, HelpFormatter, Action, _MutuallyExclusiveGroup
from typing import Iterable, Optional, TypeVar, Tuple, Dict, NoReturn
from colorama import Fore


import sys
import os

# define abstract generic types for typing 
_D = TypeVar("_D")
_V = TypeVar("_V")

class CrisprHawkArgumentParser(ArgumentParser):

    class CrisprHawkHelpFormatter(HelpFormatter):

        def add_usage(self, usage: str, actions: Iterable[Action], groups: Iterable[_MutuallyExclusiveGroup], prefix: Optional[str] = None) -> None:
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
            f"{Fore.RED}\nERROR: {error}.{Fore.RESET}" +
            f"\nRun {COMMAND} -h for usage\n\n"
        )
        sys.stderr.write(errormsg)  # write error to stderr
        sys.exit(os.EX_USAGE)  # exit execution -> usage error
        

    def error_noargs(self) -> None:
        self.print_help()  # if no input argument, print help
        sys.exit(os.EX_NOINPUT)  # exit with no input code 
