"""
"""

from crisprhawk_version import __version__  

from argparse import SUPPRESS, ArgumentParser, HelpFormatter, Action, _MutuallyExclusiveGroup
from typing import Iterable, Optional, Generic, TypeVar, Tuple, Dict, NoReturn
from colorama import Fore


import sys
import os

# define abstract generic types for typing 
_D = TypeVar("_D")
_V = TypeVar("_V")

class CRISPRHAWKArgumentParser(ArgumentParser):

    class CRISPRHAWKHelpFormatter(HelpFormatter):

        def add_usage(self, usage: str, actions: Iterable[Action], groups: Iterable[_MutuallyExclusiveGroup], prefix: Optional[str] = None) -> None:
            # add usage description for help only if the set action is not to
            # suppress the display of the help formatter
            if usage != SUPPRESS:  
                args = (usage, actions, groups, "")
                self._add_item(self._format_usage, args)  # initialize the formatter
        
    def __init__(self, *args: Tuple[Generic[_D]], **kwargs: Dict[Generic[_D], Generic[_V]]) -> None:
        # set custom help formatter defined as 
        kwargs["formatter_class"] = self.CRISPRHAWKHelpFormatter   
        # replace the default version display in usage help with a custom 
        # version display formatter
        kwargs["usage"] = kwargs["usage"].replace("{version}", __version__) 
        # initialize argument parser object with input parameters for
        # usage display
        super().__init__(*args, **kwargs)

    def error(self, error: str) -> NoReturn:
        # display error messages raised by argparse in red 
        error = Fore.RED + f"\nERROR: {error}." + Fore.RESET

    def error_noargs(self) -> None:
        self.print_help()  # if no input argument, print help
        sys.exit(os.EX_NOINPUT)  # exit with no input code 
