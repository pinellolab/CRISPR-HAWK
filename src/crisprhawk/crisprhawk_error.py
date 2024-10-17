""" 
"""

class CrisprHawkError(Exception):
    def __init__(self, value: str):
        # initialize exception object when raised
        self._value = value  # error message or error related info

    def __str__(self):
        return repr(self._value)  # string representation for the exception

 