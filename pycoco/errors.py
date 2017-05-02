"""
Error definitions
"""

import numpy as np
import warnings

__all__ = ["CustomValueError", "PathError", "FilterMismatchError", "TableReadError", "StringWarning"]

##------------------------------------##
##  ERROR DEFS                        ##
##------------------------------------##


class CustomValueError(ValueError):
    """
    Raise when....
    """
    def __init__(self, *args, **kwargs):
        ValueError.__init__(self, *args, **kwargs)


class PathError(Exception):
    """
    Raise when a path is found to be invalid
    """
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


class FilterMismatchError(ValueError):
    """
    Raise when a Filter from filename doesn't match the one in the photfile
    """
    def __init__(self, *args, **kwargs):
        ValueError.__init__(self, *args, **kwargs)


class TableReadError(ValueError):
    """
    Raise when something goes wrong with the table I/O
    """
    def __init__(self, *args, **kwargs):
        ValueError.__init__(self, *args, **kwargs)


def StringWarning(path):
    """

    """
    if type(path) is not str and type(path) is not np.string_:
        warnings.warn("WARNING: You passed something that was " + str(type(path)) + "This might go wrong.",
                      stacklevel = 2)

    else:
        pass
