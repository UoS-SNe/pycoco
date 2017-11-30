#!/usr/bin env python

"""
This is the module for the CoCo python tools.

author: Rob Firth, Southampton
date: 06-12-2016
"""

import sys

if __name__ is not '__main__':

    __author__ = "RobFirth"
    __name__ = 'pycocosn'

try:
    __file__

except NameError:

    __file__ = sys.argv[0]

from . import extinction
from . import colours
from . import utils
from . import errors
from . import kcorr
from . import functions
from . import defaults
from . import classes
from . import models
from . import litdata
from . import coco

# from . import test
# test.test("foo")

if __name__ == "__main__":
    pass
