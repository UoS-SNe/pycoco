#!/usr/bin env python

"""
This is the module for the CoCo python tools.

author: Rob Firth, Southampton
date: 06-12-2016
"""

import sys

if __name__ is not '__main__':

    __name__ = 'pycoco'

try:
    __file__

except NameError:

    __file__ = sys.argv[0]

# from .extinction import *
# from .colours import *
# from .utils import *
# from .errors import *
# from .kcorr import *
# from .functions import *
# from .defaults import *
# from .classes import *
# from .models import *
# from .litdata import *
# from .coco_calls import *

## TODO - implement following structure.
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
from . import coco_calls

# from . import test
# test.test("foo")

if __name__ == "__main__":
    pass
