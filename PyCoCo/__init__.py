'''
This is the module for the PyCoCo python tools

author: Rob Firth, Southampton

date: 06-12-2016
'''

from __future__ import print_function ## Force python3-like printing

__name__ = 'PyCoCo'
__version__ = 0.1

try:
    __file__
except NameError:
    __file__ = sys.argv[0]

import os
import warnings

##----------------------------------------------------------------------------##
##                                   TOOLS                                    ##
##----------------------------------------------------------------------------##

##------------------------------------##
##  DUMMY CODE                        ##
##------------------------------------##

class CustomValueError(ValueError):
	"""
	Raise when....
	"""
	def __init__(self, *args, **kwargs):
		ValueError.__init__(self, *args, **kwargs)

class DummyClass():
    '''
    Quick dummy class.

    Contains a test class variable and test class method that prints the
    variable.

    RF
    '''

    def __init__(self):
        self.dummy_string = 'Hello, World!'

    def print_dummy_string(self):
        print(self.test_string)

def dummy_function(*args, **kwargs):
    '''
    Quick dummy function.

    Prints supplied **args and **kwargs
    Issues warnings if nothing passed
    
    RF
    '''

    warnings.simplefilter('always')
    print(args)
    print(kwargs)


    # warnings.warn("WARNING")

    if not args and not kwargs:
        warnings.warn( "You didn't pass any *args or **kwargs", RuntimeWarning)

    else:
        if args:
            for i, arg in enumerate(args):
                print('an arg passed via *args: ', repr(arg))
        else:
            warnings.warn( "You didn't pass any *args", RuntimeWarning)

        if kwargs:
            for key, value in kwargs.iteritems():
                print('a **kwarg: ', repr(key), ' == ' , repr(value))
        else:
            warnings.warn( "You didn't pass any **kwargs", RuntimeWarning)
    pass

##----------------------------------------------------------------------------##
##  CODE                                                                      ##
##----------------------------------------------------------------------------##

##------------------------------------##
##  ERROR DEFS                        ##
##------------------------------------##

class CustomValueError(ValueError):
	"""
	Raise when....
	"""
	def __init__(self, *args, **kwargs):
		ValueError.__init__(self, *args, **kwargs)

##------------------------------------##
##                                    ##
##------------------------------------##

class SNClass():
    """docstring for SNClass."""

    def __init__(self):
        self.data_directory = self.get_data_directory()

    def get_data_directory(self):
        """

        """

        return os.environ.get('PYCOCO_DATA_DIR', os.path.abspath(os.path.join(__file__, os.pardir)) + '/testdata/')




##----------------------------------------------------------------------------##
##  /CODE                                                                     ##
##----------------------------------------------------------------------------##
