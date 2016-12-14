'''This is the module for the PyCoCo python tools.

author: Rob Firth, Southampton
date: 06-12-2016'''

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
##  TESTING                           ##
##------------------------------------##



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

class PathError(StandardError):
	"""
	Raise when a path is found to be invalid
	"""
	def __init__(self, *args, **kwargs):
		StandardError.__init__(self, *args, **kwargs)

##------------------------------------##
##                                    ##
##------------------------------------##

class SNClass():
    """docstring for SNClass."""

    def __init__(self, verbose = False):
        """

        """

        ## Initialise the class variables
        self._default_data_dir_path = os.path.abspath(os.path.join(__file__, os.pardir, os.pardir) + '/testdata/')

        ## Initialise using class methods
        self.set_data_directory(self._get_data_directory())

    def _get_data_directory(self):
        """
        Looks for the data data directory set as environment variable
        $PYCOCO_DATA_DIR. if not found, returns

        returns: Absolute path in environment variable $PYCOCO_DATA_DIR, or
                 default datalocation: '../testdata/'.
        """

        return os.environ.get('PYCOCO_DATA_DIR', self._default_data_dir_path)

    def set_data_directory(self, data_dir_path = '', verbose = False):
        """
        Enables the data directory to be changed by the user.

        """
        try:
            if os.path.isdir(os.path.abspath(data_dir_path)):
                self.data_directory = os.path.abspath(data_dir_path)
                pass
            else:
                warnings.warn(os.path.abspath(data_dir_path) +
                " is not a valid directory. Restoring default path: " +
                self._default_data_dir_path, UserWarning)
                self.data_directory = self._default_data_dir_path

                if not os.path.isdir(self.data_directory):
                    if verbose: print(os.path.isdir(self.data_directory))
                    raise PathError("The default data directory '" + self.data_directory
                     + "' doesn't exist. Or isn't a directory. Or can't be located.")
                else:
                    pass
        except:
            if verbose: print("foo")
            raise PathError("The default data directory '" + self._default_data_dir_path
             + "' doesn't exist. Or isn't a directory. Or can't be located. Have"
             + " you messed with _default_data_dir_path?")
            pass






##----------------------------------------------------------------------------##
##  /CODE                                                                     ##
##----------------------------------------------------------------------------##
