'''
This is the module for the PyCoCo python tools

author: Rob Firth, Southampton

date: 06-12-2016
'''

from __future__ import print_function ## Force python3-like printing

##---------------------------------------------------------------------------##
##                                   TOOLS                                   ##
##---------------------------------------------------------------------------##

##------------------------------------##
##  TEST CODE                         ##
##------------------------------------##

class TestClass():
    '''
    Quick test class.

    Contains a test class variable and test class method that prints the
    variable.

    RF
    '''

    def __init__(self):
        self.test_string = 'Hello, World!'

    def print_test_string(self):
        print(self.test_string)

def test_function(*args, **kwargs):
    '''
    Quick test function.

    Prints supplied **args and **kwargs

    RF
    '''

    if args is not None:
        for i, arg in enumerate(args):
            print('an arg passed via *args: ', repr(arg))
    if kwargs is not None:
        for key, value in kwargs.iteritems():
            print('a **kwarg: ', repr(key), ' == ' , repr(value))
    pass

##------------------------------------##
##  CODE                              ##
##------------------------------------##

class SNClass(object):
    """docstring for sn_obj."""
    def __init__(self, arg):
        super(sn_obj, self).__init__()
        self.arg = arg
