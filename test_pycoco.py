"""
Test suite for pycoco
"""

import os
import unittest

# from pycoco import *
import pycoco as pcc

##------------------------------------##
##  TESTING                           ##
##------------------------------------##

class TestClass(unittest.TestCase):
    """
    Class for testing pycoco
    """

    def test_SNClass_get_data_dir_returns_default(self):
        x = pcc.SNClass()
        self.assertEqual(x._get_data_directory(), pcc._default_data_dir_path)

    def test_load_all_phot_returns_PathError_for_None(self):
        self.assertRaises(pcc.PathError, pcc.load_all_phot, None)

    def test_find_phot_finds_SN2005bf_B(self):
        directory_path_to_search = "/Users/berto/Code/verbose-enigma/testdata/"
        phot_path = pcc.find_phot(directory_path_to_search, snname = "SN2005bf", verbose = False)[0]
        phot_filename = phot_path.split('/')[-1]
        self.assertEqual(phot_filename, u'SN2005bf_B.dat')

    def test_find_phot_finds_no_SN2011fe_data(self):
        directory_path_to_search = "/Users/berto/Code/verbose-enigma/testdata/"
        self.assertEqual(len(pcc.find_phot(directory_path_to_search, snname = "SN2011fe", verbose = False)),
                         0)

    def test_find_phot_throws_path_error_for_None(self):
        self.assertRaises(pcc.PathError, pcc.find_phot, None)

    def test_find_phot_returns_False_for_zoidberg(self):
        self.assertEqual(pcc.find_phot('Zoidberg!'), False)

    def test_check_dir_path_finds_pycoco_dir(self):
        self.assertEqual(pcc.check_dir_path(pcc._default_data_dir_path), True)

    def test_check_dir_path_raises_PathError_for_None(self):
        self.assertRaises(pcc.PathError, pcc.check_dir_path, None)

    def test_check_dir_path_returns_False_for_file(self):
        self.assertEqual(pcc.check_dir_path(__file__), False)

    def test_check_file_path_finds_SN2005bf_B(self):
        self.assertEqual(pcc.check_file_path(os.path.join(pcc._default_data_dir_path, 'SN2005bf_B.dat')), True)

    def test_check_file_path_raises_PathError_for_None(self):
        self.assertRaises(pcc.PathError, pcc.check_file_path, None)

    def test_check_file_path_returns_False_for_dir(self):
        self.assertEqual(pcc.check_file_path(pcc._default_data_dir_path), False)


if __name__ is '__main__':

    test = False
    test = True

    if test:

        print("Running test suite:\n")

        suite = unittest.TestLoader().loadTestsFromTestCase(TestClass)
        unittest.TextTestRunner(verbosity=2).run(suite)

else:

    pass
