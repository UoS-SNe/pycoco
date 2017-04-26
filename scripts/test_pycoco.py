"""
Test suite for pycoco
"""

import os
import unittest

# from pycoco import *
import pycoco as pcc
import astropy.units as u

try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3

reload(pcc)
##------------------------------------##
##  TESTING                           ##
##------------------------------------##

class TestClass(unittest.TestCase):
    """
    Class for testing pycoco
    """

    def test_environment_variables_exist(self):
        self.assertTrue("COCO_ROOT_DIR" in os.environ)

    # def test_SNClass_get_data_dir_returns_default(self):
    #     x = pcc.SNClass()
    #     self.assertEqual(x._get_data_directory(), pcc._default_data_dir_path)

    def test_load_all_phot_returns_PathError_for_None(self):
        self.assertRaises(pcc.PathError, pcc.load_all_phot, None)

    def test_find_filter_phot_finds_5_SN2005bf(self):
        directory_path_to_search = os.path.abspath(os.path.join(pcc._default_data_dir_path, "lc"))
        phot_path = pcc.find_filter_phot(directory_path_to_search, snname = "SN2005bf", verbose = False, prefix = "")
        # phot_filename = phot_path.split('/')[-1]
        # self.assertEqual(phot_filename, 'SN2005bf.dat')
        self.assertEqual(len(phot_path), 5)

    def test_find_formatted_phot_finds_SN2005bf(self):
        directory_path_to_search = os.path.abspath(os.path.join(pcc._default_data_dir_path, "lc"))
        phot_path = pcc.find_formatted_phot(directory_path_to_search, snname = "SN2005bf", verbose = False, prefix = "")[0]
        phot_filename = phot_path.split('/')[-1]
        self.assertEqual(phot_filename, 'SN2005bf.dat')

    def test_find_filter_phot_finds_no_SN2011fe_data(self):
        directory_path_to_search = os.path.abspath(os.path.join(pcc._default_data_dir_path, "lc"))
        self.assertEqual(len(pcc.find_filter_phot(directory_path_to_search, snname = "SN2011fe", verbose = False)), 0)

    def test_find_formatted_phot_finds_no_SN2011fe_data(self):
        directory_path_to_search = os.path.abspath(os.path.join(pcc._default_data_dir_path, "lc"))
        self.assertEqual(len(pcc.find_formatted_phot(directory_path_to_search, snname = "SN2011fe", verbose = False)), 0)

    def test_find_formatted_phot_throws_path_error_for_None(self):
        self.assertRaises(pcc.PathError, pcc.find_formatted_phot, None)

    def test_find_formatted_phot_throws_path_error_for_None(self):
        self.assertRaises(pcc.PathError, pcc.find_formatted_phot, None)

    def test_find_phot_returns_PathError_for_zoidberg(self):
        # self.assertEqual(pcc.find_filter_phot('Zoidberg!'), False)
        self.assertRaises(pcc.PathError, pcc.find_filter_phot, "Zoidberg!")

    def test_check_dir_path_finds_pycoco_dir(self):
        self.assertEqual(pcc.check_dir_path(pcc._default_data_dir_path), True)

    def test_check_dir_path_raises_PathError_for_None(self):
        self.assertRaises(pcc.PathError, pcc.check_dir_path, None)

    def test_check_dir_path_raises_PathError_for_file(self):
        # self.assertEqual(pcc.check_dir_path(__file__), False)
        self.assertRaises(pcc.PathError, pcc.check_dir_path, __file__)

    def test_check_file_path_finds_SN2005bf(self):
        self.assertEqual(pcc.check_file_path(os.path.join(pcc._default_data_dir_path, 'lc/SN2005bf.dat')), True)

    def test_check_file_path_raises_PathError_for_None(self):
        self.assertRaises(pcc.PathError, pcc.check_file_path, None)

    def test_check_file_path_raises_PathError_for_dir(self):
        # self.assertEqual(pcc.check_file_path(pcc._default_data_dir_path), False)
        self.assertRaises(pcc.PathError, pcc.check_file_path, pcc._default_data_dir_path)

    ## CLASS TESTS ##

    # PhotometryClass

    def test_PhotometryClass_get_data_dir_returns_default(self):
        x = pcc.PhotometryClass()
        self.assertEqual(os.path.abspath(os.path.join(x._get_data_directory(), os.pardir)), os.path.abspath(pcc._default_data_dir_path))

    # def LOAD FROM FILE

    # BaseSpectrumClass

    def test_BaseSpectrumClass_get_list_dir_returns_default(self):
        x = pcc.BaseSpectrumClass()
        listpath = os.path.abspath(os.path.join(pcc._default_coco_dir_path, "lists"))
        self.assertEqual(os.path.abspath(x.list_directory), listpath)

    # SpectrumClass

    def test_SpectrumClass_get_data_dir_returns_default(self):
        x = pcc.SpectrumClass()
        self.assertEqual(os.path.abspath(os.path.join(x._get_data_directory(), os.pardir)), os.path.abspath(pcc._default_data_dir_path))

    def test_SpectrumClass_get_data_dir_returns_default_spec_dir(self):
        x = pcc.SpectrumClass()
        self.assertEqual(os.path.abspath(x._get_data_directory()), os.path.abspath(os.path.join(pcc._default_data_dir_path, "spec/")))

    # specfitClass

    # def test_specfitClass

    ## FilterClass tests

    def test_filter_name_parsed_OK(self):
        filter_filename = "BessellB.dat"
        path_to_filter = os.path.join(os.path.abspath(pcc._default_filter_dir_path), filter_filename)
        B = pcc.load_filter(path_to_filter)
        self.assertEqual(B.filter_name, "BessellB")

    def test_filter_default_colour(self):
        filter_filename = "BessellB.dat"
        path_to_filter = os.path.join(os.path.abspath(pcc._default_filter_dir_path), filter_filename)
        B = pcc.load_filter(path_to_filter)
        self.assertEqual(B._plot_colour, "#0000ff")

    def test_filter_default_wavelength_units(self):
        filter_filename = "BessellV.dat"
        path_to_filter = os.path.join(os.path.abspath(pcc._default_filter_dir_path), filter_filename)
        V = pcc.load_filter(path_to_filter)
        self.assertEqual(V._wavelength_units, u.Angstrom)

    def test_filter_default_frequency_units(self):
        filter_filename = "BessellV.dat"
        path_to_filter = os.path.join(os.path.abspath(pcc._default_filter_dir_path), filter_filename)
        V = pcc.load_filter(path_to_filter)
        self.assertEqual(V._frequency_units, u.Hertz)

    def test_filter_edge_calc(self):
        filter_filename = "BessellB.dat"
        path_to_filter = os.path.join(os.path.abspath(pcc._default_filter_dir_path), filter_filename)
        B = pcc.load_filter(path_to_filter)
        self.assertEqual(round(float(B._upper_edge), 2), 5203.28)
        self.assertEqual(round(float(B._lower_edge), 2), 3784.99)

    def test_filter_LSST_filters_load(self):
        LSST_filter_list = ["LSST_u",
                            "LSST_g",
                            "LSST_r",
                            "LSST_i",
                            "LSST_z",
                            "LSST_y"]
        success = []

        for LSST_filter_name in LSST_filter_list:
            path_to_filter = os.path.join(os.path.abspath(pcc._default_filter_dir_path), LSST_filter_name + ".dat")
            try:
                F = pcc.load_filter(path_to_filter)
                success.append(1)
            except:
                success.append(0)

        self.assertEqual(sum(success), 6)

    ##


if __name__ is '__main__':

    test = False
    test = True

    if test:

        print("Running test suite:\n")

        suite = unittest.TestLoader().loadTestsFromTestCase(TestClass)
        unittest.TextTestRunner(verbosity=2).run(suite)

else:

    pass
