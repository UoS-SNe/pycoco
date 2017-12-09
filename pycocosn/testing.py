"""
Test suite for pycoco
"""
try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3

import os
import unittest

import pycoco as pcc
import astropy.units as u
import matplotlib as mpl
import numpy as np


#  #------------------------------------#  #
#  #  TESTING                           #  #
#  #------------------------------------#  #

class TestClass(unittest.TestCase):
    """
    Class for testing pycocosn
    """

    ## Environment Variable Tests
    def test_COCO_ROOT_DIR_environment_variable_exists(self):
        self.assertTrue("COCO_ROOT_DIR" in os.environ)

    # def test_PYCOCO_FILTER_DIR_environment_variable_exists(self):
    #     self.assertTrue("PYCOCO_FILTER_DIR" in os.environ)
    #
    # def test_PYCOCO_DATA_DIR_environment_variable_exists(self):
    #     self.assertTrue("PYCOCO_DATA_DIR" in os.environ)

    def test_SFD_DIR_environment_variable_exists(self):
        self.assertTrue("SFD_DIR" in os.environ)

    def test_LSST_THROUGHPUTS_environment_variable_exists(self):
        self.assertTrue("LSST_THROUGHPUTS" in os.environ)

    # def test_LSST_THROUGHPUTS_BASELINE_environment_variable_exists(self):
    #     self.assertTrue("LSST_THROUGHPUTS_BASELINE" in os.environ)

    ##

    def test_COCO_ROOT_DIR_exists(self):
        self.assertTrue(os.path.isdir(pcc.defaults._default_coco_dir_path))

    def test_PYCOCO_FILTER_DIR_exists(self):
        self.assertTrue(os.path.isdir(pcc.defaults._default_filter_dir_path))

    def test_PYCOCO_DATA_DIR_exists(self):
        self.assertTrue(os.path.isdir(pcc.defaults._default_data_dir_path))

    def test_SFD_DIR_exists(self):
        self.assertTrue(os.path.isdir(pcc.defaults._default_dust_dir))

    def test_LSST_THROUGHPUTS_DIR_exists(self):
        self.assertTrue(os.path.isdir(pcc.defaults._default_lsst_throughputs_path))

    def test_LSST_THROUGHPUTS_BASELINE_DIR_exists(self):
        self.assertTrue(os.path.isdir(os.path.abspath(os.environ["LSST_THROUGHPUTS_BASELINE"])))

    ## Tests

    # def test_SNClass_get_data_dir_returns_default(self):
    #     x = pcc.classes.SNClass()
    #     self.assertEqual(x._get_data_directory(), pcc.defaults._default_data_dir_path)

    def test_load_all_phot_returns_PathError_for_None(self):
        self.assertRaises(pcc.errors.PathError, pcc.functions.load_all_phot, None)

    # def test_find_filter_phot_finds_5_SN2005bf(self):
    #     directory_path_to_search = os.path.abspath(os.path.join(pcc.defaults._default_data_dir_path, "lc"))
    #     phot_path = pcc.find_filter_phot(directory_path_to_search, snname = "SN2005bf", verbose = False, prefix = "")
    #     # phot_filename = phot_path.split('/')[-1]
    #     # self.assertEqual(phot_filename, 'SN2005bf.dat')
    #     self.assertEqual(len(phot_path), 5)

    def test_find_formatted_phot_finds_SN2005bf(self):
        directory_path_to_search = os.path.abspath(os.path.join(pcc.defaults._default_data_dir_path, "lc"))
        phot_path = \
        pcc.functions.find_formatted_phot(directory_path_to_search, snname="SN2005bf", verbose=False, prefix="")[0]
        phot_filename = phot_path.split('/')[-1]
        self.assertEqual(phot_filename, 'SN2005bf.dat')

    def test_find_filter_phot_finds_no_SN2011fe_data(self):
        directory_path_to_search = os.path.abspath(os.path.join(pcc.defaults._default_data_dir_path, "lc"))
        self.assertEqual(
            len(pcc.functions.find_filter_phot(directory_path_to_search, snname="SN2011fe", verbose=False)), 0)

    def test_find_formatted_phot_finds_no_SN2011fe_data(self):
        directory_path_to_search = os.path.abspath(os.path.join(pcc.defaults._default_data_dir_path, "lc"))
        self.assertEqual(
            len(pcc.functions.find_formatted_phot(directory_path_to_search, snname="SN2011fe", verbose=False)), 0)

    def test_find_formatted_phot_throws_path_error_for_None(self):
        self.assertRaises(pcc.errors.PathError, pcc.functions.find_formatted_phot, None)

    def test_find_formatted_phot_throws_path_error_for_None(self):
        self.assertRaises(pcc.errors.PathError, pcc.functions.find_formatted_phot, None)

    def test_find_phot_returns_PathError_for_zoidberg(self):
        # self.assertEqual(pcc.find_filter_phot('Zoidberg!'), False)
        self.assertRaises(pcc.errors.PathError, pcc.functions.find_filter_phot, "Zoidberg!")

    def test_check_dir_path_finds_pycoco_dir(self):
        self.assertEqual(pcc.utils.check_dir_path(pcc.defaults._default_data_dir_path), True)

    def test_check_dir_path_raises_PathError_for_None(self):
        self.assertRaises(pcc.errors.PathError, pcc.utils.check_dir_path, None)

    def test_check_dir_path_raises_PathError_for_file(self):
        # self.assertEqual(pcc.utils.check_dir_path(__file__), False)
        self.assertRaises(pcc.errors.PathError, pcc.utils.check_dir_path, __file__)

    def test_check_file_path_finds_SN2005bf(self):
        self.assertEqual(
            pcc.utils.check_file_path(os.path.join(pcc.defaults._default_data_dir_path, 'lc/SN2005bf.dat')), True)

    def test_check_file_path_raises_PathError_for_None(self):
        self.assertRaises(pcc.errors.PathError, pcc.utils.check_file_path, None)

    def test_check_file_path_raises_PathError_for_dir(self):
        self.assertRaises(pcc.errors.PathError, pcc.utils.check_file_path, pcc.defaults._default_data_dir_path)

    def test_find_specphase_spec_raises_PathError_for_None(self):
        self.assertRaises(pcc.errors.PathError, pcc.classes.find_specphase_spec, None)

    def test_find_specphase_spec_returns_empty_array_for_zoidberg(self):
        self.assertEqual(pcc.classes.find_specphase_spec("Zoidberg"), [])

    def test_find_specphase_spec_returns_19_for_SN2006aj(self):
        self.assertEqual(len(pcc.classes.find_specphase_spec("SN2006aj")), 19)

    def test_load_info_finds_and_loads_default(self):
        i = pcc.functions.load_info()
        self.assertTrue(i.table.meta["success"])

    def test_load_info_finds_default(self):
        i = pcc.functions.load_info()
        self.assertEqual(len(i.table), 28)

    ## CLASS TESTS ##
    # classes.PhotometryClass

    def test_PhotometryClass_get_data_dir_returns_default(self):
        x = pcc.classes.PhotometryClass()
        self.assertEqual(os.path.abspath(os.path.join(x._get_data_directory(), os.pardir)),
                         os.path.abspath(pcc.defaults._default_data_dir_path))

    def test_PhotometryClass_get_1993J(self):
        x = pcc.classes.PhotometryClass()
        x.load(os.path.join(pcc.defaults._default_data_dir_path, "lc/SN1993J.dat"))

        self.assertTrue(hasattr(x, "phot"))

    def test_PhotometryClass_1993J_phot_size(self):
        x = pcc.classes.PhotometryClass()
        x.load(os.path.join(pcc.defaults._default_data_dir_path, "lc/SN1993J.dat"))

        self.assertEqual(len(x.phot), 607)

    def test_PhotometryClass_get_and_plot_1993J(self):
        x = pcc.classes.PhotometryClass()
        x.load(os.path.join(pcc.defaults._default_data_dir_path, "lc/SN1993J.dat"))
        fig = x.plot(return_figure=True)
        self.assertIsInstance(fig, mpl.figure.Figure)

    # def LOAD FROM FILE

    # classes.SpectrumClass

    def test_BaseSpectrumClass_get_list_dir_returns_default(self):
        x = pcc.classes.SpectrumClass()
        listpath = os.path.abspath(os.path.join(pcc.defaults._default_coco_dir_path, "lists"))
        self.assertEqual(os.path.abspath(x.list_directory), listpath)

    # SpectrumClass

    def test_SpectrumClass_get_data_dir_returns_default(self):
        x = pcc.classes.SpectrumClass()
        self.assertEqual(os.path.abspath(os.path.join(x._get_data_directory(), os.pardir)),
                         os.path.abspath(pcc.defaults._default_data_dir_path))

    # def test_SpectrumClass_get_data_dir_returns_default_spec_dir(self):
    #     x = pcc.classes.SpectrumClass()
    #     self.assertEqual(os.path.abspath(x._get_data_directory()),
    #                      os.path.abspath(os.path.join(pcc.defaults._default_data_dir_path, "spec/")))

    # def test_SpectrumClass_get_specphot_works_more_than_one_overlapping_filter(self):


    def test_SpectrumClass_get_specphot_works_with_one_requested_and_one_overlapping_filter(self):
        snname = "SN1993J"

        listfile = os.path.join(pcc.defaults._default_list_dir_path, snname + ".list")

        sn = pcc.classes.SNClass(snname)
        sn.load_phot(path=os.path.join(pcc.defaults._default_data_dir_path, "lc/" + snname + ".dat"))
        sn.load_list(listfile)
        sn.load_spec()
        sn.check_overlaps()
        S = sn.spec["1993J_-11.0.txt"]
        S.get_specphot(sn.phot.data_filters["BessellV"], verbose=False)
        self.assertAlmostEqual(1.373e-17, S.specphot["flux"][np.where(S.specphot["filter"] == "BessellV")[0][0]], 3)

    def test_SpectrumClass_get_specphot_works_with_one_requested_and_many_overlapping_filters(self):
        snname = "SN1993J"

        listfile = os.path.join(pcc.defaults._default_list_dir_path, snname + ".list")

        sn = pcc.classes.SNClass(snname)
        sn.load_phot(path=os.path.join(pcc.defaults._default_data_dir_path, "lc/" + snname + ".dat"))
        sn.load_list(listfile)
        sn.load_spec()
        sn.check_overlaps()
        S = sn.spec["1993J_-3.0.txt"]

        S.get_specphot(sn.phot.data_filters["BessellV"], verbose=False)
        self.assertAlmostEqual(1.907e-17, S.specphot["flux"][np.where(S.specphot["filter"] == "BessellV")[0][0]], 3)

    def test_SpectrumClass_get_specphot_works_with_many_requested_and_many_overlapping_filters(self):
        snname = "SN2007uy"

        listfile = os.path.join(pcc.defaults._default_list_dir_path, snname + ".list")

        sn = pcc.classes.SNClass(snname)
        sn.load_phot(path=os.path.join(pcc.defaults._default_data_dir_path, "lc/" + snname + ".dat"))
        sn.load_list(listfile)
        sn.load_spec()
        sn.check_overlaps()
        S = sn.spec["2007uy_-5.06.txt"]

        S.get_specphot(sn.phot.data_filters, verbose=False)

        expected = [1.31306328e-17, 1.48349665e-17, 1.18406322e-17, 7.79000154e-18]
        for i in zip(expected, np.array(S.specphot["flux"])):
            self.assertAlmostEqual(i[0], i[1])
            # self.assertAlmostEqual(expected, np.array(S.specphot["flux"]), 3)

    # specfitClass

    # def test_specfitClass

    # def test_compare_spec()

    ## FilterClass tests

    def test_filter_name_parsed_OK(self):
        filter_filename = "BessellB.dat"
        path_to_filter = os.path.join(os.path.abspath(pcc.defaults._default_filter_dir_path), filter_filename)
        B = pcc.functions.load_filter(path_to_filter)
        self.assertEqual(B.filter_name, "BessellB")

    def test_filter_default_colour(self):
        filter_filename = "BessellB.dat"
        path_to_filter = os.path.join(os.path.abspath(pcc.defaults._default_filter_dir_path), filter_filename)
        B = pcc.functions.load_filter(path_to_filter)
        self.assertEqual(B._plot_colour, "#0000ff")

    def test_filter_default_wavelength_units(self):
        filter_filename = "BessellV.dat"
        path_to_filter = os.path.join(os.path.abspath(pcc.defaults._default_filter_dir_path), filter_filename)
        V = pcc.functions.load_filter(path_to_filter)
        self.assertEqual(V._wavelength_units, u.Angstrom)

    def test_filter_default_frequency_units(self):
        filter_filename = "BessellV.dat"
        path_to_filter = os.path.join(os.path.abspath(pcc.defaults._default_filter_dir_path), filter_filename)
        V = pcc.functions.load_filter(path_to_filter)
        self.assertEqual(V._frequency_units, u.Hertz)

    def test_filter_edge_calc(self):
        filter_filename = "BessellB.dat"
        path_to_filter = os.path.join(os.path.abspath(pcc.defaults._default_filter_dir_path), filter_filename)
        B = pcc.functions.load_filter(path_to_filter)
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
            print("loading", LSST_filter_name)
            path_to_filter = os.path.join(os.path.abspath(pcc.defaults._default_filter_dir_path),
                                          LSST_filter_name + ".dat")
            print("from", path_to_filter)
            print("file exists? ", os.path.isfile(path_to_filter))
            try:
                F = pcc.functions.load_filter(path_to_filter)
                success.append(1)
            except:
                success.append(0)

        self.assertEqual(sum(success), 6)

    ## kcorr tests

    def test_kcorr_load_AB_pseudospectrum(self):
        ABspec = pcc.kcorr.load_AB()
        self.assertEqual(ABspec.success, True)

    def test_kcorr_load_vega_pseudospectrum(self):
        Vegaspec = pcc.kcorr.load_vega()
        self.assertEqual(Vegaspec.success, True)

    def test_dark_sky_spectrum_exists(self):
        dark_sky_path = os.path.join(os.environ["LSST_THROUGHPUTS_BASELINE"], "darksky.dat")
        print(dark_sky_path)
        print("file exists?", os.path.isfile(dark_sky_path))
        self.assertTrue(os.path.isfile(dark_sky_path))

    def test_kcorr_load_dark_sky_spectrum(self):
        dark_sky_path = os.path.join(os.environ["LSST_THROUGHPUTS_BASELINE"], "darksky.dat")

        darksky = pcc.classes.SpectrumClass()
        darksky.load(dark_sky_path, wavelength_u=u.nm, fmt="ascii.commented_header", abspath=True)
        self.assertEqual(darksky.success, True)

    def test_kcorr_calc_m_darksky(self):
        m_darkskyV = pcc.kcorr.calc_m_darksky("BessellV", abspath=True)
        self.assertAlmostEqual(21.717677839340244, m_darkskyV, 2)

    ###

    ## test CoCo stuff

    # def test_coco_calls_test_LCfit(self):

    def test_coco_lists_phase_is_monotonic(self):
        checklist = pcc.utils.check_all_lists(pcc.defaults._default_list_dir_path)
        self.assertTrue(all(checklist))


def runtests(test = True):
    """

    :param test:
    :return:
    """
    if test:
        print("Running test suite:\n")

        suite = unittest.TestLoader().loadTestsFromTestCase(TestClass)
        unittest.TextTestRunner(verbosity=2).run(suite)
