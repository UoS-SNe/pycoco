from __future__ import print_function  ## Force python3-like printing

import os
from collections import OrderedDict

import numpy as np
from astropy.constants import c
from astropy.io import ascii
from astropy.table import Table, Column, vstack, MaskedColumn
from astropy.time import Time

# from .classes import *
# from .functions import *
from . import classes
from . import functions

__all__ = ["LitLightCurveClass",
           "JD_to_MJD",
           "write_dict",
           "load_bianco_phot_table",
           "load_bianco_extinction",
           "get_EBV",
           "make_bianco_AB_phot",
           "translate_filter_names",
           "make_new_CfA_phot_table",
           "load_CfA_phot_table",
           "SN2009jf_read_lit",
           "SN2009jf_read_ap"]


class LitLightCurveClass(classes.BaseLightCurveClass):
    pass

translation_dict = OrderedDict()

translation_dict = {"U" : "BessellU",
                    "B" : "BessellB",
                    "V" : "BessellV",
                    "R" : "BessellR",
                    "I" : "BessellI",
                    "u'": "SDSS_u",
                    "g'": "SDSS_g",
                    "r'": "SDSS_r",
                    "i'": "SDSS_i",
                    "z'": "SDSS_z",
                    "u": "SDSS_u",
                    "g": "SDSS_g",
                    "r": "SDSS_r",
                    "i": "SDSS_i",
                    "z": "SDSS_z",
                    "J" : "2MASS_J",
                    "H" : "2MASS_H",
                    "Ks": "2MASS_Ks"}

def JD_to_MJD(JD_column):
    """
    convert astropy column object of JD floats into MJD's
    """
    MJD_col = Column(Time(JD_column, format = "jd").mjd, name = "MJD")
    return MJD_col


def write_dict(out_dir, sn_dict):
    """
    Each dict entry should be an astropy table. Taking this intermediate step in
    order to make it easier to get the formatting correct
    """
    for key in sn_dict:
        # print(key)
        out_path = os.path.join(out_dir, key + ".dat")
        print(out_path)
        sn_dict[key].write(out_path, format = "ascii.fast_commented_header")
    pass


def load_bianco_phot_table(verbose = False):
    """
    from http://iopscience.iop.org/article/10.1088/0067-0049/213/2/19
    Parameters
    ----------
    Returns
    -------
    """

    fname = '/Users/berto/data/CoreCollapse/phot/apjs496616t6_mrt.txt'
    return ascii.read(fname, data_start = 18)


def load_bianco_extinction(verbose = False):
    """
    from http://iopscience.iop.org/article/10.1088/0067-0049/213/2/19
    Parameters
    ----------
    Returns
    -------
    """
    fname = "/Users/berto/data/CoreCollapse/phot/apjs496616t2_mrt.txt"
    data = ascii.read(fname, format = "commented_header", delimiter = "\t")

    zcol = Column(data["cz"]*1000./c.value, name = "z") ## cz is in km s-1

    data.add_column(zcol)
    return data


def get_EBV(snname):
    data = load_bianco_extinction()
    wsn = np.where(data["snname"] == snname)
    return data["EBV"][wsn]


def make_bianco_AB_phot():
    data = load_bianco_phot_table()

    for sn_id in np.unique(data["SN"]):
        snname = "SN" + sn_id

        print(snname)


def translate_filter_names(data, verbose = False):

    translation_dict = OrderedDict()

    translation_dict = {
    "U" : "BessellU",
    "B" : "BessellB",
    "V" : "BessellV",
    "R" : "BessellR",
    "I" : "BessellI",
    "u'": "SDSS_u",
    "g'": "SDSS_g",
    "r'": "SDSS_r",
    "i'": "SDSS_i",
    "z'": "SDSS_z",
    "u": "SDSS_u",
    "g": "SDSS_g",
    "r": "SDSS_r",
    "i": "SDSS_i",
    "z": "SDSS_z",
    "J" : "2MASS_J",
    "H" : "2MASS_H",
    "Ks": "2MASS_Ks"
    }

    translated_arr = np.array([])
    for i, entry in enumerate(data):
        if verbose: print(i, entry)
        # data["NewFilter"][i] = translation_dict[entry["Filter"]]
        translated_arr = np.append(translated_arr, translation_dict[entry["Filter"]])
    newfilter_Column = Column(name = "Filter", data = translated_arr)
    print(newfilter_Column)
    data.remove_column("Filter")
    data["Filter"] = newfilter_Column

    return data


def make_new_CfA_phot_table(verbose = False):
    """
    from http://iopscience.iop.org/article/10.1088/0067-0049/213/2/19

    Parameters
    ----------
    Returns
    -------
    """

    data = load_bianco_phot_table(verbose = verbose)

    data = translate_filter_names(data, verbose = verbose)

    extinction_table = load_bianco_extinction(verbose = verbose)

    EBV_col = data['mag']*0
    EBV_col.name = "EBV"

    for entry in extinction_table:
        if verbose: print(entry)
        w = np.where(data["SN"] == entry["snname"])
        if verbose: print(w)
        EBV_col[w] = entry["EBV"]

    data.add_column(EBV_col)
    return data


def load_CfA_phot_table(verbose = False):
    fname = "/Users/berto/data/CoreCollapse/phot/apjs496616_combined_translated_table_RF.txt"
    data = ascii.read(fname)
    return data


def SN2009jf_read_lit():
    """
    Reads the tables I have collected for SN2009jf
    """

    sn2009jf_dict = OrderedDict()

    ## Get the Phot from Bianco et al.
    data = load_CfA_phot_table()
    sn2009jf_dict['sn2009jf_CfA'] = data[np.where(data['SN'] == '2009jf')]

    ## Get Valenti et al. 2001
    freadme = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/J_MNRAS_416_3138/ReadMe"

    fnamec1 = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/J_MNRAS_416_3138/tablec1.dat"
    fnamec2 = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/J_MNRAS_416_3138/tablec2.dat"
    """
    These next few break due to the 'limits' flags in the table. Instead I have downloaded the tables
    in tab-separated-variable VO tables.
    One caveat to this is the the format and the dashed header row get read in too, even with
    data_start set.
    """
    # fnamec3 = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/J_MNRAS_416_3138/tablec3.dat"
    # fnamec4 = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/J_MNRAS_416_3138/tablec4.dat"
    # fnamec5 = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/J_MNRAS_416_3138/tablec5.dat"
    # fnamec6 = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/J_MNRAS_416_3138/tablec6.dat"

    fnamec3 = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/tablec3.tsv"
    fnamec4 = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/tablec4.tsv"
    fnamec5 = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/tablec5.tsv"
    fnamec6 = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/tablec6.tsv"

    sn2009jf_dict['standards_UBVRI'] = ascii.read(fnamec1, readme = freadme)
    sn2009jf_dict['standards_ugriz'] = ascii.read(fnamec2, readme = freadme)

    sn2009jf_dict['sn2009jf_UBVRI'] = ascii.read(fnamec3)
    sn2009jf_dict['sn2009jf_UBVRI'].remove_rows([0,1])
    sn2009jf_dict['sn2009jf_ugriz'] = ascii.read(fnamec4)
    sn2009jf_dict['sn2009jf_ugriz'].remove_rows([0,1])
    sn2009jf_dict['sn2009jf_JHK'] = ascii.read(fnamec5)
    sn2009jf_dict['sn2009jf_JHK'].remove_rows([0,1])
    sn2009jf_dict['sn2009jf_swift'] = ascii.read(fnamec6)
    sn2009jf_dict['sn2009jf_swift'].remove_rows([0,1])

    return sn2009jf_dict


def SN2009jf_read_ap(format = "ascii", names = ('MJD', 'flux', 'flux_err', 'filter')):
    """
    Reads the output from write_dict(SN2009jf_read_lit) in order to unpack and
    correct the photometry
    """

    snjf = LitLightCurveClass()
    AB_zp = 48.60

    jfdict = OrderedDict()

    ## Vega Landolt
    sn2009jf_UBVRI = ascii.read("/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/astropy_tables/sn2009jf_UBVRI.dat")
    ## MJDs being calculated twice :@:@
    sn2009jf_UBVRI["mjd"] = Time(sn2009jf_UBVRI["JD"], format = "jd").mjd
    ## Unclear if limits are 3sigma or 5sigma - mask
    filter_names = ["U", "B", "V", "R", "I"]

    for filter_name in filter_names:
        print(filter_name + "mag")

        if filter_name in translation_dict:
            full_filter_name = translation_dict[filter_name]

        magmask = np.logical_not(sn2009jf_UBVRI[filter_name + "mag"].mask)
        vega_mag = sn2009jf_UBVRI[filter_name + "mag"]

        if "l_" + filter_name + "mag" in sn2009jf_UBVRI.keys():
            ulim_flag = sn2009jf_UBVRI["l_" + filter_name + "mag"]
            ## See above sigmas comment!!
            magmask = np.logical_and(ulim_flag.mask, magmask)

        AB_offset = kcorr.calc_offset_AB_minus_Vega(full_filter_name)## This is AB - Vega

        mjd_col = JD_to_MJD(sn2009jf_UBVRI["JD"])

        AB_mag = MaskedColumn(vega_mag + AB_offset, name = "mag")
        # e_AB_mag = Column(vega_mag, name = "e_" + AB_mag.name)
        e_AB_mag = MaskedColumn(sn2009jf_UBVRI["e_" + filter_name + "mag"], name = "emag")
        # print(vega_mag[magmask], AB_mag[magmask])
        # print(AB_mag + AB_zp)
        f_nu = MaskedColumn(np.power(10., (AB_mag + AB_zp)/-2.5), name = "f_nu")
        ef_nu = MaskedColumn(f_nu/(1.086/e_AB_mag), name = "ef_nu")
        # print(f_nu)

        filter_column = Column(np.array([full_filter_name for i in mjd_col]), name = "filter")

        phot_table = Table([mjd_col[magmask], AB_mag[magmask], e_AB_mag[magmask], f_nu[magmask], ef_nu[magmask], filter_column[magmask]])

        jfdict[full_filter_name] = phot_table



    sn2009jf_ugriz = ascii.read("/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/astropy_tables/sn2009jf_ugriz.dat")
    sn2009jf_ugriz["mjd"] = Time(sn2009jf_ugriz["JD"], format = "jd").mjd
    filter_names = ["u", "g", "r", "r", "i"]

    for filter_name in filter_names:
        print(filter_name + "mag")

        if filter_name in translation_dict:
            full_filter_name = translation_dict[filter_name]

        ## No ulims for ugriz phot
    #     ulim_flag = sn2009jf_ugriz["l_" + filter_name + "mag"]

        magmask = np.logical_not(sn2009jf_ugriz[filter_name + "mag"].mask)
        mjd_col = JD_to_MJD(sn2009jf_ugriz["JD"])

        AB_mag = MaskedColumn(sn2009jf_ugriz[filter_name + "mag"], name = "mag")
        e_AB_mag = MaskedColumn(sn2009jf_ugriz["e_" + filter_name + "mag"], name = "emag")

        f_nu = MaskedColumn(np.power(10., (AB_mag + AB_zp)/-2.5), name = "f_nu")
        ef_nu = MaskedColumn(f_nu/(1.086/e_AB_mag), name = "ef_nu")

        filter_column = Column(np.array([full_filter_name for i in mjd_col]), name = "filter")
        print(len(mjd_col), len(AB_mag), len(e_AB_mag), len(filter_column))
        phot_table = Table([mjd_col[magmask], AB_mag[magmask], e_AB_mag[magmask], f_nu[magmask], ef_nu[magmask], filter_column[magmask]])
        jfdict[full_filter_name] = phot_table

    # sn2009jf_JHK = ascii.read("/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/astropy_tables/sn2009jf_JHK.dat")
    # sn2009jf_JHK["mjd"] = Time(sn2009jf_JHK["JD"], format = "jd").mjd
    # sn2009jf_swift = ascii.read("/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/astropy_tables/sn2009jf_swift.dat")
    # sn2009jf_swift["mjd"] = Time(sn2009jf_swift["JD"], format = "jd").mjd

    filter_dir = '/Users/berto/Code/CoCo/data/filters/'


    # return jfdict
    for i, phot_filter in enumerate(jfdict.keys()):
        filter_path = os.path.join(filter_dir, phot_filter + ".dat")
        FilterObj = functions.load_filter(filter_path)
        FilterObj.calculate_frequency()
        FilterObj.calculate_effective_frequency()

        ## REMEBER lambda*f_lambda = nu*f_nu
        f_lambda = MaskedColumn((jfdict[phot_filter]["f_nu"]*FilterObj.nu_effective)/FilterObj.lambda_effective, name = "flux")
        ef_lambda = MaskedColumn((jfdict[phot_filter]["ef_nu"]*FilterObj.nu_effective)/FilterObj.lambda_effective, name = "flux_err")

        jfdict[phot_filter].add_column(f_lambda)
        jfdict[phot_filter].add_column(ef_lambda)

        if i == 0:
            full_phot = jfdict[phot_filter]
        else:
            full_phot = vstack([full_phot, jfdict[phot_filter]])

    # flux = np.power(10., (full_phot["mag"] + AB_zp)/-2.5)
    # return full_phot
    # return jfdict
    ## Write Combined Table
    final_phot_table = Table([full_phot["MJD"].filled(), full_phot["flux"].filled(), full_phot["flux_err"].filled(), full_phot["filter"].filled()], masked = False)
    outpath = "/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/final/SN2009jf.dat"
    final_phot_table.write(outpath, format = "ascii.fast_commented_header")
    pass


#
# if __name__ == "__main__":
#     # print("Foo")
#     # print(get_EBV("2009jf"))
#     # make_bianco_AB_phot()
#
#     ## SN2009jf
#     data = load_bianco_phot_table()
#     snjf = data[np.where(data['SN'] == '2009jf')]
#
# else:
#     pass
