from __future__ import print_function ## Force python3-like printing

import numpy as np
import rfutils as rfu
import rfcolours as rfc

from matplotlib import pyplot as plt

import os
from astropy.io import ascii
from astropy.time import Time
from astropy.table import Table, Column
from astropy.constants import c
from collections import OrderedDict

import pycoco as pcc
import pycoco.kcorr as kcorr

class LitLightCurveClass(pcc.BaseLightCurveClass):
    pass


def JD_to_MJD():
    """
    """

    return mjd

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

    ## Vega Landolt
    sn2009jf_UBVRI = ascii.read("/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/astropy_tables/sn2009jf_UBVRI.dat")
    sn2009jf_UBVRI["mjd"] = Time(sn2009jf_UBVRI["JD"], format = "jd").mjd

    filter_names = ["U", "B", "V", "R", "I"]
    for filter_name in filter_names:
        print(filter_name + "mag")

        if filter_name in translation_dict:
            full_filter_name = translation_dict[filter_name]

        vega_mag = sn2009jf_UBVRI[filter_name + "mag"]

        magmask = np.logical_not(sn2009jf_UBVRI[filter_name + "mag"].mask)
        AB_offset = kcorr.calc_offset_AB_minus_Vega(full_filter_name)## This is AB - Vega

        AB_mag = vega_mag + AB_offset
        phot_table = Table

        # print(vega_mag[magmask], AB_mag[magmask])
        # print(AB_mag + AB_zp)
        # f_nu = np.power(10., (AB_mag + AB_zp)/-2.5)
        # print(f_nu)

    sn2009jf_ugriz = ascii.read("/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/astropy_tables/sn2009jf_ugriz.dat")
    sn2009jf_ugriz["mjd"] = Time(sn2009jf_ugriz["JD"], format = "jd").mjd
    sn2009jf_JHK = ascii.read("/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/astropy_tables/sn2009jf_JHK.dat")
    sn2009jf_JHK["mjd"] = Time(sn2009jf_JHK["JD"], format = "jd").mjd
    sn2009jf_swift = ascii.read("/Users/berto/data/CoreCollapse/phot/rf/SN2009jf/astropy_tables/sn2009jf_swift.dat")
    sn2009jf_swift["mjd"] = Time(sn2009jf_swift["JD"], format = "jd").mjd



if __name__ == "__main__":
    # print("Foo")
    # print(get_EBV("2009jf"))
    # make_bianco_AB_phot()

    ## SN2009jf
    data = load_bianco_phot_table()
    snjf = data[np.where(data['SN'] == '2009jf')]

else:
    pass
