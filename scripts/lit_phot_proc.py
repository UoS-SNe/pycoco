from __future__ import print_function ## Force python3-like printing

import numpy as np
import rfutils as rfu
import rfcolours as rfc

from matplotlib import pyplot as plt

import os
from astropy.io import ascii
from astropy.table import Table, Column
from astropy.constants import c
from collections import OrderedDict

import pycoco as pcc


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

def SN2009jf():
    get_EBV("2009jf")
    pass

if __name__ == "__main__":
    # print("Foo")
    # print(get_EBV("2009jf"))
    # make_bianco_AB_phot()

    ## SN2009jf
    data = load_bianco_phot_table()
    snjf = data[np.where(data['SN'] == '2009jf')]

else:
    pass
