from __future__ import print_function ## Force python3-like printing

import numpy as np
import rfutils as rfu
import rfcolours as rfc

from matplotlib import pyplot as plt

import os
from astropy.io import ascii
from astropy.table import Table
from astropy.constants import c

import pycoco as pcc


def load_bianco_phot_table():
    fname = '/Users/berto/data/CoreCollapse/phot/apjs496616t6_mrt.txt'
    return ascii.read(fname, data_start = 18)

def load_bianco_extinction():
    fname = "/Users/berto/data/CoreCollapse/phot/apjs496616t2_mrt.txt"
    data = ascii.read(fname, format = "commented_header", delimiter = "\t")
    return data

def get_EBV(snname):
    data = load_bianco_extinction()
    wsn = np.where(data["snname"] == snname)
    return data["EBV"][wsn]

def make_bianco_AB_phot():
    data = load_bianco_table()

    for sn_id in np.unique(data["SN"]):
        snname = "SN" + sn_id

        print(snname)

if __name__ == "__main__":
    # print("Foo")
    print(get_EBV("2009jf"))
    # make_bianco_AB_phot()

else:
    pass
