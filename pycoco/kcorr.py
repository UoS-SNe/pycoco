"""

author: Rob Firth, Southampton
date: 01-2017

'offset' Based on Blanton et al. 2007: 2007AJ....133..734B
The UBRVI filters are those of Bessell (1990). The ugriz filters are those
determined by M. Doi, D. Eisenstein, and J. Gunn and are available on the
SDSS DR4 Web site ( http://www.sdss.org/dr4/ ). The JHKs filters are those from
Cohen et al. (2003).
"""

from __future__ import print_function

import os
import sys

from numpy import log10
from scipy.integrate import simps
from astropy import units as u

from .classes import *
from .functions import *
from .defaults import *

__all__ = ["offset",
            # "convert_AB_to_Vega",
            # "convert_Vega_to_AB",
            "calc_AB_zp",
            "calc_vega_zp",
            "load_dark_sky_spectrum",
            "calc_spectrum_filter_flux",
            "load_atmosphere"]

## offset is calculated as m_AB - m_vega
offset = {
    "U" : 0.79,
    "B" : 0.09,
    "V" : 0.02,
    "R" : 0.21,
    "I" : 0.45,
    "u" : 0.91,
    "g" : 0.08,
    "r" : 0.16,
    "i" : 0.37,
    "z" : 0.54,
    "J" : 0.91,
    "H" : 1.39,
    "Ks" : 1.85 ,
    "0:1u" : 1.25 ,
    "0:1g" : 0.01 ,
    "0:1r" : 0.04 ,
    "0:1i" : 0.27 ,
    "0:1z" : 0.46,
}

# def convert_Vega_to_AB(phot_table, filters = False):
#     """
#     Parameters
#     ----------
#     Returns
#     -------
#     """
#
#     return phot_table
#

# def convert_AB_to_Vega():
#     """
#     Parameters
#     ----------
#     Returns
#     -------
#     """
#
#     return phot_table


def load_vega(path = os.path.join(_default_kcorr_data_path, "alpha_lyr_stis_002.dat")):
    """
    returns spectrum of Vega as a SpectrumClass instance
    """
    vega = SpectrumClass()
    vega.load(path)

    return vega


def load_AB(path = os.path.join(_default_kcorr_data_path, "AB_pseudospectrum.dat")):
    """
    returns 'spectrum' as a SpectrumClass instance
    """
    vega = SpectrumClass()
    vega.load(path)

    return vega


def load_atmosphere(path = os.path.join(_default_lsst_throughputs_path, "baseline/atmos_std.dat")):
    """
    reads in atmosphere from LSST_THROUGHPUTS, default is at airmass 1.2
    """

    atmos = BaseFilterClass()
    atmos.load(path, wavelength_u = u.nm, fmt = "ascii.commented_header")

    return atmos


def calc_filter_area(filter_name):
    filter_object = load_filter(_default_filter_dir_path + filter_name + ".dat")
    filter_area = simps(filter_object.throughput, filter_object.wavelength)
    return filter_area


def calc_spectrum_filter_flux(filter_name, SpecClass):
    filter_object = load_filter(_default_filter_dir_path + filter_name + ".dat")
    filter_object.resample_response(new_wavelength = SpecClass.wavelength)
    filter_area = simps(filter_object.throughput, filter_object.wavelength)

    transmitted_spec = filter_object.throughput * SpecClass.flux

    integrated_flux = simps(transmitted_spec, SpecClass.wavelength)

    return  integrated_flux/filter_area


def calc_AB_flux(filter_name):

    AB = load_AB()

    filter_object = load_filter("/Users/berto/Code/CoCo/data/filters/" + filter_name + ".dat")
    filter_object.resample_response(new_wavelength = AB.wavelength)

    transmitted_spec = filter_object.throughput * AB.flux

    integrated_flux = simps(transmitted_spec, AB.wavelength)

    return integrated_flux


def calc_AB_zp(filter_name):
    """

    """

    integrated_flux = calc_AB_flux(filter_name)
    area_corr_integrated_flux = integrated_flux / calc_filter_area(filter_name)

    return -2.5 * log10(area_corr_integrated_flux)


def calc_vega_flux(filter_name, filter_object = False,):
    """

    """

    vega = load_vega()

    if not filter_object:
        filter_object = load_filter("/Users/berto/Code/CoCo/data/filters/" + filter_name + ".dat")
    # else if hasattr(filter_object, "wavelength"):

    filter_object.resample_response(new_wavelength = vega.wavelength)

    transmitted_spec = filter_object.throughput * vega.flux

    integrated_flux = simps(transmitted_spec, vega.wavelength)

    return integrated_flux


def calc_vega_zp(filter_name, filter_object = False, vega_Vmag = 0.03):
    """

    """

    if not filter_object:
        filter_object = load_filter("/Users/berto/Code/CoCo/data/filters/" + filter_name + ".dat")

    integrated_flux = calc_vega_flux(filter_name)
    area_corr_integrated_flux = integrated_flux / calc_filter_area(filter_name)

    # return -2.5 * log10(integrated_V_flux) - vega_Vmag
    return -2.5 * log10(area_corr_integrated_flux)


# def calc_vega_mag(filter_name):
#     """
#
#     """
#     zp = calc_vega_zp(filter_name)
#     flux = calc_vega_flux(filter_name)
#
#     mag = -2.5 * log10(flux) - zp
#     return mag


def load_dark_sky_spectrum():
    """
    requires https://github.com/lsst/throughputs/ and environment vars LSST_THROUGHPUTS
    and LSST_THROUGHPUTS_BASELINE.

    e.g.:
    setenv LSST_THROUGHPUTS ${HOME}/projects/LSST/throughputs
    setenv LSST_THROUGHPUTS_BASELINE ${LSST_THROUGHPUTS}/baseline
    """
    dark_sky_path = os.path.join(os.environ["LSST_THROUGHPUTS_BASELINE"],"darksky.dat")
    darksky = SpectrumClass()
    darksky.load(dark_sky_path, wavelength_u = u.nm, fmt = "ascii.commented_header",
                 wmin = 3500*u.angstrom, wmax = 11000*u.angstrom)

    darksky.success = True

    return darksky




def convert_f_nu_to_f_lambda():
    pass
