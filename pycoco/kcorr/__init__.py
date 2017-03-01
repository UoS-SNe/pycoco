"""
Based on Blanton et al. 2007: 2007AJ....133..734B
The UBRVI filters are those of Bessell (1990). The ugriz filters are those
determined by M. Doi, D. Eisenstein, and J. Gunn and are available on the
SDSS DR4 Web site ( http://www.sdss.org/dr4/ ). The JHKs filters are those from
Cohen et al. (2003).
"""

from __future__ import print_function

from numpy import log10
from scipy.integrate import simps

import pycoco

__all__ = ['offset', 'convert_AB_to_Vega', 'convert_Vega_to_AB']

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

def convert_Vega_to_AB(phot_table, filters = False):
    """
    Parameters
    ----------
    Returns
    -------
    """

    return phot_table


def convert_AB_to_Vega():
    """
    Parameters
    ----------
    Returns
    -------
    """

    return phot_table


def load_vega(path = "/Users/berto/Code/verbose-enigma/pycoco/kcorr/data/alpha_lyr_stis_002.dat"):
    vega = pycoco.SpectrumClass()
    vega.load(path)

    return vega


def load_AB(path = "/Users/berto/Code/verbose-enigma/pycoco/kcorr/data/AB_pseudospectrum.dat"):
    vega = pycoco.SpectrumClass()
    vega.load(path)

    return vega


def calc_AB_flux(filter_name):

    AB = load_AB()

    filter_object = pycoco.load_filter("/Users/berto/Code/CoCo/data/filters/" + filter_name + ".dat")
    filter_object.resample_response(new_wavelength = AB.wavelength)

    transmitted_spec = filter_object.throughput * AB.flux

    integrated_flux = simps(transmitted_spec, AB.wavelength)

    return integrated_flux


def calc_AB_zp():

    integrated_flux = calc_AB_flux("BessellV")

    return -2.5 * log10(integrated_flux)


def calc_vega_flux(filter_name):

    vega = load_vega()

    filter_object = pycoco.load_filter("/Users/berto/Code/CoCo/data/filters/" + filter_name + ".dat")
    filter_object.resample_response(new_wavelength = vega.wavelength)

    transmitted_spec = filter_object.throughput * vega.flux

    integrated_flux = simps(transmitted_spec, vega.wavelength)

    return integrated_flux


def calc_vega_zp(filter_name, vega_Vmag = 0.03):

    integrated_V_flux = calc_vega_flux("BessellV")
    integrated_V_flux = calc_vega_flux("BessellV")

    return -2.5 * log10(integrated_flux) - vega_Vmag

def calc_vega_mag(filter_name):

    zp = calc_vega_zp()
    flux = calc_vega_flux(filter_name)

    mag = -2.5 * log10(flux) - zp
    return mag

# def calc_ve

def calc_offset_AB_minus_Vega(filter_name, vega_Vmag = 0.03):
    """

    """

    integrated_Vega_flux = calc_vega_flux(filter_name)
    integrated_AB_flux = calc_AB_flux(filter_name)

    AB  = -2.5*log10(integrated_AB_flux)
    Vega = -2.5*log10(integrated_Vega_flux) - vega_Vmag

    return AB - Vega

def convert_f_nu_to_f_lambda():
    pass
