"""
This contains the extinction sub-module for pycoco
"""
from __future__ import print_function

import numpy as np

__all__ = ['extinction_law', 'unred']

def extinction_law(a, b, Rv = 3.1):
    """Eqn 1 from Cardelli 1989"""

    a_lam_aV = a + b / Rv

    return a_lam_aV

def coeffs(x, wavl = False, verbose = False, cardelli = True, odonnell = False):
    """ x = 1 / lambda """

    # wavl = np.array(wavl)
    # x = np.array(1./wavl)
    x = np.array(x)
    if verbose: print(x)

    a = np.array([])
    b = np.array([])

    pass_x = np.array([])

    for i in enumerate(x):

        ## IR
        if i[1] >= 0.3 and i[1] <= 1.1:
            if verbose: print('IR')
            a = np.append(a, 0.574*np.power(i[1], 1.61))
            b = np.append(b, -0.527*np.power(i[1], 1.61))

            pass_x= np.append(pass_x, i[1])

        ## Near-IR/OPT
        if i[1] > 1.1 and i[1] <= 3.3:
            if verbose: print('Near-IR/Optical')
            y = i[1] - 1.82
            if cardelli:
                a_coefficients = np.array([1., 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530, 0.32999][::-1]) ## use [::-1] to reverse
                b_coefficients = np.array([0., 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260, -2.09002][::-1])
            if odonnell:
                a_coefficients = np.array([1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505][::-1])        ##from O'Donnell
                b_coefficients = np.array([0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347][::-1])

            apoly = np.poly1d(a_coefficients)
            bpoly = np.poly1d(b_coefficients)

            a = np.append(a, apoly(y))
            b = np.append(b, bpoly(y))
            # a_i = 1. + 0.17699 - 0.50447 - 0.02427 + 0.72085 + 0.01979 - 0.77530 + 0.32999
            pass_x= np.append(pass_x, i[1])

    return {'a': a, 'b': b, 'x': pass_x}

def angstrom_to_inv_micron(wavl):
    wavl_m = wavl*1e-10
    wavl_mu = wavl_m/1e-6
    inv_micron = 1./wavl_mu
    return inv_micron

def inv_micron_to_angstrom(x):
    return 10000./x

def unred(wave, flux, wav_in_m = False, r_v = 3.1, EBV_MW = False, EBV_host = False):
    """

    """

    if wav_in_m:
        wav_inv = 1./(1e6*wave)
    else:
        wav_inv = angstrom_to_inv_micron(wave)

    EBV = EBV_MW + EBV_host

    A_V = r_v*EBV
    vals = coeffs(wav_inv)

    A_lambda = A_V *  (vals['a'] + vals['b']/r_v)
    funred = flux * np.power(10.,(0.4 * A_lambda))

    return funred

def ccm(wave, wav_in_m = False, r_v = 3.1, EBV_MW = False, EBV_host = False, return_frac = True):
    """

    """

    if wav_in_m:
        wav_inv = 1./(1e6*wave)
    else:
        wav_inv = angstrom_to_inv_micron(wave)

    EBV = EBV_MW + EBV_host

    A_V = r_v*EBV
    vals = coeffs(wav_inv)

    A_lambda = A_V*(vals['a'] + vals['b']/r_v)
    A_lambda_A_v =  (vals['a'] + vals['b']/r_v)
    if return_frac:
        return A_lambda_A_v, inv_micron_to_angstrom(vals['x'])
    else:
        return A_lambda
