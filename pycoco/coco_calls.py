"""

"""

import os
import sys
import warnings
import subprocess

from numpy import asarray

from .classes import *
from .defaults import *
from .utils import *

__all__ = ["test_LCfit",
           "run_LCfit",
           "test_specfit",
           "run_specfit",
           "specfit_sn",
           "run_LCfit_fileinput"
           ]

#  #------------------------------------#  #
#  # CoCo Functions                     #  #
#  #------------------------------------#  #


def test_LCfit(snname, coco_dir = _default_coco_dir_path,
               verbose = True):
    """
    Check to see if a fit has been done. Does this by
    looking for reconstructed LC files
    Parameters
    ----------
    Returns
    -------
    """

    # try:
    #     if not coco_dir:
    #         coco_dir = _default_coco_dir_path
    #
    # except:
    #     warnings.warn("Something funky with your input")

    check_dir_path(coco_dir)

    if verbose: print(coco_dir)

    try:
        path_to_test_dat = os.path.join(coco_dir, 'recon', snname + '.dat')
        path_to_test_stat = os.path.join(coco_dir, 'recon', snname + '.stat')

        for path in [path_to_test_stat, path_to_test_dat]:

            if os.path.isfile(os.path.abspath(path)):
                if verbose: print("Looks like you have done a fit, I found ", path )
                boolflag = True
            else:
                warnings.warn(os.path.abspath(path) +
                " not found. Have you done a fit?")
                boolflag = False

    except:

        warnings.warn("Failing gracefully. Can't find the droids you are looking for.")
        boolflag = False

    return boolflag


def run_LCfit(path, coco_dir = _default_coco_dir_path, verbose = True,):
    """
    Parameters
    ----------
    Returns
    -------
    """
    check_file_path(path)
    relist() ## Check filter file is up to date

    if verbose: print("Running CoCo lcfit on " + path)
    cwd = os.getcwd()
    os.chdir(coco_dir)
    subprocess.call([os.path.join(_default_coco_dir_path, "lcfit"), path])
    os.chdir(cwd)
    if verbose: print("Fit complete")
    pass


def run_LCfit_fileinput(listfile_path, coco_dir = _default_coco_dir_path, data_dir = _default_data_dir_path, verbose = True):
    """

    :param listfile:
    :param verbose:
    :return:
    """

    check_file_path(listfile_path)

    if verbose: print("Reading ", listfile_path)

    file_list = []

    with  open(listfile_path) as infile:
        for line in infile:
            if verbose: print(os.path.join(data_dir, os.pardir, line.strip("\n")))
            file_list.append(os.path.join(data_dir, os.pardir, line.strip("\n")))

    file_list = asarray(file_list)

    for lc_path in file_list:
        run_LCfit(lc_path, coco_dir=coco_dir, verbose=verbose)

    pass

def run_all_SNe_LCfit(verbose = True):
    """

    :param verbose:
    :return:
    """
    pass

def test_specfit(snname, coco_dir = False,
               verbose = True):
    """
    Check to see if a fit has been done. Does this by
    looking for reconstructed .spec filess
    Parameters
    ----------
    Returns
    -------
    """

    try:
        if not coco_dir:
            coco_dir = _default_coco_dir_path

    except:
        warnings.warn("Something funky with your input")

    check_dir_path(coco_dir)

    if verbose: print(coco_dir)

    try:
        ##
        path_to_test_dat = os.path.join(coco_dir, 'recon', snname + '.dat')
        path_to_test_stat = os.path.join(coco_dir, 'recon', snname + '.stat')
        ## NEED TO THINK OF THE BEST WAY TO DO THIS

        for path in [path_to_test_stat, path_to_test_dat]:

            if os.path.isfile(os.path.abspath(path)):
                if verbose: print("Looks like you have done a fit, I found ", path )
                boolflag = True
            else:
                warnings.warn(os.path.abspath(path) +
                " not found. Have you done a fit?")
                boolflag = False

    except:

        warnings.warn("Failing gracefully. Can't find the droids you are looking for.")
        boolflag = False

    return boolflag


def run_specfit(path, verbose = True):
    """
    runs CoCo specfit on the listfile supplied in path

    Parameters
    ----------
    Returns
    -------
    """
    check_file_path(path)
    relist() ## Check filter file is up to date

    if verbose: print("Running CoCo specfit on " + path)
    subprocess.call([os.path.join(_default_coco_dir_path, "./specfit"), path])

    pass


def specfit_sn(snname, verbose = True):
    """
    runs CoCo specfit on the listfile supplied in path.

    Parameters
    ----------

    snname -

    Returns
    -------

    """

    ## Need to look for the recon lc files for snname
    # sn = SNClass(snname)
    lcfit = LCfitClass()

    path = os.path.join(lcfit.recon_directory, snname+".dat")
    lcfit.load_formatted_phot(path)
    lcfit.unpack()
    lcfit._sort_phot()
    lcfit.get_fit_splines()

    ## Need to make new recon lc files for mangling - no overlaps
    # lcfit.
    manglefilters = [i for i in lcfit.filter_names]

    if "BessellR" and "SDSS_r" in manglefilters:
        ## Only use SDSS_r - less hassle
        manglefilters.remove("BessellR")
    if "BessellI" and "SDSS_i" in manglefilters:
        ## Only use SDSS_r - less hassle
        manglefilters.remove("BessellI")

    filename = snname + "_m.dat"
    outpath = lcfit.recon_directory

    lcfit.save(filename, path = outpath, filters = manglefilters, squash = True)

    # lcfit.
    ## Need to change the listfile to one that has snname matches the new lc file
    listpath = os.path.join(_default_coco_dir_path, "lists", snname + ".list")

    origlist = read_list_file(listpath)
    origlist.rename_column('snname', 'snname_nomangle')
    origlist["snname"] = [j+"_m" for j in origlist["snname_nomangle"]]

    newlistpath = os.path.join(_default_coco_dir_path, "lists", snname + "_m.list")
    newlist = origlist["spec_path", "snname", "mjd_obs", "z"]
    # newlist.write(newlistpath, format = "ascii.fast_commented_header", overwrite = True)
    newlist.write(newlistpath, format = "ascii.fast_commented_header")


    ## Need to call run_specfit with path of new list file

    run_specfit(newlistpath)

    pass
