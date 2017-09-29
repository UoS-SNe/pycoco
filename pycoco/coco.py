"""

"""

import os
import sys
import warnings
import subprocess
import copy

import numpy as np

from . import classes
from . import defaults
from . import functions
from . import utils
from . import kcorr

__all__ = [
            "test_LCfit",
            "run_LCfit",
            "run_LCfit_fileinput",
            "run_all_SNe_LCfit",
            "test_specfit",
            "run_specfit",
            "run_cocospecfit",
            "get_all_spec_lists",
            "specfit_all",
            "specfit_sn",
            "run_specphase",
            "check_specphase",
           ]


#  #------------------------------------#  #
#  # CoCo Functions                     #  #
#  #------------------------------------#  #


def test_LCfit(snname, coco_dir = defaults._default_coco_dir_path,
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
    #         coco_dir = defaults._default_coco_dir_path
    #
    # except:
    #     warnings.warn("Something funky with your input")

    utils.check_dir_path(coco_dir)

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


def run_LCfit(path, coco_dir = defaults._default_coco_dir_path, model = False,
              verbose = True,):
    """

    :param path:
    :param coco_dir:
    :param model:
    :param verbose:
    :return:
    """

    utils.check_file_path(path)
    utils.relist() ## Check filter file is up to date

    if model:
        models = np.unique([i.split(".")[0] for i in os.listdir(os.path.join(defaults._default_coco_dir_path, "src/models"))])

        if model not in models:
            return False
        callargs = [os.path.join(defaults._default_coco_dir_path, "lcfit"), path, "-m", model]
        print("running with", model)
    else:
        print("No Model supplied - running with default")
        callargs = [os.path.join(defaults._default_coco_dir_path, "lcfit"), path]
    if verbose: print("Running CoCo lcfit on " + path)
    if verbose: print("callargs are ", callargs)

    cwd = os.getcwd()
    os.chdir(coco_dir)
    subprocess.call(callargs)
    os.chdir(cwd)
    if verbose: print("Fit complete")
    pass


def run_LCfit_fileinput(listfile_path, coco_dir = defaults._default_coco_dir_path, data_dir = defaults._default_data_dir_path,
                        verbose = True):
    """

    :param listfile_path:
    :param coco_dir:
    :param data_dir:
    :param verbose:
    :return:
    """

    utils.check_file_path(listfile_path)

    if verbose: print("Reading ", listfile_path)

    file_list = []

    with  open(listfile_path) as infile:
        for line in infile:
            if verbose: print(os.path.join(data_dir, os.pardir, line.strip("\n")))
            file_list.append(os.path.join(data_dir, os.pardir, line.strip("\n")))

    file_list = np.asarray(file_list)

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
            coco_dir = defaults._default_coco_dir_path

    except:
        warnings.warn("Something funky with your input")

    utils.check_dir_path(coco_dir)

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


def run_specfit(SNObject, wantedfilters=False, anchor_distance=1000, save=True, plot = False, coco_dir=defaults._default_coco_dir_path, verbose = True):
    """
    replacement for `run_cocospecfit`. Mangles the spectra in the listfiles. Built for comfort, not speed.

    Parameters
    ----------
    Returns
    -------
    """

    if not wantedfilters:
        wantedfilters = SNObject.phot.filter_names.data

    outfile_log = []
    if hasattr(SNObject, "spec") and hasattr(SNObject, "lcfit"):
        for name, mS in SNObject.spec.items():
            #     print(name, mS)

            S = copy.deepcopy(mS)
            fit_dict = kcorr.mangle(SNObject, mS, mS.mjd_obs, wantedfilters, anchor_distance=anchor_distance)
            if plot:
                functions.plot_mangledata(S, fit_dict["data_table"], mS=fit_dict["SpectrumObject"], spl=fit_dict["final_spl"],
                                    show_linear_extrap=True, normalise=True)

            outfile = SNObject.name + "_" + str(fit_dict["SpectrumObject"].mjd_obs).ljust(12, "0") + ".spec"

            while outfile in outfile_log:
                j = 1
                outfile = SNObject.name + "_" + str(fit_dict["SpectrumObject"].mjd_obs + j * 0.00001).ljust(12, "0") + ".spec"
                j += 1
            print(outfile)
            outfile_log.append(outfile)

            if save:
                kcorr.save_mangle(fit_dict["SpectrumObject"], outfile, fit_dict["SpectrumObject"].infile)
    else:
        print("SNObject needs lcfit and spectra")
    # utils.check_file_path(path)
    # utils.relist() ## Check filter file is up to date
    # cwd = os.getcwd()
    # os.chdir(coco_dir)
    # if verbose: print("Running CoCo specfit on " + path)
    # subprocess.call([os.path.join(defaults._default_coco_dir_path, "./specfit"), path])
    # os.chdir(cwd)
    pass


def run_cocospecfit(path, coco_dir=defaults._default_coco_dir_path, verbose = True):
    """
    runs CoCo specfit on the listfile supplied in path

    Parameters
    ----------
    Returns
    -------
    """
    utils.check_file_path(path)
    utils.relist() ## Check filter file is up to date
    cwd = os.getcwd()
    os.chdir(coco_dir)
    if verbose: print("Running CoCo specfit on " + path)
    subprocess.call([os.path.join(defaults._default_coco_dir_path, "./specfit"), path])
    os.chdir(cwd)
    pass


def get_all_spec_lists(dirpath = defaults._default_list_dir_path, verbose=False):
    ignore = [".DS_Store", "master.list", "lightcurves.list"]

    utils.check_dir_path(dirpath)

    if verbose: print(dirpath)

    listfile_list = [i for i in os.listdir(dirpath) if i not in ignore]

    fullpath_list = [os.path.join(dirpath, j) for j in listfile_list]

    return fullpath_list


def specfit_all(verbose=True, dirpath=defaults._default_list_dir_path):

    fullpath_list = get_all_spec_lists(dirpath)

    for i, path in enumerate(fullpath_list):

        if verbose: print(i, path)

        run_specfit(path, verbose=verbose)

        if verbose: print("Done")

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
    # sn = classes.SNClass(snname)
    lcfit = classes.LCfitClass()

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
    listpath = os.path.join(defaults._default_coco_dir_path, "lists", snname + ".list")

    origlist = utils.read_list_file(listpath)
    origlist.rename_column('snname', 'snname_nomangle')
    origlist["snname"] = [j+"_m" for j in origlist["snname_nomangle"]]

    newlistpath = os.path.join(defaults._default_coco_dir_path, "lists", snname + "_m.list")
    newlist = origlist["spec_path", "snname", "mjd_obs", "z"]
    # newlist.write(newlistpath, format = "ascii.fast_commented_header", overwrite = True)
    newlist.write(newlistpath, format = "ascii.fast_commented_header")


    ## Need to call run_specfit with path of new list file

    run_specfit(newlistpath)

    pass


def run_specphase(filtername, phase_path, filetype=".dat", coco_dir=defaults._default_coco_dir_path, verbose = True):
    """
    runs CoCo specphase.

    Parameters
    ----------
    Returns
    -------
    """
    filters = utils._get_current_filter_registry()

    if filtername+filetype not in filters:
        warnings.warn("Filtername not recognised")
        return False

    utils.check_file_path(phase_path)

    utils.relist() ## Check filter file is up to date
    cwd = os.getcwd()
    os.chdir(coco_dir)
    # if verbose: print("Running CoCo specfit on " + path)
    subprocess.call([os.path.join(defaults._default_coco_dir_path, "./specphase"), phase_path, filtername])
    os.chdir(cwd)
    pass


def check_specphase(snname, spectra_dir="spectra/", coco_dir=defaults._default_coco_dir_path, absolute_path=False, verbose = False):
    """

    :param verbose:
    :param snname:
    :param spectra_dir:
    :param coco_dir:
    :param absolute_path:
    :return:
    """
    if not absolute_path:
        spectra_dir = os.path.join(coco_dir, spectra_dir)

    utils.check_dir_path(spectra_dir)

    dir_list = [spec for spec in os.listdir(spectra_dir) if spec[:len(snname)] == snname]
    if verbose:
        print(len(dir_list), "files found matching", snname)
        print(dir_list)
    return dir_list

