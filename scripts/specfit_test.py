import pycoco as pcc
import pycoco.classes as classes
import pycoco.defaults as defaults
import os



verbose = True
fitpath = False
snname = "SN2007uy"
photpath = os.path.join(pcc.defaults._default_data_dir_path, "lc/" + snname + ".dat")


if snname:
    SNobject = classes.SNClass(snname)
    if not photpath:
        photpath = os.path.join(defaults._default_data_dir_path, "lc/" + snname + ".dat")

    SNobject.load_phot(path = photpath, verbose=verbose)

if not hasattr(SNobject, "lcfit"):
    if not fitpath:
        fitpath = os.path.join(defaults._default_recon_dir_path, snname + ".dat")
    SNobject.get_lcfit(fitpath, verbose=verbose)

manglefilters = [i for i in SNobject.lcfit.filter_names]
print(manglefilters)
if "BessellR" in manglefilters and "SDSS_r" in manglefilters:
    if verbose:
        print("Has both BessellR and SDSS_r - using SDSS_r")
    ## Only use SDSS_r - less hassle
    manglefilters.remove("BessellR")
if "BessellI" in manglefilters and "SDSS_i" in manglefilters:
    if verbose:
        print("Has both BessellI and SDSS_i - using SDSS_i")
    ## Only use SDSS_i - less hassle
    manglefilters.remove("BessellI")

filename = snname + "_m.dat"
outpath = SNobject.lcfit.recon_directory

if verbose: print(outpath, filename)
