from __future__ import print_function ## Force python3-like printing

try:
    from importlib import reload
except:
    pass

from matplotlib import pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

import os
import numpy as np
from astropy.table import Table,Column

import pycocosn as pcc

reload(pcc) ## FOR DEV
reload(pcc.defaults)
reload(pcc.functions)
reload(pcc.classes)

import pyCoCo as pccsims

def convert_column_string_encoding(column):
    column = Column([pcc.utils.b(x) for x in column.data], name = column.name)
    return column

def get_mjdmax_BessellV(sn):
    v = sn.lcfit.spline["BessellV"]
    mjd_spline = np.arange(np.nanmin(sn.phot.data["BessellV"]["MJD"]),
                 np.nanmax(sn.phot.data["BessellV"]["MJD"]),
                 0.001)
    w = np.where(v(mjd_spline) == np.nanmax(v(mjd_spline)))

    mjdmax = mjd_spline[w]

    return mjdmax

if __name__ == "__main__":
    print(pccsims.__file__)

    filter_path = "/Users/berto/Code/CoCo/data/filters"
    coco_root_path = "/Users/berto/Code/CoCo"
    coco = pccsims.pyCoCo(pcc.utils.b(filter_path), pcc.utils.b(coco_root_path))

    # snname = "SN2007uy"
    snname = "SN2009jf"

    sn = pcc.SNClass(snname)

    phot_path = os.path.join(coco_root_path, "data/lc/", snname + ".dat")
    speclist_path = os.path.join(str(coco_root_path),"lists/" + snname + ".list")
    recon_filename = os.path.abspath(os.path.join(str(coco_root_path), "recon/", snname + ".dat"))

    print(phot_path)
    sn.load_phot(path = phot_path)
    # sn.phot.plot()
    sn.get_lcfit(recon_filename)

    sn.load_list(path = speclist_path)
    sn.load_spec()
    # sn.load_mangledspec()
    # sn.plot_spec()
    # sn.plot_mangledspec()
    sn.plot_lc(multiplot = False)

    sn.lcfit.get_fit_splines()

    plt.close()

    plt.plot(sn.phot.data["BessellV"]["MJD"], sn.lcfit.spline["BessellV"](sn.phot.data["BessellV"]["MJD"]), label = r"$\textnormal{Spline}$")
    plt.scatter(sn.phot.data["BessellV"]["MJD"], sn.phot.data["BessellV"]["flux"], label = r"$\textnormal{Photometry}$")
    plt.plot(sn.lcfit.data["BessellV"]["MJD"], sn.lcfit.data["BessellV"]["flux"], label = r"$\textnormal{Fit}$")
    plt.legend()

    mjdmax = get_mjdmax_BessellV(sn)[0]

    filters_to_sim = convert_column_string_encoding(sn.phot.phot["filter"]).data
    mjd_to_sim = sn.phot.phot["MJD"].data

    print(mjdmax)
    print(mjd_to_sim)
    print(filters_to_sim)

    tablepath = "/Users/berto/Code/verbose-enigma/testdata/info/info.dat"
    info = Table.read(tablepath, format = "ascii.commented_header")

    z_obs = info[np.where(info["snname"] == "SN2009jf")]["z_obs"].data[0]
    print(z_obs)

    flux, flux_err = coco.simulate(b"SN2009jf",
                    0.008, 0.0, 0.00001, 0.00001, 3.1,
                    mjdmax, mjd_to_sim,
                    filters_to_sim)

    p = pcc.PhotometryClass()

    p.load_table(pcc.utils.simulate_out_to_ap_table(mjd_to_sim, flux, flux_err, filters_to_sim))

    p.plot(legend=True)

else:
    pass
