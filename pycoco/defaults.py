"""

"""

import os

import astropy.units as u
import matplotlib.pyplot as plt

__all__ = ["_default_data_dir_path",
           "_default_filter_dir_path",
           "_default_coco_dir_path",
           "_default_recon_dir_path",
           "_default_specphase_dir_path",
           "_default_sn_dist_path",
           "_default_lcsim_path",
           "_default_list_dir_path",
           "_colourmap_name",
           "_spec_colourmap_name",
           "spec_colourmap",
           "_colour_upper_lambda_limit",
           "_colour_lower_lambda_limit",
           "_default_info_path",
           "_default_kcorr_data_path",
           "_default_lsst_throughputs_path"]

## Important variables

COCO_ROOT_DIR = os.environ["COCO_ROOT_DIR"]
LSST_THROUGHPUTS_ROOT = os.environ["LSST_THROUGHPUTS"]
SFD_DIR = os.environ["SFD_DIR"]

_default_data_dir_path = os.path.join(COCO_ROOT_DIR, "data/")
_default_filter_dir_path = os.path.join(COCO_ROOT_DIR, "data/filters/")
_default_list_dir_path = os.path.join(COCO_ROOT_DIR, "lists/")
_default_coco_dir_path = os.path.join(COCO_ROOT_DIR)
_default_recon_dir_path = os.path.join(COCO_ROOT_DIR, "recon/")
_default_specphase_dir_path = os.path.join(COCO_ROOT_DIR, "spectra/")
_default_sn_dist_path = os.path.join(COCO_ROOT_DIR, "sndist.list")
_default_lcsim_path = os.path.join(COCO_ROOT_DIR, "sim/")
_default_info_path = os.path.join(_default_data_dir_path, "info/info.dat")
_default_kcorr_data_path = os.path.join(os.path.abspath(os.path.join(__file__, os.pardir)), 'kcorr_data/')
_default_lsst_throughputs_path = os.path.abspath(LSST_THROUGHPUTS_ROOT)
_default_dust_dir = os.path.abspath(SFD_DIR)

# _colormap_name = 'jet'
# _colourmap_name = 'rainbow'
_spec_colourmap_name = 'viridis'
# _spec_colourmap_name = 'plasma'
# _spec_colourmap_name = 'jet'
_colourmap_name = 'plasma'

colourmap = plt.get_cmap(_colourmap_name)
spec_colourmap = plt.get_cmap(_spec_colourmap_name)

_colour_upper_lambda_limit = 11000 * u.angstrom
_colour_lower_lambda_limit = 3500 * u.angstrom
