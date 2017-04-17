"""

"""

import os
import warnings
import re

from .defaults import *
from .errors import *
from .filters import FilterClass
from .utils import *

# from .filters import FilterClass
# from .
##------------------------------------##
##  Functions                         ##
##------------------------------------##

def load_filter(path, cmap = False, verbose = False):
    """
    Loads a filter response into FilterClass and returns it.

    Parameters
    ----------
    Returns
    -------
    """

    if check_file_path(os.path.abspath(path)):
        filter_object = FilterClass()
        filter_object.read_filter_file(os.path.abspath(path), verbose = verbose)

        if cmap:
            filter_object.calculate_plot_colour(verbose = verbose)

        return filter_object
    else:
        warnings.warn("Couldn't load the filter")
        return None


def get_filter_from_filename(path, snname, file_type):
    """
    Parameters
    ----------
    Returns
    -------
    """
    filename_from_path = path.split("/")[-1]
    filename_no_extension = filename_from_path.replace(file_type, "")
    filter_string = filename_no_extension.replace(snname+"_", "")

    # phot_file.replace(file_type, '').split('_')[-1]

    return filter_string


def _get_filter_directory():
    """
    Get the defaul path to the filter directory.

    Looks for the filter directory set as environment variable
    $PYCOCO_FILTER_DIR. if not found, returns default.

    returns: Absolute path in environment variable $PYCOCO_DATA_DIR, or
             default datalocation: '../testdata/'.
    """

    return os.environ.get('PYCOCO_FILTER_DIR', _default_filter_dir_path)


def load_phot(path, names = ('MJD', 'flux', 'flux_err', 'filter'),
              format = 'ascii', verbose = True):
    """
    Loads a single photometry file.

    Parameters
    ----------
    Returns
    -------
    """

    StringWarning(path)

    # phot_table = ap.table.Table.read(path, format = format, names = names)
    phot_table = Table.read(path, format = format, names = names)

    phot_table.replace_column("MJD", Time(phot_table["MJD"], format = 'mjd'))

    phot_table["flux"].unit = u.cgs.erg / u.si.angstrom / u.si.cm ** 2 / u.si.s
    phot_table["flux_err"].unit =  phot_table["flux"].unit


    return phot_table


def load_formatted_phot(path, format = "ascii", names = False,
                        verbose = True):
    """
    Loads a single photometry file.

    Parameters
    ----------
    Returns
    -------
    """

    StringWarning(path)

    if names:
        phot_table = Table.read(path, format = format, names = names)
    else:
        phot_table = Table.read(path, format = format)

    phot_table.meta = {"filename" : path}

    phot_table["MJD"].unit = u.day
    phot_table["flux"].unit = u.cgs.erg / u.si.angstrom / u.si.cm ** 2 / u.si.s
    phot_table["flux_err"].unit =  phot_table["flux"].unit

    return phot_table


def load(path, format = "ascii", verbose = True):
    pc = PhotometryClass()
    pc.phot = load_formatted_phot(path, format = format, verbose = verbose)
    pc.unpack(verbose = verbose)
    return pc


def load_all_phot(path = _default_data_dir_path, format = "ascii", verbose = True):
    """
    loads photometry into AstroPy Table.

    returns: AstroPy table
    """
    ## Do i need the following? Errors should be handled in find_phot?
    # try:
    #     if os.path.isdir(os.path.abspath(path)):
    #         pass
    #     else:
    #         warnings.warn(os.path.abspath(data_dir_path) +
    #         " is not a valid directory. Returning 'False'.")
    # except:
    #     raise PathError("The data directory '" + path + "' doesn't exist.")
    phot_list = find_filter_phot(path = path)

    if len(phot_list) > 0:
        # phot_table = Table()
        phot_table = ap.table.Table()

        for phot_file in phot_list:
            print(phot_file)
            print(phot_table.read(phot_file, format = format))

        return phot_table
    else:
        warning.warn("Couldn't find any photometry")


def find_filter_phot(path = _default_data_dir_path, snname = False,
              prefix = 'SN', file_type = '.dat',
              verbose = True):
    """
    Tries to find photometry in the supplied directory.

    Looks in a directory for things that match SN*.dat. Uses regex via `re` -
    probably overkill.

    Parameters
    ----------

    path :

    snname :

    prefix :

    file_type :


    Returns
    -------

    phot_list :

    """
    # regex = re.compile("^SN.*.dat")

    StringWarning(path)
    if not check_dir_path(path):
        return False

    try:
        if snname:
            match_string = "^" + str(snname) + ".*" + '.dat'
        else:
            match_string = "^" + str(prefix) + ".*" + '.dat'
    except:
        raise TypeError

    regex = re.compile(match_string)

    ls = os.listdir(path)

    phot_list = [os.path.abspath(os.path.join(path, match.group(0))) for file_name in ls for match in [regex.search(file_name)] if match]

    if os.path.join(path, snname + file_type) in phot_list:
        phot_list.remove(os.path.join(path,snname + file_type))
        warnings.warn("Found " + os.path.join(path,snname + file_type) + " - you could just read that in.")

    if verbose:
        print("searching for", match_string)
        print("Found: ")
        print(ls)
        print("Matched:")
        print(phot_list)
    if len(phot_list) is 0:
        warnings.warn("No matches found.")
    return phot_list


def find_formatted_phot(path = _default_data_dir_path, snname = False,
              prefix = 'SN', file_type = '.dat',
              verbose = True):
    """
    Tries to find photometry in the supplied directory.

    Looks in a directory for things that match SNname.

    Originally did something else - this is a bit hacky.
    Parameters
    ----------

    path :

    snname :

    prefix :

    file_type :


    Returns
    -------

    phot_list :

    """
    # regex = re.compile("^SN.*.dat")

    StringWarning(path)
    if not check_dir_path(path):
        return False

    try:
        if snname:
            match_string = str(snname) + file_type
        else:
            warnings.warn("No SN name given")
            return 0
    except:
        raise TypeError
    #
    # regex = re.compile(match_string)

    ls = os.listdir(path)
    if verbose: print(ls)
    # phot_list = [os.path.abspath(os.path.join(path, match.group(0))) for file_name in ls for match in [regex.search(file_name)] if match]
    phot_list = [os.path.abspath(os.path.join(path, file_name)) for file_name in ls if file_name == match_string]

    if os.path.join(path, snname + file_type) in phot_list:
        # phot_list.remove(os.path.join(path,snname + file_type))
        warnings.warn("Found " + os.path.join(path,snname + file_type) + " - you could just read that in.")

    if verbose:
        print("searching for", match_string)
        print("Found: ")
        print(ls)
        print("Matched:")
        print(phot_list)

    if len(phot_list) is 0:
        warnings.warn("No matches found.")
    return phot_list


def find_recon_spec(snname, dir_path = _default_recon_dir_path, verbose = False):
    """

    Parameters
    ----------

    Returns
    -------
    """
    file_type = ".spec"
    StringWarning(dir_path)
    if not check_dir_path(dir_path):
        return False

    try:
        ls = np.array(os.listdir(dir_path))

        wspec = np.where(np.char.find(ls, file_type, start = -len(file_type)) > -1)
        spec_list = ls[wspec]

        ## The last 18 chars are for the MJD and file_type
        wsn = np.where([i[:-18] == snname for i in spec_list])
        snmatch_list = spec_list[wsn]

        if verbose:
            print("Found: ")
            print(ls)
            print("Spec:")
            print(spec_list)
            print("Matched:")
            print(snmatch_list)
        if len(snmatch_list) is 0:
            warnings.warn("No matches found.")
        return snmatch_list

    except:
        warnings.warn("Something went wrong")
        return False


def find_specphase_spec(snname, dir_path = _default_specphase_dir_path, file_type = ".spec", verbose = False):
    """

    Parameters
    ----------

    Returns
    -------
    """

    StringWarning(dir_path)
    if not check_dir_path(dir_path):
        return False

    try:
        ls = np.array(os.listdir(dir_path))

        # wspec = np.where(np.char.find(ls, file_type, start = -len(file_type)) > -1)
        # spec_list = ls[wspec]
        speclist = [i for i in ls if i[-5:] == ".spec"]
        ## The last 18 chars are for the MJD and file_type
        # wsn = np.where([i[:-18] == snname for i in spec_list])
        # snmatch_list = spec_list[wsn]
        snmatch_list = [i for i in speclist if i[:len(snname)] == snname ]

        if verbose:
            print("Found: ")
            print(ls)
            print("Spec:")
            print(spec_list)
            print("Matched:")
            print(snmatch_list)
        if len(snmatch_list) is 0:
            warnings.warn("No matches found.")
        return snmatch_list

    except:
        warnings.warn("Something went wrong")
        return False


# def check_url_status(url):
#     """
#     Snippet from http://stackoverflow.com/questions/6471275 .
#
#     Checks the status of a website - a status flag of < 400 means the site
#     is up.
#
#     """
#     p = urlparse(url)
#     conn = httplib.HTTPConnection(p.netloc)
#     conn.request('HEAD', p.path)
#     resp = conn.getresponse()
#
#     return resp.status
#
#
# def check_url(url):
#     """
#     Wrapper for check_url_status - considers the status, True if < 400.
#     """
#     return check_url_status(url) < 400


def setup_plot_defaults():
    """

    """

    plt.rcParams['ps.useafm'] = True
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 14
    plt.rcParams['figure.subplot.hspace'] = 0.1
    plt.rc('font', family='sans-serif')
    plt.rc('font', serif='Helvetica')
    pass


def read_list_file(path, names = ('spec_path', 'snname', 'mjd_obs', 'z'), verbose = True):
    """
    Parameters
    ----------
    Returns
    -------
    """
    check_file_path(path)
    #
    # ifile = open(path, 'r')
    #
    # for line in ifile:
    #     if verbose: print(line.strip('\n'))
    # ifile.close()
    data = Table.read(path, names = names, format = 'ascii')
    return data


def load_specfit(path):
    """
    Parameters
    ----------
    Returns
    -------
    """

    specfit = specfitClass()

    return specfit


def compare_spec(orig_spec, specfit,
                 xminorticks = 250, legend = True, verbose = True,
                 *args, **kwargs):
        """
        Parameters
        ----------
        Returns
        -------
        """

        if hasattr(orig_spec, "data") and hasattr(specfit, "data"):

            setup_plot_defaults()

            fig = plt.figure(figsize=[8, 4])
            fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                                right = 0.99, hspace=0, wspace = 0)

            ax1 = fig.add_subplot(111)

            if verbose: print(np.nanmean(specfit.flux), np.nanmean(orig_spec.flux))

            plot_label_string = r'$\rm{' + orig_spec.data.meta["filename"].split('/')[-1].replace('_', '\_') + '}$'

            ax1.plot(orig_spec.data['wavelength'], orig_spec.flux/np.nanmean(orig_spec.flux), lw = 2,
                         label = plot_label_string, color = 'Red',
                         *args, **kwargs)

            plot_label_string = r'$\rm{' + specfit.data.meta["filename"].split('/')[-1].replace('_', '\_') + '}$'


            ax1.plot(specfit.data['wavelength'], specfit.flux/np.nanmean(specfit.flux), lw = 2,
                         label = plot_label_string, color = 'Blue',
                         *args, **kwargs)

            maxplotydata = np.nanmax([specfit.flux/np.nanmean(specfit.flux), orig_spec.flux/np.nanmean(orig_spec.flux)])
            minplotydata = np.nanmin([specfit.flux/np.nanmean(specfit.flux), orig_spec.flux/np.nanmean(orig_spec.flux)])

            if legend:

                plot_legend = ax1.legend(loc = [1.,0.0], scatterpoints = 1,
                                      numpoints = 1, frameon = False, fontsize = 12)

            ax1.set_ylim(0, maxplotydata*1.02)

            ## Label the axes
            xaxis_label_string = r'$\textnormal{Wavelength (\AA)}$'
            yaxis_label_string = r'$\textnormal{Flux, erg s}^{-1}\textnormal{cm}^{-2}$'

            ax1.set_xlabel(xaxis_label_string)
            ax1.set_ylabel(yaxis_label_string)

            xminorLocator = MultipleLocator(xminorticks)
            ax1.xaxis.set_minor_locator(xminorLocator)

            plt.show()
        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")
        pass


# def load_stat(stats_path = '/Users/berto/Code/CoCo/chains/SN2011dh_Bessell/BessellB-stats.dat'):
#     verbose = False
#     fitparams = None
#     j = None
#     stat = None
#     modenumber = 0
#     modes = OrderedDict()
#
#     del j
#     del stat
#     del fitparams
#
#     # stats_path = '/Users/berto/Code/CoCo/chains/SN2006aj/B-stats.dat'
#
#
#     stat =  OrderedDict()
#     nparams = 8
#
#     key_string = 'Params'
#
#     ifile = open(stats_path, 'r')
#
#     for i, line in enumerate(ifile):
#         if line != '\n':
#             if verbose: print(i, line)
#             if line[:17] == "Total Modes Found":
#                 nmodes = int(line.split()[-1])
#
#             if line[:4] == "Mode":
#                 modenumber = modenumber + 1
#                 fitparams = OrderedDict()
#             if line[:8] == 'Strictly':
#                 try:
#                     stat["SLLE"] = np.float64(line.split()[3])
#                     stat["SLLE_err"] = np.float64(line.split()[5])
#                 except:
#                     if verbose: print("NOPE1")
#                 if verbose: print("BAR")
#             if line[:9] == "Local Log":
#                 try:
#                     stat["LLE"] = np.float64(line.split()[2])
#                     stat["LLE_err"] = np.float64(line.split()[4])
#                 except:
#                     if verbose: print("NOPE2")
#                 if verbose: print("SPAM")
#
#             if line[:3] == "MAP":
#                 key_string = key_string + "-MAP"
#                 fitparams = OrderedDict()
#
#             if line[:7] == "Maximum":
#                 key_string = key_string + "-ML"
#                 fitparams = OrderedDict()
#
#             if line[:7] == 'Dim No.':
#                 j = i + nparams+1
#                 if verbose: print("FOO")
#                 try:
#                     fitparams["Dim No."] = np.array([])
#                     fitparams["Mean"] = np.array([])
#                     if len(line.split()) > 2:
#                         fitparams["Sigma"] = np.array([])
#                 except:
#                     if verbose: print("NOPE3")
#             try:
#                 if verbose: print(i, j)
#                 if i<j and line[:7] != 'Dim No.' :
#
#                     if verbose: print("READ")
#                     if verbose: print(len(line.split()) ,line.split())
#                     fitparams["Dim No."] = np.append(fitparams["Dim No."], int(line.split()[0]))
#                     fitparams["Mean"] = np.append(fitparams["Mean"], np.float64(line.split()[1]))
#                     if len(line.split()) > 2:
#                         if verbose: print("GT 2")
#                         if verbose: print(np.float64(line.split()[2]))
#                         fitparams["Sigma"] = np.append(fitparams["Sigma"], np.float64(line.split()[2]))
#                         if verbose: print("fitparams sigma", fitparams["Sigma"])
#                 if i > j:
#     #                 print("SAVE")
#                     if verbose: print("SAVE")
#                     if verbose: print(key_string)
#                     stat[key_string] = fitparams
#
#                     if key_string[-3:] ==  "MAP":
#                         modes["Mode"+str(modenumber)] = stat
#                     key_string = 'Params'
#
#             except:
#                 if verbose: print("NOPE4")
#                 pass
#     #         if line
#
#     ifile.close()
#     print(nmodes)
#     print(stat.keys())
#     return modes


def filter_within_spec(filter_obj, spec_obj):
    """
    returns true if filter_edges are within spectrum, False otherwise

    Parameters
    ----------

    Returns
    -------
    """
    try:
        if hasattr(filter_obj, "_lower_edge") and hasattr(filter_obj, "_upper_edge") and hasattr(spec_obj, "data"):
            blue_bool = filter_obj._lower_edge > spec_obj.min_wavelength
            red_bool = filter_obj._upper_edge < spec_obj.max_wavelength

            if blue_bool and red_bool:
                return True
            else:
                return False
        else:
            warnings.warn("Filter object has no edges or spectrum object has no data")
            return False
    except:
        raise StandardError


def list_lcfits(verbose = False):
    """

    Parameters
    ----------

    Returns
    -------
    """
    fits_in_recon = [i for i in os.listdir(_default_recon_dir_path) if i[-4:] == ".dat"]
    if verbose: print(fits_in_recon)
    return fits_in_recon


def list_lcs(verbose = False):
    """

    Parameters
    ----------

    Returns
    -------
    """
    lcs_in_data_dir = [i for i in os.listdir(os.path.join(_default_data_dir_path, "lc")) if i[-4:] == ".dat"]
    if verbose: print(lcs_in_data_dir)
    return lcs_in_data_dir


def read_sndist_file(path = _default_sn_dist_path, format = "ascii"):
    """

    Parameters
    ----------

    Returns
    -------
    """
    check_file_path(path)
    table = Table.read(path, format = format)
    return table


def load_sndist(snname, *args, **kwargs):
    """

    Parameters
    ----------

    Returns
    -------
    """
    sndistlist = read_sndist_file(*args, **kwargs)

    try:
        w = np.where(sndistlist["snname"] == snname)
        row = sndistlist[w]
    except:
        warnings.warn("Failed to find distance info for " + snname + ". is it in the list?")
    return row
