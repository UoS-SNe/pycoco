'''
This is the module for the PyCoCo python tools.

author: Rob Firth, Southampton
date: 06-12-2016
'''

from __future__ import print_function ## Force python3-like printing

if __name__ is not '__main__':

    __name__ = 'pycoco'
    __version__ = 0.1

try:
    __file__

except NameError:

    __file__ = sys.argv[0]

import os
import warnings
import unittest
import httplib
import re
from urlparse import urlparse
from collections import OrderedDict

import astropy as ap
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d as interp1d


##----------------------------------------------------------------------------##
##                                   TOOLS                                    ##
##----------------------------------------------------------------------------##

##------------------------------------##
##  DUMMY CODE                        ##
##------------------------------------##

class CustomValueError(ValueError):
	"""
	Raise when....
	"""


	def __init__(self, *args, **kwargs):
		ValueError.__init__(self, *args, **kwargs)


class DummyClass():
    '''
    Quick dummy class.

    Contains a test class variable and test class method that prints the
    variable.

    RF
    '''


    def __init__(self):
        self.dummy_string = 'Hello, World!'


    def print_dummy_string(self):
        print(self.test_string)


def dummy_function(*args, **kwargs):
    '''
    Quick dummy function.

    Prints supplied **args and **kwargs
    Issues warnings if nothing passed

    RF
    '''

    warnings.simplefilter('always')
    print(args)
    print(kwargs)


    # warnings.warn("WARNING")

    if not args and not kwargs:
        warnings.warn( "You didn't pass any *args or **kwargs", RuntimeWarning)

    else:
        if args:
            for i, arg in enumerate(args):
                print('an arg passed via *args: ', repr(arg))
        else:
            warnings.warn( "You didn't pass any *args", RuntimeWarning)

        if kwargs:
            for key, value in kwargs.iteritems():
                print('a **kwarg: ', repr(key), ' == ' , repr(value))
        else:
            warnings.warn( "You didn't pass any **kwargs", RuntimeWarning)
    pass


_somevar = 'Foo'


##----------------------------------------------------------------------------##
##  CODE                                                                      ##
##----------------------------------------------------------------------------##

__all__ = ["_default_data_dir_path", "_default_filter_dir_path", "_colourmap_name", "_colour_upper_lambda_limit", "_colour_lower_lambda_limit"]

## Important variables

_default_data_dir_path = os.path.abspath(os.path.join(__file__, os.pardir, os.pardir) + '/testdata/')
_default_filter_dir_path = "/Users/berto/Code/CoCo/data/filters/"
# _colormap_name = 'jet'
_colourmap_name = 'rainbow'
colourmap = plt.get_cmap(_colourmap_name)

_colour_upper_lambda_limit = 10000 * u.angstrom
_colour_lower_lambda_limit = 3000 * u.angstrom

##------------------------------------##
##  ERROR DEFS                        ##
##------------------------------------##


class CustomValueError(ValueError):
	"""
	Raise when....
	"""
	def __init__(self, *args, **kwargs):
		ValueError.__init__(self, *args, **kwargs)


class PathError(StandardError):
	"""
	Raise when a path is found to be invalid
	"""
	def __init__(self, *args, **kwargs):
		StandardError.__init__(self, *args, **kwargs)


class FilterMismatchError(ValueError):
	"""
	Raise when a Filter from filename doesn't match the one in the photfile
	"""
	def __init__(self, *args, **kwargs):
		ValueError.__init__(self, *args, **kwargs)


def StringWarning(path):
    """

    """
    if type(path) is not str and type(path) is not unicode:
        warnings.warn("WARNING: You passed something that was " + str(type(path)) + "This might go wrong.",
                      stacklevel = 2)

    else:
        pass



##------------------------------------##
##                                    ##
##------------------------------------##

class PhotometryClass():
    """
    Probably also overkill - but should be easier to store metadata etc. Hopefully
    flexible enough to just be a wrapper for AP tables of phot.

    PhotometryClass should have a FilterClass method describing the observations.
    """

    def __init__(self, verbose = False):
        """

        """

        ## Initialise the class variables
        self._default_data_dir_path = _default_data_dir_path
        self._default_filter_dir_path = _default_filter_dir_path
        self.data = OrderedDict()
        self.data_filters = OrderedDict()

        ## Initialise using class methods
        self.set_data_directory(self._get_data_directory())
        self.set_filter_directory(self._get_filter_directory())


    def _get_filter_directory(self):
        """
        Get the default path to the filter directory.

        Looks for the filter data directory set as environment variable
        $PYCOCO_FILTER_DIR. if not found, returns default.

        returns: Absolute path in environment variable $PYCOCO_FILTER_DIR, or
                 default datalocation: '/Users/berto/Code/CoCo/data/filters/'.
        """
        return os.path.abspath(os.environ.get('PYCOCO_FILTER_DIR', self._default_filter_dir_path))


    def _get_data_directory(self):
        """
        Get the default path to the data directory.

        Looks for the data data directory set as environment variable
        $PYCOCO_DATA_DIR. if not found, returns default.

        returns: Absolute path in environment variable $PYCOCO_DATA_DIR, or
                 default datalocation: '../testdata/'.
        """

        return os.path.abspath(os.environ.get('PYCOCO_DATA_DIR', self._default_data_dir_path))


    def set_filter_directory(self, filter_dir_path = '', verbose = False):
        """
        Set a new filter directory path.

        Enables the data directory to be changed by the user.

        """
        try:
            if os.path.isdir(os.path.abspath(filter_dir_path)):
                self.filter_directory = os.path.abspath(filter_dir_path)
                pass
            else:
                warnings.warn(os.path.abspath(filter_dir_path) +
                " is not a valid directory. Restoring default path: " +
                self._default_filter_dir_path, UserWarning)
                self.data_directory = self._default_data_dir_path

                if not os.path.isdir(self.filter_directory):
                    if verbose: print(os.path.isdir(self.filter_directory))
                    raise PathError("The default data directory '" + self.filter_directory
                     + "' doesn't exist. Or isn't a directory. Or can't be located.")
                else:
                    pass
        except:
            if verbose: print("foo")
            raise PathError("The default filter directory '" + self._default_filter_dir_path
             + "' doesn't exist. Or isn't a directory. Or can't be located. Have"
             + " you messed with _default_filter_dir_path?")
            pass


    def set_data_directory(self, data_dir_path = '', verbose = False):
        """
        Set a new data directory path.

        Enables the data directory to be changed by the user.

        """
        try:
            if os.path.isdir(os.path.abspath(data_dir_path)):
                self.data_directory = os.path.abspath(data_dir_path)
                pass
            else:
                warnings.warn(os.path.abspath(data_dir_path) +
                " is not a valid directory. Restoring default path: " +
                self._default_data_dir_path, UserWarning)
                self.data_directory = self._default_data_dir_path

                if not os.path.isdir(self.data_directory):
                    if verbose: print(os.path.isdir(self.data_directory))
                    raise PathError("The default data directory '" + self.data_directory
                     + "' doesn't exist. Or isn't a directory. Or can't be located.")
                else:
                    pass
        except:
            if verbose: print("foo")
            raise PathError("The default data directory '" + self._default_data_dir_path
             + "' doesn't exist. Or isn't a directory. Or can't be located. Have"
             + " you messed with _default_data_dir_path?")
            pass


    def load_phot_from_file(self, path, names = ('MJD', 'flux', 'flux_err', 'filter'),
                  format = 'ascii', verbose = True):
        """

        """
        phot_table = load_phot(path, names = names, format = format, verbose = verbose)
        self.data[np.unique(phot_table["filter"])[0]] = phot_table
        # self.photometry
        # phot_list = find_phot(self.data_directory)
        # print(phot_list)

        pass


    def load_phot_ap_tables(self):
        """

        """

        pass


    def load(self, path = _default_data_dir_path, snname = False, prefix = 'SN',
             file_type = '.dat', names = ('MJD', 'flux', 'flux_err', 'filter'),
             format = 'ascii', filter_file_type = '.dat', verbose = True):
        """
        Finds and loads in data (from file) into phot objects.

        Parameters
        ----------

        Returns
        -------

        """

        if snname:
            ## Find matching photometry
            phot_list = find_phot(path = path, snname = snname, prefix = prefix,
                                  file_type = file_type, verbose = verbose)

            ## Loop over files (shouldn't be that many really)
            if len(phot_list) > 0:

                for phot_file in phot_list:

                    if verbose: print(phot_file)
                    phot_table = Table.read(phot_file, names = names, format = format)

                    filter_string = get_filter_from_filename(phot_file, snname, file_type)
                    phot_table.meta = {"filename" : phot_file,
                                       "filter" : filter_string,
                                       "filter_filename": filter_string + filter_file_type}

                    ## Sort out units
                    phot_table.replace_column("MJD", Time(phot_table["MJD"], format = 'mjd'))

                    phot_table["flux"].unit = u.cgs.erg / u.si.angstrom / u.si.cm ** 2 / u.si.s
                    phot_table["flux_err"].unit =  phot_table["flux"].unit

                    ## Put in dictionary - use filter from the file
                    filter_key = np.unique(phot_table["filter"])[0]
                    if verbose: print(len(np.unique(phot_table["filter"])) , phot_table.meta["filter"], filter_key)

                    if len(np.unique(phot_table["filter"])) > 1 or filter_key != phot_table.meta["filter"]:
                        raise FilterMismatchError("There is a mismatch between the filter filename and that in the "
                                                   + "photometry file")

                    self.data[filter_key] = phot_table

                    path_to_filter = os.path.join(self.filter_directory, phot_table.meta['filter_filename'])
                    self.data_filters[filter_key] = load_filter(path_to_filter)

                # return phot_table
            else:
                warning.warn("Couldn't find any photometry")
        else:
            warnings.warn("Provide a SN name")

        pass


    def _combine_phot(self, verbose = True):
        """

        """

        if hasattr(self, "data"):
            if verbose: print(self.data.keys())

            for i, phot_filter in enumerate(self.data.keys()):
                print(i, phot_filter)



        else:
            warnings.warn("Cant find self.data")



        pass


    def save(self, verbose = True):
        """
        Outputs

        Parameters
        ----------
        Returns
        -------
        """

        if hasattr(self, "data"):
            if verbose: print("has data")


        pass


    def plot(self):
        """

        """


        pass


    def save_phot(self, path):
        """
        Output the photometry loaded into the SNClass via self.load_phot* into a format
        and location recognised by CoCo.

        Parameters
        ----------

        Returns
        -------

        """


        pass


class FilterClass():
    """Docstring for FilterClass"""

    def __init__(self, verbose = True):
        self._wavelength_units = u.Angstrom
        self._wavelength_units._format['latex'] = r'\rm{\AA}'
        pass


    def read_filter_file(self, path):
        """
        Assumes Response function is fractional rather than %.
        """
        if check_file_path(os.path.abspath(path)):
            self.wavelength, self.throughput = np.loadtxt(path).T

            self._filter_file_path = path

            filename = path.split('/')[-1]
            filename_no_extension = filename.split('.')[0]
            self.filter_name = filename_no_extension

            self.calculate_effective_wavelength()

        else:
            warnings.warn("Foo")


    def calculate_effective_wavelength(self):
        """
        Well, what are you expecting something called `calculate_effective_wavelength`
         to do?
        """

        spline_rev = interp1d((np.cumsum(self.wavelength*self.throughput)/np.sum(self.wavelength*self.throughput)), self.wavelength)
        lambda_eff = spline_rev(0.5)

        self.lambda_effective = lambda_eff * self._wavelength_units


    def plot(self, xminorticks = 250, yminorticks = 0.1, *args, **kwargs):
        """
        Plots filter throughput, so you can double check it.

        Parameters
        ----------

        Returns
        -------
        """

        ## Check if there is something in the class to plot
        if hasattr(self, "wavelength") and hasattr(self, "throughput"):

            setup_plot_defaults()
            xaxis_label_string = r'$\textnormal{Wavelength, ' + self._wavelength_units.name + ' (}' + self._wavelength_units._format['latex'] +')$'
            yaxis_label_string = r'$\textnormal{Fractional Throughput}$'

            plot_label_string = r'$' + self.filter_name + '$'

            yminorLocator = MultipleLocator(yminorticks)
            xminorLocator = MultipleLocator(xminorticks)

            fig = plt.figure(figsize=[8, 4])
            fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                                right = 0.99, hspace=0, wspace = 0)

            ax1 = fig.add_subplot(111)

            if hasattr(self, "_plot_colour"):
                ax1.plot(self.wavelength, self.throughput, color = self._plot_colour,
                         lw = 2, label = plot_label_string)
            else:
                ax1.plot(self.wavelength, self.throughput, lw = 2, label = plot_label_string)

            ax1.set_xlabel(xaxis_label_string)
            ax1.set_ylabel(yaxis_label_string)

            ax1.yaxis.set_minor_locator(yminorLocator)
            ax1.xaxis.set_minor_locator(xminorLocator)

            ax1.legend(loc = 0)

            plt.show()
            pass
        else:
            warning.warn("Doesn't look like you have loaded a filter into the object")


    def resample_response(self, new_wavelength):
        """
        Bit dodgy - spline has weird results for poorly sampled filters

        Parameters
        ----------

        Returns
        -------
        """

        if hasattr(self, "wavelength") and hasattr(self, "throughput"):
            self._wavelength_orig = self.wavelength
            self._throughput_orig = self.throughput

            self.wavelength = np.concatenate(([0,1], self._wavelength_orig, [24999,25000]))
            self.throughput = np.concatenate(([0,0], self._throughput_orig, [0,0]))

            interp_func = InterpolatedUnivariateSpline(self.wavelength, self.throughput)
            self.throughput = interp_func(new_wavelength)
            self.wavelength = new_wavelength

            self.throughput[np.where(self.throughput < 0.0)] = 0.0
        else:
            warning.warn("Doesn't look like you have loaded a filter into the object")


    def calculate_plot_colour(self, colourmap = colourmap, verbose = True):
        """


        Parameters
        ----------

        Returns
        -------
        """

        if hasattr(self, 'lambda_effective'):

            relative_lambda = self.lambda_effective - _colour_lower_lambda_limit
            relative_lambda = relative_lambda / _colour_lower_lambda_limit

            if verbose: print("relative_lambda = ", relative_lambda)

            self._plot_colour = colourmap(relative_lambda)

        else:
            warnings.warn("No self.lambda_effective set.")


class SpectrumClass():
    """

    """

    def __init__(self):
        pass


class SNClass():
    """docstring for SNClass."""

    def __init__(self):
        pass

##------------------------------------##
##                                    ##
##------------------------------------##

def load_filter(path, verbose = True):
    """
    Loads a filter response into FilterClass and returns it.

    Parameters
    ----------
    Returns
    -------
    """
    if check_file_path(os.path.abspath(path)):
        filter_object = FilterClass()
        filter_object.read_filter_file(os.path.abspath(path))

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


def _get_filter_directory(self):
    """
    Get the defaul path to the filter directory.

    Looks for the filter directory set as environment variable
    $PYCOCO_FILTER_DIR. if not found, returns default.

    returns: Absolute path in environment variable $PYCOCO_DATA_DIR, or
             default datalocation: '../testdata/'.
    """

    return os.environ.get('PYCOCO_FILTER_DIR', self._default_filter_dir_path)


def check_dir_path(path, verbose = False):
    """
    Parameters
    ----------
    Returns
    -------
    """
    try:
        if os.path.isdir(os.path.abspath(path)):
            if verbose: print("foo")
            return True
        else:
        #     if verbose: print("bar")
            warnings.warn(os.path.abspath(path) +
            " is not a valid directory. Returning 'False'.")
            return False
    except:
        raise PathError("The path '" + str(path) + "'is not a directory or doesn't exist.")
        return False


def check_file_path(path, verbose = True):
    """

    """
    try:
        if os.path.isfile(os.path.abspath(str(path))):
            if verbose: print("foo")
            return True
        else:
            warnings.warn(os.path.abspath(path) +
            " is not a valid file. Returning 'False'.")
            return False
    except:
        raise PathError("The data file '" + str(path) + "' doesn't exist or is a directory.")
        return False


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


def load_all_phot(path = _default_data_dir_path, format = 'ascii', verbose = True):
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
    phot_list = find_phot(path = path)

    if len(phot_list) > 0:
        # phot_table = Table()
        phot_table = ap.table.Table()

        for phot_file in phot_list:
            print(phot_file)
            print(phot_table.read(phot_file, format = format))

        return phot_table
    else:
        warning.warn("Couldn't find any photometry")


def find_phot(path = _default_data_dir_path, snname = False,
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

    if verbose:
        print("Found: ")
        print(ls)
        print("Matched:")
        print(phot_list)
    if len(phot_list) is 0:
        warnings.warn("No matches found.")
    return phot_list


def check_url_status(url):
    """
    Snippet from http://stackoverflow.com/questions/6471275 .

    Checks the status of a website - a status flag of < 400 means the site
    is up.

    """
    p = urlparse(url)
    conn = httplib.HTTPConnection(p.netloc)
    conn.request('HEAD', p.path)
    resp = conn.getresponse()

    return resp.status


def check_url(url):
    """
    Wrapper for check_url_status - considers the status, True if < 400.
    """
    return check_url_status(url) < 400


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



##----------------------------------------------------------------------------##
##  /CODE                                                                     ##
##----------------------------------------------------------------------------##
