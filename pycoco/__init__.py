'''
This is the module for the pycoco python tools.

author: Rob Firth, Southampton
date: 06-12-2016
'''

from __future__ import print_function ## Force python3-like printing

if __name__ is not '__main__':

    __name__ = 'pycoco'
    __version__ = 0.2

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
from astropy.table import Table, vstack
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d as interp1d

from .extinction import *
from .colours import *

warnings.simplefilter("error") ## Turn warnings into erros - good for debugging

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

__all__ = ["_default_data_dir_path",
           "_default_filter_dir_path",
           "_default_coco_dir_path",
           "_colourmap_name",
           "_spec_colourmap_name",
           "_colour_upper_lambda_limit",
           "_colour_lower_lambda_limit"]

## Important variables

_default_data_dir_path = os.path.abspath(os.path.join(__file__, os.pardir, os.pardir) + '/testdata/')
_default_filter_dir_path = os.path.abspath("/Users/berto/Code/CoCo/data/filters/")
_default_coco_dir_path = os.path.abspath("/Users/berto/Code/CoCo/")
# _colormap_name = 'jet'
_colourmap_name = 'rainbow'
_spec_colourmap_name = 'viridis'
_colourmap_name = 'plasma'

colourmap = plt.get_cmap(_colourmap_name)
spec_colourmap = plt.get_cmap(_spec_colourmap_name)

_colour_upper_lambda_limit = 11000 * u.angstrom
_colour_lower_lambda_limit = 3500 * u.angstrom

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


class TableReadError(ValueError):
    """
    Raise when something goes wrong with the table I/O
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
##  Classes                           ##
##------------------------------------##

class PhotometryClass():
    """
    Probably also overkill - but should be easier to store metadata etc. Hopefully
    flexible enough to just be a wrapper for AP tables of phot.

    Photometry stored in PhotometryClass.data should have a FilterClass method
    describing the observations stored in PhotometryClass.data_filters.

    ## NOTE should I use properties instead of get/set? http://www.python-course.eu/python3_properties.php
    looks like only python3?
    """

    def __init__(self, verbose = False):
        """

        """

        ## Initialise the class variables
        self._default_data_dir_path = os.path.join(_default_data_dir_path, 'lc/')
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
                 default datalocation: '../testdata/', with '/lc/' appended.
        """

        return os.path.join(os.path.abspath(os.environ.get('PYCOCO_DATA_DIR', os.path.join(self._default_data_dir_path, os.pardir))), "lc/")


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
            if verbose: print(data_dir_path, self._default_data_dir_path)
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


    def unpack(self, filter_file_type = '.dat', verbose = False):
        """
        If loading from preformatted file, then unpack the table into self.data
        OrderedDict and load FilterClass objects into self.data_filters OrderedDict

        Parameters
        ----------

        Returns
        -------

        """
        if hasattr(self, "phot"):
            filter_names = np.unique(self.phot["filter"])
            self.phot.add_index('filter', unique = True)


            for filter_name in filter_names:

                phot_table = self.phot.loc["filter", filter_name]
                filter_filename = filter_name + filter_file_type
                phot_table.meta = {"filter_filename": filter_filename}

                if verbose: print(phot_table)
                indices = phot_table.argsort("MJD")
                # for column_name in phot_table.colnames:
                #     phot_table[column_name] = phot_table[column_name][indices]
                sorted_phot_table = Table([phot_table[column_name][indices] for column_name in phot_table.colnames])
                filter_key = np.unique(phot_table["filter"])[0]

                if len(np.unique(phot_table["filter"])) > 1 or filter_key != filter_name:

                    raise FilterMismatchError("There is a more than one filterdata in here! or there is a mismatch with filename")

                path_to_filter = os.path.join(self.filter_directory, phot_table.meta['filter_filename'])

                self.data_filters[filter_key] = load_filter(path_to_filter, verbose = verbose)
                self.data[filter_name] = sorted_phot_table

        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")

        pass


    def load(self, path, names = ('MJD', 'flux', 'flux_err', 'filter'),
                  format = 'ascii.commented_header', verbose = True):
        """
        Loads a single photometry file.

        Parameters
        ----------
        Returns
        -------
        """
        StringWarning(path)
        try:
            phot_table = self._load_formatted_phot(path, names = names, format = format, verbose = verbose)
            self.phot = phot_table
            self.unpack()

            ## Sort the OrderedDict
            self._sort_phot()
        except:
            raise StandardError


    def _load_formatted_phot(self, path, format = "ascii", names = False,
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


    def load_phot_from_file(self, path, names = ('MJD', 'flux', 'flux_err', 'filter'),
                  format = 'ascii', verbose = True):
        """

        """
        StringWarning(path)
        try:
            phot_table = load_phot(path, names = names, format = format, verbose = verbose)
            self.data[np.unique(phot_table["filter"])[0]] = phot_table

            ## Sort the OrderedDict
            self._sort_phot()
        except:
            raise StandardError

        pass


    def load_phot_ap_tables(self):
        """

        """

        pass


    def load_phot_from_files(self, path = False, snname = False, prefix = 'SN',
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
            if not path:
                path = self._default_data_dir_path
            ## Find matching photometry
            phot_list = find_phot(path = path, snname = snname, prefix = prefix,
                              file_type = file_type, verbose = verbose)

            full_phot_table = Table()

            ## Loop over files (shouldn't be that many really)
            if len(phot_list) > 0:

                for phot_file in phot_list:

                    if verbose: print(phot_file)
                    phot_table = Table.read(phot_file, names = names, format = format)

                    ## NOTE astropy vstack does not support mixin columns http://docs.astropy.org/en/stable/table/mixin_columns.html
                    # This means I might have problems joining the tables together if I don't add together as I go along.

                    full_phot_table = vstack([full_phot_table, phot_table])

                    filter_string = get_filter_from_filename(phot_file, snname, file_type)
                    phot_table.meta = {"filename" : phot_file,
                                       "filter" : filter_string,
                                       "filter_filename": filter_string + filter_file_type}

                    ## Sort out units
                    phot_table.sort("MJD")
                    phot_table["t"] = Time(phot_table["MJD"], format = 'mjd')

                    phot_table["MJD"].unit = u.day
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


                ## NOTE doing it this way because vstack doesn't like mixin columns (see above comment)
                full_phot_table.sort("MJD")
                # full_phot_table["t"] = Time(full_phot_table["MJD"], format = 'mjd')
                full_phot_table["MJD"].unit = u.day

                full_phot_table["flux"].unit = u.cgs.erg / u.si.angstrom / u.si.cm ** 2 / u.si.s
                full_phot_table["flux_err"].unit =  full_phot_table["flux"].unit

                self.phot = full_phot_table

                ## Sort the OrderedDict
                self._sort_phot()
            else:
                warning.warn("Couldn't find any photometry")
        else:
            warnings.warn("Provide a SN name")

        pass


    def _sort_phot(self):
        """
        resorts the photometry according to effective wavelength of the filter.

        Parameters
        ----------

        Returns
        -------

        """
        if hasattr(self, "data") and hasattr(self, "data_filters"):
            ## This looks fugly.
            newkeys = np.array(self.data_filters.keys())[np.argsort([self.data_filters[i].lambda_effective.value for i in self.data_filters])]

            sorted_data = OrderedDict()
            sorted_data_filters = OrderedDict()

            for newkey in newkeys:
                sorted_data[newkey] = self.data[newkey]
                sorted_data_filters[newkey] = self.data_filters[newkey]

            self.data = sorted_data
            self.data_filters = sorted_data_filters

        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")
        pass


    def _combine_phot(self, verbose = True):
        """

        """

        if hasattr(self, "data"):
            if verbose: print(self.data.keys())

            for i, phot_filter in enumerate(self.data.keys()):

                if verbose: print(i, phot_filter)

                if i == 0:

                    full_phot = self.data[phot_filter]

                else:

                    full_phot = vstack([full_phot, self.data[phot_filter]])

                    pass

            self.data['full'] = full_phot

        else:
            warnings.warn("Cant find self.data")

        pass


    def save(self, filename, path = False,
             squash = False, verbose = True, *args, **kwargs):
        """
        Output the photometry loaded into the SNClass via self.load_phot* into a format
        and location recognised by CoCo.

        Parameters
        ----------
        Returns
        -------
        """

        if hasattr(self, "data"):
            if verbose: print("has data")
            if not path:
                if verbose: print("No directory specified, assuming " + self._default_data_dir_path)
                path = self._default_data_dir_path
            else:
                StringWarning(path)

            outpath = os.path.join(path, filename)

            check_dir_path(path)

            if os.path.isfile(outpath):
                warnings.warn("Found existing file matching " + path + ". Run with squash = True to overwrite")
                if squash:
                    print("Overwriting " + outpath)
                    self._phot_format_for_save().write(outpath, format = "ascii.fast_commented_header")


            else:
                    print("Writing " + outpath)
                    self._phot_format_for_save().write(outpath, format = "ascii")

        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")
        pass


    def _phot_format_for_save(self):
        """
        This is hacky - clear it up!

        Parameters
        ----------
        Returns
        -------
        """

        save_table = self.phot
        save_table['MJD'].format = "5.5f"
        save_table['flux'].format = "5.5e"
        save_table['flux_err'].format = "5.5e"

        return save_table


    def plot(self, legend = True, xminorticks = 5,
             verbose = False, *args, **kwargs):
        """
        Plots phot.

        Parameters
        ----------

        Returns
        -------
        """

        if hasattr(self, "data"):

            setup_plot_defaults()

            fig = plt.figure(figsize=[8, 4])
            fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                                right = 0.99, hspace=0, wspace = 0)

            ax1 = fig.add_subplot(111)

            for i, filter_key in enumerate(self.data_filters):
                if verbose: print(i, self.data[filter_key].__dict__)
                plot_label_string = r'$\rm{' + self.data_filters[filter_key].filter_name.replace('_', '\\_') + '}$'
                if filter_key in hex.keys():
                    self.data_filters[filter_key]._plot_colour = hex[filter_key]

                ax1.errorbar(self.data[filter_key]['MJD'], self.data[filter_key]['flux'],
                             yerr = self.data[filter_key]['flux_err'],
                             capsize = 0, fmt = 'o', color = self.data_filters[filter_key]._plot_colour,
                             label = plot_label_string, ecolor = hex['batman'],
                             *args, **kwargs)

            if legend:

                plot_legend = ax1.legend(loc = [1.,0.0], scatterpoints = 1,
                                      numpoints = 1, frameon = False, fontsize = 12)

            ## Use ap table groups instead? - can't; no support for mixin columns.
            ax1.set_ylim(np.nanmin(self.phot['flux']), np.nanmax(self.phot['flux']))

            ## Label the axes
            xaxis_label_string = r'$\textnormal{Time, MJD (days)}$'
            yaxis_label_string = r'$\textnormal{Flux, erg s}^{-1}\textnormal{\AA}^{-1}\textnormal{cm}^{-2}$'

            ax1.set_xlabel(xaxis_label_string)
            ax1.set_ylabel(yaxis_label_string)

            xminorLocator = MultipleLocator(xminorticks)
            ax1.xaxis.set_minor_locator(xminorLocator)

            plt.show()
        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")
        pass


    def plot_filters(self, xminorticks = 250, yminorticks = 0.1,
                     legend = True, use_cmap = False, verbose = False):
        """
        Plots filters.

        Parameters
        ----------

        Returns
        -------
        """
        if hasattr(self, "data_filters"):

            setup_plot_defaults()
            xaxis_label_string = r'$\textnormal{Wavelength, (\AA)}$'
            yaxis_label_string = r'$\textnormal{Fractional Throughput}$'
            yminorLocator = MultipleLocator(yminorticks)
            xminorLocator = MultipleLocator(xminorticks)

            fig = plt.figure(figsize=[8, 4])
            fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                                right = 0.99, hspace=0, wspace = 0)
            ax1 = fig.add_subplot(111)

            ## Plot the throughput for each filter
            for i, filter_key in enumerate(self.data_filters):
                if verbose: print(i, self.data_filters[filter_key].__dict__)
                plot_label_string = r'$\rm{' + self.data_filters[filter_key].filter_name.replace('_', '\\_') + '}$'
                if hasattr(self.data_filters[filter_key], "_plot_colour") and use_cmap:
                    ax1.plot((self.data_filters[filter_key].wavelength_u).to(u.angstrom),
                             self.data_filters[filter_key].throughput,
                             color = self.data_filters[filter_key]._plot_colour,
                             lw = 2, label = plot_label_string)
                else:
                    ax1.plot((self.data_filters[filter_key].wavelength_u).to(u.angstrom),
                             self.data_filters[filter_key].throughput,
                             lw = 2, label = plot_label_string)
            # if hasattr(self, "_plot_colour"):
            #     ax1.plot(self.wavelength, self.throughput, color = self._plot_colour,
            #              lw = 2, label = plot_label_string)
            # else:
            #     ax1.plot(self.wavelength, self.throughput, lw = 2, label = plot_label_string)

            ax1.set_xlabel(xaxis_label_string)
            ax1.set_ylabel(yaxis_label_string)

            ax1.yaxis.set_minor_locator(yminorLocator)
            ax1.xaxis.set_minor_locator(xminorLocator)

            if legend:
                plot_legend = ax1.legend(loc = [1.,0.0], scatterpoints = 1,
                                      numpoints = 1, frameon = False, fontsize = 12)

            plt.show()
        else:
            warnings.warn("Doesn't seem to be any filters here (empty self.filter_data)")
        pass


class FilterClass():
    """Docstring for FilterClass"""

    def __init__(self, verbose = True):
        self._wavelength_units = u.Angstrom
        self._wavelength_units._format['latex'] = r'\rm{\AA}'
        pass


    def read_filter_file(self, path, wavelength_units = u.angstrom, verbose = False):
        """
        Assumes Response function is fractional rather than %.
        """
        if check_file_path(os.path.abspath(path), verbose = verbose):
            self.wavelength, self.throughput = np.loadtxt(path).T
            self.wavelength_u = self.wavelength * wavelength_units
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


    def calculate_plot_colour(self, colourmap = colourmap, verbose = False):
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
    Class for handling Spectra.
    """

    def __init__(self):
        """

        """

        ## Initialise the class variables
        self._default_data_dir_path = os.path.abspath(os.path.join(_default_data_dir_path, "spec/"))
        self._default_list_dir_path = self._default_data_dir_path

        ## Initialise using class methods
        self.set_data_directory(self._get_data_directory())

        pass


    def _get_data_directory(self):
        """
        Get the default path to the data directory.

        Looks for the data data directory set as environment variable
        $PYCOCO_DATA_DIR. if not found, returns default.

        returns: Absolute path in environment variable $PYCOCO_DATA_DIR, or
                 default datalocation: '../testdata/', with '/spec/' appended.
        """

        return os.path.join(os.path.abspath(os.environ.get('PYCOCO_DATA_DIR', os.path.join(self._default_data_dir_path, os.pardir))), "spec/")


    def set_data_directory(self, data_dir_path = '', verbose = False):
        """
        Set a new data directory path.

        Enables the data directory to be changed by the user.

        """
        try:
            if verbose: print(data_dir_path, self._default_data_dir_path)
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


    def load(self, filename, directory = False, fmt = "ascii",
             names = ("wavelength", "flux"), wavelength_u = u.angstrom,
             flux_u = u.cgs.erg / u.si.cm ** 2 / u.si.s, verbose = True):
        """
        Parameters
        ----------

        Returns
        -------
        """


        StringWarning(filename)

        if not directory:
            path = os.path.abspath(os.path.join(self.data_directory, filename))
            if verbose: print("You didn't supply a directory, so using self.data_directory")
        else:
            StringWarning(directory)
            check_dir_path(directory)

            path = os.path.abspath(os.path.join(directory, filename))

        if os.path.isfile(path):

            ## Some might have three columns, deal with laters - this is untidy
            try:
                spec_table = Table.read(path, format = fmt, names = names)

            except:
                names = names + ("flux_err",)
                spec_table = Table.read(path, format = fmt, names = names)

            if verbose:print("Reading " + path)

            spec_table.meta = {"filename": path}

            spec_table['wavelength'].unit = wavelength_u
            spec_table['flux'].unit = flux_u

            if "flux_err" in spec_table.colnames:
                spec_table["flux_err"].unit = flux_u

            self.data = spec_table
            self.wavelength = spec_table["wavelength"]
            self.flux = spec_table["flux"]

        else:
            warnings.warn(path + " is not a valid file path")


    def plot(self, xminorticks = 250, legend = True,
             verbose = True, compare_red = True,
             *args, **kwargs):
        """
        Plots spec.

        Parameters
        ----------

        Returns
        -------
        """

        if hasattr(self, "data"):

            setup_plot_defaults()

            fig = plt.figure(figsize=[8, 4])
            fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                                right = 0.99, hspace=0, wspace = 0)

            ax1 = fig.add_subplot(111)


            if verbose: print(self.data.__dict__)
            plot_label_string = r'$\rm{' + self.data.meta["filename"] + '}$'

            ax1.plot(self.data['wavelength'], self.flux, lw = 2,
                         label = plot_label_string, color = 'Red',
                         *args, **kwargs)

            maxplotydata = np.nanmax(self.flux)
            minplotydata = np.nanmin(self.flux)

            if hasattr(self, 'flux_dered') and compare_red:
                ax1.plot(self.data['wavelength'], self.data['flux_dered'], lw = 2,
                             label = plot_label_string, color = 'Blue',
                             *args, **kwargs)
                maxplotydata = np.nanmax(np.append(maxplotydata, np.nanmax(self.data['flux_dered'])))
                minplotydata = np.nanmin(np.append(minplotydata, np.nanmin(self.data['flux_dered'])))
            if legend:

                plot_legend = ax1.legend(loc = [1.,0.0], scatterpoints = 1,
                                      numpoints = 1, frameon = False, fontsize = 12)

            ax1.set_ylim(minplotydata*0.98, maxplotydata*1.02)

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


    def get_MJD_obs(self, list_filename, list_dir = False, verbose = True):
        """
        Retrieve the MJD of the observation from a '.list' file.

        Parameters
        ----------

        Returns
        -------
        """

        try:

            if not list_dir:
                list_dir = self._default_list_dir_path

            check_dir_path(list_dir)
            list_path = os.path.abspath(os.path.join(list_dir, list_filename))
            check_file_path(list_path)
        except:

            raise PathError("The data file '" + str(path) + "' doesn't exist or is a directory.")

            return False

        data = np.genfromtxt(list_path, dtype = np.str)
        short_filenames = [f.split('/')[-1] for f in  data.T[0]]
        filename = self.data.meta['filename'].split('/')[-1]
        print(filename)
        if verbose: print(data.T[0])

        if filename in short_filenames:
            print("Foo")

        # pass
        return data


    def set_MJD_obs(self, mjd):
        """
        Log MJD of the observation.

        Parameters
        ----------

        Returns
        -------
        """
        self.mjd_obs = mjd

        pass


    def set_EBV(self, EBV):
        """
        Parameters
        ----------

        Returns
        -------
        """
        self.EBV = EBV


    def deredden(self, verbose = True):
        """
        Parameters
        ----------

        Returns
        -------
        """

        if hasattr(self, "EBV") and hasattr(self, "data"):
            if verbose: print("Foo")

            self.flux_dered = unred(self.wavelength, self.flux, EBV_MW = self.EBV)
            self.data["flux_dered"] = self.flux_dered

        else:
            warnings.warn("No extinction value set")
        pass


    def use_flux_dered(self):
        """
        Parameters
        ----------

        Returns
        -------
        """

        if hasattr(self, "data"):
            self.flux_red = self.flux
            self.flux = self.data['flux_dered']
        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")
        pass


    def _spec_format_for_save(self):
        """
        Parameters
        ----------
        Returns
        -------
        """

        save_table = Table()

        save_table['wavelength'] = self.wavelength
        save_table['flux'] = self.flux

        save_table['wavelength'].format = "5.5f"
        save_table['flux'].format = "5.5e"

        return save_table


    def save(self, filename, path = False,
             squash = False, verbose = True, *args, **kwargs):
        """
        Output the spectrum loaded into the Class via self.load into a format
        and location recognised by CoCo.

        Parameters
        ----------
        Returns
        -------
        """

        if hasattr(self, "data"):
            if verbose: print("has data")
            if not path:
                if verbose: print("No directory specified, assuming " + self._default_data_dir_path)
                path = self._default_data_dir_path
            else:
                StringWarning(path)

            outpath = os.path.join(path, filename)

            check_dir_path(path)

            if os.path.isfile(outpath):
                warnings.warn("Found existing file matching " + path + ". Run with squash = True to overwrite")
                if squash:
                    print("Overwriting " + outpath)
                    self._spec_format_for_save().write(outpath, format = "ascii.fast_commented_header")


            else:
                    print("Writing " + outpath)
                    self._spec_format_for_save().write(outpath, format = "ascii")

        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")
        pass


class SNClass():
    """docstring for SNClass."""

    def __init__(self, snname):
        """
        Parameters
        ----------

        Returns
        -------
        """
        ## Initialise
        self.spec = OrderedDict()
        # self.spec = SpectrumClass()
        self.phot = PhotometryClass()

        self.coco_directory = self._get_coco_directory()

        self.name = snname
        pass


    def _get_coco_directory(self):
        """
        Get the default path to the data directory.

        Looks for the CoCo home directory set as environment variable
        $COCO_ROOT_DIR. if not found, returns default.

        returns: Absolute path in environment variable $COCO_ROOT_DIR, or
                 default CoCo location: '~/Code/CoCo/', with appended.
        """

        return os.path.abspath(os.environ.get('COCO_ROOT_DIR', os.path.abspath(_default_coco_dir_path)))


    def load_phot(self, snname = False, path = False, file_type = '.dat',
                  verbose = True):
        """
        Parameters
        ----------

        Returns
        -------
        """

        if not snname:
            snname = self.name
        if not path:
            path = os.path.abspath(os.path.join(self.phot._default_data_dir_path, snname + file_type))
        if verbose: print(path)
        self.phot.load(path, verbose = verbose)

        pass


    def load_list(self, path, verbose = True):
        """
        Parameters
        ----------
        Returns
        -------
        """
        listdata = read_list_file(path, verbose = verbose)
        listdata.sort('mjd_obs')
        self.list  = listdata


    def load_spec(self, snname = False, spec_dir_path = False, verbose = False):
        """
        Parameters
        ----------

        Returns
        -------
        """


        # if not snname:
        #     snname = self.name
        #
        # if not spec_dir_path:
        #     spec_dir_path = os.path.abspath(os.path.join(self._default_spec_data_dir_path, snname))
        #
        # if verbose: print("Loading spectra from: ", spec_dir_path)

        # spec_dir_path =


        if hasattr(self, 'coco_directory') and hasattr(self, 'list'):
            for i, path in enumerate(self.list['spec_path']):
                spec_fullpath = os.path.abspath(os.path.join(self.coco_directory, path))
                spec_filename = path.split('/')[-1]
                spec_dir_path = spec_fullpath.replace(spec_filename, '')
                if verbose: print(spec_dir_path, spec_filename)

                self.spec[spec_filename] = SpectrumClass()
                self.spec[spec_filename].load(spec_filename, directory = spec_dir_path, verbose = verbose)
                self.spec[spec_filename].set_MJD_obs(self.list['mjd_obs'][i])
                # self.spec[spec_filename].data.add_index('wavelength')

        else:
            warnings.warn("no coco or no listfile")
        pass


    def plot_lc(self, filters = False, legend = True, xminorticks = 5, mark_spectra = True,
                fit = True,
                verbose = False, *args, **kwargs):
        """
        Parameters
        ----------

        Returns
        -------
        """
        if hasattr(self.phot, "data"):

            if not filters:
                filters = self.phot.data_filters

            setup_plot_defaults()

            fig = plt.figure(figsize=[8, 4])
            fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                                right = 0.99, hspace=0, wspace = 0)

            ax1 = fig.add_subplot(111)

            for i, filter_key in enumerate(filters):
                if filter_key in self.phot.data:
                    if verbose: print(i, self.phot.data[filter_key].__dict__)
                    plot_label_string = r'$\rm{' + self.phot.data_filters[filter_key].filter_name.replace('_', '\\_') + '}$'
                    if filter_key in hex.keys():
                        self.phot.data_filters[filter_key]._plot_colour = hex[filter_key]

                    ax1.errorbar(self.phot.data[filter_key]['MJD'], self.phot.data[filter_key]['flux'],
                                 yerr = self.phot.data[filter_key]['flux_err'],
                                 capsize = 0, fmt = 'o', color = self.phot.data_filters[filter_key]._plot_colour,
                                 label = plot_label_string, ecolor = hex['batman'],
                                 *args, **kwargs)

                    if fit and hasattr(self, 'fit'):
                        ax1.fill_between(self.fit.data[filter_key]['MJD'], self.fit.data[filter_key]['flux_upper'], self.fit.data[filter_key]['flux_lower'],
                                         color = self.phot.data_filters[filter_key]._plot_colour,
                                         alpha = 0.8, zorder = 0,
                                         *args, **kwargs)
                else:
                    warnings.warn("Filter '" + filter_key + "' not found")

            if mark_spectra:
                for spec_key in self.spec:
                    plt.plot([self.spec[spec_key].mjd_obs, self.spec[spec_key].mjd_obs],
                             [0.0, np.nanmax(self.phot.phot['flux'])*1.5],
                             ls = ':', color = hex['batman'], zorder = 0)

            if legend:

                plot_legend = ax1.legend(loc = [1.,0.0], scatterpoints = 1,
                                      numpoints = 1, frameon = False, fontsize = 12)

            ## Use ap table groups instead? - can't; no support for mixin columns.
            ax1.set_ylim(np.nanmin(self.phot.phot['flux']), np.nanmax(self.phot.phot['flux']))

            ## Label the axes
            xaxis_label_string = r'$\textnormal{Time, MJD (days)}$'
            yaxis_label_string = r'$\textnormal{Flux, erg s}^{-1}\textnormal{\AA}^{-1}\textnormal{cm}^{-2}$'

            ax1.set_xlabel(xaxis_label_string)
            ax1.set_ylabel(yaxis_label_string)

            xminorLocator = MultipleLocator(xminorticks)
            ax1.xaxis.set_minor_locator(xminorLocator)

            plt.show()
        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")
        pass


    def plot_spec(self, xminorticks = 250, legend = True,
                  verbose = False, add_mjd = True,
                  *args, **kwargs):
        """
        Parameters
        ----------

        Returns
        -------
        """
        if hasattr(self, "spec"):

            setup_plot_defaults()

            fig = plt.figure(figsize=[8, 10])
            fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                                right = 0.99, hspace=0, wspace = 0)

            ax1 = fig.add_subplot(111)

            cmap_indices = np.linspace(0,1, len(self.spec))

            j = 0
            for i, spec_key in enumerate(self.spec):
                # if verbose: print(self.spec[spec_key].data.__dict__)

                plot_label_string = r'$\rm{' + self.spec[spec_key].data.meta["filename"].split('/')[-1].replace('_', '\_') + '}$'


                v_eff = 5436.87 ##Angstrom
                w = np.logical_and(self.spec[spec_key].data['wavelength'] > (v_eff-100.),self.spec[spec_key].data['wavelength'] < v_eff+100.)

                if verbose: print(i, len(w[np.where(w == True)]), spec_key, len(self.spec[spec_key].data['wavelength']), len(self.spec[spec_key].data['flux']), len(self.spec[spec_key].flux))
                if len(w[np.where(w == True)]) > 0:
                    if verbose: print(len(w), 'Foo')
                    flux_norm = self.spec[spec_key].flux / np.nanmean(self.spec[spec_key].flux[w])

                    ax1.plot(self.spec[spec_key].data['wavelength'], flux_norm - 0.5*j, lw = 2,
                                 label = plot_label_string, color = spec_colourmap(cmap_indices[i]),
                                 *args, **kwargs)

                    maxspecxdata = np.nanmax(self.spec[spec_key].data['wavelength'])
                    minspecxdata = np.nanmin(self.spec[spec_key].data['wavelength'])

                    yatmaxspecxdata = np.nanmean((flux_norm - 0.5*j)[-10:-1])
                    yatminspecxdata = np.nanmean((flux_norm - 0.5*j)[0:10])

                    if i == 0:
                        maxplotydata = np.nanmax(flux_norm - 0.5*j)
                        minplotydata = np.nanmin(flux_norm - 0.5*j)

                        maxplotxdata = maxspecxdata
                        minplotxdata = np.nanmin(self.spec[spec_key].data['wavelength'])
                    else:
                        maxplotydata = np.nanmax(np.append(maxplotydata, flux_norm - 0.5*j))
                        minplotydata = np.nanmin(np.append(minplotydata, flux_norm - 0.5*j))

                        maxplotxdata = np.nanmax(np.append(maxplotxdata, np.nanmax(self.spec[spec_key].data['wavelength'])))
                        minplotxdata = np.nanmin(np.append(minplotxdata, np.nanmin(self.spec[spec_key].data['wavelength'])))
                    if add_mjd:
                        # ax1.plot([maxspecxdata, 11000],[1 - 0.5*j, 1 - 0.5*j], ls = '--', color = hex['batman'])
                        # ax1.plot([maxspecxdata, 11000],[yatmaxspecxdata, yatmaxspecxdata], ls = '--', color = hex['batman'])
                        ax1.plot([1000, minspecxdata],[yatminspecxdata, yatminspecxdata], ls = '--', color = hex['batman'])
                        txt = ax1.text(1000, yatminspecxdata, r'$' + str(self.spec[spec_key].mjd_obs) + '$',
                                       horizontalalignment = 'right', verticalalignment = 'center')
                        # ax1.text(1000, 1 - 0.5*j, r'$' + str(self.spec[spec_key].mjd_obs) + '$', horizontalalignment = 'right')
                    j = j + 1
                else:
                    if verbose: print("Not enough data to normalise")
            if legend:

                plot_legend = ax1.legend(loc = [1.,0.0], scatterpoints = 1,
                                      numpoints = 1, frameon = False, fontsize = 12)

            ax1.set_ylim(minplotydata - 0.5, maxplotydata + 0.5)
            if verbose: print(minplotydata, maxplotydata)
            ## Label the axes
            xaxis_label_string = r'$\textnormal{Wavelength (\AA)}$'
            yaxis_label_string = r'$\textnormal{Flux, Arbitrary}$'

            ax1.set_xlabel(xaxis_label_string)
            ax1.set_ylabel(yaxis_label_string)

            ax1.set_yticklabels('')

            xminorLocator = MultipleLocator(xminorticks)
            ax1.xaxis.set_minor_locator(xminorLocator)

            plt.show()
        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")
        pass


    def get_fit(self, path):
        StringWarning(path)
        self.fit = LCfitClass()
        self.fit.load_formatted_phot(path)
        self.fit.unpack()
        self.fit._sort_phot()
        pass


class LCfitClass():
    """
    Small class to hold the output from CoCo LCfit
    """

    def __init__(self):

        ## Initialise the class variables
        self._default_recon_dir_path = os.path.join(_default_coco_dir_path, "recon/")
        self._default_filter_dir_path = _default_filter_dir_path

        ## Initialise using class methods
        self.set_recon_directory(self._get_recon_directory())
        self.set_filter_directory(self._get_filter_directory())

        ## Initialise some other stuff
        self.data = OrderedDict()
        self.data_filters = OrderedDict()

        pass


    def _get_recon_directory(self):
        """
        Get the default path to the data directory.

        Looks for the CoCo home directory set as environment variable
        $COCO_ROOT_DIR. if not found, returns default.

        returns: Absolute path in environment variable $COCO_ROOT_DIR, or
                 default CoCo location: '~/Code/CoCo/', with 'recon/' appended.
        """

        return os.path.join(os.path.abspath(os.environ.get('COCO_ROOT_DIR', os.path.join(self._default_recon_dir_path, os.path.pardir))), "recon/")


    def set_recon_directory(self, recon_dir_path = '', verbose = False):
        """
        Set a new recon directory path.

        Enables the recon directory to be changed by the user.

        """
        try:
            if verbose: print(recon_dir_path, self._default_recon_dir_path)
            if os.path.isdir(os.path.abspath(recon_dir_path)):
                self.recon_directory = os.path.abspath(recon_dir_path)
                pass
            else:
                warnings.warn(os.path.abspath(recon_dir_path) +
                " is not a valid directory. Restoring default path: " +
                self._default_recon_dir_path, UserWarning)
                self.recon_directory = self._default_recon_dir_path

                if not os.path.isdir(self.recon_directory):
                    if verbose: print(os.path.isdir(self.recon_directory))
                    raise PathError("The default recon directory '" + self.recon_directory
                     + "' doesn't exist. Or isn't a directory. Or can't be located.")
                else:
                    pass
        except:
            if verbose: print("foo")
            raise PathError("The default recon directory '" + self._default_recon_dir_path
             + "' doesn't exist. Or isn't a directory. Or can't be located. Have"
             + " you messed with _default_recon_dir_path?")
            pass


    def load_formatted_phot(self, path, names = ('MJD', 'flux', 'flux_err', 'filter'),
                  format = 'ascii', verbose = True):
        """

        """
        StringWarning(path)

        try:
            phot_table = load_formatted_phot(path, format = format, names = names,
                                             verbose = verbose)
            self.phot = phot_table

            self.phot['flux_upper'] = phot_table['flux'] + phot_table['flux_err']
            self.phot['flux_lower'] = phot_table['flux'] - phot_table['flux_err']

        except:
            raise StandardError

        pass


    def unpack(self, filter_file_type = '.dat', verbose = False):
        """
        If loading from preformatted file, then unpack the table into self.data
        OrderedDict and load FilterClass objects into self.data_filters OrderedDict

        Parameters
        ----------

        Returns
        -------

        """

        if hasattr(self, "phot"):
            filter_names = np.unique(self.phot["filter"])
            self.phot.add_index('filter', unique = True)


            for filter_name in filter_names:
                phot_table = self.phot.loc["filter", filter_name]
                filter_filename = filter_name + filter_file_type
                phot_table.meta = {"filter_filename": filter_filename}

                if verbose: print(phot_table)
                indices = phot_table.argsort("MJD")
                # for column_name in phot_table.colnames:
                #     phot_table[column_name] = phot_table[column_name][indices]
                sorted_phot_table = Table([phot_table[column_name][indices] for column_name in phot_table.colnames])
                filter_key = np.unique(phot_table["filter"])[0]

                if len(np.unique(phot_table["filter"])) > 1 or filter_key != filter_name:
                    raise FilterMismatchError("There is a more than one filterdata in here! or there is a mismatch with filename")
                path_to_filter = os.path.join(self.filter_directory, phot_table.meta['filter_filename'])

                self.data_filters[filter_key] = load_filter(path_to_filter)
                self.data[filter_name] = sorted_phot_table
        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")

        pass


    def _get_filter_directory(self):
        """
        Get the default path to the filter directory.

        Looks for the filter data directory set as environment variable
        $PYCOCO_FILTER_DIR. if not found, returns default.

        returns: Absolute path in environment variable $PYCOCO_FILTER_DIR, or
                 default datalocation: '/Users/berto/Code/CoCo/data/filters/'.
        """
        return os.path.abspath(os.environ.get('PYCOCO_FILTER_DIR', self._default_filter_dir_path))


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


    def plot(self, legend = True, xminorticks = 5,
             verbose = False, *args, **kwargs):
        """
        Plots phot.

        Parameters
        ----------

        Returns
        -------
        """

        if hasattr(self, "data"):

            setup_plot_defaults()

            fig = plt.figure(figsize=[8, 4])
            fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                                right = 0.99, hspace=0, wspace = 0)

            ax1 = fig.add_subplot(111)

            for i, filter_key in enumerate(self.data_filters):
                if verbose: print(i, self.data[filter_key].__dict__)
                plot_label_string = r'$\rm{' + self.data_filters[filter_key].filter_name.replace('_', '\\_') + '}$'

                # ax1.errorbar(self.data[filter_key]['MJD'], self.data[filter_key]['flux'],
                #              yerr = self.data[filter_key]['flux_err'],
                #              capsize = 0, fmt = 'o',
                #              label = plot_label_string,
                #              *args, **kwargs)

                # ## Best Fit
                # ax1.plot(self.data[filter_key]['MJD'], self.data[filter_key]['flux'],
                #          lw = 2, label = plot_label_string,
                #           *args, **kwargs)

                ## With error
                ax1.fill_between(self.data[filter_key]['MJD'], self.data[filter_key]['flux_upper'], self.data[filter_key]['flux_lower'],
                                 label = plot_label_string, color = self.data_filters[filter_key]._plot_colour,
                                 alpha = 0.8,
                                 *args, **kwargs)
            if legend:

                plot_legend = ax1.legend(loc = [1.,0.0], scatterpoints = 1,
                                      numpoints = 1, frameon = False, fontsize = 12)

            ## Use ap table groups instead? - can't; no support for mixin columns.
            ax1.set_ylim(np.nanmin(self.phot['flux']), np.nanmax(self.phot['flux']))

            ## Label the axes
            xaxis_label_string = r'$\textnormal{Time, MJD (days)}$'
            yaxis_label_string = r'$\textnormal{Flux, erg s}^{-1}\textnormal{\AA}^{-1}\textnormal{cm}^{-2}$'

            ax1.set_xlabel(xaxis_label_string)
            ax1.set_ylabel(yaxis_label_string)

            xminorLocator = MultipleLocator(xminorticks)
            ax1.xaxis.set_minor_locator(xminorLocator)

            plt.show()
        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")
        pass


    def _sort_phot(self):
        """
        resorts the photometry according to effective wavelength of the filter.

        Parameters
        ----------

        Returns
        -------

        """
        if hasattr(self, "data") and hasattr(self, "data_filters"):
            ## This looks fugly.
            newkeys = np.array(self.data_filters.keys())[np.argsort([self.data_filters[i].lambda_effective.value for i in self.data_filters])]

            sorted_data = OrderedDict()
            sorted_data_filters = OrderedDict()

            for newkey in newkeys:
                sorted_data[newkey] = self.data[newkey]
                sorted_data_filters[newkey] = self.data_filters[newkey]

            self.data = sorted_data
            self.data_filters = sorted_data_filters

        else:
            warnings.warn("Doesn't seem to be any data here (empty self.data)")
        pass


##------------------------------------##
##                                    ##
##------------------------------------##

def load_filter(path, verbose = False):
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


def check_file_path(path, verbose = False):
    """

    """
    try:
        if os.path.isfile(os.path.abspath(str(path))):
            if verbose: print("bar")
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

    if os.path.join(path, snname + file_type) in phot_list:
        phot_list.remove(os.path.join(path,snname + file_type))
        warnings.warn("Found " + os.path.join(path,snname + file_type) + " - you could just read that in.")

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



##------------------------------------##
## CoCo                               ##
##------------------------------------##

def run_LCfitClass():
    """
    Parameters
    ----------
    Returns
    -------
    """
    check_file_path(path)
    if verbose: print("Running CoCo lcfit on " + path)
    subprocess.call(["./lcfit", path])

    pass

##----------------------------------------------------------------------------##
##  /CODE                                                                     ##
##----------------------------------------------------------------------------##
