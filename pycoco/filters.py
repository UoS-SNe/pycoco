"""

"""

import os

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d as interp1d

from .colours import *
from .utils import check_file_path


__all__ = ["FilterClass",]

class FilterClass():
    """Docstring for FilterClass"""

    def __init__(self, verbose = True):
        self._wavelength_units = u.Angstrom
        self._wavelength_units._format['latex'] = r'\rm{\AA}'
        self._frequency_units = u.Hertz
        # self.calculate_frequency()
        # self.calculate_effective_frequency()
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

            self.set_plot_colour(verbose = verbose)
            # self.
            self.calculate_effective_wavelength()
            self.calculate_edges()

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
        pass


    def calculate_frequency(self):
        nu = c/self.wavelength_u
        self.frequency_u = nu.to(self._frequency_units)
        self.frequency = self.frequency_u.value


    def calculate_effective_frequency(self):
        """

        """

        if hasattr(self, "frequency"):
            spline_rev = interp1d((np.cumsum(self.frequency*self.throughput)/np.sum(self.frequency*self.throughput)), self.frequency)
            nu_eff = spline_rev(0.5)

            self.nu_effective = nu_eff * self._frequency_units
        pass


    def calculate_edges_zero(self, verbose = False):
        """
        calculates the first and last wavelength that has non-zero and steps one
         away

        Parameters
        ----------

        Returns
        -------
        """

        ## calculates the first and last wavelength that has non-zero
        # w = np.where(self.throughput > 0)[0]
        # if verbose: print(w)
        # self._upper_edge = self.wavelength[w[-1]]
        # self._lower_edge = self.wavelength[w[0]]

        w = np.where(self.throughput > 0)[0]
        if verbose: print(w)
        if w[0] - 1 < 0:
            w_low = 0
        else:
            w_low =  w[0] - 1

        if w[-1] + 1 == len(self.throughput):
            w_high = w[-1]
        else:
            w_high = w[-1] + 1

        self._upper_edge = self.wavelength[w_high]
        self._lower_edge = self.wavelength[w_low]


    def calculate_edges(self, pc = 3., verbose = True):
        """
        calculates edges by defining the region that contains (100 - pc)% of the
        flux.

        Parameters
        ----------

        Returns
        -------
        """
        self._cumulative_throughput = np.cumsum(self.throughput)/np.sum(self.throughput)
        self._cumulative_throughput_spline = interp1d(self._cumulative_throughput, self.wavelength)

        self._upper_edge = self._cumulative_throughput_spline(1.0 - 0.5*(0.01*pc))
        self._lower_edge = self._cumulative_throughput_spline(0.0 + 0.5*(0.01*pc))

        pass


    def plot(self, xminorticks = 250, yminorticks = 0.1,
             show_lims = False, small = False, cumulative = False,
             *args, **kwargs):
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

            plot_label_string = r'$\textnormal{' + self.filter_name.replace('_', '\\_') + '}$'

            yminorLocator = MultipleLocator(yminorticks)
            xminorLocator = MultipleLocator(xminorticks)

            if not small:
                fig = plt.figure(figsize=[8, 4])
            else:
                fig = plt.figure(figsize=[4, 2])
                plt.rcParams['font.size'] = 10

            fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                                right = 0.99, hspace=0, wspace = 0)

            ax1 = fig.add_subplot(111)

            if cumulative:
                throughput = np.cumsum(self.throughput)/np.sum(self.throughput)
                yaxis_label_string = r'$\textnormal{Cumulative Throughput}$'

            else:
                throughput = self.throughput
                yaxis_label_string = r'$\textnormal{Fractional Throughput}$'


            if hasattr(self, "_plot_colour"):
                ax1.plot(self.wavelength, throughput, color = self._plot_colour,
                         lw = 2, label = plot_label_string)
            else:
                ax1.plot(self.wavelength, throughput, lw = 2, label = plot_label_string)

            if show_lims:
                try:
                    ax1.plot([self._upper_edge, self._upper_edge], [0,1] ,
                             lw = 1.5, alpha = 0.5, ls = ':',
                             color = hex['batman'], zorder = 0, )
                    ax1.plot([self._lower_edge, self._lower_edge], [0,1] ,
                             lw = 1.5, alpha = 0.5, ls = ':',
                             color = hex['batman'], zorder = 0, )
                except:
                    print("Failed")

            ax1.spines['top'].set_visible(True)

            ax1.set_xlabel(xaxis_label_string)
            ax1.set_ylabel(yaxis_label_string)

            ax1.yaxis.set_minor_locator(yminorLocator)
            ax1.xaxis.set_minor_locator(xminorLocator)

            ax1.legend(loc = 0)

            plt.show()
            pass
        else:
            warning.warn("Doesn't look like you have loaded a filter into the object")


    def resample_response(self, new_wavelength = False, k = 1,
                          *args, **kwargs):
        """
        Bit dodgy - spline has weird results for poorly sampled filters.
        Now the order is by default 1, seems to be less likely to introduce artifacts

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

            interp_func = InterpolatedUnivariateSpline(self.wavelength, self.throughput, k = k,
                                                       *args, **kwargs)
            self.throughput = interp_func(new_wavelength)
            self.wavelength = new_wavelength

            self.throughput[np.where(self.throughput < 0.0)] = 0.0
        else:
            warning.warn("Doesn't look like you have loaded a filter into the object")


    def calculate_plot_colour(self, colourmap_name = "plasma", verbose = False):
        """


        Parameters
        ----------

        Returns
        -------
        """

        if not hasattr(self, "_colourmap"):
            self._colourmap = plt.get_cmap(_colourmap_name)

        if hasattr(self, 'lambda_effective'):

            relative_lambda = self.lambda_effective - _colour_lower_lambda_limit
            relative_lambda = relative_lambda / _colour_lower_lambda_limit

            if verbose: print("relative_lambda = ", relative_lambda)

            self._plot_colour = self._colourmap(relative_lambda)

        else:
            warnings.warn("No self.lambda_effective set.")


    def set_plot_colour(self, colour = False, verbose = False):
        """


        Parameters
        ----------

        Returns
        -------
        """
        if colour:
            self._plot_colour = colour

        else:
            if verbose: print(hex[self.filter_name])
            try:
                self._plot_colour = hex[self.filter_name]
            except:
                if verbose: print("Nope")
                self.calculate_plot_colour(verbose = verbose)

        pass
