#!/usr/bin env python
'''
This is the utilities sub-module for the pycoco python tools.

author: Rob Firth, Southampton
date: 28-02-2016
'''

from __future__ import print_function

import os
import pycoco as pcc
from numpy import loadtxt, savetxt, array, array_equal

__all__ = ["relist"]


def _get_filters():
    """
    Parameters
    ----------

    Returns
    -------
    """
    filter_dir = pcc._get_filter_directory()
    file_list = os.listdir(filter_dir)

    for filter_file in file_list:
        if not os.path.isfile(os.path.join(filter_dir, filter_file)):
            file_list.remove(filter_file)
        elif filter_file[0] == ".":
            file_list.remove(filter_file)
        elif filter_file == "list.txt":
            file_list.remove(filter_file)

    return array(file_list)

def _check_filters():
    """
    Parameters
    ----------

    Returns
    -------
    """
    filter_dir = pcc._get_filter_directory()
    path = os.path.join(filter_dir, "list.txt")

    current_arr = loadtxt(path, dtype = str)

    filter_arr = _get_filters()

    return array_equal(current_arr, filter_arr)


def make_list_dot_txt():
    """
    Parameters
    ----------

    Returns
    -------
    """
    filter_dir = pcc._get_filter_directory()
    outpath = os.path.join(filter_dir, "list.txt")
    new_list = _get_filters()
    savetxt(outpath, new_list, fmt = "%s")
    pass

def relist(force = False):
    """
    Parameters
    ----------

    Returns
    -------
    """
    if force or not _check_filters:
        make_list_dot_txt()
    else:
        print("current list.txt is up to date. re run with force = True to force.")
    pass
