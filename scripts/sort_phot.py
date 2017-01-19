from __future__ import print_function ## Force python3-like printing

from matplotlib import pyplot as plt

import os
import numpy as np
from astropy.table import Table
import shutil

path = u'/Users/berto/data/CoreCollapse/phot/natasha_photometry_forCODE'
outpath = u'/Users/berto/data/CoreCollapse/phot/sn'
dirlist = os.listdir(path)

for dir in dirlist:
    print(os.path.abspath(os.path.join(path, dir)))
    if os.path.isdir(os.path.abspath(os.path.join(path, dir))):
        print(os.path.abspath(os.path.join(path, dir)))
        filelist = os.listdir(os.path.abspath(os.path.join(path, dir)))
        print(filelist)
        for file in filelist:
            if file != '.DS_Store':
                # print(os.path.abspath(os.path.join(path, dir, file)))
                print(dir)
                data = Table.read(os.path.abspath(os.path.join(path, dir, file)),
                                  format = 'ascii', names = ('MJD', 'flux', 'flux_err'))
                charar = np.chararray(len(data['MJD']), itemsize = len(dir))
                charar[:] = dir
                data['filter'] = charar
                print(charar)
                snname = file.split('.')[0]
                print(snname)
                if not os.path.isdir(os.path.abspath(os.path.join(outpath, snname))):
                    os.makedirs(os.path.abspath(os.path.join(outpath, snname)))

                data.write(os.path.abspath(os.path.join(outpath, snname, snname + '_' + dir+ '.dat')), format = 'ascii.commented_header')
