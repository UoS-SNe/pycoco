import os
import shutil
import pycoco as pcc

dir_list = os.listdir('/Users/berto/data/CoreCollapse/phot/sn')

for sn in dir_list:
    if sn != '.DS_Store':
        photlist = os.listdir(os.path.join('/Users/berto/data/CoreCollapse/phot/sn', sn))
        for photfile in photlist:
            print photfile

                # shutil.move(os.path.join('/Users/berto/data/CoreCollapse/phot/sn', sn, photfile), os.path.join('/Users/berto/data/CoreCollapse/phot/sn', sn, newphotfile))
        P = pcc.PhotometryClass()
        P.load_phot_from_files(path = '/Users/berto/data/CoreCollapse/phot/sn/'+sn,
                               snname = sn)
        P.save(sn + '.dat', squash = True)
