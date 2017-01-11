import pycoco as pcc

P = pcc.PhotometryClass()
P.load_phot_from_files(path = '/Users/berto/data/CoreCollapse/phot/sn/2011dh', snname = '2011dh')

P.save('SN2011dh.dat', squash = True )
