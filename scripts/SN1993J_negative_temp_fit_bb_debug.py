"""

"""
import os

import pycoco as pcc

snname = "SN1993J"

sn = pcc.classes.SNClass(snname)
sn.load_phot(path=os.path.join(pcc.defaults._default_data_dir_path, "lc/" + snname + ".dat"))
sn.load_list(os.path.join(pcc.defaults._default_list_dir_path, "SN1993J.list"))

sn.load_spec()
sn.get_lcfit(os.path.join(pcc.defaults._default_recon_dir_path, snname + ".dat"))
sn.check_overlaps()

spec_key = "1993J_16.0.txt"

print(spec_key)
S = sn.spec[spec_key]
S.get_specphot(sn.phot.data_filters, verbose=True)

new_spec = pcc.classes.SpectrumClass()

new_spec.load_table(pcc.kcorr.fit_bb(S, filter_dict=sn.phot.data_filters, return_table=True, verbose=True))

new_spec.plot()