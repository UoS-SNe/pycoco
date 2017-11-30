"""
incorporated into test_pycoco.py as
test_SpectrumClass_get_specphot_works_with_one_overlapping_filter
"""

import os
import pycocosn as pcc

snname = "SN1993J"
all_listfiles = pcc.coco.get_all_spec_lists()
all_listfiles
listfile = '/Users/berto/Code/CoCo/lists/SN1993J.list'
list_table = pcc.utils.read_list_file(listfile)
snname = list_table["snname"][0]
snname
sn = pcc.classes.SNClass(snname)
sn.load_phot(path=os.path.join(pcc.defaults._default_data_dir_path, "lc/" + snname + ".dat"))
sn.load_list(listfile)
sn.load_spec()
sn.check_overlaps()
S = sn.spec["1993J_-11.0.txt"]

S.get_specphot(sn.phot.data_filters["BessellV"], verbose=True)
