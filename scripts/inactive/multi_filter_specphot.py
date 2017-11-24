"""
adapted from single_filter_specphot.py
"""

import os
import pycoco as pcc

# snname = "SN1993J"
snname = "SN2007uy"
all_listfiles = pcc.coco.get_all_spec_lists()
all_listfiles
listfile = "/Users/berto/Code/CoCo/lists/"+snname+".list"
list_table = pcc.utils.read_list_file(listfile)
snname = list_table["snname"][0]
snname
sn = pcc.classes.SNClass(snname)
sn.load_phot(path=os.path.join(pcc.defaults._default_data_dir_path, "lc/" + snname + ".dat"))
sn.load_list(listfile)
sn.load_spec()
sn.check_overlaps()
S = sn.spec["2007uy_-5.06.txt"]
S.get_specphot(sn.phot.data_filters, verbose=True)


