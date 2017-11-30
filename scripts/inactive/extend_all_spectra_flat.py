"""

"""
import os

import pycoco as pcc

all_listfiles = pcc.coco.get_all_spec_lists()

for listfile in all_listfiles:
    list_table = pcc.utils.read_list_file(listfile)

    snname = list_table["snname"][0]
    print(snname)

    sn = pcc.classes.SNClass(snname)
    sn.load_phot(path=os.path.join(pcc.defaults._default_data_dir_path, "lc/" + snname + ".dat"))
    sn.load_list(listfile)

    sn.load_spec()
    sn.get_lcfit(os.path.join(pcc.defaults._default_recon_dir_path, snname + ".dat"))
    sn.check_overlaps()

    for spec_key in sn.spec:
        print(spec_key)
        S = sn.spec[spec_key]
        S.get_specphot(sn.phot.data_filters, verbose=True)

        new_spec = pcc.classes.SpectrumClass()

        new_spec.load_table(pcc.kcorr.flat_extend(S, return_table=True))

        new_spec.save(filename = spec_key.replace(".txt", ".flat.txt"), path = os.path.join(pcc.defaults._default_data_dir_path, "spec_extended"))
