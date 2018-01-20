"""

"""
import os

import pycoco as pcc
from astropy.table import Column, Table
import numpy as np

all_listfiles = pcc.coco.get_all_spec_lists()

for listfile in all_listfiles:
    print(listfile)

    list_table = pcc.utils.read_list_file(listfile)

    snname = list_table["snname"][0]
    print(snname)

    print(list_table)

    sn = pcc.classes.SNClass(snname)

    sn.load_phot(path=os.path.join(pcc.defaults._default_data_dir_path, "lc/" + snname + ".dat"))
    sn.load_list(listfile)

    # make_new_listfile = True
    make_new_listfile = False

    if make_new_listfile:
        # newpath = Column(
        #     ["/Users/berto/data/CoreCollapse/Blue_extension/spectra_blue_extended/" + snname + "/" + path.split("/")[-1] for
        #      path in sn.list["spec_path"]], name="spec_path")
        basedir = "/Users/berto/data/CoreCollapse/Blue_extension/spectra_blue_extended/"
        newpath = Column(
            [basedir + snname + "/" + i for i in os.listdir(os.path.join(basedir, snname))], name="spec_path")

        newmjdobs = Column([np.float(specname.split("_")[-1].replace(".spec", "")) for specname in newpath], name="mjd_obs")
        newsnnames = Column([snname for i in newmjdobs], name="snname")
        new_z = Column([list_table["z"][0] for i in newmjdobs])

        list_table_orig = sn.list
        # new_list_table = Table([newpath, list_table["snname"], list_table["mjd_obs"], list_table["z"]])
        new_list_table = Table([newpath, newsnnames, newmjdobs, new_z])
        sn.list = new_list_table

        print(new_list_table)
        outpath = os.path.join( "/Users/berto/data/CoreCollapse/Blue_extension/lists", snname+".list")

        new_list_table.write(outpath, format = "ascii.fast_commented_header", overwrite=True)

    listfile = os.path.join("/Users/berto/data/CoreCollapse/Blue_extension/lists", snname+".list")

    sn.load_list(listfile)
    print(sn.list)

    sn.load_spec()
    sn.check_overlaps()

    sn.get_lcfit(os.path.join(pcc.defaults._default_recon_dir_path, snname + ".dat"))
    sn.check_overlaps()

    # print(63)
    # print(len(sn.spec))

    for spec_key in sn.spec:
        print(spec_key)
        S = sn.spec[spec_key]
        S.get_specphot(sn.phot.data_filters, verbose=True)

        # new_spec = pcc.classes.SpectrumClass()
        new_spec = pcc.kcorr.donor_extend(S, sn, verbose=True)
        # new_spec.load_table(pcc.kcorr.linear_extend(S, return_table=True))
        # new_spec.load_table(pcc.kcorr.donor_extend(S, snobject=sn))

        # new_spec.save(filename = spec_key.replace(".spec", ".donor.spec"), path = os.path.join(pcc.defaults._default_data_dir_path, "spec_extended"), verbose=True)
        new_spec.save(filename = spec_key, path = os.path.join(pcc.defaults._default_data_dir_path, "spec_extended"), verbose=True)

