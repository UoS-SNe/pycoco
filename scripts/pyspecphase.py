#!/usr/bin/env python

"""

"""
import os
import pycoco as pcc


if __name__ == "__main__":
    phase_path = os.path.join(pcc.defaults._default_coco_dir_path, "examples/phase.list")
    pcc.coco.run_specphase("BessellV", phase_path)
else:
    pass