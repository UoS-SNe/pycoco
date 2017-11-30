try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os
import re

packageName = 'pycocosn'
packageDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          packageName)
versionFile = os.path.join(packageDir, 'version.py')

with open(versionFile, 'r') as f:
      s = f.read()

# Look up the string value assigned to __version__ in version.py using regexp
versionRegExp = re.compile("__version__ = \"(.*?)\"")

# Assign to __version__
__version__ =  versionRegExp.findall(s)[0]

## Compile list of data files to be installed
# data_files = [("lc", [ os.path.join(packageDir, "data/lc", i) for i in os.listdir(os.path.join(packageDir, "data/lc"))] ),
#               ("filters", [ os.path.join(packageDir, "data/filters", i) for i in os.listdir(os.path.join(packageDir, "data/filters"))] )
#              ]
#
# for i in os.listdir(os.path.join(packageDir, "data/spec")):
#     path = os.path.join("spec/", i)
#     datalist = []
#     for j in os.listdir(os.path.join(packageDir, "data/", path)):
#         datalist.append(os.path.join(packageDir, "data", path, j))
#
#     datatuple = (path, datalist)
#     data_files.append(datatuple)

setup(# package information
      name=packageName,
      version=__version__,
      description='Python tools for the CoCo templates',
      long_description=''' ''',
      # What code to include as packages
    #   packages=[packageName,
    #             packageName+".classes",
    #             packageName+".colours",
    #             packageName+".defaults",
    #             packageName+".extinction",
    #             packageName+".errors",
    #             packageName+".functions",
    #             packageName+".kcorr",
    #             packageName+".sims",
    #             packageName+".utils"
    #             ],
          packages=[packageName,
                    ],
      package_dir={packageName:'pycocosn'},
      package_data={packageName:['kcorr_data/*']},
      install_requires=['numpy', 'matplotlib', 'pandas', 'sqlalchemy', 'astropy', 'sfdmap', 'lmfit']  # ,

    #   include_package_data = True,
    #   data_files = data_files
      )
