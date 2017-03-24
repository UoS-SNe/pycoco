from distutils.core import setup
import sys
import os
import re

packageName = 'pycoco'
packageDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          packageName)
versionFile = os.path.join(packageDir, 'version.py')



with open(versionFile, 'r') as f:
      s = f.read()

# Look up the string value assigned to __version__ in version.py using regexp
versionRegExp = re.compile("__version__ = \"(.*?)\"")

# Assign to __version__
__version__ =  versionRegExp.findall(s)[0]
setup(# package information
      name=packageName,
      version=__version__,
      description='Python tools for the CoCo templates',
      long_description=''' ''',
      # What code to include as packages
      packages=[packageName,
                packageName+".extinction",
                packageName+".kcorr",
                packageName+".sims",
                packageName+".utils",
                ],
      package_dir={packageName:'pycoco'}

      # What data to include as packages
      )
