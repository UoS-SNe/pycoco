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

setup(# package information
      name="pycocosn",
      version=__version__,
      author="Rob Firth",
      author_email="science@robertfirth.co.uk",
      url="https://github.com/RobFirth/pycoco",
      description='Python tools for the CoCo templates',
      long_description='''Python tools for the CoCo templates''',
      packages=["pycoco",],
      package_dir={"pycoco":"pycocosn"},
      package_data={"pycoco":['kcorr_data/*']},
      install_requires=['numpy', 'matplotlib', 'pandas', 'sqlalchemy', 'astropy>=2.0', 'sfdmap', 'lmfit']
      )
