# **`pycoco`**
___

## v0.9.0
___
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/) [![Build Status](https://travis-ci.org/RobFirth/pycoco.svg?branch=master)](https://travis-ci.org/RobFirth/pycoco)[![DOI](https://zenodo.org/badge/74136059.svg)](https://zenodo.org/badge/latestdoi/74136059)
___
This is the development repo for the python frontend for the core-collapse SNe template code 'CoCo':  https://github.com/UoS-SNe/CoCo
(my fork is currently https://github.com/RobFirth/CoCo).

CoCo was originally started by Natasha Karpenka, and is currently being updated and maintained by S. Prajs (https://github.com/SzymonPrajs).

A paper, Firth et. al. 2017, is currently in prep.
___
 
 * Now available on PyPi and via `pip` - packagen available to install as `pycocosn`
 
 * Extending templates now possible with Black Body spectrum, flat and linear

 * calling CoCo (or reproduced CoCo functions) now more straightforward

 * Added more test cases to test_pycoco

 * Travis-CI implemented and Master and dev are passing

 * Dependancies now handled in setup.py

 * Calling all parts of CoCo from pycoco now operational

 * Migrated to better, and clearer code structure and sub-modules

 * Mangling now done within python, rather than C++ in CoCo

 * File I/O and interaction with `CoCo LCfit` output now operational - 11/01/17

 * Usage of a `SN` class now solid

 * Implementing calls to `CoCo LCfit` now

 * Calls to `CoCo specfit` now implemented

 * Mangling now stable

 * Adding tools for calculating magnitude offsets to make it easier to import new data

 * regeneration of filter list file now possible

 * added colours and bandpasses for LSST filters

 * Less dependence on environment variables

 * installation through `setup.py` now possible

 * better handing of CoCo sim outputs

 * SN position and mu now stored in infofile `./testdata/info/info.dat`

 * dark sky calculations and integration with LSST Throughputs now done

 * Improved stability

 * input from `astropy tables` now more straightforward - better integration with `coco.simulate`

 * can now batch fit light curves and spectra from within python
___

To install:

```
pip install pycocosn

```

To install from source:

```
git clone https://github.com/RobFirth/verbose-enigma.git
```

then:

```
cd verbose-enigma
python setup.py install --user
```

(The --user argument only installs current user only, omitting flag will install for all users on the system if there are appropriate permissions)

_NOTE: make sure that the python used to install is the one that you will use with pycoco_
___

Ideally set the following environment variables:

`COCO_ROOT_DIR` (my default is `~/Code/CoCo/`)
`PYCOCO_FILTER_DIR`(my default is `~/Code/CoCo/data/filters/`)
`PYCOCO_DATA_DIR` (my default is `~/Code/CoCo/data/`)
`SFD_DIR` (my default is `~/data/Dust/sfddata-master/`; see below)
`LSST_THROUGHPUTS` (my default is `${HOME}/projects/LSST/throughputs`)
`LSST_THROUGHPUTS_BASELINE` (my default is `${LSST_THROUGHPUTS}/baseline`)

also `pycoco/` and `CoCo` need to be in your path and pythonpath, i.e.:

 ```
 setenv PATH /Users/berto/Code/pycoco/:$PATH

 setenv PYTHONPATH "/Users/berto/Code/pycoco/:$PYTHONPATH
 ```

___

# Requirements
## python packages

* matplotlib
* numpy
* scipy
* astropy
* lmfit
* sfdmap (https://github.com/kbarbary/sfdmap)

# additionally

* lsst throughputs (https://github.com/lsst/throughputs; `$LSST_THROUGHPUTS` points to this directory)
___


for `sfdmap`, the environment variable `SFD_DIR` needs to point at the path to the parent directory of the appropriate dust map files. See the installation instructions here: https://github.com/kbarbary/sfdmap
___  
