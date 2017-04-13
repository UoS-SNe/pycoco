# verbose-enigma
___

#v0.5.3#
___

This is the development repo for the python frontend for the core-collapse SNe template code 'CoCo':  https://github.com/UoS-SNe/CoCo
(my fork is currently https://github.com/RobFirth/CoCo).

CoCo was originally started by Natasha Karpenka, and is currently being updated and maintained by S. Prajs (https://github.com/SzymonPrajs).

A paper, Firth et. al. 2017, is currently in prep.
___


 * File I/O and interaction with `CoCo LCfit` output now operational - 11/01/17

 * Usage of a `SN` class now solid

 * Implementing calls to `CoCo LCfit` now

 * Calls to `CoCo specfit` now implemented

 * Mangling now stable

 * Adding tools for calculating magnitude offsets to make it easier to import new data

 * regeneration of filter list file now possible
___

Ideally set the following environment variables:

`COCO_ROOT_DIR` (my default is `~/Code/CoCo/`)
`PYCOCO_FILTER_DIR`(my default is `~/Code/CoCo/data/filters/`)
`PYCOCO_DATA_DIR` (my default is `~/Code/CoCo/data/`)

also `pycoco/` and `CoCo` need to be in your path and pythonpath, i.e.:

 `setenv PATH /Users/berto/Code/pycoco/:$PATH`
 `setenv PYTHONPATH "/Users/berto/Code/pycoco/:$PYTHONPATH`

___

for `sfdmap`, the environment variable `SFD_DIR` needs to point at the path to the parent directory of the appropriate dust map files. See the installation instructions here: https://github.com/kbarbary/sfdmap

___  
