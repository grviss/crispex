# CRISPEX
[![Latest](https://img.shields.io/badge/Latest_release-v1.7.4-blue.svg?style=flat)](https://github.com/grviss/crispex/releases/tag/v1.7.4)

`CRISPEX` ([Vissers and Rouppe van der Voort (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...750...22V/abstract), [Löfdahl et al. (2018)](https://arxiv.org/abs/1804.03030)) is a graphical user interface written in the Interactive Data Language (IDL) for multi-dimensional (solar) data. 
It has been developed at the Institute of Theoretical Astrophysics (University of Oslo) and the Institute for Solar Physics (Stockholm University) in order to facilitate browsing and analysis of large data sets obtained with the Swedish 1-m Solar Telescope, and in particular the data obtained with the CRisp Imaging SpectroPolarimeter (CRISP) instrument. 

Since version 1.7 it has been extended to accomodate FITS files with data from the *Interface Region Imaging Spectrograph* (IRIS), shifting the focus to browsing and analysis of multi-instrument data sets. 
`CRISPEX` can handle any data cube regardless of its origin, provided it has been formatted according to certain rules.

## Installation
Simply clone the repository and execute `setup.sh` in the crispex directory.
Note that `setup.sh` should be re-run after every `git pull`.
