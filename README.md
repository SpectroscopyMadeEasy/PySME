![Python application](https://github.com/SpectroscopyMadeEasy/PySME/workflows/Python%20application/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/pysme-astro/badge/?version=latest)](https://pysme-astro.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5547527.svg)](https://doi.org/10.5281/zenodo.5547527)

# PySME

Spectroscopy Made Easy (SME) is a software tool that fits an observed
spectrum of a star with a model spectrum. Since its initial release in
[1996](http://adsabs.harvard.edu/abs/1996A%26AS..118..595V), SME has been a
suite of IDL routines that call a dynamically linked library, which is
compiled from C++ and fortran. This classic IDL version of SME is available
for [download](http://www.stsci.edu/~valenti/sme.html).

In 2018, we began began reimplmenting the IDL part of SME in python 3,
adopting an object oriented paradigm and continuous itegration practices
(code repository, build automation, self-testing, frequent builds).

# Installation

Stable release (recommended):

```bash
pip install pysme-astro
```

Install latest source from GitHub:

```bash
git clone --recurse-submodules https://github.com/SpectroscopyMadeEasy/PySME.git
cd PySME
pip install .
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

For full installation options (virtual environments, specific tags, uninstall), see:
https://pysme-astro.readthedocs.io/en/latest/getting_started/installation.html

# Poster

A poster about PySME can be found here: [Poster](http://sme.astro.uu.se/poster.html)

# GUI

A GUI for PySME is available in its own repository [PySME-GUI](https://github.com/SpectroscopyMadeEasy/PySME-GUI).

# Windows

PySME is supported on Windows via WSL2 (install and run PySME inside the Linux subsystem).
