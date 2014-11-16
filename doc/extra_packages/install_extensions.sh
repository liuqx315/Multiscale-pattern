#!/bin/bash
# ----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ----------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# ----------------------------------------------------------------

# Sphinx
echo "installing Sphinx"
tar -zxf Sphinx-1.2b1.tgz
cd Sphinx-1.2b1
python ./setup.py install
cd -

# Sphinx fortran domain
echo "installing Sphinx fortran domain"
tar -zxf sphinx-fortran-extension.tgz
cd sphinx-fortran-extension
python ./setup.py install
cd -

# Sphinx bootstrap theme
echo "installing Sphinx bootstrap theme"
tar -zxf sphinx-bootstrap-theme.tgz
cd sphinx-bootstrap-theme
python ./setup.py install
cd -

echo "finished"
