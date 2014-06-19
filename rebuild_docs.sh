#!/bin/bash
# ----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ----------------------------------------------------------------
# Copyright (c) 2014, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# ----------------------------------------------------------------


# go to the repository directory
cd /home/dreynolds/www/arkode_dev

# update the ARKode_dev source code
hg pull
hg up

# Rebuild all ARKode documents
cd doc/guide
make clean
make html latexpdf
cd -

cd doc/examples
make clean
make html latexpdf
cd -

