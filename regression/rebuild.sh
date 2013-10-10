#!/bin/bash
#------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
#------------------------------------------------------------

# arkode regression tests
cd arkode
make
cd ..

# cvode regression tests
cd cvode
make
cd ..

# ida regression tests
cd ida
make
cd ..

# sundials regression tests
cd sundials
make
cd ..

