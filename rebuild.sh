#!/bin/bash

# update the ARKode_dev source code
hg pull
hg up

# reconfigure and rebuild ARKode library -- requires 
# sundials to have previously been configured with cmake, 
# with directory structure:
#    <arkode>/sundials -- source directory (already exists)
#    <arkode>/build    -- build directory
#    <arkode>/install  -- installation directory
cd build
cmake ../sundials/
make
make install
