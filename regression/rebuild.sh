#!/bin/bash

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

