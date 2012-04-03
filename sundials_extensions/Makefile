###################################################################
#  Daniel R. Reynolds, SMU Mathematics
#  Copyright 2010; all rights reserved
###################################################################

# Make.machine is used for setting machine-dependent compilation
# variables: e.g. compilers, libraries, F90-to-C conversion vars
include Make.machine

# append problem and interface preprocessor defs to compilation defs.
FFLAGS_ = $(FFLAGS) $(SUND_DEFS)
CFLAGS_ = $(CFLAGS) $(SUND_DEFS)


# extension directories
EXTEN = nvector


# build all library source code and tests
all : 
	for dir in ${EXTEN}; do \
	  $(MAKE) -C $$dir; \
	done
	ar -r sundials_extensions.a */*.o
	make -C tests

clean : 
	for dir in ${EXTEN}; do \
	  cd $$dir && $(MAKE) clean; \
	done
	cd tests && $(MAKE) clean


####### End of Makefile #######
