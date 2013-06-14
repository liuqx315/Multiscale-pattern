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
EXTEN = nvector \
        arkode \
        arkode/fcmix \
        arkode/pcg \
        arkode/spfgmr

# lib directories
LIBFILES = nvector/*.o \
           arkode/*.o \
           arkode/fcmix/*.o \
           arkode/pcg/*.o \
           arkode/spfgmr/*.o

# executable directories
EXEC = examples/arkode/serial \
       examples/arkode/parallel \
       examples/arkode/fcmix_serial \
       examples/cvode/serial \
       examples/ida/serial \
       regression/arkode \
       regression/cvode \
       regression/ida \
       regression/sundials

# build all library source code and tests
all : 
	for dir in ${EXTEN}; do \
	  $(MAKE) -C $$dir; \
	done
	ar -cr sundials_extensions.a ${LIBFILES}
	for dir in ${EXEC}; do \
	  $(MAKE) -C $$dir; \
	done

clean : 
	for dir in ${EXTEN}; do \
	  cd $$dir && $(MAKE) clean && cd -; \
	done
	for dir in ${EXEC}; do \
	  cd $$dir && $(MAKE) clean && cd -; \
	done


####### End of Makefile #######
