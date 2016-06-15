# configure generated Makefile
# configure generated Makefile
#
# Makefile for LovoAlign: Read the comments if you have some
#                         problem while compiling.
#
# You may use the ./configure script to search automatically for
# some fortran compiler.
#
# This make file will try to compile LovoAlign with the default
# fortran compiler, defined by the FC directive. For doing this,
# just type
#
#          make 
#
FORTRAN = /usr/bin/gfortran
#
FLAGS = -O3 -ffast-math 
#
###################################################################
#                                                                 #
# Generally no modifications are required after this.             #
#                                                                 #
###################################################################
#
# Get the default fortran compiler
#
ifeq ($(FORTRAN),AUTO)
FORTRAN = $(FC)
endif 
FILES = lovoalign.f90 title.f90 flashsort.f jacobi.f optim.f 
BIN=../bin
lovoalign : $(FILES)
	@echo " ------------------------------------------------------ " 
	@echo " Compiling LovoAlign with $(FORTRAN) " 
	@echo " Flags: $(FLAGS) " 
	@echo " ------------------------------------------------------ " 
	@$(FORTRAN) $(FLAGS) $(FILES) -o $(BIN)/lovoalign
	@echo " LovoAlign succesfully built." 
	@echo " ------------------------------------------------------ " 
devel: $(FILES)
	@echo " ------------------------------------------------------ " 
	@echo " Compiling LovoAlign with $(FORTRAN) " 
	@echo " Flags: -Wall -fbounds-check" 
	@echo " ------------------------------------------------------ " 
	@$(FORTRAN) -Wall -fbounds-check $(FILES) -o $(BIN)/lovoalign
	@echo " LovoAlign succesfully built." 
	@echo " ------------------------------------------------------ " 
clean:
	\rm -f *.mod *.o
cleanall:
	\rm -f $(BIN)/lovoalign *.mod *.o
