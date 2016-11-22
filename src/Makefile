#
# Makefile for LovoAlign: Read the comments if you have some
#                         problem while compiling.
#
# If you want, set a different fortran compiler here:
#
FORTRAN = gfortran
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
MODULES = sizes.f90             \
          inputpars.f90         \
          file_operations.f90   \
          initrandom.f90        \
          ioformat.f90 
FILES = lovoalign.f90     \
        computegdt.f90    \
        get_nearest.f90   \
        getpars.f90       \
        getrmsd2.f90      \
        getrmsd.f90       \
        help.f90          \
        initial.f90       \
        moveprot.f90      \
        nonbscore.f90     \
        orderpdb.f90      \
        orprot.f90        \
        printdata.f90     \
        procrustes.f90    \
        prodin.f90        \
        protall.f90       \
        pseudoprot.f90    \
        readfile.f90      \
        readlist.f90      \
        structal.f90      \
        title.f90         \
        tmscore.f90       \
        tocm.f90          \
        triang.f90        \
        writebije.f90     \
        writenonbije.f90  \
        writepdb.f90      \
        writermsf.f90     \
        dist.f90          \
        dval.f90          \
        randomini.f90     \
        ic.f90            \
        ival.f90          \
        length.f90        \
        letter.f90        \
        jacobi.f          \
        flashsort.f90     \
        optim.f 
BIN=../bin
lovoalign : $(MODULES) $(FILES)
	@echo " ------------------------------------------------------ " 
	@echo " Compiling LovoAlign with $(FORTRAN) " 
	@echo " Flags: $(FLAGS) " 
	@echo " ------------------------------------------------------ " 
	@$(FORTRAN) $(FLAGS) $(MODULES) $(FILES) -o $(BIN)/lovoalign
	@echo " LovoAlign succesfully built." 
	@echo " ------------------------------------------------------ " 
devel: $(MODULES) $(FILES)
	@echo " ------------------------------------------------------ " 
	@echo " Compiling LovoAlign with $(FORTRAN) " 
	@echo " Flags: -Wall -fbounds-check" 
	@echo " ------------------------------------------------------ " 
	@$(FORTRAN) -Wall -fbounds-check $(MODULES) $(FILES) -o $(BIN)/lovoalign
	@echo " LovoAlign succesfully built." 
	@echo " ------------------------------------------------------ " 
clean:
	\rm -f *.mod *.o
cleanall:
	\rm -f $(BIN)/lovoalign *.mod *.o
