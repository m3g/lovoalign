#include(FetchContent)
#FetchContent_Declare(
#  lapack
#  GIT_REPOSITORY https://github.com/Reference-LAPACK/lapack 
#  GIT_TAG        04b044e020a3560ccfa9988c8a80a1fb7083fc2e # release 3.12.0
#)
#FetchContent_MakeAvailable(lapack)
#find_package(lapack)

cmake_minimum_required(VERSION 2.8.12)
project(lovoalign Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    add_compile_options(
        -Wall "$<$<CONFIG:Debug>:-Werror>"
        -llapack -lblas
        -fmax-stack-var-size=1296000000
    )
endif()

# Build the executable
add_executable(lovoalign
        src/sizes.f90
        src/bijetype.f90
        src/inputpars.f90
        src/file_operations.f90
        src/initrandom.f90
        src/warnings.f90
        src/ioformat.f90
        src/lovoalign.f90
        src/computegdt.f90
        src/get_nearest.f90
        src/getpars.f90
        src/getrmsd2.f90
        src/getrmsd.f90
        src/help.f90
        src/initial.f90
        src/moveprot.f90
        src/nonbscore.f90
        src/orderpdb.f90
        src/orprot.f90
        src/printdata.f90
        src/procrustes.f90
        src/prodin.f90
        src/protall.f90
        src/pseudoprot.f90
        src/readfile.f90
        src/readlist.f90
        src/structal.f90
        src/title.f90
        src/tmscore.f90
        src/tocm.f90
        src/triang.f90
        src/writebije.f90
        src/writenonbije.f90
        src/writepdb.f90
        src/writermsf.f90
        src/readfasta.f90
        src/dist.f90
        src/dval.f90
        src/randomini.f90
        src/ic.f90
        src/ival.f90
        src/length.f90
        src/letter.f90
        src/flashsort.f90
        src/optim.f 
)

# Installation directive
install(TARGETS lovoalign DESTINATION bin)
