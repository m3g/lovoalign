### LovoAlign ###

Lovoalign is a structural alignment package. The methods used for structural alignment are based on Low Order Value Optimization (LOVO) theory. The use of LOVO theory led to the development of fast convergent algorithms that provide very robust optimization of scoring functions.

Numerical experiments show that the LOVO algorithms implemented here provide the most reliable optimization of the STRUCTAL alignment while being very fast.

Simple input parameters can be used to align two structures, a single structure to a whole database, or to perform an all-on-all database structural alignment.

See https://m3g.github.io/lovoalign

The present site is contains current and previous versions of the package.

The alignment can be highly customized, since any subatom selection of the structures, main chains or not, can be aligned by using the "-beta and -ocup" input options.

### How to install and use ###

Currently the package is mostly developed for linux/unix users. Under those systems, instalation is done by:

1 - Download the lovoalign.DATE.tar.gz file.

2 - Unpack it and compile it:

```
    tar -xzvf lovoalign.DATE.tar.gz
    cd lovoalign/src
    make
```

(You must have the gfortran compiler and Lapack installed [[see these notes]](#installing-lapack)). This will create a "lovoalign" executable file in the "lovoalign/bin" directory. Add this directory to your path. 

The directory lovoalign/input contains an example input file for the "lovoalign.sh"
script, which is mostly self-explicative, but required the Visual Molecular Dynamics (VMD) package. With this input file, run lovoalign.sh with:

` lovoalign.sh lovoalign.inp `

Otherwise, the "lovoalign" executable can be run as a standalone program, as indicated by the instructions in the lovoalign main site.

### References ###

**Primary reference, please cite this work when using lovoalign:**

> L. Martinez, R. Andreani, J. M. Martinez.
> Convergent algorithms for protein structural alignment
> BMC Bioinformatics, 8:306 2007
> [Full Text](http://www.biomedcentral.com/1471-2105/8/306/abstract)

**Related work:**

> R. Andreani, J. M. Martinez, L. Martinez, F. Yano.
> Continuous Optimization Methods for Structural Alignment.
> Mathematical Pogramming, 112:93-124, 2008
> [article](http://www.springerlink.com/content/hv537728r6k15qu4/)

> R. Andreani, J. M. Martinez, L. Martinez
> Trust region superposition methods for protein alignment
> IMA Journal on Numerical Analysis, 28, 690-710, 2008.
> [article](http://imajna.oxfordjournals.org/cgi/content/abstract/drn021)

> R. Andreani, J. M. Martinez, L. Martinez, F. Yano.
> Low Order Value Optimization and Applications.
> Journal of Global Optimization, 43, 1-22, 2009.
> [article](http://www.springerlink.com/content/hw66734228844619/)

For full text articles go to http://www.ime.unicamp.br/~martinez/lovoalign

### Installing lapack

To compile `lovoalign`, you need first a fortran compiler. Typically the `gfortran` compiler can be installed using the package
manager of any linux system. On debian-based systems use, for example,
```
sudo apt install gfortran
```

You also need the `lapack` library. In debian-based systems, it can be installed with:
```
sudo apt install liblapack-dev
```
with that, just follow the compilation instructions provided above.

Other distributions probably have a lapack library available, but the package name might vary. Alternativelly, you can 
install and use local version of `lapack`, with:

1. Download and compile lapack:
   ```bash
   git clone https://github.com/Reference-LAPACK/lapack
   cd lapack
   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_LIBDIR=$HOME/.local/lapack ..
   cmake --build . -j --target install
   ```
2. Download lovoalign and edit the `lovoalign/src/Makefile` file:
   Change the `FLAGS` line to:
   ```
   FLAGS = -O3 -ffast-math -L$HOME/.local/lapack -llapack -lblas
   ```

3. Compile `lovoalign` as usual:
   ```
   cd lovoalign/src
   make
   ```

This should build the `lovoalign` executable in the directory `/lovoalign/bin`. 



