RASPA2
======

This software is a general purpose classical simulation package. It has been developed at
Northwestern University (Evanston, USA) during 2006-2011 in active collaboration
with University Pablo de Olavide (Seville, Spain), and the Technical University of Delft.

It can be used for the simulation of molecules in gases, fluids, zeolites, aluminosilicates,
metal-organic frameworks, carbon nanotubes and external fields.

Documentation: https://iraspa.org/raspa/

References
==========
* D. Dubbeldam, S. Calero, D.E. Ellis, and R.Q. Snurr, "RASPA: Molecular Simulation Software for Adsorption and Diffusion in Flexible Nanoporous Materials",
   Mol. Sim., 42 (2), 81-101, 2016.
   [link to article](https://www.tandfonline.com/doi/full/10.1080/08927022.2015.1010082)
* D. Dubbeldam, A. Torres-Knoop, and K.S. Walton,  "On the Inner Workings of Monte Carlo Codes",
   Mol. Sim. (special issue on Monte Carlo) 39(14-15), 1253-1292, 2013.
   [link to article](http://www.tandfonline.com/doi/full/10.1080/08927022.2013.819102)
* D. Dubbeldam and R.Q. Snurr, "Recent developments in the molecular modeling of diffusion in nanoporous materials",
   Mol. Sim., 33 (4-5), 305-325, 2007.
   [link to article](http://www.tandfonline.com/doi/abs/10.1080/08927020601156418)
* D. Dubbeldam, K.S. Walton, T.J.H. Vlugt, and S. Calero, "Design, Parameterization, and Implementation of Atomic Force Fields for Adsorption in Nanoporous Materials",
   Adv. Theory Simulat., 2(11), 1900135, 2019.
   [link to article](https://onlinelibrary.wiley.com/doi/full/10.1002/adts.201900135)

Acknowledgements
================
Development of this work is partially supported by the U.S. Department 
of Energy, Office of Basic Energy Sciences, Division of Chemical 
Sciences, Geosciences and Biosciences through the Nanoporous Materials 
Genome Center under award DE-FG02-17ER16362.

Installation
============

```
rm -rf autom4te.cache  
mkdir m4  
aclocal  
autoreconf -i  
automake --add-missing  
autoconf  

./configure --prefix=${RASPA_DIR}  
# or ./scripts/CompileScript/make-gcc-local  

make  
make install  
```

Conda-forge
-----------
raspa2 can be installed through the [conda](https://docs.conda.io/) package manager on Linux and MacOS (done by Leopold Talirz):

```
conda install -c conda-forge raspa2
export RASPA_DIR=/path/to/environment/root
# e.g. export RASPA_DIR=/Users/leopold/Applications/miniconda3/envs/raspa-env
```

Documentation
-------------

Build the latest RASPA manual using:
```
cd Docs
pdflatex raspa.tex
```

