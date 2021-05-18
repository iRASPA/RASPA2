RASPA2
======

This software is a general purpose classical simulation package. It has been developed at
Northwestern University (Evanston, USA) during 2006-2011 in active collaboration
with University Pablo de Olavide (Seville, Spain), and the Technical University of Delft.

It can be used for the simulation of molecules in gases, fluids, zeolites, aluminosilicates,
metal-organic frameworks, carbon nanotubes and external fields.

Documentation: https://iraspa.org/raspa/

References
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

Force fields used as examples in RASPA
======================================

Note: all force field-, structure-, and molecule-files are meant as examples to illustrate the input and workings of the code.
For simulations on your own systems, you will have to construct your own classical force field and validate it.
The force fields used in the examples are based on the following force fields taken from literature:

small adsorbates
------------------
* CO2 using a shifted LJ potential (can be used both in Monte Carlo and in Molecular Dynamics).\
  "Transferable Force Field for Carbon Dioxide Adsorption in Zeolites",
  A. Garcia-Sanchez, C. O. Ania, J. B. Parra, D. Dubbeldam, T. J. H. Vlugt, R. Krishna, S. Calero, J. Phys. Chem. C 2009, 113, 8814-8820.
  [link to article](https://pubs.acs.org/doi/abs/10.1021/jp810871f)
* N2, O2, and Ar using shifted LJ potentials (can be used both in Monte Carlo and in Molecular Dynamics).\
  "Effect of air humidity on the removal of carbon tetrachloride from air using Cu–BTC metal–organic framework",
  A. Martin-Calvo, E. Garcia-Perez, A. Garcia-Sanchez, R. Bueno- Perez, S. Hamad, S. Calero, Phys. Chem. Chem. Phys. 2011, 13, 11165-11174.
  [link to article](https://pubs.rsc.org/en/content/articlelanding/2011/cp/c1cp20168a)
* CH4 using shifted LJ potentials (can be used both in Monte Carlo and in Molecular Dynamics).\
  "Effect of pressure, membrane thickness, and placement of control volumes on the flux of methane through thin silicalite membranes: A dual
control volume grand canonical molecular dynamics study",
  M. G. Martin, A. P. Thompson, T. M. Nenoff, J. Chem. Phys. 2001, 114, 7174-7181.
  [link to article](https://aip.scitation.org/doi/10.1063/1.1360256)
* Helium potential for computing void-fractions.\
  J.O. Hirschfelder, C.F. Curtiss, R.B. Bird, Molecular Theory of Gases and Liquids, Wiley, New York, 1954, p. 1114.
  [reference in article](https://www.sciencedirect.com/science/article/abs/pii/S0927775701006288)
* Linear/Branched alkane model.\
  "United atom force field for alkanes in nanoporous materials",
  D. Dubbeldam, S. Calero, T.J.H. Vlugt, R. Krishna, T.L.M. Maesen, B. Smit, J. Phys. Chem. B 2004, 108(33), 12301-12313
  [link to article](https://pubs.acs.org/doi/abs/10.1021/jp0376727)

rigid MOFs
------------
* "DREIDING: a generic force field for molecular simulations",
  S.L. Mayo, B.D. Olafson, W.A. Goddard, J. Phys. Chem. 1990, 94, 8897-8909
  [link to article](https://pubs.acs.org/doi/10.1021/j100389a010)
* "UFF, a full periodic table force field for molecular mechanics and molecular dynamics simulations",
  A.K. Rappé, C.J. Casewit, K.S. Colwell, W.A. Goddard, W.M. Skiff, J. Am. Chem. Soc. 1992, 114, 10024-10035.
  [link to article](https://pubs.acs.org/doi/10.1021/ja00051a040)

rigid zeolites
----------------
* "TraPPE-zeo: Transferable Potentials for Phase Equilibria Force Field for All-Silica Zeolites",
   P. Bai, M. Tsapatsis, J. I. Siepmann, J. Phys. Chem. C 2013, 117, 24375-24387.
  [link to article](https://pubs.acs.org/doi/10.1021/jp4074224)

flexible MOFs
---------------
* Flexible force fields for IRMOF-1, IRMOF-10, and IRMOF-16.
  "Exceptional Negative Thermal Expansion in Isoreticular Metal–Organic Frameworks",
  D. Dubbeldam, K.S. Walton, D.E. Ellis, R.Q. Snurr, Angew. Chem. Int. Ed. 2007, 46, 4496-4499.
  [link to article](https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.200700218)

flexible zeolites
-------------------
* "Molecular dynamics studies on zeolites. 3. Dehydrated zeolite A",
  P. Demontis, G. B. Suffritti, S. Quartieri, E. S. Fois, A. Gamba, J. Phys. Chem. 1988, 92, 867-871.
  [link to article](https://pubs.acs.org/doi/10.1021/j100315a003)
* "Molecular modeling of zeolite structure. 2. Structure and dynamics of silica sodalite and silicate force field",
  J. B. Nicholas, A. J. Hopfinger, F. R. Trouw, L. E. Iton, J. Am. Chem. Soc. 1991, 113, 4792-4800.
  [link to article](https://pubs.acs.org/doi/10.1021/ja00013a012)
* "Potential Functions for Silica and Zeolite Catalysts Based on ab Initio Calculations. 3. A Shell Model Ion Pair Potential for Silica and Aluminosilicates",
  K.-P. Schröder, J. Sauer, J. Phys. Chem. 1996, 110, 11043-11049.
  [link to article](https://pubs.acs.org/doi/abs/10.1021/jp953405s)

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

