
n.n.n / 2020-12-19
==================

  * Merge pull request #15 from ltalirz/cleanup
  * Updated VTK visualization to VTK 9
  * update version number
  * repository cleanup

v2.0.39 / 2020-12-05
====================

  * Bugfix: fixed the computation of the error bar (they were too high). 
  * Fixed compilation under macOS Big Sur.

v2.0.38 / 2020-05-27
====================

  * Updated the submit-scripts for large batch simulations in the 'scripts' directory
  * Updated partition factors for ReactionEnsembleAmmonia example
  * Added examples from a forthcoming bookchapter
  * Changed CO2 and N2 to shifted potentials in GenericZeolites and GenericMOFs
  * Modified RXMC to avoid under/overflow in acceptance rule. 
  * Added option 'LnPartitionFunction' to enter ln(q/V) in units of ln(A^(-3)).
  * Added reference to the README
  * Fixed type: "gr/gr" should be "mg/g" 
  * Added ITQ-29-2x2x2 
  * Added carbon monoxide
  * Incorporated changes to make it possible to change the volume of the unit cell along a particular direction/directions.
  * Do not check for energy-drift for MuPT-, MuPTPR-, and MuVT-ensembles
  * Fixed typo from 'J. Phys. Chem. B 1999, 103, 4508-4517' and changed torsion parameters 335.03 to 355.03 (as in 'J. Phys. Chem. B 1998, 102, 2569-2577')

v2.0.37 / 2019-04-27
====================

  * Added FrameworkChangeMove CFCMC compatibility
  * Prepared for upload to github

v2.0.36 / 2019-01-06
====================

  * Prepared for upload to github

v2.0.35 / 2018-07-22
====================

  * Added CFWidomProbability move; Fixed binary-restart to be backwards-compatible (from 2.0.35 onwards)

v2.0.34 / 2018-07-05
====================

  * Set pseudo-atom to charged if a framework-atom from a CIF-file is charged, even though the pseudo-atom was defined without charge
  * Fixed compiler error on some compilers (reported by N. Scott Bobbitt)

v2.0.33 / 2018-06-29
====================

  * Updated bookkeeping for CFCMC and CB/CFCMC with grids

v2.0.32 / 2018-06-29
====================

  * PDB-movies and CFCMC: Lambda is printed as the occupancy 
  * Added CRYST1-line to allcomponents pdb 
  * Modified VTK-files from 'SetAAFrames' to 'SetMultiSamples' (needed on some nvidia cards on linux)
  * Python version bump

v2.0.31 / 2018-05-08
====================

  * Added LTA dcTST example
  * Minor fixes in documentation and example molecule definitions

v2.0.30 / 2018-03-25
====================

  * Fixed degrees-of-freedom in GibbsIdentityChangeAdsorbateMove/GibbsIdentityChangeCationMove

v2.0.29 / 2018-03-14
====================

  * Fixed Intra 1-4 bug for molecules (Reported by Hirad S. Salehi)

v2.0.28 / 2018-03-05
====================

  * Changed criteria to detect identical framework-atoms when expanding spacegroup

v2.0.27 / 2017-12-22
====================

  * Added examples for MuPT, MuPTPR and MuVT ensemble MD simulations.
  * Implemented MuVT, MuPT, and MuPTPR ensembles for MD simulations.

v2.0.26 / 2017-09-27
====================

  * Fixed adsorbate/cations Fourier energy for CFCMC

v2.0.25 / 2017-09-06
====================

  * Added host-guest and guest-guest to SampleEnergyHistogram

v2.0.24 / 2017-07-17
====================

  * Fixed error in 'AddRealMatrix6x6' and 'SubtractRealMatrix6x6' (they are now symmetric)

v2.0.23 / 2017-07-12
====================

  * Removed erroneous debug-stop in 'SelectRandomMoleculeOfType'

v2.0.22 / 2017-06-26
====================

  * Report chemical potentials for individual components when using multiple phases in the new CFCMC GE

v2.0.20 / 2017-06-18
====================

  * Fixed order of deletion in CFCRXMCLambdaChangeMove (previous fix only fixed half)
  * Added Tutorial in Examples and Docs

v2.0.19 / 2017-05-09
====================

  * Changed condition for overlap for spacegroup-expanding back to previous definition

v2.0.18 / 2017-04-04
====================

  * CFGibbsSwapFractionalMoleculeToOtherBoxMove: extended to handle flexible molecules

v2.0.17 / 2017-04-01
====================

  * Extended CFGibbsFractionalToIntegerMove to multiple components

v2.0.16 / 2017-04-01
====================

  * Fixed order of deletion in CFCRXMCLambdaChangeMove

v2.0.15 / 2017-03-21
====================

  * Fixed errors in Hessian cation-cation minimization

v2.0.14 / 2017-02-21
====================

  * Fixed memory allocation-error in core-shells (reported by Ismael)

v2.0.13 / 2017-02-06
====================

  * Corrected CFMoleculePresent for new Gibbs CFCMC

v2.0.12 / 2016-12-18
====================

  * Added check for enough reaction molecules in Gibbs-swap move
  * Omit CFCRXMCLambdaChangeMove when no reactions defined

v2.0.11 / 2016-12-05
====================

  * Added check for enough reaction molecules in identity-change move
  * Fixed segmentation-fault in 'shell'-output
  * Added 'shell'-output for frameworks

v2.0.10 / 2016-08-29
====================

  * Error fixed (by Ozgur Yazaydin) in BoxShapeChangeMove
  * Gibbs CFCMC example added using the new methodology of Poursaeidesfahani et al. (J. Chem. Theory Comput., 2016, 12 (4), pp 1481-1490)
  * Small corrections for spectra, and cases without adsorbates
  * Merge branch 'master' of ssh://itsaio.science.uva.nl/~raspa2-dev/RASPA-2.0
  * Added ortho C11, C12, C13 strains
  * Added option to measure properties for CFCMC only when lambda is small ('MeasureLambdaBelow')

v2.0.9 / 2016-07-04
===================

  * Added option to measure properties for CFCMC only when lambda is small ('MeasureLambdaBelow')

v2.0.8 / 2016-05-30
===================

  * Fixes NPH-PR energy drift (reported by Ozgur Yazaydin
  * Fixed limit length of file-name in 'output.c'
  * Patch by Jure Varlec on configure.ac to allow specification different versions of blas, lapack and fftw

v2.0.7 / 2016-04-14
===================

  * Bug fixed by Ozgur Yazaydin in routine 'CalculateEnergyDifferenceFrameworkMoveCharge'
  * Added labeling of pseudo-atoms when reading cif-files from the chemical element when no explicit label (_atom_site_label) is given

v2.0.6 / 2016-03-05
===================

  * Added modified TraPPE potential (Salvador Rodriguez Gomez and Ismael Matito Martos)
  * Fixed error in binary restart for core-shell flexible framework

v2.0.5 / 2016-02-24
===================

  * Error in tail-corrections introduced in 2.0.2 fixed (reported by Ismael Matito Martos)
  * Added lambda-extrapolation to estimate the chemical potential for grand-canonical CFCMC

v2.0.4 / 2016-02-12
===================

  * Added lambda-extrapolation to estimate the chemical potential for grand-canonical CFCMC

v2.0.3 / 2016-02-11
===================

  * Fixed extrapolation in new Gibbs-CFCMC method and corrected inverse-density in ideal-gas term

v2.0.2 / 2016-02-11
===================

  * Fixed error in non-default number of trial-moves for CBMC-Gibbs; Added new CFCMC-method for Gibbs (Ali Poursaeidesfahani, Ariana Torres-Knoop, David Dubbeldam, Thijs J.H. Vlugt).
  * Fixed spaces in pdb-output

v2.0.1 / 2015-11-26
===================

  * Fixed enthalpy of adsorption for mixture and cations
  * Corrected an small error in the ideal-gas part of the fluctuation formula for the enthalpy of adsorption (at 300K it was about 0.5 kJ/mol off)
  * Added MIE potential with cutoff
  * Rewrite of Numerical-Recipes matrix-inverse as zero-based
  * Fixed error in cation-cation intra Coulomb total energy (reported by José Manuel Vicent Luna)
  * Added: enthalpy of adsorption corrected for CFMC
  * Added: total heat of adsorption from components of the mixture and their measured mol-fractions
  * Fixed binary-restart 'exchange the CF fractional particle'
  * Added: Enthalpy of adsorption for mixtures (Ali Poursaeidesfahani)
  * Add torsion potential: TRAPPE_DIHEDRAL_EXTENDED [5-parameters, p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)]
  * Added MC-move: to exchange the CF fractional particle of a random component with a random integer molecule of the same component
  * Uploaded to PyPI, can now be installed with 'pip install RASPA2'
  * ln -s COPYING LICENSE.txt
  * Updated README with links to public wiki
  * Recommiting Yamil Colon's warnings.c/h changes to new repo
  * Change default overlap-criteria in FreeEnergyProfile3D to remove artifacts in VTK pictures
  * Added PrintEnergyStatusToStdOut and changed DebugEnergyStatus to print energy in Kelvin
  * Added example 'UmbrellaSampling'
  * Added example 'minimization of MIL-53 with water'
  * Added pdf of the RASPA Molecular Simulation article
  * Added averages to SampleProjectedLengthsDistributionFunction()
  * Added SampleCOMDensityProfile3DVTKGrid to output the center-of-mass of molecule
  * Added pdf 'On the Inner workings of Monte Carlo codes'

v2.0 / 2015-02-03
=================

  * Initial Release RASPA-2.0
