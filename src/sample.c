/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'rigid.c' is part of RASPA-2.0

    RASPA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RASPA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *************************************************************************************************************/

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/stat.h>
//#include <fftw3.h>
#include "vector.h"
#include "potentials.h"
#include "simulation.h"
#include "molecule.h"
#include "cbmc.h"
#include "utils.h"
#include "inter_energy.h"
#include "internal_force.h"
#include "inter_force.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "internal_energy.h"
#include "ewald.h"
#include "cubic_spline_1d.h"
#include "grids.h"
#include "sample.h"
#include "rigid.h"
#include "molecule_properties.h"
#include "cubic_spline_1d.h"
#include "recrossing.h"
#include "thermo_baro_stats.h"
#include "spacegroup.h"
#include "statistics.h"
#include "integrate.h"
#include "integration.h"
#include "mc_moves.h"
#include "movies.h"

// An external var for saving PSD simulation output
extern char *PORE_SIZE_DISTRIBUTION_OUTPUT;
extern size_t PORE_SIZE_DISTRIBUTION_OUTPUT_SIZE;
extern bool STREAM;

// sampling the radial distribution function (RDF)
int *ComputeRDF;
int *WriteRDFEvery;
static REAL *CountRDF;
int *RDFHistogramSize;
REAL *RDFRange;
static REAL ****RadialDistributionFunction;
static REAL ***RadialDistributionFunctionWithFramework;

int *ComputeProjectedLengths;
int *WriteProjectedLengthsEvery;
static REAL **CountProjectedLengths;
int *ProjectedLengthsHistogramSize;
REAL *ProjectedLengthsRange;
static VECTOR ***ProjectedLengthsDistributionFunction;
static VECTOR **ProjectedLengthsAverage;

int *ComputeProjectedAngles;
int *WriteProjectedAnglesEvery;
static REAL **CountProjectedAngles;
int *ProjectedAnglesHistogramSize;
REAL *ProjectedAnglesRange;
static VECTOR ***ProjectedAnglesDistributionFunction;

//----------------------------------------------------------------------------------------
// CFC-RXMC : sampling lambda histogram
//----------------------------------------------------------------------------------------

int *ComputeCFCRXMCLambdaHistogram;
int *WriteCFCRXMCLambdaHistogramEvery;
int *CFCRXMCLambdaHistogramBins;
static REAL ***CFCRXMCLambdaHistogram;

//----------------------------------------------------------------------------------------


// sampling the number-of-molecules histogram
int *ComputeNumberOfMoleculesHistogram;
int *WriteNumberOfMoleculesHistogramEvery;
int *NumberOfMoleculesHistogramSize;
REAL *NumberOfMoleculesRange;
static REAL ***NumberOfMoleculesHistogram;

// sampling position histograms/free energies
int *ComputePositionHistogram;
int *WritePositionHistogramEvery;
int *PositionHistogramSize;
int *PositionHistogramMappingType;
static VECTOR ****PositionABCHistogram;
static VECTOR ****Position2DDiagonalHistogram;
static VECTOR ****Position2DDiagonalHistogram2;
static VECTOR4 ****Position3DDiagonalHistogram;

// samples the free energy profiles in a,b,c directions
int *ComputeFreeEnergyProfile;
int *WriteFreeEnergyProfileEvery;
int *FreeEnergyHistogramSize;
int *FreeEnergyMappingType;
static REAL (***RosenBinSum)[13];
static REAL (***RosenBinSumSquared)[13];
static REAL (***RosenBinCount)[13];

// sampling the pore-size distribution (PSD)
int *ComputePSDHistogram;
int *PSDHistogramSize;
REAL *PSDRange;
int *WritePSDHistogramEvery;
static REAL **PoreSizeDistributionHistogram;

// sampling the end-to-end histograms
int *ComputeEndToEndDistanceHistogram;
int *WriteEndToEndDistanceHistogramEvery;
int *EndToEndHistogramSize;
REAL *EndToEndRange;
static REAL ***EndToEndDistanceHistogram;

// sampling the energy histogram
int *ComputeEnergyHistogram;
int *WriteEnergyHistogramEvery;
int *EnergyHistogramSize;
REAL *EnergyHistogramLowerLimit;
REAL *EnergyHistogramUpperLimit;
static REAL ***EnergyHistogram;

// sampling the thermodynamic factor
int *ComputeThermoDynamicFactor;
int *WriteThermoDynamicFactorEvery;
static REAL ***ThermoDynamicFactorNumberOfMoleculesCrossTerm;
static REAL **ThermoDynamicFactorNumberOfMolecules;
static REAL *ThermoDynamicFactorNumberOfSamples;

// sampling the inter-framework spacing histogram
int *ComputeFrameworkSpacingHistogram;
int *WriteFrameworkSpacingHistogramEvery;
int *FrameworkSpacingHistogramSize;
REAL *FrameworkSpacingRange;
REAL ***OriginalFrameworkShiftDirAvg;
VECTOR ***OriginalFrameworkShift;
static REAL *****FrameworkDistanceHistogram;

// sampling histograms of the residence times
int *ComputeResidenceTimes;      // whether to compute the residence times or not
int *WriteResidenceTimesEvery;   // writes the output every 'WriteResidenceTimesEvery' times
REAL **ResidenceTimesHistogram;  // the data for the histogram
int *ResidenceTimesHistogramSize; // the number of elements of the histogram
REAL *RangeResidenceTimes;     // the maximum range of the histogram
static long long **ResidenceOriginAdsorbates;
static long long **ResidenceOriginCations;
static int **ResidenceStatusAdsorbates;
static int **ResidenceStatusCations;
static REAL (**ResidenceTimesFractionAdsorbates)[NR_BLOCKS];
static REAL (**ResidenceTimesFractionCations)[NR_BLOCKS];
static REAL *ResidenceTimesFractionCounts;

// sampling histograms of the distance between 2 selected atoms
int *ComputeDistanceHistograms;               // whether to compute the distance histograms or not
int *WriteDistanceHistogramsEvery;            // writes the output every 'WriteDistanceHistogramsEvery' times
int NumberOfElementsDistanceHistogram;        // the number of elements of the histogram
REAL MaxRangeDistanceHistogram;               // the maximum range of the histogram
int *NumberOfDistanceHistogramDefinitions;
int (**DistanceHistogramDefinitions)[2][3];
ATOM* (**DistanceHistogramPairs)[2];
REAL ***DistanceHistograms;

// sampling histograms of the bend angle between 3 selected atoms
int *ComputeBendAngleHistograms;               // whether to compute the angle histograms or not
int *WriteBendAngleHistogramsEvery;            // writes the output every 'WriteAngleHistogramsEvery' times
int NumberOfElementsBendAngleHistogram;        // the number of elements of the histogram
REAL MaxRangeBendAngleHistogram;               // the maximum range of the histogram
int *NumberOfBendAngleHistogramDefinitions;
int (**BendAngleHistogramDefinitions)[3][3];
ATOM* (**BendAngleHistogramPairs)[3];
REAL ***BendAngleHistograms;

// sampling histograms of the dihedral angle between 4 selected atoms
int *ComputeDihedralAngleHistograms;               // whether to compute the angle histograms or not
int *WriteDihedralAngleHistogramsEvery;            // writes the output every 'WriteDihedralAngleHistogramsEvery' times
int NumberOfElementsDihedralAngleHistogram;        // the number of elements of the histogram
REAL MaxRangeDihedralAngleHistogram;               // the maximum range of the histogram
int *NumberOfDihedralAngleHistogramDefinitions;
int (**DihedralAngleHistogramDefinitions)[4][3];
ATOM* (**DihedralAngleHistogramPairs)[4];
REAL ***DihedralAngleHistograms;

// sampling histograms of the angle between two planes (each formed by 3 chosen atoms)
int *ComputeAngleBetweenPlanesHistograms;
int *WriteAngleBetweenPlanesHistogramsEvery;
int NumberOfElementsAngleBetweenPlanesHistogram;
REAL MaxRangeAngleBetweenPlanesHistogram;
int *NumberOfAngleBetweenPlanesHistogramDefinitions;
int (**AngleBetweenPlanesHistogramDefinitions)[6][3];
ATOM* (**AngleBetweenPlanesHistogramPairs)[6];
REAL ***AngleBetweenPlanesHistograms;

// sampling molecular properties (bond distance, bend angle, dihedral angle)
int *ComputeMoleculeProperties;
int *WriteMoleculePropertiesEvery;
int *BondLengthHistogramSize;
REAL *BondLengthRange;
int *BendAngleHistogramSize;
REAL *BendAngleRange;
int *DihedralHistogramSize;
REAL *DihedralRange;
REAL ****BondLengthHistogram;
REAL ****UreyBradleyLengthHistogram;
REAL ****BendAngleHistogram;
REAL ****TorsionAngleHistogram;
REAL ****TorsionConformationHistogram;
REAL ***FrameworkBondLengthHistogram;
REAL ***FrameworkUreyBradleyLengthHistogram;
REAL ***FrameworkBendAngleHistogram;
REAL ***FrameworkTorsionAngleHistogram;
REAL **FrameworkAverageBondLength;
REAL **FrameworkBondLengthCount;
REAL **FrameworkAverageBendAngle;
REAL **FrameworkBendAngleCount;
REAL **FrameworkAverageTorsionAngle;
REAL **FrameworkTorsionAngleCount;

// sampling the IR spectra (spacings: 2048, 4196, 8192, 16384, 32768 points)
int *ComputeInfraRedSpectra;
int *WriteInfraRedSpectraEvery;
int SampleEveryInfraRed;
REAL ****Spectrum;
REAL ****SpectrumAverage;
REAL ****UnweightedSpectrum;
REAL ****UnweightedSpectrumAverage;
REAL *****SpectrumPseudoAtoms;
REAL *****SpectrumPseudoAtomsAverage;
REAL *sumw;
REAL **SpectrumCount;

// sampling the mean-squared displacement using a modified order-N algorithm
int *ComputeMSDOrderN;                           // whether or not to compute the msd
int *SampleMSDOrderNEvery;                       // the sample frequency
int *WriteMSDOrderNEvery;                        // write output every 'WriteMSDOrderNEvery' steps
int *NumberOfSitesMSDOrderN;                     // the number of dfferent sites
int NumberOfBlockElementsMSDOrderN;              // the number of elements per block
int MaxNumberOfBlocksMSDOrderN;                  // the maxmimum amount of blocks (data beyond this block is ignored)
int ComputeIndividualMSDOrderN;                  // whether or not to compute (self-)msd's for individual molecules
int ComputeSiteTypeMSDOrderN;                    // whether or not to compute (self-)msd's for individual molecules
int ComputeMSDOrderNPerPseudoAtom;               // whether or not to compute (self-)msd's for (pseudo-)atoms
static int *CountMSDOrderN;                      // counter for the amount of msd's measured
static int *NumberOfBlocksMSDOrderN;             // the current number of blocks in use
static int **BlockLengthMSDOrderN;               // the current length of the blocks
static VECTOR ****BlockDataMSDOrderN;            // array for all blocks containing all the molecule positions per component
static VECTOR ****BlockDataMSDOrderNOnsager;     // array for all blocks containing the summed molecule positions per component
static VECTOR ****MsdOrderN;                     // the self-msd in x,y,z per component
static REAL ****MsdOrderNDirAvg;                 // the msd directionally averaged per component
static VECTOR ****MsdOrderNPerMolecule;          // the self-msd in x,y,z per molecule
static REAL ****MsdOrderNPerMoleculeDirAvg;      // the self-msd directionally averaged per molecule
static REAL ****MsdOrderNCount;                  // counter for the amount of self msd-samples
static REAL ****MsdOrderNCountPerMolecule;       // counter for the amount of self msd-samples per molecule
static VECTOR *****MsdOrderNOnsager;             // the Onsager-msd in x,y,z per component
static REAL *****MsdOrderNOnsagerDirAvg;         // the Onsager-msd directionally averaged per component
static REAL ****MsdOrderNOnsagerCount;           // counter for the amount of Onsager msd-samples
static VECTOR ***MsdOrderNTotalOnsager;          // the Onsager-msd in x,y,z for the fluid
static REAL ***MsdOrderNTotalOnsagerCount;       // the Onsager-msd directionally for the fluid
static REAL ***MsdOrderNTotalOnsagerDirAvg;      // counter for the amount of Onsager fluid msd-samples
static REAL ****MsdOrderNCountPerSiteType;       // counter for the amount of self msd-samples per site-type
static VECTOR ****MsdOrderNPerSiteType;          // the self-msd in x,y,z per site-type
static REAL ****MsdOrderNPerSiteTypeDirAvg;      // the self-msd directionally averaged per site-type
static int ****BlockDataSiteTypeMSDOrderN;       // array for all blocks containing all the molecule positions per component

// sampling the velocity autocorrelation function using a modified order-N algorithm
int *ComputeVACFOrderN;                          // whether or not to compute the vacf
int *SampleVACFOrderNEvery;                      // the sample frequency
int *WriteVACFOrderNEvery;                       // write output every 'WriteVACFOrderNEvery' steps
int NumberOfBlockElementsVACFOrderN;             // the number of elements per block
int MaxNumberOfBlocksVACFOrderN;                 // the maxmimum amount of blocks (data beyond this block is ignored)
int ComputeIndividualVACFOrderN;                 // whether or not to compute (self-)vacf's for individual molecules
int ComputeVACFOrderNPerPseudoAtom;              // whether or not to compute (self-)vacf's for (pseudo-)atoms
static int *CountVACFOrderN;                     // counter for the amount of vacf's measured
static int *NumberOfBlocksVACFOrderN;            // the current number of blocks in use
static int **BlockLengthVACFOrderN;              // the current length of the blocks
static VECTOR ****BlockDataVACFOrderN;           // array for all blocks containing all the molecule velocities per component
static VECTOR ****BlockDataVACFOrderNOnsager;    // array for all blocks containing the summed molecule velocities per component
static VECTOR ****VacfOrderN;                    // the self-vacf in x,y,z per component
static REAL ****VacfOrderNDirAvg;                // the vacf directionally averaged per component
static VECTOR ****VacfOrderNPerMolecule;         // the self-vacf in x,y,z per molecule
static REAL ****VacfOrderNPerMoleculeDirAvg;     // the self-vacf directionally averaged per molecule
static REAL ****VacfOrderNCount;                 // counter for the amount of self vacf-samples
static REAL ****VacfOrderNCountPerMolecule;      // counter for the amount of self vacf-samples per molecule
static VECTOR *****VacfOrderNOnsager;            // the Onsager-vacf in x,y,z per component
static REAL *****VacfOrderNOnsagerDirAvg;        // the Onsager-vacf directionally averaged per component
static REAL ****VacfOrderNOnsagerCount;          // counter for the amount of Onsager vacf-samples
static VECTOR ***VacfOrderNTotalOnsager;         // the Onsager-msd in x,y,z for the fluid
static REAL ***VacfOrderNTotalOnsagerCount;      // the Onsager-msd directionally for the fluid
static REAL ***VacfOrderNTotalOnsagerDirAvg;     // counter for the amount of Onsager fluid msd-samples

// sampling of the rotational velocity autocorrelation function using a modified order-N algorithm
int *ComputeRVACFOrderN;                          // whether or not to compute the vacf
int *SampleRVACFOrderNEvery;                      // the sample frequency
int *WriteRVACFOrderNEvery;                       // write output every 'WriteVACFOrderNEvery' steps
int NumberOfBlockElementsRVACFOrderN;             // the number of elements per block
int MaxNumberOfBlocksRVACFOrderN;                 // the maxmimum amount of blocks (data beyond this block is ignored)
int ComputeIndividualRVACFOrderN;                 // whether or not to compute (self-)vacf's for individual molecules
int ComputeRVACFOrderNPerPseudoAtom;              // whether or not to compute (self-)vacf's for (pseudo-)atoms
static int *CountRVACFOrderN;                     // counter for the amount of vacf's measured
static int *NumberOfBlocksRVACFOrderN;            // the current number of blocks in use
static int **BlockLengthRVACFOrderN;              // the current length of the blocks
static VECTOR ****BlockDataRVACFOrderN;           // array for all blocks containing all the molecule velocities per component
static VECTOR ****BlockDataRVACFOrderNOnsager;    // array for all blocks containing the summed molecule velocities per component
static VECTOR ****RvacfOrderN;                    // the self-vacf in x,y,z per component
static REAL ****RvacfOrderNDirAvg;                // the vacf directionally averaged per component
static VECTOR ****RvacfOrderNPerMolecule;         // the self-vacf in x,y,z per molecule
static REAL ****RvacfOrderNPerMoleculeDirAvg;     // the self-vacf directionally averaged per molecule
static REAL ****RvacfOrderNCount;                 // counter for the amount of self vacf-samples
static REAL ****RvacfOrderNCountPerMolecule;      // counter for the amount of self vacf-samples per molecule
static VECTOR *****RvacfOrderNOnsager;            // the Onsager-vacf in x,y,z per component
static REAL *****RvacfOrderNOnsagerDirAvg;        // the Onsager-vacf directionally averaged per component
static REAL ****RvacfOrderNOnsagerCount;         // counter for the amount of Onsager vacf-samples
static VECTOR ***RvacfOrderNTotalOnsager;         // the Onsager-msd in x,y,z for the fluid
static REAL ***RvacfOrderNTotalOnsagerCount;      // the Onsager-msd directionally for the fluid
static REAL ***RvacfOrderNTotalOnsagerDirAvg;     // counter for the amount of Onsager fluid msd-samples

// sampling of the molecular orientation autocorrelation function using a modified order-N algorithm
int *ComputeMolecularOrientationOrderN;                          // whether or not to compute the vacf
int *SampleMolecularOrientationOrderNEvery;                      // the sample frequency
int *WriteMolecularOrientationOrderNEvery;                       // write output every 'CountMolecularOrientationOrderN' steps
int MolecularOrientationType;
VECTOR MolecularOrientationVector;
int MolecularOrientationGroup;
int NumberOfBlockElementsMolecularOrientationOrderN;             // the number of elements per block
int MaxNumberOfBlocksMolecularOrientationOrderN;                 // the maxmimum amount of blocks (data beyond this block is ignored)
int ComputeIndividualMolecularOrientationOrderN;                 // whether or not to compute (self-)vacf's for individual molecules
int ComputeMolecularOrientationOrderNPerPseudoAtom;              // whether or not to compute (self-)vacf's for (pseudo-)atoms
static int *CountMolecularOrientationOrderN;                     // counter for the amount of vacf's measured
static int *NumberOfBlocksMolecularOrientationOrderN;            // the current number of blocks in use
static int **BlockLengthMolecularOrientationOrderN;              // the current length of the blocks
static VECTOR ****BlockDataMolecularOrientationOrderN;           // array for all blocks containing all the molecule velocities per component
static VECTOR ****MolecularOrientationOrderN;                    // the self-vacf in x,y,z per component
static REAL ****MolecularOrientationOrderNDirAvg;                // the vacf directionally averaged per component
static REAL ****MolecularOrientationOrderNCount;                 // counter for the amount of self vacf-samples

// sampling of the bond orientation autocorrelation function using a modified order-N algorithm
int *ComputeBondOrientationOrderN;                          // whether or not to compute the vacf
int *SampleBondOrientationOrderNEvery;                      // the sample frequency
int *WriteBondOrientationOrderNEvery;                       // write output every 'WriteBondOrientationOrderNEvery' steps
int NumberOfBlockElementsBondOrientationOrderN;             // the number of elements per block
int MaxNumberOfBlocksBondOrientationOrderN;                 // the maxmimum amount of blocks (data beyond this block is ignored)
int ComputeIndividualBondOrientationOrderN;                 // whether or not to compute (self-)oacf's for individual molecules
int ComputeBondOrientationOrderNPerPseudoAtom;              // whether or not to compute (self-)oacf's for (pseudo-)atoms
static int **CountBondOrientationOrderN;                    // counter for the amount of oacf's measured
static int **NumberOfBlocksBondOrientationOrderN;           // the current number of blocks in use
static int ***BlockLengthBondOrientationOrderN;             // the current length of the blocks
static VECTOR ******BlockDataBondOrientationOrderN;         // array for all blocks containing all the molecule velocities per component
static VECTOR *****BondOrientationOrderN;                   // the self-vacf in x,y,z per component
static REAL *****BondOrientationOrderNDirAvg;               // the oacf directionally averaged per component
static REAL *****BondOrientationOrderNCount;                // counter for the amount of self oacf-samples

int **NumberOfOrientationFrameworkBonds;
char (***OrientationFrameworkBonds)[2][256];
int ***OrientationFrameworkBondTypes;

int ***NumberOfOrientationFrameworkBondPairs;
PAIR ****OrientationFrameworkBondPairs;
VECTOR ****BondOrientationAngleDistributionFunction;
int *BondOrientationAngleHistogramSize;

// sampling the mean-square displacement function using a conventional algorithm
int *ComputeMSD;                                     // whether or not to compute the vacf
int *SampleMSDEvery;                                 // the sample frequency
int *WriteMSDEvery;                                  // write output every 'WriteMSDEvery' steps
int NumberOfBuffersMSD;                              // the number of overlapping buffers
int BufferLengthMSD;                                 // the length of the buffer-arrays
int SampleEveryMSD;                                  // the sample frequency
static VECTOR ***OriginMSD;                          // the stored origins in x,y,z
static VECTOR ***OriginOnsagerMSD;                   // the stored summed origins in x,y,z
static VECTOR ****AcfMSD;                            // the self-msd buffers in x,y,z
static VECTOR *****AcfOnsagerMSD;                    // the Onsager-msd buffers in x,y,z
static REAL ****AcfDirAvgMSD;                        // the self-msd buffers directionally averaged
static REAL *****AcfOnsagerDirAvgMSD;                // the Onsager-msd buffers directionally averaged
static VECTOR ***AccumulatedAcfMSD;                  // the accumulated self-msd's
static VECTOR ****AccumulatedAcfOnsagerMSD;          // the accumulated Onsager-msd's
static REAL ***AccumulatedAcfDirAvgMSD;              // the accumulated self-msd's directionally averaged
static REAL ****AccumulatedAcfOnsagerDirAvgMSD;      // the accumulated Onager-msd's directionally averaged
static int **CountMSD;                               // counter for the msd per component
static int *CountAccumulatedMSD;                     // counter for the accumulated msd's

// sampling the velocity autocorrelation function using a conventional algorithm
int *ComputeVACF;                                    // whether or not to compute the vacf
int *SampleVACFEvery;                                // the sample frequency
int *WriteVACFEvery;                                 // write output every 'WriteVACFEvery' steps
int NumberOfBuffersVACF;                             // the number of overlapping buffers
int BufferLengthVACF;                                // the length of the buffer-arrays
static VECTOR ***OriginVACF;                         // the stored origins in x,y,z
static VECTOR ***OriginOnsagerVACF;                  // the stored summed origins in x,y,z
static VECTOR ****AcfVACF;                           // the self-vacf buffers in x,y,z
static VECTOR *****AcfOnsagerVACF;                   // the Onsager-vacf buffers in x,y,z
static REAL ****AcfDirAvgVACF;                       // the self-vacf buffers directionally averaged
static REAL *****AcfOnsagerDirAvgVACF;               // the Onsager-vacf buffers directionally averaged
static VECTOR ***AccumulatedAcfVACF;                 // the accumulated self-vacf's
static VECTOR ****AccumulatedAcfOnsagerVACF;         // the accumulated Onsager-vacf's
static REAL ***AccumulatedAcfDirAvgVACF;             // the accumulated self-vacf's directionally averaged
static REAL ****AccumulatedAcfOnsagerDirAvgVACF;     // the accumulated Onager-vacf's directionally averaged
static int **CountVACF;                              // counter for the vacf per component
static int *CountAccumulatedVACF;                    // counter for the accumulated vacf's

// sampling the 3D histograms of position (i.e. 3D free energy)
// Calculates a free-energy profile in 3d and outputs in VTK-format
// Used for visualization of a zeolite.
int *ComputeDensityProfile3DVTKGrid;
static REAL ***DensityProfile3D;
static REAL ***COMDensityProfile3D;
int *WriteDensityProfile3DVTKGridEvery;
INT_VECTOR3 DensityProfile3DVTKGridPoints;

// samples the cation sites and adsorption sites
int *ComputeCationAndAdsorptionSites;
int *WriteCationAndAdsorptionSitesEvery;

// samples initial configurations for the tranmission coefficient (dcTST)
int *WritedcTSTSnapShotsToFile;
int *WritedcTSTSnapShotsEvery;

// samples the principle moments of inertia
int ComputePrincipleMomentsOfInertia;

int *ComputeMolecularPressure;

/*********************************************************************************************************
 * Name       | SampleRadialDistributionFunction                                                         *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Sampling the radial distribution function (RDF)                                          *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleRadialDistributionFunction(int Switch)
{
  int f1,f2;
  int i,j,k,l;
  int TypeA,TypeB;
  VECTOR posA,posB,dr;
  REAL r,deltaR;
  FILE *FilePtr;
  REAL normalization;
  char buffer[256];

  switch(Switch)
  {
    case ALLOCATE:
      // allocate mememory for the rdf's, a 4D array
      RadialDistributionFunction=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      RadialDistributionFunctionWithFramework=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      CountRDF=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeRDF[i])
        {
          RadialDistributionFunctionWithFramework[i]=(REAL**)calloc(NumberOfPseudoAtoms,sizeof(REAL*));
          RadialDistributionFunction[i]=(REAL***)calloc(NumberOfPseudoAtoms,sizeof(REAL**));
          for(j=0;j<NumberOfPseudoAtoms;j++)
          {
            RadialDistributionFunctionWithFramework[i][j]=(REAL*)calloc(RDFHistogramSize[i],sizeof(REAL));
            RadialDistributionFunction[i][j]=(REAL**)calloc(NumberOfPseudoAtoms,sizeof(REAL*));
            for(k=0;k<NumberOfPseudoAtoms;k++)
              RadialDistributionFunction[i][j][k]=(REAL*)calloc(RDFHistogramSize[i],sizeof(REAL));
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      // return if the rdf does not has to be calculated for this system
      if(!ComputeRDF[CurrentSystem]) return;

      deltaR=RDFRange[CurrentSystem]/RDFHistogramSize[CurrentSystem];

      // contributions from intra-frameworks
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
        {
          for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
          {
            TypeA=Framework[CurrentSystem].Atoms[f1][i].Type;
            if(PseudoAtoms[TypeA].PrintToPDB)
            {
              posA=Framework[CurrentSystem].Atoms[f1][i].Position;

              for(j=i+1;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
              {
                if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],0))
                {
                  TypeB=Framework[CurrentSystem].Atoms[f1][j].Type;
                  if(PseudoAtoms[TypeB].PrintToPDB)
                  {
                    posB=Framework[CurrentSystem].Atoms[f1][j].Position;

                    dr.x=posA.x-posB.x;
                    dr.y=posA.y-posB.y;
                    dr.z=posA.z-posB.z;
                    dr=ApplyBoundaryCondition(dr);
                    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                    if(r<RDFRange[CurrentSystem])
                    {
                      RadialDistributionFunction[CurrentSystem][TypeA][TypeB][(int)(r/deltaR)]+=1.0;
                      RadialDistributionFunction[CurrentSystem][TypeB][TypeA][(int)(r/deltaR)]+=1.0;
                    }
                  }
                }
              }
            }
          }
        }
      }

      // contributions from interactions between the frameworks
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(f2=f1+1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
        {
          for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
          {
            TypeA=Framework[CurrentSystem].Atoms[f1][i].Type;
            if(PseudoAtoms[TypeA].PrintToPDB)
            {
              posA=Framework[CurrentSystem].Atoms[f1][i].Position;

              for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
              {
                TypeB=Framework[CurrentSystem].Atoms[f2][j].Type;
                if(PseudoAtoms[TypeB].PrintToPDB)
                {
                  posB=Framework[CurrentSystem].Atoms[f2][j].Position;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                  if(r<RDFRange[CurrentSystem])
                  {
                    RadialDistributionFunction[CurrentSystem][TypeA][TypeB][(int)(r/deltaR)]+=1.0;
                    RadialDistributionFunction[CurrentSystem][TypeB][TypeA][(int)(r/deltaR)]+=1.0;
                  }
                }
              }
            }
          }
        }
      }


      // outer loop over adsorbates
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
        {
          TypeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
          if(PseudoAtoms[TypeA].PrintToPDB)
          {
            posA=Adsorbates[CurrentSystem][i].Atoms[k].Position;

            // inner loop over framework atoms
            for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
            {
              for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
              {
                TypeB=Framework[CurrentSystem].Atoms[f1][j].Type;
                if(PseudoAtoms[TypeB].PrintToPDB)
                {
                  posB=Framework[CurrentSystem].Atoms[f1][j].Position;
                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                  if(r<RDFRange[CurrentSystem])
                  {
                    RadialDistributionFunction[CurrentSystem][TypeA][TypeB][(int)(r/deltaR)]+=1.0;
                    RadialDistributionFunction[CurrentSystem][TypeB][TypeA][(int)(r/deltaR)]+=1.0;

                    RadialDistributionFunctionWithFramework[CurrentSystem][TypeA][(int)(r/deltaR)]+=1.0;
                  }
                }
              }
            }

            // inner loop over adsorbate molecules
            for(j=i+1;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
            {
              for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
              {
                TypeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
                if(PseudoAtoms[TypeB].PrintToPDB)
                {
                  posB=Adsorbates[CurrentSystem][j].Atoms[l].Position;
                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

                  if(r<RDFRange[CurrentSystem])
                  {
                    RadialDistributionFunction[CurrentSystem][TypeA][TypeB][(int)(r/deltaR)]+=1.0;
                    RadialDistributionFunction[CurrentSystem][TypeB][TypeA][(int)(r/deltaR)]+=1.0;
                  }
                }
              }
            }

            // inner loop over cation molecules
            for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
            {
              for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
              {
                TypeB=Cations[CurrentSystem][j].Atoms[l].Type;
                if(PseudoAtoms[TypeB].PrintToPDB)
                {
                  posB=Cations[CurrentSystem][j].Atoms[l].Position;
                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

                  if(r<RDFRange[CurrentSystem])
                  {
                    RadialDistributionFunction[CurrentSystem][TypeA][TypeB][(int)(r/deltaR)]+=1.0;
                    RadialDistributionFunction[CurrentSystem][TypeB][TypeA][(int)(r/deltaR)]+=1.0;
                  }
                }
              }
            }
          }
        }
      }

      // outer loop over cations
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
        {
          TypeA=Cations[CurrentSystem][i].Atoms[k].Type;
          if(PseudoAtoms[TypeA].PrintToPDB)
          {
            posA=Cations[CurrentSystem][i].Atoms[k].Position;

            // inner loop over framework atoms
            for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
            {
              for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
              {
                TypeB=Framework[CurrentSystem].Atoms[f1][j].Type;
                if(PseudoAtoms[TypeB].PrintToPDB)
                {
                  posB=Framework[CurrentSystem].Atoms[f1][j].Position;
                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                  if(r<RDFRange[CurrentSystem])
                  {
                    RadialDistributionFunction[CurrentSystem][TypeA][TypeB][(int)(r/deltaR)]+=1.0;
                    RadialDistributionFunction[CurrentSystem][TypeB][TypeA][(int)(r/deltaR)]+=1.0;

                    RadialDistributionFunctionWithFramework[CurrentSystem][TypeA][(int)(r/deltaR)]+=1.0;
                  }
                }
              }
            }

            // inner loop over cation molecules
            for(j=i+1;j<NumberOfCationMolecules[CurrentSystem];j++)
            {
              for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
              {
                TypeB=Cations[CurrentSystem][j].Atoms[l].Type;
                if(PseudoAtoms[TypeB].PrintToPDB)
                {
                  posB=Cations[CurrentSystem][j].Atoms[l].Position;
                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

                  if(r<RDFRange[CurrentSystem])
                  {
                    RadialDistributionFunction[CurrentSystem][TypeA][TypeB][(int)(r/deltaR)]+=1.0;
                    RadialDistributionFunction[CurrentSystem][TypeB][TypeA][(int)(r/deltaR)]+=1.0;
                  }
                }
              }
            }
          }
        }
      }
      CountRDF[CurrentSystem]+=2.0;
      break;
    case PRINT:
      // return if the rdf does not has to be calculated for this system
      if((!ComputeRDF[CurrentSystem])||(CurrentCycle%WriteRDFEvery[CurrentSystem]!=0)) return;

      // make the output directory
      mkdir("RadialDistributionFunctions",S_IRWXU);

      // make the system directory in the output directory
      sprintf(buffer,"RadialDistributionFunctions/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      deltaR=RDFRange[CurrentSystem]/RDFHistogramSize[CurrentSystem];

      for(i=0;i<NumberOfPseudoAtoms;i++)
      {
        if(PseudoAtoms[i].PrintToPDB)
        {
          normalization=Volume[CurrentSystem]/(2.0*M_PI*pow(deltaR,(REAL)3.0)*
                  (REAL)(NumberOfPseudoAtomsType[CurrentSystem][i]*Framework[CurrentSystem].TotalNumberOfAtoms)*CountRDF[CurrentSystem]);

          if(NumberOfPseudoAtomsType[CurrentSystem][i]>0.0)
          {
            sprintf(buffer,"RadialDistributionFunctions/System_%d/RDF_Framework_%s%s.dat",
                  CurrentSystem,PseudoAtoms[i].Name,FileNameAppend);
            FilePtr=fopen(buffer,"w");
            fprintf(FilePtr,"# column 1: index\n");
            fprintf(FilePtr,"# column 2: distance [A]\n");
            fprintf(FilePtr,"# column 3: RDF histogram\n");
            fprintf(FilePtr,"# column 4: unnormalized distance histogram\n");
            for(k=0;k<RDFHistogramSize[CurrentSystem];k++)
            {
              fprintf(FilePtr,"%d %lf %lf %lf\n",
                k,(double)((k+0.5)*deltaR),(double)(RadialDistributionFunctionWithFramework[CurrentSystem][i][k]*normalization/SQR(k+0.5)),
                                           (double)(RadialDistributionFunctionWithFramework[CurrentSystem][i][k]/CountRDF[CurrentSystem]));
            }
            fclose(FilePtr);
          }
          

          for(j=i;j<NumberOfPseudoAtoms;j++)
          {
            if(PseudoAtoms[j].PrintToPDB)
            {
              normalization=Volume[CurrentSystem]/(2.0*M_PI*pow(deltaR,(REAL)3.0)*
                  (REAL)((NumberOfPseudoAtomsType[CurrentSystem][i]-(i==j?1.0:0.0))*NumberOfPseudoAtomsType[CurrentSystem][j])*CountRDF[CurrentSystem]);

              if((NumberOfPseudoAtomsType[CurrentSystem][i]>0.0)&&(NumberOfPseudoAtomsType[CurrentSystem][j]>0.0))
              {
                sprintf(buffer,"RadialDistributionFunctions/System_%d/RDF_%s_%s%s.dat",
                      CurrentSystem,PseudoAtoms[i].Name,PseudoAtoms[j].Name,FileNameAppend);
                FilePtr=fopen(buffer,"w");
                fprintf(FilePtr,"# column 1: index\n");
                fprintf(FilePtr,"# column 2: distance [A]\n");
                fprintf(FilePtr,"# column 3: RDF histogram\n");
                fprintf(FilePtr,"# column 4: unnormalized distance histogram\n");
                for(k=0;k<RDFHistogramSize[CurrentSystem];k++)
                {
                  fprintf(FilePtr,"%d %lf %lf %lf\n",
                    k,(double)((k+0.5)*deltaR),(double)(RadialDistributionFunction[CurrentSystem][i][j][k]*normalization/SQR(k+0.5)),
                                               (double)(RadialDistributionFunction[CurrentSystem][i][j][k]/CountRDF[CurrentSystem]));
                }
                fclose(FilePtr);
              }
            }
          }
        }
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeRDF[i])
        {
          for(j=0;j<NumberOfPseudoAtoms;j++)
          {
            for(k=0;k<NumberOfPseudoAtoms;k++)
              free(RadialDistributionFunction[i][j][k]);
            free(RadialDistributionFunction[i][j]);
          }
          free(RadialDistributionFunction[i]);
        }
      }
      free(RadialDistributionFunction);
      free(CountRDF);
      break;
  }
}

/*********************************************************************************************************
 * Name       | SampleProjectedLengthsDistributionFunction                                               *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Sampling the distribution function of projected lengths                                  *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleProjectedLengthsDistributionFunction(int Switch)
{
  int i,j,k;
  int TypeA,Type;
  VECTOR posA;
  REAL deltaR;
  FILE *FilePtr;
  REAL normalization;
  char buffer[256];
  VECTOR min,max;

  switch(Switch)
  {
    case ALLOCATE:
      // allocate mememory for the storage, a 3D array
      ProjectedLengthsDistributionFunction=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      ProjectedLengthsAverage=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
      CountProjectedLengths=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeProjectedLengths[i])
        {
          ProjectedLengthsDistributionFunction[i]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
          ProjectedLengthsAverage[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
          CountProjectedLengths[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
          for(j=0;j<NumberOfComponents;j++)
          {
            ProjectedLengthsDistributionFunction[i][j]=(VECTOR*)calloc(ProjectedLengthsHistogramSize[i],sizeof(VECTOR));
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      // return if the rdf does not has to be calculated for this system
      if(!ComputeProjectedLengths[CurrentSystem]) return;

      deltaR=ProjectedLengthsRange[CurrentSystem]/ProjectedLengthsHistogramSize[CurrentSystem];

      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        Type=Adsorbates[CurrentSystem][i].Type;

        min.x=DBL_MAX; max.x=-DBL_MAX;
        min.y=DBL_MAX; max.y=-DBL_MAX;
        min.z=DBL_MAX; max.z=-DBL_MAX;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          TypeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          if(PseudoAtoms[TypeA].PrintToPDB)
          {
            posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
            if(posA.x<min.x) min.x=posA.x;
            if(posA.x>max.x) max.x=posA.x;
            if(posA.y<min.y) min.y=posA.y;
            if(posA.y>max.y) max.y=posA.y;
            if(posA.z<min.z) min.z=posA.z;
            if(posA.z>max.z) max.z=posA.z;
          }
        }

        ProjectedLengthsDistributionFunction[CurrentSystem][Type][(int)(fabs(max.x-min.x)/deltaR)].x+=1.0;
        ProjectedLengthsDistributionFunction[CurrentSystem][Type][(int)(fabs(max.y-min.y)/deltaR)].y+=1.0;
        ProjectedLengthsDistributionFunction[CurrentSystem][Type][(int)(fabs(max.z-min.z)/deltaR)].z+=1.0;
        ProjectedLengthsAverage[CurrentSystem][Type].x+=fabs(max.x-min.x);
        ProjectedLengthsAverage[CurrentSystem][Type].y+=fabs(max.y-min.y);
        ProjectedLengthsAverage[CurrentSystem][Type].z+=fabs(max.z-min.z);
        CountProjectedLengths[CurrentSystem][Type]++;
      }
      break;
    case PRINT:
      // return if the rdf does not has to be calculated for this system
      if((!ComputeProjectedLengths[CurrentSystem])||(CurrentCycle%WriteProjectedLengthsEvery[CurrentSystem]!=0)) return;

      // make the output directory
      mkdir("ProjectedLengthsDistributionFunctions",S_IRWXU);

      // make the system directory in the output directory
      sprintf(buffer,"ProjectedLengthsDistributionFunctions/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      deltaR=ProjectedLengthsRange[CurrentSystem]/ProjectedLengthsHistogramSize[CurrentSystem];

      for(i=0;i<NumberOfComponents;i++)
      {
        sprintf(buffer,"ProjectedLengthsDistributionFunctions/System_%d/ProjectedLength_%s%s.dat",
              CurrentSystem,Components[i].Name,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        normalization=CountProjectedLengths[CurrentSystem][i];

        fprintf(FilePtr,"# average x: %18.10f\n",ProjectedLengthsAverage[CurrentSystem][i].x/normalization);
        fprintf(FilePtr,"# average y: %18.10f\n",ProjectedLengthsAverage[CurrentSystem][i].y/normalization);
        fprintf(FilePtr,"# average z: %18.10f\n",ProjectedLengthsAverage[CurrentSystem][i].z/normalization);

        for(k=0;k<ProjectedLengthsHistogramSize[CurrentSystem];k++)
        {
          fprintf(FilePtr,"%d %lf %lf %lf %lf\n",
            k,
            (double)((k+0.5)*deltaR),
            (double)(ProjectedLengthsDistributionFunction[CurrentSystem][i][k].x/normalization),
            (double)(ProjectedLengthsDistributionFunction[CurrentSystem][i][k].y/normalization),
            (double)(ProjectedLengthsDistributionFunction[CurrentSystem][i][k].z/normalization));

        }
        fclose(FilePtr);
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeProjectedLengths[i])
        {
          for(j=0;j<NumberOfComponents;j++)
            free(ProjectedLengthsDistributionFunction[i][j]);
          free(ProjectedLengthsDistributionFunction[i]);
          free(ProjectedLengthsAverage[i]);
          free(CountProjectedLengths[i]);
        }
      }
      free(ProjectedLengthsDistributionFunction);
      free(ProjectedLengthsAverage);
      free(CountProjectedLengths);
      break;
  }
}

/*********************************************************************************************************
 * Name       | SampleProjectedAnglesDistributionFunction                                                *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Sampling the distribution function of projected angles                                   *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleProjectedAnglesDistributionFunction(int Switch)
{
  int i,j,k;
  int Type;
  VECTOR posA,posB,posC;
  REAL deltaR;
  FILE *FilePtr;
  VECTOR normalization;
  REAL length;
  char buffer[256];
  VECTOR v,w;
  int A,B,C;
  VECTOR dr1,dr2;
  REAL angle;

  switch(Switch)
  {
    case ALLOCATE:
      // allocate mememory for the storage, a 3D array
      ProjectedAnglesDistributionFunction=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      CountProjectedAngles=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeProjectedAngles[i])
        {
          ProjectedAnglesDistributionFunction[i]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
          CountProjectedAngles[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
          for(j=0;j<NumberOfComponents;j++)
          {
            ProjectedAnglesDistributionFunction[i][j]=(VECTOR*)calloc(ProjectedAnglesHistogramSize[i],sizeof(VECTOR));
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      // return if the rdf does not has to be calculated for this system
      if(!ComputeProjectedAngles[CurrentSystem]) return;

      deltaR=ProjectedAnglesRange[CurrentSystem]/ProjectedAnglesHistogramSize[CurrentSystem];

      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        Type=Adsorbates[CurrentSystem][i].Type;

        A=Components[Type].orientation.A;
        B=Components[Type].orientation.B;
        C=Components[Type].orientation.C;

        posA=Adsorbates[CurrentSystem][i].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][i].Atoms[B].Position;
        posC=Adsorbates[CurrentSystem][i].Atoms[C].Position;

        dr1.x=posB.x-posA.x;     dr2.x=posB.x-posC.x;
        dr1.y=posB.y-posA.y;     dr2.y=posB.y-posC.y;
        dr1.z=posB.z-posA.z;     dr2.z=posB.z-posC.z;

        v=CrossProduct(dr1,dr2);
        length=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
        v.x/=length;
        v.y/=length;
        v.z/=length;

        w.x=acos(v.x*1.0+v.y*0.0+v.z*0.0);
        w.y=acos(v.x*0.0+v.y*1.0+v.z*0.0);
        w.z=acos(v.x*0.0+v.y*0.0+v.z*1.0);

        if(w.x<0.0) w.x+=2.0*M_PI;
        if(w.y<0.0) w.y+=2.0*M_PI;
        if(w.z<0.0) w.z+=2.0*M_PI;

        ProjectedAnglesDistributionFunction[CurrentSystem][Type][(int)(w.x*RAD2DEG/deltaR)].x+=1.0;
        ProjectedAnglesDistributionFunction[CurrentSystem][Type][(int)(w.y*RAD2DEG/deltaR)].y+=1.0;
        ProjectedAnglesDistributionFunction[CurrentSystem][Type][(int)(w.z*RAD2DEG/deltaR)].z+=1.0;
        CountProjectedAngles[CurrentSystem][Type]++;
      }
      break;
    case PRINT:
      // return if the rdf does not has to be calculated for this system
      if((!ComputeProjectedAngles[CurrentSystem])||(CurrentCycle%WriteProjectedAnglesEvery[CurrentSystem]!=0)) return;

      // make the output directory
      mkdir("ProjectedAnglesDistributionFunctions",S_IRWXU);

      // make the system directory in the output directory
      sprintf(buffer,"ProjectedAnglesDistributionFunctions/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      deltaR=ProjectedAnglesRange[CurrentSystem]/ProjectedAnglesHistogramSize[CurrentSystem];

      for(i=0;i<NumberOfComponents;i++)
      {
        sprintf(buffer,"ProjectedAnglesDistributionFunctions/System_%d/ProjectedAngle_%s%s.dat",
              CurrentSystem,Components[i].Name,FileNameAppend);
        FilePtr=fopen(buffer,"w");

        normalization.x=normalization.y=normalization.z=0.0;
        for(k=0;k<ProjectedAnglesHistogramSize[CurrentSystem];k++)
        {
          normalization.x+=ProjectedAnglesDistributionFunction[CurrentSystem][i][k].x;
          normalization.y+=ProjectedAnglesDistributionFunction[CurrentSystem][i][k].y;
          normalization.z+=ProjectedAnglesDistributionFunction[CurrentSystem][i][k].z;
        }

        for(k=0;k<ProjectedAnglesHistogramSize[CurrentSystem];k++)
        {
          angle=(k+0.5)*deltaR;
          fprintf(FilePtr,"%d %lf %lf %lf %lf %lf %lf %lf\n",
            k,
            (double)angle,
            (double)(ProjectedAnglesDistributionFunction[CurrentSystem][i][k].x*ProjectedAnglesHistogramSize[CurrentSystem]*2.0/(M_PI*normalization.x)),
            (double)(ProjectedAnglesDistributionFunction[CurrentSystem][i][k].y*ProjectedAnglesHistogramSize[CurrentSystem]*2.0/(M_PI*normalization.y)),
            (double)(ProjectedAnglesDistributionFunction[CurrentSystem][i][k].z*ProjectedAnglesHistogramSize[CurrentSystem]*2.0/(M_PI*normalization.z)),
            (double)(ProjectedAnglesDistributionFunction[CurrentSystem][i][k].x*ProjectedAnglesHistogramSize[CurrentSystem]*2.0/(M_PI*normalization.x*sin(angle*DEG2RAD))),
            (double)(ProjectedAnglesDistributionFunction[CurrentSystem][i][k].y*ProjectedAnglesHistogramSize[CurrentSystem]*2.0/(M_PI*normalization.z*sin(angle*DEG2RAD))),
            (double)(ProjectedAnglesDistributionFunction[CurrentSystem][i][k].z*ProjectedAnglesHistogramSize[CurrentSystem]*2.0/(M_PI*normalization.z*sin(angle*DEG2RAD))));


        }
        fclose(FilePtr);
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeProjectedAngles[i])
        {
          for(j=0;j<NumberOfComponents;j++)
            free(ProjectedAnglesDistributionFunction[i][j]);
          free(ProjectedAnglesDistributionFunction[i]);
          free(CountProjectedAngles[i]);
        }
      }
      free(ProjectedAnglesDistributionFunction);
      free(CountProjectedAngles);
      break;
  }
}



/********************************************************************************************************
* Name       | SampleCFCRXMCLambdaHistogram                                                             *
* ----------------------------------------------------------------------------------------------------- *
* Function   | Samples the lambda histograms                                                            *
* Parameters | -                                                                                        *
*********************************************************************************************************/

void SampleCFCRXMCLambdaHistogram(int Switch)
{
  int i,j,k,index;
  REAL norm,r,delta;
  char buffer[256];
  FILE *FilePtr;

  switch(Switch)
  {
    case ALLOCATE:
      CFCRXMCLambdaHistogram=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
         if(ComputeCFCRXMCLambdaHistogram[i])
         {
           CFCRXMCLambdaHistogram[i]=(REAL**)calloc(NumberOfReactions,sizeof(REAL*));
           for(j=0;j<NumberOfReactions;j++)
              CFCRXMCLambdaHistogram[i][j]=(REAL*)calloc(CFCRXMCLambdaHistogramBins[i]+1,sizeof(REAL));
         }
      }
      break;
    case INITIALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
         if(ComputeCFCRXMCLambdaHistogram[i])
         {
           for(j=0;j<NumberOfReactions;j++)
              for(k=0;k<=CFCRXMCLambdaHistogramBins[i];k++)
                 CFCRXMCLambdaHistogram[i][j][k]=0.0;
         }
      }
      break;
    case SAMPLE:
      if(!ComputeCFCRXMCLambdaHistogram[CurrentSystem]) return;

      for(j=0;j<NumberOfReactions;j++)
      {
        for(i=0;i<NumberOfReactions;i++)
        {
           r=CFCRXMCLambda[j][i];
           delta=CFCRXMCLambdaHistogramBins[CurrentSystem];
           index=(int)(r*delta);
           CFCRXMCLambdaHistogram[CurrentSystem][i][index]+=1.0;
        }
      }
      break;
    case PRINT:
      if((!ComputeCFCRXMCLambdaHistogram[CurrentSystem])||(CurrentCycle%WriteCFCRXMCLambdaHistogramEvery[CurrentSystem]!=0)) return;
      mkdir("CFCRXMCLambdaHistograms",S_IRWXU);
      sprintf(buffer,"CFCRXMCLambdaHistograms/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      delta=CFCRXMCLambdaHistogramBins[CurrentSystem];
      for(i=0;i<NumberOfReactions;i++)
      {
        norm=0.0;
        for(k=0;k<=CFCRXMCLambdaHistogramBins[CurrentSystem];k++)
          norm+=CFCRXMCLambdaHistogram[CurrentSystem][i][k];
        norm/=delta;

        sprintf(buffer,"CFCRXMCLambdaHistograms/System_%d/Histogram_Reaction_%d%s.dat",CurrentSystem,i,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        for(k=0;k<=CFCRXMCLambdaHistogramBins[CurrentSystem];k++)
        {
          r=(REAL)k/delta;
          if(CFCRXMCLambdaHistogram[CurrentSystem][i][k]>0.0)
            fprintf(FilePtr,"%g %g\n",(double)r,(double)(CFCRXMCLambdaHistogram[CurrentSystem][i][k]/norm));
        }
        fclose(FilePtr);
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
         if(ComputeCFCRXMCLambdaHistogram[i])
         {
           for(j=0;j<NumberOfReactions;j++)
             free(CFCRXMCLambdaHistogram[i][j]);
           free(CFCRXMCLambdaHistogram[i]);
         }
      }
      free(CFCRXMCLambdaHistogram);
      break;
  }
}



/*********************************************************************************************************
 * Name       | SampleNumberOfMoleculesHistogram                                                         *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the histograms of the number of molecules.                                       *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleNumberOfMoleculesHistogram(int Switch)
{
  int i,j,k,index;
  REAL norm,r,delta;
  char buffer[256];
  FILE *FilePtr;

  switch(Switch)
  {
    case ALLOCATE:
      NumberOfMoleculesHistogram=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeNumberOfMoleculesHistogram[i])
        {
          NumberOfMoleculesHistogram[i]=(REAL**)calloc(NumberOfComponents+1,sizeof(REAL*));
          for(j=0;j<NumberOfComponents+1;j++)
            NumberOfMoleculesHistogram[i][j]=(REAL*)calloc(NumberOfMoleculesHistogramSize[i],sizeof(REAL));
        }
      }
      break;
    case INITIALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeNumberOfMoleculesHistogram[i])
        {
          for(j=0;j<NumberOfComponents+1;j++)
            for(k=0;k<NumberOfMoleculesHistogramSize[i];k++)
              NumberOfMoleculesHistogram[i][j][k]=0.0;
        }
      }
      break;
    case SAMPLE:
      if(!ComputeNumberOfMoleculesHistogram[CurrentSystem]) return;

      delta=NumberOfMoleculesHistogramSize[CurrentSystem]/NumberOfMoleculesRange[CurrentSystem];
      for(i=0;i<NumberOfComponents;i++)
      {
        r=Components[i].NumberOfMolecules[CurrentSystem];
        index=(int)(r*delta);
        if(index>=0&&index<NumberOfMoleculesHistogramSize[CurrentSystem])
          NumberOfMoleculesHistogram[CurrentSystem][i][index]+=1.0;
      }
      break;
    case PRINT:
      if((!ComputeNumberOfMoleculesHistogram[CurrentSystem])||(CurrentCycle%WriteNumberOfMoleculesHistogramEvery[CurrentSystem]!=0)) return;
      mkdir("NumberOfMoleculesHistograms",S_IRWXU);
      sprintf(buffer,"NumberOfMoleculesHistograms/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      delta=NumberOfMoleculesHistogramSize[CurrentSystem]/NumberOfMoleculesRange[CurrentSystem];
      for(i=0;i<NumberOfComponents;i++)
      {
        norm=0.0;
        for(k=0;k<NumberOfMoleculesHistogramSize[CurrentSystem];k++)
          norm+=NumberOfMoleculesHistogram[CurrentSystem][i][k];
        norm/=delta;

        sprintf(buffer,"NumberOfMoleculesHistograms/System_%d/Histogram_%s_%d_Density_%g%s.dat",
                CurrentSystem,Components[i].Name,i,
                (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                FileNameAppend);
        FilePtr=fopen(buffer,"w");
        for(k=0;k<NumberOfMoleculesHistogramSize[CurrentSystem];k++)
        {
          r=(REAL)k/delta;
          if(NumberOfMoleculesHistogram[CurrentSystem][i][k]>0.0)
            fprintf(FilePtr,"%g %g\n",(double)r,(double)(NumberOfMoleculesHistogram[CurrentSystem][i][k]/norm));
        }
        fclose(FilePtr);
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeNumberOfMoleculesHistogram[i])
        {
          for(j=0;j<NumberOfComponents+1;j++)
            free(NumberOfMoleculesHistogram[i][j]);
          free(NumberOfMoleculesHistogram[i]);
        }
      }
      free(NumberOfMoleculesHistogram);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SamplePositionHistogram                                                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the position histograms and free energy profiles in a,b,c directions.            *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SamplePositionHistogram(int Switch)
{
  int i,j,k;
  int index,type;
  REAL F;
  VECTOR pos,s,fcom,drift;
  FILE *FilePtr;
  char buffer[256];
  VECTOR norm;
  VECTOR4 norm4;
  VECTOR dr;
  REAL q;

  switch(Switch)
  {
    case ALLOCATE:
      PositionABCHistogram=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      Position2DDiagonalHistogram=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      Position2DDiagonalHistogram2=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      Position3DDiagonalHistogram=(VECTOR4****)calloc(NumberOfSystems,sizeof(VECTOR4***));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputePositionHistogram[i])
        {
          PositionABCHistogram[i]=(VECTOR***)calloc(NumberOfComponents,sizeof(VECTOR**));
          Position2DDiagonalHistogram[i]=(VECTOR***)calloc(NumberOfComponents,sizeof(VECTOR**));
          Position2DDiagonalHistogram2[i]=(VECTOR***)calloc(NumberOfComponents,sizeof(VECTOR**));
          Position3DDiagonalHistogram[i]=(VECTOR4***)calloc(NumberOfComponents,sizeof(VECTOR4**));
          for(j=0;j<NumberOfComponents;j++)
          {
            PositionABCHistogram[i][j]=(VECTOR**)calloc(Components[j].NumberOfAtoms+2,sizeof(VECTOR*));
            Position2DDiagonalHistogram[i][j]=(VECTOR**)calloc(Components[j].NumberOfAtoms+2,sizeof(VECTOR*));
            Position2DDiagonalHistogram2[i][j]=(VECTOR**)calloc(Components[j].NumberOfAtoms+2,sizeof(VECTOR*));
            Position3DDiagonalHistogram[i][j]=(VECTOR4**)calloc(Components[j].NumberOfAtoms+2,sizeof(VECTOR4*));
            for(k=0;k<Components[j].NumberOfAtoms+2;k++)
            {
              PositionABCHistogram[i][j][k]=(VECTOR*)calloc(PositionHistogramSize[i],sizeof(VECTOR));
              Position2DDiagonalHistogram[i][j][k]=(VECTOR*)calloc(PositionHistogramSize[i],sizeof(VECTOR));
              Position2DDiagonalHistogram2[i][j][k]=(VECTOR*)calloc(PositionHistogramSize[i],sizeof(VECTOR));
              Position3DDiagonalHistogram[i][j][k]=(VECTOR4*)calloc(PositionHistogramSize[i],sizeof(VECTOR4));
            }
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      if(!ComputePositionHistogram[CurrentSystem]) return;

      fcom=GetFrameworkCenterOfMass();
      drift.x=fcom.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x;
      drift.y=fcom.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y;
      drift.z=fcom.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z;

      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        type=Adsorbates[CurrentSystem][i].Type;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          pos=Adsorbates[CurrentSystem][i].Atoms[j].Position;

          //pos.x-=drift.x;
          //pos.y-=drift.y;
          //pos.z-=drift.z;

          s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
          s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
          s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;

          // apply boundary condition
          s.x-=(REAL)NINT(s.x);
          s.y-=(REAL)NINT(s.y);
          s.z-=(REAL)NINT(s.z);

          // 's' the fractional position within the first unit cell
          if(s.x<0.0) s.x+=1.0;
          if(s.y<0.0) s.y+=1.0;
          if(s.z<0.0) s.z+=1.0;

          // map the fractional position to a reaction coordinate
          switch(PositionHistogramMappingType[CurrentSystem])
          {
            case NO_MAPPING:
              break;
            case A_MAPPING:
              q=s.x;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                PositionABCHistogram[CurrentSystem][type][j][index].x+=1.0;
              break;
            case B_MAPPING:
              q=s.y;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                PositionABCHistogram[CurrentSystem][type][j][index].y+=1.0;
              break;
            case C_MAPPING:
              q=s.z;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                PositionABCHistogram[CurrentSystem][type][j][index].z+=1.0;
              break;
            case ABC_MAPPING:
              q=s.x;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                PositionABCHistogram[CurrentSystem][type][j][index].x+=1.0;
              q=s.y;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                PositionABCHistogram[CurrentSystem][type][j][index].y+=1.0;
              q=s.z;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                PositionABCHistogram[CurrentSystem][type][j][index].z+=1.0;
              break;
            case MAP_AB_DIAGONAL:   // e.g. DDR
              dr.x=M_SQRT1_2;
              dr.y=-M_SQRT1_2;
              dr.z=0.0;
              q=((1.0-s.x)*dr.x+(-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_2;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                Position2DDiagonalHistogram[CurrentSystem][type][j][index].x+=1.0;
              break;
            case MAP_AC_DIAGONAL:
              dr.x=M_SQRT1_2;
              dr.y=0.0;
              dr.z=-M_SQRT1_2;
              q=((1.0-s.x)*dr.x+(-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_2;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                Position2DDiagonalHistogram[CurrentSystem][type][j][index].y+=1.0;
              break;
            case MAP_BC_DIAGONAL:
              dr.x=0.0;
              dr.y=M_SQRT1_2;
              dr.z=-M_SQRT1_2;
              q=((-s.x)*dr.x+(1.0-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_2;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                Position2DDiagonalHistogram[CurrentSystem][type][j][index].z+=1.0;
              break;
            case MAP_O_AB_DIAGONAL:
              dr.x=M_SQRT1_2;
              dr.y=M_SQRT1_2;
              dr.z=0.0;
              q=(s.x*dr.x+s.y*dr.y+s.z*dr.z)*M_SQRT1_2;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                Position2DDiagonalHistogram2[CurrentSystem][type][j][index].x+=1.0;
              break;
            case MAP_O_AC_DIAGONAL:
              dr.x=M_SQRT1_2;
              dr.y=0.0;
              dr.z=M_SQRT1_2;
              q=(s.x*dr.x+s.y*dr.y+s.z*dr.z)*M_SQRT1_2;

              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                Position2DDiagonalHistogram2[CurrentSystem][type][j][index].y+=1.0;
              break;
            case MAP_O_BC_DIAGONAL:
              dr.x=0.0;
              dr.y=M_SQRT1_2;
              dr.z=M_SQRT1_2;
              q=(s.x*dr.x+s.y*dr.y+s.z*dr.z)*M_SQRT1_2;

              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                Position2DDiagonalHistogram2[CurrentSystem][type][j][index].z+=1.0;
              break;
            case MAP_A_BC_DIAGONAL:
              dr.x=M_SQRT1_3;
              dr.y=-M_SQRT1_3;
              dr.z=-M_SQRT1_3;
              q=((1.0-s.x)*dr.x+(-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_3;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                Position3DDiagonalHistogram[CurrentSystem][type][j][index].r+=1.0;
              break;
            case MAP_B_AC_DIAGONAL:
              dr.x=-M_SQRT1_3;
              dr.y=M_SQRT1_3;
              dr.z=-M_SQRT1_3;
              q=((-s.x)*dr.x+(1.0-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_3;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                Position3DDiagonalHistogram[CurrentSystem][type][j][index].i+=1.0;
              break;
            case MAP_C_AB_DIAGONAL:
              dr.x=-M_SQRT1_3;
              dr.y=-M_SQRT1_3;
              dr.z=M_SQRT1_3;
              q=((-s.x)*dr.x+(-s.y)*dr.y+(1.0-s.z)*dr.z)*M_SQRT1_3;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                Position3DDiagonalHistogram[CurrentSystem][type][j][index].j+=1.0;
              break;
            case MAP_O_ABC_DIAGONAL:
              dr.x=-M_SQRT1_3;
              dr.y=-M_SQRT1_3;
              dr.z=-M_SQRT1_3;
              q=((-s.x)*dr.x+(-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_3;
              index=(int)(q*(REAL)PositionHistogramSize[CurrentSystem]);
              if(index>=0&&index<PositionHistogramSize[CurrentSystem])
                Position3DDiagonalHistogram[CurrentSystem][type][j][index].k+=1.0;
              break;
          }
        }
      }

      break;
    case PRINT:
      if((!ComputePositionHistogram[CurrentSystem])||(CurrentCycle%WritePositionHistogramEvery[CurrentSystem]!=0)) return;

      for(i=0;i<NumberOfComponents;i++)
      {
        for(j=0;j<Components[i].NumberOfAtoms;j++)
        {
          // mappings onto a,b,c-directions
          // ------------------------------

          norm.x=norm.y=norm.z=0.0;
          for(k=0;k<PositionHistogramSize[CurrentSystem];k++)
          {
            norm.x+=PositionABCHistogram[CurrentSystem][i][j][k].x;
            norm.y+=PositionABCHistogram[CurrentSystem][i][j][k].y;
            norm.z+=PositionABCHistogram[CurrentSystem][i][j][k].z;
          }

          mkdir("Histograms",S_IRWXU);

          sprintf(buffer,"Histograms/System_%d",CurrentSystem);
          mkdir(buffer,S_IRWXU);

          switch(PositionHistogramMappingType[CurrentSystem])
          {
            case NO_MAPPING:
              break;
            case A_MAPPING:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_A",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. a in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].x/norm.x),
                    (double)(-log(PositionABCHistogram[CurrentSystem][i][j][index].x/norm.x)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==A_MAPPING)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].x/norm.x),
                    (double)F,
                    (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].x/norm.x)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_A%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. a in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==A_MAPPING))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].x/norm.x)));
              }
              fclose(FilePtr);
              break;
            case B_MAPPING:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_B",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. b in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].y/norm.y),
                    (double)(-log(PositionABCHistogram[CurrentSystem][i][j][index].y/norm.y)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==B_MAPPING)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].y/norm.y),
                    (double)F,
                    (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].y/norm.y)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_B%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. a in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==B_MAPPING))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].y/norm.y)));
              }
              fclose(FilePtr);
              break;
            case C_MAPPING:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_C",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. c in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].z/norm.z),
                    (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].z/norm.z)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==C_MAPPING)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].z/norm.z),
                    (double)F,
                    (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].z/norm.z)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_C%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. c in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==C_MAPPING))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].z/norm.z)));
              }
              fclose(FilePtr);
              break;
            case ABC_MAPPING:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_A",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. a in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].x/norm.x),
                    (double)(-log(PositionABCHistogram[CurrentSystem][i][j][index].x/norm.x)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==A_MAPPING)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].x/norm.x),
                    (double)F,
                    (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].x/norm.x)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_A%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. a in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==A_MAPPING))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].x/norm.x)));
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_B",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. b in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].y/norm.y),
                    (double)(-log(PositionABCHistogram[CurrentSystem][i][j][index].y/norm.y)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==B_MAPPING)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].y/norm.y),
                    (double)F,
                    (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].y/norm.y)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_B%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. a in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==B_MAPPING))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].y/norm.y)));
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_C",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. c in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].z/norm.z),
                    (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].z/norm.z)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==C_MAPPING)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(PositionABCHistogram[CurrentSystem][i][j][index].z/norm.z),
                    (double)F,
                    (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].z/norm.z)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_C%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. c in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==C_MAPPING))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(PositionABCHistogram[CurrentSystem][i][j][index].z/norm.z)));
              }
              fclose(FilePtr);
              break;
          }

          // mappings onto 2D diagonal-directions
          // ------------------------------------

          norm.x=norm.y=norm.z=0.0;
          for(k=0;k<PositionHistogramSize[CurrentSystem];k++)
          {
            norm.x+=Position2DDiagonalHistogram[CurrentSystem][i][j][k].x;
            norm.y+=Position2DDiagonalHistogram[CurrentSystem][i][j][k].y;
            norm.z+=Position2DDiagonalHistogram[CurrentSystem][i][j][k].z;
          }

          mkdir("Histograms",S_IRWXU);

          sprintf(buffer,"Histograms/System_%d",CurrentSystem);
          mkdir(buffer,S_IRWXU);
          switch(PositionHistogramMappingType[CurrentSystem])
          {
            case MAP_AB_DIAGONAL:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_AB_2D",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D AB-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram[CurrentSystem][i][j][index].x/norm.x),
                    (double)(-log(Position2DDiagonalHistogram[CurrentSystem][i][j][index].x/norm.x)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==MAP_AB_DIAGONAL)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram[CurrentSystem][i][j][index].x/norm.x),
                    (double)F,
                    (double)(F-log(Position2DDiagonalHistogram[CurrentSystem][i][j][index].x/norm.x)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_AB_2D%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D AB-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==MAP_AB_DIAGONAL))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(Position2DDiagonalHistogram[CurrentSystem][i][j][index].x/norm.x)));
              }
              fclose(FilePtr);
              break;
            case MAP_AC_DIAGONAL:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_AC_2D",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D AC-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram[CurrentSystem][i][j][index].y/norm.y),
                    (double)(-log(Position2DDiagonalHistogram[CurrentSystem][i][j][index].y/norm.y)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==MAP_AC_DIAGONAL)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram[CurrentSystem][i][j][index].y/norm.y),
                    (double)F,
                    (double)(F-log(Position2DDiagonalHistogram[CurrentSystem][i][j][index].y/norm.y)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_AC_2D%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D AC-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==MAP_AC_DIAGONAL))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(Position2DDiagonalHistogram[CurrentSystem][i][j][index].y/norm.y)));
              }
              fclose(FilePtr);
              break;
            case MAP_BC_DIAGONAL:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_BC_2D",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D BC-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram[CurrentSystem][i][j][index].z/norm.z),
                    (double)(-log(Position2DDiagonalHistogram[CurrentSystem][i][j][index].z/norm.z)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==MAP_BC_DIAGONAL)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram[CurrentSystem][i][j][index].z/norm.z),
                    (double)F,
                    (double)(F-log(Position2DDiagonalHistogram[CurrentSystem][i][j][index].z/norm.z)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_BC_2D%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D BC-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==MAP_BC_DIAGONAL))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(Position2DDiagonalHistogram[CurrentSystem][i][j][index].z/norm.z)));
              }
              fclose(FilePtr);
              break;
          }


          norm.x=norm.y=norm.z=0.0;
          for(k=0;k<PositionHistogramSize[CurrentSystem];k++)
          {
            norm.x+=Position2DDiagonalHistogram2[CurrentSystem][i][j][k].x;
            norm.y+=Position2DDiagonalHistogram2[CurrentSystem][i][j][k].y;
            norm.z+=Position2DDiagonalHistogram2[CurrentSystem][i][j][k].z;
          }

          mkdir("Histograms",S_IRWXU);

          sprintf(buffer,"Histograms/System_%d",CurrentSystem);
          mkdir(buffer,S_IRWXU);
          switch(PositionHistogramMappingType[CurrentSystem])
          {
            case MAP_O_AB_DIAGONAL:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_O_AB_2D",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D O-A+B-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].x/norm.x),
                    (double)(-log(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].x/norm.x)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==MAP_AB_DIAGONAL)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].x/norm.x),
                    (double)F,
                    (double)(F-log(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].x/norm.x)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_O_AB_2D%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D O-A+B-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==MAP_AB_DIAGONAL))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].x/norm.x)));
              }
              fclose(FilePtr);
              break;
            case MAP_O_AC_DIAGONAL:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_O_AC_2D",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D O-A+C-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].y/norm.y),
                    (double)(-log(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].y/norm.y)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==MAP_AC_DIAGONAL)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].y/norm.y),
                    (double)F,
                    (double)(F-log(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].y/norm.y)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_O_AC_2D%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D O-A+C-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==MAP_AC_DIAGONAL))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].y/norm.y)));
              }
              fclose(FilePtr);
              break;
            case MAP_O_BC_DIAGONAL:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_O_BC_2D",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D O-B+C-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].z/norm.z),
                    (double)(-log(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].z/norm.z)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==MAP_BC_DIAGONAL)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].z/norm.z),
                    (double)F,
                    (double)(F-log(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].z/norm.z)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_O_BC_2D%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D O-B+C-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==MAP_BC_DIAGONAL))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(Position2DDiagonalHistogram2[CurrentSystem][i][j][index].z/norm.z)));
              }
              fclose(FilePtr);
              break;
          }


          // mappings onto 3D diagonal-directions
          // ------------------------------------

          norm4.r=norm4.i=norm4.j=norm4.k=0.0;
          for(k=0;k<PositionHistogramSize[CurrentSystem];k++)
          {
            norm4.r+=Position3DDiagonalHistogram[CurrentSystem][i][j][k].r;
            norm4.i+=Position3DDiagonalHistogram[CurrentSystem][i][j][k].i;
            norm4.j+=Position3DDiagonalHistogram[CurrentSystem][i][j][k].j;
            norm4.k+=Position3DDiagonalHistogram[CurrentSystem][i][j][k].k;
          }

          mkdir("Histograms",S_IRWXU);

          sprintf(buffer,"Histograms/System_%d",CurrentSystem);
          mkdir(buffer,S_IRWXU);
          switch(PositionHistogramMappingType[CurrentSystem])
          {
            case MAP_A_BC_DIAGONAL:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_A_BC_3D",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D AB-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position3DDiagonalHistogram[CurrentSystem][i][j][index].r/norm4.r),
                    (double)(-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].r/norm4.r)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==MAP_A_BC_DIAGONAL)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position3DDiagonalHistogram[CurrentSystem][i][j][index].r/norm4.r),
                    (double)F,
                    (double)(F-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].r/norm4.r)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_A_BC_3D%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D AB-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==MAP_A_BC_DIAGONAL))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].r/norm4.r)));
              }
              fclose(FilePtr);
              break;
            case MAP_B_AC_DIAGONAL:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_B_AC_3D",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D AC-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position3DDiagonalHistogram[CurrentSystem][i][j][index].i/norm4.i),
                    (double)(-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].i/norm4.i)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==MAP_B_AC_DIAGONAL)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position3DDiagonalHistogram[CurrentSystem][i][j][index].i/norm4.i),
                    (double)F,
                    (double)(F-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].i/norm4.i)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_B_AC_3D%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D AC-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==MAP_B_AC_DIAGONAL))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].i/norm4.i)));
              }
              fclose(FilePtr);
              break;
            case MAP_C_AB_DIAGONAL:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_C_AB_3D",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D BC-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position3DDiagonalHistogram[CurrentSystem][i][j][index].j/norm4.j),
                    (double)(-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].j/norm4.j)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==MAP_C_AB_DIAGONAL)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position3DDiagonalHistogram[CurrentSystem][i][j][index].j/norm4.j),
                    (double)F,
                    (double)(F-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].j/norm4.j)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_C_AB_3D%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D BC-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==MAP_C_AB_DIAGONAL))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].j/norm4.j)));
              }
              fclose(FilePtr);
              break;
            case MAP_O_ABC_DIAGONAL:
              sprintf(buffer,"Histograms/System_%d/Histogram_%s_%d_Bead%d%s.dat_O_ABC_3D",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# Biased %s\n",(Components[i].Biased!=NO_BIASING)?"yes":"no");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D BC-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: histogram of the q-position\n");
              if(Components[i].Biased!=NO_BIASING)
              {
                fprintf(FilePtr,"# column 3: the value of the biasing spline at q (optional)\n");
                fprintf(FilePtr,"# column 4: the free energy profile [k_BT]\n");
              }
              else
                fprintf(FilePtr,"# column 3: the free energy profile [k_BT]\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                if(Components[i].Biased==NO_BIASING)
                {
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position3DDiagonalHistogram[CurrentSystem][i][j][index].k/norm4.k),
                    (double)(-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].k/norm4.k)));
                }
                else
                {
                  F=0.0;
                  if(Components[i].BiasingDirection==MAP_O_ABC_DIAGONAL)
                    F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                  fprintf(FilePtr,"%g %g %g %g\n",
                    (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                    (double)(Position3DDiagonalHistogram[CurrentSystem][i][j][index].k/norm4.k),
                    (double)F,
                    (double)(F-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].k/norm4.k)));
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"Histograms/System_%d/FreeEnergy_%s_%d_Bead%d_O_ABC_3D%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
              FilePtr=fopen(buffer,"w");
              fprintf(FilePtr,"# column 1: the reaction coordinate q, i.e. the 2D BC-diagonal in fractional units\n");
              fprintf(FilePtr,"# column 2: the free energy profile [k_BT], (including a possible biasing free energy)\n");
              fprintf(FilePtr,"# column 3: error (not computed yet, fixed at 0.2)\n");
              for(k=0;k<=PositionHistogramSize[CurrentSystem];k++)
              {
                index=k%PositionHistogramSize[CurrentSystem];

                F=0.0;
                if((Components[i].Biased!=NO_BIASING)&&(Components[i].BiasingDirection==MAP_O_ABC_DIAGONAL))
                  F=BiasingPotentialUmbrellaQ(i,(REAL)k/(REAL)PositionHistogramSize[CurrentSystem]);

                fprintf(FilePtr,"%g %g 0.2\n",
                  (double)k/(REAL)PositionHistogramSize[CurrentSystem],
                  (double)(F-log(Position3DDiagonalHistogram[CurrentSystem][i][j][index].k/norm4.k)));
              }
              fclose(FilePtr);
              break;
          }
        }
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputePositionHistogram[i])
        {
          for(j=0;j<NumberOfComponents;j++)
          {
            for(k=0;k<Components[j].NumberOfAtoms+2;k++)
            {
              free(PositionABCHistogram[i][j][k]);
              free(Position2DDiagonalHistogram[i][j][k]);
              free(Position2DDiagonalHistogram2[i][j][k]);
              free(Position3DDiagonalHistogram[i][j][k]);
            }
            free(PositionABCHistogram[i][j]);
            free(Position2DDiagonalHistogram[i][j]);
            free(Position2DDiagonalHistogram2[i][j]);
            free(Position3DDiagonalHistogram[i][j]);
          }
          free(PositionABCHistogram[i]);
          free(Position2DDiagonalHistogram[i]);
          free(Position2DDiagonalHistogram2[i]);
          free(Position3DDiagonalHistogram[i]);
        }
      }
      free(PositionABCHistogram);
      free(Position2DDiagonalHistogram);
      free(Position2DDiagonalHistogram2);
      free(Position3DDiagonalHistogram);
      break;
    default:
      break;
  }
}

POINT GenerateRandomPointInCartesianUnitCellSpace(void)
{
  VECTOR SizeGrid;
  VECTOR ShiftGrid;
  POINT P,s;

  SizeGrid.x=SizeGrid.y=SizeGrid.z=0.0;
  ShiftGrid.x=ShiftGrid.y=ShiftGrid.z=0.0;

  SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].ax);
  SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].ay);
  SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].az);
  if(UnitCellBox[CurrentSystem].ax<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].ax;
  if(UnitCellBox[CurrentSystem].ay<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].ay;
  if(UnitCellBox[CurrentSystem].az<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].az;

  SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].bx);
  SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].by);
  SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].bz);
  if(UnitCellBox[CurrentSystem].bx<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].bx;
  if(UnitCellBox[CurrentSystem].by<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].by;
  if(UnitCellBox[CurrentSystem].bz<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].bz;

  SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].cx);
  SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].cy);
  SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].cz);
  if(UnitCellBox[CurrentSystem].cx<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].cx;
  if(UnitCellBox[CurrentSystem].cy<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].cy;
  if(UnitCellBox[CurrentSystem].cz<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].cz;

  do
  {
    P.x=RandomNumber()*SizeGrid.x+ShiftGrid.x;
    P.y=RandomNumber()*SizeGrid.y+ShiftGrid.y;
    P.z=RandomNumber()*SizeGrid.z+ShiftGrid.z;

    s.x=InverseUnitCellBox[CurrentSystem].ax*P.x+InverseUnitCellBox[CurrentSystem].bx*P.y+InverseUnitCellBox[CurrentSystem].cx*P.z;
    s.y=InverseUnitCellBox[CurrentSystem].ay*P.x+InverseUnitCellBox[CurrentSystem].by*P.y+InverseUnitCellBox[CurrentSystem].cy*P.z;
    s.z=InverseUnitCellBox[CurrentSystem].az*P.x+InverseUnitCellBox[CurrentSystem].bz*P.y+InverseUnitCellBox[CurrentSystem].cz*P.z;
  } while(!((s.x>=0.0)&&(s.x<1.0)&&(s.y>=0.0)&&(s.y<1.0)&&(s.z>=0.0)&&(s.z<1.0)));

  return s;
}


/*********************************************************************************************************
 * Name       | SampleFreeEnergyProfile                                                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the free energy profiles in a,b,c directions using Widom insertion.              *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleFreeEnergyProfile(int Switch)
{
  int i,j,k,index,StartingBead;
  REAL value,IdealGasRosenBluth,error;
  VECTOR A,C,drift,fcom;
  char buffer[1024];
  FILE *FilePtr;
  VECTOR dr;
  REAL q;
  int StoredNumberOfTrialPositions;
  int StoredNumberOfTrialPositionsFirstBead;
  REAL UDeltaPolarization;

  switch(Switch)
  {
    case ALLOCATE:
      RosenBinSum=(REAL(***)[13])calloc(NumberOfSystems,sizeof(REAL(**)[13]));
      RosenBinSumSquared=(REAL(***)[13])calloc(NumberOfSystems,sizeof(REAL(**)[13]));
      RosenBinCount=(REAL(***)[13])calloc(NumberOfSystems,sizeof(REAL(**)[13]));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeFreeEnergyProfile[i])
        {
          RosenBinSum[i]=(REAL(**)[13])calloc(NumberOfComponents,sizeof(REAL(*)[13]));
          RosenBinSumSquared[i]=(REAL(**)[13])calloc(NumberOfComponents,sizeof(REAL(*)[13]));
          RosenBinCount[i]=(REAL(**)[13])calloc(NumberOfComponents,sizeof(REAL(*)[13]));
          for(j=0;j<NumberOfComponents;j++)
          {                                       // +1 to include index 0 also as the last value
            RosenBinSum[i][j]=(REAL(*)[13])calloc(FreeEnergyHistogramSize[i]+1,sizeof(REAL[13]));
            RosenBinSumSquared[i][j]=(REAL(*)[13])calloc(FreeEnergyHistogramSize[i]+1,sizeof(REAL[13]));
            RosenBinCount[i][j]=(REAL(*)[13])calloc(FreeEnergyHistogramSize[i]+1,sizeof(REAL[13]));
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      if(!ComputeFreeEnergyProfile[CurrentSystem]) return;

      fcom=GetFrameworkCenterOfMass();
      drift.x=fcom.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x;
      drift.y=fcom.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y;
      drift.z=fcom.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z;

      StoredNumberOfTrialPositions=NumberOfTrialPositions;
      StoredNumberOfTrialPositionsFirstBead=NumberOfTrialPositionsForTheFirstBead;

      for(i=0;i<MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem]);i++)
      {
        for(CurrentComponent=0;CurrentComponent<NumberOfComponents;CurrentComponent++)
        {
          if(Components[CurrentComponent].ComputeFreeEnergyProfile[CurrentSystem])
          {
            IdealGasRosenBluth=Components[CurrentComponent].IdealGasRosenbluthWeight[CurrentSystem];
            CurrentAdsorbateMolecule=NumberOfAdsorbateMolecules[CurrentSystem];
            CurrentCationMolecule=NumberOfCationMolecules[CurrentSystem];
            StartingBead=Components[CurrentComponent].StartingBead;

            // compute point 'A' as a random fractional position within the first unit cell
            // compute point 'C' as the corresponding Cartesian position
            // the fractional position is used for the 'a','b','c'-mapping
            // the Cartesian position is used for the 'diagonal'-mappings

            do
            {
              // generate random fractional point between 0.0 and 1.0 and try again when it is not inside the specified range
              A=GenerateRandomPointInCartesianUnitCellSpace();
            }while(!ValidFractionalPoint(CurrentComponent,A));

            // compute the Cartesian position 'C' in the first unit cell
            C.x=UnitCellBox[CurrentSystem].ax*A.x+UnitCellBox[CurrentSystem].bx*A.y+UnitCellBox[CurrentSystem].cx*A.z;
            C.y=UnitCellBox[CurrentSystem].ay*A.x+UnitCellBox[CurrentSystem].by*A.y+UnitCellBox[CurrentSystem].cy*A.z;
            C.z=UnitCellBox[CurrentSystem].az*A.x+UnitCellBox[CurrentSystem].bz*A.y+UnitCellBox[CurrentSystem].cz*A.z;

            // grow a trial molecule at position 'C' from scratch, the new Rosenbluth weight is stored in variable 'value'
            NewPosition[CurrentSystem][StartingBead]=C;
            NumberOfBeadsAlreadyPlaced=0;

            for(k=0;k<Components[CurrentComponent].NumberOfAtoms;k++)
            {
              CFVDWScaling[k]=1.0;
              CFChargeScaling[k]=1.0;
            }

            NumberOfTrialPositions=NumberOfTrialPositionsWidom;
            NumberOfTrialPositionsForTheFirstBead=NumberOfTrialPositionsForTheFirstBeadWidom;
            value=GrowMolecule(CBMC_PARTIAL_INSERTION);
            NumberOfTrialPositions=StoredNumberOfTrialPositions;
            NumberOfTrialPositionsForTheFirstBead=StoredNumberOfTrialPositionsFirstBead;

            if(OVERLAP||BlockedPocket(NewPosition[CurrentSystem][StartingBead]))
              value=0.0; // set the Rosenbluth weight at zero if an overlap occurs or when the position is 'blocked'
            else
            {
             // correct for tail-corrections
              value*=exp(-Beta[CurrentSystem]*TailMolecularEnergyDifferenceAdd());

              // correct for the long-range part of the Ewald-summation (Fourier part)
              if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
              {
                CalculateEwaldFourierAdsorbate(TRUE,FALSE,NumberOfAdsorbateMolecules[CurrentSystem],0);
                value*=exp(-Beta[CurrentSystem]*(
                    UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                    UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                    UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
              }
/*
              if(ComputePolarization)
              {
                ComputeNewPolarizationEnergy(TRUE,NumberOfAdsorbateMolecules[CurrentSystem],-1);
                UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                                   (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                                   UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                                   UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                                   (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                                   UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
                value*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
              }
*/

            }


            switch(FreeEnergyMappingType[CurrentSystem])
            {
              case NO_MAPPING:
                break;
              case A_MAPPING:
                index=(int)(A.x*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][0]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][0]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][0]+=1.0;
                }
                break;
              case B_MAPPING:
                index=(int)(A.y*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][1]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][1]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][1]+=1.0;
                }
                break;
              case C_MAPPING:
                index=(int)(A.z*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][2]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][2]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][2]+=1.0;
                }
                break;
              case ABC_MAPPING:
                index=(int)(A.x*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][0]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][0]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][0]+=1.0;
                }
                index=(int)(A.y*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][1]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][1]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][1]+=1.0;
                }
                index=(int)(A.z*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][2]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][2]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][2]+=1.0;
                }
                break;
              case MAP_AB_DIAGONAL:   // e.g. DDR
                dr.x=M_SQRT1_2;
                dr.y=-M_SQRT1_2;
                dr.z=0.0;
                q=((1.0-A.x)*dr.x+(-A.y)*dr.y+(-A.z)*dr.z)*M_SQRT1_2;
                index=(int)(q*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][3]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][3]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][3]+=1.0;
                }
                break;
              case MAP_AC_DIAGONAL:
                dr.x=M_SQRT1_2;
                dr.y=0.0;
                dr.z=-M_SQRT1_2;
                q=((1.0-A.x)*dr.x+(-A.y)*dr.y+(-A.z)*dr.z)*M_SQRT1_2;
                index=(int)(q*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][4]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][4]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][4]+=1.0;
                }
                break;
              case MAP_BC_DIAGONAL:
                dr.x=0.0;
                dr.y=M_SQRT1_2;
                dr.z=-M_SQRT1_2;
                q=((-A.x)*dr.x+(1.0-A.y)*dr.y+(-A.z)*dr.z)*M_SQRT1_2;
                index=(int)(q*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][5]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][5]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][5]+=1.0;
                }
                break;
              case MAP_O_AB_DIAGONAL:
                dr.x=-M_SQRT1_2;
                dr.y=-M_SQRT1_2;
                dr.z=0.0;
                q=((-A.x)*dr.x+(-A.y)*dr.y+(-A.z)*dr.z)*M_SQRT1_2;
                index=(int)(q*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][6]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][6]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][6]+=1.0;
                }
                break;
              case MAP_O_AC_DIAGONAL:
                dr.x=-M_SQRT1_2;
                dr.y=0.0;
                dr.z=-M_SQRT1_2;
                q=((-A.x)*dr.x+(-A.y)*dr.y+(-A.z)*dr.z)*M_SQRT1_2;
                index=(int)(q*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][7]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][7]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][7]+=1.0;
                }
                break;
              case MAP_O_BC_DIAGONAL:
                dr.x=0.0;
                dr.y=-M_SQRT1_2;
                dr.z=-M_SQRT1_2;
                q=((-A.x)*dr.x+(-A.y)*dr.y+(-A.z)*dr.z)*M_SQRT1_2;
                index=(int)(q*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][8]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][8]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][8]+=1.0;
                }
                break;
              case MAP_A_BC_DIAGONAL:
                dr.x=M_SQRT1_3;
                dr.y=-M_SQRT1_3;
                dr.z=-M_SQRT1_3;
                q=((1.0-A.x)*dr.x+(-A.y)*dr.y+(-A.z)*dr.z)*M_SQRT1_3;
                index=(int)(q*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][9]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][9]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][9]+=1.0;
                }
                break;
              case MAP_B_AC_DIAGONAL:
                dr.x=-M_SQRT1_3;
                dr.y=M_SQRT1_3;
                dr.z=-M_SQRT1_3;
                q=((-A.x)*dr.x+(1.0-A.y)*dr.y+(-A.z)*dr.z)*M_SQRT1_3;
                index=(int)(q*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][10]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][10]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][10]+=1.0;
                }
                break;
              case MAP_C_AB_DIAGONAL:
                dr.x=-M_SQRT1_3;
                dr.y=-M_SQRT1_3;
                dr.z=M_SQRT1_3;
                q=((-A.x)*dr.x+(-A.y)*dr.y+(1.0-A.z)*dr.z)*M_SQRT1_3;
                index=(int)(q*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][11]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][11]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][11]+=1.0;
                }
                break;
              case MAP_O_ABC_DIAGONAL:
                dr.x=-M_SQRT1_3;
                dr.y=-M_SQRT1_3;
                dr.z=-M_SQRT1_3;
                q=((-A.x)*dr.x+(-A.y)*dr.y+(-A.z)*dr.z)*M_SQRT1_3;
                index=(int)(q*(REAL)FreeEnergyHistogramSize[CurrentSystem]);
                if(index>=0&&index<FreeEnergyHistogramSize[CurrentSystem])
                {
                  RosenBinSum[CurrentSystem][CurrentComponent][index][12]+=value;
                  RosenBinSumSquared[CurrentSystem][CurrentComponent][index][12]+=SQR(value);
                  RosenBinCount[CurrentSystem][CurrentComponent][index][12]+=1.0;
                }
                break;
            }


          }
        }
      }
      break;
    case PRINT:
      // return if the rdf does not has to be calculated for this system
      if((!ComputeFreeEnergyProfile[CurrentSystem])||(CurrentCycle%WriteFreeEnergyProfileEvery[CurrentSystem]!=0)) return;

      // make the output directory
      mkdir("FreeEnergyProfile",S_IRWXU);

      // make the system directory in the output directory
      sprintf(buffer,"FreeEnergyProfile/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(CurrentComponent=0;CurrentComponent<NumberOfComponents;CurrentComponent++)
      {
        if(Components[CurrentComponent].ComputeFreeEnergyProfile[CurrentSystem])
        {
          switch(FreeEnergyMappingType[CurrentSystem])
          {
            case NO_MAPPING:
              break;
            case A_MAPPING:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_A",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)   // note the '<='
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];        // the last element of index is actually taken to be the value at 0
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][0]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][0]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][0]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][0]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][0]*(RosenBinCount[CurrentSystem][CurrentComponent][index][0]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][0]/RosenBinCount[CurrentSystem][CurrentComponent][index][0]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][0]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][0]/RosenBinCount[CurrentSystem][CurrentComponent][index][0])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case B_MAPPING:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_B",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][1]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][1]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][1]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][1]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][1]*(RosenBinCount[CurrentSystem][CurrentComponent][index][1]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][1]/RosenBinCount[CurrentSystem][CurrentComponent][index][1]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][1]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][1]/RosenBinCount[CurrentSystem][CurrentComponent][index][1])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case C_MAPPING:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_C",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][2]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][2]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][2]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][2]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][2]*(RosenBinCount[CurrentSystem][CurrentComponent][index][2]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][2]/RosenBinCount[CurrentSystem][CurrentComponent][index][2]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][2]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][2]/RosenBinCount[CurrentSystem][CurrentComponent][index][2])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case ABC_MAPPING:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_A",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)   // note the '<='
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];        // the last element of index is actually taken to be the value at 0
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][0]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][0]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][0]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][0]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][0]*(RosenBinCount[CurrentSystem][CurrentComponent][index][0]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][0]/RosenBinCount[CurrentSystem][CurrentComponent][index][0]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][0]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][0]/RosenBinCount[CurrentSystem][CurrentComponent][index][0])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_B",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][1]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][1]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][1]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][1]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][1]*(RosenBinCount[CurrentSystem][CurrentComponent][index][1]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][1]/RosenBinCount[CurrentSystem][CurrentComponent][index][1]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][1]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][1]/RosenBinCount[CurrentSystem][CurrentComponent][index][1])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_C",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][2]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][2]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][2]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][2]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][2]*(RosenBinCount[CurrentSystem][CurrentComponent][index][2]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][2]/RosenBinCount[CurrentSystem][CurrentComponent][index][2]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][2]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][2]/RosenBinCount[CurrentSystem][CurrentComponent][index][2])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case MAP_AB_DIAGONAL:   // e.g. DDR
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_AB_2D",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][3]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][3]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][3]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][3]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][3]*(RosenBinCount[CurrentSystem][CurrentComponent][index][3]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][3]/RosenBinCount[CurrentSystem][CurrentComponent][index][3]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][3]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][3]/RosenBinCount[CurrentSystem][CurrentComponent][index][3])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case MAP_AC_DIAGONAL:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_AC_2D",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][4]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][4]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][4]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][4]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][4]*(RosenBinCount[CurrentSystem][CurrentComponent][index][4]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][4]/RosenBinCount[CurrentSystem][CurrentComponent][index][4]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][4]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][4]/RosenBinCount[CurrentSystem][CurrentComponent][index][4])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case MAP_BC_DIAGONAL:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_BC_2D",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][5]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][5]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][5]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][5]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][5]*(RosenBinCount[CurrentSystem][CurrentComponent][index][5]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][5]/RosenBinCount[CurrentSystem][CurrentComponent][index][5]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][5]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][5]/RosenBinCount[CurrentSystem][CurrentComponent][index][5])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case MAP_O_AB_DIAGONAL:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_O_AB_2D",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][6]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][6]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][6]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][6]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][6]*(RosenBinCount[CurrentSystem][CurrentComponent][index][6]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][6]/RosenBinCount[CurrentSystem][CurrentComponent][index][6]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][6]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][6]/RosenBinCount[CurrentSystem][CurrentComponent][index][6])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case MAP_O_AC_DIAGONAL:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_O_AC_2D",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][7]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][7]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][7]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][7]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][7]*(RosenBinCount[CurrentSystem][CurrentComponent][index][7]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][7]/RosenBinCount[CurrentSystem][CurrentComponent][index][7]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][7]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][7]/RosenBinCount[CurrentSystem][CurrentComponent][index][7])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case MAP_O_BC_DIAGONAL:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_O_BC_2D",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][8]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][8]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][8]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][8]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][8]*(RosenBinCount[CurrentSystem][CurrentComponent][index][8]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][8]/RosenBinCount[CurrentSystem][CurrentComponent][index][8]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][8]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][8]/RosenBinCount[CurrentSystem][CurrentComponent][index][8])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case MAP_A_BC_DIAGONAL:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_A_BC_3D",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][9]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][9]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][9]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][9]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][9]*(RosenBinCount[CurrentSystem][CurrentComponent][index][9]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][9]/RosenBinCount[CurrentSystem][CurrentComponent][index][9]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][9]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][9]/RosenBinCount[CurrentSystem][CurrentComponent][index][9])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case MAP_B_AC_DIAGONAL:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_B_AC_3D",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][10]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][10]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][10]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][10]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][10]*(RosenBinCount[CurrentSystem][CurrentComponent][index][10]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][10]/RosenBinCount[CurrentSystem][CurrentComponent][index][10]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][10]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][10]/RosenBinCount[CurrentSystem][CurrentComponent][index][10])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case MAP_C_AB_DIAGONAL:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_C_BC_3D",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][11]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][11]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][11]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][11]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][11]*(RosenBinCount[CurrentSystem][CurrentComponent][index][11]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][11]/RosenBinCount[CurrentSystem][CurrentComponent][index][11]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][11]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][11]/RosenBinCount[CurrentSystem][CurrentComponent][index][11])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
            case MAP_O_ABC_DIAGONAL:
              sprintf(buffer,"FreeEnergyProfile/System_%d/FreeEnergy_%d_%s%s.dat_%d_O_ABC_3D",CurrentSystem,CurrentComponent,Components[CurrentComponent].Name,FileNameAppend,
                      Components[CurrentComponent].StartingBead);
              FilePtr=fopen(buffer,"w");
              for(j=0;j<=FreeEnergyHistogramSize[CurrentSystem];j++)
              {
                index=j%FreeEnergyHistogramSize[CurrentSystem];
                if(RosenBinCount[CurrentSystem][CurrentComponent][index][12]>0.5)
                {
                  error=sqrt((RosenBinCount[CurrentSystem][CurrentComponent][index][12]*RosenBinSumSquared[CurrentSystem][CurrentComponent][index][12]-
                       SQR(RosenBinSum[CurrentSystem][CurrentComponent][index][12]))/
                      (RosenBinCount[CurrentSystem][CurrentComponent][index][12]*(RosenBinCount[CurrentSystem][CurrentComponent][index][12]-1.0)));
                  error/=(RosenBinSum[CurrentSystem][CurrentComponent][index][12]/RosenBinCount[CurrentSystem][CurrentComponent][index][12]);
                  error/=sqrt(RosenBinCount[CurrentSystem][CurrentComponent][index][12]);
                  error*=2.0;  // 95 % confidence interval is 2*sigma
                  fprintf(FilePtr,"%g %g %g\n",
                    (double)((REAL)j/(REAL)FreeEnergyHistogramSize[CurrentSystem]),
                    (double)(-log(RosenBinSum[CurrentSystem][CurrentComponent][index][12]/RosenBinCount[CurrentSystem][CurrentComponent][index][12])),
                    (double)error);
                }
              }
              fclose(FilePtr);
              break;
          }
        }
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeFreeEnergyProfile[i])
        {
          for(j=0;j<NumberOfComponents;j++)
          {                     
            free(RosenBinSum[i][j]);
            free(RosenBinSumSquared[i][j]);
            free(RosenBinCount[i][j]);
          }
          free(RosenBinSum[i]);
          free(RosenBinSumSquared[i]);
          free(RosenBinCount[i]);
        }
      }
      free(RosenBinSum);
      free(RosenBinSumSquared);
      free(RosenBinCount);
      break;
  }
}


/*********************************************************************************************************
 * Name       | CheckSphereOverlap                                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the largest free radius of a sphere that does not overlap with the framework    *
 *            | given that the center of the sphere is 'posA'.                                           *
 * Parameters | -                                                                                        *
 * Used in    | Auxiliary function used in 'SamplePoreSizeDistribution'                                  *
 *********************************************************************************************************/

int CheckSphereOverlap(VECTOR posA,REAL *SmallestRadius)
{
  int i,j;
  int typeB,f1;
  VECTOR dr,posB;
  REAL rr,well_depth_factor;
  REAL Radius;

  well_depth_factor=Framework[CurrentSystem].PoreSizeDistributionProbeDistance;

  *SmallestRadius=DBL_MAX;

  if(BlockedPocket(posA))
    return TRUE;

  // check overlap with framework
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
    {
      posB=Framework[CurrentSystem].Atoms[f1][j].Position;
      typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      if(rr<SQR(0.5*well_depth_factor*PotentialParms[typeB][typeB][1]))
        return TRUE;

      Radius=sqrt(rr)-0.5*well_depth_factor*PotentialParms[typeB][typeB][1];
      if(Radius<*SmallestRadius) *SmallestRadius=Radius;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeB=Cations[CurrentSystem][i].Atoms[j].Type;
      posB=Cations[CurrentSystem][i].Atoms[j].Position;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      if(rr<SQR(0.5*well_depth_factor*PotentialParms[typeB][typeB][1]))
        return TRUE;

      Radius=sqrt(rr)-0.5*well_depth_factor*PotentialParms[typeB][typeB][1];
      if(Radius<*SmallestRadius) *SmallestRadius=Radius;
    }
  }
  return FALSE;
}


/*********************************************************************************************************
 * Name       | SamplePoreSizeDistribution                                                               *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the Pore-Size Distribution function (PSD).                                       *
 * Parameters | -                                                                                        *
 * Note       | Based on code of Lev Sarkisov who wrote a PSD program after discussions with Lev Gelb.   *
 *********************************************************************************************************/

void SamplePoreSizeDistribution(int Switch)
{
  int i,j,k;
  int index;
  VECTOR posA,posB,dr;
  REAL Radius,LargestRadius,rr,r;
  REAL deltaR;
  FILE *FilePtr;
  REAL temp,deriv;
  REAL left,mid,right;
  char buffer[256];


  switch(Switch)
  {
    case ALLOCATE:
      PoreSizeDistributionHistogram=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputePSDHistogram[i])
          PoreSizeDistributionHistogram[i]=(REAL*)calloc(PSDHistogramSize[i],sizeof(REAL));
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      if(!ComputePSDHistogram[CurrentSystem]) return;

      deltaR=PSDRange[CurrentSystem]/PSDHistogramSize[CurrentSystem];

      // choose a random point in the simulation box
      posA.x=RandomNumber();
      posA.y=RandomNumber();
      posA.z=RandomNumber();
      posA=ConvertFromABCtoXYZ(posA);

      CurrentComponent=0;

      // if overlap with framework then return
      if(CheckSphereOverlap(posA,&Radius)) return;

      LargestRadius=-10.0;

      for(j=0;j<20*NumberOfCycles;j++)
      {
        // choose a random center of a sphere in the simulation box
        posB.x=RandomNumber();
        posB.y=RandomNumber();
        posB.z=RandomNumber();
        posB=ConvertFromABCtoXYZ(posB);

        // Computes the largest free radius of a sphere that does not overlap with the framework
        // given that the center of the sphere is 'posB'.
        if(CheckSphereOverlap(posB,&Radius)) continue;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        // check that point A is inside sphere of point B
        if(r>Radius) continue;

        if(Radius>LargestRadius) LargestRadius=Radius;
      }

      index=(int)(LargestRadius/deltaR);
      if((index>=0)&&(index<PSDHistogramSize[CurrentSystem]))
        for(k=0;k<=index;k++)
          PoreSizeDistributionHistogram[CurrentSystem][k]+=1.0;
      break;
    case PRINT:
      if((!ComputePSDHistogram[CurrentSystem])||(CurrentCycle%WritePSDHistogramEvery[CurrentSystem]!=0)) return;

      deltaR=PSDRange[CurrentSystem]/PSDHistogramSize[CurrentSystem];

      if (STREAM)
      {
#ifdef __unix__
        // Loads PSD output contents into a global
        FilePtr=open_memstream(&PORE_SIZE_DISTRIBUTION_OUTPUT,
                               &PORE_SIZE_DISTRIBUTION_OUTPUT_SIZE);
#else
        fprintf(stderr, "Streaming only allowed on POSIX systems (for now)\n.");
        exit(1);
#endif
      }
      else
      {
        mkdir("PoreSizeDistributionHistogram",S_IRWXU);
        sprintf(buffer,"PoreSizeDistributionHistogram/System_%d",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        sprintf(buffer,"PoreSizeDistributionHistogram/System_%d/HistogramPoreSizeDistribution_%s_%d.%d.%d_%lf_%lf%s.dat",
                CurrentSystem,
                Framework[CurrentSystem].Name[0],
                NumberOfUnitCells[CurrentSystem].x,
                NumberOfUnitCells[CurrentSystem].y,
                NumberOfUnitCells[CurrentSystem].z,
                (double)therm_baro_stats.ExternalTemperature[CurrentSystem],
                (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                FileNameAppend);
        FilePtr=fopen(buffer,"w");
      }

      fprintf(FilePtr,"# column 1: dameter d [A]\n");
      fprintf(FilePtr,"# column 2: Connoly distribution\n");
      fprintf(FilePtr,"# column 3: PSD\n");
      fprintf(FilePtr,"# value at d=0 is void-fraction\n");

      for(k=0;k<PSDHistogramSize[CurrentSystem];k++)
      {
        if(k==0)
        {
          mid=PoreSizeDistributionHistogram[CurrentSystem][k]/CurrentCycle;
          right=PoreSizeDistributionHistogram[CurrentSystem][k+1]/CurrentCycle;
          deriv=(right-mid)/(deltaR);
        }
        else if(k==(PSDHistogramSize[CurrentSystem]-1))
        {
          left=PoreSizeDistributionHistogram[CurrentSystem][k-1]/CurrentCycle;
          mid=PoreSizeDistributionHistogram[CurrentSystem][k]/CurrentCycle;
          deriv=(mid-left)/(deltaR);
        }
        else
        {
          left=PoreSizeDistributionHistogram[CurrentSystem][k-1]/CurrentCycle;
          mid=PoreSizeDistributionHistogram[CurrentSystem][k]/CurrentCycle;
          right=PoreSizeDistributionHistogram[CurrentSystem][k+1]/CurrentCycle;
          deriv=(right-left)/(2.0*deltaR);
        }

        fprintf(FilePtr,"%g %g %g\n",2.0*deltaR*(k+0.5),mid,-deriv);
      }
      fclose(FilePtr);
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputePSDHistogram[i])
          free(PoreSizeDistributionHistogram[i]);
      }
      free(PoreSizeDistributionHistogram);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleEndToEndDistanceHistogram                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the end-to-end distance histograms.                                              *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleEndToEndDistanceHistogram(int Switch)
{
  int i,j,k,index;
  int A,B,Type;
  VECTOR dr;
  REAL r,norm,delta;
  FILE *FilePtr;
  char buffer[256];


  switch(Switch)
  {
    case ALLOCATE:
      EndToEndDistanceHistogram=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeEndToEndDistanceHistogram[i])
        {
          EndToEndDistanceHistogram[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
          for(j=0;j<NumberOfComponents;j++)
            EndToEndDistanceHistogram[i][j]=(REAL*)calloc(EndToEndHistogramSize[i],sizeof(REAL));
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      if(!ComputeEndToEndDistanceHistogram[CurrentSystem]) return;

      delta=EndToEndHistogramSize[CurrentSystem]/EndToEndRange[CurrentSystem];

      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        Type=Adsorbates[CurrentSystem][i].Type;
        A=0;
        B=Components[Type].NumberOfAtoms-1;
        if(B>0) // no end-to-end for single atom
        {
          dr.x=Adsorbates[CurrentSystem][i].Atoms[A].Position.x-Adsorbates[CurrentSystem][i].Atoms[B].Position.x;
          dr.y=Adsorbates[CurrentSystem][i].Atoms[A].Position.y-Adsorbates[CurrentSystem][i].Atoms[B].Position.y;
          dr.z=Adsorbates[CurrentSystem][i].Atoms[A].Position.z-Adsorbates[CurrentSystem][i].Atoms[B].Position.z;
          r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
          index=(int)(r*delta);
          if(index>=0&&index<EndToEndHistogramSize[CurrentSystem])
            EndToEndDistanceHistogram[CurrentSystem][Type][index]+=1.0;
        }
      }
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        Type=Cations[CurrentSystem][i].Type;
        A=0;
        B=Components[Type].NumberOfAtoms-1;
        if(B>0) // no end-to-end for single atom
        {
          dr.x=Cations[CurrentSystem][i].Atoms[A].Position.x-Cations[CurrentSystem][i].Atoms[B].Position.x;
          dr.y=Cations[CurrentSystem][i].Atoms[A].Position.y-Cations[CurrentSystem][i].Atoms[B].Position.y;
          dr.z=Cations[CurrentSystem][i].Atoms[A].Position.z-Cations[CurrentSystem][i].Atoms[B].Position.z;
          r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
          index=(int)(r*delta);
          if(index>=0&&index<EndToEndHistogramSize[CurrentSystem])
            EndToEndDistanceHistogram[CurrentSystem][Type][index]+=1.0;
        }
      }
      break;
    case PRINT:
      if((!ComputeEndToEndDistanceHistogram[CurrentSystem])||(CurrentCycle%WriteEndToEndDistanceHistogramEvery[CurrentSystem]!=0)) return;

      delta=EndToEndHistogramSize[CurrentSystem]/EndToEndRange[CurrentSystem];

      mkdir("EndToEndDistanceHistograms",S_IRWXU);
      sprintf(buffer,"EndToEndDistanceHistograms/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(i=0;i<NumberOfComponents;i++)
      {
        norm=0.0;
        for(k=0;k<EndToEndHistogramSize[CurrentSystem];k++)
          norm+=EndToEndDistanceHistogram[CurrentSystem][i][k];
        norm/=delta;

        sprintf(buffer,"EndToEndDistanceHistograms/System_%d/Histogram_%s_%d_Density_%g%s.dat",
                CurrentSystem,Components[i].Name,i,
                (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                FileNameAppend);
        FilePtr=fopen(buffer,"w");
        for(k=0;k<EndToEndHistogramSize[CurrentSystem];k++)
        {
          r=(REAL)k/delta;
          if(EndToEndDistanceHistogram[CurrentSystem][i][k]>0.0)
            fprintf(FilePtr,"%g %g\n",(double)r,(double)(EndToEndDistanceHistogram[CurrentSystem][i][k]/norm));
        }
        fclose(FilePtr);
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeEndToEndDistanceHistogram[i])
        {
          for(j=0;j<NumberOfComponents;j++)
            free(EndToEndDistanceHistogram[i][j]);
          free(EndToEndDistanceHistogram[i]);
        }
      }
      free(EndToEndDistanceHistogram);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleEnergyHistogram                                                                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the histogram of the energy.                                                     *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleEnergyHistogram(int Switch)
{
  int i,j,k,index;
  REAL norm,U,Range;
  char buffer[256];
  FILE *FilePtr;

  switch(Switch)
  {
    case ALLOCATE:
      EnergyHistogram=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
        EnergyHistogram[i]=(REAL**)calloc(6,sizeof(REAL*));
        for(j=0;j<6;j++)
          EnergyHistogram[i][j]=(REAL*)calloc(EnergyHistogramSize[i],sizeof(REAL));
      }
      break;
    case INITIALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        for(j=0;j<6;j++)
          for(k=0;k<EnergyHistogramSize[i];k++)
            EnergyHistogram[i][j][k]=0.0;
      }
      break;
    case SAMPLE:
      if(!ComputeEnergyHistogram[CurrentSystem]) return;

      Range=fabs(EnergyHistogramUpperLimit[CurrentSystem]-EnergyHistogramLowerLimit[CurrentSystem]);

      U=UTotal[CurrentSystem]*ENERGY_TO_KELVIN;
      index=(int)((U-EnergyHistogramLowerLimit[CurrentSystem])*EnergyHistogramSize[CurrentSystem]/Range);
      if(index>=0&&index<EnergyHistogramSize[CurrentSystem])
        EnergyHistogram[CurrentSystem][0][index]+=1.0;

      U=(UHostHostCoulomb[CurrentSystem]+UHostAdsorbateCoulomb[CurrentSystem]+UHostCationCoulomb[CurrentSystem]+
         UAdsorbateAdsorbateCoulomb[CurrentSystem]+UCationCationCoulomb[CurrentSystem]+UAdsorbateCationCoulomb[CurrentSystem])*ENERGY_TO_KELVIN;
      index=(int)((U-EnergyHistogramLowerLimit[CurrentSystem])*EnergyHistogramSize[CurrentSystem]/Range);
      if(index>=0&&index<EnergyHistogramSize[CurrentSystem])
        EnergyHistogram[CurrentSystem][1][index]+=1.0;

      U=(UHostHostVDW[CurrentSystem]+UHostAdsorbateVDW[CurrentSystem]+UHostCationVDW[CurrentSystem]+
         UAdsorbateAdsorbateVDW[CurrentSystem]+UCationCationVDW[CurrentSystem]+UAdsorbateCationVDW[CurrentSystem])*ENERGY_TO_KELVIN;
      index=(int)((U-EnergyHistogramLowerLimit[CurrentSystem])*EnergyHistogramSize[CurrentSystem]/Range);
      if(index>=0&&index<EnergyHistogramSize[CurrentSystem])
        EnergyHistogram[CurrentSystem][2][index]+=1.0;

      U=(UHostPolarization[CurrentSystem]+UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem])*ENERGY_TO_KELVIN;
      index=(int)((U-EnergyHistogramLowerLimit[CurrentSystem])*EnergyHistogramSize[CurrentSystem]/Range);
      if(index>=0&&index<EnergyHistogramSize[CurrentSystem])
        EnergyHistogram[CurrentSystem][3][index]+=1.0;

      // Energy host-guest
      U=(UHostAdsorbate[CurrentSystem]+UHostCation[CurrentSystem])*ENERGY_TO_KELVIN;
      index=(int)((U-EnergyHistogramLowerLimit[CurrentSystem])*EnergyHistogramSize[CurrentSystem]/Range);
      if(index>=0&&index<EnergyHistogramSize[CurrentSystem])
        EnergyHistogram[CurrentSystem][4][index]+=1.0;

      // Energy guest-guest
      U=(UAdsorbateAdsorbate[CurrentSystem]+UAdsorbateCation[CurrentSystem]+UCationCation[CurrentSystem])*ENERGY_TO_KELVIN;
      index=(int)((U-EnergyHistogramLowerLimit[CurrentSystem])*EnergyHistogramSize[CurrentSystem]/Range);
      if(index>=0&&index<EnergyHistogramSize[CurrentSystem])
        EnergyHistogram[CurrentSystem][5][index]+=1.0;
      break;
    case PRINT:
      if((!ComputeEnergyHistogram[CurrentSystem])||(CurrentCycle%WriteEnergyHistogramEvery[CurrentSystem]!=0)) return;
      mkdir("EnergyHistograms",S_IRWXU);
      sprintf(buffer,"EnergyHistograms/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      Range=fabs(EnergyHistogramUpperLimit[CurrentSystem]-EnergyHistogramLowerLimit[CurrentSystem]);

      norm=0.0;
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
        norm+=EnergyHistogram[CurrentSystem][0][k];
      norm*=Range/(REAL)EnergyHistogramSize[CurrentSystem];

      sprintf(buffer,"EnergyHistograms/System_%d/Histogram_Total_Energy_%g%s.dat",CurrentSystem,
              (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
              FileNameAppend);
      FilePtr=fopen(buffer,"w");
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
      {
        U=(REAL)k*(Range/EnergyHistogramSize[CurrentSystem])+EnergyHistogramLowerLimit[CurrentSystem];
        if(EnergyHistogram[CurrentSystem][0][k]>0.0)
          fprintf(FilePtr,"%g %g\n",(double)U,(double)(EnergyHistogram[CurrentSystem][0][k]/norm));
      }
      fclose(FilePtr);

      norm=0.0;
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
        norm+=EnergyHistogram[CurrentSystem][1][k];
      norm*=Range/(REAL)EnergyHistogramSize[CurrentSystem];

      sprintf(buffer,"EnergyHistograms/System_%d/Histogram_Coulomb_Energy_%g%s.dat",CurrentSystem,
              (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
              FileNameAppend);
      FilePtr=fopen(buffer,"w");
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
      {
        U=(REAL)k*(Range/EnergyHistogramSize[CurrentSystem])+EnergyHistogramLowerLimit[CurrentSystem];
        if(EnergyHistogram[CurrentSystem][1][k]>0.0)
          fprintf(FilePtr,"%g %g\n",(double)U,(double)(EnergyHistogram[CurrentSystem][1][k]/norm));
      }
      fclose(FilePtr);

      norm=0.0;
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
        norm+=EnergyHistogram[CurrentSystem][2][k];
      norm*=Range/(REAL)EnergyHistogramSize[CurrentSystem];

      sprintf(buffer,"EnergyHistograms/System_%d/Histogram_VDW_Energy_%g%s.dat",CurrentSystem,
              (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
              FileNameAppend);
      FilePtr=fopen(buffer,"w");
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
      {
        U=(REAL)k*(Range/EnergyHistogramSize[CurrentSystem])+EnergyHistogramLowerLimit[CurrentSystem];
        if(EnergyHistogram[CurrentSystem][2][k]>0.0)
          fprintf(FilePtr,"%g %g\n",(double)U,(double)(EnergyHistogram[CurrentSystem][2][k]/norm));
      }
      fclose(FilePtr);

      norm=0.0;
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
        norm+=EnergyHistogram[CurrentSystem][3][k];
      norm*=Range/(REAL)EnergyHistogramSize[CurrentSystem];

      sprintf(buffer,"EnergyHistograms/System_%d/Histogram_Polarization_Energy_%g%s.dat",CurrentSystem,
              (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
              FileNameAppend);
      FilePtr=fopen(buffer,"w");
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
      {
        U=(REAL)k*(Range/EnergyHistogramSize[CurrentSystem])+EnergyHistogramLowerLimit[CurrentSystem];
        if(EnergyHistogram[CurrentSystem][3][k]>0.0)
          fprintf(FilePtr,"%g %g\n",(double)U,(double)(EnergyHistogram[CurrentSystem][3][k]/norm));
      }
      fclose(FilePtr);

      norm=0.0;
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
        norm+=EnergyHistogram[CurrentSystem][4][k];
      norm*=Range/(REAL)EnergyHistogramSize[CurrentSystem];

      sprintf(buffer,"EnergyHistograms/System_%d/Histogram_HostGuest_Energy_%g%s.dat",CurrentSystem,
              (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
              FileNameAppend);
      FilePtr=fopen(buffer,"w");
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
      {
        U=(REAL)k*(Range/EnergyHistogramSize[CurrentSystem])+EnergyHistogramLowerLimit[CurrentSystem];
        if(EnergyHistogram[CurrentSystem][4][k]>0.0)
          fprintf(FilePtr,"%g %g\n",(double)U,(double)(EnergyHistogram[CurrentSystem][4][k]/norm));
      }
      fclose(FilePtr);

      norm=0.0;
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
        norm+=EnergyHistogram[CurrentSystem][5][k];
      norm*=Range/(REAL)EnergyHistogramSize[CurrentSystem];

      sprintf(buffer,"EnergyHistograms/System_%d/Histogram_GuestGuest_Energy_%g%s.dat",CurrentSystem,
              (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
              FileNameAppend);
      FilePtr=fopen(buffer,"w");
      for(k=0;k<EnergyHistogramSize[CurrentSystem];k++)
      {
        U=(REAL)k*(Range/EnergyHistogramSize[CurrentSystem])+EnergyHistogramLowerLimit[CurrentSystem];
        if(EnergyHistogram[CurrentSystem][5][k]>0.0)
          fprintf(FilePtr,"%g %g\n",(double)U,(double)(EnergyHistogram[CurrentSystem][5][k]/norm));
      }
      fclose(FilePtr);
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        for(j=0;j<6;j++)
          free(EnergyHistogram[i][j]);
        free(EnergyHistogram[i]);
      }
      free(EnergyHistogram);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleThermoDynamicsFactor                                                               *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the thermodynamic factor.                                                        *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleThermoDynamicsFactor(int Choice)
{
  int i,j,k;
  REAL temp1,temp2,temp3;
  char buffer[256];
  FILE *FilePtr;

  switch(Choice)
  {
    case ALLOCATE:
      ThermoDynamicFactorNumberOfMoleculesCrossTerm=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      ThermoDynamicFactorNumberOfMolecules=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
      ThermoDynamicFactorNumberOfSamples=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

      for(k=0;k<NumberOfSystems;k++)
      {
        if(ComputeThermoDynamicFactor[k])
        {
          ThermoDynamicFactorNumberOfMoleculesCrossTerm[k]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
          ThermoDynamicFactorNumberOfMolecules[k]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
          for(i=0;i<NumberOfComponents;i++)
            ThermoDynamicFactorNumberOfMoleculesCrossTerm[k][i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
        }
      }
      break;
    case INITIALIZE:
      for(k=0;k<NumberOfSystems;k++)
      {
        if(ComputeThermoDynamicFactor[k])
        {
          ThermoDynamicFactorNumberOfSamples[k]=0.0;
          for(i=0;i<NumberOfComponents;i++)
          {
            ThermoDynamicFactorNumberOfMolecules[k][i]=0;
            for(j=0;j<NumberOfComponents;j++)
              ThermoDynamicFactorNumberOfMoleculesCrossTerm[k][i][j]=0.0;
          }
        }
      }
      break;
    case SAMPLE:
      if(!ComputeThermoDynamicFactor[CurrentSystem]) return;
      for(i=0;i<NumberOfComponents;i++)
      {
        ThermoDynamicFactorNumberOfMolecules[CurrentSystem][i]+=Components[i].NumberOfMolecules[CurrentSystem];
        for(j=0;j<NumberOfComponents;j++)
          ThermoDynamicFactorNumberOfMoleculesCrossTerm[CurrentSystem][i][j]+=Components[i].NumberOfMolecules[CurrentSystem]*Components[j].NumberOfMolecules[CurrentSystem];
      }
      ThermoDynamicFactorNumberOfSamples[CurrentSystem]+=1.0;
      break;
    case PRINT:
      if((!ComputeThermoDynamicFactor[CurrentSystem])||(CurrentCycle%WriteThermoDynamicFactorEvery[CurrentSystem]!=0)) return;
      mkdir("ThermoDynamicFactor",S_IRWXU);
      sprintf(buffer,"ThermoDynamicFactor/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      sprintf(buffer,"ThermoDynamicFactor/System_%d/ThermoDynamicFactor%s.dat",CurrentSystem,FileNameAppend);
      FilePtr=fopen(buffer,"w");
      for(i=0;i<NumberOfComponents;i++)
        fprintf(FilePtr,"# component %d: %s\n",i,Components[i].Name);
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"# <N_i N_j>-<N_i><N_j>\n");
      for(i=0;i<NumberOfComponents;i++)
      {
        for(j=0;j<NumberOfComponents;j++)
        {
          temp1=ThermoDynamicFactorNumberOfMoleculesCrossTerm[CurrentSystem][i][j]/ThermoDynamicFactorNumberOfSamples[CurrentSystem];
          temp2=ThermoDynamicFactorNumberOfMolecules[CurrentSystem][i]/ThermoDynamicFactorNumberOfSamples[CurrentSystem];
          temp3=ThermoDynamicFactorNumberOfMolecules[CurrentSystem][j]/ThermoDynamicFactorNumberOfSamples[CurrentSystem];
          fprintf(FilePtr,"%g ",temp1-temp2*temp3);
        }
        fprintf(FilePtr,"\n");
      }
      fprintf(FilePtr,"\n");
      fclose(FilePtr);
      break;
    case FINALIZE:
      for(k=0;k<NumberOfSystems;k++)
      {
        if(ComputeThermoDynamicFactor[k])
        {
          for(i=0;i<NumberOfComponents;i++)
            free(ThermoDynamicFactorNumberOfMoleculesCrossTerm[k][i]);
          free(ThermoDynamicFactorNumberOfMolecules[k]);
          free(ThermoDynamicFactorNumberOfMoleculesCrossTerm[k]);
        }
      }
      free(ThermoDynamicFactorNumberOfSamples);
      free(ThermoDynamicFactorNumberOfMolecules);
      free(ThermoDynamicFactorNumberOfMoleculesCrossTerm);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleFrameworkSpacingHistogram                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the inter-framework spacing histogram.                                           *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleFrameworkSpacingHistogram(int Choice)
{
  int i,j,k,l,f1,f2;
  int index;
  VECTOR dr;
  POINT posf1,posf2;
  REAL r,norm[4],delta;
  char buffer[256];
  FILE *FilePtr;


  switch(Choice)
  {
    case ALLOCATE:
      FrameworkDistanceHistogram=(REAL*****)calloc(NumberOfSystems,sizeof(REAL****));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeFrameworkSpacingHistogram[i])
        {
          FrameworkDistanceHistogram[i]=(REAL****)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL***));
          for(j=0;j<Framework[i].NumberOfFrameworks;j++)
          {
            FrameworkDistanceHistogram[i][j]=(REAL***)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL**));
            for(k=0;k<Framework[i].NumberOfFrameworks;k++)
            {
              FrameworkDistanceHistogram[i][j][k]=(REAL**)calloc(4,sizeof(REAL*));
              for(l=0;l<4;l++)
                FrameworkDistanceHistogram[i][j][k][l]=(REAL*)calloc(FrameworkSpacingHistogramSize[i],sizeof(REAL));
            }
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      if(!ComputeFrameworkSpacingHistogram[CurrentSystem]) return;

      delta=FrameworkSpacingHistogramSize[CurrentSystem]/FrameworkSpacingRange[CurrentSystem];

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        CurrentFramework=f1;
        posf1=GetFrameworkCenterOfMassPosition();
        for(f2=f1+1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
        {
          CurrentFramework=f2;
          posf2=GetFrameworkCenterOfMassPosition();
          dr.x=posf1.x-posf2.x;
          dr.y=posf1.y-posf2.y;
          dr.z=posf1.z-posf2.z;
          dr=ApplyBoundaryCondition(dr);
          r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

          // xyz
          index=(int)(r*delta);
          if(index>=0&&index<FrameworkSpacingHistogramSize[CurrentSystem])
            FrameworkDistanceHistogram[CurrentSystem][f1][f2][0][index]+=1.0;

          // x
          index=(int)(fabs(dr.x)*delta);
          if(index>=0&&index<FrameworkSpacingHistogramSize[CurrentSystem])
            FrameworkDistanceHistogram[CurrentSystem][f1][f2][1][index]+=1.0;

          // y
          index=(int)(fabs(dr.y)*delta);
          if(index>=0&&index<FrameworkSpacingHistogramSize[CurrentSystem])
            FrameworkDistanceHistogram[CurrentSystem][f1][f2][2][index]+=1.0;

          // z
          index=(int)(fabs(dr.z)*delta);
          if(index>=0&&index<FrameworkSpacingHistogramSize[CurrentSystem])
            FrameworkDistanceHistogram[CurrentSystem][f1][f2][3][index]+=1.0;
        }
      }
      break;
    case PRINT:
      if((!ComputeFrameworkSpacingHistogram[CurrentSystem])||(CurrentCycle%WriteFrameworkSpacingHistogramEvery[CurrentSystem]!=0)) return;

      delta=FrameworkSpacingHistogramSize[CurrentSystem]/FrameworkSpacingRange[CurrentSystem];

      mkdir("FrameworkDistanceHistogram",S_IRWXU);
      sprintf(buffer,"FrameworkDistanceHistogram/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(f2=f1+1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
        {
          for(i=0;i<4;i++)
          {
            norm[i]=0.0;
            for(k=0;k<FrameworkSpacingHistogramSize[CurrentSystem];k++)
              norm[i]+=FrameworkDistanceHistogram[CurrentSystem][f1][f2][i][k];
            norm[i]/=delta;
          }

          sprintf(buffer,"FrameworkDistanceHistogram/System_%d/Histogram_%d_%d_FrameworkDistance_%g%s.dat",
              CurrentSystem,f1,f2,therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR,
              FileNameAppend);
          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# xyz, x, y, z distances of the frameworks in Angstrom\n");
          fprintf(FilePtr,"# original distances: %g %g %g %g\n",
               OriginalFrameworkShiftDirAvg[CurrentSystem][f1][f2],
               OriginalFrameworkShift[CurrentSystem][f1][f2].x,
               OriginalFrameworkShift[CurrentSystem][f1][f2].y,
               OriginalFrameworkShift[CurrentSystem][f1][f2].z);
          for(k=0;k<FrameworkSpacingHistogramSize[CurrentSystem];k++)
          {
            r=(REAL)k/delta;
            fprintf(FilePtr,"%g %g %g %g %g\n",(double)r,
              (double)(norm[0]>0.0?(FrameworkDistanceHistogram[CurrentSystem][f1][f2][0][k]/norm[0]):0.0),
              (double)(norm[1]>0.0?(FrameworkDistanceHistogram[CurrentSystem][f1][f2][1][k]/norm[1]):0.0),
              (double)(norm[2]>0.0?(FrameworkDistanceHistogram[CurrentSystem][f1][f2][2][k]/norm[2]):0.0),
              (double)(norm[3]>0.0?(FrameworkDistanceHistogram[CurrentSystem][f1][f2][3][k]/norm[3]):0.0));
          }
          fclose(FilePtr);
        }
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeFrameworkSpacingHistogram[i])
        {
          for(j=0;j<Framework[i].NumberOfFrameworks;j++)
          {
            for(k=0;k<Framework[i].NumberOfFrameworks;k++)
            {
              for(l=0;l<4;l++)
                free(FrameworkDistanceHistogram[i][j][k][l]);
              free(FrameworkDistanceHistogram[i][j][k]);
            }
            free(FrameworkDistanceHistogram[i][j]);
          }
          free(FrameworkDistanceHistogram[i]);
        }
      }
      free(FrameworkDistanceHistogram);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleResidenceTimes                                                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the histograms of the residence times.                                           *
 * Parameters | -                                                                                        *
 * Note       | Samples the distribution of residence-times of molecules inside 'blocking-spheres'       *
 *            | (usually used to block inaccessible regions in Monte-Carlo). Also computes the fraction  *
 *            | each molecule spends inside a pocket. For a well equilibrated MD simulation these should *
 *            | be equal for the same type of molecule.                                                  *
 *********************************************************************************************************/

void SampleResidenceTimes(int Switch)
{
  int i,j,k,index;
  int StartingBead,Type;
  REAL time_difference,deltaR,temp;
  char buffer[256];
  FILE *FilePtr;

  switch(Switch)
  {
    case ALLOCATE:
      ResidenceOriginAdsorbates=(long long **)calloc(NumberOfSystems,sizeof(long long*));
      ResidenceOriginCations=(long long **)calloc(NumberOfSystems,sizeof(long long*));
      ResidenceStatusAdsorbates=(int **)calloc(NumberOfSystems,sizeof(int*));
      ResidenceStatusCations=(int**)calloc(NumberOfSystems,sizeof(int*));
      ResidenceTimesHistogram=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
      ResidenceTimesFractionAdsorbates=(REAL(**)[NR_BLOCKS])calloc(NumberOfSystems,sizeof(REAL(*)[NR_BLOCKS]));
      ResidenceTimesFractionCations=(REAL(**)[NR_BLOCKS])calloc(NumberOfSystems,sizeof(REAL(*)[NR_BLOCKS]));
      ResidenceTimesFractionCounts=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

      for(i=0;i<NumberOfSystems;i++)
      {
        ResidenceOriginAdsorbates[i]=(long long*)calloc(NumberOfAdsorbateMolecules[i],sizeof(long long));
        ResidenceOriginCations[i]=(long long*)calloc(NumberOfCationMolecules[i],sizeof(long long));
        ResidenceStatusAdsorbates[i]=(int*)calloc(NumberOfAdsorbateMolecules[i],sizeof(int));
        ResidenceStatusCations[i]=(int*)calloc(NumberOfCationMolecules[i],sizeof(int));
        ResidenceTimesHistogram[i]=(REAL*)calloc(ResidenceTimesHistogramSize[i],sizeof(REAL));
        ResidenceTimesFractionAdsorbates[i]=(REAL(*)[NR_BLOCKS])calloc(NumberOfAdsorbateMolecules[i],sizeof(REAL[NR_BLOCKS]));
        ResidenceTimesFractionCations[i]=(REAL(*)[NR_BLOCKS])calloc(NumberOfCationMolecules[i],sizeof(REAL[NR_BLOCKS]));
      }
      break;
    case INITIALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeResidenceTimes[i])
        {
          ResidenceTimesFractionCounts[i]=0;
          for(j=0;j<NumberOfAdsorbateMolecules[i];j++)
          {
            Type=Adsorbates[i][j].Type;
            StartingBead=Components[Type].StartingBead;

            for(k=0;k<NR_BLOCKS;k++)
              ResidenceTimesFractionAdsorbates[i][j][k]=0.0;

            if(BlockedPocket(Adsorbates[i][j].Atoms[StartingBead].Position))
              ResidenceStatusAdsorbates[i][j]=TRUE;
            else
              ResidenceStatusAdsorbates[i][j]=FALSE;
          }
          for(j=0;j<NumberOfCationMolecules[i];j++)
          {
            Type=Cations[i][j].Type;
            StartingBead=Components[Type].StartingBead;

            for(k=0;k<NR_BLOCKS;k++)
              ResidenceTimesFractionCations[i][j][k]=0.0;

            if(BlockedPocket(Cations[i][j].Atoms[StartingBead].Position))
              ResidenceStatusCations[i][j]=TRUE;
            else
              ResidenceStatusCations[i][j]=FALSE;
          }
        }
      }
      break;
    case SAMPLE:
      if(!ComputeResidenceTimes[CurrentSystem]) return;

      ResidenceTimesFractionCounts[CurrentSystem]+=1.0;
      for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
      {
        CurrentComponent=Adsorbates[CurrentSystem][j].Type;
        StartingBead=Components[CurrentComponent].StartingBead;

        if(BlockedPocket(Adsorbates[CurrentSystem][j].Atoms[StartingBead].Position))
          ResidenceTimesFractionAdsorbates[CurrentSystem][j][Block]+=1.0;

        // outside blocking area and residence-status shows 'true' -> leaving a pocket
        if((!BlockedPocket(Adsorbates[CurrentSystem][j].Atoms[StartingBead].Position))&&(ResidenceStatusAdsorbates[CurrentSystem][j]==TRUE))
        {
          ResidenceStatusAdsorbates[CurrentSystem][j]=FALSE;
          time_difference=(REAL)(CurrentCycle-ResidenceOriginAdsorbates[CurrentSystem][j])*DeltaT;

          deltaR=RangeResidenceTimes[CurrentSystem]/ResidenceTimesHistogramSize[CurrentSystem];
          index=(int)(time_difference/deltaR);
          if((index>=0)&&(index<ResidenceTimesHistogramSize[CurrentSystem]))
            ResidenceTimesHistogram[CurrentSystem][index]+=1.0;
        }
        // inside blocking area and residence-status shows 'false' -> entering a pocket
        if((BlockedPocket(Adsorbates[CurrentSystem][j].Atoms[StartingBead].Position))&&(ResidenceStatusAdsorbates[CurrentSystem][j]==FALSE))
        {
          ResidenceStatusAdsorbates[CurrentSystem][j]=TRUE;
          ResidenceOriginAdsorbates[CurrentSystem][j]=CurrentCycle;
        }
      }

      for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
      {
        CurrentComponent=Cations[CurrentSystem][j].Type;
        StartingBead=Components[CurrentComponent].StartingBead;

        if(BlockedPocket(Cations[CurrentSystem][j].Atoms[StartingBead].Position))
          ResidenceTimesFractionCations[CurrentSystem][j][Block]+=1.0;

        // outside blocking area and residence-status shows 'true' -> leaving a pocket
        if((!BlockedPocket(Cations[CurrentSystem][j].Atoms[StartingBead].Position))&&(ResidenceStatusCations[CurrentSystem][j]==TRUE))
        {
          ResidenceStatusCations[CurrentSystem][j]=FALSE;
          time_difference=(REAL)(CurrentCycle-ResidenceOriginCations[CurrentSystem][j])*DeltaT;

          deltaR=RangeResidenceTimes[CurrentSystem]/ResidenceTimesHistogramSize[CurrentSystem];
          index=(int)(time_difference/deltaR);
          if((index>=0)&&(index<ResidenceTimesHistogramSize[CurrentSystem]))
            ResidenceTimesHistogram[CurrentSystem][index]+=1.0;
        }
        // inside blocking area and residence-status shows 'false' -> entering a pocket
        if((BlockedPocket(Cations[CurrentSystem][j].Atoms[StartingBead].Position))&&(ResidenceStatusCations[CurrentSystem][j]==FALSE))
        {
          ResidenceStatusCations[CurrentSystem][j]=TRUE;
          ResidenceOriginCations[CurrentSystem][j]=CurrentCycle;
        }
      }

      break;
    case PRINT:
      if((!ComputeResidenceTimes[CurrentSystem])||(CurrentCycle%WriteResidenceTimesEvery[CurrentSystem]!=0)) return;

      deltaR=RangeResidenceTimes[CurrentSystem]/ResidenceTimesHistogramSize[CurrentSystem];
      mkdir("ResidenceTimesHistogram",S_IRWXU);
      sprintf(buffer,"ResidenceTimesHistogram/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      sprintf(buffer,"ResidenceTimesHistogram/System_%d/HistogramResidenceTimes%s.dat",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"# column 1: t [ps]\n");
      fprintf(FilePtr,"# column 2: Residence Times\n");
      for(i=0;i<ResidenceTimesHistogramSize[CurrentSystem];i++)
        fprintf(FilePtr,"%g %g\n",deltaR*(i+0.5),ResidenceTimesHistogram[CurrentSystem][i]);
      fclose(FilePtr);

      sprintf(buffer,"ResidenceTimesHistogram/System_%d/FractionResidenceTimes%s.dat",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"# fraction of time spent inside the 'blocking-pockets'\n");
      if(NumberOfAdsorbateMolecules[CurrentSystem]>0)
        fprintf(FilePtr,"# adsorbates\n");
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        temp=0.0;
        for(k=0;k<NR_BLOCKS;k++)
          temp+=ResidenceTimesFractionAdsorbates[CurrentSystem][i][k];
        fprintf(FilePtr,"%d %g\n",i,temp/ResidenceTimesFractionCounts[CurrentSystem]);
      }

      if(NumberOfCationMolecules[CurrentSystem]>0)
        fprintf(FilePtr,"# cations\n");
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        temp=0.0;
        for(k=0;k<NR_BLOCKS;k++)
          temp+=ResidenceTimesFractionCations[CurrentSystem][i][k];
        fprintf(FilePtr,"%d %g\n",i,temp/ResidenceTimesFractionCounts[CurrentSystem]);
      }
      fclose(FilePtr);
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        free(ResidenceOriginAdsorbates[i]);
        free(ResidenceOriginCations[i]);
        free(ResidenceStatusAdsorbates[i]);
        free(ResidenceStatusCations[i]);
        free(ResidenceTimesHistogram[i]);
        free(ResidenceTimesFractionAdsorbates[i]);
        free(ResidenceTimesFractionCations[i]);
      }
      free(ResidenceOriginAdsorbates);
      free(ResidenceOriginCations);
      free(ResidenceStatusAdsorbates);
      free(ResidenceStatusCations);
      free(ResidenceTimesHistogram);
      free(ResidenceTimesFractionAdsorbates);
      free(ResidenceTimesFractionCations);
      free(ResidenceTimesFractionCounts);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleDistanceHistogram                                                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the histograms of the distance between 2 selected atoms.                         *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleDistanceHistogram(int Switch)
{
  int i,j,k,index;
  VECTOR posA,posB,dr;
  REAL r,norm;
  FILE *FilePtr;
  char buffer[256];

  switch(Switch)
  {
    case ALLOCATE:
      DistanceHistograms=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
        DistanceHistograms[i]=(REAL**)calloc(NumberOfDistanceHistogramDefinitions[i],sizeof(REAL*));
        for(j=0;j<NumberOfDistanceHistogramDefinitions[i];j++)
          DistanceHistograms[i][j]=(REAL*)calloc(NumberOfElementsDistanceHistogram,sizeof(REAL));
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      if(!ComputeDistanceHistograms[CurrentSystem]) return;

      for(i=0;i<NumberOfDistanceHistogramDefinitions[CurrentSystem];i++)
      {
        posA=DistanceHistogramPairs[CurrentSystem][i][0]->Position;
        posB=DistanceHistogramPairs[CurrentSystem][i][1]->Position;
        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

        index=(int)(r*(REAL)NumberOfElementsDistanceHistogram/MaxRangeDistanceHistogram);
        if(index>=0&&index<NumberOfElementsDistanceHistogram)
          DistanceHistograms[CurrentSystem][i][index]+=1.0;
      }
      break;
    case PRINT:
      if((!ComputeDistanceHistograms[CurrentSystem])||(CurrentCycle%WriteDistanceHistogramsEvery[CurrentSystem]!=0)) return;

      mkdir("DistanceHistograms",S_IRWXU);
      sprintf(buffer,"DistanceHistograms/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(i=0;i<NumberOfDistanceHistogramDefinitions[CurrentSystem];i++)
      {
        norm=0.0;
        for(k=0;k<NumberOfElementsDistanceHistogram;k++)
          norm+=DistanceHistograms[CurrentSystem][i][k];
        norm*=(REAL)MaxRangeDistanceHistogram/(REAL)NumberOfElementsDistanceHistogram;

        sprintf(buffer,"DistanceHistograms/System_%d/Histogram_%d%s.dat",CurrentSystem,i,FileNameAppend);
        FilePtr=fopen(buffer,"w");

        switch(DistanceHistogramDefinitions[CurrentSystem][i][0][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Atom A: Framework %d atom %d\n",DistanceHistogramDefinitions[CurrentSystem][i][0][1],DistanceHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Atom A: Adsorbate molecule %d atom %d\n",DistanceHistogramDefinitions[CurrentSystem][i][0][1],DistanceHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Atom A: Cation molecule %d atom %d\n",DistanceHistogramDefinitions[CurrentSystem][i][0][1],DistanceHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          default:
            break;
        }
        switch(DistanceHistogramDefinitions[CurrentSystem][i][1][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Atom B: Framework %d atom %d\n",DistanceHistogramDefinitions[CurrentSystem][i][1][1],DistanceHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Atom B: Adsorbate molecule %d atom %d\n",DistanceHistogramDefinitions[CurrentSystem][i][1][1],DistanceHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Atom B: Cation molecule %d atom %d\n",DistanceHistogramDefinitions[CurrentSystem][i][1][1],DistanceHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          default:
            break;
        }


        for(k=0;k<NumberOfElementsDistanceHistogram;k++)
        {
          r=(REAL)k*((REAL)MaxRangeDistanceHistogram/(REAL)NumberOfElementsDistanceHistogram);
          if(DistanceHistograms[CurrentSystem][i][k]>0.0)
            fprintf(FilePtr,"%g %g\n",(double)r,(double)(DistanceHistograms[CurrentSystem][i][k]/norm));
        }
        fclose(FilePtr);
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        for(j=0;j<NumberOfDistanceHistogramDefinitions[i];j++)
          free(DistanceHistograms[i][j]);
        free(DistanceHistograms[i]);
      }
      free(DistanceHistograms);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleBendAngleHistogram                                                                 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the histograms of the bend angle between 3 selected atoms.                       *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleBendAngleHistogram(int Switch)
{
  int i,j,k,index;
  VECTOR posA,posB,posC;
  REAL theta,norm;
  FILE *FilePtr;
  char buffer[256];

  switch(Switch)
  {
    case ALLOCATE:
      BendAngleHistograms=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
        BendAngleHistograms[i]=(REAL**)calloc(NumberOfBendAngleHistogramDefinitions[i],sizeof(REAL*));
        for(j=0;j<NumberOfBendAngleHistogramDefinitions[i];j++)
          BendAngleHistograms[i][j]=(REAL*)calloc(NumberOfElementsBendAngleHistogram,sizeof(REAL));
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      if(!ComputeBendAngleHistograms[CurrentSystem]) return;

      for(i=0;i<NumberOfBendAngleHistogramDefinitions[CurrentSystem];i++)
      {
        posA=BendAngleHistogramPairs[CurrentSystem][i][0]->Position;
        posB=BendAngleHistogramPairs[CurrentSystem][i][1]->Position;
        posC=BendAngleHistogramPairs[CurrentSystem][i][2]->Position;
        theta=ReturnBendAngle(posA,posB,posC)*RAD2DEG;

        index=(int)(theta*(REAL)NumberOfElementsBendAngleHistogram/MaxRangeBendAngleHistogram);
        if(index>=0&&index<NumberOfElementsBendAngleHistogram)
          BendAngleHistograms[CurrentSystem][i][index]+=1.0;
      }
      break;
    case PRINT:
      if((!ComputeBendAngleHistograms[CurrentSystem])||(CurrentCycle%WriteBendAngleHistogramsEvery[CurrentSystem]!=0)) return;

      mkdir("BendAngleHistograms",S_IRWXU);
      sprintf(buffer,"BendAngleHistograms/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(i=0;i<NumberOfBendAngleHistogramDefinitions[CurrentSystem];i++)
      {
        norm=0.0;
        for(k=0;k<NumberOfElementsBendAngleHistogram;k++)
          norm+=BendAngleHistograms[CurrentSystem][i][k];
        norm*=(REAL)MaxRangeBendAngleHistogram/(REAL)NumberOfElementsBendAngleHistogram;

        sprintf(buffer,"BendAngleHistograms/System_%d/Histogram_%d%s.dat",CurrentSystem,i,FileNameAppend);
        FilePtr=fopen(buffer,"w");

        switch(BendAngleHistogramDefinitions[CurrentSystem][i][0][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Atom A: Framework %d atom %d\n",BendAngleHistogramDefinitions[CurrentSystem][i][0][1],BendAngleHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Atom A: Adsorbate molecule %d atom %d\n",BendAngleHistogramDefinitions[CurrentSystem][i][0][1],BendAngleHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Atom A: Cation molecule %d atom %d\n",BendAngleHistogramDefinitions[CurrentSystem][i][0][1],BendAngleHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          default:
            break;
        }
        switch(BendAngleHistogramDefinitions[CurrentSystem][i][1][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Atom B: Framework %d atom %d\n",BendAngleHistogramDefinitions[CurrentSystem][i][1][1],BendAngleHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Atom B: Adsorbate molecule %d atom %d\n",BendAngleHistogramDefinitions[CurrentSystem][i][1][1],BendAngleHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Atom B: Cation molecule %d atom %d\n",BendAngleHistogramDefinitions[CurrentSystem][i][1][1],BendAngleHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          default:
            break;
        }
        switch(BendAngleHistogramDefinitions[CurrentSystem][i][2][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Atom C: Framework %d atom %d\n",BendAngleHistogramDefinitions[CurrentSystem][i][2][1],BendAngleHistogramDefinitions[CurrentSystem][i][2][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Atom C: Adsorbate molecule %d atom %d\n",BendAngleHistogramDefinitions[CurrentSystem][i][2][1],BendAngleHistogramDefinitions[CurrentSystem][i][2][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Atom C: Cation molecule %d atom %d\n",BendAngleHistogramDefinitions[CurrentSystem][i][2][1],BendAngleHistogramDefinitions[CurrentSystem][i][2][2]);
            break;
          default:
            break;
        }

        for(k=0;k<NumberOfElementsBendAngleHistogram;k++)
        {
          theta=(REAL)k*((REAL)MaxRangeBendAngleHistogram/(REAL)NumberOfElementsBendAngleHistogram);
          if(BendAngleHistograms[CurrentSystem][i][k]>0.0)
            fprintf(FilePtr,"%g %g\n",(double)theta,(double)(BendAngleHistograms[CurrentSystem][i][k]/norm));
        }
        fclose(FilePtr);
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        for(j=0;j<NumberOfBendAngleHistogramDefinitions[i];j++)
          free(BendAngleHistograms[i][j]);
        free(BendAngleHistograms[i]);
      }
      free(BendAngleHistograms);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleDihedralAngleHistogram                                                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the histograms of the dihedral angle between 4 selected atoms.                   *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleDihedralAngleHistogram(int Switch)
{
  int i,j,k,index;
  VECTOR posA,posB,posC,posD;
  REAL phi,norm;
  FILE *FilePtr;
  char buffer[256];

  switch(Switch)
  {
    case ALLOCATE:
      DihedralAngleHistograms=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
        DihedralAngleHistograms[i]=(REAL**)calloc(NumberOfDihedralAngleHistogramDefinitions[i],sizeof(REAL*));
        for(j=0;j<NumberOfDihedralAngleHistogramDefinitions[i];j++)
          DihedralAngleHistograms[i][j]=(REAL*)calloc(NumberOfElementsDihedralAngleHistogram,sizeof(REAL));
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      if(!ComputeDihedralAngleHistograms[CurrentSystem]) return;

      for(i=0;i<NumberOfDihedralAngleHistogramDefinitions[CurrentSystem];i++)
      {
        posA=DihedralAngleHistogramPairs[CurrentSystem][i][0]->Position;
        posB=DihedralAngleHistogramPairs[CurrentSystem][i][1]->Position;
        posC=DihedralAngleHistogramPairs[CurrentSystem][i][2]->Position;
        posD=DihedralAngleHistogramPairs[CurrentSystem][i][3]->Position;
        phi=ReturnDihedralAngle(posA,posB,posC,posD)*RAD2DEG;
        if(phi<0.0) phi+=360.0;

        index=(int)(phi*(REAL)NumberOfElementsDihedralAngleHistogram/MaxRangeDihedralAngleHistogram);
        if(index>=0&&index<NumberOfElementsDihedralAngleHistogram)
          DihedralAngleHistograms[CurrentSystem][i][index]+=1.0;
      }
      break;
    case PRINT:
      if((!ComputeDihedralAngleHistograms[CurrentSystem])||(CurrentCycle%WriteDihedralAngleHistogramsEvery[CurrentSystem]!=0)) return;

      mkdir("DihedralAngleHistograms",S_IRWXU);
      sprintf(buffer,"DihedralAngleHistograms/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(i=0;i<NumberOfDihedralAngleHistogramDefinitions[CurrentSystem];i++)
      {
        norm=0.0;
        for(k=0;k<NumberOfElementsDihedralAngleHistogram;k++)
          norm+=DihedralAngleHistograms[CurrentSystem][i][k];
        norm*=(REAL)MaxRangeDihedralAngleHistogram/(REAL)NumberOfElementsDihedralAngleHistogram;

        sprintf(buffer,"DihedralAngleHistograms/System_%d/Histogram_%d%s.dat",CurrentSystem,i,FileNameAppend);
        FilePtr=fopen(buffer,"w");

        switch(DihedralAngleHistogramDefinitions[CurrentSystem][i][0][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Atom A: Framework %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][0][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Atom A: Adsorbate molecule %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][0][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Atom A: Cation molecule %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][0][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          default:
            break;
        }
        switch(DihedralAngleHistogramDefinitions[CurrentSystem][i][1][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Atom B: Framework %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][1][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Atom B: Adsorbate molecule %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][1][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Atom B: Cation molecule %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][1][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          default:
            break;
        }
        switch(DihedralAngleHistogramDefinitions[CurrentSystem][i][2][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Atom C: Framework %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][2][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][2][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Atom C: Adsorbate molecule %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][2][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][2][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Atom C: Cation molecule %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][2][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][2][2]);
            break;
          default:
            break;
        }
        switch(DihedralAngleHistogramDefinitions[CurrentSystem][i][3][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Atom D: Framework %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][3][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][3][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Atom D: Adsorbate molecule %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][3][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][3][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Atom D: Cation molecule %d atom %d\n",DihedralAngleHistogramDefinitions[CurrentSystem][i][3][1],DihedralAngleHistogramDefinitions[CurrentSystem][i][3][2]);
            break;
          default:
            break;
        }


        for(k=0;k<NumberOfElementsDihedralAngleHistogram;k++)
        {
          phi=(REAL)k*((REAL)MaxRangeDihedralAngleHistogram/(REAL)NumberOfElementsDihedralAngleHistogram);
          if(DihedralAngleHistograms[CurrentSystem][i][k]>0.0)
            fprintf(FilePtr,"%g %g\n",(double)phi,(double)(DihedralAngleHistograms[CurrentSystem][i][k]/norm));
        }
        fclose(FilePtr);
      }

      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        for(j=0;j<NumberOfDihedralAngleHistogramDefinitions[i];j++)
          free(DihedralAngleHistograms[i][j]);
        free(DihedralAngleHistograms[i]);
      }
      free(DihedralAngleHistograms);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleAngleBetweenPlanesHistogram                                                        *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the histograms of the angle between two planes (each formed by 3 chosen atoms).  *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleAngleBetweenPlanesHistogram(int Switch)
{
  int i,j,k,index;
  VECTOR posA1,posB1,posC1;
  VECTOR posA2,posB2,posC2;
  VECTOR Rba,Rca;
  REAL theta,norm,length,CosTheta;
  VECTOR cross_product1,cross_product2;
  FILE *FilePtr;
  char buffer[256];

  switch(Switch)
  {
    case ALLOCATE:
      AngleBetweenPlanesHistograms=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
        AngleBetweenPlanesHistograms[i]=(REAL**)calloc(NumberOfAngleBetweenPlanesHistogramDefinitions[i],sizeof(REAL*));
        for(j=0;j<NumberOfAngleBetweenPlanesHistogramDefinitions[i];j++)
          AngleBetweenPlanesHistograms[i][j]=(REAL*)calloc(NumberOfElementsAngleBetweenPlanesHistogram,sizeof(REAL));
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      if(!ComputeAngleBetweenPlanesHistograms[CurrentSystem]) return;

      for(i=0;i<NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem];i++)
      {
        posA1=AngleBetweenPlanesHistogramPairs[CurrentSystem][i][0]->Position;
        posB1=AngleBetweenPlanesHistogramPairs[CurrentSystem][i][1]->Position;
        posC1=AngleBetweenPlanesHistogramPairs[CurrentSystem][i][2]->Position;

        posA2=AngleBetweenPlanesHistogramPairs[CurrentSystem][i][3]->Position;
        posB2=AngleBetweenPlanesHistogramPairs[CurrentSystem][i][4]->Position;
        posC2=AngleBetweenPlanesHistogramPairs[CurrentSystem][i][5]->Position;

/*
        Va=RandomNumberOnUnitSphere();
        Vb=Va;

        dot_product=Va.x*(posB1.x-posA1.x)+Va.y*(posB1.y-posA1.y)+Va.z*(posB1.z-posA1.z);
        Va.x-=dot_product*(posB1.x-posA1.x)/(SQR(posB1.x-posA1.x)+SQR(posB1.y-posA1.y)+SQR(posB1.z-posA1.z));
        Va.y-=dot_product*(posB1.y-posA1.y)/(SQR(posB1.x-posA1.x)+SQR(posB1.y-posA1.y)+SQR(posB1.z-posA1.z));
        Va.z-=dot_product*(posB1.z-posA1.z)/(SQR(posB1.x-posA1.x)+SQR(posB1.y-posA1.y)+SQR(posB1.z-posA1.z));

        dot_product=Vb.x*(posB2.x-posA2.x)+Vb.y*(posB2.y-posA2.y)+Vb.z*(posB2.z-posA2.z);
        Vb.x-=dot_product*(posB2.x-posA2.x)/(SQR(posB2.x-posA2.x)+SQR(posB2.y-posA1.y)+SQR(posB2.z-posA2.z));
        Vb.y-=dot_product*(posB2.y-posA2.y)/(SQR(posB2.x-posA2.x)+SQR(posB2.y-posA1.y)+SQR(posB2.z-posA2.z));
        Vb.z-=dot_product*(posB2.z-posA2.z)/(SQR(posB2.x-posA2.x)+SQR(posB2.y-posA1.y)+SQR(posB2.z-posA2.z));

        dot_product=Va.x*(posC1.x-posA1.x)+Va.y*(posC1.y-posA1.y)+Va.z*(posC1.z-posA1.z);
        Va.x-=dot_product*(posC1.x-posA1.x)/(SQR(posC1.x-posA1.x)+SQR(posC1.y-posA1.y)+SQR(posC1.z-posA1.z));
        Va.y-=dot_product*(posC1.y-posA1.y)/(SQR(posC1.x-posA1.x)+SQR(posC1.y-posA1.y)+SQR(posC1.z-posA1.z));
        Va.z-=dot_product*(posC1.z-posA1.z)/(SQR(posC1.x-posA1.x)+SQR(posC1.y-posA1.y)+SQR(posC1.z-posA1.z));

        dot_product=Vb.x*(posC2.x-posA2.x)+Vb.y*(posC2.y-posA2.y)+Vb.z*(posC2.z-posA2.z);
        Vb.x-=dot_product*(posC2.x-posA2.x)/(SQR(posC2.x-posA2.x)+SQR(posC2.y-posA2.y)+SQR(posC2.z-posA2.z));
        Vb.y-=dot_product*(posC2.y-posA2.y)/(SQR(posC2.x-posA2.x)+SQR(posC2.y-posA2.y)+SQR(posC2.z-posA2.z));
        Vb.z-=dot_product*(posC2.z-posA2.z)/(SQR(posC2.x-posA2.x)+SQR(posC2.y-posA2.y)+SQR(posC2.z-posA2.z));

        theta=acos((Va.x*Vb.x+Va.y*Vb.y+Va.z*Vb.z)/(sqrt(SQR(Va.x)+SQR(Va.y)+SQR(Va.z))*sqrt(SQR(Vb.x)+SQR(Vb.y)+SQR(Vb.z))))*RAD2DEG;
*/

        Rba.x=(posB1.x-posA1.x);
        Rba.y=(posB1.y-posA1.y);
        Rba.z=(posB1.z-posA1.z);
        Rca.x=(posC1.x-posA1.x);
        Rca.y=(posC1.y-posA1.y);
        Rca.z=(posC1.z-posA1.z);
        cross_product1.x=Rba.y*Rca.z-Rba.z*Rca.y;
        cross_product1.y=Rba.z*Rca.x-Rba.x*Rca.z;
        cross_product1.z=Rba.x*Rca.y-Rba.y*Rca.x;
        length=sqrt(SQR(cross_product1.x)+SQR(cross_product1.y)+SQR(cross_product1.z));
        cross_product1.x/=length;
        cross_product1.y/=length;
        cross_product1.z/=length;

        Rba.x=(posB2.x-posA2.x);
        Rba.y=(posB2.y-posA2.y);
        Rba.z=(posB2.z-posA2.z);
        Rca.x=(posC2.x-posA2.x);
        Rca.y=(posC2.y-posA2.y);
        Rca.z=(posC2.z-posA2.z);
        cross_product2.x=Rba.y*Rca.z-Rba.z*Rca.y;
        cross_product2.y=Rba.z*Rca.x-Rba.x*Rca.z;
        cross_product2.z=Rba.x*Rca.y-Rba.y*Rca.x;
        length=sqrt(SQR(cross_product2.x)+SQR(cross_product2.y)+SQR(cross_product2.z));
        cross_product2.x/=length;
        cross_product2.y/=length;
        cross_product2.z/=length;

        CosTheta=cross_product1.x*cross_product2.x+cross_product1.y*cross_product2.y+cross_product1.z*cross_product2.z;
        theta=SIGN(acos(CosTheta)*RAD2DEG,(posB1.x-posA1.x)*cross_product2.x+(posB1.y-posA1.y)*cross_product2.y+(posB1.z-posA1.z)*cross_product2.z);
        if(theta<0.0) theta+=360.0;

        index=(int)(theta*(REAL)NumberOfElementsAngleBetweenPlanesHistogram/MaxRangeAngleBetweenPlanesHistogram);
        if(index>=0&&index<NumberOfElementsAngleBetweenPlanesHistogram)
          AngleBetweenPlanesHistograms[CurrentSystem][i][index]+=1.0;
      }
      break;
    case PRINT:
      if((!ComputeAngleBetweenPlanesHistograms[CurrentSystem])||(CurrentCycle%WriteAngleBetweenPlanesHistogramsEvery[CurrentSystem]!=0)) return;

      mkdir("AngleBetweenPlanesHistograms",S_IRWXU);
      sprintf(buffer,"AngleBetweenPlanesHistograms/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(i=0;i<NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem];i++)
      {
        norm=0.0;
        for(k=0;k<NumberOfElementsAngleBetweenPlanesHistogram;k++)
          norm+=AngleBetweenPlanesHistograms[CurrentSystem][i][k];
        norm*=(REAL)MaxRangeAngleBetweenPlanesHistogram/(REAL)NumberOfElementsAngleBetweenPlanesHistogram;

        sprintf(buffer,"AngleBetweenPlanesHistograms/System_%d/Histogram_%d%s.dat",CurrentSystem,i,FileNameAppend);
        FilePtr=fopen(buffer,"w");

        switch(AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][0][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Plane 1 Atom A: Framework %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][0][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Plane 1 Atom A: Adsorbate molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][0][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Plane 1 Atom A: Cation molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][0][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][0][2]);
            break;
          default:
            break;
        }
        switch(AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][1][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Plane 1 Atom B: Framework %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][1][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Plane 1 Atom B: Adsorbate molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][1][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Plane 1 Atom B: Cation molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][1][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][1][2]);
            break;
          default:
            break;
        }
        switch(AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][2][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Plane 1 Atom C: Framework %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][2][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][2][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Plane 1 Atom C: Adsorbate molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][2][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][2][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Plane 1 Atom C: Cation molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][2][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][2][2]);
            break;
          default:
            break;
        }

        switch(AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][3][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Plane 2 Atom A: Framework %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][3][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][3][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Plane 2 Atom A: Adsorbate molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][3][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][3][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Plane 2 Atom A: Cation molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][3][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][3][2]);
            break;
          default:
            break;
        }
        switch(AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][4][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Plane 2 Atom B: Framework %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][4][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][4][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Plane 2 Atom B: Adsorbate molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][4][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][4][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Plane 2 Atom B: Cation molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][4][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][4][2]);
            break;
          default:
            break;
        }
        switch(AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][5][0])
        {
          case FRAMEWORK:
            fprintf(FilePtr,"# Plane 2 Atom C: Framework %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][5][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][5][2]);
            break;
          case ADSORBATE:
            fprintf(FilePtr,"# Plane 2 Atom C: Adsorbate molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][5][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][5][2]);
            break;
          case CATION:
            fprintf(FilePtr,"# Plane 2 Atom C: Cation molecule %d atom %d\n",AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][5][1],AngleBetweenPlanesHistogramDefinitions[CurrentSystem][i][5][2]);
            break;
          default:
            break;
        }
        for(k=0;k<NumberOfElementsAngleBetweenPlanesHistogram;k++)
        {
          theta=(REAL)k*((REAL)MaxRangeAngleBetweenPlanesHistogram/(REAL)NumberOfElementsAngleBetweenPlanesHistogram);
          if(AngleBetweenPlanesHistograms[CurrentSystem][i][k]>0.0)
            fprintf(FilePtr,"%g %g\n",(double)theta,(double)(AngleBetweenPlanesHistograms[CurrentSystem][i][k]/norm));
        }
        fclose(FilePtr);
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        for(j=0;j<NumberOfAngleBetweenPlanesHistogramDefinitions[i];j++)
          free(AngleBetweenPlanesHistograms[i][j]);
        free(AngleBetweenPlanesHistograms[i]);
      }
      free(AngleBetweenPlanesHistograms);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleMoleculePropertyHistogram                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the molecular properties (bond distance, bend angle, dihedral angle).            *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleMoleculePropertyHistogram(int Switch)
{
  int i,j,k,index,type;
  int NumberOfBonds,NumberOfUreyBradleys,NumberOfBends,NumberOfTorsions,Type;
  REAL r,theta,phi,norm;
  FILE *FilePtr;
  char buffer[256];

  switch(Switch)
  {
    case ALLOCATE:
      BondLengthHistogram=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      UreyBradleyLengthHistogram=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      BendAngleHistogram=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      TorsionAngleHistogram=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      TorsionConformationHistogram=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));

      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeMoleculeProperties[i])
        {
          BondLengthHistogram[i]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));
          UreyBradleyLengthHistogram[i]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));
          BendAngleHistogram[i]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));
          TorsionAngleHistogram[i]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));
          TorsionConformationHistogram[i]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));
          for(j=0;j<NumberOfComponents;j++)
          {
            BondLengthHistogram[i][j]=(REAL**)calloc(MaxNumberOfBonds,sizeof(REAL*));
            for(k=0;k<MaxNumberOfBonds;k++)
              BondLengthHistogram[i][j][k]=(REAL*)calloc(BondLengthHistogramSize[i],sizeof(REAL));

            UreyBradleyLengthHistogram[i][j]=(REAL**)calloc(MaxNumberOfUreyBradleys,sizeof(REAL*));
            for(k=0;k<MaxNumberOfUreyBradleys;k++)
              UreyBradleyLengthHistogram[i][j][k]=(REAL*)calloc(BondLengthHistogramSize[i],sizeof(REAL));

            BendAngleHistogram[i][j]=(REAL**)calloc(MaxNumberOfBends,sizeof(REAL*));
            for(k=0;k<MaxNumberOfBends;k++)
              BendAngleHistogram[i][j][k]=(REAL*)calloc(BendAngleHistogramSize[i],sizeof(REAL));

            TorsionAngleHistogram[i][j]=(REAL**)calloc(MaxNumberOfTorsions,sizeof(REAL*));
            for(k=0;k<MaxNumberOfTorsions;k++)
              TorsionAngleHistogram[i][j][k]=(REAL*)calloc(DihedralHistogramSize[i],sizeof(REAL));

            TorsionConformationHistogram[i][j]=(REAL**)calloc(MaxNumberOfTorsions,sizeof(REAL*));
            for(k=0;k<MaxNumberOfTorsions;k++)
              TorsionConformationHistogram[i][j][k]=(REAL*)calloc(6,sizeof(REAL));
          }
        }
      }

      FrameworkBondLengthHistogram=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      FrameworkUreyBradleyLengthHistogram=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      FrameworkBendAngleHistogram=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      FrameworkTorsionAngleHistogram=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

      FrameworkAverageBondLength=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
      FrameworkBondLengthCount=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
      FrameworkAverageBendAngle=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
      FrameworkBendAngleCount=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
      FrameworkAverageTorsionAngle=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
      FrameworkTorsionAngleCount=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeMoleculeProperties[i]&&(Framework[i].FrameworkModel==FLEXIBLE))
        {
          if(Framework[i].FrameworkModel==FLEXIBLE)
          FrameworkBondLengthHistogram[i]=(REAL**)calloc(Framework[i].NumberOfBonds[0],sizeof(REAL*));
          for(j=0;j<Framework[i].NumberOfBonds[0];j++)
            FrameworkBondLengthHistogram[i][j]=(REAL*)calloc(BondLengthHistogramSize[i],sizeof(REAL));

          FrameworkUreyBradleyLengthHistogram[i]=(REAL**)calloc(Framework[i].NumberOfUreyBradleys[0],sizeof(REAL*));
          for(j=0;j<Framework[i].NumberOfUreyBradleys[0];j++)
            FrameworkUreyBradleyLengthHistogram[i][j]=(REAL*)calloc(BondLengthHistogramSize[i],sizeof(REAL));

          FrameworkBendAngleHistogram[i]=(REAL**)calloc(Framework[i].NumberOfBends[0],sizeof(REAL*));
          for(j=0;j<Framework[i].NumberOfBends[0];j++)
            FrameworkBendAngleHistogram[i][j]=(REAL*)calloc(BendAngleHistogramSize[i],sizeof(REAL));

          FrameworkTorsionAngleHistogram[i]=(REAL**)calloc(Framework[i].NumberOfTorsions[0],sizeof(REAL*));
          for(j=0;j<Framework[i].NumberOfTorsions[0];j++)
            FrameworkTorsionAngleHistogram[i][j]=(REAL*)calloc(DihedralHistogramSize[i],sizeof(REAL));

          FrameworkAverageBondLength[i]=(REAL*)calloc(Framework[i].NumberOfBonds[0],sizeof(REAL));
          FrameworkBondLengthCount[i]=(REAL*)calloc(Framework[i].NumberOfBonds[0],sizeof(REAL));
          FrameworkAverageBendAngle[i]=(REAL*)calloc(Framework[i].NumberOfBends[0],sizeof(REAL));
          FrameworkBendAngleCount[i]=(REAL*)calloc(Framework[i].NumberOfBends[0],sizeof(REAL));
          FrameworkAverageTorsionAngle[i]=(REAL*)calloc(Framework[i].NumberOfTorsions[0],sizeof(REAL));
          FrameworkTorsionAngleCount[i]=(REAL*)calloc(Framework[i].NumberOfTorsions[0],sizeof(REAL));
        }
      }
      break;
    case INITIALIZE:
     break;
    case SAMPLE:
      if(!ComputeMoleculeProperties[CurrentSystem]) return;

      CurrentFramework=0;

      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfBonds[CurrentFramework];i++)
        {
          r=ComputeBondDistanceFramework(i)/BondLengthRange[CurrentSystem];

          index=(int)(r*(REAL)BondLengthHistogramSize[CurrentSystem]);
          if(index>=0&&index<BondLengthHistogramSize[CurrentSystem])
          {
            FrameworkBondLengthHistogram[CurrentSystem][i][index]+=1.0;
            FrameworkAverageBondLength[CurrentSystem][i]+=r;
            FrameworkBondLengthCount[CurrentSystem][i]+=1.0;
          }
        }

        for(i=0;i<Framework[CurrentSystem].NumberOfBends[CurrentFramework];i++)
        {
          theta=ComputeBendAngleFramework(i)/BendAngleRange[CurrentSystem];

          index=(int)(theta*(REAL)BendAngleHistogramSize[CurrentSystem]);
          if(index>=0&&index<BendAngleHistogramSize[CurrentSystem])
          {
            FrameworkBendAngleHistogram[CurrentSystem][i][index]+=1.0;
            FrameworkAverageBendAngle[CurrentSystem][i]+=theta;
            FrameworkBendAngleCount[CurrentSystem][i]+=1.0;
          }
        }

        for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[CurrentFramework];i++)
        {
          phi=ComputeTorsionAngleFramework(i)/DihedralRange[CurrentSystem];

          index=(int)(phi*(REAL)DihedralHistogramSize[CurrentSystem]);
          if(index>=0&&index<DihedralHistogramSize[CurrentSystem])
          {
            FrameworkTorsionAngleHistogram[CurrentSystem][i][index]+=1.0;
            FrameworkAverageTorsionAngle[CurrentSystem][i]+=phi;
            FrameworkTorsionAngleCount[CurrentSystem][i]+=1.0;
          }
        }
      }

      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        Type=Adsorbates[CurrentSystem][i].Type;
        NumberOfBonds=Components[Type].NumberOfBonds;
        NumberOfUreyBradleys=Components[Type].NumberOfUreyBradleys;
        NumberOfBends=Components[Type].NumberOfBends;
        NumberOfTorsions=Components[Type].NumberOfTorsions;

        for(j=0;j<NumberOfBonds;j++)
        {
          r=ComputeBondDistanceAdsorbate(i,j)/BondLengthRange[CurrentSystem];
          index=(int)(r*(REAL)BondLengthHistogramSize[CurrentSystem]);
          if(index>=0&&index<BondLengthHistogramSize[CurrentSystem])
            BondLengthHistogram[CurrentSystem][Type][j][index]+=1.0;
        }

        for(j=0;j<NumberOfUreyBradleys;j++)
        {
          r=ComputeUreyBradleyDistanceAdsorbate(i,j)/BondLengthRange[CurrentSystem];
          index=(int)(r*(REAL)BondLengthHistogramSize[CurrentSystem]);
          if(index>=0&&index<BondLengthHistogramSize[CurrentSystem])
            UreyBradleyLengthHistogram[CurrentSystem][Type][j][index]+=1.0;
        }

        for(j=0;j<NumberOfBends;j++)
        {
          theta=ComputeBendAngleAdsorbate(i,j)/BendAngleRange[CurrentSystem];
          index=(int)(theta*(REAL)BendAngleHistogramSize[CurrentSystem]);
          if(index>=0&&index<BendAngleHistogramSize[CurrentSystem])
            BendAngleHistogram[CurrentSystem][Type][j][index]+=1.0;
        }

        for(j=0;j<NumberOfTorsions;j++)
        {
          theta=ComputeTorsionAngleAdsorbate(i,j);
          type=ReturnTorsionConformation(theta);
          TorsionConformationHistogram[CurrentSystem][Type][j][type]+=1.0;
          index=(int)((theta/DihedralRange[CurrentSystem])*(REAL)DihedralHistogramSize[CurrentSystem]);
          if(index>=0&&index<DihedralHistogramSize[CurrentSystem])
            TorsionAngleHistogram[CurrentSystem][Type][j][index]+=1.0;
        }
      }

      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        Type=Cations[CurrentSystem][i].Type;
        NumberOfBonds=Components[Type].NumberOfBonds;
        NumberOfUreyBradleys=Components[Type].NumberOfUreyBradleys;
        NumberOfBends=Components[Type].NumberOfBends;
        NumberOfTorsions=Components[Type].NumberOfTorsions;

        for(j=0;j<NumberOfBonds;j++)
        {
          r=ComputeBondDistanceCation(i,j)/BondLengthRange[CurrentSystem];
          index=(int)(r*(REAL)BondLengthHistogramSize[CurrentSystem]);
          if(index>=0&&index<BondLengthHistogramSize[CurrentSystem])
            BondLengthHistogram[CurrentSystem][Type][j][index]+=1.0;
        }

        for(j=0;j<NumberOfUreyBradleys;j++)
        {
          r=ComputeUreyBradleyDistanceCation(i,j)/BondLengthRange[CurrentSystem];
          index=(int)(r*(REAL)BondLengthHistogramSize[CurrentSystem]);
          if(index>=0&&index<BondLengthHistogramSize[CurrentSystem])
            UreyBradleyLengthHistogram[CurrentSystem][Type][j][index]+=1.0;
        }

        for(j=0;j<NumberOfBends;j++)
        {
          theta=ComputeBendAngleCation(i,j)/BendAngleRange[CurrentSystem];
          index=(int)(theta*(REAL)BendAngleHistogramSize[CurrentSystem]);
          if(index>=0&&index<BendAngleHistogramSize[CurrentSystem])
            BendAngleHistogram[CurrentSystem][Type][j][index]+=1.0;
        }
        for(j=0;j<NumberOfTorsions;j++)
        {
          theta=ComputeTorsionAngleCation(i,j);
          type=ReturnTorsionConformation(theta);
          TorsionConformationHistogram[CurrentSystem][Type][j][type]+=1.0;
          index=(int)((theta/DihedralRange[CurrentSystem])*(REAL)DihedralHistogramSize[CurrentSystem]);
          if(index>=0&&index<DihedralHistogramSize[CurrentSystem])
            TorsionAngleHistogram[CurrentSystem][Type][j][index]+=1.0;
        }
      }
      break;
    case PRINT:
      if((!ComputeMoleculeProperties[CurrentSystem])||(CurrentCycle%WriteMoleculePropertiesEvery[CurrentSystem]!=0)) return;
      mkdir("MoleculeProperties",S_IRWXU);
      sprintf(buffer,"MoleculeProperties/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfBonds[0];j++)
        {
          norm=0.0;
          for(k=0;k<BondLengthHistogramSize[CurrentSystem];k++)
            norm+=FrameworkBondLengthHistogram[CurrentSystem][j][k];
          norm*=BondLengthRange[CurrentSystem]/BondLengthHistogramSize[CurrentSystem];

          sprintf(buffer,"MoleculeProperties/System_%d/Histogram_Framework_BondLength_%d_%d_%d%s.dat",CurrentSystem,j,
              Framework[CurrentSystem].Bonds[0][j].A,
              Framework[CurrentSystem].Bonds[0][j].B,
              FileNameAppend);
          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# average bond length: %lf [A]\n",
             (double)(FrameworkAverageBondLength[CurrentSystem][j]/FrameworkBondLengthCount[CurrentSystem][j]));
          for(k=0;k<BondLengthHistogramSize[CurrentSystem];k++)
          {
            r=(REAL)k/(REAL)BondLengthHistogramSize[CurrentSystem];
            if(FrameworkBondLengthHistogram[CurrentSystem][j][k]>0.0)
              fprintf(FilePtr,"%lf %lf\n",(double)(r*BondLengthRange[CurrentSystem]),(double)(FrameworkBondLengthHistogram[CurrentSystem][j][k]/norm));
          }
          fclose(FilePtr);
        }

        for(j=0;j<Framework[CurrentSystem].NumberOfBends[0];j++)
        {
          norm=0.0;
          for(k=0;k<BendAngleHistogramSize[CurrentSystem];k++)
            norm+=FrameworkBendAngleHistogram[CurrentSystem][j][k];
          norm*=BendAngleRange[CurrentSystem]/BendAngleHistogramSize[CurrentSystem];

          sprintf(buffer,"MoleculeProperties/System_%d/Histogram_Framework_BendAngle_%d_%d_%d_%d%s.dat",CurrentSystem,j,
             Framework[CurrentSystem].Bends[0][j].A,
             Framework[CurrentSystem].Bends[0][j].B,
             Framework[CurrentSystem].Bends[0][j].C,
             FileNameAppend);
          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# average bend angle: %lf [degrees]\n",
            (double)(FrameworkAverageBendAngle[CurrentSystem][j]/FrameworkBendAngleCount[CurrentSystem][j]));
          for(k=0;k<BendAngleHistogramSize[CurrentSystem];k++)
          {
            theta=(REAL)k/(REAL)BendAngleHistogramSize[CurrentSystem];
            if(FrameworkBendAngleHistogram[CurrentSystem][j][k]>0.0)
              fprintf(FilePtr,"%lf %lf\n",(double)(theta*BendAngleRange[CurrentSystem]),(double)(FrameworkBendAngleHistogram[CurrentSystem][j][k]/norm));
          }
          fclose(FilePtr);
        }

        for(j=0;j<Framework[CurrentSystem].NumberOfTorsions[0];j++)
        {
          norm=0.0;
          for(k=0;k<DihedralHistogramSize[CurrentSystem];k++)
            norm+=FrameworkTorsionAngleHistogram[CurrentSystem][j][k];
          norm*=DihedralRange[CurrentSystem]/DihedralHistogramSize[CurrentSystem];

          sprintf(buffer,"MoleculeProperties/System_%d/Histogram_Framework_TorsionAngle_%d_%d_%d_%d_%d%s.dat",CurrentSystem,j,
             Framework[CurrentSystem].Torsions[0][j].A,
             Framework[CurrentSystem].Torsions[0][j].B,
             Framework[CurrentSystem].Torsions[0][j].C,
             Framework[CurrentSystem].Torsions[0][j].D,
             FileNameAppend);
          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# average torsion angle: %lf [degrees]\n",
            (double)(FrameworkAverageTorsionAngle[CurrentSystem][j]/FrameworkTorsionAngleCount[CurrentSystem][j]));
          for(k=0;k<DihedralHistogramSize[CurrentSystem];k++)
          {
            phi=(REAL)k/(REAL)DihedralHistogramSize[CurrentSystem];
            if(FrameworkTorsionAngleHistogram[CurrentSystem][j][k]>0.0)
              fprintf(FilePtr,"%lf %lf\n",(double)(phi*DihedralRange[CurrentSystem]),(double)(FrameworkTorsionAngleHistogram[CurrentSystem][j][k]/norm));
          }
          fclose(FilePtr);
        }

      }

      for(i=0;i<NumberOfComponents;i++)
      {
        NumberOfBonds=Components[i].NumberOfBonds;
        for(j=0;j<NumberOfBonds;j++)
        {
          norm=0.0;
          for(k=0;k<BondLengthHistogramSize[CurrentSystem];k++)
            norm+=BondLengthHistogram[CurrentSystem][i][j][k];
          norm*=BondLengthRange[CurrentSystem]/BondLengthHistogramSize[CurrentSystem];

          sprintf(buffer,"MoleculeProperties/System_%d/Histogram_%s_%d_BondLength%d%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
          FilePtr=fopen(buffer,"w");
          for(k=0;k<BondLengthHistogramSize[CurrentSystem];k++)
          {
            r=(REAL)k/(REAL)BondLengthHistogramSize[CurrentSystem];
            if(BondLengthHistogram[CurrentSystem][i][j][k]>0.0)
              fprintf(FilePtr,"%g %g\n",(double)(r*BondLengthRange[CurrentSystem]),(double)(BondLengthHistogram[CurrentSystem][i][j][k]/norm));
          }
          fclose(FilePtr);
        }

        NumberOfUreyBradleys=Components[i].NumberOfUreyBradleys;
        for(j=0;j<NumberOfUreyBradleys;j++)
        {
          norm=0.0;
          for(k=0;k<BondLengthHistogramSize[CurrentSystem];k++)
            norm+=UreyBradleyLengthHistogram[CurrentSystem][i][j][k];
          norm*=BondLengthRange[CurrentSystem]/BondLengthHistogramSize[CurrentSystem];

          sprintf(buffer,"MoleculeProperties/System_%d/Histogram_%s_%d_UreyBradley%d%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
          FilePtr=fopen(buffer,"w");
          for(k=0;k<BondLengthHistogramSize[CurrentSystem];k++)
          {
            r=(REAL)k/(REAL)BondLengthHistogramSize[CurrentSystem];
            if(UreyBradleyLengthHistogram[CurrentSystem][i][j][k]>0.0)
              fprintf(FilePtr,"%g %g\n",(double)(r*BondLengthRange[CurrentSystem]),(double)(UreyBradleyLengthHistogram[CurrentSystem][i][j][k]/norm));
          }
          fclose(FilePtr);
        }

        NumberOfBends=Components[i].NumberOfBends;
        for(j=0;j<NumberOfBends;j++)
        {
          norm=0.0;
          for(k=0;k<BendAngleHistogramSize[CurrentSystem];k++)
            norm+=BendAngleHistogram[CurrentSystem][i][j][k];
          norm*=BendAngleRange[CurrentSystem]/BendAngleHistogramSize[CurrentSystem];

          sprintf(buffer,"MoleculeProperties/System_%d/Histogram_%s_%d_BendAngle%d%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
          FilePtr=fopen(buffer,"w");
          for(k=0;k<BendAngleHistogramSize[CurrentSystem];k++)
          {
            theta=(REAL)k/(REAL)BendAngleHistogramSize[CurrentSystem];
            if(BendAngleHistogram[CurrentSystem][i][j][k]>0.0)
              fprintf(FilePtr,"%g %g\n",(double)(theta*BendAngleRange[CurrentSystem]),(double)(BendAngleHistogram[CurrentSystem][i][j][k]/norm));
          }
          fclose(FilePtr);
        }

        NumberOfTorsions=Components[i].NumberOfTorsions;
        for(j=0;j<NumberOfTorsions;j++)
        {
          sprintf(buffer,"MoleculeProperties/System_%d/Histogram_%s_%d_TorsionAngle%d%s.dat",CurrentSystem,Components[i].Name,i,j,FileNameAppend);
          FilePtr=fopen(buffer,"w");

          norm=TorsionConformationHistogram[CurrentSystem][i][j][SYNPERIPLANAR]+
               TorsionConformationHistogram[CurrentSystem][i][j][SYNCLINAL_PLUS]+
               TorsionConformationHistogram[CurrentSystem][i][j][ANTICLINAL_PLUS]+
               TorsionConformationHistogram[CurrentSystem][i][j][ANTIPERIPLANAR_PLUS]+
               TorsionConformationHistogram[CurrentSystem][i][j][ANTICLINAL_MIN]+
               TorsionConformationHistogram[CurrentSystem][i][j][SYNCLINAL_MIN];

          fprintf(FilePtr,"#Percentage conformations found in torsional forms:\n");
          fprintf(FilePtr,"#--------------------------------------------------\n");
          fprintf(FilePtr,"# synperiplanar form         [-30:30,330:360]: %lf %%\n",
                  (double)(TorsionConformationHistogram[CurrentSystem][i][j][SYNPERIPLANAR]*100.0/norm));
          fprintf(FilePtr,"# synclinal+ (gauche+) form           [30:90]: %lf %%\n",
                  (double)(TorsionConformationHistogram[CurrentSystem][i][j][SYNCLINAL_PLUS]*100.0/norm));
          fprintf(FilePtr,"# anticlinal+ form                   [90:150]: %lf %%\n",
                  (double)(TorsionConformationHistogram[CurrentSystem][i][j][ANTICLINAL_PLUS]*100.0/norm));
          fprintf(FilePtr,"# antiperiplanar+ (trans,anti) form [150:210]: %lf %%\n",
                  (double)(TorsionConformationHistogram[CurrentSystem][i][j][ANTIPERIPLANAR_PLUS]*100.0/norm));
          fprintf(FilePtr,"# anticlinal- form                  [210:270]: %lf %%\n",
                  (double)(TorsionConformationHistogram[CurrentSystem][i][j][ANTICLINAL_MIN]*100.0/norm));
          fprintf(FilePtr,"# synclinal- (gauche-) form         [270:330]: %lf %%\n",
                  (double)(TorsionConformationHistogram[CurrentSystem][i][j][SYNCLINAL_MIN]*100.0/norm));

          norm=0.0;
          for(k=0;k<DihedralHistogramSize[CurrentSystem];k++)
            norm+=TorsionAngleHistogram[CurrentSystem][i][j][k];
          norm*=DihedralRange[CurrentSystem]/DihedralHistogramSize[CurrentSystem];

          for(k=0;k<DihedralHistogramSize[CurrentSystem];k++)
          {
            theta=(REAL)k/(REAL)DihedralHistogramSize[CurrentSystem];
            if(TorsionAngleHistogram[CurrentSystem][i][j][k]>0.0)
              fprintf(FilePtr,"%g %g\n",(double)(theta*DihedralRange[CurrentSystem]),(double)(TorsionAngleHistogram[CurrentSystem][i][j][k]/norm));
          }
          fclose(FilePtr);
        }
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeMoleculeProperties[i])
        {
          for(j=0;j<NumberOfComponents;j++)
          {
            for(k=0;k<MaxNumberOfBonds;k++)
              free(BondLengthHistogram[i][j][k]);
            free(BondLengthHistogram[i][j]);

            for(k=0;k<MaxNumberOfUreyBradleys;k++)
              free(UreyBradleyLengthHistogram[i][j][k]);
            free(UreyBradleyLengthHistogram[i][j]);

            for(k=0;k<MaxNumberOfBends;k++)
              free(BendAngleHistogram[i][j][k]);
            free(BendAngleHistogram[i][j]);

            for(k=0;k<MaxNumberOfTorsions;k++)
              free(TorsionAngleHistogram[i][j][k]);
            free(TorsionAngleHistogram[i][j]);

            for(k=0;k<MaxNumberOfTorsions;k++)
              free(TorsionConformationHistogram[i][j][k]);
            free(TorsionConformationHistogram[i][j]);
          }
          free(BondLengthHistogram[i]);
          free(UreyBradleyLengthHistogram[i]);
          free(BendAngleHistogram[i]);
          free(TorsionAngleHistogram[i]);
          free(TorsionConformationHistogram[i]);
        }
      }

      free(BondLengthHistogram);
      free(UreyBradleyLengthHistogram);
      free(BendAngleHistogram);
      free(TorsionAngleHistogram);
      free(TorsionConformationHistogram);


      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeMoleculeProperties[i]&&(Framework[i].FrameworkModel==FLEXIBLE))
        {

          for(j=0;j<Framework[i].NumberOfBonds[0];j++)
            free(FrameworkBondLengthHistogram[i][j]);
          free(FrameworkBondLengthHistogram[i]);

          for(j=0;j<Framework[i].NumberOfUreyBradleys[0];j++)
            free(FrameworkUreyBradleyLengthHistogram[i][j]);
          free(FrameworkUreyBradleyLengthHistogram[i]);

          for(j=0;j<Framework[i].NumberOfBends[0];j++)
            free(FrameworkBendAngleHistogram[i][j]);
          free(FrameworkBendAngleHistogram[i]);

          for(j=0;j<Framework[i].NumberOfTorsions[0];j++)
            free(FrameworkTorsionAngleHistogram[i][j]);
          free(FrameworkTorsionAngleHistogram[i]);

          free(FrameworkAverageBondLength[i]);
          free(FrameworkBondLengthCount[i]);
          free(FrameworkAverageBendAngle[i]);
          free(FrameworkBendAngleCount[i]);
          free(FrameworkAverageTorsionAngle[i]);
          free(FrameworkTorsionAngleCount[i]);
        }
      }
      free(FrameworkBondLengthHistogram);
      free(FrameworkUreyBradleyLengthHistogram);
      free(FrameworkBendAngleHistogram);
      free(FrameworkTorsionAngleHistogram);
      free(FrameworkAverageBondLength);
      free(FrameworkBondLengthCount);
      free(FrameworkAverageBendAngle);
      free(FrameworkBendAngleCount);
      free(FrameworkAverageTorsionAngle);
      free(FrameworkTorsionAngleCount);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleInfraRedSpectra                                                                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the IR spectra (spacings: 2048, 4196, 8192, 16384, 32768 points).                *
 * Parameters | -                                                                                        *
 * Note       | zero-padding and/or windowing gives poor results                                         *
 *********************************************************************************************************/


#define TRIANGULAR(j,m) (1.0-fabs(((j))-(m))/(m))            /* Bartlett window */
#define SQUARE(j,m) 1.0                                        /* Square window */
#define WELTCH(j,m) (1.0-SQR(((j)-(m))/(m)))             /* Weltch window */
#define HANNING(j,m) (0.5*(1.0+cos(((j)-(m))*M_PI/(m))))   /* Hanning window */
#define HAMMING(j,m) (0.54+0.46*cos((((j)-(m)))*M_PI/(m))) /* Hamming window */

#define FFT_WINDOW(j,m) SQUARE((double)j,(double)m)    /* use Weltch window */

void SampleInfraRedSpectra(int Switch)
{
  int i,j,k,l,m,n;
  int index,Type;
  REAL Charge,w;
  REAL corrvacf,corrfvacf,corravacf,corrcvacf;
  REAL corr,corrf,corra,corrc;
  VECTOR vel,drift_system;
  FILE *FilePtr;
  char buffer[256];
  static REAL (*SpectrumBuffer)[2];
  REAL conversion_factor;
  VECTOR corrvacf_all,corrvacf_framework,corrvacf_adsorbates,corrvacf_cations;
  VECTOR corr_all,corr_framework,corr_adsorbates,corr_cations;
  static VECTOR *corr_pseudoatoms=NULL,*corrvacf_pseudoatoms=NULL;
  static REAL *corr_pseudoatoms_isotropic=NULL,*corrvacf_pseudoatoms_isotropic=NULL;
  int fft_size,half_fft_size;

  switch(Switch)
  {
    case ALLOCATE:
      // allocate all the needed memory
      corr_pseudoatoms=(VECTOR*)calloc(NumberOfPseudoAtoms,sizeof(VECTOR));
      corrvacf_pseudoatoms=(VECTOR*)calloc(NumberOfPseudoAtoms,sizeof(VECTOR));
      corr_pseudoatoms_isotropic=(REAL*)calloc(NumberOfPseudoAtoms,sizeof(REAL));
      corrvacf_pseudoatoms_isotropic=(REAL*)calloc(NumberOfPseudoAtoms,sizeof(REAL));
      Spectrum=(REAL ****)calloc(NumberOfSystems,sizeof(REAL***));
      SpectrumAverage=(REAL ****)calloc(NumberOfSystems,sizeof(REAL***));
      UnweightedSpectrum=(REAL ****)calloc(NumberOfSystems,sizeof(REAL***));
      UnweightedSpectrumAverage=(REAL ****)calloc(NumberOfSystems,sizeof(REAL***));
      SpectrumPseudoAtoms=(REAL *****)calloc(NumberOfSystems,sizeof(REAL****));
      SpectrumPseudoAtomsAverage=(REAL *****)calloc(NumberOfSystems,sizeof(REAL****));
      SpectrumCount=(REAL **)calloc(NumberOfSystems,sizeof(REAL*));
      sumw=(REAL*)calloc(5,sizeof(REAL));
      SpectrumBuffer=(REAL(*)[2])calloc(5*32768,sizeof(REAL[2]));
      for(k=0;k<NumberOfSystems;k++)
      {
        if(ComputeInfraRedSpectra[k])
        {
          Spectrum[k]=(REAL ***)calloc(4,sizeof(REAL**));
          SpectrumAverage[k]=(REAL ***)calloc(4,sizeof(REAL**));
          UnweightedSpectrum[k]=(REAL ***)calloc(4,sizeof(REAL**));
          UnweightedSpectrumAverage[k]=(REAL ***)calloc(4,sizeof(REAL**));
          SpectrumPseudoAtoms[k]=(REAL ****)calloc(2,sizeof(REAL***));
          SpectrumPseudoAtomsAverage[k]=(REAL ****)calloc(2,sizeof(REAL***));
          SpectrumCount[k]=(REAL *)calloc(5,sizeof(REAL));
          for(i=0;i<4;i++)
          {
            Spectrum[k][i]=(REAL **)malloc(5*sizeof(REAL*));
            Spectrum[k][i][0]=(REAL *)calloc(4*2048,sizeof(REAL));
            Spectrum[k][i][1]=(REAL *)calloc(4*4096,sizeof(REAL));
            Spectrum[k][i][2]=(REAL *)calloc(4*8192,sizeof(REAL));
            Spectrum[k][i][3]=(REAL *)calloc(4*16384,sizeof(REAL));
            Spectrum[k][i][4]=(REAL *)calloc(4*32768,sizeof(REAL));

            SpectrumAverage[k][i]=(REAL **)malloc(5*sizeof(REAL*));
            SpectrumAverage[k][i][0]=(REAL *)calloc(4*2048+1,sizeof(REAL));
            SpectrumAverage[k][i][1]=(REAL *)calloc(4*4096+1,sizeof(REAL));
            SpectrumAverage[k][i][2]=(REAL *)calloc(4*8192+1,sizeof(REAL));
            SpectrumAverage[k][i][3]=(REAL *)calloc(4*16384+1,sizeof(REAL));
            SpectrumAverage[k][i][4]=(REAL *)calloc(4*32768+1,sizeof(REAL));

            UnweightedSpectrum[k][i]=(REAL **)malloc(5*sizeof(REAL*));
            UnweightedSpectrum[k][i][0]=(REAL *)calloc(4*2048+1,sizeof(REAL));
            UnweightedSpectrum[k][i][1]=(REAL *)calloc(4*4096+1,sizeof(REAL));
            UnweightedSpectrum[k][i][2]=(REAL *)calloc(4*8192+1,sizeof(REAL));
            UnweightedSpectrum[k][i][3]=(REAL *)calloc(4*16384+1,sizeof(REAL));
            UnweightedSpectrum[k][i][4]=(REAL *)calloc(4*32768+1,sizeof(REAL));

            UnweightedSpectrumAverage[k][i]=(REAL **)malloc(5*sizeof(REAL*));
            UnweightedSpectrumAverage[k][i][0]=(REAL *)calloc(4*2048+1,sizeof(REAL));
            UnweightedSpectrumAverage[k][i][1]=(REAL *)calloc(4*4096+1,sizeof(REAL));
            UnweightedSpectrumAverage[k][i][2]=(REAL *)calloc(4*8192+1,sizeof(REAL));
            UnweightedSpectrumAverage[k][i][3]=(REAL *)calloc(4*16384+1,sizeof(REAL));
            UnweightedSpectrumAverage[k][i][4]=(REAL *)calloc(4*32768+1,sizeof(REAL));
          }
          for(i=0;i<2;i++)
          {
            SpectrumPseudoAtoms[k][i]=(REAL ***)malloc(NumberOfPseudoAtoms*sizeof(REAL**));
            SpectrumPseudoAtomsAverage[k][i]=(REAL ***)malloc(NumberOfPseudoAtoms*sizeof(REAL**));
            for(j=0;j<NumberOfPseudoAtoms;j++)
            {
               if(NumberOfPseudoAtomsType[k][j]>0)
               {
                 SpectrumPseudoAtoms[k][i][j]=(REAL **)malloc(5*sizeof(REAL*));
                 SpectrumPseudoAtoms[k][i][j][0]=(REAL *)calloc(4*2048+1,sizeof(REAL));
                 SpectrumPseudoAtoms[k][i][j][1]=(REAL *)calloc(4*4096+1,sizeof(REAL));
                 SpectrumPseudoAtoms[k][i][j][2]=(REAL *)calloc(4*8192+1,sizeof(REAL));
                 SpectrumPseudoAtoms[k][i][j][3]=(REAL *)calloc(4*16384+1,sizeof(REAL));
                 SpectrumPseudoAtoms[k][i][j][4]=(REAL *)calloc(4*32768+1,sizeof(REAL));

                 SpectrumPseudoAtomsAverage[k][i][j]=(REAL **)malloc(5*sizeof(REAL*));
                 SpectrumPseudoAtomsAverage[k][i][j][0]=(REAL *)calloc(4*2048+1,sizeof(REAL));
                 SpectrumPseudoAtomsAverage[k][i][j][1]=(REAL *)calloc(4*4096+1,sizeof(REAL));
                 SpectrumPseudoAtomsAverage[k][i][j][2]=(REAL *)calloc(4*8192+1,sizeof(REAL));
                 SpectrumPseudoAtomsAverage[k][i][j][3]=(REAL *)calloc(4*16384+1,sizeof(REAL));
                 SpectrumPseudoAtomsAverage[k][i][j][4]=(REAL *)calloc(4*32768+1,sizeof(REAL));
               }
             }
          }
        }
      }
      break;
    case INITIALIZE:
      for(j=0;j<NumberOfSystems;j++)
      {
        if(ComputeInfraRedSpectra[j])
        {

          m=2048;
          for(k=0;k<5;k++)
          {
            fft_size=4*m;
            half_fft_size=2*m;
            SpectrumCount[j][k]=0.0;
            sumw[k]=0.0;
            for(i=0;i<fft_size;i++)
              sumw[k]+=SQR(FFT_WINDOW(i,half_fft_size));
            m*=2;
          }
        }
      }
      break;
    case SAMPLE:
      if((!ComputeInfraRedSpectra[CurrentSystem])||(CurrentCycle%SampleEveryInfraRed!=0)) return;

      drift_system=GetCenterOfMassVelocityCurrentSystem();

      corrvacf_all.x=0.0;
      corrvacf_all.y=0.0;
      corrvacf_all.z=0.0;
      corr_all.x=0.0;
      corr_all.y=0.0;
      corr_all.z=0.0;

      corrvacf_framework.x=0.0;
      corrvacf_framework.y=0.0;
      corrvacf_framework.z=0.0;
      corr_framework.x=0.0;
      corr_framework.y=0.0;
      corr_framework.z=0.0;

      for(n=0;n<NumberOfPseudoAtoms;n++)
      {
        corrvacf_pseudoatoms[n].x=corrvacf_pseudoatoms[n].y=corrvacf_pseudoatoms[n].z=0.0;
        corr_pseudoatoms[n].x=corr_pseudoatoms[n].y=corr_pseudoatoms[n].z=0.0;
        corr_pseudoatoms_isotropic[n]=corrvacf_pseudoatoms_isotropic[n]=0.0;
      }

      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
        {
          Type=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
          vel=Framework[CurrentSystem].Atoms[CurrentFramework][i].Velocity;
          vel.x-=drift_system.x;
          vel.y-=drift_system.y;
          vel.z-=drift_system.z;
          Charge=Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge;

          corrvacf_all.x+=vel.x;
          corrvacf_all.y+=vel.y;
          corrvacf_all.z+=vel.z;
          corrvacf_framework.x+=vel.x;
          corrvacf_framework.y+=vel.y;
          corrvacf_framework.z+=vel.z;
          corrvacf_pseudoatoms[Type].x+=vel.x;
          corrvacf_pseudoatoms[Type].y+=vel.y;
          corrvacf_pseudoatoms[Type].z+=vel.z;

          corr_framework.x+=Charge*vel.x;
          corr_framework.y+=Charge*vel.y;
          corr_framework.z+=Charge*vel.z;
          corr_all.x+=Charge*vel.x;
          corr_all.y+=Charge*vel.y;
          corr_all.z+=Charge*vel.z;
          corr_pseudoatoms[Type].x+=Charge*vel.x;
          corr_pseudoatoms[Type].y+=Charge*vel.y;
          corr_pseudoatoms[Type].z+=Charge*vel.z;
        }
      }

      corrvacf_adsorbates.x=0.0;
      corrvacf_adsorbates.y=0.0;
      corrvacf_adsorbates.z=0.0;
      corr_adsorbates.x=0.0;
      corr_adsorbates.y=0.0;
      corr_adsorbates.z=0.0;
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          vel=Adsorbates[CurrentSystem][i].Atoms[j].Velocity;
          vel.x-=drift_system.x;
          vel.y-=drift_system.y;
          vel.z-=drift_system.z;
          Charge=Adsorbates[CurrentSystem][i].Atoms[j].Charge;

          corrvacf_all.x+=vel.x;
          corrvacf_all.y+=vel.y;
          corrvacf_all.z+=vel.z;
          corrvacf_adsorbates.x+=vel.x;
          corrvacf_adsorbates.y+=vel.y;
          corrvacf_adsorbates.z+=vel.z;
          corrvacf_pseudoatoms[Type].x+=vel.x;
          corrvacf_pseudoatoms[Type].y+=vel.y;
          corrvacf_pseudoatoms[Type].z+=vel.z;

          corr_adsorbates.x+=Charge*vel.x;
          corr_adsorbates.y+=Charge*vel.y;
          corr_adsorbates.z+=Charge*vel.z;
          corr_all.x+=Charge*vel.x;
          corr_all.y+=Charge*vel.y;
          corr_all.z+=Charge*vel.z;
          corr_pseudoatoms[Type].x+=Charge*vel.x;
          corr_pseudoatoms[Type].y+=Charge*vel.y;
          corr_pseudoatoms[Type].z+=Charge*vel.z;
        }
      }

      corrvacf_cations.x=0.0;
      corrvacf_cations.y=0.0;
      corrvacf_cations.z=0.0;
      corr_cations.x=0.0;
      corr_cations.y=0.0;
      corr_cations.z=0.0;
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          Type=Cations[CurrentSystem][i].Atoms[j].Type;
          vel=Cations[CurrentSystem][i].Atoms[j].Velocity;
          vel.x-=drift_system.x;
          vel.y-=drift_system.y;
          vel.z-=drift_system.z;
          Charge=Cations[CurrentSystem][i].Atoms[j].Charge;

          corrvacf_all.x+=vel.x;
          corrvacf_all.y+=vel.y;
          corrvacf_all.z+=vel.z;
          corrvacf_cations.x+=vel.x;
          corrvacf_cations.y+=vel.y;
          corrvacf_cations.z+=vel.z;
          corrvacf_pseudoatoms[Type].x+=vel.x;
          corrvacf_pseudoatoms[Type].y+=vel.y;
          corrvacf_pseudoatoms[Type].z+=vel.z;

          corr_cations.x+=Charge*vel.x;
          corr_cations.y+=Charge*vel.y;
          corr_cations.z+=Charge*vel.z;
          corr_all.x+=Charge*vel.x;
          corr_all.y+=Charge*vel.y;
          corr_all.z+=Charge*vel.z;
          corr_pseudoatoms[Type].x+=Charge*vel.x;
          corr_pseudoatoms[Type].y+=Charge*vel.y;
          corr_pseudoatoms[Type].z+=Charge*vel.z;
        }
      }
      corrvacf=(corrvacf_all.x+corrvacf_all.y+corrvacf_all.z)/3.0;
      corrfvacf=(corrvacf_framework.x+corrvacf_framework.y+corrvacf_framework.z)/3.0;
      corravacf=(corrvacf_adsorbates.x+corrvacf_adsorbates.y+corrvacf_adsorbates.z)/3.0;
      corrcvacf=(corrvacf_cations.x+corrvacf_cations.y+corrvacf_cations.z)/3.0;

      corr=(corr_all.x+corr_all.y+corr_all.z)/3.0;
      corrf=(corr_framework.x+corr_framework.y+corr_framework.z)/3.0;
      corra=(corr_adsorbates.x+corr_adsorbates.y+corr_adsorbates.z)/3.0;
      corrc=(corr_cations.x+corr_cations.y+corr_cations.z)/3.0;

      for(n=0;n<NumberOfPseudoAtoms;n++)
      {
        if(NumberOfPseudoAtomsType[CurrentSystem][n]>0)
        {
          corrvacf_pseudoatoms_isotropic[n]=(corrvacf_pseudoatoms[n].x+corrvacf_pseudoatoms[n].y+corrvacf_pseudoatoms[n].z)/3.0;
          corr_pseudoatoms_isotropic[n]=(corr_pseudoatoms[n].x+corr_pseudoatoms[n].y+corr_pseudoatoms[n].z)/3.0;
        }
      }

      m=2048;
      for(k=0;k<5;k++)
      {
        if(CurrentCycle<4*m)
        {
          UnweightedSpectrum[CurrentSystem][0][k][CurrentCycle]=corrvacf;
          UnweightedSpectrum[CurrentSystem][1][k][CurrentCycle]=corrfvacf;
          UnweightedSpectrum[CurrentSystem][2][k][CurrentCycle]=corravacf;
          UnweightedSpectrum[CurrentSystem][3][k][CurrentCycle]=corrcvacf;

          Spectrum[CurrentSystem][0][k][CurrentCycle]=corr;
          Spectrum[CurrentSystem][1][k][CurrentCycle]=corrf;
          Spectrum[CurrentSystem][2][k][CurrentCycle]=corra;
          Spectrum[CurrentSystem][3][k][CurrentCycle]=corrc;

          for(n=0;n<NumberOfPseudoAtoms;n++)
          {
            if(NumberOfPseudoAtomsType[CurrentSystem][n]>0)
            {
              SpectrumPseudoAtoms[CurrentSystem][0][n][k][CurrentCycle]=corr_pseudoatoms_isotropic[n];
              SpectrumPseudoAtoms[CurrentSystem][1][n][k][CurrentCycle]=corrvacf_pseudoatoms_isotropic[n];
            }
          }
        }
        else
        {
          if((index=(CurrentCycle%(2*m)))!=0)
          {
            UnweightedSpectrum[CurrentSystem][0][k][index+2*m]=corrvacf;
            UnweightedSpectrum[CurrentSystem][1][k][index+2*m]=corrfvacf;
            UnweightedSpectrum[CurrentSystem][2][k][index+2*m]=corravacf;
            UnweightedSpectrum[CurrentSystem][3][k][index+2*m]=corrcvacf;

            Spectrum[CurrentSystem][0][k][index+2*m]=corr;
            Spectrum[CurrentSystem][1][k][index+2*m]=corrf;
            Spectrum[CurrentSystem][2][k][index+2*m]=corra;
            Spectrum[CurrentSystem][3][k][index+2*m]=corrc;

            for(n=0;n<NumberOfPseudoAtoms;n++)
            {
              if(NumberOfPseudoAtomsType[CurrentSystem][n]>0)
              {
                SpectrumPseudoAtoms[CurrentSystem][0][n][k][index+2*m]=corr_pseudoatoms_isotropic[n];
                SpectrumPseudoAtoms[CurrentSystem][1][n][k][index+2*m]=corrvacf_pseudoatoms_isotropic[n];
              }
            }
          }
          else // periodogram full, do FFT and shift down
          {
            SpectrumCount[CurrentSystem][k]+=1.0;
            for(l=0;l<4;l++)
            {
              //================================================================
              // FFT of the velocity autocorrelation function (unweighted)
              //================================================================
              fft_size=4*m;
              half_fft_size=2*m;

              // fill even (real) part
              for (i=0;i<fft_size;i++)
              {
                SpectrumBuffer[i][0]=UnweightedSpectrum[CurrentSystem][l][k][i];
                SpectrumBuffer[i][1]=0.0;
              }

              // "window" the data
              for(i=0;i<fft_size;i++)
              {
                w=FFT_WINDOW(i,half_fft_size);
                SpectrumBuffer[i][0]*=w;
              }

              // FFT the frame
              FastFourierTransform(SpectrumBuffer,fft_size,1);
              for(i=0;i<2*(fft_size/2+1);i++)
                UnweightedSpectrumAverage[CurrentSystem][l][k][i]+=SQR(SpectrumBuffer[i][0])+SQR(SpectrumBuffer[i][1]);

              // shift data down by half_fft_size
              for(i=0;i<2*m;i++)
                UnweightedSpectrum[CurrentSystem][l][k][i]=UnweightedSpectrum[CurrentSystem][l][k][i+2*m];

              //================================================================
              // charge weighted FFT of the velocity autocorrelation function
              //================================================================

              // fill even (real) part
              for (i=0;i<fft_size;i++)
              {
                SpectrumBuffer[i][0]=Spectrum[CurrentSystem][l][k][i];
                SpectrumBuffer[i][1]=0.0;
              }

              // "window" the data
              for(i=0;i<fft_size;i++)
              {
                w=FFT_WINDOW(i,half_fft_size);
                SpectrumBuffer[i][0]*=w;
              }

              // FFT the frame
              FastFourierTransform(SpectrumBuffer,fft_size,1);
              for(i=0;i<2*(fft_size/2+1);i++)
                SpectrumAverage[CurrentSystem][l][k][i]+=SQR(SpectrumBuffer[i][0])+SQR(SpectrumBuffer[i][1]);

              for(i=0;i<2*m;i++)
                Spectrum[CurrentSystem][l][k][i]=Spectrum[CurrentSystem][l][k][i+2*m];
            }

            //================================================================
            // FFT of the velocity autocorrelation function (unweighted)
            //================================================================

            for(n=0;n<NumberOfPseudoAtoms;n++)
            {
              if(NumberOfPseudoAtomsType[CurrentSystem][n]>0)
              {
                // fill even (real) part
                for (i=0;i<fft_size;i++)
                {
                  SpectrumBuffer[i][0]=SpectrumPseudoAtoms[CurrentSystem][0][n][k][i];
                  SpectrumBuffer[i][1]=0.0;
                }

                // "window" the data
                for(i=0;i<fft_size;i++)
                {
                  w=FFT_WINDOW(i,half_fft_size);
                  SpectrumBuffer[i][0]*=w;
                }

                // FFT the frame
                FastFourierTransform(SpectrumBuffer,fft_size,1);
                for(i=0;i<2*(fft_size/2+1);i++)
                  SpectrumPseudoAtomsAverage[CurrentSystem][0][n][k][i]+=SQR(SpectrumBuffer[i][0])+SQR(SpectrumBuffer[i][1]);

                for(i=0;i<2*m;i++)
                  SpectrumPseudoAtoms[CurrentSystem][0][n][k][i]=SpectrumPseudoAtoms[CurrentSystem][0][n][k][i+2*m];
              }
            }

            //================================================================
            // charge weighted FFT of the velocity autocorrelation function
            //================================================================

            for(n=0;n<NumberOfPseudoAtoms;n++)
            {
              if(NumberOfPseudoAtomsType[CurrentSystem][n]>0)
              {
                // fill even (real) part
                for (i=0;i<fft_size;i++)
                {
                  SpectrumBuffer[i][0]=SpectrumPseudoAtoms[CurrentSystem][1][n][k][i];
                  SpectrumBuffer[i][1]=0.0;
                }

                // "window" the data
                for(i=0;i<fft_size;i++)
                {
                  w=FFT_WINDOW(i,half_fft_size);
                  SpectrumBuffer[i][0]*=w;
                }

                // FFT the frame
                FastFourierTransform(SpectrumBuffer,fft_size,1);
                for(i=0;i<2*(fft_size/2+1);i++)
                  SpectrumPseudoAtomsAverage[CurrentSystem][1][n][k][i]+=SQR(SpectrumBuffer[i][0])+SQR(SpectrumBuffer[i][1]);

                // shift down by
                for(i=0;i<2*m;i++)
                  SpectrumPseudoAtoms[CurrentSystem][1][n][k][i]=SpectrumPseudoAtoms[CurrentSystem][1][n][k][i+2*m];
              }
            }
          }
        }
        m*=2;
      }
      break;
    case PRINT:
      if((!ComputeInfraRedSpectra[CurrentSystem])||(CurrentCycle%WriteInfraRedSpectraEvery[CurrentSystem]!=0)) return;

      mkdir("MD_Spectra",S_IRWXU);
      sprintf(buffer,"MD_Spectra/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);
      sprintf(buffer,"MD_Spectra/System_%d/ChargeWeightedVACF_FFT",CurrentSystem);
      mkdir(buffer,S_IRWXU);
      sprintf(buffer,"MD_Spectra/System_%d/VACF_FFT",CurrentSystem);
      mkdir(buffer,S_IRWXU);
      sprintf(buffer,"MD_Spectra/System_%d/ChargeWeightedVACF_FFT_PseudoAtoms",CurrentSystem);
      mkdir(buffer,S_IRWXU);
      sprintf(buffer,"MD_Spectra/System_%d/VACF_FFT_PseudoAtoms",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      m=2048;
      for(k=0;k<5;k++)
      {
        conversion_factor=TO_WAVENUMBERS*M_PI/(2.0*m*DeltaT);

        fft_size=4*m;
        half_fft_size=2*m;

        //================================================================
        // FFT of the velocity autocorrelation function (unweighted)
        //================================================================

        sprintf(buffer,"MD_Spectra/System_%d/VACF_FFT/Spectrum_%d_all%s.dat",CurrentSystem,m,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# %lld samples (half overlapping)\n",(long long)SpectrumCount[CurrentSystem][k]);
        fprintf(FilePtr,"# column 1 frequencies in cm^-1\n");
        fprintf(FilePtr,"# column 2 intensity\n");
        for(i=0;i<half_fft_size;i++)
          if(SpectrumCount[CurrentSystem][k]>0.0)
            fprintf(FilePtr,"%g %g\n",i*conversion_factor,UnweightedSpectrumAverage[CurrentSystem][0][k][i]/(m*SpectrumCount[CurrentSystem][k]));
        fclose(FilePtr);

        sprintf(buffer,"MD_Spectra/System_%d/VACF_FFT/Spectrum_%d_framework%s.dat",CurrentSystem,m,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# %lld samples (half overlapping)\n",(long long)SpectrumCount[CurrentSystem][k]);
        fprintf(FilePtr,"# column 1 frequencies in cm^-1\n");
        fprintf(FilePtr,"# column 2 intensity\n");
        for(i=0;i<half_fft_size;i++)
          if(SpectrumCount[CurrentSystem][k]>0.0)
            fprintf(FilePtr,"%g %g\n",i*conversion_factor,UnweightedSpectrumAverage[CurrentSystem][1][k][i]/(m*SpectrumCount[CurrentSystem][k]));
        fclose(FilePtr);

        sprintf(buffer,"MD_Spectra/System_%d/VACF_FFT/Spectrum_%d_adsorbates%s.dat",CurrentSystem,m,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# %lld samples (half overlapping)\n",(long long)SpectrumCount[CurrentSystem][k]);
        fprintf(FilePtr,"# column 1 frequencies in cm^-1\n");
        fprintf(FilePtr,"# column 2 intensity\n");
        for(i=0;i<half_fft_size;i++)
          if(SpectrumCount[CurrentSystem][k]>0.0)
            fprintf(FilePtr,"%g %g\n",i*conversion_factor,UnweightedSpectrumAverage[CurrentSystem][2][k][i]/(m*SpectrumCount[CurrentSystem][k]));
        fclose(FilePtr);

        sprintf(buffer,"MD_Spectra/System_%d/VACF_FFT/Spectrum_%d_cations%s.dat",CurrentSystem,m,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# %lld samples (half overlapping)\n",(long long)SpectrumCount[CurrentSystem][k]);
        fprintf(FilePtr,"# column 1 frequencies in cm^-1\n");
        fprintf(FilePtr,"# column 2 intensity\n");
        for(i=0;i<half_fft_size;i++)
          if(SpectrumCount[CurrentSystem][k]>0.0)
            fprintf(FilePtr,"%g %g\n",i*conversion_factor,UnweightedSpectrumAverage[CurrentSystem][3][k][i]/(m*SpectrumCount[CurrentSystem][k]));
        fclose(FilePtr);


        //================================================================
        // charge weighted FFT of the velocity autocorrelation function
        //================================================================

        sprintf(buffer,"MD_Spectra/System_%d/ChargeWeightedVACF_FFT/Spectrum_%d_all%s.dat",CurrentSystem,m,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# %lld samples (half overlapping)\n",(long long)SpectrumCount[CurrentSystem][k]);
        fprintf(FilePtr,"# column 1 frequencies in cm^-1\n");
        fprintf(FilePtr,"# column 2 intensity\n");
        for(i=0;i<half_fft_size;i++)
          if(SpectrumCount[CurrentSystem][k]>0.0)
            fprintf(FilePtr,"%g %g\n",i*conversion_factor,SpectrumAverage[CurrentSystem][0][k][i]/(m*SpectrumCount[CurrentSystem][k]));
        fclose(FilePtr);

        sprintf(buffer,"MD_Spectra/System_%d/ChargeWeightedVACF_FFT/Spectrum_%d_framework%s.dat",CurrentSystem,m,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# %lld samples (half overlapping)\n",(long long)SpectrumCount[CurrentSystem][k]);
        fprintf(FilePtr,"# column 1 frequencies in cm^-1\n");
        fprintf(FilePtr,"# column 2 intensity\n");
        for(i=0;i<half_fft_size;i++)
          if(SpectrumCount[CurrentSystem][k]>0.0)
            fprintf(FilePtr,"%g %g\n",i*conversion_factor,SpectrumAverage[CurrentSystem][1][k][i]/(m*SpectrumCount[CurrentSystem][k]));
        fclose(FilePtr);

        sprintf(buffer,"MD_Spectra/System_%d/ChargeWeightedVACF_FFT/Spectrum_%d_adsorbates%s.dat",CurrentSystem,m,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# %lld samples (half overlapping)\n",(long long)SpectrumCount[CurrentSystem][k]);
        fprintf(FilePtr,"# column 1 frequencies in cm^-1\n");
        fprintf(FilePtr,"# column 2 intensity\n");
        for(i=0;i<half_fft_size;i++)
          if(SpectrumCount[CurrentSystem][k]>0.0)
            fprintf(FilePtr,"%g %g\n",i*conversion_factor,SpectrumAverage[CurrentSystem][2][k][i]/(m*SpectrumCount[CurrentSystem][k]));
        fclose(FilePtr);

        sprintf(buffer,"MD_Spectra/System_%d/ChargeWeightedVACF_FFT/Spectrum_%d_cations%s.dat",CurrentSystem,m,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# %lld samples (half overlapping)\n",(long long)SpectrumCount[CurrentSystem][k]);
        fprintf(FilePtr,"# column 1 frequencies in cm^-1\n");
        fprintf(FilePtr,"# column 2 intensity\n");
        for(i=0;i<half_fft_size;i++)
          if(SpectrumCount[CurrentSystem][k]>0.0)
            fprintf(FilePtr,"%g %g\n",i*conversion_factor,SpectrumAverage[CurrentSystem][3][k][i]/(m*SpectrumCount[CurrentSystem][k]));
        fclose(FilePtr);

        for(n=0;n<NumberOfPseudoAtoms;n++)
        {
          if(NumberOfPseudoAtomsType[CurrentSystem][n]>0)
          {
            sprintf(buffer,"MD_Spectra/System_%d/ChargeWeightedVACF_FFT_PseudoAtoms/Spectrum_%s_%d%s.dat",CurrentSystem,PseudoAtoms[n].Name,m,FileNameAppend);
            FilePtr=fopen(buffer,"w");
            fprintf(FilePtr,"# %lld samples (half overlapping)\n",(long long)SpectrumCount[CurrentSystem][k]);
            fprintf(FilePtr,"# column 1 frequencies in cm^-1\n");
            fprintf(FilePtr,"# column 2 intensity\n");
            for(i=0;i<half_fft_size;i++)
              if(SpectrumCount[CurrentSystem][k]>0.0)
                fprintf(FilePtr,"%g %g\n",i*conversion_factor,SpectrumPseudoAtomsAverage[CurrentSystem][0][n][k][i]/(m*SpectrumCount[CurrentSystem][k]));
            fclose(FilePtr);
          }
        }

        for(n=0;n<NumberOfPseudoAtoms;n++)
        {
          if(NumberOfPseudoAtomsType[CurrentSystem][n]>0)
          {
            sprintf(buffer,"MD_Spectra/System_%d/VACF_FFT_PseudoAtoms/Spectrum_%s_%d%s.dat",CurrentSystem,PseudoAtoms[n].Name,m,FileNameAppend);
            FilePtr=fopen(buffer,"w");
            fprintf(FilePtr,"# %lld samples (half overlapping)\n",(long long)SpectrumCount[CurrentSystem][k]);
            fprintf(FilePtr,"# column 1 frequencies in cm^-1\n");
            fprintf(FilePtr,"# column 2 intensity\n");
            for(i=0;i<half_fft_size;i++)
              if(SpectrumCount[CurrentSystem][k]>0.0)
                fprintf(FilePtr,"%g %g\n",i*conversion_factor,SpectrumPseudoAtomsAverage[CurrentSystem][1][n][k][i]/(m*SpectrumCount[CurrentSystem][k]));
            fclose(FilePtr);
          }
        }
        m*=2;
      }
      break;
    case FINALIZE:
      // allocate all the needed memory
      for(k=0;k<NumberOfSystems;k++)
      {
        if(ComputeInfraRedSpectra[k])
        {
          for(i=0;i<4;i++)
          {
            free(Spectrum[k][i][0]);
            free(Spectrum[k][i][1]);
            free(Spectrum[k][i][2]);
            free(Spectrum[k][i][3]);
            free(Spectrum[k][i][4]);

            free(SpectrumAverage[k][i][0]);
            free(SpectrumAverage[k][i][1]);
            free(SpectrumAverage[k][i][2]);
            free(SpectrumAverage[k][i][3]);
            free(SpectrumAverage[k][i][4]);

            free(UnweightedSpectrum[k][i][0]);
            free(UnweightedSpectrum[k][i][1]);
            free(UnweightedSpectrum[k][i][2]);
            free(UnweightedSpectrum[k][i][3]);
            free(UnweightedSpectrum[k][i][4]);

            free(UnweightedSpectrumAverage[k][i][0]);
            free(UnweightedSpectrumAverage[k][i][1]);
            free(UnweightedSpectrumAverage[k][i][2]);
            free(UnweightedSpectrumAverage[k][i][3]);
            free(UnweightedSpectrumAverage[k][i][4]);

            free(Spectrum[k][i]);
            free(SpectrumAverage[k][i]);
            free(UnweightedSpectrum[k][i]);
            free(UnweightedSpectrumAverage[k][i]);
          }
          for(i=0;i<2;i++)
          {
            for(j=0;j<NumberOfPseudoAtoms;j++)
            {
              if(NumberOfPseudoAtomsType[k][j]>0)
              {
                free(SpectrumPseudoAtoms[k][i][j][0]);
                free(SpectrumPseudoAtoms[k][i][j][1]);
                free(SpectrumPseudoAtoms[k][i][j][2]);
                free(SpectrumPseudoAtoms[k][i][j][3]);
                free(SpectrumPseudoAtoms[k][i][j][4]);

                free(SpectrumPseudoAtomsAverage[k][i][j][0]);
                free(SpectrumPseudoAtomsAverage[k][i][j][1]);
                free(SpectrumPseudoAtomsAverage[k][i][j][2]);
                free(SpectrumPseudoAtomsAverage[k][i][j][3]);
                free(SpectrumPseudoAtomsAverage[k][i][j][4]);

                free(SpectrumPseudoAtoms[k][i][j]);
                free(SpectrumPseudoAtomsAverage[k][i][j]);
              }
            }
            free(SpectrumPseudoAtoms[k][i]);
            free(SpectrumPseudoAtomsAverage[k][i]);
          }
          free(Spectrum[k]);
          free(SpectrumAverage[k]);
          free(UnweightedSpectrum[k]);
          free(UnweightedSpectrumAverage[k]);
          free(SpectrumPseudoAtoms[k]);
          free(SpectrumPseudoAtomsAverage[k]);
          free(SpectrumCount[k]);
        }
      }
      free(corr_pseudoatoms);
      free(corrvacf_pseudoatoms);
      free(corr_pseudoatoms_isotropic);
      free(corrvacf_pseudoatoms_isotropic);
      free(Spectrum);
      free(SpectrumAverage);
      free(UnweightedSpectrum);
      free(UnweightedSpectrumAverage);
      free(SpectrumPseudoAtoms);
      free(SpectrumPseudoAtomsAverage);
      free(SpectrumCount);
      free(sumw);
      free(SpectrumBuffer);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleMeanSquaredDisplacementOrderN                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the mean-squared displacement using a modified order-N algorithm.                *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

int CalculateSiteType(VECTOR pos)
{
  int type;
  int closest;
  REAL distance;

  if (Framework[CurrentSystem].NumberOfIons==0) return 0;

  ClosestCrystallographicPosition2(pos,&closest,&distance);

  type = Framework[CurrentSystem].Ions[closest].AssymetricType;

  if(distance>CutOffIons) return 0;

  return type+1;
}

void SampleMeanSquaredDisplacementOrderN(int Switch)
{
  int i,j,k,l,CurrentBlock,index;
  int CurrentBlocklength,type,shift;
  VECTOR value,drift,com;
  int site_type;
  static VECTOR *value_onsager;
  FILE *FilePtr;
  char buffer[256];
  REAL fac,dt,value_dir_avg,count;

  switch(Switch)
  {
    case ALLOCATE:
      CountMSDOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
      NumberOfBlocksMSDOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
      value_onsager=(VECTOR*)calloc(NumberOfComponents+1,sizeof(VECTOR));

      BlockLengthMSDOrderN=(int**)calloc(NumberOfSystems,sizeof(int*));
      MsdOrderNCount=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      MsdOrderNOnsagerCount=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      MsdOrderNCountPerMolecule=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      MsdOrderNDirAvg=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      MsdOrderNOnsagerDirAvg=(REAL*****)calloc(NumberOfSystems,sizeof(REAL****));
      MsdOrderNPerMoleculeDirAvg=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      MsdOrderN=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      MsdOrderNOnsager=(VECTOR*****)calloc(NumberOfSystems,sizeof(VECTOR****));
      MsdOrderNPerMolecule=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      BlockDataMSDOrderN=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      BlockDataMSDOrderNOnsager=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      MsdOrderNTotalOnsager=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      MsdOrderNTotalOnsagerDirAvg=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      MsdOrderNTotalOnsagerCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

      MsdOrderNCountPerSiteType=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      MsdOrderNPerSiteType=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      MsdOrderNPerSiteTypeDirAvg=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      BlockDataSiteTypeMSDOrderN=(int****)calloc(NumberOfSystems,sizeof(int***));

      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeMSDOrderN[i])
        {
          BlockLengthMSDOrderN[i]=(int*)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(int));
          MsdOrderNCount[i]=(REAL***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(REAL**));
          MsdOrderNOnsagerCount[i]=(REAL***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(REAL**));
          MsdOrderNDirAvg[i]=(REAL***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(REAL**));
          MsdOrderNOnsagerDirAvg[i]=(REAL****)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(REAL***));
          MsdOrderN[i]=(VECTOR***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(VECTOR**));
          MsdOrderNOnsager[i]=(VECTOR****)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(VECTOR***));
          BlockDataMSDOrderN[i]=(VECTOR***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(VECTOR**));
          BlockDataMSDOrderNOnsager[i]=(VECTOR***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(VECTOR**));
          MsdOrderNTotalOnsager[i]=(VECTOR**)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(VECTOR*));
          MsdOrderNTotalOnsagerDirAvg[i]=(REAL**)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(REAL*));
          MsdOrderNTotalOnsagerCount[i]=(REAL**)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(REAL*));

          if(ComputeIndividualMSDOrderN)
          {
            MsdOrderNCountPerMolecule[i]=(REAL***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(REAL**));
            MsdOrderNPerMoleculeDirAvg[i]=(REAL***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(REAL**));
            MsdOrderNPerMolecule[i]=(VECTOR***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(VECTOR**));
          }
          if(ComputeSiteTypeMSDOrderN)
          {
            MsdOrderNCountPerSiteType[i]=(REAL***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(REAL**));
            MsdOrderNPerSiteTypeDirAvg[i]=(REAL***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(REAL**));
            MsdOrderNPerSiteType[i]=(VECTOR***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(VECTOR**));
            BlockDataSiteTypeMSDOrderN[i]=(int***)calloc(MaxNumberOfBlocksMSDOrderN,sizeof(int**));
          }

          for(j=0;j<MaxNumberOfBlocksMSDOrderN;j++)
          {
            BlockDataMSDOrderN[i][j]=(VECTOR**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR*));
            for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              BlockDataMSDOrderN[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(VECTOR));

            if(ComputeSiteTypeMSDOrderN)
            {
              BlockDataSiteTypeMSDOrderN[i][j]=(int**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(int*));
              for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
                BlockDataSiteTypeMSDOrderN[i][j][k]=(int*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(int));
            }

            MsdOrderNCount[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            MsdOrderNDirAvg[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            MsdOrderN[i][j]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
            for(k=0;k<NumberOfComponents;k++)
            {
              MsdOrderNCount[i][j][k]=(REAL*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(REAL));
              MsdOrderNDirAvg[i][j][k]=(REAL*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(REAL));
              MsdOrderN[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(VECTOR));
            }

            if(ComputeIndividualMSDOrderN)
            {
              MsdOrderNCountPerMolecule[i][j]=(REAL**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(REAL*));
              MsdOrderNPerMoleculeDirAvg[i][j]=(REAL**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(REAL*));
              MsdOrderNPerMolecule[i][j]=(VECTOR**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR*));
              for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              {
                MsdOrderNCountPerMolecule[i][j][k]=(REAL*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(REAL));
                MsdOrderNPerMolecule[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(VECTOR));
                MsdOrderNPerMoleculeDirAvg[i][j][k]=(REAL*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(REAL));
              }
            }

            if(ComputeSiteTypeMSDOrderN)
            {
              MsdOrderNCountPerSiteType[i][j]=(REAL**)calloc(NumberOfSitesMSDOrderN[i],sizeof(REAL*));
              MsdOrderNPerSiteTypeDirAvg[i][j]=(REAL**)calloc(NumberOfSitesMSDOrderN[i],sizeof(REAL*));
              MsdOrderNPerSiteType[i][j]=(VECTOR**)calloc(NumberOfSitesMSDOrderN[i],sizeof(VECTOR*));
              for(k=0;k<NumberOfSitesMSDOrderN[i];k++)
              {
                MsdOrderNCountPerSiteType[i][j][k]=(REAL*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(REAL));
                MsdOrderNPerSiteType[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(VECTOR));
                MsdOrderNPerSiteTypeDirAvg[i][j][k]=(REAL*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(REAL));
              }
            }

            // Onsager data
            BlockDataMSDOrderNOnsager[i][j]=(VECTOR**)calloc(NumberOfComponents+1,sizeof(VECTOR*));
            MsdOrderNTotalOnsager[i][j]=(VECTOR*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(VECTOR));
            MsdOrderNTotalOnsagerDirAvg[i][j]=(REAL*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(REAL));
            MsdOrderNTotalOnsagerCount[i][j]=(REAL*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(REAL));
            MsdOrderNOnsager[i][j]=(VECTOR***)calloc(NumberOfComponents,sizeof(VECTOR**));
            MsdOrderNOnsagerDirAvg[i][j]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));
            MsdOrderNOnsagerCount[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

            for(k=0;k<=NumberOfComponents;k++)
              BlockDataMSDOrderNOnsager[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(VECTOR));

            for(k=0;k<NumberOfComponents;k++)
            {
              MsdOrderNOnsagerCount[i][j][k]=(REAL*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(REAL));
              MsdOrderNOnsagerDirAvg[i][j][k]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
              MsdOrderNOnsager[i][j][k]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
              for(l=0;l<NumberOfComponents;l++)
              {
                MsdOrderNOnsager[i][j][k][l]=(VECTOR*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(VECTOR));
                MsdOrderNOnsagerDirAvg[i][j][k][l]=(REAL*)calloc(NumberOfBlockElementsMSDOrderN,sizeof(REAL));
              }
            }
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      // return if the msd does not has to be calculated for this system
      if(!(ComputeMSDOrderN[CurrentSystem]&&(CurrentCycle%SampleMSDOrderNEvery[CurrentSystem]==0))) return;

      // compute drift of the system
      if(Framework[CurrentSystem].FrameworkModel==NONE)
      {
        com=GetCenterOfMassCurrentSystem();
        drift.x=com.x-IntialCenterOfMassPosition[CurrentSystem].x;
        drift.y=com.y-IntialCenterOfMassPosition[CurrentSystem].y;
        drift.z=com.z-IntialCenterOfMassPosition[CurrentSystem].z;
      }
      else if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        com=GetFrameworkCenterOfMass();
        drift.x=com.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x;
        drift.y=com.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y;
        drift.z=com.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z;
      }
      else
      {
        drift.x=0.0;
        drift.y=0.0;
        drift.z=0.0;
      }

      // determine current number of blocks
      NumberOfBlocksMSDOrderN[CurrentSystem]=1;
      i=CountMSDOrderN[CurrentSystem]/NumberOfBlockElementsMSDOrderN;
      while(i!=0)
      {
        NumberOfBlocksMSDOrderN[CurrentSystem]++;
        i/=NumberOfBlockElementsMSDOrderN;
      }

      // ignore everything beyond the last block
      NumberOfBlocksMSDOrderN[CurrentSystem]=MIN2(NumberOfBlocksMSDOrderN[CurrentSystem],MaxNumberOfBlocksMSDOrderN);

      for(CurrentBlock=0;CurrentBlock<NumberOfBlocksMSDOrderN[CurrentSystem];CurrentBlock++)
      {
        // test for blocking operation: CountMSDOrderN is a multiple of NumberOfBlockElementsMSDOrderN^CurrentBlock
        if((CountMSDOrderN[CurrentSystem])%((int)pow((REAL)NumberOfBlockElementsMSDOrderN,CurrentBlock))==0)
        {
          // increase the current block-length
          BlockLengthMSDOrderN[CurrentSystem][CurrentBlock]++;

          // limit length to NumberOfBlockElementsMSDOrderN
          CurrentBlocklength=MIN2(BlockLengthMSDOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsMSDOrderN);

          // Self diffusion
          // ========================================================================================================================
          shift=NumberOfAdsorbateMolecules[CurrentSystem];
          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];i++)
          {
            if(i<NumberOfAdsorbateMolecules[CurrentSystem])
              value=GetAdsorbateCenterOfMass(i);
            else
              value=GetCationCenterOfMass(i-shift);
            value.x-=drift.x;
            value.y-=drift.y;
            value.z-=drift.z;

            for(j=CurrentBlocklength-1;j>0;j--)
              BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j]=BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j-1];
            BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][0]=value;

            if(ComputeSiteTypeMSDOrderN)
            {

              for(j=CurrentBlocklength-1;j>0;j--)
                BlockDataSiteTypeMSDOrderN[CurrentSystem][CurrentBlock][i][j]=BlockDataSiteTypeMSDOrderN[CurrentSystem][CurrentBlock][i][j-1];

              // set the type of the site for this measurement-point
              site_type=CalculateSiteType(value);
              BlockDataSiteTypeMSDOrderN[CurrentSystem][CurrentBlock][i][0]=site_type;
            }

            // get the type of the molecule
            if(i<NumberOfAdsorbateMolecules[CurrentSystem])
              type=Adsorbates[CurrentSystem][i].Type;
            else
              type=Cations[CurrentSystem][i-shift].Type;

            for(j=0;j<CurrentBlocklength;j++)
            {
              // msd for each component
              MsdOrderNCount[CurrentSystem][CurrentBlock][type][j]+=1.0;
              MsdOrderN[CurrentSystem][CurrentBlock][type][j].x+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].x-value.x);
              MsdOrderN[CurrentSystem][CurrentBlock][type][j].y+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].y-value.y);
              MsdOrderN[CurrentSystem][CurrentBlock][type][j].z+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].z-value.z);
              MsdOrderNDirAvg[CurrentSystem][CurrentBlock][type][j]+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].x-value.x)+
                                                                     SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].y-value.y)+
                                                                     SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].z-value.z);

              if(ComputeIndividualMSDOrderN)
              {
                MsdOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j]+=1.0;
                MsdOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].x+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].x-value.x);
                MsdOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].y+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].y-value.y);
                MsdOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].z+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].z-value.z);
                MsdOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][i][j]+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].x-value.x)+
                                                                               SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].y-value.y)+
                                                                               SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].z-value.z);
              }

              if(ComputeSiteTypeMSDOrderN)
              {
                if(site_type>=NumberOfSitesMSDOrderN[CurrentSystem])
                {
                  printf("ERROR: site-type %d > maximum number of sites %d in routine 'SampleMeanSquaredDisplacementOrderN'\n",site_type,NumberOfSitesMSDOrderN[CurrentSystem]);
                  exit(0);
                }

                // only count when the site is currently the same as it was at the acf-origin
                if (BlockDataSiteTypeMSDOrderN[CurrentSystem][CurrentBlock][i][j] == site_type)
                {
                  MsdOrderNCountPerSiteType[CurrentSystem][CurrentBlock][site_type][j]+=1.0;
                  MsdOrderNPerSiteType[CurrentSystem][CurrentBlock][site_type][j].x+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].x-value.x);
                  MsdOrderNPerSiteType[CurrentSystem][CurrentBlock][site_type][j].y+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].y-value.y);
                  MsdOrderNPerSiteType[CurrentSystem][CurrentBlock][site_type][j].z+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].z-value.z);
                  MsdOrderNPerSiteTypeDirAvg[CurrentSystem][CurrentBlock][site_type][j]+=SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].x-value.x)+
                                                                                         SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].y-value.y)+
                                                                                         SQR(BlockDataMSDOrderN[CurrentSystem][CurrentBlock][i][j].z-value.z);
                }
              }
            }
          }

          // Onsager diffusion
          // ========================================================================================================================

          for(i=0;i<=NumberOfComponents;i++)
            value_onsager[i].x=value_onsager[i].y=value_onsager[i].z=0.0;

          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
          {
            type=Adsorbates[CurrentSystem][i].Type;
            value=GetAdsorbateCenterOfMass(i);
            value.x-=drift.x;
            value.y-=drift.y;
            value.z-=drift.z;

            // add per component
            value_onsager[type].x+=value.x;
            value_onsager[type].y+=value.y;
            value_onsager[type].z+=value.z;

            // add to total
            value_onsager[NumberOfComponents].x+=value.x;
            value_onsager[NumberOfComponents].y+=value.y;
            value_onsager[NumberOfComponents].z+=value.z;
          }

          for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
          {
            type=Cations[CurrentSystem][i].Type;
            value=GetCationCenterOfMass(i);
            value.x-=drift.x;
            value.y-=drift.y;
            value.z-=drift.z;

            // add per component
            value_onsager[type].x+=value.x;
            value_onsager[type].y+=value.y;
            value_onsager[type].z+=value.z;

            // add to total
            value_onsager[NumberOfComponents].x+=value.x;
            value_onsager[NumberOfComponents].y+=value.y;
            value_onsager[NumberOfComponents].z+=value.z;
          }


          for(i=0;i<=NumberOfComponents;i++)
          {
            for(j=CurrentBlocklength-1;j>0;j--)
              BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][i][j]=BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][i][j-1];
            BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][i][0]=value_onsager[i];
          }

          // msd for each component
          for(k=0;k<CurrentBlocklength;k++)
          {
            for(i=0;i<NumberOfComponents;i++)
            {
              MsdOrderNOnsagerCount[CurrentSystem][CurrentBlock][i][k]+=1.0;
              for(j=0;j<NumberOfComponents;j++)
              {
                MsdOrderNOnsager[CurrentSystem][CurrentBlock][i][j][k].x+=(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][i][k].x-value_onsager[i].x)*
                                                                          (BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][j][k].x-value_onsager[j].x);
                MsdOrderNOnsager[CurrentSystem][CurrentBlock][i][j][k].y+=(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][i][k].y-value_onsager[i].y)*
                                                                          (BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][j][k].y-value_onsager[j].y);
                MsdOrderNOnsager[CurrentSystem][CurrentBlock][i][j][k].z+=(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][i][k].z-value_onsager[i].z)*
                                                                          (BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][j][k].z-value_onsager[j].z);
                MsdOrderNOnsagerDirAvg[CurrentSystem][CurrentBlock][i][j][k]+=(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][i][k].x-value_onsager[i].x)*
                                                                              (BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][j][k].x-value_onsager[j].x)
                                                                             +(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][i][k].y-value_onsager[i].y)*
                                                                              (BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][j][k].y-value_onsager[j].y)
                                                                             +(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][i][k].z-value_onsager[i].z)*
                                                                              (BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][j][k].z-value_onsager[j].z);
              }
            }
          }
          // msd for the total fluid
          for(k=0;k<CurrentBlocklength;k++)
          {
            MsdOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][k]+=1.0;
            MsdOrderNTotalOnsager[CurrentSystem][CurrentBlock][k].x+=SQR(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].x-
                                                                         value_onsager[NumberOfComponents].x);
            MsdOrderNTotalOnsager[CurrentSystem][CurrentBlock][k].y+=SQR(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].y-
                                                                         value_onsager[NumberOfComponents].y);
            MsdOrderNTotalOnsager[CurrentSystem][CurrentBlock][k].z+=SQR(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].z-
                                                                         value_onsager[NumberOfComponents].z);
            MsdOrderNTotalOnsagerDirAvg[CurrentSystem][CurrentBlock][k]+=
               SQR(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].x-value_onsager[NumberOfComponents].x)+
               SQR(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].y-value_onsager[NumberOfComponents].y)+
               SQR(BlockDataMSDOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].z-value_onsager[NumberOfComponents].z);
          }
        }
      }
      // CountMSDOrderN the current sampling
      CountMSDOrderN[CurrentSystem]++;
      break;
    case PRINT:
      // return if the msd does not has to be calculated for this system
      if((!ComputeMSDOrderN[CurrentSystem])||(CurrentCycle%WriteMSDOrderNEvery[CurrentSystem]!=0)) return;

      mkdir("MSDOrderN",S_IRWXU);

      sprintf(buffer,"MSDOrderN/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      // Self diffusion
      // ========================================================================================================================

      sprintf(buffer,"MSDOrderN/System_%d/msd_self_total%s.dat",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"# column 1: time [ps]\n");
      fprintf(FilePtr,"# column 2: msd xyz [A^2]\n");
      fprintf(FilePtr,"# column 3: msd x [A^2]\n");
      fprintf(FilePtr,"# column 4: msd y [A^2]\n");
      fprintf(FilePtr,"# column 5: msd z [A^2]\n");
      fprintf(FilePtr,"# column 6: number of samples [-]\n");

      for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksMSDOrderN,NumberOfBlocksMSDOrderN[CurrentSystem]);CurrentBlock++)
      {
        CurrentBlocklength=MIN2(BlockLengthMSDOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsMSDOrderN);
        dt=SampleMSDOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsMSDOrderN,CurrentBlock);
        for(j=1;j<CurrentBlocklength;j++)
        {
          count=0.0;
          value_dir_avg=0.0;
          value.x=value.y=value.z=0.0;
          for(i=0;i<NumberOfComponents;i++)
          {
            value.x+=MsdOrderN[CurrentSystem][CurrentBlock][i][j].x;
            value.y+=MsdOrderN[CurrentSystem][CurrentBlock][i][j].y;
            value.z+=MsdOrderN[CurrentSystem][CurrentBlock][i][j].z;
            value_dir_avg+=MsdOrderNDirAvg[CurrentSystem][CurrentBlock][i][j];
            count+=MsdOrderNCount[CurrentSystem][CurrentBlock][i][j];
          }

          if(count>0.0)
          {
            fac=1.0/count;

            fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                (double)(j*dt),
                (double)(fac*value_dir_avg),
                (double)(fac*value.x),
                (double)(fac*value.y),
                (double)(fac*value.z),
                (double)count);
          }
        }
      }

      fclose(FilePtr);


      // print out msd per component
      for(i=0;i<NumberOfComponents;i++)
      {
        sprintf(buffer,"MSDOrderN/System_%d/msd_self_%s_%d%s.dat",
                CurrentSystem,Components[i].Name,i,FileNameAppend);

        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# column 1: time [ps]\n");
        fprintf(FilePtr,"# column 2: msd xyz [A^2]\n");
        fprintf(FilePtr,"# column 3: msd x [A^2]\n");
        fprintf(FilePtr,"# column 4: msd y [A^2]\n");
        fprintf(FilePtr,"# column 5: msd z [A^2]\n");
        fprintf(FilePtr,"# column 6: number of samples [-]\n");


        // write msd results averaged per component to file
        for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksMSDOrderN,NumberOfBlocksMSDOrderN[CurrentSystem]);CurrentBlock++)
        {
          CurrentBlocklength=MIN2(BlockLengthMSDOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsMSDOrderN);
          dt=SampleMSDOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsMSDOrderN,CurrentBlock);
          for(j=1;j<CurrentBlocklength;j++)
          {
            if(MsdOrderNCount[CurrentSystem][CurrentBlock][i][j]>0.0)
            {
              fac=1.0/MsdOrderNCount[CurrentSystem][CurrentBlock][i][j];

              fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                  (double)(j*dt),
                  (double)(fac*MsdOrderNDirAvg[CurrentSystem][CurrentBlock][i][j]),
                  (double)(fac*MsdOrderN[CurrentSystem][CurrentBlock][i][j].x),
                  (double)(fac*MsdOrderN[CurrentSystem][CurrentBlock][i][j].y),
                  (double)(fac*MsdOrderN[CurrentSystem][CurrentBlock][i][j].z),
                  (double)MsdOrderNCount[CurrentSystem][CurrentBlock][i][j]);
            }
          }
        }

        fclose(FilePtr);
      }

      if(ComputeIndividualMSDOrderN)
      {
        sprintf(buffer,"MSDOrderN/System_%d/MSDOrderN_per_adsorbate",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        // print out msd per adsorbate
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          sprintf(buffer,"MSDOrderN/System_%d/MSDOrderN_per_adsorbate/msd_self_%d%s.dat",
                  CurrentSystem,i,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: msd xyz [A^2]\n");
          fprintf(FilePtr,"# column 3: msd x [A^2]\n");
          fprintf(FilePtr,"# column 4: msd y [A^2]\n");
          fprintf(FilePtr,"# column 5: msd z [A^2]\n");
          fprintf(FilePtr,"# column 6: number of samples [-]\n");

          // write results to file
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksMSDOrderN,NumberOfBlocksMSDOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthMSDOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsMSDOrderN);
            dt=SampleMSDOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsMSDOrderN,CurrentBlock);
            for(j=1;j<CurrentBlocklength;j++)
            {
              if(MsdOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j]>0.0)
              {
                fac=1.0/MsdOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j];
                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(j*dt),
                    (double)(fac*MsdOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][i][j]),
                    (double)(fac*MsdOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].x),
                    (double)(fac*MsdOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].y),
                    (double)(fac*MsdOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].z),
                    (double)MsdOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j]);
              }
            }
          }

          fclose(FilePtr);
        }

        sprintf(buffer,"MSDOrderN/System_%d/MSDOrderN_per_cation",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        // print out msd per cation
        shift=NumberOfAdsorbateMolecules[CurrentSystem];
        for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
        {
          index=i+shift;
          sprintf(buffer,"MSDOrderN/System_%d/MSDOrderN_per_cation/msd_self_%d%s.dat",
                  CurrentSystem,i,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: msd xyz [A^2]\n");
          fprintf(FilePtr,"# column 3: msd x [A^2]\n");
          fprintf(FilePtr,"# column 4: msd y [A^2]\n");
          fprintf(FilePtr,"# column 5: msd z [A^2]\n");
          fprintf(FilePtr,"# column 6: number of samples [-]\n");

          // write results to file
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksMSDOrderN,NumberOfBlocksMSDOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthMSDOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsMSDOrderN);
            dt=SampleMSDOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsMSDOrderN,CurrentBlock);
            for(j=1;j<CurrentBlocklength;j++)
            {
              if(MsdOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][j]>0.0)
              {
                fac=1.0/MsdOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][j];
                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(j*dt),
                    (double)(fac*MsdOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][index][j]),
                    (double)(fac*MsdOrderNPerMolecule[CurrentSystem][CurrentBlock][index][j].x),
                    (double)(fac*MsdOrderNPerMolecule[CurrentSystem][CurrentBlock][index][j].y),
                    (double)(fac*MsdOrderNPerMolecule[CurrentSystem][CurrentBlock][index][j].z),
                    (double)MsdOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][j]);
              }
            }
          }

          fclose(FilePtr);
        }

      }

      if(ComputeSiteTypeMSDOrderN)
      {
        sprintf(buffer,"MSDOrderN/System_%d/MSDOrderN_per_site",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        // print out msd per adsorbate
        for(i=0;i<NumberOfSitesMSDOrderN[CurrentSystem];i++)
        {
          sprintf(buffer,"MSDOrderN/System_%d/MSDOrderN_per_site/msd_self_%d%s.dat",
                  CurrentSystem,i,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: msd xyz [A^2]\n");
          fprintf(FilePtr,"# column 3: msd x [A^2]\n");
          fprintf(FilePtr,"# column 4: msd y [A^2]\n");
          fprintf(FilePtr,"# column 5: msd z [A^2]\n");
          fprintf(FilePtr,"# column 6: number of samples [-]\n");

          // write results to file
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksMSDOrderN,NumberOfBlocksMSDOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthMSDOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsMSDOrderN);
            dt=SampleMSDOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsMSDOrderN,CurrentBlock);
            for(j=1;j<CurrentBlocklength;j++)
            {
              if(MsdOrderNCountPerSiteType[CurrentSystem][CurrentBlock][i][j]>0.0)
              {
                fac=1.0/MsdOrderNCountPerSiteType[CurrentSystem][CurrentBlock][i][j];
                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(j*dt),
                    (double)(fac*MsdOrderNPerSiteTypeDirAvg[CurrentSystem][CurrentBlock][i][j]),
                    (double)(fac*MsdOrderNPerSiteType[CurrentSystem][CurrentBlock][i][j].x),
                    (double)(fac*MsdOrderNPerSiteType[CurrentSystem][CurrentBlock][i][j].y),
                    (double)(fac*MsdOrderNPerSiteType[CurrentSystem][CurrentBlock][i][j].z),
                    (double)MsdOrderNCountPerSiteType[CurrentSystem][CurrentBlock][i][j]);
              }
            }
          }

          fclose(FilePtr);
        }
      }

      // Onsager diffusion
      // ========================================================================================================================

      // write msd results for the fluid to file
      sprintf(buffer,"MSDOrderN/System_%d/msd_onsager_total%s.dat",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"# column 1: time [ps]\n");
      fprintf(FilePtr,"# column 2: msd xyz [A^2]\n");
      fprintf(FilePtr,"# column 3: msd x [A^2]\n");
      fprintf(FilePtr,"# column 4: msd y [A^2]\n");
      fprintf(FilePtr,"# column 5: msd z [A^2]\n");
      fprintf(FilePtr,"# column 6: number of samples [-]\n");

      for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksMSDOrderN,NumberOfBlocksMSDOrderN[CurrentSystem]);CurrentBlock++)
      {
        CurrentBlocklength=MIN2(BlockLengthMSDOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsMSDOrderN);
        dt=SampleMSDOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsMSDOrderN,CurrentBlock);
        for(l=1;l<CurrentBlocklength;l++)
        {
          if(MsdOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][l]>0.0)
          {
             fac=1.0/((NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem])*MsdOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][l]);
             fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                (double)(l*dt),
                (double)(fac*MsdOrderNTotalOnsagerDirAvg[CurrentSystem][CurrentBlock][l]),
                (double)(fac*MsdOrderNTotalOnsager[CurrentSystem][CurrentBlock][l].x),
                (double)(fac*MsdOrderNTotalOnsager[CurrentSystem][CurrentBlock][l].y),
                (double)(fac*MsdOrderNTotalOnsager[CurrentSystem][CurrentBlock][l].z),
                (double)MsdOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][l]);
          }
        }
      }
      fclose(FilePtr);

      // write the msd per component
      for(j=0;j<NumberOfComponents;j++)
      {
        for(k=0;k<NumberOfComponents;k++)
        {
          sprintf(buffer,"MSDOrderN/System_%d/msd_onsager_%s_%d_%s_%d%s.dat",
                  CurrentSystem,Components[j].Name,j,Components[k].Name,k,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: msd xyz [A^2]\n");
          fprintf(FilePtr,"# column 3: msd x [A^2]\n");
          fprintf(FilePtr,"# column 4: msd y [A^2]\n");
          fprintf(FilePtr,"# column 5: msd z [A^2]\n");
          fprintf(FilePtr,"# column 6: number of samples [-]\n");

          // write msd results averaged per component to file
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksMSDOrderN,NumberOfBlocksMSDOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthMSDOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsMSDOrderN);
            dt=SampleMSDOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsMSDOrderN,CurrentBlock);
            for(l=1;l<CurrentBlocklength;l++)
            {
              if(MsdOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][l]>0.0)
              {
                fac=1.0/(Components[j].NumberOfMolecules[CurrentSystem]*MsdOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][l]);
                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(l*dt),
                    (double)(fac*MsdOrderNOnsagerDirAvg[CurrentSystem][CurrentBlock][j][k][l]),
                    (double)(fac*MsdOrderNOnsager[CurrentSystem][CurrentBlock][j][k][l].x),
                    (double)(fac*MsdOrderNOnsager[CurrentSystem][CurrentBlock][j][k][l].y),
                    (double)(fac*MsdOrderNOnsager[CurrentSystem][CurrentBlock][j][k][l].z),
                    (double)MsdOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][l]);
              }
            }
          }

          fclose(FilePtr);
        }
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeMSDOrderN[i])
        {
          for(j=0;j<MaxNumberOfBlocksMSDOrderN;j++)
          {
            for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              free(BlockDataMSDOrderN[i][j][k]);

            free(BlockDataMSDOrderN[i][j]);

            for(k=0;k<NumberOfComponents;k++)
            {
              free(MsdOrderNCount[i][j][k]);
              free(MsdOrderNDirAvg[i][j][k]);
              free(MsdOrderN[i][j][k]);
            }
            free(MsdOrderNCount[i][j]);
            free(MsdOrderNDirAvg[i][j]);
            free(MsdOrderN[i][j]);

            if(ComputeIndividualMSDOrderN)
            {
              for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              {
                free(MsdOrderNCountPerMolecule[i][j][k]);
                free(MsdOrderNPerMolecule[i][j][k]);
                free(MsdOrderNPerMoleculeDirAvg[i][j][k]);
              }
              free(MsdOrderNCountPerMolecule[i][j]);
              free(MsdOrderNPerMoleculeDirAvg[i][j]);
              free(MsdOrderNPerMolecule[i][j]);
            }

            if(ComputeSiteTypeMSDOrderN)
            {
              for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
                free(BlockDataSiteTypeMSDOrderN[i][j][k]);

              free(BlockDataSiteTypeMSDOrderN[i][j]);

              for(k=0;k<NumberOfSitesMSDOrderN[i];k++)
              {
                free(MsdOrderNCountPerSiteType[i][j][k]);
                free(MsdOrderNPerSiteType[i][j][k]);
                free(MsdOrderNPerSiteTypeDirAvg[i][j][k]);
              }
              free(MsdOrderNCountPerSiteType[i][j]);
              free(MsdOrderNPerSiteTypeDirAvg[i][j]);
              free(MsdOrderNPerSiteType[i][j]);
            }

            // Onsager data
            for(k=0;k<NumberOfComponents;k++)
            {
              for(l=0;l<NumberOfComponents;l++)
              {
                free(MsdOrderNOnsager[i][j][k][l]);
                free(MsdOrderNOnsagerDirAvg[i][j][k][l]);
              }
              free(MsdOrderNOnsagerCount[i][j][k]);
              free(MsdOrderNOnsagerDirAvg[i][j][k]);
              free(MsdOrderNOnsager[i][j][k]);
            }
            for(k=0;k<=NumberOfComponents;k++)
              free(BlockDataMSDOrderNOnsager[i][j][k]);

            free(BlockDataMSDOrderNOnsager[i][j]);
            free(MsdOrderNTotalOnsager[i][j]);
            free(MsdOrderNTotalOnsagerDirAvg[i][j]);
            free(MsdOrderNTotalOnsagerCount[i][j]);
            free(MsdOrderNOnsager[i][j]);
            free(MsdOrderNOnsagerDirAvg[i][j]);
            free(MsdOrderNOnsagerCount[i][j]);
          }

          if(ComputeIndividualMSDOrderN)
          {
            free(MsdOrderNCountPerMolecule[i]);
            free(MsdOrderNPerMoleculeDirAvg[i]);
            free(MsdOrderNPerMolecule[i]);
          }

          if(ComputeIndividualMSDOrderN)
          {
            free(MsdOrderNCountPerSiteType[i]);
            free(MsdOrderNPerSiteTypeDirAvg[i]);
            free(MsdOrderNPerSiteType[i]);
          }

          free(BlockLengthMSDOrderN[i]);
          free(MsdOrderNCount[i]);
          free(MsdOrderNOnsagerCount[i]);
          free(MsdOrderNDirAvg[i]);
          free(MsdOrderNOnsagerDirAvg[i]);
          free(MsdOrderN[i]);
          free(MsdOrderNOnsager[i]);
          free(BlockDataMSDOrderN[i]);
          free(BlockDataSiteTypeMSDOrderN[i]);
          free(BlockDataMSDOrderNOnsager[i]);
          free(MsdOrderNTotalOnsager[i]);
          free(MsdOrderNTotalOnsagerDirAvg[i]);
          free(MsdOrderNTotalOnsagerCount[i]);
        }
      }
      free(CountMSDOrderN);
      free(NumberOfBlocksMSDOrderN);
      free(value_onsager);
      free(BlockLengthMSDOrderN);
      free(MsdOrderNCount);
      free(MsdOrderNOnsagerCount);
      free(MsdOrderNCountPerMolecule);
      free(MsdOrderNCountPerSiteType);
      free(MsdOrderNDirAvg);
      free(MsdOrderNOnsagerDirAvg);
      free(MsdOrderNPerMoleculeDirAvg);
      free(MsdOrderNPerSiteTypeDirAvg);
      free(MsdOrderN);
      free(MsdOrderNOnsager);
      free(MsdOrderNPerMolecule);
      free(MsdOrderNPerSiteType);
      free(BlockDataMSDOrderN);
      free(BlockDataSiteTypeMSDOrderN);
      free(BlockDataMSDOrderNOnsager);
      free(MsdOrderNTotalOnsager);
      free(MsdOrderNTotalOnsagerDirAvg);
      free(MsdOrderNTotalOnsagerCount);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleVelocityAutoCorrelationFunctionOrderN                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the velocity autocorrelation function using a modified order-N algorithm.        *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleVelocityAutoCorrelationFunctionOrderN(int Switch)
{
  int i,j,k,l,CurrentBlock,index;
  int CurrentBlocklength,type,shift;
  VECTOR value,drift;
  static VECTOR *value_onsager;
  FILE *FilePtr;
  char buffer[256];
  REAL fac,Integral,D,dt,count,value_dir_avg;

  switch(Switch)
  {
    case ALLOCATE:
      CountVACFOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
      NumberOfBlocksVACFOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
      value_onsager=(VECTOR*)calloc(NumberOfComponents+1,sizeof(VECTOR));

      BlockLengthVACFOrderN=(int**)calloc(NumberOfSystems,sizeof(int*));
      VacfOrderNCount=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      VacfOrderNOnsagerCount=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      VacfOrderNCountPerMolecule=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      VacfOrderNDirAvg=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      VacfOrderNOnsagerDirAvg=(REAL*****)calloc(NumberOfSystems,sizeof(REAL****));
      VacfOrderNPerMoleculeDirAvg=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      VacfOrderN=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      VacfOrderNOnsager=(VECTOR*****)calloc(NumberOfSystems,sizeof(VECTOR****));
      VacfOrderNPerMolecule=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      BlockDataVACFOrderN=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      BlockDataVACFOrderNOnsager=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      VacfOrderNTotalOnsager=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      VacfOrderNTotalOnsagerDirAvg=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      VacfOrderNTotalOnsagerCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeVACFOrderN[i])
        {
          BlockLengthVACFOrderN[i]=(int*)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(int));
          VacfOrderNCount[i]=(REAL***)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(REAL**));
          VacfOrderNOnsagerCount[i]=(REAL***)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(REAL**));
          VacfOrderNDirAvg[i]=(REAL***)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(REAL**));
          VacfOrderNOnsagerDirAvg[i]=(REAL****)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(REAL***));
          VacfOrderN[i]=(VECTOR***)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(VECTOR**));
          VacfOrderNOnsager[i]=(VECTOR****)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(VECTOR***));
          BlockDataVACFOrderN[i]=(VECTOR***)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(VECTOR**));
          BlockDataVACFOrderNOnsager[i]=(VECTOR***)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(VECTOR**));
          VacfOrderNTotalOnsager[i]=(VECTOR**)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(VECTOR*));
          VacfOrderNTotalOnsagerDirAvg[i]=(REAL**)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(REAL*));
          VacfOrderNTotalOnsagerCount[i]=(REAL**)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(REAL*));

          if(ComputeIndividualVACFOrderN)
          {
            VacfOrderNCountPerMolecule[i]=(REAL***)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(REAL**));
            VacfOrderNPerMoleculeDirAvg[i]=(REAL***)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(REAL**));
            VacfOrderNPerMolecule[i]=(VECTOR***)calloc(MaxNumberOfBlocksVACFOrderN,sizeof(VECTOR**));
          }

          for(j=0;j<MaxNumberOfBlocksVACFOrderN;j++)
          {
            BlockDataVACFOrderN[i][j]=(VECTOR**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR*));
            for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              BlockDataVACFOrderN[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(VECTOR));

            VacfOrderNCount[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            VacfOrderNDirAvg[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            VacfOrderN[i][j]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
            for(k=0;k<NumberOfComponents;k++)
            {
              VacfOrderNCount[i][j][k]=(REAL*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(REAL));
              VacfOrderNDirAvg[i][j][k]=(REAL*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(REAL));
              VacfOrderN[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(VECTOR));
            }

            if(ComputeIndividualVACFOrderN)
            {
              VacfOrderNCountPerMolecule[i][j]=(REAL**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(REAL*));
              VacfOrderNPerMoleculeDirAvg[i][j]=(REAL**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(REAL*));
              VacfOrderNPerMolecule[i][j]=(VECTOR**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR*));
              for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              {
                VacfOrderNCountPerMolecule[i][j][k]=(REAL*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(REAL));
                VacfOrderNPerMolecule[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(VECTOR));
                VacfOrderNPerMoleculeDirAvg[i][j][k]=(REAL*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(REAL));
              }
            }

            // Onsager data
            BlockDataVACFOrderNOnsager[i][j]=(VECTOR**)calloc(NumberOfComponents+1,sizeof(VECTOR*));
            VacfOrderNTotalOnsager[i][j]=(VECTOR*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(VECTOR));
            VacfOrderNTotalOnsagerDirAvg[i][j]=(REAL*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(REAL));
            VacfOrderNTotalOnsagerCount[i][j]=(REAL*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(REAL));
            VacfOrderNOnsagerCount[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            VacfOrderNOnsagerDirAvg[i][j]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));
            VacfOrderNOnsager[i][j]=(VECTOR***)calloc(NumberOfComponents,sizeof(VECTOR**));

            for(k=0;k<=NumberOfComponents;k++)
              BlockDataVACFOrderNOnsager[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(VECTOR));

            for(k=0;k<NumberOfComponents;k++)
            {
              VacfOrderNOnsagerCount[i][j][k]=(REAL*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(REAL));
              VacfOrderNOnsagerDirAvg[i][j][k]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
              VacfOrderNOnsager[i][j][k]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
              for(l=0;l<NumberOfComponents;l++)
              {
                VacfOrderNOnsager[i][j][k][l]=(VECTOR*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(VECTOR));
                VacfOrderNOnsagerDirAvg[i][j][k][l]=(REAL*)calloc(NumberOfBlockElementsVACFOrderN,sizeof(REAL));
              }
            }
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      // return if the vacf does not has to be calculated for this system
      if(!(ComputeVACFOrderN[CurrentSystem]&&(CurrentCycle%SampleVACFOrderNEvery[CurrentSystem]==0))) return;

      // compute drift of the system
      if(Framework[CurrentSystem].FrameworkModel==NONE)
        drift=GetCenterOfMassVelocityCurrentSystem();
      else if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
        drift=GetFrameworkCenterOfMassVelocity();
      else
      {
        drift.x=0.0;
        drift.y=0.0;
        drift.z=0.0;
      }

      // determine current number of blocks
      NumberOfBlocksVACFOrderN[CurrentSystem]=1;
      i=CountVACFOrderN[CurrentSystem]/NumberOfBlockElementsVACFOrderN;
      while(i!=0)
      {
        NumberOfBlocksVACFOrderN[CurrentSystem]++;
        i/=NumberOfBlockElementsVACFOrderN;
      }

      // ignore everything beyond the last block
      NumberOfBlocksVACFOrderN[CurrentSystem]=MIN2(NumberOfBlocksVACFOrderN[CurrentSystem],MaxNumberOfBlocksVACFOrderN);

      for(CurrentBlock=0;CurrentBlock<NumberOfBlocksVACFOrderN[CurrentSystem];CurrentBlock++)
      {
        // test for blocking operation: CountVACFOrderN is a multiple of NumberOfBlockElementsVACFOrderN^CurrentBlock
        if ((CountVACFOrderN[CurrentSystem])%((int)pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock))==0)
        {
          // increase the current block-length
          BlockLengthVACFOrderN[CurrentSystem][CurrentBlock]++;

          // limit length to NumberOfBlockElementsVACFOrderN
          CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);

          // Self diffusion
          // ========================================================================================================================

          shift=NumberOfAdsorbateMolecules[CurrentSystem];
          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];i++)
          {
            if(i<NumberOfAdsorbateMolecules[CurrentSystem])
              value=GetAdsorbateCenterOfMassVelocity(i);
            else
              value=GetCationCenterOfMassVelocity(i-shift);
            value.x-=drift.x;
            value.y-=drift.y;
            value.z-=drift.z;

            for(j=CurrentBlocklength-1;j>0;j--)
              BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j]=BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j-1];
            BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][0]=value;


            // get the type of the molecule
            if(i<NumberOfAdsorbateMolecules[CurrentSystem])
              type=Adsorbates[CurrentSystem][i].Type;
            else
              type=Cations[CurrentSystem][i-shift].Type;

            for(j=0;j<CurrentBlocklength;j++)
            {
              // vacf for each component
              VacfOrderNCount[CurrentSystem][CurrentBlock][type][j]+=1.0;
              VacfOrderN[CurrentSystem][CurrentBlock][type][j].x+=(BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].x*value.x);
              VacfOrderN[CurrentSystem][CurrentBlock][type][j].y+=(BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].y*value.y);
              VacfOrderN[CurrentSystem][CurrentBlock][type][j].z+=(BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].z*value.z);
              VacfOrderNDirAvg[CurrentSystem][CurrentBlock][type][j]+=(BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].x*value.x)+
                                                                      (BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].y*value.y)+
                                                                      (BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].z*value.z);

              if(ComputeIndividualVACFOrderN)
              {
                VacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j]+=1.0;
                VacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].x+=(BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].x*value.x);
                VacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].y+=(BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].y*value.y);
                VacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].z+=(BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].z*value.z);
                VacfOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][i][j]+=(BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].x*value.x)+
                                                                                (BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].y*value.y)+
                                                                                (BlockDataVACFOrderN[CurrentSystem][CurrentBlock][i][j].z*value.z);
              }
            }
          }

          // Onsager diffusion
          // ========================================================================================================================

          for(i=0;i<=NumberOfComponents;i++)
            value_onsager[i].x=value_onsager[i].y=value_onsager[i].z=0.0;

          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
          {
            type=Adsorbates[CurrentSystem][i].Type;

           if(i<NumberOfAdsorbateMolecules[CurrentSystem])
              value=GetAdsorbateCenterOfMassVelocity(i);
            else
              value=GetCationCenterOfMassVelocity(i-shift);
            value.x-=drift.x;
            value.y-=drift.y;
            value.z-=drift.z;

            value_onsager[type].x+=value.x;
            value_onsager[type].y+=value.y;
            value_onsager[type].z+=value.z;

            value_onsager[NumberOfComponents].x+=value.x;
            value_onsager[NumberOfComponents].y+=value.y;
            value_onsager[NumberOfComponents].z+=value.z;
          }

          for(i=0;i<=NumberOfComponents;i++)
          {
            for(j=CurrentBlocklength-1;j>0;j--)
              BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][j]=BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][j-1];
            BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][0]=value_onsager[i];
          }

          // vacf for each component
          for(k=0;k<CurrentBlocklength;k++)
          {
            for(i=0;i<NumberOfComponents;i++)
            {
              VacfOrderNOnsagerCount[CurrentSystem][CurrentBlock][i][k]+=1.0;
              for(j=0;j<NumberOfComponents;j++)
              {
                VacfOrderNOnsager[CurrentSystem][CurrentBlock][i][j][k].x+=(BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].x*value_onsager[j].x);
                VacfOrderNOnsager[CurrentSystem][CurrentBlock][i][j][k].y+=(BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].y*value_onsager[j].y);
                VacfOrderNOnsager[CurrentSystem][CurrentBlock][i][j][k].z+=(BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].z*value_onsager[j].z);
                VacfOrderNOnsagerDirAvg[CurrentSystem][CurrentBlock][i][j][k]+=(BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].x*value_onsager[j].x)+
                                                                               (BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].y*value_onsager[j].y)+
                                                                               (BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].z*value_onsager[j].z);
              }
            }
          }
          // vacf for the total fluid
          for(k=0;k<CurrentBlocklength;k++)
          {
            VacfOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][k]+=1.0;
            VacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][k].x+=
              (BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].x*value_onsager[NumberOfComponents].x);
            VacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][k].y+=
              (BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].x*value_onsager[NumberOfComponents].y);
            VacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][k].z+=
              (BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].x*value_onsager[NumberOfComponents].z);
            VacfOrderNTotalOnsagerDirAvg[CurrentSystem][CurrentBlock][k]+=
                  (BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].x*value_onsager[NumberOfComponents].x)+
                  (BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].y*value_onsager[NumberOfComponents].y)+
                  (BlockDataVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].z*value_onsager[NumberOfComponents].z);
          }
        }
      }
      // CountVACFOrderN the current sampling
      CountVACFOrderN[CurrentSystem]++;
      break;
    case PRINT:
      // return if the vacf does not has to be calculated for this system
      if((!ComputeVACFOrderN[CurrentSystem])||(CurrentCycle%WriteVACFOrderNEvery[CurrentSystem]!=0)) return;

      mkdir("VACFOrderN",S_IRWXU);

      sprintf(buffer,"VACFOrderN/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      // Self diffusion
      // ========================================================================================================================

      sprintf(buffer,"VACFOrderN/System_%d/vacf_self_total%s.dat",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"# column 1: time [ps]\n");
      fprintf(FilePtr,"# column 2: vacf xyz\n");
      fprintf(FilePtr,"# column 3: vacf x\n");
      fprintf(FilePtr,"# column 4: vacf y\n");
      fprintf(FilePtr,"# column 5: vacf z\n");
      fprintf(FilePtr,"# column 6: number of samples [-]\n");

      for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksVACFOrderN,NumberOfBlocksVACFOrderN[CurrentSystem]);CurrentBlock++)
      {
        CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);
        dt=SampleVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock);
        for(j=1;j<CurrentBlocklength;j++)
        {
          count=0.0;
          value_dir_avg=0.0;
          value.x=value.y=value.z=0.0;
          for(i=0;i<NumberOfComponents;i++)
          {
            value.x+=VacfOrderN[CurrentSystem][CurrentBlock][i][j].x;
            value.y+=VacfOrderN[CurrentSystem][CurrentBlock][i][j].y;
            value.z+=VacfOrderN[CurrentSystem][CurrentBlock][i][j].z;
            value_dir_avg+=VacfOrderNDirAvg[CurrentSystem][CurrentBlock][i][j];
            count+=VacfOrderNCount[CurrentSystem][CurrentBlock][i][j];
          }

          if(count>0.0)
          {
            fac=1.0/count;

            fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                (double)(j*dt),
                (double)(fac*value_dir_avg),
                (double)(fac*value.x),
                (double)(fac*value.y),
                (double)(fac*value.z),
                (double)count);
          }
        }
      }

      fclose(FilePtr);


      // print out vacf per component
      for(i=0;i<NumberOfComponents;i++)
      {
        sprintf(buffer,"VACFOrderN/System_%d/vacf_self_%s_%d%s.dat",
                CurrentSystem,Components[i].Name,i,FileNameAppend);

        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# column 1: time [ps]\n");
        fprintf(FilePtr,"# column 2: vacf xyz\n");
        fprintf(FilePtr,"# column 3: vacf x\n");
        fprintf(FilePtr,"# column 4: vacf y\n");
        fprintf(FilePtr,"# column 5: vacf z\n");
        fprintf(FilePtr,"# column 6: number of samples [-]\n");

        // write vacf results averaged per component to file
        D=0.0;
        for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksVACFOrderN,NumberOfBlocksVACFOrderN[CurrentSystem]);CurrentBlock++)
        {
          CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);
          dt=SampleVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock);
          if(CurrentBlocklength>0)
          {
            fac=DIFFUSION_CONVERSION_FACTOR/3.0;
            Integral=fac*IntegrateGeneralizedSimpsonOrderN(&VacfOrderNDirAvg[CurrentSystem][CurrentBlock][i][1],
                        &VacfOrderNCount[CurrentSystem][CurrentBlock][i][1],dt,CurrentBlocklength-1);
            D+=Integral;
            fprintf(FilePtr,"# Block: %d D %g [m^2/s] contribution of this block: %g [m^2/s]\n",CurrentBlock,D,Integral);
          }
        }

        for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksVACFOrderN,NumberOfBlocksVACFOrderN[CurrentSystem]);CurrentBlock++)
        {
          CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);
          dt=SampleVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock);
          for(j=1;j<CurrentBlocklength;j++)
          {
            if(VacfOrderNCount[CurrentSystem][CurrentBlock][i][j]>0.0)
            {
              fac=1.0/VacfOrderNCount[CurrentSystem][CurrentBlock][i][j];

              fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                  (double)(j*dt),
                  (double)(fac*VacfOrderNDirAvg[CurrentSystem][CurrentBlock][i][j]),
                  (double)(fac*VacfOrderN[CurrentSystem][CurrentBlock][i][j].x),
                  (double)(fac*VacfOrderN[CurrentSystem][CurrentBlock][i][j].y),
                  (double)(fac*VacfOrderN[CurrentSystem][CurrentBlock][i][j].z),
                  (double)VacfOrderNCount[CurrentSystem][CurrentBlock][i][j]);
            }
          }
        }

        fclose(FilePtr);
      }

      if(ComputeIndividualVACFOrderN)
      {
        sprintf(buffer,"VACFOrderN/System_%d/VACFOrderN_per_adsorbate",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        // print out vacf per adsorbate
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          sprintf(buffer,"VACFOrderN/System_%d/VACFOrderN_per_adsorbate/vacf_self_%d_%s.dat",
                  CurrentSystem,i,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: vacf xyz\n");
          fprintf(FilePtr,"# column 3: vacf x\n");
          fprintf(FilePtr,"# column 4: vacf y\n");
          fprintf(FilePtr,"# column 5: vacf z\n");
          fprintf(FilePtr,"# column 6: number of samples [-]\n");

          D=0.0;
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksVACFOrderN,NumberOfBlocksVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);
            dt=SampleVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock);
            fac=DIFFUSION_CONVERSION_FACTOR/3.0;
            Integral=fac*IntegrateGeneralizedSimpsonOrderN(&VacfOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][i][1],
                        &VacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][1],dt,CurrentBlocklength-1);
            D+=Integral;
            fprintf(FilePtr,"# Block: %d D %g [m^2/s] contribution of this block: %g [m^2/s]\n",CurrentBlock,D,Integral);
          }


          // write results to file
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksVACFOrderN,NumberOfBlocksVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);
            dt=SampleVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock);
            for(j=1;j<CurrentBlocklength;j++)
            {
              if(VacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j]>0.0)
              {
                fac=1.0/VacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j];
                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(j*dt),
                    (double)(fac*VacfOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][i][j]),
                    (double)(fac*VacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].x),
                    (double)(fac*VacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].y),
                    (double)(fac*VacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].z),
                    (double)VacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j]);
              }
            }
          }

          fclose(FilePtr);
        }

        sprintf(buffer,"VACFOrderN/System_%d/VACFOrderN_per_cation",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        // print out vacf per cation
        shift=NumberOfAdsorbateMolecules[CurrentSystem];
        for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
        {
          index=i+shift;
          sprintf(buffer,"VACFOrderN/System_%d/VACFOrderN_per_cation/vacf_self_%d%s.dat",
                  CurrentSystem,i,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: vacf xyz\n");
          fprintf(FilePtr,"# column 3: vacf x\n");
          fprintf(FilePtr,"# column 4: vacf y\n");
          fprintf(FilePtr,"# column 5: vacf z\n");
          fprintf(FilePtr,"# column 6: number of samples [-]\n");

          // write results to file

          D=0.0;
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksVACFOrderN,NumberOfBlocksVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);
            dt=SampleVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock);
            fac=DIFFUSION_CONVERSION_FACTOR/3.0;
            Integral=fac*IntegrateGeneralizedSimpsonOrderN(&VacfOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][index][1],
                        &VacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][1],dt,CurrentBlocklength-1);
            D+=Integral;
            fprintf(FilePtr,"# Block: %d D %g [m^2/s] contribution of this block: %g [m^2/s]\n",CurrentBlock,D,Integral);
          }

          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksVACFOrderN,NumberOfBlocksVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);
            dt=SampleVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock);
            for(j=1;j<CurrentBlocklength;j++)
            {
              if(VacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][j]>0.0)
              {
                fac=1.0/VacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][j];
                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(j*dt),
                    (double)(fac*VacfOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][index][j]),
                    (double)(fac*VacfOrderNPerMolecule[CurrentSystem][CurrentBlock][index][j].x),
                    (double)(fac*VacfOrderNPerMolecule[CurrentSystem][CurrentBlock][index][j].y),
                    (double)(fac*VacfOrderNPerMolecule[CurrentSystem][CurrentBlock][index][j].z),
                    (double)VacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][j]);
              }
            }
          }

          fclose(FilePtr);
        }

      }

      // Onsager diffusion
      // ========================================================================================================================

      // write vacf results for the fluid to file
      sprintf(buffer,"VACFOrderN/System_%d/vacf_onsager_total%s.dat",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"# column 1: time [ps]\n");
      fprintf(FilePtr,"# column 2: vacf xyz\n");
      fprintf(FilePtr,"# column 3: vacf x\n");
      fprintf(FilePtr,"# column 4: vacf y\n");
      fprintf(FilePtr,"# column 5: vacf z\n");
      fprintf(FilePtr,"# column 6: number of samples [-]\n");

      for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksVACFOrderN,NumberOfBlocksVACFOrderN[CurrentSystem]);CurrentBlock++)
      {
        CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);
        dt=SampleVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock);
        for(l=1;l<CurrentBlocklength;l++)
        {
          if(VacfOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][l]>0.0)
          {
             fac=1.0/((NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem])*VacfOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][l]);
             fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                (double)(l*dt),
                (double)(fac*VacfOrderNTotalOnsagerDirAvg[CurrentSystem][CurrentBlock][l]),
                (double)(fac*VacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][l].x),
                (double)(fac*VacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][l].y),
                (double)(fac*VacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][l].z),
                (double)VacfOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][l]);
          }
        }
      }
      fclose(FilePtr);


      // write the vacf per component
      for(j=0;j<NumberOfComponents;j++)
      {
        for(k=0;k<NumberOfComponents;k++)
        {
          sprintf(buffer,"VACFOrderN/System_%d/vacf_onsager_%s_%d_%s_%d%s.dat",
                  CurrentSystem,Components[j].Name,j,Components[k].Name,k,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: vacf xyz\n");
          fprintf(FilePtr,"# column 3: vacf x\n");
          fprintf(FilePtr,"# column 4: vacf y\n");
          fprintf(FilePtr,"# column 5: vacf z\n");
          fprintf(FilePtr,"# column 6: number of samples [-]\n");

          // write vacf results averaged per component to file
          D=0.0;
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksVACFOrderN,NumberOfBlocksVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);
            dt=SampleVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock);
            fac=DIFFUSION_CONVERSION_FACTOR/(3.0*Components[j].NumberOfMolecules[CurrentSystem]);
            Integral=fac*IntegrateGeneralizedSimpsonOrderN(&VacfOrderNOnsagerDirAvg[CurrentSystem][CurrentBlock][j][k][1],
                            &VacfOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][1],dt,CurrentBlocklength-1);
            D+=Integral;
            fprintf(FilePtr,"# Block: %d D %g [m^2/s] contribution of this block: %g [m^2/s]\n",CurrentBlock,D,Integral);
          }

          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksVACFOrderN,NumberOfBlocksVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsVACFOrderN);
            dt=SampleVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsVACFOrderN,CurrentBlock);
            for(l=1;l<CurrentBlocklength;l++)
            {
              if(VacfOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][l]>0.0)
              {
                fac=1.0/(Components[j].NumberOfMolecules[CurrentSystem]*VacfOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][l]);
                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(l*dt),
                    (double)(fac*VacfOrderNOnsagerDirAvg[CurrentSystem][CurrentBlock][j][k][l]),
                    (double)(fac*VacfOrderNOnsager[CurrentSystem][CurrentBlock][j][k][l].x),
                    (double)(fac*VacfOrderNOnsager[CurrentSystem][CurrentBlock][j][k][l].y),
                    (double)(fac*VacfOrderNOnsager[CurrentSystem][CurrentBlock][j][k][l].z),
                    (double)VacfOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][l]);
              }
            }
          }

          fclose(FilePtr);
        }
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeVACFOrderN[i])
        {

          for(j=0;j<MaxNumberOfBlocksVACFOrderN;j++)
          {
            for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              free(BlockDataVACFOrderN[i][j][k]);
            free(BlockDataVACFOrderN[i][j]);

            for(k=0;k<NumberOfComponents;k++)
            {
              free(VacfOrderNCount[i][j][k]);
              free(VacfOrderNDirAvg[i][j][k]);
              free(VacfOrderN[i][j][k]);
            }
            free(VacfOrderNCount[i][j]);
            free(VacfOrderNDirAvg[i][j]);
            free(VacfOrderN[i][j]);

            if(ComputeIndividualVACFOrderN)
            {
              for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              {
                free(VacfOrderNCountPerMolecule[i][j][k]);
                free(VacfOrderNPerMolecule[i][j][k]);
                free(VacfOrderNPerMoleculeDirAvg[i][j][k]);
              }
              free(VacfOrderNCountPerMolecule[i][j]);
              free(VacfOrderNPerMoleculeDirAvg[i][j]);
              free(VacfOrderNPerMolecule[i][j]);
            }

            // Onsager data
            for(k=0;k<NumberOfComponents;k++)
            {
              for(l=0;l<NumberOfComponents;l++)
              {
                free(VacfOrderNOnsager[i][j][k][l]);
                free(VacfOrderNOnsagerDirAvg[i][j][k][l]);
              }
              free(VacfOrderNOnsagerCount[i][j][k]);
              free(VacfOrderNOnsagerDirAvg[i][j][k]);
              free(VacfOrderNOnsager[i][j][k]);
            }

            for(k=0;k<=NumberOfComponents;k++)
              free(BlockDataVACFOrderNOnsager[i][j][k]);
            free(BlockDataVACFOrderNOnsager[i][j]);
            free(VacfOrderNTotalOnsager[i][j]);
            free(VacfOrderNTotalOnsagerDirAvg[i][j]);
            free(VacfOrderNTotalOnsagerCount[i][j]);
            free(VacfOrderNOnsagerCount[i][j]);
            free(VacfOrderNOnsagerDirAvg[i][j]);
            free(VacfOrderNOnsager[i][j]);
          }

          if(ComputeIndividualVACFOrderN)
          {
            free(VacfOrderNCountPerMolecule[i]);
            free(VacfOrderNPerMoleculeDirAvg[i]);
            free(VacfOrderNPerMolecule[i]);
          }
          free(BlockLengthVACFOrderN[i]);
          free(VacfOrderNCount[i]);
          free(VacfOrderNOnsagerCount[i]);
          free(VacfOrderNDirAvg[i]);
          free(VacfOrderNOnsagerDirAvg[i]);
          free(VacfOrderN[i]);
          free(VacfOrderNOnsager[i]);
          free(BlockDataVACFOrderN[i]);
          free(BlockDataVACFOrderNOnsager[i]);
          free(VacfOrderNTotalOnsager[i]);
          free(VacfOrderNTotalOnsagerDirAvg[i]);
          free(VacfOrderNTotalOnsagerCount[i]);

        }
      }
      free(CountVACFOrderN);
      free(NumberOfBlocksVACFOrderN);
      free(value_onsager);
      free(BlockLengthVACFOrderN);
      free(VacfOrderNCount);
      free(VacfOrderNOnsagerCount);
      free(VacfOrderNCountPerMolecule);
      free(VacfOrderNDirAvg);
      free(VacfOrderNOnsagerDirAvg);
      free(VacfOrderNPerMoleculeDirAvg);
      free(VacfOrderN);
      free(VacfOrderNOnsager);
      free(VacfOrderNPerMolecule);
      free(BlockDataVACFOrderN);
      free(BlockDataVACFOrderNOnsager);
      free(VacfOrderNTotalOnsager);
      free(VacfOrderNTotalOnsagerDirAvg);
      free(VacfOrderNTotalOnsagerCount);
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleRotationalVelocityAutoCorrelationFunctionOrderN                                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the rotational velocity autocorrelation function using a modified order-N alg.   *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleRotationalVelocityAutoCorrelationFunctionOrderN(int Switch)
{
  int i,j,k,l,CurrentBlock,index;
  int CurrentBlocklength,type,shift;
  VECTOR value,drift;
  static VECTOR *value_onsager;
  FILE *FilePtr;
  char buffer[256];
  REAL fac,Integral,D,dt,count,value_dir_avg;

  switch(Switch)
  {
    case ALLOCATE:
      CountRVACFOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
      NumberOfBlocksRVACFOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
      value_onsager=(VECTOR*)calloc(NumberOfComponents+1,sizeof(VECTOR));

      BlockLengthRVACFOrderN=(int**)calloc(NumberOfSystems,sizeof(int*));
      RvacfOrderNCount=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      RvacfOrderNOnsagerCount=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      RvacfOrderNCountPerMolecule=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      RvacfOrderNDirAvg=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      RvacfOrderNOnsagerDirAvg=(REAL*****)calloc(NumberOfSystems,sizeof(REAL****));
      RvacfOrderNPerMoleculeDirAvg=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      RvacfOrderN=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      RvacfOrderNOnsager=(VECTOR*****)calloc(NumberOfSystems,sizeof(VECTOR****));
      RvacfOrderNPerMolecule=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      BlockDataRVACFOrderN=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      BlockDataRVACFOrderNOnsager=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      RvacfOrderNTotalOnsager=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      RvacfOrderNTotalOnsagerDirAvg=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      RvacfOrderNTotalOnsagerCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeRVACFOrderN[i])
        {
          BlockLengthRVACFOrderN[i]=(int*)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(int));
          RvacfOrderNCount[i]=(REAL***)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(REAL**));
          RvacfOrderNOnsagerCount[i]=(REAL***)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(REAL**));
          RvacfOrderNDirAvg[i]=(REAL***)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(REAL**));
          RvacfOrderNOnsagerDirAvg[i]=(REAL****)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(REAL***));
          RvacfOrderN[i]=(VECTOR***)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(VECTOR**));
          RvacfOrderNOnsager[i]=(VECTOR****)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(VECTOR***));
          BlockDataRVACFOrderN[i]=(VECTOR***)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(VECTOR**));
          BlockDataRVACFOrderNOnsager[i]=(VECTOR***)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(VECTOR**));
          RvacfOrderNTotalOnsager[i]=(VECTOR**)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(VECTOR*));
          RvacfOrderNTotalOnsagerDirAvg[i]=(REAL**)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(REAL*));
          RvacfOrderNTotalOnsagerCount[i]=(REAL**)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(REAL*));

          if(ComputeIndividualRVACFOrderN)
          {
            RvacfOrderNCountPerMolecule[i]=(REAL***)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(REAL**));
            RvacfOrderNPerMoleculeDirAvg[i]=(REAL***)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(REAL**));
            RvacfOrderNPerMolecule[i]=(VECTOR***)calloc(MaxNumberOfBlocksRVACFOrderN,sizeof(VECTOR**));
          }

          for(j=0;j<MaxNumberOfBlocksRVACFOrderN;j++)
          {
            BlockDataRVACFOrderN[i][j]=(VECTOR**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR*));
            for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              BlockDataRVACFOrderN[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR));

            RvacfOrderNCount[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            RvacfOrderNDirAvg[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            RvacfOrderN[i][j]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
            for(k=0;k<NumberOfComponents;k++)
            {
              RvacfOrderNCount[i][j][k]=(REAL*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(REAL));
              RvacfOrderNDirAvg[i][j][k]=(REAL*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(REAL));
              RvacfOrderN[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR));
            }

            if(ComputeIndividualRVACFOrderN)
            {
              RvacfOrderNCountPerMolecule[i][j]=(REAL**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(REAL*));
              RvacfOrderNPerMoleculeDirAvg[i][j]=(REAL**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(REAL*));
              RvacfOrderNPerMolecule[i][j]=(VECTOR**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR*));
              for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              {
                RvacfOrderNCountPerMolecule[i][j][k]=(REAL*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(REAL));
                RvacfOrderNPerMolecule[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR));
                RvacfOrderNPerMoleculeDirAvg[i][j][k]=(REAL*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(REAL));
              }
            }

            // Onsager data
            BlockDataRVACFOrderNOnsager[i][j]=(VECTOR**)calloc(NumberOfComponents+1,sizeof(VECTOR*));
            RvacfOrderNTotalOnsager[i][j]=(VECTOR*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR));
            RvacfOrderNTotalOnsagerDirAvg[i][j]=(REAL*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(REAL));
            RvacfOrderNTotalOnsagerCount[i][j]=(REAL*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(REAL));
            RvacfOrderNOnsagerCount[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            RvacfOrderNOnsagerDirAvg[i][j]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));
            RvacfOrderNOnsager[i][j]=(VECTOR***)calloc(NumberOfComponents,sizeof(VECTOR**));

            for(k=0;k<=NumberOfComponents;k++)
              BlockDataRVACFOrderNOnsager[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR));

            for(k=0;k<NumberOfComponents;k++)
            {
              RvacfOrderNOnsagerCount[i][j][k]=(REAL*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(REAL));
              RvacfOrderNOnsagerDirAvg[i][j][k]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
              RvacfOrderNOnsager[i][j][k]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
              for(l=0;l<NumberOfComponents;l++)
              {
                RvacfOrderNOnsager[i][j][k][l]=(VECTOR*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR));
                RvacfOrderNOnsagerDirAvg[i][j][k][l]=(REAL*)calloc(NumberOfBlockElementsRVACFOrderN,sizeof(REAL));
              }
            }
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      // return if the vacf does not has to be calculated for this system
      if(!(ComputeRVACFOrderN[CurrentSystem]&&(CurrentCycle%SampleRVACFOrderNEvery[CurrentSystem]==0))) return;

      // compute drift of the system
      if(Framework[CurrentSystem].FrameworkModel==NONE)
        drift=GetCenterOfMassVelocityCurrentSystem();
      else if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
        drift=GetFrameworkCenterOfMassVelocity();
      else
      {
        drift.x=0.0;
        drift.y=0.0;
        drift.z=0.0;
      }

      // determine current number of blocks
      NumberOfBlocksRVACFOrderN[CurrentSystem]=1;
      i=CountRVACFOrderN[CurrentSystem]/NumberOfBlockElementsRVACFOrderN;
      while(i!=0)
      {
        NumberOfBlocksRVACFOrderN[CurrentSystem]++;
        i/=NumberOfBlockElementsRVACFOrderN;
      }

      // ignore everything beyond the last block
      NumberOfBlocksRVACFOrderN[CurrentSystem]=MIN2(NumberOfBlocksRVACFOrderN[CurrentSystem],MaxNumberOfBlocksRVACFOrderN);

      for(CurrentBlock=0;CurrentBlock<NumberOfBlocksRVACFOrderN[CurrentSystem];CurrentBlock++)
      {
        // test for blocking operation: CountRVACFOrderN is a multiple of NumberOfBlockElementsRVACFOrderN^CurrentBlock
        if ((CountRVACFOrderN[CurrentSystem])%((int)pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock))==0)
        {
          // increase the current block-length
          BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock]++;

          // limit length to NumberOfBlockElementsRVACFOrderN
          CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);

          // Self diffusion
          // ========================================================================================================================

          shift=NumberOfAdsorbateMolecules[CurrentSystem];
          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];i++)
          {
            if(i<NumberOfAdsorbateMolecules[CurrentSystem])
              value=QuaternionMomentumToAngularVelocityAdsorbates(i,0);
            else
              value=QuaternionMomentumToAngularVelocityCations(i-shift,0);

            for(j=CurrentBlocklength-1;j>0;j--)
              BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j]=BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j-1];
            BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][0]=value;


            // get the type of the molecule
            if(i<NumberOfAdsorbateMolecules[CurrentSystem])
              type=Adsorbates[CurrentSystem][i].Type;
            else
              type=Cations[CurrentSystem][i-shift].Type;

            for(j=0;j<CurrentBlocklength;j++)
            {
              // vacf for each component
              RvacfOrderNCount[CurrentSystem][CurrentBlock][type][j]+=1.0;
              RvacfOrderN[CurrentSystem][CurrentBlock][type][j].x+=(BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].x*value.x);
              RvacfOrderN[CurrentSystem][CurrentBlock][type][j].y+=(BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].y*value.y);
              RvacfOrderN[CurrentSystem][CurrentBlock][type][j].z+=(BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].z*value.z);
              RvacfOrderNDirAvg[CurrentSystem][CurrentBlock][type][j]+=(BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].x*value.x)+
                                                                      (BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].y*value.y)+
                                                                      (BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].z*value.z);

              if(ComputeIndividualRVACFOrderN)
              {
                RvacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j]+=1.0;
                RvacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].x+=(BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].x*value.x);
                RvacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].y+=(BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].y*value.y);
                RvacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].z+=(BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].z*value.z);
                RvacfOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][i][j]+=(BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].x*value.x)+
                                                                                (BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].y*value.y)+
                                                                                (BlockDataRVACFOrderN[CurrentSystem][CurrentBlock][i][j].z*value.z);
              }
            }
          }

          // Onsager diffusion
          // ========================================================================================================================

          for(i=0;i<=NumberOfComponents;i++)
            value_onsager[i].x=value_onsager[i].y=value_onsager[i].z=0.0;

          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
          {
            type=Adsorbates[CurrentSystem][i].Type;

           if(i<NumberOfAdsorbateMolecules[CurrentSystem])
              value=QuaternionMomentumToAngularVelocityAdsorbates(i,0);
            else
              value=QuaternionMomentumToAngularVelocityCations(i-shift,0);

            value_onsager[type].x+=value.x;
            value_onsager[type].y+=value.y;
            value_onsager[type].z+=value.z;

            value_onsager[NumberOfComponents].x+=value.x;
            value_onsager[NumberOfComponents].y+=value.y;
            value_onsager[NumberOfComponents].z+=value.z;
          }

          for(i=0;i<=NumberOfComponents;i++)
          {
            for(j=CurrentBlocklength-1;j>0;j--)
              BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][j]=BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][j-1];
            BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][0]=value_onsager[i];
          }

          // vacf for each component
          for(k=0;k<CurrentBlocklength;k++)
          {
            for(i=0;i<NumberOfComponents;i++)
            {
              RvacfOrderNOnsagerCount[CurrentSystem][CurrentBlock][i][k]+=1.0;
              for(j=0;j<NumberOfComponents;j++)
              {
                RvacfOrderNOnsager[CurrentSystem][CurrentBlock][i][j][k].x+=(BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].x*value_onsager[j].x);
                RvacfOrderNOnsager[CurrentSystem][CurrentBlock][i][j][k].y+=(BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].y*value_onsager[j].y);
                RvacfOrderNOnsager[CurrentSystem][CurrentBlock][i][j][k].z+=(BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].z*value_onsager[j].z);
                RvacfOrderNOnsagerDirAvg[CurrentSystem][CurrentBlock][i][j][k]+=(BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].x*value_onsager[j].x)+
                                                                               (BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].y*value_onsager[j].y)+
                                                                               (BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][i][k].z*value_onsager[j].z);
              }
            }
          }
          // vacf for the total fluid
          for(k=0;k<CurrentBlocklength;k++)
          {
            RvacfOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][k]+=1.0;
            RvacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][k].x+=
              (BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].x*value_onsager[NumberOfComponents].x);
            RvacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][k].y+=
              (BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].x*value_onsager[NumberOfComponents].y);
            RvacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][k].z+=
              (BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].x*value_onsager[NumberOfComponents].z);
            RvacfOrderNTotalOnsagerDirAvg[CurrentSystem][CurrentBlock][k]+=
                  (BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].x*value_onsager[NumberOfComponents].x)+
                  (BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].y*value_onsager[NumberOfComponents].y)+
                  (BlockDataRVACFOrderNOnsager[CurrentSystem][CurrentBlock][NumberOfComponents][k].z*value_onsager[NumberOfComponents].z);
          }
        }
      }
      // CountRVACFOrderN the current sampling
      CountRVACFOrderN[CurrentSystem]++;
      break;
    case PRINT:
      // return if the vacf does not has to be calculated for this system
      if((!ComputeRVACFOrderN[CurrentSystem])||(CurrentCycle%WriteRVACFOrderNEvery[CurrentSystem]!=0)) return;

      mkdir("RVACFOrderN",S_IRWXU);

      sprintf(buffer,"RVACFOrderN/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      // Self diffusion
      // ========================================================================================================================

      sprintf(buffer,"RVACFOrderN/System_%d/vacf_self_total%s.dat",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"# column 1: time [ps]\n");
      fprintf(FilePtr,"# column 2: rvacf xyz\n");
      fprintf(FilePtr,"# column 3: rvacf x\n");
      fprintf(FilePtr,"# column 4: rvacf y\n");
      fprintf(FilePtr,"# column 5: rvacf z\n");
      fprintf(FilePtr,"# column 6: number of samples [-]\n");

      for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksRVACFOrderN,NumberOfBlocksRVACFOrderN[CurrentSystem]);CurrentBlock++)
      {
        CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);
        dt=SampleRVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock);
        for(j=1;j<CurrentBlocklength;j++)
        {
          count=0.0;
          value_dir_avg=0.0;
          value.x=value.y=value.z=0.0;
          for(i=0;i<NumberOfComponents;i++)
          {
            value.x+=RvacfOrderN[CurrentSystem][CurrentBlock][i][j].x;
            value.y+=RvacfOrderN[CurrentSystem][CurrentBlock][i][j].y;
            value.z+=RvacfOrderN[CurrentSystem][CurrentBlock][i][j].z;
            value_dir_avg+=RvacfOrderNDirAvg[CurrentSystem][CurrentBlock][i][j];
            count+=RvacfOrderNCount[CurrentSystem][CurrentBlock][i][j];
          }

          if(count>0.0)
          {
            fac=1.0/count;

            fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                (double)(j*dt),
                (double)(fac*value_dir_avg),
                (double)(fac*value.x),
                (double)(fac*value.y),
                (double)(fac*value.z),
                (double)count);
          }
        }
      }

      fclose(FilePtr);


      // print out vacf per component
      for(i=0;i<NumberOfComponents;i++)
      {
        sprintf(buffer,"RVACFOrderN/System_%d/vacf_self_%s_%d%s.dat",
                CurrentSystem,Components[i].Name,i,FileNameAppend);

        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# column 1: time [ps]\n");
        fprintf(FilePtr,"# column 2: rvacf xyz\n");
        fprintf(FilePtr,"# column 3: rvacf x\n");
        fprintf(FilePtr,"# column 4: rvacf y\n");
        fprintf(FilePtr,"# column 5: rvacf z\n");
        fprintf(FilePtr,"# column 6: number of samples [-]\n");


        // write vacf results averaged per component to file
        D=0.0;
        for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksRVACFOrderN,NumberOfBlocksRVACFOrderN[CurrentSystem]);CurrentBlock++)
        {
          CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);
          dt=SampleRVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock);
          if(CurrentBlocklength>0)
          {
            fac=DIFFUSION_CONVERSION_FACTOR/3.0;
            Integral=fac*IntegrateGeneralizedSimpsonOrderN(&RvacfOrderNDirAvg[CurrentSystem][CurrentBlock][i][1],
                        &RvacfOrderNCount[CurrentSystem][CurrentBlock][i][1],dt,CurrentBlocklength-1);
            D+=Integral;
            fprintf(FilePtr,"# Block: %d D %g [m^2/s] contribution of this block: %g [m^2/s]\n",CurrentBlock,D,Integral);
          }
        }

        for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksRVACFOrderN,NumberOfBlocksRVACFOrderN[CurrentSystem]);CurrentBlock++)
        {
          CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);
          dt=SampleRVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock);
          for(j=1;j<CurrentBlocklength;j++)
          {
            if(RvacfOrderNCount[CurrentSystem][CurrentBlock][i][j]>0.0)
            {
              fac=1.0/RvacfOrderNCount[CurrentSystem][CurrentBlock][i][j];

              fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                  (double)(j*dt),
                  (double)(fac*RvacfOrderNDirAvg[CurrentSystem][CurrentBlock][i][j]),
                  (double)(fac*RvacfOrderN[CurrentSystem][CurrentBlock][i][j].x),
                  (double)(fac*RvacfOrderN[CurrentSystem][CurrentBlock][i][j].y),
                  (double)(fac*RvacfOrderN[CurrentSystem][CurrentBlock][i][j].z),
                  (double)RvacfOrderNCount[CurrentSystem][CurrentBlock][i][j]);
            }
          }
        }

        fclose(FilePtr);
      }

      if(ComputeIndividualRVACFOrderN)
      {
        sprintf(buffer,"RVACFOrderN/System_%d/RVACFOrderN_per_adsorbate",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        // print out vacf per adsorbate
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          sprintf(buffer,"RVACFOrderN/System_%d/RVACFOrderN_per_adsorbate/vacf_self_%d%s.dat",
                  CurrentSystem,i,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: rvacf xyz\n");
          fprintf(FilePtr,"# column 3: rvacf x\n");
          fprintf(FilePtr,"# column 4: rvacf y\n");
          fprintf(FilePtr,"# column 5: rvacf z\n");
          fprintf(FilePtr,"# column 6: number of samples [-]\n");

          D=0.0;
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksRVACFOrderN,NumberOfBlocksRVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);
            dt=SampleRVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock);
            fac=DIFFUSION_CONVERSION_FACTOR/3.0;
            Integral=fac*IntegrateGeneralizedSimpsonOrderN(&RvacfOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][i][1],
                        &RvacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][1],dt,CurrentBlocklength-1);
            D+=Integral;
            fprintf(FilePtr,"# Block: %d D %g [m^2/s] contribution of this block: %g [m^2/s]\n",CurrentBlock,D,Integral);
          }


          // write results to file
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksRVACFOrderN,NumberOfBlocksRVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);
            dt=SampleRVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock);
            for(j=1;j<CurrentBlocklength;j++)
            {
              if(RvacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j]>0.0)
              {
                fac=1.0/RvacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j];
                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(j*dt),
                    (double)(fac*RvacfOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][i][j]),
                    (double)(fac*RvacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].x),
                    (double)(fac*RvacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].y),
                    (double)(fac*RvacfOrderNPerMolecule[CurrentSystem][CurrentBlock][i][j].z),
                    (double)RvacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][i][j]);
              }
            }
          }

          fclose(FilePtr);
        }

        sprintf(buffer,"RVACFOrderN/System_%d/RVACFOrderN_per_cation",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        // print out vacf per cation
        shift=NumberOfAdsorbateMolecules[CurrentSystem];
        for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
        {
          index=i+shift;
          sprintf(buffer,"RVACFOrderN/System_%d/RVACFOrderN_per_cation/vacf_self_%d%s.dat",
                  CurrentSystem,i,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: rvacf xyz\n");
          fprintf(FilePtr,"# column 3: rvacf x\n");
          fprintf(FilePtr,"# column 4: rvacf y\n");
          fprintf(FilePtr,"# column 5: rvacf z\n");
          fprintf(FilePtr,"# column 6: number of samples [-]\n");

          // write results to file

          D=0.0;
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksRVACFOrderN,NumberOfBlocksRVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);
            dt=SampleRVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock);
            fac=DIFFUSION_CONVERSION_FACTOR/3.0;
            Integral=fac*IntegrateGeneralizedSimpsonOrderN(&RvacfOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][index][1],
                        &RvacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][1],dt,CurrentBlocklength-1);
            D+=Integral;
            fprintf(FilePtr,"# Block: %d D %g [m^2/s] contribution of this block: %g [m^2/s]\n",CurrentBlock,D,Integral);
          }

          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksRVACFOrderN,NumberOfBlocksRVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);
            dt=SampleRVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock);
            for(j=1;j<CurrentBlocklength;j++)
            {
              if(RvacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][j]>0.0)
              {
                fac=1.0/RvacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][j];
                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(j*dt),
                    (double)(fac*RvacfOrderNPerMoleculeDirAvg[CurrentSystem][CurrentBlock][index][j]),
                    (double)(fac*RvacfOrderNPerMolecule[CurrentSystem][CurrentBlock][index][j].x),
                    (double)(fac*RvacfOrderNPerMolecule[CurrentSystem][CurrentBlock][index][j].y),
                    (double)(fac*RvacfOrderNPerMolecule[CurrentSystem][CurrentBlock][index][j].z),
                    (double)RvacfOrderNCountPerMolecule[CurrentSystem][CurrentBlock][index][j]);
              }
            }
          }

          fclose(FilePtr);
        }

      }

      // Onsager diffusion
      // ========================================================================================================================

      // write vacf results for the fluid to file
      sprintf(buffer,"RVACFOrderN/System_%d/vacf_onsager_total%s.dat",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"# column 1: time [ps]\n");
      fprintf(FilePtr,"# column 2: rvacf xyz\n");
      fprintf(FilePtr,"# column 3: rvacf x\n");
      fprintf(FilePtr,"# column 4: rvacf y\n");
      fprintf(FilePtr,"# column 5: rvacf z\n");
      fprintf(FilePtr,"# column 6: number of samples [-]\n");

      for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksRVACFOrderN,NumberOfBlocksRVACFOrderN[CurrentSystem]);CurrentBlock++)
      {
        CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);
        dt=SampleRVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock);
        for(l=1;l<CurrentBlocklength;l++)
        {
          if(RvacfOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][l]>0.0)
          {
             fac=1.0/((NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem])*RvacfOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][l]);
             fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                (double)(l*dt),
                (double)(fac*RvacfOrderNTotalOnsagerDirAvg[CurrentSystem][CurrentBlock][l]),
                (double)(fac*RvacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][l].x),
                (double)(fac*RvacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][l].y),
                (double)(fac*RvacfOrderNTotalOnsager[CurrentSystem][CurrentBlock][l].z),
                (double)RvacfOrderNTotalOnsagerCount[CurrentSystem][CurrentBlock][l]);
          }
        }
      }
      fclose(FilePtr);


      // write the vacf per component
      for(j=0;j<NumberOfComponents;j++)
      {
        for(k=0;k<NumberOfComponents;k++)
        {
          sprintf(buffer,"RVACFOrderN/System_%d/vacf_onsager_%s_%d_%s_%d%s.dat",
                  CurrentSystem,Components[j].Name,j,Components[k].Name,k,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: rvacf xyz\n");
          fprintf(FilePtr,"# column 3: rvacf x\n");
          fprintf(FilePtr,"# column 4: rvacf y\n");
          fprintf(FilePtr,"# column 5: rvacf z\n");
          fprintf(FilePtr,"# column 6: number of samples [-]\n");

          // write vacf results averaged per component to file
          D=0.0;
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksRVACFOrderN,NumberOfBlocksRVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);
            dt=SampleRVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock);
            fac=DIFFUSION_CONVERSION_FACTOR/(3.0*Components[j].NumberOfMolecules[CurrentSystem]);
            Integral=fac*IntegrateGeneralizedSimpsonOrderN(&RvacfOrderNOnsagerDirAvg[CurrentSystem][CurrentBlock][j][k][1],
                            &RvacfOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][1],dt,CurrentBlocklength-1);
            D+=Integral;
            fprintf(FilePtr,"# Block: %d D %g [m^2/s] contribution of this block: %g [m^2/s]\n",CurrentBlock,D,Integral);
          }

          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksRVACFOrderN,NumberOfBlocksRVACFOrderN[CurrentSystem]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthRVACFOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsRVACFOrderN);
            dt=SampleRVACFOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsRVACFOrderN,CurrentBlock);
            for(l=1;l<CurrentBlocklength;l++)
            {
              if(RvacfOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][l]>0.0)
              {
                fac=1.0/(Components[j].NumberOfMolecules[CurrentSystem]*RvacfOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][l]);
                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(l*dt),
                    (double)(fac*RvacfOrderNOnsagerDirAvg[CurrentSystem][CurrentBlock][j][k][l]),
                    (double)(fac*RvacfOrderNOnsager[CurrentSystem][CurrentBlock][j][k][l].x),
                    (double)(fac*RvacfOrderNOnsager[CurrentSystem][CurrentBlock][j][k][l].y),
                    (double)(fac*RvacfOrderNOnsager[CurrentSystem][CurrentBlock][j][k][l].z),
                    (double)RvacfOrderNOnsagerCount[CurrentSystem][CurrentBlock][j][l]);
              }
            }
          }

          fclose(FilePtr);
        }
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeRVACFOrderN[i])
        {
          for(j=0;j<MaxNumberOfBlocksRVACFOrderN;j++)
          {
            for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              free(BlockDataRVACFOrderN[i][j][k]);
            free(BlockDataRVACFOrderN[i][j]);

            for(k=0;k<NumberOfComponents;k++)
            {
              free(RvacfOrderNCount[i][j][k]);
              free(RvacfOrderNDirAvg[i][j][k]);
              free(RvacfOrderN[i][j][k]);
            }
            free(RvacfOrderNCount[i][j]);
            free(RvacfOrderNDirAvg[i][j]);
            free(RvacfOrderN[i][j]);

            if(ComputeIndividualRVACFOrderN)
            {
              for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              {
                free(RvacfOrderNCountPerMolecule[i][j][k]);
                free(RvacfOrderNPerMolecule[i][j][k]);
                free(RvacfOrderNPerMoleculeDirAvg[i][j][k]);
              }
              free(RvacfOrderNCountPerMolecule[i][j]);
              free(RvacfOrderNPerMoleculeDirAvg[i][j]);
              free(RvacfOrderNPerMolecule[i][j]);
            }

            // Onsager data
            for(k=0;k<=NumberOfComponents;k++)
              free(BlockDataRVACFOrderNOnsager[i][j][k]);

            for(k=0;k<NumberOfComponents;k++)
            {
              for(l=0;l<NumberOfComponents;l++)
              {
                free(RvacfOrderNOnsager[i][j][k][l]);
                free(RvacfOrderNOnsagerDirAvg[i][j][k][l]);
              }
              free(RvacfOrderNOnsagerCount[i][j][k]);
              free(RvacfOrderNOnsagerDirAvg[i][j][k]);
              free(RvacfOrderNOnsager[i][j][k]);
            }

            free(BlockDataRVACFOrderNOnsager[i][j]);
            free(RvacfOrderNTotalOnsager[i][j]);
            free(RvacfOrderNTotalOnsagerDirAvg[i][j]);
            free(RvacfOrderNTotalOnsagerCount[i][j]);
            free(RvacfOrderNOnsagerCount[i][j]);
            free(RvacfOrderNOnsagerDirAvg[i][j]);
            free(RvacfOrderNOnsager[i][j]);

          }

          if(ComputeIndividualRVACFOrderN)
          {
            free(RvacfOrderNCountPerMolecule[i]);
            free(RvacfOrderNPerMoleculeDirAvg[i]);
            free(RvacfOrderNPerMolecule[i]);
          }
          free(BlockLengthRVACFOrderN[i]);
          free(RvacfOrderNCount[i]);
          free(RvacfOrderNOnsagerCount[i]);
          free(RvacfOrderNDirAvg[i]);
          free(RvacfOrderNOnsagerDirAvg[i]);
          free(RvacfOrderN[i]);
          free(RvacfOrderNOnsager[i]);
          free(BlockDataRVACFOrderN[i]);
          free(BlockDataRVACFOrderNOnsager[i]);
          free(RvacfOrderNTotalOnsager[i]);
          free(RvacfOrderNTotalOnsagerDirAvg[i]);
          free(RvacfOrderNTotalOnsagerCount[i]);
        }
      }
      free(CountRVACFOrderN);
      free(NumberOfBlocksRVACFOrderN);
      free(value_onsager);

      free(BlockLengthRVACFOrderN);
      free(RvacfOrderNCount);
      free(RvacfOrderNOnsagerCount);
      free(RvacfOrderNCountPerMolecule);
      free(RvacfOrderNDirAvg);
      free(RvacfOrderNOnsagerDirAvg);
      free(RvacfOrderNPerMoleculeDirAvg);
      free(RvacfOrderN);
      free(RvacfOrderNOnsager);
      free(RvacfOrderNPerMolecule);
      free(BlockDataRVACFOrderN);
      free(BlockDataRVACFOrderNOnsager);
      free(RvacfOrderNTotalOnsager);
      free(RvacfOrderNTotalOnsagerDirAvg);
      free(RvacfOrderNTotalOnsagerCount);
      break;
  }
}

VECTOR GetOrientationalVectorAdsorbates(int i)
{
  int last;
  VECTOR orientation;
  VECTOR com,t;
  REAL r;
  REAL_MATRIX3x3 M;

  orientation.x=0.0;
  orientation.y=1.0;
  orientation.z=0.0;
  switch(MolecularOrientationType)
  {
    case END_TO_END_VECTOR:
      last=Adsorbates[CurrentSystem][i].NumberOfAtoms;
      orientation.x=Adsorbates[CurrentSystem][i].Atoms[0].Position.x-Adsorbates[CurrentSystem][i].Atoms[last-1].Position.x;
      orientation.y=Adsorbates[CurrentSystem][i].Atoms[0].Position.y-Adsorbates[CurrentSystem][i].Atoms[last-1].Position.y;
      orientation.z=Adsorbates[CurrentSystem][i].Atoms[0].Position.z-Adsorbates[CurrentSystem][i].Atoms[last-1].Position.z;
      r=sqrt(SQR(orientation.x)+SQR(orientation.y)+SQR(orientation.z));
      orientation.x/=r;
      orientation.y/=r;
      orientation.z/=r;
      return orientation;
    default:
    case MOLECULAR_VECTOR:
      BuildRotationMatrixInverse(&M,Adsorbates[CurrentSystem][i].Groups[MolecularOrientationGroup].Quaternion);
      orientation.x=M.ax*MolecularOrientationVector.x+M.bx*MolecularOrientationVector.y+M.cx*MolecularOrientationVector.z;
      orientation.y=M.ay*MolecularOrientationVector.x+M.by*MolecularOrientationVector.y+M.cy*MolecularOrientationVector.z;
      orientation.z=M.az*MolecularOrientationVector.x+M.bz*MolecularOrientationVector.y+M.cz*MolecularOrientationVector.z;
      r=sqrt(SQR(orientation.x)+SQR(orientation.y)+SQR(orientation.z));
      orientation.x/=r;
      orientation.y/=r;
      orientation.z/=r;
      return orientation;
  }
  return orientation;
}

VECTOR GetOrientationalVectorCations(int i)
{
  int last;
  VECTOR orientation;
  VECTOR com,t;
  REAL r;
  REAL_MATRIX3x3 M;

  orientation.x=0.0;
  orientation.y=1.0;
  orientation.z=0.0;
  switch(MolecularOrientationType)
  {
    case END_TO_END_VECTOR:
      last=Cations[CurrentSystem][i].NumberOfAtoms;
      orientation.x=Cations[CurrentSystem][i].Atoms[0].Position.x-Cations[CurrentSystem][i].Atoms[last-1].Position.x;
      orientation.y=Cations[CurrentSystem][i].Atoms[0].Position.y-Cations[CurrentSystem][i].Atoms[last-1].Position.y;
      orientation.z=Cations[CurrentSystem][i].Atoms[0].Position.z-Cations[CurrentSystem][i].Atoms[last-1].Position.z;
      r=sqrt(SQR(orientation.x)+SQR(orientation.y)+SQR(orientation.z));
      orientation.x/=r;
      orientation.y/=r;
      orientation.z/=r;
      return orientation;
    default:
    case MOLECULAR_VECTOR:
      BuildRotationMatrixInverse(&M,Cations[CurrentSystem][i].Groups[MolecularOrientationGroup].Quaternion);
      orientation.x=M.ax*MolecularOrientationVector.x+M.bx*MolecularOrientationVector.y+M.cx*MolecularOrientationVector.z;
      orientation.y=M.ay*MolecularOrientationVector.x+M.by*MolecularOrientationVector.y+M.cy*MolecularOrientationVector.z;
      orientation.z=M.az*MolecularOrientationVector.x+M.bz*MolecularOrientationVector.y+M.cz*MolecularOrientationVector.z;
      r=sqrt(SQR(orientation.x)+SQR(orientation.y)+SQR(orientation.z));
      orientation.x/=r;
      orientation.y/=r;
      orientation.z/=r;
      return orientation;
  }
  return orientation;
}

/*********************************************************************************************************
 * Name       | SampleMolecularOrientationAutoCorrelationFunctionOrderN                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the molecular orientation autocorrelation function using a modified order-N alg. *
 * Parameters | -                                                                                        *
 * Note       | the autocorrelation function of the orientation dependent factor of the dipolar          *
 *            | interaction: A(t)=sqrt(5/8)*(3*cos(theta(t))^2-1)                                        *
 *********************************************************************************************************/

void SampleMolecularOrientationAutoCorrelationFunctionOrderN(int Switch)
{
  int i,j,k,CurrentBlock;
  int CurrentBlocklength,type,shift;
  VECTOR value,temp;
  FILE *FilePtr;
  char buffer[256];
  REAL fac,dt,count,value_dir_avg;

  switch(Switch)
  {
    case ALLOCATE:
      CountMolecularOrientationOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
      NumberOfBlocksMolecularOrientationOrderN=(int*)calloc(NumberOfSystems,sizeof(int));

      BlockLengthMolecularOrientationOrderN=(int**)calloc(NumberOfSystems,sizeof(int*));
      MolecularOrientationOrderNCount=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      MolecularOrientationOrderNDirAvg=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      MolecularOrientationOrderN=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      BlockDataMolecularOrientationOrderN=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeMolecularOrientationOrderN[i])
        {
          BlockLengthMolecularOrientationOrderN[i]=(int*)calloc(MaxNumberOfBlocksMolecularOrientationOrderN,sizeof(int));
          MolecularOrientationOrderNCount[i]=(REAL***)calloc(MaxNumberOfBlocksMolecularOrientationOrderN,sizeof(REAL**));
          MolecularOrientationOrderNDirAvg[i]=(REAL***)calloc(MaxNumberOfBlocksMolecularOrientationOrderN,sizeof(REAL**));
          MolecularOrientationOrderN[i]=(VECTOR***)calloc(MaxNumberOfBlocksMolecularOrientationOrderN,sizeof(VECTOR**));
          BlockDataMolecularOrientationOrderN[i]=(VECTOR***)calloc(MaxNumberOfBlocksMolecularOrientationOrderN,sizeof(VECTOR**));

          for(j=0;j<MaxNumberOfBlocksMolecularOrientationOrderN;j++)
          {
            BlockDataMolecularOrientationOrderN[i][j]=(VECTOR**)calloc(NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR*));
            for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
              BlockDataMolecularOrientationOrderN[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsMolecularOrientationOrderN,sizeof(VECTOR));

            MolecularOrientationOrderNCount[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            MolecularOrientationOrderNDirAvg[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            MolecularOrientationOrderN[i][j]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
            for(k=0;k<NumberOfComponents;k++)
            {
              MolecularOrientationOrderNCount[i][j][k]=(REAL*)calloc(NumberOfBlockElementsMolecularOrientationOrderN,sizeof(REAL));
              MolecularOrientationOrderNDirAvg[i][j][k]=(REAL*)calloc(NumberOfBlockElementsMolecularOrientationOrderN,sizeof(REAL));
              MolecularOrientationOrderN[i][j][k]=(VECTOR*)calloc(NumberOfBlockElementsMolecularOrientationOrderN,sizeof(VECTOR));
            }
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      // return if the vacf does not has to be calculated for this system
      if(!(ComputeMolecularOrientationOrderN[CurrentSystem]&&(CurrentCycle%SampleMolecularOrientationOrderNEvery[CurrentSystem]==0))) return;

      // determine current number of blocks
      NumberOfBlocksMolecularOrientationOrderN[CurrentSystem]=1;
      i=CountMolecularOrientationOrderN[CurrentSystem]/NumberOfBlockElementsMolecularOrientationOrderN;
      while(i!=0)
      {
        NumberOfBlocksMolecularOrientationOrderN[CurrentSystem]++;
        i/=NumberOfBlockElementsMolecularOrientationOrderN;
      }

      // ignore everything beyond the last block
      NumberOfBlocksMolecularOrientationOrderN[CurrentSystem]=MIN2(NumberOfBlocksMolecularOrientationOrderN[CurrentSystem],
                   MaxNumberOfBlocksMolecularOrientationOrderN);

      for(CurrentBlock=0;CurrentBlock<NumberOfBlocksMolecularOrientationOrderN[CurrentSystem];CurrentBlock++)
      {
        // test for blocking operation: CountRVACFOrderN is a multiple of NumberOfBlockElementsRVACFOrderN^CurrentBlock
        if ((CountMolecularOrientationOrderN[CurrentSystem])%((int)pow((REAL)NumberOfBlockElementsMolecularOrientationOrderN,CurrentBlock))==0)
        {
          // increase the current block-length
          BlockLengthMolecularOrientationOrderN[CurrentSystem][CurrentBlock]++;

          // limit length to NumberOfBlockElementsRVACFOrderN
          CurrentBlocklength=MIN2(BlockLengthMolecularOrientationOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsMolecularOrientationOrderN);

          // Self diffusion
          // ========================================================================================================================

          shift=NumberOfAdsorbateMolecules[CurrentSystem];
          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];i++)
          {
            if(i<NumberOfAdsorbateMolecules[CurrentSystem])
              value=GetOrientationalVectorAdsorbates(i);
            else
              value=GetOrientationalVectorCations(i-shift);

            for(j=CurrentBlocklength-1;j>0;j--)
              BlockDataMolecularOrientationOrderN[CurrentSystem][CurrentBlock][i][j]=BlockDataMolecularOrientationOrderN[CurrentSystem][CurrentBlock][i][j-1];
            BlockDataMolecularOrientationOrderN[CurrentSystem][CurrentBlock][i][0]=value;


            // get the type of the molecule
            if(i<NumberOfAdsorbateMolecules[CurrentSystem])
              type=Adsorbates[CurrentSystem][i].Type;
            else
              type=Cations[CurrentSystem][i-shift].Type;

            for(j=0;j<CurrentBlocklength;j++)
            {
              // vacf for each component
              MolecularOrientationOrderNCount[CurrentSystem][CurrentBlock][type][j]+=1.0;
              temp.x=BlockDataMolecularOrientationOrderN[CurrentSystem][CurrentBlock][i][j].x*value.x;
              temp.y=BlockDataMolecularOrientationOrderN[CurrentSystem][CurrentBlock][i][j].y*value.y;
              temp.z=BlockDataMolecularOrientationOrderN[CurrentSystem][CurrentBlock][i][j].z*value.z;

              MolecularOrientationOrderN[CurrentSystem][CurrentBlock][type][j].x+=temp.x;
              MolecularOrientationOrderN[CurrentSystem][CurrentBlock][type][j].y+=temp.y;
              MolecularOrientationOrderN[CurrentSystem][CurrentBlock][type][j].z+=temp.z;
              MolecularOrientationOrderNDirAvg[CurrentSystem][CurrentBlock][type][j]+=temp.x+temp.y+temp.z;
            }
          }
        }
      }
      // CountRVACFOrderN the current sampling
      CountMolecularOrientationOrderN[CurrentSystem]++;
      break;
    case PRINT:
      // return if the vacf does not has to be calculated for this system
      if((!ComputeMolecularOrientationOrderN[CurrentSystem])||(CurrentCycle%WriteMolecularOrientationOrderNEvery[CurrentSystem]!=0)) return;

      mkdir("MolecularOrientationOrderN",S_IRWXU);

      sprintf(buffer,"MolecularOrientationOrderN/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      // Self diffusion
      // ========================================================================================================================

      sprintf(buffer,"MolecularOrientationOrderN/System_%d/molecular_orientation_self_total%s.dat",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"# column 1: time [ps]\n");
      fprintf(FilePtr,"# column 2: autocorrelation of the bond averaged over x,y, and z\n");
      fprintf(FilePtr,"# column 3: autocorrelation of the bond in the x-direction\n");
      fprintf(FilePtr,"# column 4: autocorrelation of the bond in the y-direction\n");
      fprintf(FilePtr,"# column 5: autocorrelation of the bond in the z-direction\n");
      fprintf(FilePtr,"# column 6: the number of measurements (per point)\n");


      for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksMolecularOrientationOrderN,NumberOfBlocksMolecularOrientationOrderN[CurrentSystem]);CurrentBlock++)
      {
        CurrentBlocklength=MIN2(BlockLengthMolecularOrientationOrderN[CurrentSystem][CurrentBlock],NumberOfBlockElementsMolecularOrientationOrderN);
        dt=SampleMolecularOrientationOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsMolecularOrientationOrderN,CurrentBlock);
        for(j=1;j<CurrentBlocklength;j++)
        {
          count=0.0;
          value_dir_avg=0.0;
          value.x=value.y=value.z=0.0;
          for(i=0;i<NumberOfComponents;i++)
          {
            value.x+=MolecularOrientationOrderN[CurrentSystem][CurrentBlock][i][j].x;
            value.y+=MolecularOrientationOrderN[CurrentSystem][CurrentBlock][i][j].y;
            value.z+=MolecularOrientationOrderN[CurrentSystem][CurrentBlock][i][j].z;
            value_dir_avg+=MolecularOrientationOrderNDirAvg[CurrentSystem][CurrentBlock][i][j];
            count+=MolecularOrientationOrderNCount[CurrentSystem][CurrentBlock][i][j];
          }

          if(count>0.0)
          {
            fac=1.0/count;

            fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                (double)(j*dt),
                (double)(fac*value_dir_avg),
                (double)(fac*value.x),
                (double)(fac*value.y),
                (double)(fac*value.z),
                (double)count);
          }
        }
      }

      fclose(FilePtr);
      break;
    case FINALIZE:
      break;
  }
}

VECTOR BondOrientation(int system,int f1,int i,int j)
{
  int A,B;
  VECTOR posA,posB,dr;
  VECTOR CosTheta,factor;
  VECTOR orientation;
  VECTOR e1,e2,e3;
  REAL rr,r;

  A=OrientationFrameworkBondPairs[system][f1][i][j].A;
  B=OrientationFrameworkBondPairs[system][f1][i][j].B;

  posA=Framework[system].Atoms[f1][A].Position;
  posB=Framework[system].Atoms[f1][B].Position;

  orientation.x=posA.x-posB.x;
  orientation.y=posA.y-posB.y;
  orientation.z=posA.z-posB.z;
  orientation=ApplyBoundaryCondition(orientation);
  rr=SQR(orientation.x)+SQR(orientation.y)+SQR(orientation.z);
  r=sqrt(rr);

  orientation.x/=r;
  orientation.y/=r;
  orientation.z/=r;

  return orientation;
}



/*********************************************************************************************************
 * Name       | SampleBondOrientationAutoCorrelationFunctionOrderN                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the molecular orientation autocorrelation function using a modified order-N alg. *
 * Parameters | -                                                                                        *
 * Note       | the autocorrelation function of the orientation dependent factor of the dipolar          *
 *            | interaction: A(t)=sqrt(5/8)*(3*cos(theta(t))^2-1)                                        *
 *********************************************************************************************************/

// fit to f(x)=A*exp(-(t/tau_a)**beta)
// The relaxation time tau_a, the exponent beta, and the non-ergodicity factor A are fitting parameters that depend on temperature T and density rho

void SampleBondOrientationAutoCorrelationFunctionOrderN(int Switch)
{
  int i,j,k,l,CurrentBlock;
  int CurrentBlocklength;
  VECTOR value,norm;
  FILE *FilePtr;
  char buffer[256];
  REAL fac,dt,count,value_dir_avg;
  int A,B,f1,index;
  int typeA,typeB;
  int NumberOfTypes;
  REAL angle;

  switch(Switch)
  {
    case ALLOCATE:
      CountBondOrientationOrderN=(int**)calloc(NumberOfSystems,sizeof(int*));
      NumberOfBlocksBondOrientationOrderN=(int**)calloc(NumberOfSystems,sizeof(int*));

      BlockLengthBondOrientationOrderN=(int***)calloc(NumberOfSystems,sizeof(int**));
      BondOrientationOrderNCount=(REAL*****)calloc(NumberOfSystems,sizeof(REAL****));
      BondOrientationOrderNDirAvg=(REAL*****)calloc(NumberOfSystems,sizeof(REAL****));
      BondOrientationOrderN=(VECTOR*****)calloc(NumberOfSystems,sizeof(VECTOR****));
      BlockDataBondOrientationOrderN=(VECTOR******)calloc(NumberOfSystems,sizeof(VECTOR*****));

      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeBondOrientationOrderN[i])
        {
          CountBondOrientationOrderN[i]=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
          NumberOfBlocksBondOrientationOrderN[i]=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));

          BlockLengthBondOrientationOrderN[i]=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
          BondOrientationOrderNCount[i]=(REAL****)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL***));
          BondOrientationOrderNDirAvg[i]=(REAL****)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL***));
          BondOrientationOrderN[i]=(VECTOR****)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR***));
          BlockDataBondOrientationOrderN[i]=(VECTOR*****)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR****));

          for(f1=0;f1<Framework[i].NumberOfFrameworks;f1++)
          {
            BlockLengthBondOrientationOrderN[i][f1]=(int*)calloc(MaxNumberOfBlocksBondOrientationOrderN,sizeof(int));
            BondOrientationOrderNCount[i][f1]=(REAL***)calloc(MaxNumberOfBlocksBondOrientationOrderN,sizeof(REAL**));
            BondOrientationOrderNDirAvg[i][f1]=(REAL***)calloc(MaxNumberOfBlocksBondOrientationOrderN,sizeof(REAL**));
            BondOrientationOrderN[i][f1]=(VECTOR***)calloc(MaxNumberOfBlocksBondOrientationOrderN,sizeof(VECTOR**));
            BlockDataBondOrientationOrderN[i][f1]=(VECTOR****)calloc(MaxNumberOfBlocksBondOrientationOrderN,sizeof(VECTOR***));

            NumberOfTypes=NumberOfOrientationFrameworkBonds[i][f1];
            BondOrientationAngleDistributionFunction[i][f1]=(VECTOR**)calloc(NumberOfTypes,sizeof(VECTOR*));
            for(k=0;k<NumberOfTypes;k++)
              BondOrientationAngleDistributionFunction[i][f1][k]=(VECTOR*)calloc(BondOrientationAngleHistogramSize[i],sizeof(VECTOR));

            for(j=0;j<MaxNumberOfBlocksBondOrientationOrderN;j++)
            {
              BlockDataBondOrientationOrderN[i][f1][j]=(VECTOR***)calloc(NumberOfTypes,sizeof(VECTOR**));
              for(k=0;k<NumberOfTypes;k++)
              {
                BlockDataBondOrientationOrderN[i][f1][j][k]=(VECTOR**)calloc(NumberOfOrientationFrameworkBondPairs[i][f1][k],sizeof(VECTOR*));

                for(l=0;l<NumberOfOrientationFrameworkBondPairs[i][f1][k];l++)
                  BlockDataBondOrientationOrderN[i][f1][j][k][l]=(VECTOR*)calloc(NumberOfBlockElementsBondOrientationOrderN,sizeof(VECTOR));
              }

              BondOrientationOrderNCount[i][f1][j]=(REAL**)calloc(NumberOfTypes,sizeof(REAL*));
              BondOrientationOrderNDirAvg[i][f1][j]=(REAL**)calloc(NumberOfTypes,sizeof(REAL*));
              BondOrientationOrderN[i][f1][j]=(VECTOR**)calloc(NumberOfTypes,sizeof(VECTOR*));
              for(k=0;k<NumberOfTypes;k++)
              {
                BondOrientationOrderNCount[i][f1][j][k]=(REAL*)calloc(NumberOfBlockElementsBondOrientationOrderN,sizeof(REAL));
                BondOrientationOrderNDirAvg[i][f1][j][k]=(REAL*)calloc(NumberOfBlockElementsBondOrientationOrderN,sizeof(REAL));
                BondOrientationOrderN[i][f1][j][k]=(VECTOR*)calloc(NumberOfBlockElementsBondOrientationOrderN,sizeof(VECTOR));
              }
            }
          }
        }
      }
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      // return if the vacf does not has to be calculated for this system
      if(!(ComputeBondOrientationOrderN[CurrentSystem]&&(CurrentCycle%SampleBondOrientationOrderNEvery[CurrentSystem]==0))) return;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<NumberOfOrientationFrameworkBonds[CurrentSystem][f1];k++)
        {
          for(l=0;l<NumberOfOrientationFrameworkBondPairs[CurrentSystem][f1][k];l++)
          {
            value=BondOrientation(CurrentSystem,f1,k,l);

            angle=atan2(value.y,value.x)*RAD2DEG;
            if(angle<0) angle+=360.0;
            index=(int)((angle/360.0)*(REAL)BondOrientationAngleHistogramSize[CurrentSystem]);
            BondOrientationAngleDistributionFunction[CurrentSystem][f1][k][index].z+=1.0;

            angle=atan2(value.z,value.x)*RAD2DEG;
            if(angle<0) angle+=360.0;
            index=(int)((angle/360.0)*(REAL)BondOrientationAngleHistogramSize[CurrentSystem]);
            BondOrientationAngleDistributionFunction[CurrentSystem][f1][k][index].y+=1.0;

            angle=atan2(value.z,value.y)*RAD2DEG;
            if(angle<0) angle+=360.0;
            index=(int)((angle/360.0)*(REAL)BondOrientationAngleHistogramSize[CurrentSystem]);
            BondOrientationAngleDistributionFunction[CurrentSystem][f1][k][index].x+=1.0;
          }
        }


        // determine current number of blocks
        NumberOfBlocksBondOrientationOrderN[CurrentSystem][f1]=1;
        i=CountBondOrientationOrderN[CurrentSystem][f1]/NumberOfBlockElementsBondOrientationOrderN;
        while(i!=0)
        {
          NumberOfBlocksBondOrientationOrderN[CurrentSystem][f1]++;
          i/=NumberOfBlockElementsBondOrientationOrderN;
        }

        // ignore everything beyond the last block
        NumberOfBlocksBondOrientationOrderN[CurrentSystem][f1]=MIN2(NumberOfBlocksBondOrientationOrderN[CurrentSystem][f1],
                     MaxNumberOfBlocksBondOrientationOrderN);

        for(CurrentBlock=0;CurrentBlock<NumberOfBlocksBondOrientationOrderN[CurrentSystem][f1];CurrentBlock++)
        {
          // test for blocking operation: CountBondOrientationOrderN is a multiple of NumberOfBlockElementsBondOrientationOrderN^CurrentBlock
          if ((CountBondOrientationOrderN[CurrentSystem][f1])%((int)pow((REAL)NumberOfBlockElementsBondOrientationOrderN,CurrentBlock))==0)
          {
            // increase the current block-length
            BlockLengthBondOrientationOrderN[CurrentSystem][f1][CurrentBlock]++;

            // limit length to NumberOfBlockElementsRVACFOrderN
            CurrentBlocklength=MIN2(BlockLengthBondOrientationOrderN[CurrentSystem][f1][CurrentBlock],NumberOfBlockElementsBondOrientationOrderN);

            // Self diffusion
            // ========================================================================================================================

            for(k=0;k<NumberOfOrientationFrameworkBonds[CurrentSystem][f1];k++)
            {
              for(l=0;l<NumberOfOrientationFrameworkBondPairs[CurrentSystem][f1][k];l++)
              {
                value=BondOrientation(CurrentSystem,f1,k,l);

                for(j=CurrentBlocklength-1;j>0;j--)
                  BlockDataBondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][l][j]=BlockDataBondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][l][j-1];
                BlockDataBondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][l][0]=value;


                for(j=0;j<CurrentBlocklength;j++)
                {
                  // vacf for each component
                  BondOrientationOrderNCount[CurrentSystem][f1][CurrentBlock][k][j]+=1.0;
                  BondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][j].x+=
                    (BlockDataBondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][l][j].x*value.x);
                  BondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][j].y+=
                    (BlockDataBondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][l][j].y*value.y);
                  BondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][j].z+=
                    (BlockDataBondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][l][j].z*value.z);
                  BondOrientationOrderNDirAvg[CurrentSystem][f1][CurrentBlock][k][j]+=
                     (BlockDataBondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][l][j].x*value.x)+
                     (BlockDataBondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][l][j].y*value.y)+
                     (BlockDataBondOrientationOrderN[CurrentSystem][f1][CurrentBlock][k][l][j].z*value.z);
                }
              }
            }
          }
        }
        // CountRVACFOrderN the current sampling
        CountBondOrientationOrderN[CurrentSystem][f1]++;
      }
      break;
    case PRINT:
      // return if the boacf does not has to be calculated for this system
      if((!ComputeBondOrientationOrderN[CurrentSystem])||(CurrentCycle%WriteBondOrientationOrderNEvery[CurrentSystem]!=0)) return;

      mkdir("BondOrientationOrderN",S_IRWXU);

      sprintf(buffer,"BondOrientationOrderN/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        sprintf(buffer,"BondOrientationOrderN/System_%d/Framework_%d",CurrentSystem,f1);
        mkdir(buffer,S_IRWXU);

        for(i=0;i<NumberOfOrientationFrameworkBonds[CurrentSystem][f1];i++)
        {
          norm.x=norm.y=norm.z=0.0;
          for(k=0;k<BondOrientationAngleHistogramSize[CurrentSystem];k++)
          {
            norm.x+=BondOrientationAngleDistributionFunction[CurrentSystem][f1][i][k].x;
            norm.y+=BondOrientationAngleDistributionFunction[CurrentSystem][f1][i][k].y;
            norm.z+=BondOrientationAngleDistributionFunction[CurrentSystem][f1][i][k].z;
          }
          norm.x*=360.0/(REAL)BondOrientationAngleHistogramSize[CurrentSystem];
          norm.y*=360.0/(REAL)BondOrientationAngleHistogramSize[CurrentSystem];
          norm.z*=360.0/(REAL)BondOrientationAngleHistogramSize[CurrentSystem];

          sprintf(buffer,"BondOrientationOrderN/System_%d/Framework_%d/Orientation_histogram%d%s.dat",CurrentSystem,f1,i,FileNameAppend);
          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: angle [degrees]\n");
          fprintf(FilePtr,"# column 2: histogram of the angle of the vector with the y-axis in the z-y plane\n");
          fprintf(FilePtr,"# column 3: histogram of the angle of the vector with the x-axis in the z-x plane\n");
          fprintf(FilePtr,"# column 4: histogram of the angle of the vector with the x-axis in the y-x plane\n");
          for(k=0;k<BondOrientationAngleHistogramSize[CurrentSystem];k++)
          {
            fprintf(FilePtr,"%g %g %g %g\n",
               k*360.0/BondOrientationAngleHistogramSize[CurrentSystem],
               BondOrientationAngleDistributionFunction[CurrentSystem][f1][i][k].x/norm.x,
               BondOrientationAngleDistributionFunction[CurrentSystem][f1][i][k].y/norm.y,
               BondOrientationAngleDistributionFunction[CurrentSystem][f1][i][k].z/norm.z);
          }
          fclose(FilePtr);

          sprintf(buffer,"BondOrientationOrderN/System_%d/Framework_%d/Orientation_correlation_type%d%s.dat",CurrentSystem,f1,i,FileNameAppend);
          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: autocorrelation of the bond averaged over x,y, and z\n");
          fprintf(FilePtr,"# column 3: autocorrelation of the bond in the x-direction\n");
          fprintf(FilePtr,"# column 4: autocorrelation of the bond in the y-direction\n");
          fprintf(FilePtr,"# column 5: autocorrelation of the bond in the z-direction\n");
          fprintf(FilePtr,"# column 6: the number of measurements (per point)\n");
          for(CurrentBlock=0;CurrentBlock<MIN2(MaxNumberOfBlocksBondOrientationOrderN,NumberOfBlocksBondOrientationOrderN[CurrentSystem][f1]);CurrentBlock++)
          {
            CurrentBlocklength=MIN2(BlockLengthBondOrientationOrderN[CurrentSystem][f1][CurrentBlock],NumberOfBlockElementsBondOrientationOrderN);
            dt=SampleBondOrientationOrderNEvery[CurrentSystem]*DeltaT*pow((REAL)NumberOfBlockElementsBondOrientationOrderN,CurrentBlock);
            for(j=1;j<CurrentBlocklength;j++)
            {
              value=BondOrientationOrderN[CurrentSystem][f1][CurrentBlock][i][j];
              value_dir_avg=BondOrientationOrderNDirAvg[CurrentSystem][f1][CurrentBlock][i][j];
              count=BondOrientationOrderNCount[CurrentSystem][f1][CurrentBlock][i][j];

              if(count>0.0)
              {
                fac=1.0/count;

                fprintf(FilePtr,"%g %g %g %g %g (count: %g)\n",
                    (double)(j*dt),
                    (double)(fac*value_dir_avg),
                    (double)(fac*value.x),
                    (double)(fac*value.y),
                    (double)(fac*value.z),
                    (double)count);
              }
            }
          }
          fclose(FilePtr);
        }
      }
      break;
    case FINALIZE:
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleMeanSquaredDisplacement                                                            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the mean-square displacement  function using a conventional algorithm.           *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleMeanSquaredDisplacement(int Switch)
{
  int i,j,k,l;
  int type,index,CurrentBuffer;
  VECTOR pos,drift,com;
  static VECTOR *sum_pos;
  REAL fac;
  FILE *FilePtr;
  char buffer[256];

  switch(Switch)
  {
    case ALLOCATE:
      OriginMSD=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      OriginOnsagerMSD=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      AcfMSD=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      AcfOnsagerMSD=(VECTOR*****)calloc(NumberOfSystems,sizeof(VECTOR****));
      AcfDirAvgMSD=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      AcfOnsagerDirAvgMSD=(REAL*****)calloc(NumberOfSystems,sizeof(REAL****));
      AccumulatedAcfMSD=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      AccumulatedAcfDirAvgMSD=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      AccumulatedAcfOnsagerMSD=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      AccumulatedAcfOnsagerDirAvgMSD=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      CountMSD=(int**)calloc(NumberOfSystems,sizeof(int*));
      CountAccumulatedMSD=(int*)calloc(NumberOfSystems,sizeof(int));
      sum_pos=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeMSD[i])
        {
          OriginMSD[i]=(VECTOR**)calloc(NumberOfBuffersMSD,sizeof(VECTOR*));
          OriginOnsagerMSD[i]=(VECTOR**)calloc(NumberOfBuffersMSD,sizeof(VECTOR*));
          CountMSD[i]=(int*)calloc(NumberOfBuffersMSD,sizeof(int));

          AcfMSD[i]=(VECTOR***)calloc(NumberOfBuffersMSD,sizeof(VECTOR**));
          AcfOnsagerMSD[i]=(VECTOR****)calloc(NumberOfBuffersMSD,sizeof(VECTOR***));
          AcfDirAvgMSD[i]=(REAL***)calloc(NumberOfBuffersMSD,sizeof(REAL**));
          AcfOnsagerDirAvgMSD[i]=(REAL****)calloc(NumberOfBuffersMSD,sizeof(REAL***));

          for(j=0;j<NumberOfBuffersMSD;j++)
          {
            AcfMSD[i][j]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
            AcfOnsagerMSD[i][j]=(VECTOR***)calloc(NumberOfComponents,sizeof(VECTOR**));
            AcfDirAvgMSD[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            AcfOnsagerDirAvgMSD[i][j]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));
            OriginMSD[i][j]=(VECTOR*)calloc(NumberOfAdsorbateMolecules[i],sizeof(VECTOR));
            OriginOnsagerMSD[i][j]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
            for(k=0;k<NumberOfComponents;k++)
            {
              AcfMSD[i][j][k]=(VECTOR*)calloc(BufferLengthMSD,sizeof(VECTOR));
              AcfDirAvgMSD[i][j][k]=(REAL*)calloc(BufferLengthMSD,sizeof(REAL));

              AcfOnsagerMSD[i][j][k]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
              AcfOnsagerDirAvgMSD[i][j][k]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
              for(l=0;l<NumberOfComponents;l++)
              {
                AcfOnsagerMSD[i][j][k][l]=(VECTOR*)calloc(BufferLengthMSD,sizeof(VECTOR));
                AcfOnsagerDirAvgMSD[i][j][k][l]=(REAL*)calloc(BufferLengthMSD,sizeof(REAL));
              }
            }
          }

          AccumulatedAcfMSD[i]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
          AccumulatedAcfDirAvgMSD[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
          AccumulatedAcfOnsagerMSD[i]=(VECTOR***)calloc(NumberOfComponents,sizeof(VECTOR**));
          AccumulatedAcfOnsagerDirAvgMSD[i]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));

          for(j=0;j<NumberOfComponents;j++)
          {
            AccumulatedAcfMSD[i][j]=(VECTOR*)calloc(BufferLengthMSD,sizeof(VECTOR));
            AccumulatedAcfDirAvgMSD[i][j]=(REAL*)calloc(BufferLengthMSD,sizeof(REAL));

            AccumulatedAcfOnsagerMSD[i][j]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
            AccumulatedAcfOnsagerDirAvgMSD[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            for(k=0;k<NumberOfComponents;k++)
            {
              AccumulatedAcfOnsagerMSD[i][j][k]=(VECTOR*)calloc(BufferLengthMSD,sizeof(VECTOR));
              AccumulatedAcfOnsagerDirAvgMSD[i][j][k]=(REAL*)calloc(BufferLengthMSD,sizeof(REAL));
            }
          }
        }
      }
      break;
    case INITIALIZE:
      // trick to space the origins evenly (see for example Rapaport 2004)
      if(ComputeMSD[CurrentSystem])
        for(CurrentBuffer=0;CurrentBuffer<NumberOfBuffersMSD;CurrentBuffer++)
          CountMSD[CurrentSystem][CurrentBuffer]=(int)(-CurrentBuffer*BufferLengthMSD/NumberOfBuffersMSD);
      break;
    case SAMPLE:
      if(!(ComputeMSD[CurrentSystem]&&CurrentCycle%SampleMSDEvery[CurrentSystem]==0)) return;

      // compute drift of the system
      if(Framework[CurrentSystem].FrameworkModel==NONE)
      {
        com=GetCenterOfMassCurrentSystem();
        drift.x=com.x-IntialCenterOfMassPosition[CurrentSystem].x;
        drift.y=com.y-IntialCenterOfMassPosition[CurrentSystem].y;
        drift.z=com.z-IntialCenterOfMassPosition[CurrentSystem].z;
      }
      else if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        com=GetFrameworkCenterOfMass();
        drift.x=com.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x;
        drift.y=com.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y;
        drift.z=com.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z;
      }
      else
      {
        drift.x=0.0;
        drift.y=0.0;
        drift.z=0.0;
      }


      for(CurrentBuffer=0;CurrentBuffer<NumberOfBuffersMSD;CurrentBuffer++)
      {
        if(CountMSD[CurrentSystem][CurrentBuffer]==0)
        {
          for(i=0;i<NumberOfComponents;i++)
            sum_pos[i].x=sum_pos[i].y=sum_pos[i].z=0.0;
          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
          {
            type=Adsorbates[CurrentSystem][i].Type;
            pos=GetAdsorbateCenterOfMass(i);
            OriginMSD[CurrentSystem][CurrentBuffer][i].x=(pos.x-drift.x);
            OriginMSD[CurrentSystem][CurrentBuffer][i].y=(pos.y-drift.y);
            OriginMSD[CurrentSystem][CurrentBuffer][i].z=(pos.z-drift.z);

            sum_pos[type].x+=(pos.x-drift.x);
            sum_pos[type].y+=(pos.y-drift.y);
            sum_pos[type].z+=(pos.z-drift.z);
          }
          for(i=0;i<NumberOfComponents;i++)
            OriginOnsagerMSD[CurrentSystem][CurrentBuffer][i]=sum_pos[i];
        }
        if(CountMSD[CurrentSystem][CurrentBuffer]>=0)
        {
          index=(int)CountMSD[CurrentSystem][CurrentBuffer];

          for(i=0;i<NumberOfComponents;i++)
          {
            AcfMSD[CurrentSystem][CurrentBuffer][i][index].x=0.0;
            AcfMSD[CurrentSystem][CurrentBuffer][i][index].y=0.0;
            AcfMSD[CurrentSystem][CurrentBuffer][i][index].z=0.0;
            AcfDirAvgMSD[CurrentSystem][CurrentBuffer][i][index]=0.0;
          }

          for(i=0;i<NumberOfComponents;i++)
            sum_pos[i].x=sum_pos[i].y=sum_pos[i].z=0.0;
          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
          {
            type=Adsorbates[CurrentSystem][i].Type;
            pos=GetAdsorbateCenterOfMass(i);
            pos.x-=drift.x;
            pos.y-=drift.y;
            pos.z-=drift.z;
            AcfMSD[CurrentSystem][CurrentBuffer][type][index].x+=SQR(pos.x-OriginMSD[CurrentSystem][CurrentBuffer][i].x);
            AcfMSD[CurrentSystem][CurrentBuffer][type][index].y+=SQR(pos.y-OriginMSD[CurrentSystem][CurrentBuffer][i].y);
            AcfMSD[CurrentSystem][CurrentBuffer][type][index].z+=SQR(pos.z-OriginMSD[CurrentSystem][CurrentBuffer][i].z);
            AcfDirAvgMSD[CurrentSystem][CurrentBuffer][type][index]+=SQR(pos.x-OriginMSD[CurrentSystem][CurrentBuffer][i].x)+
                                                                     SQR(pos.y-OriginMSD[CurrentSystem][CurrentBuffer][i].y)+
                                                                     SQR(pos.z-OriginMSD[CurrentSystem][CurrentBuffer][i].z);
            sum_pos[type].x+=pos.x;
            sum_pos[type].y+=pos.y;
            sum_pos[type].z+=pos.z;
          }
          for(i=0;i<NumberOfComponents;i++)
          {
            for(j=0;j<NumberOfComponents;j++)
            {
              AcfOnsagerMSD[CurrentSystem][CurrentBuffer][i][j][index].x=(sum_pos[i].x-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][i].x)*
                                                                         (sum_pos[j].x-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][j].x);
              AcfOnsagerMSD[CurrentSystem][CurrentBuffer][i][j][index].y=(sum_pos[i].y-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][i].y)*
                                                                         (sum_pos[j].y-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][j].y);
              AcfOnsagerMSD[CurrentSystem][CurrentBuffer][i][j][index].z=(sum_pos[i].z-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][i].z)*
                                                                         (sum_pos[j].z-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][j].z);
              AcfOnsagerDirAvgMSD[CurrentSystem][CurrentBuffer][i][j][index]=
                 (sum_pos[i].x-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][i].x)*(sum_pos[j].x-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][j].x)+
                 (sum_pos[i].y-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][i].y)*(sum_pos[j].y-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][j].y)+
                 (sum_pos[i].z-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][i].z)*(sum_pos[j].z-OriginOnsagerMSD[CurrentSystem][CurrentBuffer][j].z);
            }
          }
        }
        CountMSD[CurrentSystem][CurrentBuffer]++;
      }

      // accumulate the msd
      for(CurrentBuffer=0;CurrentBuffer<NumberOfBuffersMSD;CurrentBuffer++)
      {
        if(CountMSD[CurrentSystem][CurrentBuffer]==BufferLengthMSD)
        {
          for(k=0;k<NumberOfComponents;k++)
          {
            for(i=0;i<BufferLengthMSD;i++)
            {
              AccumulatedAcfMSD[CurrentSystem][k][i].x+=AcfMSD[CurrentSystem][CurrentBuffer][k][i].x;
              AccumulatedAcfMSD[CurrentSystem][k][i].y+=AcfMSD[CurrentSystem][CurrentBuffer][k][i].y;
              AccumulatedAcfMSD[CurrentSystem][k][i].z+=AcfMSD[CurrentSystem][CurrentBuffer][k][i].z;
              AccumulatedAcfDirAvgMSD[CurrentSystem][k][i]+=AcfDirAvgMSD[CurrentSystem][CurrentBuffer][k][i];
            }
          }
          for(k=0;k<NumberOfComponents;k++)
          {
            for(l=0;l<NumberOfComponents;l++)
            {
              for(i=0;i<BufferLengthMSD;i++)
              {
                AccumulatedAcfOnsagerMSD[CurrentSystem][k][l][i].x+=AcfOnsagerMSD[CurrentSystem][CurrentBuffer][k][l][i].x;
                AccumulatedAcfOnsagerMSD[CurrentSystem][k][l][i].y+=AcfOnsagerMSD[CurrentSystem][CurrentBuffer][k][l][i].y;
                AccumulatedAcfOnsagerMSD[CurrentSystem][k][l][i].z+=AcfOnsagerMSD[CurrentSystem][CurrentBuffer][k][l][i].z;
                AccumulatedAcfOnsagerDirAvgMSD[CurrentSystem][k][l][i]+=AcfOnsagerDirAvgMSD[CurrentSystem][CurrentBuffer][k][l][i];
              }
            }
          }
          CountMSD[CurrentSystem][CurrentBuffer]=0;
          CountAccumulatedMSD[CurrentSystem]++;
        }
      }

      break;
    case PRINT:
      if(!(ComputeMSD[CurrentSystem]&&(CurrentCycle%WriteMSDEvery[CurrentSystem]==0))) return;

      mkdir("MSD",S_IRWXU);

      sprintf(buffer,"MSD/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(k=0;k<NumberOfComponents;k++)
      {
        fac=1.0/(Components[k].NumberOfMolecules[CurrentSystem]*CountAccumulatedMSD[CurrentSystem]);

        sprintf(buffer,"MSD/System_%d/msd_self_%s_%d%s.dat",
                CurrentSystem,
                Components[k].Name,
                k,
                FileNameAppend);

        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# column 1: time [ps]\n");
        fprintf(FilePtr,"# column 2: msd xyz [A^2]\n");
        fprintf(FilePtr,"# column 3: msd x [A^2]\n");
        fprintf(FilePtr,"# column 4: msd y [A^2]\n");
        fprintf(FilePtr,"# column 5: msd z [A^2]\n");

        for(i=0;i<BufferLengthMSD;i++)
          fprintf(FilePtr,"%g %g %g %g %g\n",i*SampleMSDEvery[CurrentSystem]*DeltaT,
            fac*AccumulatedAcfDirAvgMSD[CurrentSystem][k][i],
            fac*AccumulatedAcfMSD[CurrentSystem][k][i].x,
            fac*AccumulatedAcfMSD[CurrentSystem][k][i].y,
            fac*AccumulatedAcfMSD[CurrentSystem][k][i].z);
        fclose(FilePtr);
      }

      for(k=0;k<NumberOfComponents;k++)
      {
        fac=1.0/(Components[k].NumberOfMolecules[CurrentSystem]*CountAccumulatedMSD[CurrentSystem]);
        for(l=0;l<NumberOfComponents;l++)
        {
          sprintf(buffer,"MSD/System_%d/msd_onsager_%s_%d_%s_%d%s.dat",
                  CurrentSystem,Components[k].Name,k,Components[l].Name,l,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: msd xyz [A^2]\n");
          fprintf(FilePtr,"# column 3: msd x [A^2]\n");
          fprintf(FilePtr,"# column 4: msd y [A^2]\n");
          fprintf(FilePtr,"# column 5: msd z [A^2]\n");

          for(i=0;i<BufferLengthMSD;i++)
            fprintf(FilePtr,"%g %g %g %g %g\n",i*SampleMSDEvery[CurrentSystem]*DeltaT,
              fac*AccumulatedAcfOnsagerDirAvgMSD[CurrentSystem][k][l][i],
              fac*AccumulatedAcfOnsagerMSD[CurrentSystem][k][l][i].x,
              fac*AccumulatedAcfOnsagerMSD[CurrentSystem][k][l][i].y,
              fac*AccumulatedAcfOnsagerMSD[CurrentSystem][k][l][i].z);
          fclose(FilePtr);
        }
      }
      break;
    case FINALIZE:
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleVelocityAutoCorrelationFunction                                                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the velocity autocorrelation function using a conventional algorithm.            *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleVelocityAutoCorrelationFunction(int Switch)
{
  int i,j,k,l;
  int type,index,CurrentBuffer;
  VECTOR vel;
  static VECTOR *sum_vel;
  REAL fac,D;
  FILE *FilePtr;
  char buffer[256];

  switch(Switch)
  {
    case ALLOCATE:
      OriginVACF=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      OriginOnsagerVACF=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      AcfVACF=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      AcfOnsagerVACF=(VECTOR*****)calloc(NumberOfSystems,sizeof(VECTOR****));
      AcfDirAvgVACF=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      AcfOnsagerDirAvgVACF=(REAL*****)calloc(NumberOfSystems,sizeof(REAL****));
      AccumulatedAcfVACF=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
      AccumulatedAcfDirAvgVACF=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      AccumulatedAcfOnsagerVACF=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));
      AccumulatedAcfOnsagerDirAvgVACF=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));
      CountVACF=(int**)calloc(NumberOfSystems,sizeof(int*));
      CountAccumulatedVACF=(int*)calloc(NumberOfSystems,sizeof(int));
      sum_vel=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeVACF[i])
        {
          OriginVACF[i]=(VECTOR**)calloc(NumberOfBuffersVACF,sizeof(VECTOR*));
          OriginOnsagerVACF[i]=(VECTOR**)calloc(NumberOfBuffersVACF,sizeof(VECTOR*));
          CountVACF[i]=(int*)calloc(NumberOfBuffersVACF,sizeof(int));

          AcfVACF[i]=(VECTOR***)calloc(NumberOfBuffersVACF,sizeof(VECTOR**));
          AcfOnsagerVACF[i]=(VECTOR****)calloc(NumberOfBuffersVACF,sizeof(VECTOR***));
          AcfDirAvgVACF[i]=(REAL***)calloc(NumberOfBuffersVACF,sizeof(REAL**));
          AcfOnsagerDirAvgVACF[i]=(REAL****)calloc(NumberOfBuffersVACF,sizeof(REAL***));

          for(j=0;j<NumberOfBuffersVACF;j++)
          {
            AcfVACF[i][j]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
            AcfOnsagerVACF[i][j]=(VECTOR***)calloc(NumberOfComponents,sizeof(VECTOR**));
            AcfDirAvgVACF[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            AcfOnsagerDirAvgVACF[i][j]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));
            OriginVACF[i][j]=(VECTOR*)calloc(NumberOfAdsorbateMolecules[i],sizeof(VECTOR));
            OriginOnsagerVACF[i][j]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
            for(k=0;k<NumberOfComponents;k++)
            {
              AcfVACF[i][j][k]=(VECTOR*)calloc(BufferLengthVACF,sizeof(VECTOR));
              AcfDirAvgVACF[i][j][k]=(REAL*)calloc(BufferLengthVACF,sizeof(REAL));

              AcfOnsagerVACF[i][j][k]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
              AcfOnsagerDirAvgVACF[i][j][k]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
              for(l=0;l<NumberOfComponents;l++)
              {
                AcfOnsagerVACF[i][j][k][l]=(VECTOR*)calloc(BufferLengthVACF,sizeof(VECTOR));
                AcfOnsagerDirAvgVACF[i][j][k][l]=(REAL*)calloc(BufferLengthVACF,sizeof(REAL));
              }
            }
          }

          AccumulatedAcfVACF[i]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
          AccumulatedAcfDirAvgVACF[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
          AccumulatedAcfOnsagerVACF[i]=(VECTOR***)calloc(NumberOfComponents,sizeof(VECTOR**));
          AccumulatedAcfOnsagerDirAvgVACF[i]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));

          for(j=0;j<NumberOfComponents;j++)
          {
            AccumulatedAcfVACF[i][j]=(VECTOR*)calloc(BufferLengthVACF,sizeof(VECTOR));
            AccumulatedAcfDirAvgVACF[i][j]=(REAL*)calloc(BufferLengthVACF,sizeof(REAL));

            AccumulatedAcfOnsagerVACF[i][j]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
            AccumulatedAcfOnsagerDirAvgVACF[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
            for(k=0;k<NumberOfComponents;k++)
            {
              AccumulatedAcfOnsagerVACF[i][j][k]=(VECTOR*)calloc(BufferLengthVACF,sizeof(VECTOR));
              AccumulatedAcfOnsagerDirAvgVACF[i][j][k]=(REAL*)calloc(BufferLengthVACF,sizeof(REAL));
            }
          }
        }
      }
      break;
    case INITIALIZE:
      // trick to space the origins evenly (see for example Rapaport 2004)
      if(ComputeVACF[CurrentSystem])
        for(CurrentBuffer=0;CurrentBuffer<NumberOfBuffersVACF;CurrentBuffer++)
          CountVACF[CurrentSystem][CurrentBuffer]=(int)(-CurrentBuffer*BufferLengthVACF/NumberOfBuffersVACF);
      break;
    case SAMPLE:
      if(!(ComputeVACF[CurrentSystem]&&CurrentCycle%SampleVACFEvery[CurrentSystem]==0)) return;

      for(CurrentBuffer=0;CurrentBuffer<NumberOfBuffersVACF;CurrentBuffer++)
      {
        if(CountVACF[CurrentSystem][CurrentBuffer]==0)
        {
          for(i=0;i<NumberOfComponents;i++)
            sum_vel[i].x=sum_vel[i].y=sum_vel[i].z=0.0;
          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
          {
            type=Adsorbates[CurrentSystem][i].Type;
            vel=GetAdsorbateCenterOfMassVelocity(i);
            OriginVACF[CurrentSystem][CurrentBuffer][i]=vel;

            sum_vel[type].x+=vel.x;
            sum_vel[type].y+=vel.y;
            sum_vel[type].z+=vel.z;
          }
          for(i=0;i<NumberOfComponents;i++)
            OriginOnsagerVACF[CurrentSystem][CurrentBuffer][i]=sum_vel[i];
        }
        if(CountVACF[CurrentSystem][CurrentBuffer]>=0)
        {
          index=(int)CountVACF[CurrentSystem][CurrentBuffer];

          for(i=0;i<NumberOfComponents;i++)
          {
            AcfVACF[CurrentSystem][CurrentBuffer][i][index].x=0.0;
            AcfVACF[CurrentSystem][CurrentBuffer][i][index].y=0.0;
            AcfVACF[CurrentSystem][CurrentBuffer][i][index].z=0.0;
            AcfDirAvgVACF[CurrentSystem][CurrentBuffer][i][index]=0.0;
          }

          for(i=0;i<NumberOfComponents;i++)
            sum_vel[i].x=sum_vel[i].y=sum_vel[i].z=0.0;
          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
          {
            type=Adsorbates[CurrentSystem][i].Type;
            vel=GetAdsorbateCenterOfMassVelocity(i);
            AcfVACF[CurrentSystem][CurrentBuffer][type][index].x+=vel.x*OriginVACF[CurrentSystem][CurrentBuffer][i].x;
            AcfVACF[CurrentSystem][CurrentBuffer][type][index].y+=vel.y*OriginVACF[CurrentSystem][CurrentBuffer][i].y;
            AcfVACF[CurrentSystem][CurrentBuffer][type][index].z+=vel.z*OriginVACF[CurrentSystem][CurrentBuffer][i].z;
            AcfDirAvgVACF[CurrentSystem][CurrentBuffer][type][index]+=vel.x*OriginVACF[CurrentSystem][CurrentBuffer][i].x+
                                                                      vel.y*OriginVACF[CurrentSystem][CurrentBuffer][i].y+
                                                                      vel.z*OriginVACF[CurrentSystem][CurrentBuffer][i].z;
            sum_vel[type].x+=vel.x;
            sum_vel[type].y+=vel.y;
            sum_vel[type].z+=vel.z;
          }
          for(i=0;i<NumberOfComponents;i++)
          {
            for(j=0;j<NumberOfComponents;j++)
            {
              AcfOnsagerVACF[CurrentSystem][CurrentBuffer][i][j][index].x=sum_vel[i].x*OriginOnsagerVACF[CurrentSystem][CurrentBuffer][j].x;
              AcfOnsagerVACF[CurrentSystem][CurrentBuffer][i][j][index].y=sum_vel[i].y*OriginOnsagerVACF[CurrentSystem][CurrentBuffer][j].y;
              AcfOnsagerVACF[CurrentSystem][CurrentBuffer][i][j][index].z=sum_vel[i].z*OriginOnsagerVACF[CurrentSystem][CurrentBuffer][j].z;
              AcfOnsagerDirAvgVACF[CurrentSystem][CurrentBuffer][i][j][index]=sum_vel[i].x*OriginOnsagerVACF[CurrentSystem][CurrentBuffer][j].x+
                                                                              sum_vel[i].y*OriginOnsagerVACF[CurrentSystem][CurrentBuffer][j].y+
                                                                              sum_vel[i].z*OriginOnsagerVACF[CurrentSystem][CurrentBuffer][j].z;
            }
          }
        }
        CountVACF[CurrentSystem][CurrentBuffer]++;
      }

      // accumulate the vacf
      for(CurrentBuffer=0;CurrentBuffer<NumberOfBuffersVACF;CurrentBuffer++)
      {
        if(CountVACF[CurrentSystem][CurrentBuffer]==BufferLengthVACF)
        {
          for(k=0;k<NumberOfComponents;k++)
          {
            for(i=0;i<BufferLengthVACF;i++)
            {
              AccumulatedAcfVACF[CurrentSystem][k][i].x+=AcfVACF[CurrentSystem][CurrentBuffer][k][i].x;
              AccumulatedAcfVACF[CurrentSystem][k][i].y+=AcfVACF[CurrentSystem][CurrentBuffer][k][i].y;
              AccumulatedAcfVACF[CurrentSystem][k][i].z+=AcfVACF[CurrentSystem][CurrentBuffer][k][i].z;
              AccumulatedAcfDirAvgVACF[CurrentSystem][k][i]+=AcfDirAvgVACF[CurrentSystem][CurrentBuffer][k][i];
            }
          }
          for(k=0;k<NumberOfComponents;k++)
          {
            for(l=0;l<NumberOfComponents;l++)
            {
              for(i=0;i<BufferLengthVACF;i++)
              {
                AccumulatedAcfOnsagerVACF[CurrentSystem][k][l][i].x+=AcfOnsagerVACF[CurrentSystem][CurrentBuffer][k][l][i].x;
                AccumulatedAcfOnsagerVACF[CurrentSystem][k][l][i].y+=AcfOnsagerVACF[CurrentSystem][CurrentBuffer][k][l][i].y;
                AccumulatedAcfOnsagerVACF[CurrentSystem][k][l][i].z+=AcfOnsagerVACF[CurrentSystem][CurrentBuffer][k][l][i].z;
                AccumulatedAcfOnsagerDirAvgVACF[CurrentSystem][k][l][i]+=AcfOnsagerDirAvgVACF[CurrentSystem][CurrentBuffer][k][l][i];
              }
            }
          }
          CountVACF[CurrentSystem][CurrentBuffer]=0;
          CountAccumulatedVACF[CurrentSystem]++;
        }
      }

      break;
    case PRINT:
      if(!(ComputeVACF[CurrentSystem]&&(CurrentCycle%WriteVACFEvery[CurrentSystem]==0))) return;

      mkdir("VACF",S_IRWXU);

      sprintf(buffer,"VACF/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(k=0;k<NumberOfComponents;k++)
      {
        fac=1.0/(Components[k].NumberOfMolecules[CurrentSystem]*CountAccumulatedVACF[CurrentSystem]);

        sprintf(buffer,"VACF/System_%d/vacf_self_%s_%d%s.dat",
                CurrentSystem,
                Components[k].Name,
                k,
                FileNameAppend);

        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# column 1: time [ps]\n");
        fprintf(FilePtr,"# column 2: vacf xyz\n");
        fprintf(FilePtr,"# column 3: vacf x\n");
        fprintf(FilePtr,"# column 4: vacf y\n");
        fprintf(FilePtr,"# column 5: vacf z\n");

        D=DIFFUSION_CONVERSION_FACTOR*fac*IntegrateGeneralizedSimpsonConventional(&AccumulatedAcfDirAvgVACF[CurrentSystem][k][0],DeltaT*SampleVACFEvery[CurrentSystem],BufferLengthVACF)/3.0;
        fprintf(FilePtr,"# D: %g [m^2/s]\n",D);

        D=DIFFUSION_CONVERSION_FACTOR*fac*IntegrateGeneralizedSimpsonConventionalVectorX(&AccumulatedAcfVACF[CurrentSystem][k][0],DeltaT*SampleVACFEvery[CurrentSystem],BufferLengthVACF);
        fprintf(FilePtr,"# D(x): %g [m^2/s]\n",D);

        D=DIFFUSION_CONVERSION_FACTOR*fac*IntegrateGeneralizedSimpsonConventionalVectorY(&AccumulatedAcfVACF[CurrentSystem][k][0],DeltaT*SampleVACFEvery[CurrentSystem],BufferLengthVACF);
        fprintf(FilePtr,"# D(y): %g [m^2/s]\n",D);

        D=DIFFUSION_CONVERSION_FACTOR*fac*IntegrateGeneralizedSimpsonConventionalVectorZ(&AccumulatedAcfVACF[CurrentSystem][k][0],DeltaT*SampleVACFEvery[CurrentSystem],BufferLengthVACF);
        fprintf(FilePtr,"# D(z): %g [m^2/s]\n",D);

        for(i=0;i<BufferLengthVACF;i++)
          fprintf(FilePtr,"%g %g %g %g %g\n",i*SampleVACFEvery[CurrentSystem]*DeltaT,
            fac*AccumulatedAcfDirAvgVACF[CurrentSystem][k][i],
            fac*AccumulatedAcfVACF[CurrentSystem][k][i].x,
            fac*AccumulatedAcfVACF[CurrentSystem][k][i].y,
            fac*AccumulatedAcfVACF[CurrentSystem][k][i].z);
        fclose(FilePtr);
      }

      for(k=0;k<NumberOfComponents;k++)
      {
        fac=1.0/(Components[k].NumberOfMolecules[CurrentSystem]*CountAccumulatedVACF[CurrentSystem]);
        for(l=0;l<NumberOfComponents;l++)
        {
          sprintf(buffer,"VACF/System_%d/vacf_onsager_%s_%d_%s_%d%s.dat",
                  CurrentSystem,Components[k].Name,k,Components[l].Name,l,FileNameAppend);

          FilePtr=fopen(buffer,"w");
          fprintf(FilePtr,"# column 1: time [ps]\n");
          fprintf(FilePtr,"# column 2: vacf xyz\n");
          fprintf(FilePtr,"# column 3: vacf x\n");
          fprintf(FilePtr,"# column 4: vacf y\n");
          fprintf(FilePtr,"# column 5: vacf z\n");

          D=DIFFUSION_CONVERSION_FACTOR*fac*IntegrateGeneralizedSimpsonConventional(&AccumulatedAcfOnsagerDirAvgVACF[CurrentSystem][k][l][0],
                  DeltaT*SampleVACFEvery[CurrentSystem],BufferLengthVACF)/3.0;
          fprintf(FilePtr,"# D: %g [m^2/s]\n",D);

          D=DIFFUSION_CONVERSION_FACTOR*fac*IntegrateGeneralizedSimpsonConventionalVectorX(&AccumulatedAcfOnsagerVACF[CurrentSystem][k][l][0],
                  DeltaT*SampleVACFEvery[CurrentSystem],BufferLengthVACF);
          fprintf(FilePtr,"# D(x): %g [m^2/s]\n",D);

          D=DIFFUSION_CONVERSION_FACTOR*fac*IntegrateGeneralizedSimpsonConventionalVectorY(&AccumulatedAcfOnsagerVACF[CurrentSystem][k][l][0],
                  DeltaT*SampleVACFEvery[CurrentSystem],BufferLengthVACF);
          fprintf(FilePtr,"# D(y): %g [m^2/s]\n",D);

          D=DIFFUSION_CONVERSION_FACTOR*fac*IntegrateGeneralizedSimpsonConventionalVectorZ(&AccumulatedAcfOnsagerVACF[CurrentSystem][k][l][0],
                  DeltaT*SampleVACFEvery[CurrentSystem],BufferLengthVACF);
          fprintf(FilePtr,"# D(z): %g [m^2/s]\n",D);

          for(i=0;i<BufferLengthVACF;i++)
            fprintf(FilePtr,"%g %g %g %g %g\n",i*SampleVACFEvery[CurrentSystem]*DeltaT,
              fac*AccumulatedAcfOnsagerDirAvgVACF[CurrentSystem][k][l][i],
              fac*AccumulatedAcfOnsagerVACF[CurrentSystem][k][l][i].x,
              fac*AccumulatedAcfOnsagerVACF[CurrentSystem][k][l][i].y,
              fac*AccumulatedAcfOnsagerVACF[CurrentSystem][k][l][i].z);
          fclose(FilePtr);
        }
      }
      break;
    case FINALIZE:
      break;
  }
}

/*********************************************************************************************************                                        
 * Name       | SampleCOMDensityProfile3DVTKGrid                                                         *                                         
 * ----------------------------------------------------------------------------------------------------- *                                         
 * Function   | Samples the 3D histograms of center of mass position (e.g. 3D free energy).              *                                        
 * Parameters | -                                                                                        *                                          
 *********************************************************************************************************/
void SampleCOMDensityProfile3DVTKGrid(int Switch)
{
  int i,j,k,l,index,size,A,f;
  int x,y,z,Type;
  int k1,k2,k3;
  VECTOR pos;
  char buffer[256];
  FILE *FilePtr;
  REAL max_value,min_value;
  static REAL max;
  static VECTOR shift;
  static VECTOR Size;
  VECTOR C,s,t;
  static REAL_MATRIX3x3 Cell;
  static REAL_MATRIX3x3 InverseCell;
  
  REAL Mass,TotalMass;
  VECTOR com;
  
  switch(Switch)
  {
  	case ALLOCATE:
    if((DensityProfile3DVTKGridPoints.x<=0) || (DensityProfile3DVTKGridPoints.y<=0) || (DensityProfile3DVTKGridPoints.z<=0))
    {
	fprintf(stderr, "ERROR: number of gridpoint (%d, %d, %d) in CreateDensity3DVTKProfileGrid should be >0\n",
		DensityProfile3DVTKGridPoints.x, DensityProfile3DVTKGridPoints.y, DensityProfile3DVTKGridPoints.z);
		exit(-1);
	}
	
	size=DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*DensityProfile3DVTKGridPoints.z;
	
	COMDensityProfile3D=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
	for(i=0;i<NumberOfSystems;i++)
	{
		if((ComputeDensityProfile3DVTKGrid[i]))
		{
			COMDensityProfile3D[i]=(REAL**)calloc(NumberOfSystems,sizeof(REAL**));
			for(j=0;j<NumberOfComponents;j++)
				COMDensityProfile3D[i][j]=(REAL*)calloc(size,sizeof(REAL));
		}
	}
	break;
    case INITIALIZE:
      switch(DensityAveragingTypeVTK)
      {
        case VTK_UNIT_CELL:
          Cell=UnitCellBox[CurrentSystem];
          InverseCell=InverseUnitCellBox[CurrentSystem];
          break;
        case VTK_FULL_BOX:
        default:
          Cell=Box[CurrentSystem];
          InverseCell=InverseBox[CurrentSystem];
          break;
      }
		
		Size.x=Size.y=Size.z=0.0;
		shift.x=shift.y=shift.z=0.0;
		C.x=1.0;
		C.y=0.0;
		C.z=0.0;
		pos.x=Cell.ax*C.x+Cell.bx*C.y+Cell.cx*C.z;
		pos.y=Cell.ay*C.x+Cell.by*C.y+Cell.cy*C.z;
		pos.z=Cell.az*C.x+Cell.bz*C.y+Cell.cz*C.z;
		Size.x+=fabs(pos.x);
		Size.y+=fabs(pos.y);
		Size.z+=fabs(pos.z);
		if(pos.x<0.0) shift.x+=pos.x;
		if(pos.y<0.0) shift.y+=pos.y;
		if(pos.z<0.0) shift.z+=pos.z;
		
		C.x=0.0;
		C.y=1.0;
		C.z=0.0;
		pos.x=Cell.ax*C.x+Cell.bx*C.y+Cell.cx*C.z;
		pos.y=Cell.ay*C.x+Cell.by*C.y+Cell.cy*C.z;
		pos.z=Cell.az*C.x+Cell.bz*C.y+Cell.cz*C.z;
		Size.x+=fabs(pos.x);
		Size.y+=fabs(pos.y);
		Size.z+=fabs(pos.z);

		if(pos.x<0.0) shift.x+=pos.x;
		if(pos.y<0.0) shift.y+=pos.y;
		if(pos.z<0.0) shift.z+=pos.z;
		
		C.x=0.0;
		C.y=0.0;
		C.z=1.0;
		pos.x=Cell.ax*C.x+Cell.bx*C.y+Cell.cx*C.z;
		pos.y=Cell.ay*C.x+Cell.by*C.y+Cell.cy*C.z;
		pos.z=Cell.az*C.x+Cell.bz*C.y+Cell.cz*C.z;
		Size.x+=fabs(pos.x);
		Size.y+=fabs(pos.y);
		Size.z+=fabs(pos.z);
		if(pos.x<0.0) shift.x+=pos.x;
		if(pos.y<0.0) shift.y+=pos.y;
		if(pos.z<0.0) shift.z+=pos.z;

		max=MAX2(Size.x,MAX2(Size.y,Size.z));
		break;
	case SAMPLE:
		if(!ComputeDensityProfile3DVTKGrid[CurrentSystem]) return;
		
		for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
		{
		  TotalMass=0.0;
		  com.x=com.y=com.z=0.0;
		  Type=Adsorbates[CurrentSystem][i].Type;
		  for(l=0;l<Components[Type].NumberOfGroups;l++)
		  {
		      if(Components[Type].Groups[l].Rigid)
		        {
			  Mass=Components[Type].Groups[l].Mass;
			  TotalMass+=Mass;
			  com.x+=Mass*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.x;
			  com.y+=Mass*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.y;
			  com.z+=Mass*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.z;
			}
			else
			{
			  for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
			    {
			      A=Components[Type].Groups[l].Atoms[k];
			      Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
			      TotalMass+=Mass;
			      com.x+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Position.x;
			      com.y+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Position.y;
			      com.z+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Position.z;
			    }
			}
		    }
		  com.x=com.x/TotalMass;
		  com.y=com.y/TotalMass;
		  com.z=com.z/TotalMass;
		
		// if averaged over all unit cell and using the full box
		if((AverageDensityOverUnitCellsVTK==TRUE)&&(DensityAveragingTypeVTK==VTK_FULL_BOX))
			{
				s.x=InverseUnitCellBox[CurrentSystem].ax*com.x+InverseUnitCellBox[CurrentSystem].bx*com.y+InverseUnitCellBox[CurrentSystem].cx*com.z;
				s.y=InverseUnitCellBox[CurrentSystem].ay*com.x+InverseUnitCellBox[CurrentSystem].by*com.y+InverseUnitCellBox[CurrentSystem].cy*com.z;
				s.z=InverseUnitCellBox[CurrentSystem].az*com.x+InverseUnitCellBox[CurrentSystem].bz*com.y+InverseUnitCellBox[CurrentSystem].cz*com.z;
			
				t.x=s.x-(REAL)NINT(s.x);
				t.y=s.y-(REAL)NINT(s.y);
				t.z=s.z-(REAL)NINT(s.z);
				if(t.x<0.0) t.x+=1.0;
				if(t.y<0.0) t.y+=1.0;
				if(t.z<0.0) t.z+=1.0;
			
				for(k1=0;k1<NumberOfUnitCells[CurrentSystem].x;k1++)
					for(k2=0;k2<NumberOfUnitCells[CurrentSystem].y;k2++)
						for(k3=0;k3<NumberOfUnitCells[CurrentSystem].z;k3++)
						{
							s.x=(t.x+k1)/NumberOfUnitCells[CurrentSystem].x;
							s.y=(t.y+k2)/NumberOfUnitCells[CurrentSystem].y;
							s.z=(t.z+k3)/NumberOfUnitCells[CurrentSystem].z;
							pos.x=Cell.ax*s.x+Cell.bx*s.y+Cell.cx*s.z;
							pos.y=Cell.ay*s.x+Cell.by*s.y+Cell.cy*s.z;
							pos.z=Cell.az*s.x+Cell.bz*s.y+Cell.cz*s.z;
						
							x=(int)((pos.x-shift.x)*(REAL)DensityProfile3DVTKGridPoints.x/Size.x);
							y=(int)((pos.y-shift.y)*(REAL)DensityProfile3DVTKGridPoints.y/Size.y);
							z=(int)((pos.z-shift.z)*(REAL)DensityProfile3DVTKGridPoints.z/Size.z);
						
							index=x+y*DensityProfile3DVTKGridPoints.y+z*DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y;
							if((index>0)&&(index<DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*DensityProfile3DVTKGridPoints.z))
								COMDensityProfile3D[CurrentSystem][Type][index]+=1.0;
						}
				}
		else // if not averaged over all unit-cells or just a single unit cell
		{
				s.x=InverseCell.ax*com.x+InverseCell.bx*com.y+InverseCell.cx*com.z;
				s.y=InverseCell.ay*com.x+InverseCell.by*com.y+InverseCell.cy*com.z;
				s.z=InverseCell.az*com.x+InverseCell.bz*com.y+InverseCell.cz*com.z;
				
				t.x=s.x-(REAL)NINT(s.x);
				t.y=s.y-(REAL)NINT(s.y);
				t.z=s.z-(REAL)NINT(s.z);
				if(t.x<0.0) t.x+=1.0;
				if(t.y<0.0) t.y+=1.0;
				if(t.z<0.0) t.z+=1.0;
				pos.x=Cell.ax*t.x+Cell.bx*t.y+Cell.cx*t.z;
				pos.y=Cell.ay*t.x+Cell.by*t.y+Cell.cy*t.z;
				pos.z=Cell.az*t.x+Cell.bz*t.y+Cell.cz*t.z;
				
				x=(int)((pos.x-shift.x)*(REAL)DensityProfile3DVTKGridPoints.x/Size.x);
				y=(int)((pos.y-shift.y)*(REAL)DensityProfile3DVTKGridPoints.y/Size.y);
				z=(int)((pos.z-shift.z)*(REAL)DensityProfile3DVTKGridPoints.z/Size.z);
				
				index=x+y*DensityProfile3DVTKGridPoints.y+z*DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y;
				if((index>0)&&(index<DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*DensityProfile3DVTKGridPoints.z))
					COMDensityProfile3D[CurrentSystem][Type][index]+=1.0;
			}
		}
		break;
	case PRINT:
      if((!ComputeDensityProfile3DVTKGrid[CurrentSystem])||(CurrentCycle%WriteDensityProfile3DVTKGridEvery[CurrentSystem]!=0)) return;
     // make the output directory
      mkdir("VTK",S_IRWXU);

      // make the system directory in the output directory
      sprintf(buffer,"VTK/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(i=0;i<NumberOfComponents;i++)
      {
        sprintf(buffer,"VTK/System_%d/COMDensityProfile_%s%s.vtk",CurrentSystem,Components[i].Name,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
        fprintf(FilePtr,"Energy zeolite: %s (%lf K)\n",Framework[CurrentSystem].Name[0],
          (double)therm_baro_stats.ExternalTemperature[CurrentSystem]);
        fprintf(FilePtr,"ASCII\n");
        fprintf(FilePtr,"DATASET STRUCTURED_POINTS\n");
        fprintf(FilePtr,"DIMENSIONS %d %d %d\n",DensityProfile3DVTKGridPoints.x,DensityProfile3DVTKGridPoints.y,DensityProfile3DVTKGridPoints.z);

        //fprintf(FilePtr,"ASPECT_RATIO %f %f %f\n",(double)(Size.x/(max*NumberOfUnitCells[CurrentSystem].x)),(double)(Size.y/(max*NumberOfUnitCells[CurrentSystem].y)),(double)(Size.z/(max*NumberOfUnitCells[CurrentSystem].z)));
        fprintf(FilePtr,"ASPECT_RATIO %f %f %f\n",(double)(Size.x/max),(double)(Size.y/max),(double)(Size.z/max));
        fprintf(FilePtr,"ORIGIN 0.0 0.0 0.0\n");
        fprintf(FilePtr,"\n");
        fprintf(FilePtr,"POINT_DATA %d\n",DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*DensityProfile3DVTKGridPoints.z);
        fprintf(FilePtr,"SCALARS scalars unsigned_short\n");
        fprintf(FilePtr,"LOOKUP_TABLE default\n");

        max_value=0.0;
        min_value=65536.0;
        for(z=0;z<DensityProfile3DVTKGridPoints.z;z++)
          for(y=0;y<DensityProfile3DVTKGridPoints.y;y++)
            for(x=0;x<DensityProfile3DVTKGridPoints.x;x++)
            {
              index=x+y*DensityProfile3DVTKGridPoints.y+z*DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y;
              if((COMDensityProfile3D[CurrentSystem][i][index]>max_value)&&(COMDensityProfile3D[CurrentSystem][i][index]>0.0))
                max_value=COMDensityProfile3D[CurrentSystem][i][index];
              if((COMDensityProfile3D[CurrentSystem][i][index]<min_value)&&(COMDensityProfile3D[CurrentSystem][i][index]>=0.0))
                min_value=COMDensityProfile3D[CurrentSystem][i][index];
            }
        for(z=0;z<DensityProfile3DVTKGridPoints.z;z++)
          for(y=0;y<DensityProfile3DVTKGridPoints.y;y++)
            for(x=0;x<DensityProfile3DVTKGridPoints.x;x++)
            {
              index=x+y*DensityProfile3DVTKGridPoints.y+z*DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y;
              if(COMDensityProfile3D[CurrentSystem][i][index]>0.0)
                fprintf(FilePtr,"%d\n",(int)((COMDensityProfile3D[CurrentSystem][i][index]-min_value)*65535.0/(fabs(max_value-min_value))));
              else
                fprintf(FilePtr,"%d\n",0);
            }
            fclose(FilePtr);
      }
      break;
	case FINALIZE:
		for(i=0;i<NumberOfSystems;i++)
		{
				if(ComputeDensityProfile3DVTKGrid[i])
					for(j=0;j<NumberOfComponents;j++)
						free(COMDensityProfile3D[i][j]);
		}
		break;
	}
}

/*********************************************************************************************************
 * Name       | SampleDensityProfile3DVTKGrid                                                            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the 3D histograms of position (e.g. 3D free energy).                             *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleDensityProfile3DVTKGrid(int Switch)
{
  int i,j,index,size;
  int x,y,z,Type;
  int k1,k2,k3;
  VECTOR pos;
  char buffer[256];
  FILE *FilePtr;
  REAL max_value,min_value;
  static REAL max;
  static VECTOR shift;
  static VECTOR Size;
  VECTOR C,s,t;
  static REAL_MATRIX3x3 Cell;
  static REAL_MATRIX3x3 InverseCell;

  switch(Switch)
  {
    case ALLOCATE:
      if((DensityProfile3DVTKGridPoints.x<=0)||(DensityProfile3DVTKGridPoints.y<=0)||(DensityProfile3DVTKGridPoints.z<=0))
      {
        fprintf(stderr, "ERROR: number of gridpoint (%d, %d, %d) in CreateDensity3DVTKProfileGrid should be >0\n",
          DensityProfile3DVTKGridPoints.x,DensityProfile3DVTKGridPoints.y,DensityProfile3DVTKGridPoints.z);
        exit(-1);
      }

      size=DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*DensityProfile3DVTKGridPoints.z;

      DensityProfile3D=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeDensityProfile3DVTKGrid[i])
        {
          DensityProfile3D[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
          for(j=0;j<NumberOfComponents;j++)
            DensityProfile3D[i][j]=(REAL*)calloc(size,sizeof(REAL));
        }
      }
      break;
    case INITIALIZE:
      switch(DensityAveragingTypeVTK)
      {
        case VTK_UNIT_CELL:
          Cell=UnitCellBox[CurrentSystem];
          InverseCell=InverseUnitCellBox[CurrentSystem];
          break;
        case VTK_FULL_BOX:
        default:
          Cell=Box[CurrentSystem];
          InverseCell=InverseBox[CurrentSystem];
          break;
      }

      Size.x=Size.y=Size.z=0.0;
      shift.x=shift.y=shift.z=0.0;
      C.x=1.0;
      C.y=0.0;
      C.z=0.0;
      pos.x=Cell.ax*C.x+Cell.bx*C.y+Cell.cx*C.z;
      pos.y=Cell.ay*C.x+Cell.by*C.y+Cell.cy*C.z;
      pos.z=Cell.az*C.x+Cell.bz*C.y+Cell.cz*C.z;
      Size.x+=fabs(pos.x);
      Size.y+=fabs(pos.y);
      Size.z+=fabs(pos.z);
      if(pos.x<0.0) shift.x+=pos.x;
      if(pos.y<0.0) shift.y+=pos.y;
      if(pos.z<0.0) shift.z+=pos.z;

      C.x=0.0;
      C.y=1.0;
      C.z=0.0;
      pos.x=Cell.ax*C.x+Cell.bx*C.y+Cell.cx*C.z;
      pos.y=Cell.ay*C.x+Cell.by*C.y+Cell.cy*C.z;
      pos.z=Cell.az*C.x+Cell.bz*C.y+Cell.cz*C.z;
      Size.x+=fabs(pos.x);
      Size.y+=fabs(pos.y);
      Size.z+=fabs(pos.z);
      if(pos.x<0.0) shift.x+=pos.x;
      if(pos.y<0.0) shift.y+=pos.y;
      if(pos.z<0.0) shift.z+=pos.z;

      C.x=0.0;
      C.y=0.0;
      C.z=1.0;
      pos.x=Cell.ax*C.x+Cell.bx*C.y+Cell.cx*C.z;
      pos.y=Cell.ay*C.x+Cell.by*C.y+Cell.cy*C.z;
      pos.z=Cell.az*C.x+Cell.bz*C.y+Cell.cz*C.z;
      Size.x+=fabs(pos.x);
      Size.y+=fabs(pos.y);
      Size.z+=fabs(pos.z);
      if(pos.x<0.0) shift.x+=pos.x;
      if(pos.y<0.0) shift.y+=pos.y;
      if(pos.z<0.0) shift.z+=pos.z;

      max=MAX2(Size.x,MAX2(Size.y,Size.z));
      break;
    case SAMPLE:
      if(!ComputeDensityProfile3DVTKGrid[CurrentSystem]) return;

      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        Type=Adsorbates[CurrentSystem][i].Type;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          pos.x=Adsorbates[CurrentSystem][i].Atoms[j].Position.x;
          pos.y=Adsorbates[CurrentSystem][i].Atoms[j].Position.y;
          pos.z=Adsorbates[CurrentSystem][i].Atoms[j].Position.z;

          // if averaged over all unit cell and using the full
          if((AverageDensityOverUnitCellsVTK==TRUE)&&(DensityAveragingTypeVTK==VTK_FULL_BOX))
          {
            s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
            s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
            s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;

            t.x=s.x-(REAL)NINT(s.x);
            t.y=s.y-(REAL)NINT(s.y);
            t.z=s.z-(REAL)NINT(s.z);
            if(t.x<0.0) t.x+=1.0;
            if(t.y<0.0) t.y+=1.0;
            if(t.z<0.0) t.z+=1.0;

            for(k1=0;k1<NumberOfUnitCells[CurrentSystem].x;k1++)
              for(k2=0;k2<NumberOfUnitCells[CurrentSystem].y;k2++)
                for(k3=0;k3<NumberOfUnitCells[CurrentSystem].z;k3++)
                {
                  s.x=(t.x+k1)/NumberOfUnitCells[CurrentSystem].x;
                  s.y=(t.y+k2)/NumberOfUnitCells[CurrentSystem].y;
                  s.z=(t.z+k3)/NumberOfUnitCells[CurrentSystem].z;
                  pos.x=Cell.ax*s.x+Cell.bx*s.y+Cell.cx*s.z;
                  pos.y=Cell.ay*s.x+Cell.by*s.y+Cell.cy*s.z;
                  pos.z=Cell.az*s.x+Cell.bz*s.y+Cell.cz*s.z;

                  x=(int)((pos.x-shift.x)*(REAL)DensityProfile3DVTKGridPoints.x/Size.x);
                  y=(int)((pos.y-shift.y)*(REAL)DensityProfile3DVTKGridPoints.y/Size.y);
                  z=(int)((pos.z-shift.z)*(REAL)DensityProfile3DVTKGridPoints.z/Size.z);

                  index=x+y*DensityProfile3DVTKGridPoints.y+z*DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y;
                  if((index>0)&&(index<DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*DensityProfile3DVTKGridPoints.z))
                    DensityProfile3D[CurrentSystem][Type][index]+=1.0;
                }
          }
          else // if not averaged over all unit-cells or just a single unit cell
          {
            s.x=InverseCell.ax*pos.x+InverseCell.bx*pos.y+InverseCell.cx*pos.z;
            s.y=InverseCell.ay*pos.x+InverseCell.by*pos.y+InverseCell.cy*pos.z;
            s.z=InverseCell.az*pos.x+InverseCell.bz*pos.y+InverseCell.cz*pos.z;

            t.x=s.x-(REAL)NINT(s.x);
            t.y=s.y-(REAL)NINT(s.y);
            t.z=s.z-(REAL)NINT(s.z);
            if(t.x<0.0) t.x+=1.0;
            if(t.y<0.0) t.y+=1.0;
            if(t.z<0.0) t.z+=1.0;
            pos.x=Cell.ax*t.x+Cell.bx*t.y+Cell.cx*t.z;
            pos.y=Cell.ay*t.x+Cell.by*t.y+Cell.cy*t.z;
            pos.z=Cell.az*t.x+Cell.bz*t.y+Cell.cz*t.z;

            x=(int)((pos.x-shift.x)*(REAL)DensityProfile3DVTKGridPoints.x/Size.x);
            y=(int)((pos.y-shift.y)*(REAL)DensityProfile3DVTKGridPoints.y/Size.y);
            z=(int)((pos.z-shift.z)*(REAL)DensityProfile3DVTKGridPoints.z/Size.z);

            index=x+y*DensityProfile3DVTKGridPoints.y+z*DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y;
            if((index>0)&&(index<DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*DensityProfile3DVTKGridPoints.z))
              DensityProfile3D[CurrentSystem][Type][index]+=1.0;
          }
        }
      }
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        Type=Cations[CurrentSystem][i].Type;
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          pos.x=Cations[CurrentSystem][i].Atoms[j].Position.x;
          pos.y=Cations[CurrentSystem][i].Atoms[j].Position.y;
          pos.z=Cations[CurrentSystem][i].Atoms[j].Position.z;

          // if averaged over all unit cell and using the full
          if((AverageDensityOverUnitCellsVTK==TRUE)&&(DensityAveragingTypeVTK==VTK_FULL_BOX))
          {
            s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
            s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
            s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;

            t.x=s.x-(REAL)NINT(s.x);
            t.y=s.y-(REAL)NINT(s.y);
            t.z=s.z-(REAL)NINT(s.z);
            if(t.x<0.0) t.x+=1.0;
            if(t.y<0.0) t.y+=1.0;
            if(t.z<0.0) t.z+=1.0;

            for(k1=0;k1<NumberOfUnitCells[CurrentSystem].x;k1++)
              for(k2=0;k2<NumberOfUnitCells[CurrentSystem].y;k2++)
                for(k3=0;k3<NumberOfUnitCells[CurrentSystem].z;k3++)
                {
                  s.x=(t.x+k1)/NumberOfUnitCells[CurrentSystem].x;
                  s.y=(t.y+k2)/NumberOfUnitCells[CurrentSystem].y;
                  s.z=(t.z+k3)/NumberOfUnitCells[CurrentSystem].z;
                  pos.x=Cell.ax*s.x+Cell.bx*s.y+Cell.cx*s.z;
                  pos.y=Cell.ay*s.x+Cell.by*s.y+Cell.cy*s.z;
                  pos.z=Cell.az*s.x+Cell.bz*s.y+Cell.cz*s.z;

                  x=(int)((pos.x-shift.x)*(REAL)DensityProfile3DVTKGridPoints.x/Size.x);
                  y=(int)((pos.y-shift.y)*(REAL)DensityProfile3DVTKGridPoints.y/Size.y);
                  z=(int)((pos.z-shift.z)*(REAL)DensityProfile3DVTKGridPoints.z/Size.z);

                  index=x+y*DensityProfile3DVTKGridPoints.y+z*DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y;
                  if((index>0)&&(index<DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*DensityProfile3DVTKGridPoints.z))
                    DensityProfile3D[CurrentSystem][Type][index]+=1.0;
                }
          }
          else // if not averaged over all unit-cells or just a single unit cell
          {
            s.x=InverseCell.ax*pos.x+InverseCell.bx*pos.y+InverseCell.cx*pos.z;
            s.y=InverseCell.ay*pos.x+InverseCell.by*pos.y+InverseCell.cy*pos.z;
            s.z=InverseCell.az*pos.x+InverseCell.bz*pos.y+InverseCell.cz*pos.z;

            t.x=s.x-(REAL)NINT(s.x);
            t.y=s.y-(REAL)NINT(s.y);
            t.z=s.z-(REAL)NINT(s.z);
            if(t.x<0.0) t.x+=1.0;
            if(t.y<0.0) t.y+=1.0;
            if(t.z<0.0) t.z+=1.0;
            pos.x=Cell.ax*t.x+Cell.bx*t.y+Cell.cx*t.z;
            pos.y=Cell.ay*t.x+Cell.by*t.y+Cell.cy*t.z;
            pos.z=Cell.az*t.x+Cell.bz*t.y+Cell.cz*t.z;

            x=(int)((pos.x-shift.x)*(REAL)DensityProfile3DVTKGridPoints.x/Size.x);
            y=(int)((pos.y-shift.y)*(REAL)DensityProfile3DVTKGridPoints.y/Size.y);
            z=(int)((pos.z-shift.z)*(REAL)DensityProfile3DVTKGridPoints.z/Size.z);

            index=x+y*DensityProfile3DVTKGridPoints.y+z*DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y;
            if((index>0)&&(index<DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*DensityProfile3DVTKGridPoints.z))
              DensityProfile3D[CurrentSystem][Type][index]+=1.0;
          }
        }
      }
      break;
    case PRINT:
      if((!ComputeDensityProfile3DVTKGrid[CurrentSystem])||(CurrentCycle%WriteDensityProfile3DVTKGridEvery[CurrentSystem]!=0)) return;
     // make the output directory
      mkdir("VTK",S_IRWXU);

      // make the system directory in the output directory
      sprintf(buffer,"VTK/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      for(i=0;i<NumberOfComponents;i++)
      {
        sprintf(buffer,"VTK/System_%d/DensityProfile_%s%s.vtk",CurrentSystem,Components[i].Name,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
        fprintf(FilePtr,"Energy zeolite: %s (%lf K)\n",Framework[CurrentSystem].Name[0],
          (double)therm_baro_stats.ExternalTemperature[CurrentSystem]);
        fprintf(FilePtr,"ASCII\n");
        fprintf(FilePtr,"DATASET STRUCTURED_POINTS\n");
        fprintf(FilePtr,"DIMENSIONS %d %d %d\n",DensityProfile3DVTKGridPoints.x,DensityProfile3DVTKGridPoints.y,DensityProfile3DVTKGridPoints.z);

        //fprintf(FilePtr,"ASPECT_RATIO %f %f %f\n",(double)(Size.x/(max*NumberOfUnitCells[CurrentSystem].x)),(double)(Size.y/(max*NumberOfUnitCells[CurrentSystem].y)),(double)(Size.z/(max*NumberOfUnitCells[CurrentSystem].z)));
        fprintf(FilePtr,"ASPECT_RATIO %f %f %f\n",(double)(Size.x/max),(double)(Size.y/max),(double)(Size.z/max));
        fprintf(FilePtr,"ORIGIN 0.0 0.0 0.0\n");
        fprintf(FilePtr,"\n");
        fprintf(FilePtr,"POINT_DATA %d\n",DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*DensityProfile3DVTKGridPoints.z);
        fprintf(FilePtr,"SCALARS scalars unsigned_short\n");
        fprintf(FilePtr,"LOOKUP_TABLE default\n");

        max_value=0.0;
        min_value=65536.0;
        for(z=0;z<DensityProfile3DVTKGridPoints.z;z++)
          for(y=0;y<DensityProfile3DVTKGridPoints.y;y++)
            for(x=0;x<DensityProfile3DVTKGridPoints.x;x++)
            {
              index=x+y*DensityProfile3DVTKGridPoints.y+z*DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y;
              if((DensityProfile3D[CurrentSystem][i][index]>max_value)&&(DensityProfile3D[CurrentSystem][i][index]>0.0))
                max_value=DensityProfile3D[CurrentSystem][i][index];
              if((DensityProfile3D[CurrentSystem][i][index]<min_value)&&(DensityProfile3D[CurrentSystem][i][index]>=0.0))
                min_value=DensityProfile3D[CurrentSystem][i][index];
            }
        for(z=0;z<DensityProfile3DVTKGridPoints.z;z++)
          for(y=0;y<DensityProfile3DVTKGridPoints.y;y++)
            for(x=0;x<DensityProfile3DVTKGridPoints.x;x++)
            {
              index=x+y*DensityProfile3DVTKGridPoints.y+z*DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y;
              if(DensityProfile3D[CurrentSystem][i][index]>0.0)
                fprintf(FilePtr,"%d\n",(int)((DensityProfile3D[CurrentSystem][i][index]-min_value)*65535.0/(fabs(max_value-min_value))));
              else
                fprintf(FilePtr,"%d\n",0);
            }
            fclose(FilePtr);
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(ComputeDensityProfile3DVTKGrid[i])
          for(j=0;j<NumberOfComponents;j++)
            free(DensityProfile3D[i][j]);
      }
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleCationAndAdsorptionSites                                                           *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the cation sites and adsorption sites.                                           *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleCationAndAdsorptionSites(int Switch)
{
  int i,ion_type,type;
  int closest,nr_sites,starting_bead;
  VECTOR pos,dr,posB;
  REAL distance;
  REAL count,count2;
  REAL AverageDistance;
  FILE *FilePtr;
  char buffer[1024];

  switch(Switch)
  {
    case ALLOCATE:
      break;
    case INITIALIZE:
      break;
    case SAMPLE:
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        type=Cations[CurrentSystem][i].Type;
        starting_bead=Components[type].StartingBead;

        // determine position (using the starting bead)
        pos=Cations[CurrentSystem][i].Atoms[starting_bead].Position;

        // look for the closest crystallographic adsorption site
        ClosestCrystallographicPosition2(pos,&closest,&distance);

        // if the distance is less than a certain cutoff consider them 'classified'
        if(distance<CutOffIons)
        {
          ion_type=Framework[CurrentSystem].Ions[closest].AssymetricType;

          // compute distance to closest crystallogarphic site
          dr.x=Framework[CurrentSystem].Ions[closest].Position.x-pos.x;
          dr.y=Framework[CurrentSystem].Ions[closest].Position.y-pos.y;
          dr.z=Framework[CurrentSystem].Ions[closest].Position.z-pos.z;
          dr=ApplyBoundaryCondition(dr);

          crystallographic_stats[CurrentSystem].Position[ion_type].x+=dr.x;
          crystallographic_stats[CurrentSystem].Position[ion_type].y+=dr.y;
          crystallographic_stats[CurrentSystem].Position[ion_type].z+=dr.z;
          crystallographic_stats[CurrentSystem].Distance[ion_type].x+=fabs(dr.x);
          crystallographic_stats[CurrentSystem].Distance[ion_type].y+=fabs(dr.y);
          crystallographic_stats[CurrentSystem].Distance[ion_type].z+=fabs(dr.z);
          crystallographic_stats[CurrentSystem].AverageDistance[ion_type]+=distance;

          crystallographic_stats[CurrentSystem].RelativeOccupation[ion_type]+=1.0;
          nr_sites=crystallographic_stats[CurrentSystem].NumberOfCationSites[ion_type];
          crystallographic_stats[CurrentSystem].Occupation[ion_type]+=1.0/nr_sites;
          crystallographic_stats[CurrentSystem].Count[ion_type]+=1.0;

          crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].ax+=dr.x;
          crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].ay+=dr.y;
          crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].az+=dr.z;

          crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].bx+=SQR(dr.x);
          crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].by+=SQR(dr.y);
          crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].bz+=SQR(dr.z);

          crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].cx+=dr.x*dr.y;
          crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].cy+=dr.x*dr.z;
          crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].cz+=dr.y*dr.z;
        }
        else
          crystallographic_stats[CurrentSystem].Unclassified+=1.0;

        crystallographic_stats[CurrentSystem].count+=1.0;
      }
      crystallographic_stats[CurrentSystem].count2+=1.0;
      break;
    case PRINT:
      if((!ComputeCationAndAdsorptionSites[CurrentSystem])||(CurrentCycle%WriteCationAndAdsorptionSitesEvery[CurrentSystem]!=0)) return;
      mkdir("CationAndAdsorptionSites",S_IRWXU);
      sprintf(buffer,"CationAndAdsorptionSites/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      sprintf(buffer,"CationAndAdsorptionSites/System_%d/crystalographic%s.data",CurrentSystem,FileNameAppend);
      FilePtr=fopen(buffer,"w");

      count=crystallographic_stats[CurrentSystem].count;
      count2=crystallographic_stats[CurrentSystem].count2;
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"Crystallographic data\n");
      fprintf(FilePtr,"=================================================\n\n");
      fprintf(FilePtr,"Average unclassified: %lf\n\n",
         (double)((REAL)NumberOfCationMolecules[CurrentSystem]*
         crystallographic_stats[CurrentSystem].Unclassified/
            (REAL)crystallographic_stats[CurrentSystem].count));
      for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricIons;i++)
      {
        count=crystallographic_stats[CurrentSystem].Count[i];
        pos=crystallographic_stats[CurrentSystem].Position[i];
        pos.x/=count;
        pos.y/=count;
        pos.z/=count;
        posB=ConvertFromABCtoXYZ(Framework[CurrentSystem].IonsAsymmetric[i].Position);
        pos.x+=posB.x;
        pos.y+=posB.y;
        pos.z+=posB.z;

        dr=crystallographic_stats[CurrentSystem].Distance[i];
        dr.x/=count;
        dr.y/=count;
        dr.z/=count;
        AverageDistance=crystallographic_stats[CurrentSystem].AverageDistance[i]/count;

        fprintf(FilePtr,"\nCation Position %d\n",i);
        fprintf(FilePtr,"--------------------------\n");

        fprintf(FilePtr,"Relative Occupation %10.6lf Ions %5.3lf Occupation %10.6lf\n",
            (double)(crystallographic_stats[CurrentSystem].RelativeOccupation[i]/
               (REAL)crystallographic_stats[CurrentSystem].count),
            (double)((REAL)NumberOfCationMolecules[CurrentSystem]*
                crystallographic_stats[CurrentSystem].RelativeOccupation[i]/
               (REAL)crystallographic_stats[CurrentSystem].count),
            (double)(crystallographic_stats[CurrentSystem].Occupation[i]/count2));

        posB=ConvertIonsToAsymetricUnitCell(ConvertFromABCtoXYZ(Framework[CurrentSystem].IonsAsymmetric[i].Position));
        pos=ConvertIonsToAsymetricUnitCell(pos);
        fprintf(FilePtr,"  [%10.6lf %10.6lf %10.6lf] -> [%10.6lf %10.6lf %10.6lf]\n",
          (double)posB.x,
          (double)posB.y,
          (double)posB.z,
          (double)(pos.x),
          (double)(pos.y),
          (double)(pos.z));
        fprintf(FilePtr,"  distances [x] %10.6f [y] %10.6f [z] %10.6f [xyz] %10.6f %10.6f count: %f\n",
          (double)dr.x,(double)dr.y,(double)dr.z,(double)(sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z))),(double)AverageDistance,(double)count);

        fprintf(FilePtr,"B-factor:\n");
        count=crystallographic_stats[CurrentSystem].Count[i];
        fprintf(FilePtr,"%10.6lf %10.6lf %10.6lf\n",
          (double)((crystallographic_stats[CurrentSystem].TemperatureFactor[i].bx/count)-
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ax/count)*
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ax/count)),
          (double)((crystallographic_stats[CurrentSystem].TemperatureFactor[i].cx/count)-
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ax/count)*
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ay/count)),
          (double)((crystallographic_stats[CurrentSystem].TemperatureFactor[i].by/count)-
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ax/count)*
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].az/count)));
        fprintf(FilePtr,"%10.6lf %10.6lf %10lf\n",
          (double)((crystallographic_stats[CurrentSystem].TemperatureFactor[i].cx/count)-
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ax/count)*
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ay/count)),
          (double)((crystallographic_stats[CurrentSystem].TemperatureFactor[i].bz/count)-
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ay/count)*
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ay/count)),
          (double)((crystallographic_stats[CurrentSystem].TemperatureFactor[i].cz/count)-
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ay/count)*
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].az/count)));
        fprintf(FilePtr,"%10.6lf %10.6lf %10lf\n",
          (double)((crystallographic_stats[CurrentSystem].TemperatureFactor[i].by/count)-
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ax/count)*
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].az/count)),
          (double)((crystallographic_stats[CurrentSystem].TemperatureFactor[i].cz/count)-
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].ay/count)*
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].az/count)),
          (double)((crystallographic_stats[CurrentSystem].TemperatureFactor[i].bz/count)-
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].az/count)*
          (crystallographic_stats[CurrentSystem].TemperatureFactor[i].az/count)));
      }
      fclose(FilePtr);
      break;
    case FINALIZE:
      break;
  }
}


/*********************************************************************************************************
 * Name       | SampleDcTSTConfigurationFiles                                                            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples initial configurations for the tranmission coefficient (dcTST).                  *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void SampleDcTSTConfigurationFiles(int Choice)
{
  int i,j,k;
  FILE *FilePtr;
  char buffer[1024];

  switch(Choice)
  {
    case ALLOCATE:
      break;
    case INITIALIZE:
      for(k=0;k<NumberOfSystems;k++)
      {
        if(WritedcTSTSnapShotsToFile[k])
        {
          mkdir("dcTST_starting_configurations",S_IRWXU);
          sprintf(buffer,"dcTST_starting_configurations/System_%d",k);
          mkdir(buffer,S_IRWXU);
          sprintf(buffer,"dcTST_starting_configurations/System_%d/configurations%s.data",k,FileNameAppend);
          FilePtr=fopen(buffer,"w");
          fclose(FilePtr);
        }
      }
      break;
    case SAMPLE:
      break;
    case PRINT:
      if(!(WritedcTSTSnapShotsToFile[CurrentSystem]&&(CurrentCycle%WritedcTSTSnapShotsEvery[CurrentSystem]==0))) return;

      sprintf(buffer,"dcTST_starting_configurations/System_%d/configurations%s.data",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"a");

      fprintf(FilePtr,"#NumberOfAdsorbateMolecules:%d\n",NumberOfAdsorbateMolecules[CurrentSystem]);
      fprintf(FilePtr,"#NumberOfCationMolecules:%d\n",NumberOfCationMolecules[CurrentSystem]);
      fprintf(FilePtr,"#NumberOfFrameworkAtoms:%d\n",Framework[CurrentSystem].TotalNumberOfAtoms);
      fprintf(FilePtr,"#Barrier: %f %f %f\n",(double)BarrierPosition[CurrentSystem].x,(double)BarrierPosition[CurrentSystem].y,(double)BarrierPosition[CurrentSystem].z);

      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
          fprintf(FilePtr,"%f %f %f 0.0 0.0 0.0\n",
            (double)Adsorbates[CurrentSystem][i].Atoms[j].Position.x,
            (double)Adsorbates[CurrentSystem][i].Atoms[j].Position.y,
            (double)Adsorbates[CurrentSystem][i].Atoms[j].Position.z);
        fprintf(FilePtr,"\n");
      }
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          fprintf(FilePtr,"%f %f %f 0.0 0.0 0.0\n",
            (double)Cations[CurrentSystem][i].Atoms[j].Position.x,
            (double)Cations[CurrentSystem][i].Atoms[j].Position.y,
            (double)Cations[CurrentSystem][i].Atoms[j].Position.z);
        }
        fprintf(FilePtr,"\n");
      }
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        for(CurrentFramework=0;CurrentFramework<Framework[CurrentFramework].NumberOfFrameworks;CurrentFramework++)
        {
          for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
          {
            fprintf(FilePtr,"%f %f %f 0.0 0.0 0.0\n",
               (double)Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.x,
               (double)Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.y,
               (double)Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.z);
          }
        }
        fprintf(FilePtr,"\n");
      }
      fprintf(FilePtr,"\n\n");
      fclose(FilePtr);
      break;
   case FINALIZE:
      break;
  }
}


/*********************************************************************************************************
 * Name       | MeasurePrincipleMomentsOfInertia                                                         *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Samples the principle moments of inertia.                                                *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void MeasurePrincipleMomentsOfInertia(void)
{
  int i,j;
  int Type,ComponentType;
  REAL TotalMass,Mass;
  REAL_MATRIX3x3 eigenvectors;
  VECTOR com,eigenvalues,dr;
  REAL_MATRIX3x3 InertiaTensor;
  VECTOR PrincipleMomentsOfInertia,InertiaVector,pos;
  REAL rotxyz,temp;

  PrincipleMomentsOfInertia.x=0.0;
  PrincipleMomentsOfInertia.y=0.0;
  PrincipleMomentsOfInertia.z=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    ComponentType=Adsorbates[CurrentSystem][i].Type;
    TotalMass=0.0;
    com.x=com.y=com.z=0.0;

    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      Mass=PseudoAtoms[Type].Mass;
      com.x+=Mass*Adsorbates[CurrentSystem][i].Atoms[j].Position.x;
      com.y+=Mass*Adsorbates[CurrentSystem][i].Atoms[j].Position.y;
      com.z+=Mass*Adsorbates[CurrentSystem][i].Atoms[j].Position.z;
      TotalMass+=Mass;
    }
    com.x/=TotalMass;
    com.y/=TotalMass;
    com.z/=TotalMass;


    InertiaTensor.ax=InertiaTensor.bx=InertiaTensor.cx=0.0;
    InertiaTensor.ay=InertiaTensor.by=InertiaTensor.cy=0.0;
    InertiaTensor.az=InertiaTensor.bz=InertiaTensor.cz=0.0;

    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      Mass=PseudoAtoms[Type].Mass;

      dr.x=Adsorbates[CurrentSystem][i].Atoms[j].Position.x-com.x;
      dr.y=Adsorbates[CurrentSystem][i].Atoms[j].Position.y-com.y;
      dr.z=Adsorbates[CurrentSystem][i].Atoms[j].Position.z-com.z;

      InertiaTensor.ax+=Mass*dr.x*dr.x;
      InertiaTensor.bx+=Mass*dr.y*dr.x;
      InertiaTensor.cx+=Mass*dr.z*dr.x;
      InertiaTensor.ay+=Mass*dr.x*dr.y;
      InertiaTensor.by+=Mass*dr.y*dr.y;
      InertiaTensor.cy+=Mass*dr.z*dr.y;
      InertiaTensor.az+=Mass*dr.x*dr.z;
      InertiaTensor.bz+=Mass*dr.y*dr.z;
      InertiaTensor.cz+=Mass*dr.z*dr.z;
    }

    // solve for principle moments of inertia and axes
    EigenSystem3x3(InertiaTensor,&eigenvectors,&eigenvalues);

    InertiaVector.x=0.0;
    InertiaVector.y=0.0;
    InertiaVector.z=0.0;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      Mass=PseudoAtoms[Type].Mass;

      dr.x=Adsorbates[CurrentSystem][i].Atoms[j].Position.x-com.x;
      dr.y=Adsorbates[CurrentSystem][i].Atoms[j].Position.y-com.y;
      dr.z=Adsorbates[CurrentSystem][i].Atoms[j].Position.z-com.z;

      // site positions in principal axis system
      pos.x=eigenvectors.ax*dr.x+eigenvectors.bx*dr.y+eigenvectors.cx*dr.z;
      pos.y=eigenvectors.ay*dr.x+eigenvectors.by*dr.y+eigenvectors.cy*dr.z;
      pos.z=eigenvectors.az*dr.x+eigenvectors.bz*dr.y+eigenvectors.cz*dr.z;

      if(fabs(pos.x)<1e-8) pos.x=0.0;
      if(fabs(pos.y)<1e-8) pos.y=0.0;
      if(fabs(pos.z)<1e-8) pos.z=0.0;

      // rotational inertia tensor in principal axis system
      InertiaVector.x+=Mass*(SQR(pos.y)+SQR(pos.z));
      InertiaVector.y+=Mass*(SQR(pos.x)+SQR(pos.z));
      InertiaVector.z+=Mass*(SQR(pos.x)+SQR(pos.y));
    }

    // set axis system: Ixx >= Iyy >= Izz
    rotxyz=MAX3(InertiaVector.x,InertiaVector.y,InertiaVector.z);
    if(rotxyz>=InertiaVector.x)
    {
      if(InertiaVector.y>=rotxyz)
      {
        InertiaVector.y=InertiaVector.x;
        InertiaVector.x=rotxyz;
      }
      else if(InertiaVector.z>=rotxyz)
      {
        InertiaVector.z=InertiaVector.x;
        InertiaVector.x=rotxyz;
      }
    }
    if(InertiaVector.z>InertiaVector.y)
    {
      temp=InertiaVector.z;
      InertiaVector.z=InertiaVector.y;
      InertiaVector.y=temp;
    }

    // register the principle moments of inertia per component type
    PrincipleMomentsOfInertiaAccumulated[CurrentSystem][ComponentType][Block].x+=InertiaVector.x;
    PrincipleMomentsOfInertiaAccumulated[CurrentSystem][ComponentType][Block].y+=InertiaVector.y;
    PrincipleMomentsOfInertiaAccumulated[CurrentSystem][ComponentType][Block].z+=InertiaVector.z;

    // CountMSDOrderNer needed because the amount of molecules could vary in GC-MC
    PrincipleMomentsOfInertiaCount[CurrentSystem][ComponentType][Block]+=1.0;
  }

  PrincipleMomentsOfInertia.x=0.0;
  PrincipleMomentsOfInertia.y=0.0;
  PrincipleMomentsOfInertia.z=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    ComponentType=Cations[CurrentSystem][i].Type;
    TotalMass=0.0;
    com.x=com.y=com.z=0.0;

    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][i].Atoms[j].Type;
      Mass=PseudoAtoms[Type].Mass;
      com.x+=Mass*Cations[CurrentSystem][i].Atoms[j].Position.x;
      com.y+=Mass*Cations[CurrentSystem][i].Atoms[j].Position.y;
      com.z+=Mass*Cations[CurrentSystem][i].Atoms[j].Position.z;
      TotalMass+=Mass;
    }
    com.x/=TotalMass;
    com.y/=TotalMass;
    com.z/=TotalMass;


    InertiaTensor.ax=InertiaTensor.bx=InertiaTensor.cx=0.0;
    InertiaTensor.ay=InertiaTensor.by=InertiaTensor.cy=0.0;
    InertiaTensor.az=InertiaTensor.bz=InertiaTensor.cz=0.0;

    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][i].Atoms[j].Type;
      Mass=PseudoAtoms[Type].Mass;

      dr.x=Cations[CurrentSystem][i].Atoms[j].Position.x-com.x;
      dr.y=Cations[CurrentSystem][i].Atoms[j].Position.y-com.y;
      dr.z=Cations[CurrentSystem][i].Atoms[j].Position.z-com.z;

      InertiaTensor.ax+=Mass*dr.x*dr.x;
      InertiaTensor.bx+=Mass*dr.y*dr.x;
      InertiaTensor.cx+=Mass*dr.z*dr.x;
      InertiaTensor.ay+=Mass*dr.x*dr.y;
      InertiaTensor.by+=Mass*dr.y*dr.y;
      InertiaTensor.cy+=Mass*dr.z*dr.y;
      InertiaTensor.az+=Mass*dr.x*dr.z;
      InertiaTensor.bz+=Mass*dr.y*dr.z;
      InertiaTensor.cz+=Mass*dr.z*dr.z;
    }

    // solve for principle moments of inertia and axes
    EigenSystem3x3(InertiaTensor,&eigenvectors,&eigenvalues);

    InertiaVector.x=0.0;
    InertiaVector.y=0.0;
    InertiaVector.z=0.0;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][i].Atoms[j].Type;
      Mass=PseudoAtoms[Type].Mass;

      dr.x=Cations[CurrentSystem][i].Atoms[j].Position.x-com.x;
      dr.y=Cations[CurrentSystem][i].Atoms[j].Position.y-com.y;
      dr.z=Cations[CurrentSystem][i].Atoms[j].Position.z-com.z;

      // site positions in principal axis system
      pos.x=eigenvectors.ax*dr.x+eigenvectors.bx*dr.y+eigenvectors.cx*dr.z;
      pos.y=eigenvectors.ay*dr.x+eigenvectors.by*dr.y+eigenvectors.cy*dr.z;
      pos.z=eigenvectors.az*dr.x+eigenvectors.bz*dr.y+eigenvectors.cz*dr.z;

      if(fabs(pos.x)<1e-8) pos.x=0.0;
      if(fabs(pos.y)<1e-8) pos.y=0.0;
      if(fabs(pos.z)<1e-8) pos.z=0.0;

      // rotational inertia tensor in principal axis system
      InertiaVector.x+=Mass*(SQR(pos.y)+SQR(pos.z));
      InertiaVector.y+=Mass*(SQR(pos.x)+SQR(pos.z));
      InertiaVector.z+=Mass*(SQR(pos.x)+SQR(pos.y));
    }

    // set axis system: Ixx >= Iyy >= Izz
    rotxyz=MAX3(InertiaVector.x,InertiaVector.y,InertiaVector.z);
    if(rotxyz>=InertiaVector.x)
    {
      if(InertiaVector.y>=rotxyz)
      {
        InertiaVector.y=InertiaVector.x;
        InertiaVector.x=rotxyz;
      }
      else if(InertiaVector.z>=rotxyz)
      {
        InertiaVector.z=InertiaVector.x;
        InertiaVector.x=rotxyz;
      }
    }
    if(InertiaVector.z>InertiaVector.y)
    {
      temp=InertiaVector.z;
      InertiaVector.z=InertiaVector.y;
      InertiaVector.y=temp;
    }

    // register the principle moments of inertia per component type
    PrincipleMomentsOfInertiaAccumulated[CurrentSystem][ComponentType][Block].x+=InertiaVector.x;
    PrincipleMomentsOfInertiaAccumulated[CurrentSystem][ComponentType][Block].y+=InertiaVector.y;
    PrincipleMomentsOfInertiaAccumulated[CurrentSystem][ComponentType][Block].z+=InertiaVector.z;

    // CountMSDOrderNer needed because the amount of molecules could vary in GC-MC
    PrincipleMomentsOfInertiaCount[CurrentSystem][ComponentType][Block]+=1.0;
  }
}


/*********************************************************************************************************
 * Name       | ComputeMolecularPressureTensor                                                           *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the ideal gas part of the pressure and stress, and pressure tail-correction.    *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void ComputeMolecularPressureTensor(REAL_MATRIX3x3 *Pressure_matrix, REAL* PressureIdealGas, REAL* PressureTail)
{
  // compute the ideal gas part of the pressure (related to the translational degrees of freedom, not rotational dof)
  *PressureIdealGas=(NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem])
                  *K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Volume[CurrentSystem];

  // compute the excess pressure
  CalculateMolecularExcessPressure();
  *Pressure_matrix=StrainDerivativeTensor[CurrentSystem];

  // compute the tail correction
  CalculateTailCorrection();
  *PressureTail=StrainDerivativeTailCorrection[CurrentSystem]/Volume[CurrentSystem];
}


/*********************************************************************************************************
 * Name       | WriteRestartSample                                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Write the stored data of the sampling routines in binary form to a restart-file.         *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

static int versionNumber=1;

void WriteRestartSample(FILE *FilePtr)
{
  int i,j,k,l,f1;
  REAL Check;

  fwrite(&versionNumber,sizeof(int),1,FilePtr);

  // write data for the histograms of the radial distribution function
  fwrite(ComputeRDF,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteRDFEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(RDFHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(RDFRange,NumberOfSystems,sizeof(REAL),FilePtr);
  fwrite(CountRDF,NumberOfSystems,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeRDF[i])
    {
      for(j=0;j<NumberOfPseudoAtoms;j++)
        fwrite(RadialDistributionFunctionWithFramework[i][j],RDFHistogramSize[i],sizeof(REAL),FilePtr);
      for(j=0;j<NumberOfPseudoAtoms;j++)
        for(k=0;k<NumberOfPseudoAtoms;k++)
          fwrite(RadialDistributionFunction[i][j][k],RDFHistogramSize[i],sizeof(REAL),FilePtr);
    }
  }

  // write data for the histograms of the projected lengths
  fwrite(ComputeProjectedLengths,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteProjectedLengthsEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(ProjectedLengthsHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(ProjectedLengthsRange,NumberOfSystems,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeProjectedLengths[i])
    {
      fwrite(CountProjectedLengths[i],NumberOfComponents,sizeof(REAL),FilePtr);
      fwrite(ProjectedLengthsAverage[i],NumberOfComponents,sizeof(VECTOR),FilePtr);
      for(j=0;j<NumberOfComponents;j++)
        fwrite(ProjectedLengthsDistributionFunction[i][j],ProjectedLengthsHistogramSize[i],sizeof(VECTOR),FilePtr);
    }
  }

  // write data for the histograms of the projected angles
  fwrite(ComputeProjectedAngles,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteProjectedAnglesEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(ProjectedAnglesHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(ProjectedAnglesRange,NumberOfSystems,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeProjectedAngles[i])
    {
      fwrite(CountProjectedAngles[i],NumberOfComponents,sizeof(REAL),FilePtr);
      for(j=0;j<NumberOfComponents;j++)
        fwrite(ProjectedAnglesDistributionFunction[i][j],ProjectedAnglesHistogramSize[i],sizeof(VECTOR),FilePtr);
    }
  }


  // write data for the histograms of the number of molecules
  fwrite(ComputeNumberOfMoleculesHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteNumberOfMoleculesHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(NumberOfMoleculesHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(NumberOfMoleculesRange,NumberOfSystems,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeNumberOfMoleculesHistogram[i])
    {
      for(j=0;j<NumberOfComponents+1;j++)
        fwrite(NumberOfMoleculesHistogram[i][j],NumberOfMoleculesHistogramSize[i],sizeof(REAL),FilePtr);
    }
  }

  // write data for the histograms of the position (i.e. free energy)
  fwrite(ComputePositionHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(PositionHistogramMappingType,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WritePositionHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(PositionHistogramSize,NumberOfSystems,sizeof(int),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputePositionHistogram[i])
    {
      for(j=0;j<NumberOfComponents;j++)
        for(k=0;k<Components[j].NumberOfAtoms+2;k++)
        {
          fwrite(PositionABCHistogram[i][j][k],PositionHistogramSize[i],sizeof(VECTOR),FilePtr);
          fwrite(Position2DDiagonalHistogram[i][j][k],PositionHistogramSize[i],sizeof(VECTOR),FilePtr);
          fwrite(Position3DDiagonalHistogram[i][j][k],PositionHistogramSize[i],sizeof(VECTOR4),FilePtr);
        }
    }
  }

  // write free energy profile data
  fwrite(FreeEnergyMappingType,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(ComputeFreeEnergyProfile,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(FreeEnergyHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteFreeEnergyProfileEvery,NumberOfSystems,sizeof(int),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeFreeEnergyProfile[i])
    {
      for(j=0;j<NumberOfComponents;j++)
      {
        fwrite(RosenBinSum[i][j],FreeEnergyHistogramSize[i]+1,sizeof(REAL[13]),FilePtr);
        fwrite(RosenBinCount[i][j],FreeEnergyHistogramSize[i]+1,sizeof(REAL[13]),FilePtr);
        fwrite(RosenBinSumSquared[i][j],FreeEnergyHistogramSize[i]+1,sizeof(REAL[13]),FilePtr);
      }
    }
  }

  // write the pore-size distribution (PSD)
  fwrite(ComputePSDHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WritePSDHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(PSDHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(PSDRange,NumberOfSystems,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputePSDHistogram[i])
      fwrite(PoreSizeDistributionHistogram[i],PSDHistogramSize[i],sizeof(REAL),FilePtr);
  }

  // write the end-to-end distribution
  fwrite(ComputeEndToEndDistanceHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteEndToEndDistanceHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(EndToEndHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(EndToEndRange,NumberOfSystems,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeEndToEndDistanceHistogram[i])
    {
      for(j=0;j<NumberOfComponents;j++)
        fwrite(EndToEndDistanceHistogram[i][j],EndToEndHistogramSize[i],sizeof(REAL),FilePtr);
    }
  }

  // write data for the histograms of the energy
  fwrite(ComputeEnergyHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteEnergyHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);

  fwrite(EnergyHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(EnergyHistogramLowerLimit,NumberOfSystems,sizeof(REAL),FilePtr);
  fwrite(EnergyHistogramUpperLimit,NumberOfSystems,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeEnergyHistogram[i])
    {
      for(j=0;j<4;j++)
       fwrite(EnergyHistogram[i][j],EnergyHistogramSize[i],sizeof(REAL),FilePtr);
    }
  }

  // write thermodynamic factor
  fwrite(ComputeThermoDynamicFactor,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteThermoDynamicFactorEvery,NumberOfSystems,sizeof(int),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeThermoDynamicFactor[i])
    {
      fwrite(ThermoDynamicFactorNumberOfSamples,NumberOfComponents,sizeof(REAL),FilePtr);
      fwrite(ThermoDynamicFactorNumberOfMolecules[i],NumberOfComponents,sizeof(REAL),FilePtr);
      for(k=0;k<NumberOfComponents;k++)
        fwrite(ThermoDynamicFactorNumberOfMoleculesCrossTerm[i][k],NumberOfComponents,sizeof(REAL),FilePtr);
    }
  }

  // write framework spacing data
  fwrite(ComputeFrameworkSpacingHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteFrameworkSpacingHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(FrameworkSpacingHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(FrameworkSpacingRange,NumberOfSystems,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeFrameworkSpacingHistogram[i])
    {
      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      {
        fwrite(OriginalFrameworkShiftDirAvg[i][j],Framework[i].NumberOfFrameworks,sizeof(REAL),FilePtr);
        fwrite(OriginalFrameworkShift[i][j],Framework[i].NumberOfFrameworks,sizeof(VECTOR),FilePtr);
        for(k=0;k<Framework[i].NumberOfFrameworks;k++)
          for(l=0;l<4;l++)
            fwrite(FrameworkDistanceHistogram[i][j][k][l],FrameworkSpacingHistogramSize[i],sizeof(REAL),FilePtr);
      }
    }
  }

  // write residence-times data
  fwrite(ComputeResidenceTimes,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteResidenceTimesEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(ResidenceTimesHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(RangeResidenceTimes,NumberOfSystems,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeResidenceTimes[i])
    {
      fwrite(ResidenceTimesHistogram[i],ResidenceTimesHistogramSize[i],sizeof(REAL),FilePtr);
      fwrite(ResidenceOriginAdsorbates[i],NumberOfAdsorbateMolecules[i],sizeof(long long),FilePtr);
      fwrite(ResidenceOriginCations[i],NumberOfCationMolecules[i],sizeof(long long),FilePtr);
      fwrite(ResidenceStatusAdsorbates[i],NumberOfAdsorbateMolecules[i],sizeof(int),FilePtr);
      fwrite(ResidenceStatusCations[i],NumberOfCationMolecules[i],sizeof(int),FilePtr);

      fwrite(ResidenceTimesFractionAdsorbates[i],NumberOfAdsorbateMolecules[i],sizeof(REAL[NR_BLOCKS]),FilePtr);
      fwrite(ResidenceTimesFractionCations[i],NumberOfCationMolecules[i],sizeof(REAL[NR_BLOCKS]),FilePtr);
      fwrite(&ResidenceTimesFractionCounts[i],1,sizeof(REAL),FilePtr);
    }
  }

  // write distance histograms
  fwrite(ComputeDistanceHistograms,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteDistanceHistogramsEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfElementsDistanceHistogram,1,sizeof(int),FilePtr);
  fwrite(&MaxRangeDistanceHistogram,1,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
    if(ComputeDistanceHistograms[i])
      fwrite(&NumberOfDistanceHistogramDefinitions[i],1,sizeof(int),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeDistanceHistograms[i])
    {
      fwrite(DistanceHistogramDefinitions[i],NumberOfDistanceHistogramDefinitions[i],sizeof(int[2][3]),FilePtr);
      for(j=0;j<NumberOfDistanceHistogramDefinitions[i];j++)
        fwrite(DistanceHistograms[i][j],NumberOfElementsDistanceHistogram,sizeof(REAL),FilePtr);
    }
  }

  // write bend angle histograms
  fwrite(ComputeBendAngleHistograms,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteBendAngleHistogramsEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfElementsBendAngleHistogram,1,sizeof(int),FilePtr);
  fwrite(&MaxRangeBendAngleHistogram,1,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
    if(ComputeBendAngleHistograms[i])
      fwrite(&NumberOfBendAngleHistogramDefinitions[i],1,sizeof(int),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeBendAngleHistograms[i])
    {
      fwrite(BendAngleHistogramDefinitions[i],NumberOfBendAngleHistogramDefinitions[i],sizeof(int[3][3]),FilePtr);
      for(j=0;j<NumberOfBendAngleHistogramDefinitions[i];j++)
        fwrite(BendAngleHistograms[i][j],NumberOfElementsBendAngleHistogram,sizeof(REAL),FilePtr);
    }
  }

  // write dihedral angle histograms
  fwrite(ComputeDihedralAngleHistograms,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteDihedralAngleHistogramsEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfElementsDihedralAngleHistogram,1,sizeof(int),FilePtr);
  fwrite(&MaxRangeDihedralAngleHistogram,1,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
    if(ComputeDihedralAngleHistograms[i])
      fwrite(&NumberOfDihedralAngleHistogramDefinitions[i],1,sizeof(int),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeDihedralAngleHistograms[i])
    {
      fwrite(DihedralAngleHistogramDefinitions[i],NumberOfDihedralAngleHistogramDefinitions[i],sizeof(int[4][3]),FilePtr);
      for(j=0;j<NumberOfDihedralAngleHistogramDefinitions[i];j++)
        fwrite(DihedralAngleHistograms[i][j],NumberOfElementsDihedralAngleHistogram,sizeof(REAL),FilePtr);
    }
  }

  // write angle-between planes histograms
  fwrite(ComputeAngleBetweenPlanesHistograms,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteAngleBetweenPlanesHistogramsEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfElementsAngleBetweenPlanesHistogram,1,sizeof(int),FilePtr);
  fwrite(&MaxRangeAngleBetweenPlanesHistogram,1,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
    if(ComputeAngleBetweenPlanesHistograms[i])
      fwrite(&NumberOfAngleBetweenPlanesHistogramDefinitions[i],1,sizeof(int),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeAngleBetweenPlanesHistograms[i])
    {
      fwrite(AngleBetweenPlanesHistogramDefinitions[i],NumberOfAngleBetweenPlanesHistogramDefinitions[i],sizeof(int[6][3]),FilePtr);
      for(j=0;j<NumberOfAngleBetweenPlanesHistogramDefinitions[i];j++)
        fwrite(AngleBetweenPlanesHistograms[i][j],NumberOfElementsAngleBetweenPlanesHistogram,sizeof(REAL),FilePtr);
    }
  }

  // sampling molecular properties (bond distance, bend angle, dihedral angle)
  fwrite(ComputeMoleculeProperties,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteMoleculePropertiesEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(BondLengthHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(BendAngleHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(DihedralHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(BondLengthRange,NumberOfSystems,sizeof(REAL),FilePtr);
  fwrite(BendAngleRange,NumberOfSystems,sizeof(REAL),FilePtr);
  fwrite(DihedralRange,NumberOfSystems,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeMoleculeProperties[i])
    {
      for(j=0;j<NumberOfComponents;j++)
      {
        for(k=0;k<MaxNumberOfBonds;k++)
          fwrite(BondLengthHistogram[i][j][k],BondLengthHistogramSize[i],sizeof(REAL),FilePtr);

        for(k=0;k<MaxNumberOfUreyBradleys;k++)
          fwrite(UreyBradleyLengthHistogram[i][j][k],BondLengthHistogramSize[i],sizeof(REAL),FilePtr);

        for(k=0;k<MaxNumberOfBends;k++)
          fwrite(BendAngleHistogram[i][j][k],BendAngleHistogramSize[i],sizeof(REAL),FilePtr);

        for(k=0;k<MaxNumberOfTorsions;k++)
          fwrite(TorsionAngleHistogram[i][j][k],DihedralHistogramSize[i],sizeof(REAL),FilePtr);

        for(k=0;k<MaxNumberOfTorsions;k++)
          fwrite(TorsionConformationHistogram[i][j][k],6,sizeof(REAL),FilePtr);
      }

      for(j=0;j<Framework[i].NumberOfBondsDefinitions;j++)
        fwrite(FrameworkBondLengthHistogram[i][j],BondLengthHistogramSize[i],sizeof(REAL),FilePtr);

      for(j=0;j<Framework[i].NumberOfUreyBradleyDefinitions;j++)
        fwrite(FrameworkUreyBradleyLengthHistogram[i][j],BondLengthHistogramSize[i],sizeof(REAL),FilePtr);

      for(j=0;j<Framework[i].NumberOfBendDefinitions;j++)
        fwrite(FrameworkBendAngleHistogram[i][j],BendAngleHistogramSize[i],sizeof(REAL),FilePtr);

      for(j=0;j<Framework[i].NumberOfTorsionDefinitions;j++)
        fwrite(FrameworkTorsionAngleHistogram[i][j],DihedralHistogramSize[i],sizeof(REAL),FilePtr);

      fwrite(FrameworkAverageBondLength[i],Framework[i].NumberOfBondsDefinitions,sizeof(REAL),FilePtr);
      fwrite(FrameworkBondLengthCount[i],Framework[i].NumberOfBondsDefinitions,sizeof(REAL),FilePtr);
      fwrite(FrameworkAverageBendAngle[i],Framework[i].NumberOfBendDefinitions,sizeof(REAL),FilePtr);
      fwrite(FrameworkBendAngleCount[i],Framework[i].NumberOfBendDefinitions,sizeof(REAL),FilePtr);
      fwrite(FrameworkAverageTorsionAngle[i],Framework[i].NumberOfTorsionDefinitions,sizeof(REAL),FilePtr);
      fwrite(FrameworkTorsionAngleCount[i],Framework[i].NumberOfTorsionDefinitions,sizeof(REAL),FilePtr);
    }
  }


  // write spectra data
  fwrite(ComputeInfraRedSpectra,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteInfraRedSpectraEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(sumw,5,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeInfraRedSpectra[i])
    {
      fwrite(SpectrumCount[i],5,sizeof(REAL),FilePtr);
      for(j=0;j<4;j++)
      {
        fwrite(Spectrum[i][j][0],4*2048,sizeof(REAL),FilePtr);
        fwrite(Spectrum[i][j][1],4*4096,sizeof(REAL),FilePtr);
        fwrite(Spectrum[i][j][2],4*8192,sizeof(REAL),FilePtr);
        fwrite(Spectrum[i][j][3],4*16384,sizeof(REAL),FilePtr);
        fwrite(Spectrum[i][j][4],4*32768,sizeof(REAL),FilePtr);

        fwrite(SpectrumAverage[i][j][0],4*2048,sizeof(REAL),FilePtr);
        fwrite(SpectrumAverage[i][j][1],4*4096,sizeof(REAL),FilePtr);
        fwrite(SpectrumAverage[i][j][2],4*8192,sizeof(REAL),FilePtr);
        fwrite(SpectrumAverage[i][j][3],4*16384,sizeof(REAL),FilePtr);
        fwrite(SpectrumAverage[i][j][4],4*32768,sizeof(REAL),FilePtr);

        fwrite(UnweightedSpectrum[i][j][0],4*2048,sizeof(REAL),FilePtr);
        fwrite(UnweightedSpectrum[i][j][1],4*4096,sizeof(REAL),FilePtr);
        fwrite(UnweightedSpectrum[i][j][2],4*8192,sizeof(REAL),FilePtr);
        fwrite(UnweightedSpectrum[i][j][3],4*16384,sizeof(REAL),FilePtr);
        fwrite(UnweightedSpectrum[i][j][4],4*32768,sizeof(REAL),FilePtr);

        fwrite(UnweightedSpectrumAverage[i][j][0],4*2048,sizeof(REAL),FilePtr);
        fwrite(UnweightedSpectrumAverage[i][j][1],4*4096,sizeof(REAL),FilePtr);
        fwrite(UnweightedSpectrumAverage[i][j][2],4*8192,sizeof(REAL),FilePtr);
        fwrite(UnweightedSpectrumAverage[i][j][3],4*16384,sizeof(REAL),FilePtr);
        fwrite(UnweightedSpectrumAverage[i][j][4],4*32768,sizeof(REAL),FilePtr);
      }
      for(j=0;j<2;j++)
      {
        for(k=0;k<NumberOfPseudoAtoms;k++)
        {
          if(NumberOfPseudoAtomsType[i][k]>0)
          {
            fwrite(SpectrumPseudoAtoms[i][j][k][0],4*2048,sizeof(REAL),FilePtr);
            fwrite(SpectrumPseudoAtoms[i][j][k][1],4*4096,sizeof(REAL),FilePtr);
            fwrite(SpectrumPseudoAtoms[i][j][k][2],4*8192,sizeof(REAL),FilePtr);
            fwrite(SpectrumPseudoAtoms[i][j][k][3],4*16384,sizeof(REAL),FilePtr);
            fwrite(SpectrumPseudoAtoms[i][j][k][4],4*32768,sizeof(REAL),FilePtr);
          }
        }
      }
    }
  }

  // write the mean-squared displacement using a modified order-N algorithm
  fwrite(ComputeMSDOrderN,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfBlockElementsMSDOrderN,1,sizeof(int),FilePtr);
  fwrite(&MaxNumberOfBlocksMSDOrderN,1,sizeof(int),FilePtr);
  fwrite(&ComputeIndividualMSDOrderN,1,sizeof(int),FilePtr);
  fwrite(&ComputeSiteTypeMSDOrderN,1,sizeof(int),FilePtr);
  fwrite(NumberOfSitesMSDOrderN,NumberOfSystems,sizeof(int),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeMSDOrderN[i])
    {
      fwrite(&WriteMSDOrderNEvery[i],1,sizeof(int),FilePtr);
      fwrite(&SampleMSDOrderNEvery[i],1,sizeof(int),FilePtr);
      fwrite(&CountMSDOrderN[i],1,sizeof(int),FilePtr);
      fwrite(&NumberOfBlocksMSDOrderN[i],1,sizeof(int),FilePtr);
      fwrite(BlockLengthMSDOrderN[i],MaxNumberOfBlocksMSDOrderN,sizeof(int),FilePtr);

      for(j=0;j<MaxNumberOfBlocksMSDOrderN;j++)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
        {
          fwrite(BlockDataMSDOrderN[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);

          if(ComputeSiteTypeMSDOrderN)
            fwrite(BlockDataSiteTypeMSDOrderN[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(int),FilePtr);
        }

        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(MsdOrderNCount[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          fwrite(MsdOrderNDirAvg[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          fwrite(MsdOrderN[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
        }

        // individual msd's
        if(ComputeIndividualMSDOrderN)
        {
          for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          {
            fwrite(MsdOrderNCountPerMolecule[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
            fwrite(MsdOrderNPerMolecule[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
            fwrite(MsdOrderNPerMoleculeDirAvg[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          }
        }

        if(ComputeSiteTypeMSDOrderN)
        {
          for(k=0;k<NumberOfSitesMSDOrderN[i];k++)
          {
            fwrite(MsdOrderNCountPerSiteType[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
            fwrite(MsdOrderNPerSiteType[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
            fwrite(MsdOrderNPerSiteTypeDirAvg[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          }
        }

        // Onsager data
        fwrite(MsdOrderNTotalOnsager[i][j],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
        fwrite(MsdOrderNTotalOnsagerDirAvg[i][j],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
        fwrite(MsdOrderNTotalOnsagerCount[i][j],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
        for(k=0;k<=NumberOfComponents;k++)
          fwrite(BlockDataMSDOrderNOnsager[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(MsdOrderNOnsagerCount[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          for(l=0;l<NumberOfComponents;l++)
          {
            fwrite(MsdOrderNOnsager[i][j][k][l],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
            fwrite(MsdOrderNOnsagerDirAvg[i][j][k][l],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          }
        }
      }
    }
  }

  // write the velocity autocorrelation function using a modified order-N algorithm
  fwrite(ComputeVACFOrderN,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfBlockElementsVACFOrderN,1,sizeof(int),FilePtr);
  fwrite(&MaxNumberOfBlocksVACFOrderN,1,sizeof(int),FilePtr);
  fwrite(&ComputeIndividualVACFOrderN,1,sizeof(int),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeVACFOrderN[i])
    {
      fwrite(&WriteVACFOrderNEvery[i],1,sizeof(int),FilePtr);
      fwrite(&SampleVACFOrderNEvery[i],1,sizeof(int),FilePtr);
      fwrite(&CountVACFOrderN[i],1,sizeof(int),FilePtr);
      fwrite(&NumberOfBlocksVACFOrderN[i],1,sizeof(int),FilePtr);
      fwrite(BlockLengthVACFOrderN[i],MaxNumberOfBlocksVACFOrderN,sizeof(int),FilePtr);

      for(j=0;j<MaxNumberOfBlocksVACFOrderN;j++)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          fwrite(BlockDataVACFOrderN[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);

        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(VacfOrderNCount[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
          fwrite(VacfOrderNDirAvg[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
          fwrite(VacfOrderN[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);
        }
        // individual msd's
        if(ComputeIndividualVACFOrderN)
        {
          for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          {
            fwrite(VacfOrderNCountPerMolecule[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
            fwrite(VacfOrderNPerMolecule[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);
            fwrite(VacfOrderNPerMoleculeDirAvg[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
          }
        }

        // Onsager data
        fwrite(VacfOrderNTotalOnsager[i][j],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);
        fwrite(VacfOrderNTotalOnsagerDirAvg[i][j],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
        fwrite(VacfOrderNTotalOnsagerCount[i][j],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
        for(k=0;k<=NumberOfComponents;k++)
          fwrite(BlockDataVACFOrderNOnsager[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);
        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(VacfOrderNOnsagerCount[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
          for(l=0;l<NumberOfComponents;l++)
          {
            fwrite(VacfOrderNOnsager[i][j][k][l],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);
            fwrite(VacfOrderNOnsagerDirAvg[i][j][k][l],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
          }
        }
      }
    }
  }

  // write of the rotational velocity autocorrelation function using a modified order-N algorithm
  fwrite(ComputeRVACFOrderN,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfBlockElementsRVACFOrderN,1,sizeof(int),FilePtr);
  fwrite(&MaxNumberOfBlocksRVACFOrderN,1,sizeof(int),FilePtr);
  fwrite(&ComputeIndividualRVACFOrderN,1,sizeof(int),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeRVACFOrderN[i])
    {
      fwrite(&WriteRVACFOrderNEvery[i],1,sizeof(int),FilePtr);
      fwrite(&SampleRVACFOrderNEvery[i],1,sizeof(int),FilePtr);
      fwrite(&CountRVACFOrderN[i],1,sizeof(int),FilePtr);
      fwrite(&NumberOfBlocksRVACFOrderN[i],1,sizeof(int),FilePtr);
      fwrite(BlockLengthRVACFOrderN[i],MaxNumberOfBlocksRVACFOrderN,sizeof(int),FilePtr);

      for(j=0;j<MaxNumberOfBlocksRVACFOrderN;j++)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          fwrite(BlockDataRVACFOrderN[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);

        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(RvacfOrderNCount[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
          fwrite(RvacfOrderNDirAvg[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
          fwrite(RvacfOrderN[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);
        }
        // individual msd's
        if(ComputeIndividualRVACFOrderN)
        {
          for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          {
            fwrite(RvacfOrderNCountPerMolecule[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
            fwrite(RvacfOrderNPerMolecule[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);
            fwrite(RvacfOrderNPerMoleculeDirAvg[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
          }
        }

        // Onsager data
        fwrite(RvacfOrderNTotalOnsager[i][j],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);
        fwrite(RvacfOrderNTotalOnsagerDirAvg[i][j],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
        fwrite(RvacfOrderNTotalOnsagerCount[i][j],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
        for(k=0;k<=NumberOfComponents;k++)
          fwrite(BlockDataRVACFOrderNOnsager[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);
        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(RvacfOrderNOnsagerCount[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
          for(l=0;l<NumberOfComponents;l++)
          {
            fwrite(RvacfOrderNOnsager[i][j][k][l],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);
            fwrite(RvacfOrderNOnsagerDirAvg[i][j][k][l],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
          }
        }
      }
    }
  }

  // write of the molecular orientation autocorrelation function using a modified order-N algorithm
  fwrite(ComputeMolecularOrientationOrderN,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfBlockElementsMolecularOrientationOrderN,1,sizeof(int),FilePtr);
  fwrite(&MaxNumberOfBlocksMolecularOrientationOrderN,1,sizeof(int),FilePtr);
  fwrite(&ComputeIndividualMolecularOrientationOrderN,1,sizeof(int),FilePtr);
  fwrite(&MolecularOrientationType,1,sizeof(int),FilePtr);
  fwrite(&MolecularOrientationVector,1,sizeof(VECTOR),FilePtr);
  fwrite(&MolecularOrientationGroup,1,sizeof(int),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeMolecularOrientationOrderN[i])
    {
      fwrite(&WriteMolecularOrientationOrderNEvery[i],1,sizeof(int),FilePtr);
      fwrite(&SampleMolecularOrientationOrderNEvery[i],1,sizeof(int),FilePtr);
      fwrite(&CountMolecularOrientationOrderN[i],1,sizeof(int),FilePtr);
      fwrite(&NumberOfBlocksMolecularOrientationOrderN[i],1,sizeof(int),FilePtr);
      fwrite(BlockLengthMolecularOrientationOrderN[i],MaxNumberOfBlocksMolecularOrientationOrderN,sizeof(int),FilePtr);

      for(j=0;j<MaxNumberOfBlocksMolecularOrientationOrderN;j++)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          fwrite(BlockDataMolecularOrientationOrderN[i][j][k],NumberOfBlockElementsMolecularOrientationOrderN,sizeof(VECTOR),FilePtr);

        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(MolecularOrientationOrderNCount[i][j][k],NumberOfBlockElementsMolecularOrientationOrderN,sizeof(REAL),FilePtr);
          fwrite(MolecularOrientationOrderNDirAvg[i][j][k],NumberOfBlockElementsMolecularOrientationOrderN,sizeof(REAL),FilePtr);
          fwrite(MolecularOrientationOrderN[i][j][k],NumberOfBlockElementsMolecularOrientationOrderN,sizeof(VECTOR),FilePtr);
        }
      }
    }
  }

  // write of the molecular orientation autocorrelation function using a modified order-N algorithm
  fwrite(ComputeBondOrientationOrderN,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfBlockElementsBondOrientationOrderN,1,sizeof(int),FilePtr);
  fwrite(&MaxNumberOfBlocksBondOrientationOrderN,1,sizeof(int),FilePtr);
  fwrite(&ComputeIndividualBondOrientationOrderN,1,sizeof(int),FilePtr);
  fwrite(BondOrientationAngleHistogramSize,NumberOfSystems,sizeof(int),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeBondOrientationOrderN[i])
    {
      fwrite(&WriteBondOrientationOrderNEvery[i],1,sizeof(int),FilePtr);
      fwrite(&SampleBondOrientationOrderNEvery[i],1,sizeof(int),FilePtr);

      for(f1=0;f1<Framework[i].NumberOfFrameworks;f1++)
      {

        fwrite(&NumberOfOrientationFrameworkBonds[i][f1],1,sizeof(int),FilePtr);
        fwrite(NumberOfOrientationFrameworkBondPairs[i][f1],NumberOfOrientationFrameworkBonds[i][f1],sizeof(int),FilePtr);
        fwrite(OrientationFrameworkBondTypes[i][f1],NumberOfOrientationFrameworkBonds[i][f1],sizeof(int),FilePtr);
        fwrite(OrientationFrameworkBonds[i][f1],NumberOfOrientationFrameworkBonds[i][f1],sizeof(char[2][256]),FilePtr);
      }
    }
  }

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeBondOrientationOrderN[i])
    {
      for(f1=0;f1<Framework[i].NumberOfFrameworks;f1++)
      {
        fwrite(&CountBondOrientationOrderN[i][f1],1,sizeof(int),FilePtr);
        fwrite(&NumberOfBlocksBondOrientationOrderN[i][f1],1,sizeof(int),FilePtr);
        fwrite(BlockLengthBondOrientationOrderN[i][f1],MaxNumberOfBlocksBondOrientationOrderN,sizeof(int),FilePtr);

        for(k=0;k<NumberOfOrientationFrameworkBonds[i][f1];k++)
        {
          fwrite(OrientationFrameworkBondPairs[i][f1][k],NumberOfOrientationFrameworkBondPairs[i][f1][k],sizeof(PAIR),FilePtr);
          fwrite(BondOrientationAngleDistributionFunction[i][f1][k],BondOrientationAngleHistogramSize[i],sizeof(VECTOR),FilePtr);
        }

        for(j=0;j<MaxNumberOfBlocksBondOrientationOrderN;j++)
        {
          for(k=0;k<NumberOfOrientationFrameworkBonds[i][f1];k++)
          {
            for(l=0;l<NumberOfOrientationFrameworkBondPairs[i][f1][k];l++)
              fwrite(BlockDataBondOrientationOrderN[i][f1][j][k][l],NumberOfBlockElementsBondOrientationOrderN,sizeof(VECTOR),FilePtr);

            fwrite(BondOrientationOrderNCount[i][f1][j][k],NumberOfBlockElementsBondOrientationOrderN,sizeof(REAL),FilePtr);
            fwrite(BondOrientationOrderNDirAvg[i][f1][j][k],NumberOfBlockElementsBondOrientationOrderN,sizeof(REAL),FilePtr);
            fwrite(BondOrientationOrderN[i][f1][j][k],NumberOfBlockElementsBondOrientationOrderN,sizeof(VECTOR),FilePtr);
          }
        }
      }
    }
  }


  // write the mean-square displacement  function using a conventional algorithm
  fwrite(ComputeMSD,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfBuffersMSD,1,sizeof(int),FilePtr);
  fwrite(&BufferLengthMSD,1,sizeof(int),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeMSD[i])
    {
      fwrite(&SampleMSDEvery[i],1,sizeof(int),FilePtr);
      fwrite(&WriteMSDEvery[i],1,sizeof(int),FilePtr);
      fwrite(&CountAccumulatedMSD[i],1,sizeof(int),FilePtr);
      fwrite(CountMSD[i],NumberOfBuffersMSD,sizeof(int),FilePtr);

      for(j=0;j<NumberOfBuffersMSD;j++)
      {
        fwrite(OriginMSD[i][j],NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR),FilePtr);
        fwrite(OriginOnsagerMSD[i][j],NumberOfComponents,sizeof(VECTOR),FilePtr);
        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(AcfMSD[i][j][k],BufferLengthMSD,sizeof(VECTOR),FilePtr);
          fwrite(AcfDirAvgMSD[i][j][k],BufferLengthMSD,sizeof(REAL),FilePtr);

          for(l=0;l<NumberOfComponents;l++)
          {
            fwrite(AcfOnsagerMSD[i][j][k][l],BufferLengthMSD,sizeof(VECTOR),FilePtr);
            fwrite(AcfOnsagerDirAvgMSD[i][j][k][l],BufferLengthMSD,sizeof(REAL),FilePtr);
          }
        }
      }

      for(j=0;j<NumberOfComponents;j++)
      {
        fwrite(AccumulatedAcfMSD[i][j],BufferLengthMSD,sizeof(VECTOR),FilePtr);
        fwrite(AccumulatedAcfDirAvgMSD[i][j],BufferLengthMSD,sizeof(REAL),FilePtr);

        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(AccumulatedAcfOnsagerMSD[i][j][k],BufferLengthMSD,sizeof(VECTOR),FilePtr);
          fwrite(AccumulatedAcfOnsagerDirAvgMSD[i][j][k],BufferLengthMSD,sizeof(REAL),FilePtr);
        }
      }
    }
  }

  // read the velocity autocorrelation function using a conventional algorithm
  fwrite(ComputeVACF,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&NumberOfBuffersVACF,1,sizeof(int),FilePtr);
  fwrite(&BufferLengthVACF,1,sizeof(int),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeVACF[i])
    {
      fwrite(&SampleVACFEvery[i],1,sizeof(int),FilePtr);
      fwrite(&WriteVACFEvery[i],1,sizeof(int),FilePtr);
      fwrite(&CountAccumulatedVACF[i],1,sizeof(int),FilePtr);
      fwrite(CountVACF[i],NumberOfBuffersVACF,sizeof(int),FilePtr);

      for(j=0;j<NumberOfBuffersVACF;j++)
      {
        fwrite(OriginVACF[i][j],NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR),FilePtr);
        fwrite(OriginOnsagerVACF[i][j],NumberOfComponents,sizeof(VECTOR),FilePtr);
        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(AcfVACF[i][j][k],BufferLengthVACF,sizeof(VECTOR),FilePtr);
          fwrite(AcfDirAvgVACF[i][j][k],BufferLengthVACF,sizeof(REAL),FilePtr);

          for(l=0;l<NumberOfComponents;l++)
          {
            fwrite(AcfOnsagerVACF[i][j][k][l],BufferLengthVACF,sizeof(VECTOR),FilePtr);
            fwrite(AcfOnsagerDirAvgVACF[i][j][k][l],BufferLengthVACF,sizeof(REAL),FilePtr);
          }
        }
      }

      for(j=0;j<NumberOfComponents;j++)
      {
        fwrite(AccumulatedAcfVACF[i][j],BufferLengthVACF,sizeof(VECTOR),FilePtr);
        fwrite(AccumulatedAcfDirAvgVACF[i][j],BufferLengthVACF,sizeof(REAL),FilePtr);

        for(k=0;k<NumberOfComponents;k++)
        {
          fwrite(AccumulatedAcfOnsagerVACF[i][j][k],BufferLengthVACF,sizeof(VECTOR),FilePtr);
          fwrite(AccumulatedAcfOnsagerDirAvgVACF[i][j][k],BufferLengthVACF,sizeof(REAL),FilePtr);
        }
      }
    }
  }

  // read the 3D histograms of position (i.e. 3D free energy)
  fwrite(ComputeDensityProfile3DVTKGrid,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteDensityProfile3DVTKGridEvery,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(&DensityProfile3DVTKGridPoints,1,sizeof(INT_VECTOR3),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeDensityProfile3DVTKGrid[i])
    {
      for(j=0;j<NumberOfComponents;j++)
        fwrite(DensityProfile3D[i][j],DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*
              DensityProfile3DVTKGridPoints.z,sizeof(REAL),FilePtr);
    }
  }

  fwrite(ComputeMolecularPressure,NumberOfSystems,sizeof(int),FilePtr);

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}


/*********************************************************************************************************
 * Name       | AllocateSampleMemory                                                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Allocates memory for the sampling routines.                                              *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void AllocateSampleMemory(void)
{
  int i,j;

  // sampling the radial distribution function (RDF)
  ComputeRDF=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteRDFEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  RDFHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  RDFRange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  ComputeProjectedLengths=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteProjectedLengthsEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  ProjectedLengthsHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  ProjectedLengthsRange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  ComputeProjectedAngles=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteProjectedAnglesEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  ProjectedAnglesHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  ProjectedAnglesRange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  //-----------------------------------------------------------------------------------------------------
  // CFC-RXMC : sampling lambda histogram
  //-----------------------------------------------------------------------------------------------------
  ComputeCFCRXMCLambdaHistogram=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteCFCRXMCLambdaHistogramEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  CFCRXMCLambdaHistogramBins=(int*)calloc(NumberOfSystems,sizeof(int));
  //-----------------------------------------------------------------------------------------------------

  // sampling the number-of-molecules histogram
  ComputeNumberOfMoleculesHistogram=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteNumberOfMoleculesHistogramEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfMoleculesHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfMoleculesRange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // sampling position histograms/free energies
  ComputePositionHistogram=(int*)calloc(NumberOfSystems,sizeof(int));
  PositionHistogramMappingType=(int*)calloc(NumberOfSystems,sizeof(int));
  WritePositionHistogramEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  PositionHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));

  // samples the free energy profiles in a,b,c directions
  FreeEnergyMappingType=(int*)calloc(NumberOfSystems,sizeof(int));
  ComputeFreeEnergyProfile=(int*)calloc(NumberOfSystems,sizeof(int));
  FreeEnergyHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteFreeEnergyProfileEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  // sampling the pore-size distribution (PSD)
  ComputePSDHistogram=(int*)calloc(NumberOfSystems,sizeof(int));
  WritePSDHistogramEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  PSDHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  PSDRange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // sampling the end-to-end histograms
  ComputeEndToEndDistanceHistogram=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteEndToEndDistanceHistogramEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  EndToEndHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  EndToEndRange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // sampling the energy histogram
  ComputeEnergyHistogram=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteEnergyHistogramEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  EnergyHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  EnergyHistogramLowerLimit=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  EnergyHistogramUpperLimit=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // sampling the thermodynamic factor
  ComputeThermoDynamicFactor=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteThermoDynamicFactorEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  // sampling the inter-framework spacing histogram
  ComputeFrameworkSpacingHistogram=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteFrameworkSpacingHistogramEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  FrameworkSpacingHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  FrameworkSpacingRange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  OriginalFrameworkShiftDirAvg=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  OriginalFrameworkShift=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
  for(i=0;i<NumberOfSystems;i++)
  {
    OriginalFrameworkShiftDirAvg[i]=(REAL**)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL*));
    OriginalFrameworkShift[i]=(VECTOR**)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR*));
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
    {
      OriginalFrameworkShiftDirAvg[i][j]=(REAL*)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL));
      OriginalFrameworkShift[i][j]=(VECTOR*)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR));
    }
  }

  // sampling histograms of the residence times
  ComputeResidenceTimes=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteResidenceTimesEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  ResidenceTimesHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  RangeResidenceTimes=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // sampling histograms of the distance between 2 selected atoms
  ComputeDistanceHistograms=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteDistanceHistogramsEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfDistanceHistogramDefinitions=(int*)calloc(NumberOfSystems,sizeof(int));
  DistanceHistogramDefinitions=(int (**)[2][3])calloc(NumberOfSystems,sizeof(int (*)[2][3]));
  DistanceHistogramPairs=(ATOM*(**)[2])calloc(NumberOfSystems,sizeof(ATOM*(*)[2]));

  // sampling histograms of the bend angle between 3 selected atoms
  ComputeBendAngleHistograms=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteBendAngleHistogramsEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfBendAngleHistogramDefinitions=(int*)calloc(NumberOfSystems,sizeof(int));
  BendAngleHistogramDefinitions=(int (**)[3][3])calloc(NumberOfSystems,sizeof(int (*)[3][3]));
  BendAngleHistogramPairs=(ATOM*(**)[3])calloc(NumberOfSystems,sizeof(ATOM*(*)[3]));

  // sampling histograms of the dihedral angle between 4 selected atoms
  ComputeDihedralAngleHistograms=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteDihedralAngleHistogramsEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfDihedralAngleHistogramDefinitions=(int*)calloc(NumberOfSystems,sizeof(int));
  DihedralAngleHistogramDefinitions=(int (**)[4][3])calloc(NumberOfSystems,sizeof(int (*)[4][3]));
  DihedralAngleHistogramPairs=(ATOM*(**)[4])calloc(NumberOfSystems,sizeof(ATOM*(*)[4]));

  // sampling histograms of the angle between two planes (each formed by 3 chosen atoms)
  ComputeAngleBetweenPlanesHistograms=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteAngleBetweenPlanesHistogramsEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfAngleBetweenPlanesHistogramDefinitions=(int*)calloc(NumberOfSystems,sizeof(int));
  AngleBetweenPlanesHistogramDefinitions=(int (**)[6][3])calloc(NumberOfSystems,sizeof(int (*)[6][3]));
  AngleBetweenPlanesHistogramPairs=(ATOM*(**)[6])calloc(NumberOfSystems,sizeof(ATOM*(*)[6]));

  // sampling molecular properties (bond distance, bend angle, dihedral angle)
  ComputeMoleculeProperties=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteMoleculePropertiesEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  BondLengthHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  BendAngleHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  DihedralHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  BondLengthRange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  BendAngleRange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  DihedralRange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // sampling the IR spectra (spacings: 2048, 4196, 8192, 16384, 32768 points)
  ComputeInfraRedSpectra=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteInfraRedSpectraEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  // sampling the mean-squared displacement using a modified order-N algorithm
  ComputeMSDOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
  SampleMSDOrderNEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteMSDOrderNEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfSitesMSDOrderN=(int*)calloc(NumberOfSystems,sizeof(int));

  // sampling the velocity autocorrelation function using a modified order-N algorithm
  ComputeVACFOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
  SampleVACFOrderNEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteVACFOrderNEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  // sampling of the rotational velocity autocorrelation function using a modified order-N algorithm
  ComputeRVACFOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
  SampleRVACFOrderNEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteRVACFOrderNEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  // sampling of the molecular orientation autocorrelation function using a modified order-N algorithm
  ComputeMolecularOrientationOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
  SampleMolecularOrientationOrderNEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteMolecularOrientationOrderNEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  // sampling of the bond orientation autocorrelation function using a modified order-N algorithm
  ComputeBondOrientationOrderN=(int*)calloc(NumberOfSystems,sizeof(int));
  SampleBondOrientationOrderNEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteBondOrientationOrderNEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  NumberOfOrientationFrameworkBonds=(int **)calloc(NumberOfSystems,sizeof(int*));
  OrientationFrameworkBondTypes=(int ***)calloc(NumberOfSystems,sizeof(int**));
  OrientationFrameworkBonds=(char(***)[2][256])calloc(NumberOfSystems,sizeof(char(**)[2][256]));

  NumberOfOrientationFrameworkBondPairs=(int***)calloc(NumberOfSystems,sizeof(int**));
  OrientationFrameworkBondPairs=(PAIR****)calloc(NumberOfSystems,sizeof(PAIR***));

  BondOrientationAngleHistogramSize=(int*)calloc(NumberOfSystems,sizeof(int));
  BondOrientationAngleDistributionFunction=(VECTOR****)calloc(NumberOfSystems,sizeof(VECTOR***));

  for(i=0;i<NumberOfSystems;i++)
  {
    NumberOfOrientationFrameworkBonds[i]=(int *)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
    OrientationFrameworkBondTypes[i]=(int **)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
    OrientationFrameworkBonds[i]=(char(**)[2][256])calloc(Framework[i].NumberOfFrameworks,sizeof(char(*)[2][256]));

    NumberOfOrientationFrameworkBondPairs[i]=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
    OrientationFrameworkBondPairs[i]=(PAIR***)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR**));
    BondOrientationAngleDistributionFunction[i]=(VECTOR***)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR**));
  }

  // sampling the mean-square displacement function using a conventional algorithm
  ComputeMSD=(int*)calloc(NumberOfSystems,sizeof(int));
  SampleMSDEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteMSDEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  // sampling the velocity autocorrelation function using a conventional algorithm
  ComputeVACF=(int*)calloc(NumberOfSystems,sizeof(int));
  SampleVACFEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteVACFEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  // read the 3D histograms of position (i.e. 3D free energy)
  WriteDensityProfile3DVTKGridEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  ComputeDensityProfile3DVTKGrid=(int*)calloc(NumberOfSystems,sizeof(int));

  // samples the cation sites and adsorption sites
  ComputeCationAndAdsorptionSites=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteCationAndAdsorptionSitesEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  // samples initial configurations for the tranmission coefficient (dcTST)
  WritedcTSTSnapShotsToFile=(int*)calloc(NumberOfSystems,sizeof(int));
  WritedcTSTSnapShotsEvery=(int*)calloc(NumberOfSystems,sizeof(int));

  // writes the molecular-pressure data
  ComputeMolecularPressure=(int*)calloc(NumberOfSystems,sizeof(int));
}


/*********************************************************************************************************
 * Name       | ReadRestartSample                                                                        *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Read the stored data of the sampling routines in binary form from a restart-file.        *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void ReadRestartSample(FILE *FilePtr)
{
  int i,j,k,l,f1;
  int intinput1,intinput2,intinput3,intinput4;
  REAL Check;
  int readversionNumber=0;

  // initialize and allocate memory
  AllocateSampleMemory();


  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumber)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  // read data for the histograms of the radial distribution function
  fread(ComputeRDF,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteRDFEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(RDFHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(RDFRange,NumberOfSystems,sizeof(REAL),FilePtr);

  SampleRadialDistributionFunction(ALLOCATE);

  fread(CountRDF,NumberOfSystems,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeRDF[i])
    {
      for(j=0;j<NumberOfPseudoAtoms;j++)
        fread(RadialDistributionFunctionWithFramework[i][j],RDFHistogramSize[i],sizeof(REAL),FilePtr);
      for(j=0;j<NumberOfPseudoAtoms;j++)
        for(k=0;k<NumberOfPseudoAtoms;k++)
          fread(RadialDistributionFunction[i][j][k],RDFHistogramSize[i],sizeof(REAL),FilePtr);
    }
  }

  // read data for the histograms of the projected lengths
  fread(ComputeProjectedLengths,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteProjectedLengthsEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(ProjectedLengthsHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(ProjectedLengthsRange,NumberOfSystems,sizeof(REAL),FilePtr);

  SampleProjectedLengthsDistributionFunction(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeProjectedLengths[i])
    {
      fread(CountProjectedLengths[i],NumberOfComponents,sizeof(REAL),FilePtr);
      fread(ProjectedLengthsAverage[i],NumberOfComponents,sizeof(VECTOR),FilePtr);
      for(j=0;j<NumberOfComponents;j++)
        fread(ProjectedLengthsDistributionFunction[i][j],ProjectedLengthsHistogramSize[i],sizeof(VECTOR),FilePtr);
    }
  }

  // read data for the histograms of the projected angles
  fread(ComputeProjectedAngles,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteProjectedAnglesEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(ProjectedAnglesHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(ProjectedAnglesRange,NumberOfSystems,sizeof(REAL),FilePtr);

  SampleProjectedAnglesDistributionFunction(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeProjectedAngles[i])
    {
      fread(CountProjectedAngles[i],NumberOfComponents,sizeof(REAL),FilePtr);
      for(j=0;j<NumberOfComponents;j++)
        fread(ProjectedAnglesDistributionFunction[i][j],ProjectedAnglesHistogramSize[i],sizeof(VECTOR),FilePtr);
    }
  }

  // read data for the histograms of the number of molecules
  fread(ComputeNumberOfMoleculesHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteNumberOfMoleculesHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(NumberOfMoleculesHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(NumberOfMoleculesRange,NumberOfSystems,sizeof(REAL),FilePtr);

  SampleNumberOfMoleculesHistogram(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeNumberOfMoleculesHistogram[i])
    {
      for(j=0;j<NumberOfComponents+1;j++)
        fread(NumberOfMoleculesHistogram[i][j],NumberOfMoleculesHistogramSize[i],sizeof(REAL),FilePtr);
    }
  }

  // read data for the histograms of the position (free energy)
  fread(ComputePositionHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fread(PositionHistogramMappingType,NumberOfSystems,sizeof(int),FilePtr);
  fread(WritePositionHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(PositionHistogramSize,NumberOfSystems,sizeof(int),FilePtr);

  SamplePositionHistogram(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputePositionHistogram[i])
    {
      for(j=0;j<NumberOfComponents;j++)
        for(k=0;k<Components[j].NumberOfAtoms+2;k++)
        {
          fread(PositionABCHistogram[i][j][k],PositionHistogramSize[i],sizeof(VECTOR),FilePtr);
          fread(Position2DDiagonalHistogram[i][j][k],PositionHistogramSize[i],sizeof(VECTOR),FilePtr);
          fread(Position3DDiagonalHistogram[i][j][k],PositionHistogramSize[i],sizeof(VECTOR4),FilePtr);
        }
    }
  }

  // read free energy profile data
  fread(FreeEnergyMappingType,NumberOfSystems,sizeof(int),FilePtr);
  fread(ComputeFreeEnergyProfile,NumberOfSystems,sizeof(int),FilePtr);
  fread(FreeEnergyHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteFreeEnergyProfileEvery,NumberOfSystems,sizeof(int),FilePtr);

  SampleFreeEnergyProfile(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeFreeEnergyProfile[i])
    {
      for(j=0;j<NumberOfComponents;j++)
      {
        fread(RosenBinSum[i][j],FreeEnergyHistogramSize[i]+1,sizeof(REAL[13]),FilePtr);
        fread(RosenBinCount[i][j],FreeEnergyHistogramSize[i]+1,sizeof(REAL[13]),FilePtr);
        fread(RosenBinSumSquared[i][j],FreeEnergyHistogramSize[i]+1,sizeof(REAL[13]),FilePtr);
      }
    }
  }

  // sampling the pore-size distribution (PSD)
  fread(ComputePSDHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fread(WritePSDHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(PSDHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(PSDRange,NumberOfSystems,sizeof(REAL),FilePtr);

  SamplePoreSizeDistribution(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputePSDHistogram[i])
      fread(PoreSizeDistributionHistogram[i],PSDHistogramSize[i],sizeof(REAL),FilePtr);
  }

  // read the end-to-end distribution
  fread(ComputeEndToEndDistanceHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteEndToEndDistanceHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(EndToEndHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(EndToEndRange,NumberOfSystems,sizeof(REAL),FilePtr);

  SampleEndToEndDistanceHistogram(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeEndToEndDistanceHistogram[i])
    {
      for(j=0;j<NumberOfComponents;j++)
        fread(EndToEndDistanceHistogram[i][j],EndToEndHistogramSize[i],sizeof(REAL),FilePtr);
    }
  }

  // read data for the histograms of the energy
  fread(ComputeEnergyHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteEnergyHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(EnergyHistogramSize,NumberOfSystems,sizeof(int),FilePtr);

  SampleEnergyHistogram(ALLOCATE);

  fread(EnergyHistogramLowerLimit,NumberOfSystems,sizeof(REAL),FilePtr);
  fread(EnergyHistogramUpperLimit,NumberOfSystems,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeEnergyHistogram[i])
    {
      for(j=0;j<4;j++)
        fread(EnergyHistogram[i][j],EnergyHistogramSize[i],sizeof(REAL),FilePtr);
    }
  }

  // read data for the thermodynamic factor
  fread(ComputeThermoDynamicFactor,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteThermoDynamicFactorEvery,NumberOfSystems,sizeof(int),FilePtr);

  SampleThermoDynamicsFactor(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeThermoDynamicFactor[i])
    {
      fread(ThermoDynamicFactorNumberOfSamples,NumberOfComponents,sizeof(REAL),FilePtr);
      fread(ThermoDynamicFactorNumberOfMolecules[i],NumberOfComponents,sizeof(REAL),FilePtr);
      for(k=0;k<NumberOfComponents;k++)
        fread(ThermoDynamicFactorNumberOfMoleculesCrossTerm[i][k],NumberOfComponents,sizeof(REAL),FilePtr);
    }
  }

  // read the inter-framework spacing histogram
  fread(ComputeFrameworkSpacingHistogram,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteFrameworkSpacingHistogramEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(FrameworkSpacingHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(FrameworkSpacingRange,NumberOfSystems,sizeof(REAL),FilePtr);

  SampleFrameworkSpacingHistogram(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeFrameworkSpacingHistogram[i])
    {
      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      {
        fread(OriginalFrameworkShiftDirAvg[i][j],Framework[i].NumberOfFrameworks,sizeof(REAL),FilePtr);
        fread(OriginalFrameworkShift[i][j],Framework[i].NumberOfFrameworks,sizeof(VECTOR),FilePtr);
        for(k=0;k<Framework[i].NumberOfFrameworks;k++)
          for(l=0;l<4;l++)
            fread(FrameworkDistanceHistogram[i][j][k][l],FrameworkSpacingHistogramSize[i],sizeof(REAL),FilePtr);
      }
    }
  }

  // read residence-times data
  fread(ComputeResidenceTimes,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteResidenceTimesEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(ResidenceTimesHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(RangeResidenceTimes,NumberOfSystems,sizeof(REAL),FilePtr);
  SampleResidenceTimes(ALLOCATE);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeResidenceTimes[i])
    {
      fread(ResidenceTimesHistogram[i],ResidenceTimesHistogramSize[i],sizeof(REAL),FilePtr);
      fread(ResidenceOriginAdsorbates[i],NumberOfAdsorbateMolecules[i],sizeof(long long),FilePtr);
      fread(ResidenceOriginCations[i],NumberOfCationMolecules[i],sizeof(long long),FilePtr);
      fread(ResidenceStatusAdsorbates[i],NumberOfAdsorbateMolecules[i],sizeof(int),FilePtr);
      fread(ResidenceStatusCations[i],NumberOfCationMolecules[i],sizeof(int),FilePtr);

      fread(ResidenceTimesFractionAdsorbates[i],NumberOfAdsorbateMolecules[i],sizeof(REAL[NR_BLOCKS]),FilePtr);
      fread(ResidenceTimesFractionCations[i],NumberOfCationMolecules[i],sizeof(REAL[NR_BLOCKS]),FilePtr);
      fread(&ResidenceTimesFractionCounts[i],1,sizeof(REAL),FilePtr);
    }
  }

  // read distance histograms
  fread(ComputeDistanceHistograms,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteDistanceHistogramsEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfElementsDistanceHistogram,1,sizeof(int),FilePtr);
  fread(&MaxRangeDistanceHistogram,1,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
    if(ComputeDistanceHistograms[i])
      fread(&NumberOfDistanceHistogramDefinitions[i],1,sizeof(int),FilePtr);
  SampleDistanceHistogram(ALLOCATE);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeDistanceHistograms[i])
    {
      DistanceHistogramDefinitions[i]=(int(*)[2][3])calloc(NumberOfDistanceHistogramDefinitions[i],sizeof(int[2][3]));

      fread(DistanceHistogramDefinitions[i],NumberOfDistanceHistogramDefinitions[i],sizeof(int[2][3]),FilePtr);
      for(j=0;j<NumberOfDistanceHistogramDefinitions[i];j++)
        fread(DistanceHistograms[i][j],NumberOfElementsDistanceHistogram,sizeof(REAL),FilePtr);

      // re-create the pointers to the atoms from the definitions
      if(NumberOfDistanceHistogramDefinitions[i]>0)
      {
        DistanceHistogramPairs[i]=(ATOM*(*)[2])calloc(NumberOfDistanceHistogramDefinitions[i],sizeof(ATOM*[2]));
        for(j=0;j<NumberOfDistanceHistogramDefinitions[i];j++)
        {
          intinput1=DistanceHistogramDefinitions[i][j][0][1];
          intinput2=DistanceHistogramDefinitions[i][j][0][2];
          switch(DistanceHistogramDefinitions[i][j][0][0])
          {
            case FRAMEWORK:
              DistanceHistogramPairs[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              DistanceHistogramPairs[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              DistanceHistogramPairs[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }
          intinput3=DistanceHistogramDefinitions[i][j][1][1];
          intinput4=DistanceHistogramDefinitions[i][j][1][2];
          switch(DistanceHistogramDefinitions[i][j][1][0])
          {
            case FRAMEWORK:
              DistanceHistogramPairs[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
              break;
            case ADSORBATE:
              DistanceHistogramPairs[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
              break;
            case CATION:
              DistanceHistogramPairs[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
              break;
          }
        }
      }
    }
  }

  // read bend angle histograms
  fread(ComputeBendAngleHistograms,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteBendAngleHistogramsEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfElementsBendAngleHistogram,1,sizeof(int),FilePtr);
  fread(&MaxRangeBendAngleHistogram,1,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
    if(ComputeBendAngleHistograms[i])
      fread(&NumberOfBendAngleHistogramDefinitions[i],1,sizeof(int),FilePtr);
  SampleBendAngleHistogram(ALLOCATE);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeBendAngleHistograms[i])
    {
      BendAngleHistogramDefinitions[i]=(int(*)[3][3])calloc(NumberOfBendAngleHistogramDefinitions[i],sizeof(int[3][3]));

      fread(BendAngleHistogramDefinitions[i],NumberOfBendAngleHistogramDefinitions[i],sizeof(int[3][3]),FilePtr);
      for(j=0;j<NumberOfBendAngleHistogramDefinitions[i];j++)
        fread(BendAngleHistograms[i][j],NumberOfElementsBendAngleHistogram,sizeof(REAL),FilePtr);

      // re-create the pointers to the atoms from the definitions
      if(NumberOfBendAngleHistogramDefinitions[i]>0)
      {
        BendAngleHistogramPairs[i]=(ATOM*(*)[3])calloc(NumberOfBendAngleHistogramDefinitions[i],sizeof(ATOM*[3]));
        for(j=0;j<NumberOfBendAngleHistogramDefinitions[i];j++)
        {
          intinput1=BendAngleHistogramDefinitions[i][j][0][1];
          intinput2=BendAngleHistogramDefinitions[i][j][0][2];
          switch(BendAngleHistogramDefinitions[i][j][0][0])
          {
            case FRAMEWORK:
              BendAngleHistogramPairs[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              BendAngleHistogramPairs[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              BendAngleHistogramPairs[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }
          intinput1=BendAngleHistogramDefinitions[i][j][1][1];
          intinput2=BendAngleHistogramDefinitions[i][j][1][2];
          switch(BendAngleHistogramDefinitions[i][j][1][0])
          {
            case FRAMEWORK:
              BendAngleHistogramPairs[i][j][1]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              BendAngleHistogramPairs[i][j][1]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              BendAngleHistogramPairs[i][j][1]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }
          intinput1=BendAngleHistogramDefinitions[i][j][2][1];
          intinput2=BendAngleHistogramDefinitions[i][j][2][2];
          switch(BendAngleHistogramDefinitions[i][j][2][0])
          {
            case FRAMEWORK:
              BendAngleHistogramPairs[i][j][2]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              BendAngleHistogramPairs[i][j][2]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              BendAngleHistogramPairs[i][j][2]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }
        }
      }
    }
  }

  // read dihedral angle histograms
  fread(ComputeDihedralAngleHistograms,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteDihedralAngleHistogramsEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfElementsDihedralAngleHistogram,1,sizeof(int),FilePtr);
  fread(&MaxRangeDihedralAngleHistogram,1,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
    if(ComputeDihedralAngleHistograms[i])
      fread(&NumberOfDihedralAngleHistogramDefinitions[i],1,sizeof(int),FilePtr);
  SampleDihedralAngleHistogram(ALLOCATE);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeDihedralAngleHistograms[i])
    {
      DihedralAngleHistogramDefinitions[i]=(int(*)[4][3])calloc(NumberOfDihedralAngleHistogramDefinitions[i],sizeof(int[4][3]));

      fread(DihedralAngleHistogramDefinitions[i],NumberOfDihedralAngleHistogramDefinitions[i],sizeof(int[4][3]),FilePtr);
      for(j=0;j<NumberOfDihedralAngleHistogramDefinitions[i];j++)
        fread(DihedralAngleHistograms[i][j],NumberOfElementsDihedralAngleHistogram,sizeof(REAL),FilePtr);

      // re-create the pointers to the atoms from the definitions
      if(NumberOfDihedralAngleHistogramDefinitions[i]>0)
      {
        DihedralAngleHistogramPairs[i]=(ATOM*(*)[4])calloc(NumberOfDihedralAngleHistogramDefinitions[i],sizeof(ATOM*[4]));
        for(j=0;j<NumberOfDihedralAngleHistogramDefinitions[i];j++)
        {
          intinput1=DihedralAngleHistogramDefinitions[i][j][0][1];
          intinput2=DihedralAngleHistogramDefinitions[i][j][0][2];
          switch(DihedralAngleHistogramDefinitions[i][j][0][0])
          {
            case FRAMEWORK:
              DihedralAngleHistogramPairs[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              DihedralAngleHistogramPairs[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              DihedralAngleHistogramPairs[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }

          intinput1=DihedralAngleHistogramDefinitions[i][j][1][1];
          intinput2=DihedralAngleHistogramDefinitions[i][j][1][2];
          switch(DihedralAngleHistogramDefinitions[i][j][1][0])
          {
            case FRAMEWORK:
              DihedralAngleHistogramPairs[i][j][1]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              DihedralAngleHistogramPairs[i][j][1]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              DihedralAngleHistogramPairs[i][j][1]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }

          intinput1=DihedralAngleHistogramDefinitions[i][j][2][1];
          intinput2=DihedralAngleHistogramDefinitions[i][j][2][2];
          switch(DihedralAngleHistogramDefinitions[i][j][2][0])
          {
            case FRAMEWORK:
              DihedralAngleHistogramPairs[i][j][2]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              DihedralAngleHistogramPairs[i][j][2]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              DihedralAngleHistogramPairs[i][j][2]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }

          intinput1=DihedralAngleHistogramDefinitions[i][j][3][1];
          intinput2=DihedralAngleHistogramDefinitions[i][j][3][2];
          switch(DihedralAngleHistogramDefinitions[i][j][3][0])
          {
            case FRAMEWORK:
              DihedralAngleHistogramPairs[i][j][3]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              DihedralAngleHistogramPairs[i][j][3]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              DihedralAngleHistogramPairs[i][j][3]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }
        }
      }
    }
  }

  // read angle between planes histograms
  fread(ComputeAngleBetweenPlanesHistograms,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteAngleBetweenPlanesHistogramsEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfElementsAngleBetweenPlanesHistogram,1,sizeof(int),FilePtr);
  fread(&MaxRangeAngleBetweenPlanesHistogram,1,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
    if(ComputeAngleBetweenPlanesHistograms[i])
      fread(&NumberOfAngleBetweenPlanesHistogramDefinitions[i],1,sizeof(int),FilePtr);
  SampleAngleBetweenPlanesHistogram(ALLOCATE);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeAngleBetweenPlanesHistograms[i])
    {
      AngleBetweenPlanesHistogramDefinitions[i]=(int(*)[6][3])calloc(NumberOfAngleBetweenPlanesHistogramDefinitions[i],sizeof(int[6][3]));

      fread(AngleBetweenPlanesHistogramDefinitions[i],NumberOfAngleBetweenPlanesHistogramDefinitions[i],sizeof(int[6][3]),FilePtr);
      for(j=0;j<NumberOfAngleBetweenPlanesHistogramDefinitions[i];j++)
        fread(AngleBetweenPlanesHistograms[i][j],NumberOfElementsAngleBetweenPlanesHistogram,sizeof(REAL),FilePtr);

      // re-create the pointers to the atoms from the definitions
      if(NumberOfAngleBetweenPlanesHistogramDefinitions[i]>0)
      {
        AngleBetweenPlanesHistogramPairs[i]=(ATOM*(*)[6])calloc(NumberOfAngleBetweenPlanesHistogramDefinitions[i],sizeof(ATOM*[6]));
        for(j=0;j<NumberOfAngleBetweenPlanesHistogramDefinitions[i];j++)
        {
          intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][0][1];
          intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][0][2];
          switch(AngleBetweenPlanesHistogramDefinitions[i][j][0][0])
          {
            case FRAMEWORK:
              AngleBetweenPlanesHistogramPairs[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              AngleBetweenPlanesHistogramPairs[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              AngleBetweenPlanesHistogramPairs[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }

          intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][1][1];
          intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][1][2];
          switch(AngleBetweenPlanesHistogramDefinitions[i][j][1][0])
          {
            case FRAMEWORK:
              AngleBetweenPlanesHistogramPairs[i][j][1]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              AngleBetweenPlanesHistogramPairs[i][j][1]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              AngleBetweenPlanesHistogramPairs[i][j][1]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }

          intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][2][1];
          intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][2][2];
          switch(AngleBetweenPlanesHistogramDefinitions[i][j][2][0])
          {
            case FRAMEWORK:
              AngleBetweenPlanesHistogramPairs[i][j][2]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              AngleBetweenPlanesHistogramPairs[i][j][2]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              AngleBetweenPlanesHistogramPairs[i][j][2]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }

          intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][3][1];
          intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][3][2];
          switch(AngleBetweenPlanesHistogramDefinitions[i][j][3][0])
          {
            case FRAMEWORK:
              AngleBetweenPlanesHistogramPairs[i][j][3]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              AngleBetweenPlanesHistogramPairs[i][j][3]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              AngleBetweenPlanesHistogramPairs[i][j][3]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }

          intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][4][1];
          intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][4][2];
          switch(AngleBetweenPlanesHistogramDefinitions[i][j][4][0])
          {
            case FRAMEWORK:
              AngleBetweenPlanesHistogramPairs[i][j][4]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              AngleBetweenPlanesHistogramPairs[i][j][4]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              AngleBetweenPlanesHistogramPairs[i][j][4]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }

          intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][5][1];
          intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][5][2];
          switch(AngleBetweenPlanesHistogramDefinitions[i][j][5][0])
          {
            case FRAMEWORK:
              AngleBetweenPlanesHistogramPairs[i][j][5]=&Framework[i].Atoms[intinput1][intinput2];
              break;
            case ADSORBATE:
              AngleBetweenPlanesHistogramPairs[i][j][5]=&Adsorbates[i][intinput1].Atoms[intinput2];
              break;
            case CATION:
              AngleBetweenPlanesHistogramPairs[i][j][5]=&Cations[i][intinput1].Atoms[intinput2];
              break;
          }
        }
      }
    }
  }


  // sampling molecular properties (bond distance, bend angle, dihedral angle)
  fread(ComputeMoleculeProperties,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteMoleculePropertiesEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(BondLengthHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(BendAngleHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(DihedralHistogramSize,NumberOfSystems,sizeof(int),FilePtr);
  fread(BondLengthRange,NumberOfSystems,sizeof(REAL),FilePtr);
  fread(BendAngleRange,NumberOfSystems,sizeof(REAL),FilePtr);
  fread(DihedralRange,NumberOfSystems,sizeof(REAL),FilePtr);

  SampleMoleculePropertyHistogram(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeMoleculeProperties[i])
    {
      for(j=0;j<NumberOfComponents;j++)
      {
        for(k=0;k<MaxNumberOfBonds;k++)
          fread(BondLengthHistogram[i][j][k],BondLengthHistogramSize[i],sizeof(REAL),FilePtr);

        for(k=0;k<MaxNumberOfUreyBradleys;k++)
          fread(UreyBradleyLengthHistogram[i][j][k],BondLengthHistogramSize[i],sizeof(REAL),FilePtr);

        for(k=0;k<MaxNumberOfBends;k++)
          fread(BendAngleHistogram[i][j][k],BendAngleHistogramSize[i],sizeof(REAL),FilePtr);

        for(k=0;k<MaxNumberOfTorsions;k++)
          fread(TorsionAngleHistogram[i][j][k],DihedralHistogramSize[i],sizeof(REAL),FilePtr);

        for(k=0;k<MaxNumberOfTorsions;k++)
          fread(TorsionConformationHistogram[i][j][k],6,sizeof(REAL),FilePtr);
      }

      for(j=0;j<Framework[i].NumberOfBondsDefinitions;j++)
        fread(FrameworkBondLengthHistogram[i][j],BondLengthHistogramSize[i],sizeof(REAL),FilePtr);

      for(j=0;j<Framework[i].NumberOfUreyBradleyDefinitions;j++)
        fread(FrameworkUreyBradleyLengthHistogram[i][j],BondLengthHistogramSize[i],sizeof(REAL),FilePtr);

      for(j=0;j<Framework[i].NumberOfBendDefinitions;j++)
        fread(FrameworkBendAngleHistogram[i][j],BendAngleHistogramSize[i],sizeof(REAL),FilePtr);

      for(j=0;j<Framework[i].NumberOfTorsionDefinitions;j++)
        fread(FrameworkTorsionAngleHistogram[i][j],DihedralHistogramSize[i],sizeof(REAL),FilePtr);

      fread(FrameworkAverageBondLength[i],Framework[i].NumberOfBondsDefinitions,sizeof(REAL),FilePtr);
      fread(FrameworkBondLengthCount[i],Framework[i].NumberOfBondsDefinitions,sizeof(REAL),FilePtr);
      fread(FrameworkAverageBendAngle[i],Framework[i].NumberOfBendDefinitions,sizeof(REAL),FilePtr);
      fread(FrameworkBendAngleCount[i],Framework[i].NumberOfBendDefinitions,sizeof(REAL),FilePtr);
      fread(FrameworkAverageTorsionAngle[i],Framework[i].NumberOfTorsionDefinitions,sizeof(REAL),FilePtr);
      fread(FrameworkTorsionAngleCount[i],Framework[i].NumberOfTorsionDefinitions,sizeof(REAL),FilePtr);
    }
  }


  // read spectra data
  fread(ComputeInfraRedSpectra,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteInfraRedSpectraEvery,NumberOfSystems,sizeof(int),FilePtr);

  SampleInfraRedSpectra(ALLOCATE);
  fread(sumw,5,sizeof(REAL),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeInfraRedSpectra[i])
    {
      fread(SpectrumCount[i],5,sizeof(REAL),FilePtr);
      for(j=0;j<4;j++)
      {
        fread(Spectrum[i][j][0],4*2048,sizeof(REAL),FilePtr);
        fread(Spectrum[i][j][1],4*4096,sizeof(REAL),FilePtr);
        fread(Spectrum[i][j][2],4*8192,sizeof(REAL),FilePtr);
        fread(Spectrum[i][j][3],4*16384,sizeof(REAL),FilePtr);
        fread(Spectrum[i][j][4],4*32768,sizeof(REAL),FilePtr);

        fread(SpectrumAverage[i][j][0],4*2048,sizeof(REAL),FilePtr);
        fread(SpectrumAverage[i][j][1],4*4096,sizeof(REAL),FilePtr);
        fread(SpectrumAverage[i][j][2],4*8192,sizeof(REAL),FilePtr);
        fread(SpectrumAverage[i][j][3],4*16384,sizeof(REAL),FilePtr);
        fread(SpectrumAverage[i][j][4],4*32768,sizeof(REAL),FilePtr);

        fread(UnweightedSpectrum[i][j][0],4*2048,sizeof(REAL),FilePtr);
        fread(UnweightedSpectrum[i][j][1],4*4096,sizeof(REAL),FilePtr);
        fread(UnweightedSpectrum[i][j][2],4*8192,sizeof(REAL),FilePtr);
        fread(UnweightedSpectrum[i][j][3],4*16384,sizeof(REAL),FilePtr);
        fread(UnweightedSpectrum[i][j][4],4*32768,sizeof(REAL),FilePtr);

        fread(UnweightedSpectrumAverage[i][j][0],4*2048,sizeof(REAL),FilePtr);
        fread(UnweightedSpectrumAverage[i][j][1],4*4096,sizeof(REAL),FilePtr);
        fread(UnweightedSpectrumAverage[i][j][2],4*8192,sizeof(REAL),FilePtr);
        fread(UnweightedSpectrumAverage[i][j][3],4*16384,sizeof(REAL),FilePtr);
        fread(UnweightedSpectrumAverage[i][j][4],4*32768,sizeof(REAL),FilePtr);
      }
      for(j=0;j<2;j++)
      {
        for(k=0;k<NumberOfPseudoAtoms;k++)
        {
          if(NumberOfPseudoAtomsType[i][k]>0)
          {
            fread(SpectrumPseudoAtoms[i][j][k][0],4*2048,sizeof(REAL),FilePtr);
            fread(SpectrumPseudoAtoms[i][j][k][1],4*4096,sizeof(REAL),FilePtr);
            fread(SpectrumPseudoAtoms[i][j][k][2],4*8192,sizeof(REAL),FilePtr);
            fread(SpectrumPseudoAtoms[i][j][k][3],4*16384,sizeof(REAL),FilePtr);
            fread(SpectrumPseudoAtoms[i][j][k][4],4*32768,sizeof(REAL),FilePtr);
          }
        }
      }
    }
  }

  // read the mean-squared displacement using a modified order-N algorithm
  fread(ComputeMSDOrderN,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfBlockElementsMSDOrderN,1,sizeof(int),FilePtr);
  fread(&MaxNumberOfBlocksMSDOrderN,1,sizeof(int),FilePtr);
  fread(&ComputeIndividualMSDOrderN,1,sizeof(int),FilePtr);
  fread(&ComputeSiteTypeMSDOrderN,1,sizeof(int),FilePtr);
  fread(NumberOfSitesMSDOrderN,NumberOfSystems,sizeof(int),FilePtr);

  SampleMeanSquaredDisplacementOrderN(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeMSDOrderN[i])
    {
      fread(&WriteMSDOrderNEvery[i],1,sizeof(int),FilePtr);
      fread(&SampleMSDOrderNEvery[i],1,sizeof(int),FilePtr);
      fread(&CountMSDOrderN[i],1,sizeof(int),FilePtr);
      fread(&NumberOfBlocksMSDOrderN[i],1,sizeof(int),FilePtr);
      fread(BlockLengthMSDOrderN[i],MaxNumberOfBlocksMSDOrderN,sizeof(int),FilePtr);

      for(j=0;j<MaxNumberOfBlocksMSDOrderN;j++)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
        {
          fread(BlockDataMSDOrderN[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);

          if(ComputeSiteTypeMSDOrderN)
            fread(BlockDataSiteTypeMSDOrderN[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(int),FilePtr);
        }

        for(k=0;k<NumberOfComponents;k++)
        {
          fread(MsdOrderNCount[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          fread(MsdOrderNDirAvg[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          fread(MsdOrderN[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
        }

        // individual msd's
        if(ComputeIndividualMSDOrderN)
        {
          for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          {
            fread(MsdOrderNCountPerMolecule[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
            fread(MsdOrderNPerMolecule[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
            fread(MsdOrderNPerMoleculeDirAvg[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          }
        }

        if(ComputeSiteTypeMSDOrderN)
        {
          for(k=0;k<NumberOfSitesMSDOrderN[i];k++)
          {
            fread(MsdOrderNCountPerSiteType[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
            fread(MsdOrderNPerSiteType[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
            fread(MsdOrderNPerSiteTypeDirAvg[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          }
        }


        // Onsager data
        fread(MsdOrderNTotalOnsager[i][j],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
        fread(MsdOrderNTotalOnsagerDirAvg[i][j],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
        fread(MsdOrderNTotalOnsagerCount[i][j],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
        for(k=0;k<=NumberOfComponents;k++)
          fread(BlockDataMSDOrderNOnsager[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
        for(k=0;k<NumberOfComponents;k++)
        {
          fread(MsdOrderNOnsagerCount[i][j][k],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          for(l=0;l<NumberOfComponents;l++)
          {
            fread(MsdOrderNOnsager[i][j][k][l],NumberOfBlockElementsMSDOrderN,sizeof(VECTOR),FilePtr);
            fread(MsdOrderNOnsagerDirAvg[i][j][k][l],NumberOfBlockElementsMSDOrderN,sizeof(REAL),FilePtr);
          }
        }
      }
    }
  }

  // read the velocity autocorrelation function using a modified order-N algorithm
  fread(ComputeVACFOrderN,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfBlockElementsVACFOrderN,1,sizeof(int),FilePtr);
  fread(&MaxNumberOfBlocksVACFOrderN,1,sizeof(int),FilePtr);
  fread(&ComputeIndividualVACFOrderN,1,sizeof(int),FilePtr);

  SampleVelocityAutoCorrelationFunctionOrderN(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeVACFOrderN[i])
    {
      fread(&WriteVACFOrderNEvery[i],1,sizeof(int),FilePtr);
      fread(&SampleVACFOrderNEvery[i],1,sizeof(int),FilePtr);
      fread(&CountVACFOrderN[i],1,sizeof(int),FilePtr);
      fread(&NumberOfBlocksVACFOrderN[i],1,sizeof(int),FilePtr);
      fread(BlockLengthVACFOrderN[i],MaxNumberOfBlocksVACFOrderN,sizeof(int),FilePtr);

      for(j=0;j<MaxNumberOfBlocksVACFOrderN;j++)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          fread(BlockDataVACFOrderN[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);

        for(k=0;k<NumberOfComponents;k++)
        {
          fread(VacfOrderNCount[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
          fread(VacfOrderNDirAvg[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
          fread(VacfOrderN[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);
        }
        // individual msd's
        if(ComputeIndividualVACFOrderN)
        {
          for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          {
            fread(VacfOrderNCountPerMolecule[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
            fread(VacfOrderNPerMolecule[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);
            fread(VacfOrderNPerMoleculeDirAvg[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
          }
        }

        // Onsager data
        fread(VacfOrderNTotalOnsager[i][j],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);
        fread(VacfOrderNTotalOnsagerDirAvg[i][j],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
        fread(VacfOrderNTotalOnsagerCount[i][j],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
        for(k=0;k<=NumberOfComponents;k++)
          fread(BlockDataVACFOrderNOnsager[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);
        for(k=0;k<NumberOfComponents;k++)
        {
          fread(VacfOrderNOnsagerCount[i][j][k],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
          for(l=0;l<NumberOfComponents;l++)
          {
            fread(VacfOrderNOnsager[i][j][k][l],NumberOfBlockElementsVACFOrderN,sizeof(VECTOR),FilePtr);
            fread(VacfOrderNOnsagerDirAvg[i][j][k][l],NumberOfBlockElementsVACFOrderN,sizeof(REAL),FilePtr);
          }
        }
      }
    }
  }

  // read of the rotational velocity autocorrelation function using a modified order-N algorithm
  fread(ComputeRVACFOrderN,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfBlockElementsRVACFOrderN,1,sizeof(int),FilePtr);
  fread(&MaxNumberOfBlocksRVACFOrderN,1,sizeof(int),FilePtr);
  fread(&ComputeIndividualRVACFOrderN,1,sizeof(int),FilePtr);

  SampleRotationalVelocityAutoCorrelationFunctionOrderN(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeRVACFOrderN[i])
    {
      fread(&WriteRVACFOrderNEvery[i],1,sizeof(int),FilePtr);
      fread(&SampleRVACFOrderNEvery[i],1,sizeof(int),FilePtr);
      fread(&CountRVACFOrderN[i],1,sizeof(int),FilePtr);
      fread(&NumberOfBlocksRVACFOrderN[i],1,sizeof(int),FilePtr);
      fread(BlockLengthRVACFOrderN[i],MaxNumberOfBlocksRVACFOrderN,sizeof(int),FilePtr);

      for(j=0;j<MaxNumberOfBlocksRVACFOrderN;j++)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          fread(BlockDataRVACFOrderN[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);

        for(k=0;k<NumberOfComponents;k++)
        {
          fread(RvacfOrderNCount[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
          fread(RvacfOrderNDirAvg[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
          fread(RvacfOrderN[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);
        }
        // individual msd's
        if(ComputeIndividualRVACFOrderN)
        {
          for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          {
            fread(RvacfOrderNCountPerMolecule[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
            fread(RvacfOrderNPerMolecule[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);
            fread(RvacfOrderNPerMoleculeDirAvg[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
          }
        }

        // Onsager data
        fread(RvacfOrderNTotalOnsager[i][j],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);
        fread(RvacfOrderNTotalOnsagerDirAvg[i][j],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
        fread(RvacfOrderNTotalOnsagerCount[i][j],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
        for(k=0;k<=NumberOfComponents;k++)
          fread(BlockDataRVACFOrderNOnsager[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);
        for(k=0;k<NumberOfComponents;k++)
        {
          fread(RvacfOrderNOnsagerCount[i][j][k],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
          for(l=0;l<NumberOfComponents;l++)
          {
            fread(RvacfOrderNOnsager[i][j][k][l],NumberOfBlockElementsRVACFOrderN,sizeof(VECTOR),FilePtr);
            fread(RvacfOrderNOnsagerDirAvg[i][j][k][l],NumberOfBlockElementsRVACFOrderN,sizeof(REAL),FilePtr);
          }
        }
      }
    }
  }

  // read of the molcular orientation autocorrelation function using a modified order-N algorithm
  fread(ComputeMolecularOrientationOrderN,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfBlockElementsMolecularOrientationOrderN,1,sizeof(int),FilePtr);
  fread(&MaxNumberOfBlocksMolecularOrientationOrderN,1,sizeof(int),FilePtr);
  fread(&ComputeIndividualMolecularOrientationOrderN,1,sizeof(int),FilePtr);
  fread(&MolecularOrientationType,1,sizeof(int),FilePtr);
  fread(&MolecularOrientationVector,1,sizeof(VECTOR),FilePtr);
  fread(&MolecularOrientationGroup,1,sizeof(int),FilePtr);

  SampleMolecularOrientationAutoCorrelationFunctionOrderN(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeMolecularOrientationOrderN[i])
    {
      fread(&WriteMolecularOrientationOrderNEvery[i],1,sizeof(int),FilePtr);
      fread(&SampleMolecularOrientationOrderNEvery[i],1,sizeof(int),FilePtr);
      fread(&CountMolecularOrientationOrderN[i],1,sizeof(int),FilePtr);
      fread(&NumberOfBlocksMolecularOrientationOrderN[i],1,sizeof(int),FilePtr);
      fread(BlockLengthMolecularOrientationOrderN[i],MaxNumberOfBlocksMolecularOrientationOrderN,sizeof(int),FilePtr);

      for(j=0;j<MaxNumberOfBlocksMolecularOrientationOrderN;j++)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i];k++)
          fread(BlockDataMolecularOrientationOrderN[i][j][k],NumberOfBlockElementsMolecularOrientationOrderN,sizeof(VECTOR),FilePtr);

        for(k=0;k<NumberOfComponents;k++)
        {
          fread(MolecularOrientationOrderNCount[i][j][k],NumberOfBlockElementsMolecularOrientationOrderN,sizeof(REAL),FilePtr);
          fread(MolecularOrientationOrderNDirAvg[i][j][k],NumberOfBlockElementsMolecularOrientationOrderN,sizeof(REAL),FilePtr);
          fread(MolecularOrientationOrderN[i][j][k],NumberOfBlockElementsMolecularOrientationOrderN,sizeof(VECTOR),FilePtr);
        }
      }
    }
  }

  // read of the bond orientation autocorrelation function using a modified order-N algorithm
  fread(ComputeBondOrientationOrderN,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfBlockElementsBondOrientationOrderN,1,sizeof(int),FilePtr);
  fread(&MaxNumberOfBlocksBondOrientationOrderN,1,sizeof(int),FilePtr);
  fread(&ComputeIndividualBondOrientationOrderN,1,sizeof(int),FilePtr);
  fread(BondOrientationAngleHistogramSize,NumberOfSystems,sizeof(int),FilePtr);


  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeBondOrientationOrderN[i])
    {
      fread(&WriteBondOrientationOrderNEvery[i],1,sizeof(int),FilePtr);
      fread(&SampleBondOrientationOrderNEvery[i],1,sizeof(int),FilePtr);

      for(f1=0;f1<Framework[i].NumberOfFrameworks;f1++)
      {
        fread(&NumberOfOrientationFrameworkBonds[i][f1],1,sizeof(int),FilePtr);

        OrientationFrameworkBondTypes[i][f1]=(int *)calloc(NumberOfOrientationFrameworkBonds[i][f1],sizeof(int));
        NumberOfOrientationFrameworkBondPairs[i][f1]=(int*)calloc(NumberOfOrientationFrameworkBonds[i][f1],sizeof(int));
        OrientationFrameworkBondPairs[i][f1]=(PAIR**)calloc(NumberOfOrientationFrameworkBonds[i][f1],sizeof(PAIR*));

        fread(NumberOfOrientationFrameworkBondPairs[i][f1],NumberOfOrientationFrameworkBonds[i][f1],sizeof(int),FilePtr);
        fread(OrientationFrameworkBondTypes[i][f1],NumberOfOrientationFrameworkBonds[i][f1],sizeof(int),FilePtr);

        OrientationFrameworkBonds[i][f1]=(char(*)[2][256])calloc((NumberOfOrientationFrameworkBonds[i][f1]),sizeof(char[2][256]));
        fread(OrientationFrameworkBonds[i][f1],NumberOfOrientationFrameworkBonds[i][f1],sizeof(char[2][256]),FilePtr);
      }
    }
  }

  SampleBondOrientationAutoCorrelationFunctionOrderN(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeBondOrientationOrderN[i])
    {
      for(f1=0;f1<Framework[i].NumberOfFrameworks;f1++)
      {
        fread(&CountBondOrientationOrderN[i][f1],1,sizeof(int),FilePtr);
        fread(&NumberOfBlocksBondOrientationOrderN[i][f1],1,sizeof(int),FilePtr);
        fread(BlockLengthBondOrientationOrderN[i][f1],MaxNumberOfBlocksBondOrientationOrderN,sizeof(int),FilePtr);

        for(k=0;k<NumberOfOrientationFrameworkBonds[i][f1];k++)
        {
          OrientationFrameworkBondPairs[i][f1][k]=(PAIR*)calloc(NumberOfOrientationFrameworkBondPairs[i][f1][k],sizeof(PAIR));
          fread(OrientationFrameworkBondPairs[i][f1][k],NumberOfOrientationFrameworkBondPairs[i][f1][k],sizeof(PAIR),FilePtr);
          fread(BondOrientationAngleDistributionFunction[i][f1][k],BondOrientationAngleHistogramSize[i],sizeof(VECTOR),FilePtr);
        }

        for(j=0;j<MaxNumberOfBlocksBondOrientationOrderN;j++)
        {
          for(k=0;k<NumberOfOrientationFrameworkBonds[i][f1];k++)
          {
            for(l=0;l<NumberOfOrientationFrameworkBondPairs[i][f1][k];l++)
              fread(BlockDataBondOrientationOrderN[i][f1][j][k][l],NumberOfBlockElementsBondOrientationOrderN,sizeof(VECTOR),FilePtr);

            fread(BondOrientationOrderNCount[i][f1][j][k],NumberOfBlockElementsBondOrientationOrderN,sizeof(REAL),FilePtr);
            fread(BondOrientationOrderNDirAvg[i][f1][j][k],NumberOfBlockElementsBondOrientationOrderN,sizeof(REAL),FilePtr);
            fread(BondOrientationOrderN[i][f1][j][k],NumberOfBlockElementsBondOrientationOrderN,sizeof(VECTOR),FilePtr);
          }
        }
      }
    }
  }


  // read the mean-square displacement function using a conventional algorithm
  fread(ComputeMSD,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfBuffersMSD,1,sizeof(int),FilePtr);
  fread(&BufferLengthMSD,1,sizeof(int),FilePtr);

  SampleMeanSquaredDisplacement(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeMSD[i])
    {
      fread(&SampleMSDEvery[i],1,sizeof(int),FilePtr);
      fread(&WriteMSDEvery[i],1,sizeof(int),FilePtr);
      fread(&CountAccumulatedMSD[i],1,sizeof(int),FilePtr);
      fread(CountMSD[i],NumberOfBuffersMSD,sizeof(int),FilePtr);

      for(j=0;j<NumberOfBuffersMSD;j++)
      {
        fread(OriginMSD[i][j],NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR),FilePtr);
        fread(OriginOnsagerMSD[i][j],NumberOfComponents,sizeof(VECTOR),FilePtr);
        for(k=0;k<NumberOfComponents;k++)
        {
          fread(AcfMSD[i][j][k],BufferLengthMSD,sizeof(VECTOR),FilePtr);
          fread(AcfDirAvgMSD[i][j][k],BufferLengthMSD,sizeof(REAL),FilePtr);

          for(l=0;l<NumberOfComponents;l++)
          {
            fread(AcfOnsagerMSD[i][j][k][l],BufferLengthMSD,sizeof(VECTOR),FilePtr);
            fread(AcfOnsagerDirAvgMSD[i][j][k][l],BufferLengthMSD,sizeof(REAL),FilePtr);
          }
        }
      }

      for(j=0;j<NumberOfComponents;j++)
      {
        fread(AccumulatedAcfMSD[i][j],BufferLengthMSD,sizeof(VECTOR),FilePtr);
        fread(AccumulatedAcfDirAvgMSD[i][j],BufferLengthMSD,sizeof(REAL),FilePtr);

        for(k=0;k<NumberOfComponents;k++)
        {
          fread(AccumulatedAcfOnsagerMSD[i][j][k],BufferLengthMSD,sizeof(VECTOR),FilePtr);
          fread(AccumulatedAcfOnsagerDirAvgMSD[i][j][k],BufferLengthMSD,sizeof(REAL),FilePtr);
        }
      }
    }
  }

  // read the velocity autocorrelation function using a conventional algorithm
  fread(ComputeVACF,NumberOfSystems,sizeof(int),FilePtr);
  fread(&NumberOfBuffersVACF,1,sizeof(int),FilePtr);
  fread(&BufferLengthVACF,1,sizeof(int),FilePtr);

  SampleVelocityAutoCorrelationFunction(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeVACF[i])
    {
      fread(&SampleVACFEvery[i],1,sizeof(int),FilePtr);
      fread(&WriteVACFEvery[i],1,sizeof(int),FilePtr);
      fread(&CountAccumulatedVACF[i],1,sizeof(int),FilePtr);
      fread(CountVACF[i],NumberOfBuffersVACF,sizeof(int),FilePtr);

      for(j=0;j<NumberOfBuffersVACF;j++)
      {
        fread(OriginVACF[i][j],NumberOfAdsorbateMolecules[i]+NumberOfCationMolecules[i],sizeof(VECTOR),FilePtr);
        fread(OriginOnsagerVACF[i][j],NumberOfComponents,sizeof(VECTOR),FilePtr);
        for(k=0;k<NumberOfComponents;k++)
        {
          fread(AcfVACF[i][j][k],BufferLengthVACF,sizeof(VECTOR),FilePtr);
          fread(AcfDirAvgVACF[i][j][k],BufferLengthVACF,sizeof(REAL),FilePtr);

          for(l=0;l<NumberOfComponents;l++)
          {
            fread(AcfOnsagerVACF[i][j][k][l],BufferLengthVACF,sizeof(VECTOR),FilePtr);
            fread(AcfOnsagerDirAvgVACF[i][j][k][l],BufferLengthVACF,sizeof(REAL),FilePtr);
          }
        }
      }

      for(j=0;j<NumberOfComponents;j++)
      {
        fread(AccumulatedAcfVACF[i][j],BufferLengthVACF,sizeof(VECTOR),FilePtr);
        fread(AccumulatedAcfDirAvgVACF[i][j],BufferLengthVACF,sizeof(REAL),FilePtr);

        for(k=0;k<NumberOfComponents;k++)
        {
          fread(AccumulatedAcfOnsagerVACF[i][j][k],BufferLengthVACF,sizeof(VECTOR),FilePtr);
          fread(AccumulatedAcfOnsagerDirAvgVACF[i][j][k],BufferLengthVACF,sizeof(REAL),FilePtr);
        }
      }
    }
  }

  // read the 3D histograms of position (i.e. 3D free energy)
  fread(ComputeDensityProfile3DVTKGrid,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteDensityProfile3DVTKGridEvery,NumberOfSystems,sizeof(int),FilePtr);
  fread(&DensityProfile3DVTKGridPoints,1,sizeof(INT_VECTOR3),FilePtr);

  SampleDensityProfile3DVTKGrid(ALLOCATE);
  SampleCOMDensityProfile3DVTKGrid(ALLOCATE);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(ComputeDensityProfile3DVTKGrid[i])
    {
      for(j=0;j<NumberOfComponents;j++)
        fread(DensityProfile3D[i][j],DensityProfile3DVTKGridPoints.x*DensityProfile3DVTKGridPoints.y*
              DensityProfile3DVTKGridPoints.z,sizeof(REAL),FilePtr);
    }
  }

  // read the molecular-pressure data
  fread(ComputeMolecularPressure,NumberOfSystems,sizeof(int),FilePtr);

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartSample)\n");
    ContinueAfterCrash=FALSE;
  }
}
