/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'rigid.h' is part of RASPA-2.0

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

#ifndef SAMPLE_H
#define SAMPLE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "vector.h"

#define MAX_HISTOGRAM_REWEIGHTING 500

enum{ALLOCATE,INITIALIZE,SAMPLE,PRINT,FINALIZE};
enum{VACF=0,MSDOrderN=1};

enum{NO_MAPPING,A_MAPPING,B_MAPPING,C_MAPPING,ABC_MAPPING,
     MAP_AB_DIAGONAL,MAP_AC_DIAGONAL,MAP_BC_DIAGONAL,
     MAP_O_AB_DIAGONAL,MAP_O_AC_DIAGONAL,MAP_O_BC_DIAGONAL,
     MAP_A_BC_DIAGONAL,MAP_B_AC_DIAGONAL,MAP_C_AB_DIAGONAL,MAP_O_ABC_DIAGONAL};

// sampling the radial distribution function (RDF)
void SampleRadialDistributionFunction(int);
extern int *ComputeRDF;
extern int *WriteRDFEvery;
extern int *RDFHistogramSize;
extern REAL *RDFRange;

// sampling the projected lengths distribution function 
void SampleProjectedLengthsDistributionFunction(int);
extern int *ComputeProjectedLengths;
extern int *WriteProjectedLengthsEvery;
extern int *ProjectedLengthsHistogramSize;
extern REAL *ProjectedLengthsRange;

// sampling the projected angles distribution function 
void SampleProjectedAnglesDistributionFunction(int);
extern int *ComputeProjectedAngles;
extern int *WriteProjectedAnglesEvery;
extern int *ProjectedAnglesHistogramSize;
extern REAL *ProjectedAnglesRange;

//----------------------------------------------------------------------------------------
// CFC-RXMC : sampling lambda histogram
//----------------------------------------------------------------------------------------

void SampleCFCRXMCLambdaHistogram(int);
extern int *ComputeCFCRXMCLambdaHistogram;
extern int *WriteCFCRXMCLambdaHistogramEvery;
extern int *CFCRXMCLambdaHistogramBins;

//----------------------------------------------------------------------------------------

// sampling the number-of-molecules histogram
void SampleNumberOfMoleculesHistogram(int);
extern int *ComputeNumberOfMoleculesHistogram;
extern int *WriteNumberOfMoleculesHistogramEvery;
extern int *NumberOfMoleculesHistogramSize;
extern REAL *NumberOfMoleculesRange;

// sampling position histograms/free energies
void SamplePositionHistogram(int);
extern int *ComputePositionHistogram;
extern int *WritePositionHistogramEvery;
extern int *PositionHistogramSize;
extern int *PositionHistogramMappingType;

// samples the free energy profiles in a,b,c directions
void SampleFreeEnergyProfile(int Switch);
extern int *ComputeFreeEnergyProfile;
extern int *FreeEnergyHistogramSize;
extern int *FreeEnergyMappingType;
extern int *WriteFreeEnergyProfileEvery;

// sampling the pore-size distribution (PSD)
void SamplePoreSizeDistribution(int Choice);
extern int *ComputePSDHistogram;
extern int *PSDHistogramSize;
extern REAL *PSDRange;
extern int *WritePSDHistogramEvery;

// sampling the end-to-end histograms
void SampleEndToEndDistanceHistogram(int);
extern int *ComputeEndToEndDistanceHistogram;
extern int *WriteEndToEndDistanceHistogramEvery;
extern int *EndToEndHistogramSize;
extern REAL *EndToEndRange;

// sampling the energy histogram
void SampleEnergyHistogram(int);
extern int *ComputeEnergyHistogram;
extern int *WriteEnergyHistogramEvery;
extern int *EnergyHistogramSize;
extern REAL *EnergyHistogramLowerLimit;
extern REAL *EnergyHistogramUpperLimit;

// sampling the thermodynamic factor
void SampleThermoDynamicsFactor(int);
extern int *ComputeThermoDynamicFactor;
extern int *WriteThermoDynamicFactorEvery;

// sampling the inter-framework spacing histogram
void SampleFrameworkSpacingHistogram(int);
extern int *ComputeFrameworkSpacingHistogram;
extern int *WriteFrameworkSpacingHistogramEvery;
extern int *FrameworkSpacingHistogramSize;
extern REAL *FrameworkSpacingRange;
extern REAL ***OriginalFrameworkShiftDirAvg;
extern VECTOR ***OriginalFrameworkShift;

// sampling histograms of the residence times
void SampleResidenceTimes(int Switch);
extern int *ComputeResidenceTimes;      // whether to compute the residence times or not
extern int *WriteResidenceTimesEvery;   // writes the output every 'WriteResidenceTimesEvery' times
extern REAL **ResidenceTimesHistogram;  // the data for the histogram
extern int *ResidenceTimesHistogramSize; // the number of elements of the histogram
extern REAL *RangeResidenceTimes;     // the maximum range of the histogram

// sampling histograms of the distance between 2 selected atoms
void SampleDistanceHistogram(int Switch);
extern int *ComputeDistanceHistograms;               // whether to compute the distance histograms or not
extern int *WriteDistanceHistogramsEvery;            // writes the output every 'WriteDistanceHistogramsEvery' times
extern int NumberOfElementsDistanceHistogram;        // the number of elements of the histogram
extern REAL MaxRangeDistanceHistogram;               // the maximum range of the histogram
extern int *NumberOfDistanceHistogramDefinitions;
extern int (**DistanceHistogramDefinitions)[2][3];
extern ATOM* (**DistanceHistogramPairs)[2];

// sampling histograms of the bend angle between 3 selected atoms
void SampleBendAngleHistogram(int Switch);
extern int *ComputeBendAngleHistograms;               // whether to compute the angle histograms or not
extern int *WriteBendAngleHistogramsEvery;            // writes the output every 'WriteAngleHistogramsEvery' times
extern int NumberOfElementsBendAngleHistogram;        // the number of elements of the histogram
extern REAL MaxRangeBendAngleHistogram;               // the maximum range of the histogram
extern int *NumberOfBendAngleHistogramDefinitions;
extern int (**BendAngleHistogramDefinitions)[3][3];
extern ATOM* (**BendAngleHistogramPairs)[3];

// sampling histograms of the dihedral angle between 4 selected atoms
void SampleDihedralAngleHistogram(int Switch);
extern int *ComputeDihedralAngleHistograms;               // whether to compute the angle histograms or not
extern int *WriteDihedralAngleHistogramsEvery;            // writes the output every 'WriteDihedralAngleHistogramsEvery' times
extern int NumberOfElementsDihedralAngleHistogram;        // the number of elements of the histogram
extern REAL MaxRangeDihedralAngleHistogram;               // the maximum range of the histogram
extern int *NumberOfDihedralAngleHistogramDefinitions;
extern int (**DihedralAngleHistogramDefinitions)[4][3];
extern ATOM* (**DihedralAngleHistogramPairs)[4];

// sampling histograms of the angle between two planes (each formed by 3 chosen atoms)
void SampleAngleBetweenPlanesHistogram(int Switch);
extern int *ComputeAngleBetweenPlanesHistograms;
extern int *WriteAngleBetweenPlanesHistogramsEvery;
extern int NumberOfElementsAngleBetweenPlanesHistogram;
extern REAL MaxRangeAngleBetweenPlanesHistogram;
extern int *NumberOfAngleBetweenPlanesHistogramDefinitions;
extern int (**AngleBetweenPlanesHistogramDefinitions)[6][3];
extern ATOM* (**AngleBetweenPlanesHistogramPairs)[6];

// sampling molecular properties (bond distance, bend angle, dihedral angle)
void SampleMoleculePropertyHistogram(int);
extern int *ComputeMoleculeProperties;
extern int *WriteMoleculePropertiesEvery;
extern int *BondLengthHistogramSize;
extern REAL *BondLengthRange;
extern int *BendAngleHistogramSize;
extern REAL *BendAngleRange;
extern int *DihedralHistogramSize;
extern REAL *DihedralRange;

// sampling the IR spectra (spacings: 2048, 4196, 8192, 16384, 32768 points)
void SampleInfraRedSpectra(int Switch);
extern int *ComputeInfraRedSpectra;
extern int *WriteInfraRedSpectraEvery;
extern int SampleEveryInfraRed;

// sampling the mean-squared displacement using a modified order-N algorithm
void SampleMeanSquaredDisplacementOrderN(int);
extern int *ComputeMSDOrderN;               // whether or not to compute the msd
extern int *SampleMSDOrderNEvery;           // the sample frequency
extern int *WriteMSDOrderNEvery;            // write output every 'WriteMSDOrderNEvery' steps
extern int *NumberOfSitesMSDOrderN;         // the number of sites
extern int NumberOfBlockElementsMSDOrderN;  // the number of elements per block
extern int MaxNumberOfBlocksMSDOrderN;      // the maxmimum amount of blocks (data beyond this block is ignored)
extern int ComputeIndividualMSDOrderN;      // whether or not to compute (self-)msd's for individual molecules
extern int ComputeSiteTypeMSDOrderN;        // whether or not to compute (self-)msd's for individual molecules
extern int ComputeMSDOrderNPerPseudoAtom;   // whether or not to compute (self-)msd's for (pseudo-)atoms

// sampling the velocity autocorrelation function using a modified order-N algorithm
void SampleVelocityAutoCorrelationFunctionOrderN(int);
extern int *ComputeVACFOrderN;               // whether or not to compute the vacf
extern int *SampleVACFOrderNEvery;           // the sample frequency
extern int *WriteVACFOrderNEvery;            // write output every 'WriteVACFOrderNEvery' steps
extern int NumberOfBlockElementsVACFOrderN;  // the number of elements per block
extern int MaxNumberOfBlocksVACFOrderN;      // the maxmimum amount of blocks (data beyond this block is ignored)
extern int ComputeIndividualVACFOrderN;      // whether or not to compute (self-)vacf's for individual molecules
extern int ComputeVACFOrderNPerPseudoAtom;   // whether or not to compute (self-)vacf's for (pseudo-)atoms

// sampling of the rotational velocity autocorrelation function using a modified order-N algorithm
void SampleRotationalVelocityAutoCorrelationFunctionOrderN(int);
extern int *ComputeRVACFOrderN;               // whether or not to compute the vacf
extern int *SampleRVACFOrderNEvery;           // the sample frequency
extern int *WriteRVACFOrderNEvery;            // write output every 'WriteVACFOrderNEvery' steps
extern int NumberOfBlockElementsRVACFOrderN;  // the number of elements per block
extern int MaxNumberOfBlocksRVACFOrderN;      // the maxmimum amount of blocks (data beyond this block is ignored)
extern int ComputeIndividualRVACFOrderN;      // whether or not to compute (self-)vacf's for individual molecules
extern int ComputeRVACFOrderNPerPseudoAtom;   // whether or not to compute (self-)vacf's for (pseudo-)atoms

// sampling of the molecular orientation autocorrelation function using a modified order-N algorithm
enum {END_TO_END_VECTOR,MOLECULAR_VECTOR};
void SampleMolecularOrientationAutoCorrelationFunctionOrderN(int Switch);
extern int *ComputeMolecularOrientationOrderN;               // whether or not to compute the vacf
extern int *SampleMolecularOrientationOrderNEvery;           // the sample frequency
extern int *WriteMolecularOrientationOrderNEvery;            // write output every 'WriteVACFOrderNEvery' steps
extern int MolecularOrientationType;
extern VECTOR MolecularOrientationVector;
extern int MolecularOrientationGroup;
extern int NumberOfBlockElementsMolecularOrientationOrderN;  // the number of elements per block
extern int MaxNumberOfBlocksMolecularOrientationOrderN;      // the maxmimum amount of blocks (data beyond this block is ignored)
extern int ComputeIndividualMolecularOrientationOrderN;      // whether or not to compute (self-)vacf's for individual molecules
extern int ComputeMolecularOrientationOrderNPerPseudoAtom;   // whether or not to compute (self-)vacf's for (pseudo-)atoms

extern int **NumberOfOrientationFrameworkBonds;
extern char (***OrientationFrameworkBonds)[2][256];
extern int ***OrientationFrameworkBondTypes;

extern int ***NumberOfOrientationFrameworkBondPairs;
extern PAIR ****OrientationFrameworkBondPairs;
extern VECTOR ****BondOrientationAngleDistributionFunction;
extern int *BondOrientationAngleHistogramSize;



// sampling of the bond orientation autocorrelation function using a modified order-N algorithm
void SampleBondOrientationAutoCorrelationFunctionOrderN(int Switch);
extern int *ComputeBondOrientationOrderN;               // whether or not to compute the vacf
extern int *SampleBondOrientationOrderNEvery;           // the sample frequency
extern int *WriteBondOrientationOrderNEvery;            // write output every 'WriteBondOrientationOrderNEvery' steps
extern int NumberOfBlockElementsBondOrientationOrderN;  // the number of elements per block
extern int MaxNumberOfBlocksBondOrientationOrderN;      // the maxmimum amount of blocks (data beyond this block is ignored)
extern int ComputeIndividualBondOrientationOrderN;      // whether or not to compute (self-)vacf's for individual molecules
extern int ComputeBondOrientationOrderNPerPseudoAtom;   // whether or not to compute (self-)vacf's for (pseudo-)atoms

// sampling the mean-square displacement function using a conventional algorithm
void SampleMeanSquaredDisplacement(int Switch);
extern int *ComputeMSD;                      // whether or not to compute the msd
extern int *SampleMSDEvery;                  // the sample frequency
extern int *WriteMSDEvery;                   // write output every 'WriteMSDEvery' steps
extern int NumberOfBuffersMSD;               // the number of overlapping buffers
extern int BufferLengthMSD;                  // the length of the buffer-arrays

// sampling the velocity autocorrelation function using a conventional algorithm
void SampleVelocityAutoCorrelationFunction(int Switch);
extern int *ComputeVACF;                     // whether or not to compute the vacf
extern int *SampleVACFEvery;                 // the sample frequency
extern int *WriteVACFEvery;                  // write output every 'WriteVACFEvery' steps
extern int NumberOfBuffersVACF;              // the number of overlapping buffers
extern int BufferLengthVACF;                 // the length of the buffer-arrays

// sampling the 3D histograms of position (i.e. 3D free energy)
void SampleDensityProfile3DVTKGrid(int);
void SampleCOMDensityProfile3DVTKGrid(int); // testing
extern int *WriteDensityProfile3DVTKGridEvery;
extern int *ComputeDensityProfile3DVTKGrid;
extern INT_VECTOR3 DensityProfile3DVTKGridPoints;

// samples the cation sites and adsorption sites
void SampleCationAndAdsorptionSites(int);
extern int *ComputeCationAndAdsorptionSites;
extern int *WriteCationAndAdsorptionSitesEvery;

// samples initial configurations for the tranmission coefficient (dcTST)
void SampleDcTSTConfigurationFiles(int Choice);
extern int *WritedcTSTSnapShotsToFile;
extern int *WritedcTSTSnapShotsEvery;

// samples the principle moment of inertia
void MeasurePrincipleMomentsOfInertia(void);
extern int ComputePrincipleMomentsOfInertia;

void ComputeMolecularPressureTensor(REAL_MATRIX3x3*,REAL*,REAL*);
extern int *ComputeMolecularPressure;

void WriteRestartSample(FILE *FilePtr);
void AllocateSampleMemory(void);
void ReadRestartSample(FILE *FilePtr);

#endif
