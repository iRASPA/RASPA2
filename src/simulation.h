/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'simulation.h' is part of RASPA-2.0

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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "constants.h"
#include "utils.h"
#include "vector.h"

enum {UNINITIALIZED_BOUNDARY_CONDITION,FINITE,CUBIC,RECTANGULAR,TRICLINIC};
enum {NO_CHARGE,EWALD,SMOOTHED_COULOMB,WOLF_TEST,WOLFS_METHOD,WOLFS_METHOD2,WOLFS_METHOD_DAMPED,WOLFS_METHOD_DAMPED2,
      WOLFS_METHOD_DAMPED_FG,SHIFTED_COULOMB,TRUNCATED_COULOMB,SPME};
enum {NONE,FULL,RIGID,FLEXIBLE,GRID,FORCE_GRID};
enum {FREE_ENERGY,MAKE_GRID,MAKE_ASCI_GRID,MONTE_CARLO,MOLECULAR_DYNAMICS,SPECTRA,MINIMIZATION,GLOBAL_MINIMIZATION,VISUALIZATION,PORE_SIZE_DISTRIBUTION,SURFACE_AREA,
      BARRIER_RECROSSING,BARRIER_RUIZ_MONTERO,EWALD_TEST,TEST_SPACEGROUP,NUMERICAL,USER,TEST_WOLF,PLOT_POTENTIAL,STATUS,GENERATE_FRAMEWORK,PHONON_DISPERSION};
enum {OUTPUT_TO_FILE,OUTPUT_TO_SCREEN};
enum {BIGMAC_STYLE,RASPA_STYLE,RASPA_STYLE_OLD};
enum {NVE,NVT,NPT,NPH,MuPT,MuVT,NPTPR,NPHPR,MuPTPR,GC};
enum {REGULAR,MONOCLINIC,ISOTROPIC,ANISOTROPIC,REGULAR_UPPER_TRIANGLE,MONOCLINIC_UPPER_TRIANGLE};
enum {SYNPERIPLANAR,SYNCLINAL_PLUS,ANTICLINAL_PLUS,ANTIPERIPLANAR_PLUS,ANTICLINAL_MIN,SYNCLINAL_MIN};
enum {NONLINEAR_MOLECULE,LINEAR_MOLECULE,POINT_PARTICLE};
enum {NO_CHARILITY,S_CHIRAL,R_CHIRAL};
enum {FRAMEWORK,ADSORBATE,CATION};
enum {FIXED,FREE};
enum {FORWARD,BACKWARD,NO_FORWARD_OR_BACKWARD};
enum {ENANTIOFACE_RE,ENANTIOFACE_SI};
enum {POSITION_INITIALIZATION,VELOCITY_SCALING,VELOCITY_EQUILIBRATION,CF_WANG_LANDAU_EQUILIBRATION,PRODUCTION,FINISHED};
enum {CBMC,CFCMC,GIBBS_CFCMC_MAGINN,GIBBS_CFCMC_SINGLE_FRACTIONAL_PARTICLE};
enum {XYZ_DIR,XY_DIR,XZ_DIR,YZ_DIR,X_DIR,Y_DIR,Z_DIR,ABC_DIR,AB_DIR,AC_DIR,BC_DIR,A_DIR,B_DIR,C_DIR,
      ORTHOGONAL_TO_AB_DIR,ORTHOGONAL_TO_AC_DIR,ORTHOGONAL_TO_BC_DIR,
      ORTHOGONAL_TO_O_AB_DIR,ORTHOGONAL_TO_O_AC_DIR,ORTHOGONAL_TO_O_BC_DIR,
      ORTHOGONAL_TO_A_BC_DIR,ORTHOGONAL_TO_B_AC_DIR,ORTHOGONAL_TO_C_AB_DIR,ORTHOGONAL_TO_O_ABC_DIR};


enum {MOL_FORMAT,CSSR_FORMAT,DLPOLY_FORMAT,XYZ_FORMAT,PDB_FORMAT,CIF_FORMAT};

enum {CHARGE_ATOM_FROM_PSEUDO_ATOM_DEFINITION,CHARGE_ATOM_FROM_STRUCTURE_FILE,CHARGE_NOT_DEFINED};


enum MC_ENSEMBLE
{
  CONSTANT_NUMBER_OF_FRAMEWORK_PARTICLES   = 0,   //binary value == 00000000
  CONSTANT_NUMBER_OF_PARTICLES             = 1,   //binary value == 00000001
  CONSTANT_CHEMICAL_POTENTIAL              = 2,   //binary value == 00000010
  CONSTANT_VOLUME                          = 4,   //binary value == 00000100
  CONSTANT_PRESSURE                        = 8,   //binary value == 00001000
  CONSTANT_TENSION                         = 16,  //binary value == 00010000
  GIBBS                                    = 32,  //binary value == 00100000
  CONSTANT_TEMPERATURE                     = 64,  //binary value == 01000000
  CONSTANT_ENERGY                          = 128  //binary value == 10000000
};

extern int (**CellListMap)[27];

extern int Dimension;

extern unsigned long seed;
extern char *RASPA_DIRECTORY;
extern char FileNameAppend[256];
extern char ForceField[256];

extern int SimulationType;
extern int SimulationSubType;

extern int NumberOfSystems;                 // the number of systems
extern int CurrentSystem;                   // the current system

//----------------------------------------------------------------------------------------
// CFC-RXMC Parameters
//----------------------------------------------------------------------------------------

extern int NumberOfReactions;                       // Total number of Reactions
extern REAL **CFCRXMCLambda;                        // Reaction Lambda for all reactions
extern REAL ProbabilityCFCRXMCLambdaChangeMove;
extern int **ReactantsStoichiometry;                // Reactants Stoichiometry
extern int **ProductsStoichiometry;                 // Products Stoichiometry
extern int RXMCLambdaHistogramSize;
extern REAL ***RXMCBiasingFactors;
extern REAL **CFRXMCWangLandauScalingFactor;

//----------------------------------------------------------------------------------------


extern int NumberOfPartialPressures;
extern int CurrentPartialPressure;

extern long long NumberOfCycles;
extern long long NumberOfVelocityScalingCycles;
extern long long NumberOfEquilibrationCycles;
extern long long NumberOfInitializationCycles;
extern long long CurrentCycle;
extern int SimulationStage;
extern int PrintEvery;
extern int OptimizeAcceptenceEvery;
extern int PrintPropertiesEvery;
extern int Output;
extern int WriteForcefieldToOutput;
extern int WritePseudoAtomsToOutput;
extern int WriteMoleculeDefinitionToOutput;
extern int InputFileType;
extern int ContinueAfterCrash;
extern int WriteBinaryRestartFileEvery;

extern int OmitInterMolecularInteractions;
extern int OmitInterMolecularVDWInteractions;
extern int OmitInterMolecularCoulombInteractions;
extern int OmitAdsorbateAdsorbateVDWInteractions;
extern int OmitAdsorbateAdsorbateCoulombInteractions;
extern int OmitAdsorbateCationVDWInteractions;
extern int OmitAdsorbateCationCoulombInteractions;
extern int OmitCationCationVDWInteractions;
extern int OmitCationCationCoulombInteractions;

extern int OmitAdsorbateAdsorbatePolarization;
extern int OmitCationCationPolarization;
extern int OmitAdsorbateCationPolarization;
extern int OmitIntraFrameworkPolarization;

extern int ComputeStress;
extern int ComputeHessian;
extern int ComputeBornTerm;
extern int ComputeElasticConstants;
extern int ComputeElasticConstantsThirdOrder;
extern int ComputeCrossTerm;

extern int ComputePolarization;
extern int BackPolarization;
extern int NumberOfBackPolarizationSteps;
extern int PolarizationMatrix;                      // type: isotropic, anisotropic, full

extern int ComputePrincipleMomentsOfInertia;

extern int ReinitializeVelocities;
extern int QuenchCoreShellVelocities;
extern int QuenchCoreShellVelocitiesEvery;

extern int MinimumInnerCycles;

extern REAL OverlapDistanceSquared;
extern REAL CutOffVDW;
extern REAL CutOffVDWSquared;
extern REAL CutOffVDWSwitch;
extern REAL CutOffVDWSwitchSquared;
extern REAL InverseCutOffVDW;

extern REAL *CutOffChargeCharge;
extern REAL *CutOffChargeChargeSquared;
extern REAL *CutOffChargeChargeSwitch;
extern REAL *CutOffChargeChargeSwitchSquared;
extern REAL *InverseCutOffChargeCharge;

extern REAL CutOffChargeBondDipole;
extern REAL CutOffChargeBondDipoleSquared;
extern REAL CutOffChargeBondDipoleSwitch;
extern REAL CutOffChargeBondDipoleSwitchSquared;

extern REAL CutOffBondDipoleBondDipole;
extern REAL CutOffBondDipoleBondDipoleSquared;
extern REAL CutOffBondDipoleBondDipoleSwitch;
extern REAL CutOffBondDipoleBondDipoleSwitchSquared;


extern REAL DeltaT;

extern int ChargeMethod;
extern int ChargeFromChargeEquilibration;
extern int ChargeEquilibrationMethod;
extern int SymmetrizeFrameworkCharges;
extern int UseChargesFromMOLFile;
extern int UseChargesFromCIFFile;

extern REAL_MATRIX3x3 *Box;
extern REAL_MATRIX3x3 *InverseBox;
extern REAL_MATRIX3x3 *ReplicaBox;
extern REAL_MATRIX3x3 *InverseReplicaBox;
extern INT_VECTOR3 *NumberOfReplicaCells;
extern int *TotalNumberOfReplicaCells;
extern VECTOR *ReplicaShift;
extern int *UseReplicas;
extern REAL_MATRIX3x3 *BoxProperties;
extern REAL_MATRIX3x3 *InverseBoxProperties;
extern REAL *Volume;
extern REAL *AlphaAngle;
extern REAL *BetaAngle;
extern REAL *GammaAngle;
extern int *BoundaryCondition;
extern VECTOR *IntialCenterOfMassPosition;
extern REAL *Beta;

extern INT_VECTOR3 *NumberOfCellListCells;
extern int *UseCellLists;
extern int *NumberOfCellLists;

// energies
extern REAL *UTotal;
extern REAL *ConservedEnergy;
extern REAL *Drift;
extern REAL *ReferenceEnergy;

extern REAL *UIon;
extern REAL *UDrift;
extern REAL *UTailCorrection;

extern REAL *UDistanceConstraints;
extern REAL *UAngleConstraints;
extern REAL *UDihedralConstraints;
extern REAL *UInversionBendConstraints;
extern REAL *UOutOfPlaneDistanceConstraints;
extern REAL *UExclusionConstraints;

extern REAL *UHostPolarization;
extern REAL *UAdsorbatePolarization;
extern REAL *UCationPolarization;

extern REAL *UHostBackPolarization;
extern REAL *UAdsorbateBackPolarization;
extern REAL *UCationBackPolarization;

extern REAL *UHostHost;
extern REAL *UAdsorbateAdsorbate;
extern REAL *UCationCation;
extern REAL *UHostCation;
extern REAL *UHostAdsorbate;
extern REAL *UAdsorbateCation;

extern REAL *UHostHostCoulomb;
extern REAL *UHostHostVDW;
extern REAL *UHostAdsorbateCoulomb;
extern REAL *UHostAdsorbateVDW;
extern REAL *UHostCationCoulomb;
extern REAL *UHostCationVDW;
extern REAL *UAdsorbateAdsorbateCoulomb;
extern REAL *UAdsorbateAdsorbateVDW;
extern REAL *UAdsorbateCationCoulomb;
extern REAL *UAdsorbateCationVDW;
extern REAL *UCationCationCoulomb;
extern REAL *UCationCationVDW;

extern REAL *UHostHostChargeChargeReal;
extern REAL *UHostAdsorbateChargeChargeReal;
extern REAL *UHostCationChargeChargeReal;
extern REAL *UAdsorbateAdsorbateChargeChargeReal;
extern REAL *UAdsorbateCationChargeChargeReal;
extern REAL *UCationCationChargeChargeReal;

extern REAL *UHostHostChargeBondDipoleReal;
extern REAL *UHostAdsorbateChargeBondDipoleReal;
extern REAL *UHostCationChargeBondDipoleReal;
extern REAL *UAdsorbateAdsorbateChargeBondDipoleReal;
extern REAL *UAdsorbateCationChargeBondDipoleReal;
extern REAL *UCationCationChargeBondDipoleReal;

extern REAL *UHostHostBondDipoleBondDipoleReal;
extern REAL *UHostAdsorbateBondDipoleBondDipoleReal;
extern REAL *UHostCationBondDipoleBondDipoleReal;
extern REAL *UAdsorbateAdsorbateBondDipoleBondDipoleReal;
extern REAL *UAdsorbateCationBondDipoleBondDipoleReal;
extern REAL *UCationCationBondDipoleBondDipoleReal;

extern REAL *UHostHostChargeChargeFourier;
extern REAL *UHostAdsorbateChargeChargeFourier;
extern REAL *UHostCationChargeChargeFourier;
extern REAL *UAdsorbateAdsorbateChargeChargeFourier;
extern REAL *UAdsorbateCationChargeChargeFourier;
extern REAL *UCationCationChargeChargeFourier;

extern REAL *UHostHostChargeBondDipoleFourier;
extern REAL *UHostAdsorbateChargeBondDipoleFourier;
extern REAL *UHostCationChargeBondDipoleFourier;
extern REAL *UAdsorbateAdsorbateChargeBondDipoleFourier;
extern REAL *UAdsorbateCationChargeBondDipoleFourier;
extern REAL *UCationCationChargeBondDipoleFourier;

extern REAL *UHostHostBondDipoleBondDipoleFourier;
extern REAL *UHostAdsorbateBondDipoleBondDipoleFourier;
extern REAL *UHostCationBondDipoleBondDipoleFourier;
extern REAL *UAdsorbateAdsorbateBondDipoleBondDipoleFourier;
extern REAL *UAdsorbateCationBondDipoleBondDipoleFourier;
extern REAL *UCationCationBondDipoleBondDipoleFourier;

extern REAL *UHostBond;
extern REAL *UHostUreyBradley;
extern REAL *UHostBend;
extern REAL *UHostInversionBend;
extern REAL *UHostTorsion;
extern REAL *UHostImproperTorsion;
extern REAL *UHostOutOfPlane;
extern REAL *UHostBondBond;
extern REAL *UHostBondBend;
extern REAL *UHostBendBend;
extern REAL *UHostBondTorsion;
extern REAL *UHostBendTorsion;

extern REAL *UCationBond;
extern REAL *UCationUreyBradley;
extern REAL *UCationBend;
extern REAL *UCationInversionBend;
extern REAL *UCationTorsion;
extern REAL *UCationImproperTorsion;
extern REAL *UCationOutOfPlane;
extern REAL *UCationBondBond;
extern REAL *UCationBondBend;
extern REAL *UCationBendBend;
extern REAL *UCationBondTorsion;
extern REAL *UCationBendTorsion;
extern REAL *UCationIntraVDW;
extern REAL *UCationIntraChargeCharge;
extern REAL *UCationIntraChargeBondDipole;
extern REAL *UCationIntraBondDipoleBondDipole;

extern REAL *UAdsorbateBond;
extern REAL *UAdsorbateUreyBradley;
extern REAL *UAdsorbateBend;
extern REAL *UAdsorbateInversionBend;
extern REAL *UAdsorbateTorsion;
extern REAL *UAdsorbateImproperTorsion;
extern REAL *UAdsorbateOutOfPlane;
extern REAL *UAdsorbateBondBond;
extern REAL *UAdsorbateBondBend;
extern REAL *UAdsorbateBendBend;
extern REAL *UAdsorbateBondTorsion;
extern REAL *UAdsorbateBendTorsion;
extern REAL *UAdsorbateIntraVDW;
extern REAL *UAdsorbateIntraChargeCharge;
extern REAL *UAdsorbateIntraChargeBondDipole;
extern REAL *UAdsorbateIntraBondDipoleBondDipole;


// MD

extern REAL *PressureTrace;
extern REAL *Pressure;

extern REAL *UKinetic;
extern REAL *UHostKinetic;
extern REAL *UCationKinetic;
extern REAL *UAdsorbateKinetic;

extern REAL *UAdsorbateRotationalKinetic;
extern REAL *UAdsorbateTranslationalKinetic;
extern REAL *UCationRotationalKinetic;
extern REAL *UCationTranslationalKinetic;

extern REAL *UNoseHoover;
extern REAL *UNoseHooverAdsorbates;
extern REAL *UNoseHooverCations;
extern REAL *UNoseHooverFramework;

extern REAL_MATRIX3x3 *StressTensor;
extern REAL_MATRIX3x3 *ConfigurationalStressTensor;
extern REAL_MATRIX3x3 *KineticStressTensor;
extern REAL_MATRIX3x3 *StrainDerivativeTensor;
extern REAL *StrainDerivativeTailCorrection;
extern REAL_MATRIX3x3 *StrainDerivativeRigidCorrection;
extern REAL_MATRIX3x3 *MolecularStressTensor;

extern REAL_MATRIX9x9 *BornTerm;
extern REAL_MATRIX9x9 *RelaxationTerm;

enum {VERLET,GGMT,XI_RESPA_MTS,XO_RESPA_MTS};        // integration-scheme
extern int *Ensemble;
extern int *InitEnsemble;
extern int *RunEnsemble;
extern int *NPTPRCellType;
extern int *MonoclinicAngleType;
extern int *MCEnsemble;
enum {MONOCLINIC_ALPHA_ANGLE,MONOCLINIC_BETA_ANGLE,MONOCLINIC_GAMMA_ANGLE};

extern int *DegreesOfFreedom;
extern int *DegreesOfFreedomTranslation;
extern int *DegreesOfFreedomRotation;
extern int *DegreesOfFreedomVibration;
extern int *DegreesOfFreedomConstraint;

extern int *DegreesOfFreedomFramework;

extern int *DegreesOfFreedomAdsorbates;
extern int *DegreesOfFreedomTranslationalAdsorbates;
extern int *DegreesOfFreedomRotationalAdsorbates;
extern int *DegreesOfFreedomVibrationalAdsorbates;
extern int *DegreesOfFreedomConstraintAdsorbates;

extern int *DegreesOfFreedomCations;
extern int *DegreesOfFreedomTranslationalCations;
extern int *DegreesOfFreedomRotationalCations;
extern int *DegreesOfFreedomVibrationalCations;
extern int *DegreesOfFreedomConstraintCations;

extern REAL ProbabilityParallelTemperingMove;
extern REAL ProbabilityHyperParallelTemperingMove;
extern REAL ProbabilityParallelMolFractionMove;
extern REAL ProbabilityChiralInversionMove;
extern REAL ProbabilityHybridNVEMove;
extern REAL ProbabilityHybridNPHMove;
extern REAL ProbabilityHybridNPHPRMove;
extern REAL ProbabilityVolumeChangeMove;
extern REAL ProbabilityBoxShapeChangeMove;
extern REAL ProbabilityGibbsVolumeChangeMove;
extern REAL ProbabilityFrameworkChangeMove;
extern REAL ProbabilityFrameworkShiftMove;

extern REAL CpuTimeProductionRun;
extern REAL CpuTimeInitialization;
extern REAL CpuTimeEquilibration;
extern REAL CpuTotal;

extern REAL *CpuTimeParallelTemperingMove;
extern REAL *CpuTimeHyperParallelTemperingMove;
extern REAL *CpuTimeParallelMolFractionMove;
extern REAL *CpuTimeChiralInversionMove;
extern REAL *CpuTimeHybridNVEMove;
extern REAL *CpuTimeHybridNPHMove;
extern REAL *CpuTimeHybridNPHPRMove;
extern REAL *CpuTimeVolumeChangeMove;
extern REAL *CpuTimeBoxShapeChangeMove;
extern REAL *CpuTimeGibbsVolumeChangeMove;
extern REAL *CpuTimeFrameworkChangeMove;
extern REAL *CpuTimeFrameworkShiftMove;
extern REAL *CpuTimeCFCRXMCLambdaChangeMove;

extern int OptimizeVolumeChange;
extern int OptimizeGibbsVolumeChange;
extern int OptimizeTranslation;
extern int OptimizeRotation;
extern int OptimizeFrameworkChange;
extern int OptimizeFrameworkShift;
extern int OptimizeCFLambdaChange;
extern int OptimizeCFGibbsLambdaChange;
extern int OptimizeCBCFLambdaChange;
extern int OptimizeCBCFGibbsLambdaChange;
extern int OptimizeRXMCLambdaChange;

void ScaleBornTerm(REAL r);
void AddContributionToCrossTerm(int i,REAL_MATRIX CrossTerm,REAL DDF,REAL DF,VECTOR dr);
void AddContributionToBornTerm(REAL DDF,REAL DF,VECTOR dr);
void FillInRemainingTermsInBornTerm(void);

void InitializeReplicaBox(void);
void InverseBoxMatrix(REAL_MATRIX3x3 *Box,REAL_MATRIX3x3 *InverseBox);

void WriteRestartSimulation(FILE *FilePtr);
void AllocateSimulationMemory(void);
void ReadRestartSimulation(FILE *FilePtr);

#endif
