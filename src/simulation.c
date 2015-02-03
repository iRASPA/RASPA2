/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2015 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'simulation.c' is part of RASPA-2.0

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "simulation.h"
#include "framework.h"
#include "framework_energy.h"
#include "molecule.h"
#include "recrossing.h"
#include "output.h"
#include "movies.h"
#include "utils.h"

int (**CellListMap)[27];

int Dimension;

unsigned long seed;
char *RASPA_DIRECTORY;
char FileNameAppend[256];
char ForceField[256];

int SimulationType;

int NumberOfSystems;                 // the number of systems
int CurrentSystem;                   // the current system

//----------------------------------------------------------------------------------------
// CFC-RXMC Parameters
//----------------------------------------------------------------------------------------

int NumberOfReactions;                       // Total number of Reactions
REAL **CFCRXMCLambda;                        // Reaction Lambda for all reactions
REAL ProbabilityCFCRXMCLambdaChangeMove;
int **ReactantsStoichiometry;                // Reactants Stoichiometry
int **ProductsStoichiometry;                 // Products Stoichiometry
int RXMCLambdaHistogramSize;
REAL ***RXMCBiasingFactors;
REAL **CFRXMCWangLandauScalingFactor;

//----------------------------------------------------------------------------------------

int NumberOfPartialPressures;
int CurrentPartialPressure;


long long NumberOfCycles;
long long NumberOfVelocityScalingCycles;
long long NumberOfEquilibrationCycles;
long long NumberOfInitializationCycles;
long long CurrentCycle;
int SimulationStage;
int PrintEvery;
int OptimizeAcceptenceEvery;
int PrintPropertiesEvery;
int Output;
int WriteForcefieldToOutput;
int WritePseudoAtomsToOutput;
int WriteMoleculeDefinitionToOutput;
int InputFileType;
int ContinueAfterCrash;
int WriteBinaryRestartFileEvery;


int OmitInterMolecularInteractions;
int OmitInterMolecularVDWInteractions;
int OmitInterMolecularCoulombInteractions;
int OmitAdsorbateAdsorbateVDWInteractions;
int OmitAdsorbateAdsorbateCoulombInteractions;
int OmitAdsorbateCationVDWInteractions;
int OmitAdsorbateCationCoulombInteractions;
int OmitCationCationVDWInteractions;
int OmitCationCationCoulombInteractions;

int OmitAdsorbateAdsorbatePolarization;
int OmitCationCationPolarization;
int OmitAdsorbateCationPolarization;
int OmitIntraFrameworkPolarization;

int ComputeStress;
int ComputeHessian;
int ComputeBornTerm;
int ComputeElasticConstants;
int ComputeElasticConstantsThirdOrder;
int ComputeCrossTerm;

int ComputePolarization;
int BackPolarization;
int NumberOfBackPolarizationSteps;
int PolarizationMatrix;                      // type: isotropic, anisotropic, full

int ReinitializeVelocities;
int QuenchCoreShellVelocities;
int QuenchCoreShellVelocitiesEvery;

int MinimumInnerCycles;

int ChargeMethod;
int ChargeFromChargeEquilibration;
int ChargeEquilibrationMethod;
int SymmetrizeFrameworkCharges;
int UseChargesFromMOLFile;
int UseChargesFromCIFFile;

REAL OverlapDistanceSquared;
REAL CutOffVDW;
REAL CutOffVDWSquared;
REAL CutOffVDWSwitch;
REAL CutOffVDWSwitchSquared;
REAL InverseCutOffVDW;

REAL *CutOffChargeCharge;
REAL *CutOffChargeChargeSquared;
REAL *CutOffChargeChargeSwitch;
REAL *CutOffChargeChargeSwitchSquared;
REAL *InverseCutOffChargeCharge;

REAL CutOffChargeBondDipole;
REAL CutOffChargeBondDipoleSquared;
REAL CutOffChargeBondDipoleSwitch;
REAL CutOffChargeBondDipoleSwitchSquared;

REAL CutOffBondDipoleBondDipole;
REAL CutOffBondDipoleBondDipoleSquared;
REAL CutOffBondDipoleBondDipoleSwitch;
REAL CutOffBondDipoleBondDipoleSwitchSquared;

REAL DeltaT;

REAL_MATRIX3x3 *Box;                   // the cell matrix
REAL_MATRIX3x3 *InverseBox;            // the inverse of the cell matrix
REAL_MATRIX3x3 *ReplicaBox;            // the cell matrix of the replica system
REAL_MATRIX3x3 *InverseReplicaBox;     // the inverse of the the cell matrix of the replica system
INT_VECTOR3 *NumberOfReplicaCells;     // the integere number of replicas in each direction a,b,c
int *TotalNumberOfReplicaCells;        // the total number of replica cells
VECTOR *ReplicaShift;                  // the shift in a,b,c for each replica cell
int *UseReplicas;                      // whether or not to use replicas
REAL_MATRIX3x3 *BoxProperties;         // properties of the cell matrix (i.e. perpendicular lengths)
REAL_MATRIX3x3 *InverseBoxProperties;  // properties of the inverse cell matrix
REAL *Volume;                          // the volume
REAL *AlphaAngle;                      // the alpha-angle of the cell
REAL *BetaAngle;                       // the beta-angle of the cell
REAL *GammaAngle;                      // the gamma-angle of the cell
int *BoundaryCondition;                // the boundary condition (i.e. 'RECTANGULAR' or 'TRICLINIC')
VECTOR *IntialCenterOfMassPosition;    // the intial position of the center of mass of the system

REAL *Beta;

INT_VECTOR3 *NumberOfCellListCells;
int *UseCellLists;
int *NumberOfCellLists;

// energies
REAL *UTotal;
REAL *ConservedEnergy;
REAL *Drift;
REAL *ReferenceEnergy;

REAL *UIon;
REAL *UDrift;
REAL *UTailCorrection;

REAL *UDistanceConstraints;
REAL *UAngleConstraints;
REAL *UDihedralConstraints;
REAL *UInversionBendConstraints;
REAL *UOutOfPlaneDistanceConstraints;
REAL *UExclusionConstraints;

REAL *UHostPolarization;
REAL *UAdsorbatePolarization;
REAL *UCationPolarization;

REAL *UHostBackPolarization;
REAL *UAdsorbateBackPolarization;
REAL *UCationBackPolarization;

REAL *UHostHost;
REAL *UAdsorbateAdsorbate;
REAL *UCationCation;
REAL *UHostAdsorbate;
REAL *UHostCation;
REAL *UAdsorbateCation;

REAL *UHostHostCoulomb;
REAL *UHostHostVDW;
REAL *UHostAdsorbateCoulomb;
REAL *UHostAdsorbateVDW;
REAL *UHostCationCoulomb;
REAL *UHostCationVDW;
REAL *UAdsorbateAdsorbateCoulomb;
REAL *UAdsorbateAdsorbateVDW;
REAL *UAdsorbateCationCoulomb;
REAL *UAdsorbateCationVDW;
REAL *UCationCationCoulomb;
REAL *UCationCationVDW;

REAL *UHostHostChargeChargeReal;
REAL *UHostAdsorbateChargeChargeReal;
REAL *UHostCationChargeChargeReal;
REAL *UAdsorbateAdsorbateChargeChargeReal;
REAL *UAdsorbateCationChargeChargeReal;
REAL *UCationCationChargeChargeReal;

REAL *UHostHostChargeBondDipoleReal;
REAL *UHostAdsorbateChargeBondDipoleReal;
REAL *UHostCationChargeBondDipoleReal;
REAL *UAdsorbateAdsorbateChargeBondDipoleReal;
REAL *UAdsorbateCationChargeBondDipoleReal;
REAL *UCationCationChargeBondDipoleReal;

REAL *UHostHostBondDipoleBondDipoleReal;
REAL *UHostAdsorbateBondDipoleBondDipoleReal;
REAL *UHostCationBondDipoleBondDipoleReal;
REAL *UAdsorbateAdsorbateBondDipoleBondDipoleReal;
REAL *UAdsorbateCationBondDipoleBondDipoleReal;
REAL *UCationCationBondDipoleBondDipoleReal;

REAL *UHostHostChargeChargeFourier;
REAL *UHostAdsorbateChargeChargeFourier;
REAL *UHostCationChargeChargeFourier;
REAL *UAdsorbateAdsorbateChargeChargeFourier;
REAL *UAdsorbateCationChargeChargeFourier;
REAL *UCationCationChargeChargeFourier;

REAL *UHostHostChargeBondDipoleFourier;
REAL *UHostAdsorbateChargeBondDipoleFourier;
REAL *UHostCationChargeBondDipoleFourier;
REAL *UAdsorbateAdsorbateChargeBondDipoleFourier;
REAL *UAdsorbateCationChargeBondDipoleFourier;
REAL *UCationCationChargeBondDipoleFourier;

REAL *UHostHostBondDipoleBondDipoleFourier;
REAL *UHostAdsorbateBondDipoleBondDipoleFourier;
REAL *UHostCationBondDipoleBondDipoleFourier;
REAL *UAdsorbateAdsorbateBondDipoleBondDipoleFourier;
REAL *UAdsorbateCationBondDipoleBondDipoleFourier;
REAL *UCationCationBondDipoleBondDipoleFourier;

REAL *UHostBond;
REAL *UHostUreyBradley;
REAL *UHostBend;
REAL *UHostInversionBend;
REAL *UHostTorsion;
REAL *UHostImproperTorsion;
REAL *UHostOutOfPlane;
REAL *UHostBondBond;
REAL *UHostBondBend;
REAL *UHostBendBend;
REAL *UHostBondTorsion;
REAL *UHostBendTorsion;

REAL *UCationBond;
REAL *UCationUreyBradley;
REAL *UCationBend;
REAL *UCationInversionBend;
REAL *UCationTorsion;
REAL *UCationImproperTorsion;
REAL *UCationOutOfPlane;
REAL *UCationBondBond;
REAL *UCationBondBend;
REAL *UCationBendBend;
REAL *UCationBondTorsion;
REAL *UCationBendTorsion;
REAL *UCationIntraVDW;
REAL *UCationIntraChargeCharge;
REAL *UCationIntraChargeBondDipole;
REAL *UCationIntraBondDipoleBondDipole;

REAL *UAdsorbateBond;
REAL *UAdsorbateUreyBradley;
REAL *UAdsorbateBend;
REAL *UAdsorbateInversionBend;
REAL *UAdsorbateTorsion;
REAL *UAdsorbateImproperTorsion;
REAL *UAdsorbateOutOfPlane;
REAL *UAdsorbateBondBond;
REAL *UAdsorbateBondBend;
REAL *UAdsorbateBendBend;
REAL *UAdsorbateBondTorsion;
REAL *UAdsorbateBendTorsion;
REAL *UAdsorbateIntraVDW;
REAL *UAdsorbateIntraChargeCharge;
REAL *UAdsorbateIntraChargeBondDipole;
REAL *UAdsorbateIntraBondDipoleBondDipole;

// MD

int *Ensemble;
int *InitEnsemble;
int *RunEnsemble;
int *NPTPRCellType;
int *MonoclinicAngleType;
int *MCEnsemble;

int *DegreesOfFreedom;
int *DegreesOfFreedomTranslation;
int *DegreesOfFreedomRotation;
int *DegreesOfFreedomVibration;
int *DegreesOfFreedomConstraint;

int *DegreesOfFreedomFramework;

int *DegreesOfFreedomAdsorbates;
int *DegreesOfFreedomTranslationalAdsorbates;
int *DegreesOfFreedomRotationalAdsorbates;
int *DegreesOfFreedomVibrationalAdsorbates;
int *DegreesOfFreedomConstraintAdsorbates;

int *DegreesOfFreedomCations;
int *DegreesOfFreedomTranslationalCations;
int *DegreesOfFreedomRotationalCations;
int *DegreesOfFreedomVibrationalCations;
int *DegreesOfFreedomConstraintCations;

REAL *UKinetic;
REAL *UHostKinetic;
REAL *UCationKinetic;
REAL *UAdsorbateKinetic;

REAL *UAdsorbateRotationalKinetic;
REAL *UAdsorbateTranslationalKinetic;
REAL *UCationRotationalKinetic;
REAL *UCationTranslationalKinetic;

REAL *UNoseHoover;
REAL *UNoseHooverAdsorbates;
REAL *UNoseHooverCations;
REAL *UNoseHooverFramework;

REAL *PressureTrace;
REAL *Pressure;

REAL_MATRIX3x3 *StressTensor;
REAL_MATRIX3x3 *ConfigurationalStressTensor;
REAL_MATRIX3x3 *KineticStressTensor;
REAL_MATRIX3x3 *StrainDerivativeTensor;
REAL *StrainDerivativeTailCorrection;
REAL_MATRIX3x3 *StrainDerivativeRigidCorrection;
REAL_MATRIX3x3 *MolecularStressTensor;

REAL_MATRIX9x9 *BornTerm;
REAL_MATRIX9x9 *RelaxationTerm;

REAL ProbabilityParallelTemperingMove;
REAL ProbabilityHyperParallelTemperingMove;
REAL ProbabilityParallelMolFractionMove;
REAL ProbabilityChiralInversionMove;
REAL ProbabilityHybridNVEMove;
REAL ProbabilityHybridNPHMove;
REAL ProbabilityHybridNPHPRMove;
REAL ProbabilityVolumeChangeMove;
REAL ProbabilityBoxShapeChangeMove;
REAL ProbabilityGibbsVolumeChangeMove;
REAL ProbabilityFrameworkChangeMove;
REAL ProbabilityFrameworkShiftMove;

REAL CpuTimeProductionRun;
REAL CpuTimeInitialization;
REAL CpuTimeEquilibration;
REAL CpuTotal;

REAL *CpuTimeParallelTemperingMove;
REAL *CpuTimeHyperParallelTemperingMove;
REAL *CpuTimeParallelMolFractionMove;
REAL *CpuTimeChiralInversionMove;
REAL *CpuTimeHybridNVEMove;
REAL *CpuTimeHybridNPHMove;
REAL *CpuTimeHybridNPHPRMove;
REAL *CpuTimeVolumeChangeMove;
REAL *CpuTimeBoxShapeChangeMove;
REAL *CpuTimeGibbsVolumeChangeMove;
REAL *CpuTimeFrameworkChangeMove;
REAL *CpuTimeFrameworkShiftMove;
REAL *CpuCFCRXMCLambdaChangeMove;

int OptimizeVolumeChange;
int OptimizeGibbsVolumeChange;
int OptimizeTranslation;
int OptimizeRotation;
int OptimizeFrameworkChange;
int OptimizeFrameworkShift;
int OptimizeCFLambdaChange;
int OptimizeCFGibbsLambdaChange;
int OptimizeCBCFLambdaChange;
int OptimizeCBCFGibbsLambdaChange;
int OptimizeRXMCLambdaChange;

void ScaleBornTerm(REAL r)
{
  BornTerm[CurrentSystem].xxxx*=r;
  BornTerm[CurrentSystem].yxxx*=r;
  BornTerm[CurrentSystem].zxxx*=r;
  BornTerm[CurrentSystem].xxyx*=r;
  BornTerm[CurrentSystem].yxyx*=r;
  BornTerm[CurrentSystem].zxyx*=r;
  BornTerm[CurrentSystem].xxzx*=r;
  BornTerm[CurrentSystem].yxzx*=r;
  BornTerm[CurrentSystem].zxzx*=r;

  BornTerm[CurrentSystem].xyxx*=r;
  BornTerm[CurrentSystem].yyxx*=r;
  BornTerm[CurrentSystem].zyxx*=r;
  BornTerm[CurrentSystem].xyyx*=r;
  BornTerm[CurrentSystem].yyyx*=r;
  BornTerm[CurrentSystem].zyyx*=r;
  BornTerm[CurrentSystem].xyzx*=r;
  BornTerm[CurrentSystem].yyzx*=r;
  BornTerm[CurrentSystem].zyzx*=r;

  BornTerm[CurrentSystem].xzxx*=r;
  BornTerm[CurrentSystem].yzxx*=r;
  BornTerm[CurrentSystem].zzxx*=r;
  BornTerm[CurrentSystem].xzyx*=r;
  BornTerm[CurrentSystem].yzyx*=r;
  BornTerm[CurrentSystem].zzyx*=r;
  BornTerm[CurrentSystem].xzzx*=r;
  BornTerm[CurrentSystem].yzzx*=r;
  BornTerm[CurrentSystem].zzzx*=r;

  BornTerm[CurrentSystem].xxxy*=r;
  BornTerm[CurrentSystem].yxxy*=r;
  BornTerm[CurrentSystem].zxxy*=r;
  BornTerm[CurrentSystem].xxyy*=r;
  BornTerm[CurrentSystem].yxyy*=r;
  BornTerm[CurrentSystem].zxyy*=r;
  BornTerm[CurrentSystem].xxzy*=r;
  BornTerm[CurrentSystem].yxzy*=r;
  BornTerm[CurrentSystem].zxzy*=r;

  BornTerm[CurrentSystem].xyxy*=r;
  BornTerm[CurrentSystem].yyxy*=r;
  BornTerm[CurrentSystem].zyxy*=r;
  BornTerm[CurrentSystem].xyyy*=r;
  BornTerm[CurrentSystem].yyyy*=r;
  BornTerm[CurrentSystem].zyyy*=r;
  BornTerm[CurrentSystem].xyzy*=r;
  BornTerm[CurrentSystem].yyzy*=r;
  BornTerm[CurrentSystem].zyzy*=r;

  BornTerm[CurrentSystem].xzxy*=r;
  BornTerm[CurrentSystem].yzxy*=r;
  BornTerm[CurrentSystem].zzxy*=r;
  BornTerm[CurrentSystem].xzyy*=r;
  BornTerm[CurrentSystem].yzyy*=r;
  BornTerm[CurrentSystem].zzyy*=r;
  BornTerm[CurrentSystem].xzzy*=r;
  BornTerm[CurrentSystem].yzzy*=r;
  BornTerm[CurrentSystem].zzzy*=r;

  BornTerm[CurrentSystem].xxxz*=r;
  BornTerm[CurrentSystem].yxxz*=r;
  BornTerm[CurrentSystem].zxxz*=r;
  BornTerm[CurrentSystem].xxyz*=r;
  BornTerm[CurrentSystem].yxyz*=r;
  BornTerm[CurrentSystem].zxyz*=r;
  BornTerm[CurrentSystem].xxzz*=r;
  BornTerm[CurrentSystem].yxzz*=r;
  BornTerm[CurrentSystem].zxzz*=r;

  BornTerm[CurrentSystem].xyxz*=r;
  BornTerm[CurrentSystem].yyxz*=r;
  BornTerm[CurrentSystem].zyxz*=r;
  BornTerm[CurrentSystem].xyyz*=r;
  BornTerm[CurrentSystem].yyyz*=r;
  BornTerm[CurrentSystem].zyyz*=r;
  BornTerm[CurrentSystem].xyzz*=r;
  BornTerm[CurrentSystem].yyzz*=r;
  BornTerm[CurrentSystem].zyzz*=r;

  BornTerm[CurrentSystem].xzxz*=r;
  BornTerm[CurrentSystem].yzxz*=r;
  BornTerm[CurrentSystem].zzxz*=r;
  BornTerm[CurrentSystem].xzyz*=r;
  BornTerm[CurrentSystem].yzyz*=r;
  BornTerm[CurrentSystem].zzyz*=r;
  BornTerm[CurrentSystem].xzzz*=r;
  BornTerm[CurrentSystem].yzzz*=r;
  BornTerm[CurrentSystem].zzzz*=r;
}

void AddContributionToCrossTerm(int i,REAL_MATRIX CrossTerm,REAL DDF,REAL DF,VECTOR dr)
{
  // x
  CrossTerm.element[i][0]+=DDF*dr.x*dr.x*dr.x+DF*(dr.x+dr.x);
  CrossTerm.element[i+1][0]+=DDF*dr.x*dr.y*dr.x+DF*dr.y;
  CrossTerm.element[i+2][0]+=DDF*dr.x*dr.z*dr.x+DF*dr.z;

  CrossTerm.element[i][1]+=DDF*dr.y*dr.x*dr.x+DF*dr.y;
  CrossTerm.element[i+1][1]+=DDF*dr.y*dr.y*dr.x;
  CrossTerm.element[i+2][1]+=DDF*dr.y*dr.z*dr.x;

  CrossTerm.element[i][2]+=DDF*dr.z*dr.x*dr.x+DF*dr.z;
  CrossTerm.element[i+1][2]+=DDF*dr.z*dr.y*dr.x;
  CrossTerm.element[i+2][2]+=DDF*dr.z*dr.z*dr.x;

  // y
  CrossTerm.element[i][3]+=DDF*dr.x*dr.x*dr.y;
  CrossTerm.element[i+1][3]+=DDF*dr.x*dr.y*dr.y+DF*dr.x;
  CrossTerm.element[i+2][3]+=DDF*dr.x*dr.z*dr.y;

  CrossTerm.element[i][4]+=DDF*dr.y*dr.x*dr.y+DF*dr.x;
  CrossTerm.element[i+1][4]+=DDF*dr.y*dr.y*dr.y+DF*(dr.y+dr.y);
  CrossTerm.element[i+2][4]+=DDF*dr.y*dr.z*dr.y+DF*dr.z;

  CrossTerm.element[i][5]+=DDF*dr.z*dr.x*dr.y;
  CrossTerm.element[i+1][5]+=DDF*dr.z*dr.y*dr.y+DF*dr.z;
  CrossTerm.element[i+2][5]+=DDF*dr.z*dr.z*dr.y;

  // z
  CrossTerm.element[i][6]+=DDF*dr.x*dr.x*dr.z;
  CrossTerm.element[i+1][6]+=DDF*dr.x*dr.y*dr.z;
  CrossTerm.element[i+2][6]+=DDF*dr.x*dr.z*dr.z+DF*dr.x;

  CrossTerm.element[i][7]+=DDF*dr.y*dr.x*dr.z;
  CrossTerm.element[i+1][7]+=DDF*dr.y*dr.y*dr.z;
  CrossTerm.element[i+2][7]+=DDF*dr.y*dr.z*dr.z+DF*dr.y;

  CrossTerm.element[i][8]+=DDF*dr.z*dr.x*dr.z+DF*dr.x;
  CrossTerm.element[i+1][8]+=DDF*dr.z*dr.y*dr.z+DF*dr.y;
  CrossTerm.element[i+2][8]+=DDF*dr.z*dr.z*dr.z+DF*(dr.z+dr.z);
}

void AddContributionToBornTerm(REAL DDF,REAL DF,VECTOR dr)
{
  BornTerm[CurrentSystem].xxxx+=DDF*dr.x*dr.x*dr.x*dr.x+2.0*DF*dr.x*dr.x;
  BornTerm[CurrentSystem].xxyy+=DDF*dr.x*dr.x*dr.y*dr.y;
  BornTerm[CurrentSystem].xxzz+=DDF*dr.x*dr.x*dr.z*dr.z;
  BornTerm[CurrentSystem].xxyz+=DDF*dr.x*dr.x*dr.y*dr.z;
  BornTerm[CurrentSystem].xxzx+=DDF*dr.x*dr.x*dr.z*dr.x+DF*dr.x*dr.z;
  BornTerm[CurrentSystem].xxxy+=DDF*dr.x*dr.x*dr.x*dr.y+DF*dr.x*dr.y;

  BornTerm[CurrentSystem].yyyy+=DDF*dr.y*dr.y*dr.y*dr.y+2.0*DF*dr.y*dr.y;;
  BornTerm[CurrentSystem].yyzz+=DDF*dr.y*dr.y*dr.z*dr.z;
  BornTerm[CurrentSystem].yyyz+=DDF*dr.y*dr.y*dr.y*dr.z+DF*dr.y*dr.z;
  BornTerm[CurrentSystem].yyzx+=DDF*dr.y*dr.y*dr.z*dr.x;
  BornTerm[CurrentSystem].yyxy+=DDF*dr.y*dr.y*dr.x*dr.y+DF*dr.y*dr.x;

  BornTerm[CurrentSystem].zzzz+=DDF*dr.z*dr.z*dr.z*dr.z+2.0*DF*dr.z*dr.z;
  BornTerm[CurrentSystem].zzyz+=DDF*dr.z*dr.z*dr.y*dr.z+DF*dr.z*dr.y;
  BornTerm[CurrentSystem].zzzx+=DDF*dr.z*dr.z*dr.z*dr.x+DF*dr.z*dr.x;
  BornTerm[CurrentSystem].zzxy+=DDF*dr.z*dr.z*dr.x*dr.y;

  BornTerm[CurrentSystem].yzyz+=DDF*dr.y*dr.z*dr.y*dr.z+0.5*DF*(dr.y*dr.y+dr.z*dr.z);
  BornTerm[CurrentSystem].yzzx+=DDF*dr.y*dr.z*dr.z*dr.x+0.5*DF*(dr.y*dr.x);
  BornTerm[CurrentSystem].yzxy+=DDF*dr.y*dr.z*dr.x*dr.y+0.5*DF*dr.z*dr.x;

  BornTerm[CurrentSystem].zxzx+=DDF*dr.z*dr.x*dr.z*dr.x+0.5*DF*(dr.z*dr.z+dr.x*dr.x);
  BornTerm[CurrentSystem].zxxy+=DDF*dr.z*dr.x*dr.x*dr.y+0.5*DF*dr.z*dr.y;

  BornTerm[CurrentSystem].xyxy+=DDF*dr.x*dr.y*dr.x*dr.y+0.5*DF*(dr.x*dr.x+dr.y*dr.y);
}

void FillInRemainingTermsInBornTerm(void)
{
  // use symmetry Cijkl = Cjikl = Cijlk = Cjilk = Cklij
  BornTerm[CurrentSystem].yxxx=BornTerm[CurrentSystem].xxxy;
  BornTerm[CurrentSystem].xxyx=BornTerm[CurrentSystem].xxxy;
  BornTerm[CurrentSystem].yxyx=BornTerm[CurrentSystem].xyxy;
  BornTerm[CurrentSystem].zxyx=BornTerm[CurrentSystem].zxxy;
  BornTerm[CurrentSystem].yxzx=BornTerm[CurrentSystem].zxxy;
  BornTerm[CurrentSystem].zyxx=BornTerm[CurrentSystem].xxyz;
  BornTerm[CurrentSystem].xyyx=BornTerm[CurrentSystem].xyxy;
  BornTerm[CurrentSystem].yyyx=BornTerm[CurrentSystem].yyxy;
  BornTerm[CurrentSystem].zyyx=BornTerm[CurrentSystem].yzxy;
  BornTerm[CurrentSystem].zyzx=BornTerm[CurrentSystem].yzzx;
  BornTerm[CurrentSystem].xzxx=BornTerm[CurrentSystem].xxzx;
  BornTerm[CurrentSystem].xzyx=BornTerm[CurrentSystem].zxxy;
  BornTerm[CurrentSystem].yzyx=BornTerm[CurrentSystem].yzxy;
  BornTerm[CurrentSystem].zzyx=BornTerm[CurrentSystem].zzxy;
  BornTerm[CurrentSystem].xzzx=BornTerm[CurrentSystem].zxzx;
  BornTerm[CurrentSystem].yxxy=BornTerm[CurrentSystem].xyxy;
  BornTerm[CurrentSystem].yxyy=BornTerm[CurrentSystem].yyxy;
  BornTerm[CurrentSystem].xxzy=BornTerm[CurrentSystem].xxyz;
  BornTerm[CurrentSystem].yxzy=BornTerm[CurrentSystem].yzxy;
  BornTerm[CurrentSystem].zxzy=BornTerm[CurrentSystem].yzzx;
  BornTerm[CurrentSystem].zyxy=BornTerm[CurrentSystem].yzxy;
  BornTerm[CurrentSystem].zyyy=BornTerm[CurrentSystem].yyyz;
  BornTerm[CurrentSystem].xyzy=BornTerm[CurrentSystem].yzxy;
  BornTerm[CurrentSystem].yyzy=BornTerm[CurrentSystem].yyyz;
  BornTerm[CurrentSystem].zyzy=BornTerm[CurrentSystem].yzyz;
  BornTerm[CurrentSystem].xzxy=BornTerm[CurrentSystem].zxxy;
  BornTerm[CurrentSystem].xzyy=BornTerm[CurrentSystem].yyzx;
  BornTerm[CurrentSystem].xzzy=BornTerm[CurrentSystem].yzzx;
  BornTerm[CurrentSystem].yzzy=BornTerm[CurrentSystem].yzyz;
  BornTerm[CurrentSystem].zzzy=BornTerm[CurrentSystem].zzyz;
  BornTerm[CurrentSystem].xxxz=BornTerm[CurrentSystem].xxzx;
  BornTerm[CurrentSystem].yxxz=BornTerm[CurrentSystem].zxxy;
  BornTerm[CurrentSystem].zxxz=BornTerm[CurrentSystem].zxzx;
  BornTerm[CurrentSystem].yxyz=BornTerm[CurrentSystem].yzxy;
  BornTerm[CurrentSystem].yxzz=BornTerm[CurrentSystem].zzxy;
  BornTerm[CurrentSystem].xyxz=BornTerm[CurrentSystem].zxxy;
  BornTerm[CurrentSystem].yyxz=BornTerm[CurrentSystem].yyzx;
  BornTerm[CurrentSystem].zyxz=BornTerm[CurrentSystem].yzzx;
  BornTerm[CurrentSystem].zyyz=BornTerm[CurrentSystem].yzyz;
  BornTerm[CurrentSystem].zyzz=BornTerm[CurrentSystem].zzyz;
  BornTerm[CurrentSystem].xzxz=BornTerm[CurrentSystem].zxzx;
  BornTerm[CurrentSystem].yzxz=BornTerm[CurrentSystem].yzzx;
  BornTerm[CurrentSystem].zzxz=BornTerm[CurrentSystem].zzzx;
  BornTerm[CurrentSystem].xzyz=BornTerm[CurrentSystem].yzzx;
  BornTerm[CurrentSystem].xzzz=BornTerm[CurrentSystem].zzzx;

  // use symmetry Cijkl = Cklij
  BornTerm[CurrentSystem].zxxx=BornTerm[CurrentSystem].xxzx;
  BornTerm[CurrentSystem].xyxx=BornTerm[CurrentSystem].xxxy;
  BornTerm[CurrentSystem].yyxx=BornTerm[CurrentSystem].xxyy;
  BornTerm[CurrentSystem].xyzx=BornTerm[CurrentSystem].zxxy;
  BornTerm[CurrentSystem].yzxx=BornTerm[CurrentSystem].xxyz;
  BornTerm[CurrentSystem].zzxx=BornTerm[CurrentSystem].xxzz;
  BornTerm[CurrentSystem].zxyy=BornTerm[CurrentSystem].yyzx;
  BornTerm[CurrentSystem].xyyy=BornTerm[CurrentSystem].yyxy;
  BornTerm[CurrentSystem].yzyy=BornTerm[CurrentSystem].yyyz;
  BornTerm[CurrentSystem].zzyy=BornTerm[CurrentSystem].yyzz;
  BornTerm[CurrentSystem].zxyz=BornTerm[CurrentSystem].yzzx;
  BornTerm[CurrentSystem].zxzz=BornTerm[CurrentSystem].zzzx;
  BornTerm[CurrentSystem].xyyz=BornTerm[CurrentSystem].yzxy;
  BornTerm[CurrentSystem].xyzz=BornTerm[CurrentSystem].zzxy;
  BornTerm[CurrentSystem].yzzz=BornTerm[CurrentSystem].zzyz;
}

int MapCellListId(int i,int j,int k)
{
  int imap;

  imap=((i+NumberOfCellListCells[CurrentSystem].x)%NumberOfCellListCells[CurrentSystem].x)+
       ((j+NumberOfCellListCells[CurrentSystem].y)%NumberOfCellListCells[CurrentSystem].y)*NumberOfCellListCells[CurrentSystem].x+
       ((k+NumberOfCellListCells[CurrentSystem].z)%NumberOfCellListCells[CurrentSystem].z)*
          NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;
  return imap;
}

void InverseBoxMatrix(REAL_MATRIX3x3 *Box,REAL_MATRIX3x3 *InverseBox)
{
  REAL det;

  switch(Dimension)
  {
    case 2:
      det=1.0/(Box->ax*Box->by-Box->ay*Box->bx);
      InverseBox->ax=det*Box->by;  InverseBox->bx=-det*Box->bx; InverseBox->cx=0.0;
      InverseBox->ay=-det*Box->ay; InverseBox->by=det*Box->ax;  InverseBox->cy=0.0;
      InverseBox->az=0.0;          InverseBox->cy=0.0;          InverseBox->cz=0.0;
      break;
    case 3:
      Invert3x3Matrix(Box,InverseBox,&det);
      break;
  }
}


void InitializeReplicaBox(void)
{
  int i,j,k;
  int k1,k2,k3,ncell,imap;

  if(NumberOfReplicaCells[CurrentSystem].x<NumberOfUnitCells[CurrentSystem].x) NumberOfReplicaCells[CurrentSystem].x=1;
  if(NumberOfReplicaCells[CurrentSystem].y<NumberOfUnitCells[CurrentSystem].y) NumberOfReplicaCells[CurrentSystem].y=1;
  if(NumberOfReplicaCells[CurrentSystem].z<NumberOfUnitCells[CurrentSystem].z) NumberOfReplicaCells[CurrentSystem].z=1;

  TotalNumberOfReplicaCells[CurrentSystem]=NumberOfReplicaCells[CurrentSystem].x*NumberOfReplicaCells[CurrentSystem].y*
                             NumberOfReplicaCells[CurrentSystem].z;

  ReplicaShift=(VECTOR*)calloc(TotalNumberOfReplicaCells[CurrentSystem],sizeof(VECTOR));

  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }

  ReplicaBox[CurrentSystem].ax=Box[CurrentSystem].ax*NumberOfReplicaCells[CurrentSystem].x;
  ReplicaBox[CurrentSystem].ay=Box[CurrentSystem].ay*NumberOfReplicaCells[CurrentSystem].x;
  ReplicaBox[CurrentSystem].az=Box[CurrentSystem].az*NumberOfReplicaCells[CurrentSystem].x;
  ReplicaBox[CurrentSystem].bx=Box[CurrentSystem].bx*NumberOfReplicaCells[CurrentSystem].y;
  ReplicaBox[CurrentSystem].by=Box[CurrentSystem].by*NumberOfReplicaCells[CurrentSystem].y;
  ReplicaBox[CurrentSystem].bz=Box[CurrentSystem].bz*NumberOfReplicaCells[CurrentSystem].y;
  ReplicaBox[CurrentSystem].cx=Box[CurrentSystem].cx*NumberOfReplicaCells[CurrentSystem].z;
  ReplicaBox[CurrentSystem].cy=Box[CurrentSystem].cy*NumberOfReplicaCells[CurrentSystem].z;
  ReplicaBox[CurrentSystem].cz=Box[CurrentSystem].cz*NumberOfReplicaCells[CurrentSystem].z;

  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);

  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if((NumberOfReplicaCells[CurrentSystem].x>NumberOfUnitCells[CurrentSystem].x)||
       (NumberOfReplicaCells[CurrentSystem].y>NumberOfUnitCells[CurrentSystem].y)||
       (NumberOfReplicaCells[CurrentSystem].z>NumberOfUnitCells[CurrentSystem].z))
      UseReplicas[CurrentSystem]=TRUE;
    else
      UseReplicas[CurrentSystem]=FALSE;
  }
  else
  {
    if(TotalNumberOfReplicaCells[CurrentSystem]>1) UseReplicas[CurrentSystem]=TRUE;
    else UseReplicas[CurrentSystem]=FALSE;
  }

  NumberOfCellListCells[CurrentSystem].x=(int)(BoxProperties[CurrentSystem].cx/CutOffVDW);
  NumberOfCellListCells[CurrentSystem].y=(int)(BoxProperties[CurrentSystem].cy/CutOffVDW);
  NumberOfCellListCells[CurrentSystem].z=(int)(BoxProperties[CurrentSystem].cz/CutOffVDW);
  NumberOfCellLists[CurrentSystem]=NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y*
                                   NumberOfCellListCells[CurrentSystem].z;

  // use cell lists for sufficiently large systems
  if((NumberOfCellListCells[CurrentSystem].x>3)&&(NumberOfCellListCells[CurrentSystem].y>3)&&
     (NumberOfCellListCells[CurrentSystem].z>3)&&(BoundaryCondition[CurrentSystem]==RECTANGULAR))
//  if((NumberOfCellListCells[CurrentSystem].x>3)&&(NumberOfCellListCells[CurrentSystem].y>3)&&
//     (NumberOfCellListCells[CurrentSystem].z>3))
    UseCellLists[CurrentSystem]=TRUE;
  else
    UseCellLists[CurrentSystem]=FALSE;

  UseCellLists[CurrentSystem]=FALSE;

  if(UseCellLists[CurrentSystem])
  {
    AllocateFrameworkCellList();

    CellListMap[CurrentSystem]=(int (*)[27])calloc(NumberOfCellLists[CurrentSystem],sizeof(int[27]));

    for(i=0;i<NumberOfCellListCells[CurrentSystem].x;i++)
      for(j=0;j<NumberOfCellListCells[CurrentSystem].y;j++)
        for(k=0;k<NumberOfCellListCells[CurrentSystem].z;k++)
        {
          imap=MapCellListId(i,j,k);

          CellListMap[CurrentSystem][imap][0]=MapCellListId(i-1,j-1,k-1);
          CellListMap[CurrentSystem][imap][1]=MapCellListId(i-1,j-1,k);
          CellListMap[CurrentSystem][imap][2]=MapCellListId(i-1,j-1,k+1);
          CellListMap[CurrentSystem][imap][3]=MapCellListId(i-1,j,k-1);
          CellListMap[CurrentSystem][imap][4]=MapCellListId(i-1,j,k);
          CellListMap[CurrentSystem][imap][5]=MapCellListId(i-1,j,k+1);
          CellListMap[CurrentSystem][imap][6]=MapCellListId(i-1,j+1,k-1);
          CellListMap[CurrentSystem][imap][7]=MapCellListId(i-1,j+1,k);
          CellListMap[CurrentSystem][imap][8]=MapCellListId(i-1,j+1,k+1);

          CellListMap[CurrentSystem][imap][9]=MapCellListId(i,j-1,k-1);
          CellListMap[CurrentSystem][imap][10]=MapCellListId(i,j-1,k);
          CellListMap[CurrentSystem][imap][11]=MapCellListId(i,j-1,k+1);
          CellListMap[CurrentSystem][imap][12]=MapCellListId(i,j,k-1);
          CellListMap[CurrentSystem][imap][13]=MapCellListId(i,j,k);
          CellListMap[CurrentSystem][imap][14]=MapCellListId(i,j,k+1);
          CellListMap[CurrentSystem][imap][15]=MapCellListId(i,j+1,k-1);
          CellListMap[CurrentSystem][imap][16]=MapCellListId(i,j+1,k);
          CellListMap[CurrentSystem][imap][17]=MapCellListId(i,j+1,k+1);

          CellListMap[CurrentSystem][imap][18]=MapCellListId(i+1,j-1,k-1);
          CellListMap[CurrentSystem][imap][19]=MapCellListId(i+1,j-1,k);
          CellListMap[CurrentSystem][imap][20]=MapCellListId(i+1,j-1,k+1);
          CellListMap[CurrentSystem][imap][21]=MapCellListId(i+1,j,k-1);
          CellListMap[CurrentSystem][imap][22]=MapCellListId(i+1,j,k);
          CellListMap[CurrentSystem][imap][23]=MapCellListId(i+1,j,k+1);
          CellListMap[CurrentSystem][imap][24]=MapCellListId(i+1,j+1,k-1);
          CellListMap[CurrentSystem][imap][25]=MapCellListId(i+1,j+1,k);
          CellListMap[CurrentSystem][imap][26]=MapCellListId(i+1,j+1,k+1);
        }

        MakeFrameworkCellList();
  }
}

void WriteRestartSimulation(FILE *FilePtr)
{
  int i,j;
  REAL Check;

  fwrite(&seed,sizeof(seed),1,FilePtr);
  fwrite(&ForceField,sizeof(ForceField),1,FilePtr);
  fwrite(&FileNameAppend,sizeof(FileNameAppend),1,FilePtr);

  fwrite(&Dimension,sizeof(Dimension),1,FilePtr);

  fwrite(&SimulationType,sizeof(SimulationType),1,FilePtr);

  fwrite(&NumberOfSystems,sizeof(NumberOfSystems),1,FilePtr);
  fwrite(&CurrentSystem,sizeof(CurrentSystem),1,FilePtr);

  fwrite(&NumberOfPartialPressures,sizeof(NumberOfPartialPressures),1,FilePtr);
  fwrite(&CurrentPartialPressure,sizeof(CurrentPartialPressure),1,FilePtr);

  fwrite(&NumberOfCycles,sizeof(NumberOfCycles),1,FilePtr);
  fwrite(&NumberOfVelocityScalingCycles,sizeof(NumberOfVelocityScalingCycles),1,FilePtr);
  fwrite(&NumberOfEquilibrationCycles,sizeof(NumberOfEquilibrationCycles),1,FilePtr);
  fwrite(&NumberOfInitializationCycles,sizeof(NumberOfInitializationCycles),1,FilePtr);
  fwrite(&CurrentCycle,sizeof(CurrentCycle),1,FilePtr);
  fwrite(&SimulationStage,sizeof(SimulationStage),1,FilePtr);
  fwrite(&UseReducedUnits,sizeof(UseReducedUnits),1,FilePtr);
  fwrite(&PrintEvery,sizeof(PrintEvery),1,FilePtr);
  fwrite(&OptimizeAcceptenceEvery,sizeof(OptimizeAcceptenceEvery),1,FilePtr);
  fwrite(&PrintPropertiesEvery,sizeof(PrintPropertiesEvery),1,FilePtr);
  fwrite(&Output,sizeof(Output),1,FilePtr);
  fwrite(&WriteForcefieldToOutput,sizeof(WriteForcefieldToOutput),1,FilePtr);
  fwrite(&WritePseudoAtomsToOutput,sizeof(WritePseudoAtomsToOutput),1,FilePtr);
  fwrite(&WriteMoleculeDefinitionToOutput,sizeof(WriteMoleculeDefinitionToOutput),1,FilePtr);
  fwrite(&InputFileType,sizeof(InputFileType),1,FilePtr);
  fwrite(&WriteBinaryRestartFileEvery,sizeof(WriteBinaryRestartFileEvery),1,FilePtr);

  fwrite(&OmitInterMolecularInteractions,sizeof(int),1,FilePtr);
  fwrite(&OmitInterMolecularVDWInteractions,sizeof(int),1,FilePtr);
  fwrite(&OmitInterMolecularCoulombInteractions,sizeof(int),1,FilePtr);
  fwrite(&OmitAdsorbateAdsorbateVDWInteractions,sizeof(int),1,FilePtr);
  fwrite(&OmitAdsorbateAdsorbateCoulombInteractions,sizeof(int),1,FilePtr);
  fwrite(&OmitAdsorbateCationVDWInteractions,sizeof(int),1,FilePtr);
  fwrite(&OmitAdsorbateCationCoulombInteractions,sizeof(int),1,FilePtr);
  fwrite(&OmitCationCationVDWInteractions,sizeof(int),1,FilePtr);
  fwrite(&OmitCationCationCoulombInteractions,sizeof(int),1,FilePtr);

  fwrite(&OmitAdsorbateAdsorbatePolarization,sizeof(int),1,FilePtr);
  fwrite(&OmitCationCationPolarization,sizeof(int),1,FilePtr);
  fwrite(&OmitIntraFrameworkPolarization,sizeof(int),1,FilePtr);

  fwrite(&ComputeStress,sizeof(int),1,FilePtr);
  fwrite(&ComputeHessian,sizeof(int),1,FilePtr);
  fwrite(&ComputeBornTerm,sizeof(int),1,FilePtr);
  fwrite(&ComputeCrossTerm,sizeof(int),1,FilePtr);

  fwrite(&ComputePolarization,sizeof(int),1,FilePtr);
  fwrite(&BackPolarization,sizeof(int),1,FilePtr);
  fwrite(&NumberOfBackPolarizationSteps,sizeof(int),1,FilePtr);

  fwrite(&ReinitializeVelocities,sizeof(ReinitializeVelocities),1,FilePtr);
  fwrite(&QuenchCoreShellVelocities,sizeof(QuenchCoreShellVelocities),1,FilePtr);
  fwrite(&QuenchCoreShellVelocitiesEvery,sizeof(QuenchCoreShellVelocitiesEvery),1,FilePtr);

  fwrite(&MinimumInnerCycles,sizeof(MinimumInnerCycles),1,FilePtr);

  fwrite(&ChargeMethod,sizeof(int),1,FilePtr);
  fwrite(&ChargeFromChargeEquilibration,sizeof(int),1,FilePtr);
  fwrite(&ChargeEquilibrationMethod,sizeof(int),1,FilePtr);
  fwrite(&SymmetrizeFrameworkCharges,sizeof(int),1,FilePtr);
  fwrite(&UseChargesFromMOLFile,sizeof(int),1,FilePtr);
  fwrite(&UseChargesFromCIFFile,sizeof(int),1,FilePtr);

  fwrite(&OverlapDistanceSquared,sizeof(OverlapDistanceSquared),1,FilePtr);
  fwrite(&CutOffVDW,sizeof(CutOffVDW),1,FilePtr);
  fwrite(&CutOffVDWSquared,sizeof(CutOffVDWSquared),1,FilePtr);
  fwrite(&CutOffVDWSwitch,sizeof(CutOffVDWSwitch),1,FilePtr);
  fwrite(&CutOffVDWSwitchSquared,sizeof(CutOffVDWSwitchSquared),1,FilePtr);
  fwrite(&InverseCutOffVDW,sizeof(InverseCutOffVDW),1,FilePtr);

  fwrite(CutOffChargeCharge,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CutOffChargeChargeSquared,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CutOffChargeChargeSwitch,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CutOffChargeChargeSwitchSquared,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(InverseCutOffChargeCharge,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(&CutOffChargeBondDipole,sizeof(CutOffChargeBondDipole),1,FilePtr);
  fwrite(&CutOffChargeBondDipoleSquared,sizeof(CutOffChargeBondDipoleSquared),1,FilePtr);
  fwrite(&CutOffChargeBondDipoleSwitch,sizeof(CutOffChargeBondDipoleSwitch),1,FilePtr);
  fwrite(&CutOffChargeBondDipoleSwitchSquared,sizeof(CutOffChargeBondDipoleSwitchSquared),1,FilePtr);

  fwrite(&CutOffBondDipoleBondDipole,sizeof(CutOffBondDipoleBondDipole),1,FilePtr);
  fwrite(&CutOffBondDipoleBondDipoleSquared,sizeof(CutOffBondDipoleBondDipoleSquared),1,FilePtr);
  fwrite(&CutOffBondDipoleBondDipoleSwitch,sizeof(CutOffBondDipoleBondDipoleSwitch),1,FilePtr);
  fwrite(&CutOffBondDipoleBondDipoleSwitchSquared,sizeof(CutOffBondDipoleBondDipoleSwitchSquared),1,FilePtr);

  fwrite(&DeltaT,sizeof(DeltaT),1,FilePtr);
  fwrite(Box,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(InverseBox,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(ReplicaBox,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(InverseReplicaBox,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(NumberOfReplicaCells,NumberOfSystems,sizeof(INT_VECTOR3),FilePtr);
  fwrite(TotalNumberOfReplicaCells,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(UseReplicas,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(BoxProperties,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(InverseBoxProperties,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);

  fwrite(Volume,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(AlphaAngle,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(BetaAngle,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(GammaAngle,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(BoundaryCondition,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(IntialCenterOfMassPosition,sizeof(VECTOR),NumberOfSystems,FilePtr);
  fwrite(Beta,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(StrainDerivativeTailCorrection,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(StrainDerivativeRigidCorrection,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);

  fwrite(UseCellLists,sizeof(int),NumberOfSystems,FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(UseCellLists[i])
    {
      fwrite(&NumberOfCellListCells[i].x,sizeof(int),1,FilePtr);
      fwrite(&NumberOfCellListCells[i].y,sizeof(int),1,FilePtr);
      fwrite(&NumberOfCellListCells[i].z,sizeof(int),1,FilePtr);
      fwrite(&NumberOfCellLists[i],sizeof(int),1,FilePtr);

      fwrite(CellListMap[i],sizeof(int[27]),NumberOfCellLists[i],FilePtr);
    }
  }

  fwrite(UTotal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(ConservedEnergy,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(Drift,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(ReferenceEnergy,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UIon,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UDrift,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UTailCorrection,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UDistanceConstraints,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAngleConstraints,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UDihedralConstraints,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UInversionBendConstraints,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UOutOfPlaneDistanceConstraints,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UExclusionConstraints,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostPolarization,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbatePolarization,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationPolarization,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostBackPolarization,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateBackPolarization,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationBackPolarization,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostHost,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateCation,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostHostCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostHostVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostAdsorbateCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostAdsorbateVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostCationCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostCationVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateAdsorbateCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateAdsorbateVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateCationCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateCationVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationCationCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationCationVDW,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostHostChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostAdsorbateChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostCationChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateAdsorbateChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateCationChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationCationChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostHostChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostAdsorbateChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostCationChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateAdsorbateChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateCationChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationCationChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostHostBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostAdsorbateBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostCationBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateAdsorbateBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateCationBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationCationBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostHostChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostAdsorbateChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostCationChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateAdsorbateChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateCationChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationCationChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostHostChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostAdsorbateChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostCationChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateAdsorbateChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateCationChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationCationChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostHostBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostAdsorbateBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostCationBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateAdsorbateBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateCationBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationCationBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UHostBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostUreyBradley,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostInversionBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostImproperTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostOutOfPlane,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostBondBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostBondBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostBendBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostBondTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostBendTorsion,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UCationBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationUreyBradley,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationInversionBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationImproperTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationOutOfPlane,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationBondBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationBondBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationBendBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationBondTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationBendTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationIntraVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationIntraChargeCharge,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationIntraChargeBondDipole,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationIntraBondDipoleBondDipole,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UAdsorbateBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateUreyBradley,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateInversionBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateImproperTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateOutOfPlane,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateBondBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateBondBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateBendBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateBondTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateBendTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateIntraVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateIntraChargeCharge,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateIntraChargeBondDipole,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateIntraBondDipoleBondDipole,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(Ensemble,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(InitEnsemble,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(RunEnsemble,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(NPTPRCellType,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(MonoclinicAngleType,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(MCEnsemble,sizeof(int),NumberOfSystems,FilePtr);

  fwrite(DegreesOfFreedom,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomTranslation,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomRotation,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomVibration,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomConstraint,sizeof(int),NumberOfSystems,FilePtr);

  fwrite(DegreesOfFreedomFramework,sizeof(int),NumberOfSystems,FilePtr);

  fwrite(DegreesOfFreedomAdsorbates,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomTranslationalAdsorbates,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomRotationalAdsorbates,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomVibrationalAdsorbates,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomConstraintAdsorbates,sizeof(int),NumberOfSystems,FilePtr);

  fwrite(DegreesOfFreedomCations,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomTranslationalCations,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomRotationalCations,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomVibrationalCations,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(DegreesOfFreedomConstraintCations,sizeof(int),NumberOfSystems,FilePtr);

  fwrite(UKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UHostKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationKinetic,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UAdsorbateRotationalKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UAdsorbateTranslationalKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationRotationalKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UCationTranslationalKinetic,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(UNoseHoover,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UNoseHooverAdsorbates,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UNoseHooverCations,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(UNoseHooverFramework,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(PressureTrace,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(Pressure,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(StressTensor,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(ConfigurationalStressTensor,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(KineticStressTensor,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(StrainDerivativeTensor,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(MolecularStressTensor,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);

  fwrite(BornTerm,sizeof(REAL_MATRIX9x9),NumberOfSystems,FilePtr);
  fwrite(RelaxationTerm,sizeof(REAL_MATRIX9x9),NumberOfSystems,FilePtr);

  fwrite(&ProbabilityVolumeChangeMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityBoxShapeChangeMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityGibbsVolumeChangeMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityParallelTemperingMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityHyperParallelTemperingMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityParallelMolFractionMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityChiralInversionMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityHybridNVEMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityHybridNPHMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityHybridNPHPRMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityFrameworkChangeMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityFrameworkShiftMove,sizeof(REAL),1,FilePtr);
  fwrite(&ProbabilityCFCRXMCLambdaChangeMove,sizeof(REAL),1,FilePtr);


  fwrite(&CpuTimeProductionRun,sizeof(REAL),1,FilePtr);
  fwrite(&CpuTimeInitialization,sizeof(REAL),1,FilePtr);
  fwrite(&CpuTimeEquilibration,sizeof(REAL),1,FilePtr);
  fwrite(&CpuTotal,sizeof(REAL),1,FilePtr);

  fwrite(CpuTimeParallelTemperingMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeHyperParallelTemperingMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeParallelMolFractionMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeChiralInversionMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeHybridNVEMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeHybridNPHMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeHybridNPHPRMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeVolumeChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeBoxShapeChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeGibbsVolumeChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeFrameworkChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuTimeFrameworkShiftMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(CpuCFCRXMCLambdaChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(&OptimizeVolumeChange,sizeof(int),1,FilePtr);
  fwrite(&OptimizeGibbsVolumeChange,sizeof(int),1,FilePtr);
  fwrite(&OptimizeTranslation,sizeof(int),1,FilePtr);
  fwrite(&OptimizeRotation,sizeof(int),1,FilePtr);
  fwrite(&OptimizeFrameworkChange,sizeof(int),1,FilePtr);
  fwrite(&OptimizeFrameworkShift,sizeof(int),1,FilePtr);
  fwrite(&OptimizeCFLambdaChange,sizeof(int),1,FilePtr);
  fwrite(&OptimizeCFGibbsLambdaChange,sizeof(int),1,FilePtr);
  fwrite(&OptimizeCBCFLambdaChange,sizeof(int),1,FilePtr);
  fwrite(&OptimizeCBCFGibbsLambdaChange,sizeof(int),1,FilePtr);
  fwrite(&OptimizeRXMCLambdaChange,sizeof(int),1,FilePtr);

  // CFRXMC data
  fwrite(&NumberOfReactions,sizeof(int),1,FilePtr);
  if(NumberOfReactions>0)
  {
    fwrite(&RXMCLambdaHistogramSize,sizeof(int),1,FilePtr);
    for(i=0;i<NumberOfSystems;i++)
    {
      fwrite(CFCRXMCLambda[i],sizeof(REAL),NumberOfReactions,FilePtr);
      fwrite(CFRXMCWangLandauScalingFactor[i],sizeof(REAL),NumberOfReactions,FilePtr);
      for(j=0;j<NumberOfReactions;j++)
        fwrite(RXMCBiasingFactors[i][j],sizeof(REAL),RXMCLambdaHistogramSize,FilePtr);
    }
  }


  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void AllocateSimulationMemory(void)
{
  int i,j;

  CpuTimeProductionRun=0.0;
  CpuTimeInitialization=0.0;
  CpuTimeEquilibration=0.0;
  CpuTotal=0.0;

  CpuTimeParallelTemperingMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeHyperParallelTemperingMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeParallelMolFractionMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeChiralInversionMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeHybridNVEMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeHybridNPHMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeHybridNPHPRMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeVolumeChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeBoxShapeChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeGibbsVolumeChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeFrameworkChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuTimeFrameworkShiftMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CpuCFCRXMCLambdaChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  Box=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  InverseBox=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  ReplicaBox=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  InverseReplicaBox=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  NumberOfReplicaCells=(INT_VECTOR3*)calloc(NumberOfSystems,sizeof(INT_VECTOR3));
  TotalNumberOfReplicaCells=(int*)calloc(NumberOfSystems,sizeof(int));
  UseReplicas=(int*)calloc(NumberOfSystems,sizeof(int));
  BoxProperties=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  InverseBoxProperties=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  Volume=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  AlphaAngle=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  BetaAngle=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  GammaAngle=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  BoundaryCondition=(int*)calloc(NumberOfSystems,sizeof(int));
  IntialCenterOfMassPosition=(VECTOR*)calloc(NumberOfSystems,sizeof(VECTOR));
  Beta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  StrainDerivativeTailCorrection=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  StrainDerivativeRigidCorrection=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));

  CutOffChargeCharge=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CutOffChargeChargeSquared=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CutOffChargeChargeSwitch=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  CutOffChargeChargeSwitchSquared=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  InverseCutOffChargeCharge=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // CFRXMC
  //----------------------------------------------------------------------------------------

/*
  CFCRXMCLambda=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  for(i=0;i<NumberOfSystems;i++)
     CFCRXMCLambda[i]=(REAL*)calloc(NumberOfReactions,sizeof(REAL));
  ReactantsStoichiometry=(int**)calloc(NumberOfReactions,sizeof(int*));
  ProductsStoichiometry=(int**)calloc(NumberOfReactions,sizeof(int*));
  for(i=0;i<NumberOfReactions;i++)
  {
    ReactantsStoichiometry[i]=(int*)calloc(NumberOfComponents,sizeof(int));
    ProductsStoichiometry[i]=(int*)calloc(NumberOfComponents,sizeof(int));
  }
*/

  CFCRXMCLambda=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  RXMCBiasingFactors=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  CFRXMCWangLandauScalingFactor=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  for(i=0;i<NumberOfSystems;i++)
  {
    if(NumberOfReactions>0)
    {
      CFCRXMCLambda[i]=(REAL*)calloc(NumberOfReactions,sizeof(REAL));
      RXMCBiasingFactors[i]=(REAL**)calloc(NumberOfReactions,sizeof(REAL*));
      CFRXMCWangLandauScalingFactor[i]=(REAL*)calloc(NumberOfReactions,sizeof(REAL));
      for(j=0;j<NumberOfReactions;j++)
        RXMCBiasingFactors[i][j]=(REAL*)calloc(RXMCLambdaHistogramSize,sizeof(REAL));
    }
  }


  NumberOfCellListCells=(INT_VECTOR3*)calloc(NumberOfSystems,sizeof(INT_VECTOR3));
  UseCellLists=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfCellLists=(int*)calloc(NumberOfSystems,sizeof(int));
  CellListMap=(int (**)[27])calloc(NumberOfSystems,sizeof(int(*)[27]));

  UTotal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  ConservedEnergy=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Drift=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  ReferenceEnergy=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UIon=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UDrift=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UTailCorrection=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UDistanceConstraints=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAngleConstraints=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UDihedralConstraints=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UInversionBendConstraints=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UOutOfPlaneDistanceConstraints=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UExclusionConstraints=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostPolarization=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbatePolarization=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationPolarization=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostBackPolarization=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBackPolarization=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBackPolarization=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHost=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbate=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbate=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHostCoulomb=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostHostVDW=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateCoulomb=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateVDW=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationCoulomb=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationVDW=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateCoulomb=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateVDW=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationCoulomb=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationVDW=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationCoulomb=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationVDW=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHostChargeChargeReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateChargeChargeReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationChargeChargeReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateChargeChargeReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationChargeChargeReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationChargeChargeReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHostChargeBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateChargeBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationChargeBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateChargeBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationChargeBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationChargeBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHostBondDipoleBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateBondDipoleBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationBondDipoleBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateBondDipoleBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationBondDipoleBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationBondDipoleBondDipoleReal=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHostChargeChargeFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateChargeChargeFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationChargeChargeFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateChargeChargeFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationChargeChargeFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationChargeChargeFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHostChargeBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateChargeBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationChargeBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateChargeBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationChargeBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationChargeBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHostBondDipoleBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateBondDipoleBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationBondDipoleBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateBondDipoleBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationBondDipoleBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationBondDipoleBondDipoleFourier=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostBond=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostUreyBradley=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostInversionBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostImproperTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostOutOfPlane=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBondBond=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBondBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBendBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBondTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBendTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UCationBond=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationUreyBradley=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationInversionBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationImproperTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationOutOfPlane=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBondBond=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBondBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBendBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBondTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBendTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationIntraVDW=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationIntraChargeCharge=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationIntraChargeBondDipole=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationIntraBondDipoleBondDipole=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UAdsorbateBond=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateUreyBradley=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateInversionBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateImproperTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateOutOfPlane=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBondBond=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBondBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBendBend=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBondTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBendTorsion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateIntraVDW=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateIntraChargeCharge=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateIntraChargeBondDipole=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateIntraBondDipoleBondDipole=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  Ensemble=(int*)calloc(NumberOfSystems,sizeof(int));
  InitEnsemble=(int*)calloc(NumberOfSystems,sizeof(int));
  RunEnsemble=(int*)calloc(NumberOfSystems,sizeof(int));
  NPTPRCellType=(int*)calloc(NumberOfSystems,sizeof(int));
  MonoclinicAngleType=(int*)calloc(NumberOfSystems,sizeof(int));
  MCEnsemble=(int*)calloc(NumberOfSystems,sizeof(int));

  DegreesOfFreedom=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomTranslation=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomRotation=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomVibration=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomConstraint=(int*)calloc(NumberOfSystems,sizeof(int));

  DegreesOfFreedomFramework=(int*)calloc(NumberOfSystems,sizeof(int));

  DegreesOfFreedomAdsorbates=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomTranslationalAdsorbates=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomRotationalAdsorbates=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomVibrationalAdsorbates=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomConstraintAdsorbates=(int*)calloc(NumberOfSystems,sizeof(int));

  DegreesOfFreedomCations=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomTranslationalCations=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomRotationalCations=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomVibrationalCations=(int*)calloc(NumberOfSystems,sizeof(int));
  DegreesOfFreedomConstraintCations=(int*)calloc(NumberOfSystems,sizeof(int));

  UKinetic=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostKinetic=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationKinetic=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateKinetic=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UAdsorbateRotationalKinetic=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateTranslationalKinetic=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationRotationalKinetic=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationTranslationalKinetic=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UNoseHoover=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UNoseHooverAdsorbates=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UNoseHooverCations=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UNoseHooverFramework=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  PressureTrace=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Pressure=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  StressTensor=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  ConfigurationalStressTensor=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  KineticStressTensor=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  StrainDerivativeTensor=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  MolecularStressTensor=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));

  BornTerm=(REAL_MATRIX9x9*)calloc(NumberOfSystems,sizeof(REAL_MATRIX9x9));
  RelaxationTerm=(REAL_MATRIX9x9*)calloc(NumberOfSystems,sizeof(REAL_MATRIX9x9));

  for(i=0;i<NumberOfSystems;i++)
    MonoclinicAngleType[i]=MONOCLINIC_BETA_ANGLE;
}

void ReadRestartSimulation(FILE *FilePtr)
{
  int i,j;
  REAL Check;

  fread(&seed,sizeof(seed),1,FilePtr);
  fread(&ForceField,sizeof(ForceField),1,FilePtr);
  fread(&FileNameAppend,sizeof(FileNameAppend),1,FilePtr);

  fread(&Dimension,sizeof(Dimension),1,FilePtr);

  fread(&SimulationType,sizeof(SimulationType),1,FilePtr);

  fread(&NumberOfSystems,sizeof(NumberOfSystems),1,FilePtr);
  fread(&CurrentSystem,sizeof(CurrentSystem),1,FilePtr);

  AllocateSimulationMemory();
  AllocateMovieMemory();
  AllocateRecrossingMemory();

  fread(&NumberOfPartialPressures,sizeof(NumberOfPartialPressures),1,FilePtr);
  fread(&CurrentPartialPressure,sizeof(CurrentPartialPressure),1,FilePtr);

  fread(&NumberOfCycles,sizeof(NumberOfCycles),1,FilePtr);
  fread(&NumberOfVelocityScalingCycles,sizeof(NumberOfVelocityScalingCycles),1,FilePtr);
  fread(&NumberOfEquilibrationCycles,sizeof(NumberOfEquilibrationCycles),1,FilePtr);
  fread(&NumberOfInitializationCycles,sizeof(NumberOfInitializationCycles),1,FilePtr);
  fread(&CurrentCycle,sizeof(CurrentCycle),1,FilePtr);
  fread(&SimulationStage,sizeof(SimulationStage),1,FilePtr);
  fread(&UseReducedUnits,sizeof(UseReducedUnits),1,FilePtr);
  fread(&PrintEvery,sizeof(PrintEvery),1,FilePtr);
  fread(&OptimizeAcceptenceEvery,sizeof(OptimizeAcceptenceEvery),1,FilePtr);
  fread(&PrintPropertiesEvery,sizeof(PrintPropertiesEvery),1,FilePtr);
  fread(&Output,sizeof(Output),1,FilePtr);
  fread(&WriteForcefieldToOutput,sizeof(WriteForcefieldToOutput),1,FilePtr);
  fread(&WritePseudoAtomsToOutput,sizeof(WritePseudoAtomsToOutput),1,FilePtr);
  fread(&WriteMoleculeDefinitionToOutput,sizeof(WriteMoleculeDefinitionToOutput),1,FilePtr);
  fread(&InputFileType,sizeof(InputFileType),1,FilePtr);
  fread(&WriteBinaryRestartFileEvery,sizeof(WriteBinaryRestartFileEvery),1,FilePtr);

  fread(&OmitInterMolecularInteractions,sizeof(int),1,FilePtr);
  fread(&OmitInterMolecularVDWInteractions,sizeof(int),1,FilePtr);
  fread(&OmitInterMolecularCoulombInteractions,sizeof(int),1,FilePtr);
  fread(&OmitAdsorbateAdsorbateVDWInteractions,sizeof(int),1,FilePtr);
  fread(&OmitAdsorbateAdsorbateCoulombInteractions,sizeof(int),1,FilePtr);
  fread(&OmitAdsorbateCationVDWInteractions,sizeof(int),1,FilePtr);
  fread(&OmitAdsorbateCationCoulombInteractions,sizeof(int),1,FilePtr);
  fread(&OmitCationCationVDWInteractions,sizeof(int),1,FilePtr);
  fread(&OmitCationCationCoulombInteractions,sizeof(int),1,FilePtr);

  fread(&OmitAdsorbateAdsorbatePolarization,sizeof(int),1,FilePtr);
  fread(&OmitCationCationPolarization,sizeof(int),1,FilePtr);
  fread(&OmitIntraFrameworkPolarization,sizeof(int),1,FilePtr);

  fread(&ComputeStress,sizeof(int),1,FilePtr);
  fread(&ComputeHessian,sizeof(int),1,FilePtr);
  fread(&ComputeBornTerm,sizeof(int),1,FilePtr);
  fread(&ComputeCrossTerm,sizeof(int),1,FilePtr);

  fread(&ComputePolarization,sizeof(int),1,FilePtr);
  fread(&BackPolarization,sizeof(int),1,FilePtr);
  fread(&NumberOfBackPolarizationSteps,sizeof(int),1,FilePtr);

  fread(&ReinitializeVelocities,sizeof(ReinitializeVelocities),1,FilePtr);
  fread(&QuenchCoreShellVelocities,sizeof(QuenchCoreShellVelocities),1,FilePtr);
  fread(&QuenchCoreShellVelocitiesEvery,sizeof(QuenchCoreShellVelocitiesEvery),1,FilePtr);

  fread(&MinimumInnerCycles,sizeof(MinimumInnerCycles),1,FilePtr);

  fread(&ChargeMethod,sizeof(int),1,FilePtr);
  fread(&ChargeFromChargeEquilibration,sizeof(int),1,FilePtr);
  fread(&ChargeEquilibrationMethod,sizeof(int),1,FilePtr);
  fread(&SymmetrizeFrameworkCharges,sizeof(int),1,FilePtr);
  fread(&UseChargesFromMOLFile,sizeof(int),1,FilePtr);
  fread(&UseChargesFromCIFFile,sizeof(int),1,FilePtr);

  fread(&OverlapDistanceSquared,sizeof(OverlapDistanceSquared),1,FilePtr);
  fread(&CutOffVDW,sizeof(CutOffVDW),1,FilePtr);
  fread(&CutOffVDWSquared,sizeof(CutOffVDWSquared),1,FilePtr);
  fread(&CutOffVDWSwitch,sizeof(CutOffVDWSwitch),1,FilePtr);
  fread(&CutOffVDWSwitchSquared,sizeof(CutOffVDWSwitchSquared),1,FilePtr);
  fread(&InverseCutOffVDW,sizeof(InverseCutOffVDW),1,FilePtr);

  fread(CutOffChargeCharge,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CutOffChargeChargeSquared,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CutOffChargeChargeSwitch,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CutOffChargeChargeSwitchSquared,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(InverseCutOffChargeCharge,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(&CutOffChargeBondDipole,sizeof(CutOffChargeBondDipole),1,FilePtr);
  fread(&CutOffChargeBondDipoleSquared,sizeof(CutOffChargeBondDipoleSquared),1,FilePtr);
  fread(&CutOffChargeBondDipoleSwitch,sizeof(CutOffChargeBondDipoleSwitch),1,FilePtr);
  fread(&CutOffChargeBondDipoleSwitchSquared,sizeof(CutOffChargeBondDipoleSwitchSquared),1,FilePtr);

  fread(&CutOffBondDipoleBondDipole,sizeof(CutOffBondDipoleBondDipole),1,FilePtr);
  fread(&CutOffBondDipoleBondDipoleSquared,sizeof(CutOffBondDipoleBondDipoleSquared),1,FilePtr);
  fread(&CutOffBondDipoleBondDipoleSwitch,sizeof(CutOffBondDipoleBondDipoleSwitch),1,FilePtr);
  fread(&CutOffBondDipoleBondDipoleSwitchSquared,sizeof(CutOffBondDipoleBondDipoleSwitchSquared),1,FilePtr);

  fread(&DeltaT,sizeof(DeltaT),1,FilePtr);
  fread(Box,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(InverseBox,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(ReplicaBox,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(InverseReplicaBox,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(NumberOfReplicaCells,NumberOfSystems,sizeof(INT_VECTOR3),FilePtr);
  fread(TotalNumberOfReplicaCells,sizeof(int),NumberOfSystems,FilePtr);
  fread(UseReplicas,sizeof(int),NumberOfSystems,FilePtr);
  fread(BoxProperties,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(InverseBoxProperties,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);

  fread(Volume,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(AlphaAngle,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(BetaAngle,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(GammaAngle,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(BoundaryCondition,sizeof(int),NumberOfSystems,FilePtr);
  fread(IntialCenterOfMassPosition,sizeof(VECTOR),NumberOfSystems,FilePtr);
  fread(Beta,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(StrainDerivativeTailCorrection,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(StrainDerivativeRigidCorrection,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);


  fread(UseCellLists,sizeof(int),NumberOfSystems,FilePtr);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(UseCellLists[i])
    {
      fread(&NumberOfCellListCells[i].x,sizeof(int),1,FilePtr);
      fread(&NumberOfCellListCells[i].y,sizeof(int),1,FilePtr);
      fread(&NumberOfCellListCells[i].z,sizeof(int),1,FilePtr);
      fread(&NumberOfCellLists[i],sizeof(int),1,FilePtr);

      CellListMap[i]=(int (*)[27])calloc(NumberOfCellLists[i],sizeof(int[27]));
      fread(CellListMap[i],sizeof(int[27]),NumberOfCellLists[i],FilePtr);
    }
  }

  fread(UTotal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(ConservedEnergy,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(Drift,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(ReferenceEnergy,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UIon,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UDrift,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UTailCorrection,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UDistanceConstraints,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAngleConstraints,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UDihedralConstraints,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UInversionBendConstraints,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UOutOfPlaneDistanceConstraints,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UExclusionConstraints,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostPolarization,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbatePolarization,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationPolarization,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostBackPolarization,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateBackPolarization,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationBackPolarization,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostHost,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateCation,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostHostCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostHostVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostAdsorbateCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostAdsorbateVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostCationCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostCationVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateAdsorbateCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateAdsorbateVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateCationCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateCationVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationCationCoulomb,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationCationVDW,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostHostChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostAdsorbateChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostCationChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateAdsorbateChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateCationChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationCationChargeChargeReal,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostHostChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostAdsorbateChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostCationChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateAdsorbateChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateCationChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationCationChargeBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostHostBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostAdsorbateBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostCationBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateAdsorbateBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateCationBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationCationBondDipoleBondDipoleReal,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostHostChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostAdsorbateChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostCationChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateAdsorbateChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateCationChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationCationChargeChargeFourier,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostHostChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostAdsorbateChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostCationChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateAdsorbateChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateCationChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationCationChargeBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostHostBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostAdsorbateBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostCationBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateAdsorbateBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateCationBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationCationBondDipoleBondDipoleFourier,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UHostBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostUreyBradley,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostInversionBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostImproperTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostOutOfPlane,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostBondBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostBondBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostBendBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostBondTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostBendTorsion,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UCationBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationUreyBradley,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationInversionBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationImproperTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationOutOfPlane,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationBondBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationBondBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationBendBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationBondTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationBendTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationIntraVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationIntraChargeCharge,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationIntraChargeBondDipole,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationIntraBondDipoleBondDipole,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UAdsorbateBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateUreyBradley,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateInversionBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateImproperTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateOutOfPlane,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateBondBond,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateBondBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateBendBend,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateBondTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateBendTorsion,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateIntraVDW,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateIntraChargeCharge,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateIntraChargeBondDipole,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateIntraBondDipoleBondDipole,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(Ensemble,sizeof(int),NumberOfSystems,FilePtr);
  fread(InitEnsemble,sizeof(int),NumberOfSystems,FilePtr);
  fread(RunEnsemble,sizeof(int),NumberOfSystems,FilePtr);
  fread(NPTPRCellType,sizeof(int),NumberOfSystems,FilePtr);
  fread(MonoclinicAngleType,sizeof(int),NumberOfSystems,FilePtr);
  fread(MCEnsemble,sizeof(int),NumberOfSystems,FilePtr);

  fread(DegreesOfFreedom,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomTranslation,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomRotation,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomVibration,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomConstraint,sizeof(int),NumberOfSystems,FilePtr);

  fread(DegreesOfFreedomFramework,sizeof(int),NumberOfSystems,FilePtr);

  fread(DegreesOfFreedomAdsorbates,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomTranslationalAdsorbates,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomRotationalAdsorbates,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomVibrationalAdsorbates,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomConstraintAdsorbates,sizeof(int),NumberOfSystems,FilePtr);

  fread(DegreesOfFreedomCations,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomTranslationalCations,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomRotationalCations,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomVibrationalCations,sizeof(int),NumberOfSystems,FilePtr);
  fread(DegreesOfFreedomConstraintCations,sizeof(int),NumberOfSystems,FilePtr);

  fread(UKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UHostKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationKinetic,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UAdsorbateRotationalKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UAdsorbateTranslationalKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationRotationalKinetic,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UCationTranslationalKinetic,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(UNoseHoover,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UNoseHooverAdsorbates,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UNoseHooverCations,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(UNoseHooverFramework,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(PressureTrace,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(Pressure,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(StressTensor,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(ConfigurationalStressTensor,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(KineticStressTensor,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(StrainDerivativeTensor,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(MolecularStressTensor,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);

  fread(BornTerm,sizeof(REAL_MATRIX9x9),NumberOfSystems,FilePtr);
  fread(RelaxationTerm,sizeof(REAL_MATRIX9x9),NumberOfSystems,FilePtr);

  fread(&ProbabilityVolumeChangeMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityBoxShapeChangeMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityGibbsVolumeChangeMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityParallelTemperingMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityHyperParallelTemperingMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityParallelMolFractionMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityChiralInversionMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityHybridNVEMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityHybridNPHMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityHybridNPHPRMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityFrameworkChangeMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityFrameworkShiftMove,sizeof(REAL),1,FilePtr);
  fread(&ProbabilityCFCRXMCLambdaChangeMove,sizeof(REAL),1,FilePtr);

  fread(&CpuTimeProductionRun,sizeof(REAL),1,FilePtr);
  fread(&CpuTimeInitialization,sizeof(REAL),1,FilePtr);
  fread(&CpuTimeEquilibration,sizeof(REAL),1,FilePtr);
  fread(&CpuTotal,sizeof(REAL),1,FilePtr);

  fread(CpuTimeParallelTemperingMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeHyperParallelTemperingMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeParallelMolFractionMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeChiralInversionMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeHybridNVEMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeHybridNPHMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeHybridNPHPRMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeVolumeChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeBoxShapeChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeGibbsVolumeChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeFrameworkChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuTimeFrameworkShiftMove,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(CpuCFCRXMCLambdaChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(&OptimizeVolumeChange,sizeof(int),1,FilePtr);
  fread(&OptimizeGibbsVolumeChange,sizeof(int),1,FilePtr);
  fread(&OptimizeTranslation,sizeof(int),1,FilePtr);
  fread(&OptimizeRotation,sizeof(int),1,FilePtr);
  fread(&OptimizeFrameworkChange,sizeof(int),1,FilePtr);
  fread(&OptimizeFrameworkShift,sizeof(int),1,FilePtr);
  fread(&OptimizeCFLambdaChange,sizeof(int),1,FilePtr);
  fread(&OptimizeCFGibbsLambdaChange,sizeof(int),1,FilePtr);
  fread(&OptimizeCBCFLambdaChange,sizeof(int),1,FilePtr);
  fread(&OptimizeCBCFGibbsLambdaChange,sizeof(int),1,FilePtr);
  fread(&OptimizeRXMCLambdaChange,sizeof(int),1,FilePtr);

  // CFRXMC data
  fread(&NumberOfReactions,sizeof(int),1,FilePtr);
  if(NumberOfReactions>0)
  {
    fread(&RXMCLambdaHistogramSize,sizeof(int),1,FilePtr);

    CFCRXMCLambda=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    RXMCBiasingFactors[i]=(REAL**)calloc(NumberOfReactions,sizeof(REAL*));
    CFRXMCWangLandauScalingFactor=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    for(i=0;i<NumberOfSystems;i++)
    {
      CFCRXMCLambda[i]=(REAL*)calloc(NumberOfReactions,sizeof(REAL));
      fread(CFCRXMCLambda[i],sizeof(REAL),NumberOfReactions,FilePtr);

      CFRXMCWangLandauScalingFactor[i]=(REAL*)calloc(NumberOfReactions,sizeof(REAL));
      fread(CFRXMCWangLandauScalingFactor[i],sizeof(REAL),NumberOfReactions,FilePtr);
      RXMCBiasingFactors[i]=(REAL**)calloc(NumberOfReactions,sizeof(REAL*));
      for(j=0;j<NumberOfReactions;j++)
      {
        RXMCBiasingFactors[i][j]=(REAL*)calloc(RXMCLambdaHistogramSize,sizeof(REAL));
        fread(RXMCBiasingFactors[i][j],sizeof(REAL),RXMCLambdaHistogramSize,FilePtr);
      }
    }
  }

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartSimulation)\n");
    ContinueAfterCrash=FALSE;
  }
}
