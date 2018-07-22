/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'minimization.c' is part of RASPA-2.0

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "integration.h"
#include "simulation.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "inter_force.h"
#include "internal_force.h"
#include "molecule.h"
#include "input.h"
#include "output.h"
#include "mc_moves.h"
#include "statistics.h"
#include "potentials.h"
#include "cbmc.h"
#include "grids.h"
#include "movies.h"
#include "spectra.h"
#include "sample.h"
#include "thermo_baro_stats.h"
#include "recrossing.h"
#include "rigid.h"
#include "matrix.h"
#include "minimization.h"
#include "numerical.h"
#include "spacegroup.h"
#include "ewald.h"
#include "inter_hessian.h"
#include "internal_hessian.h"
#include "framework_hessian.h"

int UseSymmetryInMinimization;

int **FrameworkFixedInitialization;
int *AdsorbateFixedInitialization;
int *CationFixedInitialization;

int **NumberOfFixedFrameworkAtoms;
int ***FixedFrameworkAtoms;
int **NumberOfFixedFrameworkAtomsX;
int ***FixedFrameworkAtomsX;
int **NumberOfFixedFrameworkAtomsY;
int ***FixedFrameworkAtomsY;
int **NumberOfFixedFrameworkAtomsZ;
int ***FixedFrameworkAtomsZ;

int **NumberOfActiveFrameworkAtoms;
int ***ActiveFrameworkAtoms;
int **NumberOfActiveFrameworkAtomsX;
int ***ActiveFrameworkAtomsX;
int **NumberOfActiveFrameworkAtomsY;
int ***ActiveFrameworkAtomsY;
int **NumberOfActiveFrameworkAtomsZ;
int ***ActiveFrameworkAtomsZ;

int *NumberOfFixedAdsorbateMolecules;
int **FixedAdsorbateMolecules;
int *NumberOfActiveAdsorbateMolecules;
int **ActiveAdsorbateMolecules;

int *NumberOfFixedAdsorbateAtoms;
PAIR **FixedAdsorbateAtoms;
int *NumberOfFixedAdsorbateAtomsX;
PAIR **FixedAdsorbateAtomsX;
int *NumberOfFixedAdsorbateAtomsY;
PAIR **FixedAdsorbateAtomsY;
int *NumberOfFixedAdsorbateAtomsZ;
PAIR **FixedAdsorbateAtomsZ;
int *NumberOfActiveAdsorbateAtoms;
PAIR **ActiveAdsorbateAtoms;
int *NumberOfActiveAdsorbateAtomsX;
PAIR **ActiveAdsorbateAtomsX;
int *NumberOfActiveAdsorbateAtomsY;
PAIR **ActiveAdsorbateAtomsY;
int *NumberOfActiveAdsorbateAtomsZ;
PAIR **ActiveAdsorbateAtomsZ;

int *NumberOfFixedAdsorbateGroups;
PAIR **FixedAdsorbateGroups;
int *NumberOfActiveAdsorbateGroups;
PAIR **ActiveAdsorbateGroups;

int *NumberOfFixedAdsorbateGroupsCenterOfMass;
PAIR **FixedAdsorbateGroupsCenterOfMass;
int *NumberOfFixedAdsorbateGroupsCenterOfMassX;
PAIR **FixedAdsorbateGroupsCenterOfMassX;
int *NumberOfFixedAdsorbateGroupsCenterOfMassY;
PAIR **FixedAdsorbateGroupsCenterOfMassY;
int *NumberOfFixedAdsorbateGroupsCenterOfMassZ;
PAIR **FixedAdsorbateGroupsCenterOfMassZ;

int *NumberOfActiveAdsorbateGroupsCenterOfMass;
PAIR **ActiveAdsorbateGroupsCenterOfMass;
int *NumberOfActiveAdsorbateGroupsCenterOfMassX;
PAIR **ActiveAdsorbateGroupsCenterOfMassX;
int *NumberOfActiveAdsorbateGroupsCenterOfMassY;
PAIR **ActiveAdsorbateGroupsCenterOfMassY;
int *NumberOfActiveAdsorbateGroupsCenterOfMassZ;
PAIR **ActiveAdsorbateGroupsCenterOfMassZ;

int *NumberOfFixedAdsorbateGroupsOrientation;
PAIR **FixedAdsorbateGroupsOrientation;
int *NumberOfFixedAdsorbateGroupsOrientationX;
PAIR **FixedAdsorbateGroupsOrientationX;
int *NumberOfFixedAdsorbateGroupsOrientationY;
PAIR **FixedAdsorbateGroupsOrientationY;
int *NumberOfFixedAdsorbateGroupsOrientationZ;
PAIR **FixedAdsorbateGroupsOrientationZ;

int *NumberOfActiveAdsorbateGroupsOrientation;
PAIR **ActiveAdsorbateGroupsOrientation;
int *NumberOfActiveAdsorbateGroupsOrientationX;
PAIR **ActiveAdsorbateGroupsOrientationX;
int *NumberOfActiveAdsorbateGroupsOrientationY;
PAIR **ActiveAdsorbateGroupsOrientationY;
int *NumberOfActiveAdsorbateGroupsOrientationZ;
PAIR **ActiveAdsorbateGroupsOrientationZ;

int *NumberOfFixedCationMolecules;
int **FixedCationMolecules;
int *NumberOfActiveCationMolecules;
int **ActiveCationMolecules;

int *NumberOfFixedCationAtoms;
PAIR **FixedCationAtoms;
int *NumberOfFixedCationAtomsX;
PAIR **FixedCationAtomsX;
int *NumberOfFixedCationAtomsY;
PAIR **FixedCationAtomsY;
int *NumberOfFixedCationAtomsZ;
PAIR **FixedCationAtomsZ;
int *NumberOfActiveCationAtoms;
PAIR **ActiveCationAtoms;
int *NumberOfActiveCationAtomsX;
PAIR **ActiveCationAtomsX;
int *NumberOfActiveCationAtomsY;
PAIR **ActiveCationAtomsY;
int *NumberOfActiveCationAtomsZ;
PAIR **ActiveCationAtomsZ;

int *NumberOfFixedCationGroups;
PAIR **FixedCationGroups;
int *NumberOfActiveCationGroups;
PAIR **ActiveCationGroups;

int *NumberOfFixedCationGroupsCenterOfMass;
PAIR **FixedCationGroupsCenterOfMass;
int *NumberOfFixedCationGroupsCenterOfMassX;
PAIR **FixedCationGroupsCenterOfMassX;
int *NumberOfFixedCationGroupsCenterOfMassY;
PAIR **FixedCationGroupsCenterOfMassY;
int *NumberOfFixedCationGroupsCenterOfMassZ;
PAIR **FixedCationGroupsCenterOfMassZ;

int *NumberOfActiveCationGroupsCenterOfMass;
PAIR **ActiveCationGroupsCenterOfMass;
int *NumberOfActiveCationGroupsCenterOfMassX;
PAIR **ActiveCationGroupsCenterOfMassX;
int *NumberOfActiveCationGroupsCenterOfMassY;
PAIR **ActiveCationGroupsCenterOfMassY;
int *NumberOfActiveCationGroupsCenterOfMassZ;
PAIR **ActiveCationGroupsCenterOfMassZ;

int *NumberOfFixedCationGroupsOrientation;
PAIR **FixedCationGroupsOrientation;
int *NumberOfFixedCationGroupsOrientationX;
PAIR **FixedCationGroupsOrientationX;
int *NumberOfFixedCationGroupsOrientationY;
PAIR **FixedCationGroupsOrientationY;
int *NumberOfFixedCationGroupsOrientationZ;
PAIR **FixedCationGroupsOrientationZ;

int *NumberOfActiveCationGroupsOrientation;
PAIR **ActiveCationGroupsOrientation;
int *NumberOfActiveCationGroupsOrientationX;
PAIR **ActiveCationGroupsOrientationX;
int *NumberOfActiveCationGroupsOrientationY;
PAIR **ActiveCationGroupsOrientationY;
int *NumberOfActiveCationGroupsOrientationZ;
PAIR **ActiveCationGroupsOrientationZ;


int NumberOfFixedAtomTypes;
char (*FixedAtomTypes)[256];
int NumberOfActiveAtomTypes;
char (*ActiveAtomTypes)[256];

// distance constraints between two arbitrary atoms
int *NumberOfDistanceConstraints;
ATOM* (**DistanceConstraints)[2];
REAL **DistanceConstraintParameter;
int (**DistanceDefinitions)[2][3];
VECTOR (**DistanceConstraintsDerivatives)[2];

// angular constraints between three arbitrary atoms
int *NumberOfAngleConstraints;
ATOM* (**AngleConstraints)[3];
REAL **AngleConstraintParameter;
int (**AngleDefinitions)[3][3];
VECTOR (**AngleConstraintsDerivatives)[3];

// dihedral constraints between four arbitrary atoms
int *NumberOfDihedralConstraints;
ATOM* (**DihedralConstraints)[4];
REAL **DihedralConstraintParameter;
int (**DihedralDefinitions)[4][3];
VECTOR (**DihedralConstraintsDerivatives)[4];

// improper dihedral constraints between four arbitrary atoms
int *NumberOfImproperDihedralConstraints;
ATOM* (**ImproperDihedralConstraints)[4];
REAL **ImproperDihedralConstraintParameter;
int (**ImproperDihedralDefinitions)[4][3];
VECTOR (**ImproperDihedralConstraintsDerivatives)[4];


// Inversion-bend constraints between four arbitrary atoms
int *NumberOfInversionBendConstraints;
ATOM* (**InversionBendConstraints)[4];
REAL **InversionBendConstraintParameter;
int (**InversionBendDefinitions)[4][3];
VECTOR (**InversionBendConstraintsDerivatives)[4];

// Out-of-plane-distance constraints between four arbitrary atoms
int *NumberOfOutOfPlaneDistanceConstraints;
ATOM* (**OutOfPlaneDistanceConstraints)[4];
REAL **OutOfPlaneDistanceConstraintParameter;
int (**OutOfPlaneDistanceDefinitions)[4][3];
VECTOR (**OutOfPlaneDistanceConstraintsDerivatives)[4];

// harmonic distance constraints between two arbitrary atoms
int *NumberOfHarmonicDistanceConstraints;
ATOM* (**HarmonicDistanceConstraints)[2];
REAL (**HarmonicDistanceConstraintParameters)[2];
int (**HarmonicDistanceDefinitions)[2][3];

// harmonic angular constraints between three arbitrary atoms
int *NumberOfHarmonicAngleConstraints;
ATOM* (**HarmonicAngleConstraints)[3];
REAL (**HarmonicAngleConstraintParameters)[2];
int (**HarmonicAngleDefinitions)[3][3];

// harmonic dihedral constraints between four arbitrary atoms
int *NumberOfHarmonicDihedralConstraints;
ATOM* (**HarmonicDihedralConstraints)[4];
REAL (**HarmonicDihedralConstraintParameters)[2];
int (**HarmonicDihedralDefinitions)[4][3];


// two point dihedral measuments between four arbitrary points, the first point is the mid-point between atom 1 and 2
int *NumberOfTwoPointDihedralDefinitions;
int (**TwoPointDihedralDefinitions)[6][3];
ATOM* (**TwoPointDihedrals)[6];

int NumberOfMinimizationVariables;                // the total number of variables used in the minimization
int NumberOfCellMinimizationVariables;            // the number of cell-variables used in the minimization (0-6)
int NumberOfCoordinatesMinimizationVariables;     // the total number of variables used for generalized coordinates in the minimization
int NumberOfPositionalMinimizationVariables;      // the total number of variables used for positions in the minimization
int NumberOfOrientationalMinimizationVariables;   // the total number of variables used for orientation in the minimization

int MinimizationMethod;
int MinimizationPotentialMethod;
int MinimizationType;
int ComputeLambda;
int MinimizationVariables;

int ShellIndex;
int ShellSize;
int CoreSize;

REAL MaximumStepLengthInput;
REAL MaximumStepLength;
REAL MinimizationConvergenceFactor;

int RemoveTranslationFromHessian;
int RemoveRotationFromHessian;

int NumberOfPositionVariables;
int NumberOfBoxVariables;

void SnymanMinimization(int np,int nb,REAL *x,int run);
void BakerMinimization(int np,int nb,REAL *x,int run);
void BakerSaddlePointSearch(int np,int nb,REAL *x,int run);
void NewtonRaphsonMinimization(int np,int nb,REAL *p,int run);
void ComputeElasticConstantsAfterMinimization(int np,int nb,REAL *x);

REAL CalculateBakerStepSize(int n,int NumberOfValidModes,int k,REAL *dx,REAL *GX,REAL *EigenValues,REAL_MATRIX HessianMatrix);
REAL CalculateBakerSaddlePointStepSize(int n,int NumberOfValidModes,int ISTEP,REAL *dx,REAL *GX,REAL *EigenValues,REAL_MATRIX HessianMatrix);
void ConvertHessianFromCartesianToFractional(int np,int nb,REAL *Gradient,REAL_MATRIX Hessian);
void FreeMinimizationLocalMemory(void);

int UseGradientInLineMinimization;

int MaximumNumberOfMinimizationSteps;
REAL RMSGradientTolerance;
REAL MaxGradientTolerance;
REAL FunctionDifferenceTolerance;

REAL_MATRIX3x3 StoredBox;
REAL_MATRIX3x3 StoredReplicaBox;
REAL_MATRIX3x3 StoredInverseBox;

VECTOR *DVecX;
VECTOR *DVecY;
VECTOR *DVecZ;
VECTOR *DDVecAX;
VECTOR *DDVecBY;
VECTOR *DDVecCZ;
VECTOR *DDVecAY;
VECTOR *DDVecAZ;
VECTOR *DDVecBZ;

REAL UnitCellDeformation;
int TransformUnitCell;
int ElasticConstantEnergyVolume;


// Minimizes the current system
void Minimization(int run)
{
  int i;
  REAL *x;
  REAL *Gradient;
  REAL_MATRIX Hessian;
  VECTOR MidPoint1,MidPoint2,dr;
  int ncell,k1,k2,k3;
  REAL delta,det;
  REAL denominator;
  REAL_MATRIX3x3 ReferenceBox,ReferenceReplicaBox,strain;
  REAL energy;
  REAL_MATRIX3x3 StrainFirstDerivative;

  AllocateMinimizationLocalMemory();

  // loop over all the framwork-atoms and fix the order and 'index' into the Hessian
  NumberOfMinimizationVariables=0;
  NumberOfCellMinimizationVariables=0;
  OrderNumberOfMinimiationVariables();

  switch(Ensemble[CurrentSystem])
  {
    case NVE:
    case NVT:
      NumberOfCellMinimizationVariables=0;
      break;
    case NPT:
    case NPH:
      NumberOfCellMinimizationVariables=1;
      NumberOfMinimizationVariables+=1;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ANISOTROPIC:
        case ISOTROPIC:
          NumberOfCellMinimizationVariables=Dimension;
          NumberOfMinimizationVariables+=Dimension;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          NumberOfCellMinimizationVariables=4;
          NumberOfMinimizationVariables+=4;
          BoundaryCondition[CurrentSystem]=TRICLINIC;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          NumberOfCellMinimizationVariables=Dimension*(Dimension+1)/2;
          NumberOfMinimizationVariables+=Dimension*(Dimension+1)/2;
          BoundaryCondition[CurrentSystem]=TRICLINIC;
          break;
      }
  }

  if(TransformUnitCell)
  {
    BoundaryCondition[CurrentSystem]=TRICLINIC;

    SaveFrameworkPositionsToReferenceValues();
    SaveAdsorbateAtomPositionsToReferenceValues();
    SaveCationAtomPositionsToReferenceValues();

    delta=UnitCellDeformation;

    ReferenceBox=Box[CurrentSystem];
    ReferenceReplicaBox=ReplicaBox[CurrentSystem];

    // Note these deformations are volume conserving
    switch(ElasticConstantEnergyVolume)
    {
      case ELASTIC_CONSTANT_ORTHO_C11:
        strain.ax=1.0;       strain.bx=0.5*delta; strain.cx=0.0;
        strain.ay=0.5*delta; strain.by=1.0;       strain.cy=0.0;
        strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0+SQR(delta)/(4.0-SQR(delta));
        break;
      case ELASTIC_CONSTANT_ORTHO_C12:
        denominator = pow(1.0-SQR(delta),1.0/3.0);
        strain.ax=(1.0+delta)/denominator; strain.bx=0.0;                     strain.cx=0.0;
        strain.ay=0.0;                     strain.by=(1.0-delta)/denominator; strain.cy=0.0;
        strain.az=0.0;                     strain.bz=0.0;                     strain.cz=1.0/denominator;
        break;
      case ELASTIC_CONSTANT_ORTHO_C13:
        denominator = pow(1.0-SQR(delta),1.0/3.0);
        strain.ax=(1.0+delta)/denominator; strain.bx=0.0;             strain.cx=0.0;
        strain.ay=0.0;                     strain.by=1.0/denominator; strain.cy=0.0;
        strain.az=0.0;                     strain.bz=0.0;             strain.cz=(1.0-delta)/denominator;
        break;
      case ELASTIC_CONSTANT_CUBIC_C44:
        strain.ax=1.0;       strain.bx=0.5*delta; strain.cx=0.0;
        strain.ay=0.5*delta; strain.by=1.0;       strain.cy=0.0;
        strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0+SQR(delta)/(4.0-SQR(delta));
        break;
      case ELASTIC_CONSTANT_CUBIC_CS:
        strain.ax=1.0+delta; strain.bx=0.0;       strain.cx=0.0;
        strain.ay=0.0;       strain.by=1.0-delta; strain.cy=0.0;
        strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0+SQR(delta)/(1.0-SQR(delta));
        break;
      case ELASTIC_CONSTANT_CUBIC_C66:
        strain.ax=1.0+delta; strain.bx=0.0;       strain.cx=0.0;
        strain.ay=0.0;       strain.by=1.0-delta; strain.cy=0.0;
        strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0+SQR(delta)/(1.0-SQR(delta));
        break;
      case ELASTIC_CONSTANT_HEXAGONAL_C44:
        strain.ax=1.0;       strain.bx=0.0;                             strain.cx=delta;
        strain.ay=0.0;       strain.by=1.0+SQR(delta)/(1.0-SQR(delta)); strain.cy=0.0;
        strain.az=delta;     strain.bz=0.0;                             strain.cz=1.0;
        break;
      case ELASTIC_CONSTANT_BULK_MODULUS:
        strain.ax=1.0+delta; strain.bx=0.0;       strain.cx=0.0;
        strain.ay=0.0;       strain.by=1.0+delta; strain.cy=0.0;
        strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0+delta;
        break;
      default:
        strain.ax=1.0+delta; strain.bx=0.0;       strain.cx=0.0;
        strain.ay=0.0;       strain.by=1.0+delta; strain.cy=0.0;
        strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0+delta;
        break;
    }
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,ReferenceBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,ReferenceReplicaBox);
    InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
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
    CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
    CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();

    SetupKVectors();
  }

  // allocate memory for the generalized coordinates, and derivatives
  x=(REAL*)calloc(NumberOfMinimizationVariables,sizeof(REAL));
  Gradient=(REAL*)calloc(NumberOfMinimizationVariables,sizeof(REAL));
  Hessian=CreateRealMatrix(NumberOfMinimizationVariables,NumberOfMinimizationVariables);

  fprintf(OutputFilePtr[CurrentSystem],"Total number of variables used in the minimization: %d\n",NumberOfMinimizationVariables);
  fprintf(OutputFilePtr[CurrentSystem],"Number of cell-variables used in the minimization: %d\n",NumberOfCellMinimizationVariables);
  fprintf(OutputFilePtr[CurrentSystem],"Total number of variables used for generalized coordinates in the minimization: %d\n",NumberOfCoordinatesMinimizationVariables);
  fprintf(OutputFilePtr[CurrentSystem],"Total number of variables used for positions in the minimization: %d\n",NumberOfPositionalMinimizationVariables);
  fprintf(OutputFilePtr[CurrentSystem],"Total number of variables used for orientation in the minimization: %d\n",NumberOfOrientationalMinimizationVariables);
  fprintf(OutputFilePtr[CurrentSystem],"\n");

  StoredBox=Box[CurrentSystem];
  StoredReplicaBox=ReplicaBox[CurrentSystem];
  StoredInverseBox=InverseBox[CurrentSystem];

  // place all the positions/orientations in a generalized position vector 'x'
  // Note: 'x' has been allocate with all elements initialized to zero, so that include possible strain
  // the routine 'CreateGeneralizedCoordinatesFromPositions' needs to be called with zero strain here
  CreateGeneralizedCoordinatesFromPositions(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x);


  EvaluateDerivatives(NumberOfCoordinatesMinimizationVariables,&energy,NULL,Hessian,&StrainFirstDerivative,FALSE,FALSE);

  fprintf(OutputFilePtr[CurrentSystem],"Minimization initial volume: %18.10f [A^3]\n",Volume[CurrentSystem]);
  fprintf(OutputFilePtr[CurrentSystem],"Minimization initial energy: %18.10f [K]\n",energy*ENERGY_TO_KELVIN);
  fprintf(OutputFilePtr[CurrentSystem],"Minimization initial strain: %18.10f [-]\n",UnitCellDeformation);
  fprintf(OutputFilePtr[CurrentSystem],"Minimization initial box: %18.10f %18.10f %18.10f[-]\n",(double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx);
  fprintf(OutputFilePtr[CurrentSystem],"                          %18.10f %18.10f %18.10f[-]\n",(double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy);
  fprintf(OutputFilePtr[CurrentSystem],"                          %18.10f %18.10f %18.10f[-]\n",(double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz);
  fprintf(OutputFilePtr[CurrentSystem],"\n");

  // select the appropriate minimization method
  switch(MinimizationMethod)
  {
    case CONJUGATE_GRADIENT_MINIMIZATION:
      //ConjugateGradient(np,nb,x,1e-5,ComputeEnergy,ComputeForce);
      break;
    case STEEPEST_DESCENT_MINIMIZATION:
      //SteepestDescent(np,nb,x,1e-5,ComputeEnergy,ComputeForce);
      break;
    case BFGS_MINIMIZATION:
      //BFGS(np,nb,x,1e-5,ComputeEnergy,ComputeForce,run);
      break;
    case BAKER_MINIMIZATION:
      BakerMinimization(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x,run);
      break;
    case NEWTON_RAPHSON_MINIMIZATION:
      NewtonRaphsonMinimization(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x,run);
      break;
   case BAKER_SADDLE_POINT:
      BakerSaddlePointSearch(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x,run);
      break;
   case SNYMAN_MINIMIZATION:
      SnymanMinimization(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x,run);
      break;
  }

  if((NumberOfDistanceConstraints[CurrentSystem]>0)||(NumberOfAngleConstraints[CurrentSystem]>0)||(NumberOfDihedralConstraints[CurrentSystem]>0)||
     (NumberOfImproperDihedralConstraints[CurrentSystem]>0)||(NumberOfInversionBendConstraints[CurrentSystem]>0)||(NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]>0))
  {
    fprintf(OutputFilePtr[CurrentSystem],"Status Hard Constraints\n");
    fprintf(OutputFilePtr[CurrentSystem],"=======================\n");
    for(i=0;i<NumberOfDistanceConstraints[CurrentSystem];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Distance constraint %d, target: %18.10f final value: %18.10f\n",i,
        DistanceConstraintParameter[CurrentSystem][i],
        ReturnBondDistance(DistanceConstraints[CurrentSystem][i][CurrentSystem]->Position,DistanceConstraints[CurrentSystem][i][1]->Position));
    }
    for(i=0;i<NumberOfAngleConstraints[CurrentSystem];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Bond-Angle constraint %d, target: %18.10f final value: %18.10f\n",i,
        AngleConstraintParameter[CurrentSystem][i]*RAD2DEG,
        ReturnBendAngle(AngleConstraints[CurrentSystem][i][CurrentSystem]->Position,AngleConstraints[CurrentSystem][i][1]->Position,
        AngleConstraints[CurrentSystem][i][2]->Position)*RAD2DEG);
    }
    for(i=0;i<NumberOfDihedralConstraints[CurrentSystem];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Dihedral angle constraint %d, target: %18.10f final value: %18.10f\n",i,
        DihedralConstraintParameter[CurrentSystem][i]*RAD2DEG,
        ReturnDihedralAngle(DihedralConstraints[CurrentSystem][i][CurrentSystem]->Position,DihedralConstraints[CurrentSystem][i][1]->Position,
             DihedralConstraints[CurrentSystem][i][2]->Position,DihedralConstraints[CurrentSystem][i][3]->Position)*RAD2DEG);
    }
    for(i=0;i<NumberOfImproperDihedralConstraints[CurrentSystem];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Improper dihedral angle constraint %d, target: %18.10f final value: %18.10f\n",i,
        ImproperDihedralConstraintParameter[CurrentSystem][i]*RAD2DEG,
        ReturnDihedralAngle(ImproperDihedralConstraints[CurrentSystem][i][CurrentSystem]->Position,ImproperDihedralConstraints[CurrentSystem][i][1]->Position,
             ImproperDihedralConstraints[CurrentSystem][i][2]->Position,ImproperDihedralConstraints[CurrentSystem][i][3]->Position)*RAD2DEG);
    }
    for(i=0;i<NumberOfInversionBendConstraints[CurrentSystem];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Inversion-bend angle constraint %d, target: %18.10f final value: %18.10f\n",i,
        InversionBendConstraintParameter[CurrentSystem][i]*RAD2DEG,
        ReturnInversionBendAngle(InversionBendConstraints[CurrentSystem][i][CurrentSystem]->Position,InversionBendConstraints[CurrentSystem][i][1]->Position,
             InversionBendConstraints[CurrentSystem][i][2]->Position,InversionBendConstraints[CurrentSystem][i][3]->Position)*RAD2DEG);
    }
    for(i=0;i<NumberOfOutOfPlaneDistanceConstraints[CurrentSystem];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Out-of-plane-distance constraint %d, target: %18.10f final value: %18.10f\n",i,
        OutOfPlaneDistanceConstraintParameter[CurrentSystem][i]*RAD2DEG,
        ReturnOutOfPlaneDistance(OutOfPlaneDistanceConstraints[CurrentSystem][i][CurrentSystem]->Position,OutOfPlaneDistanceConstraints[CurrentSystem][i][1]->Position,
             OutOfPlaneDistanceConstraints[CurrentSystem][i][2]->Position,OutOfPlaneDistanceConstraints[CurrentSystem][i][3]->Position)*RAD2DEG);
    }
    fprintf(OutputFilePtr[CurrentSystem],"\n");
  }


  if((NumberOfHarmonicDistanceConstraints[CurrentSystem]>0)||(NumberOfHarmonicAngleConstraints[CurrentSystem]>0)||(NumberOfHarmonicDihedralConstraints[CurrentSystem]>0))
  {
    fprintf(OutputFilePtr[CurrentSystem],"Status Harmonic Constraints\n");
    fprintf(OutputFilePtr[CurrentSystem],"===========================\n");
    for(i=0;i<NumberOfHarmonicDistanceConstraints[CurrentSystem];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Distance constraint %d, target: %18.10f final value: %18.10f\n",i,
        HarmonicDistanceConstraintParameters[CurrentSystem][i][1],
        ReturnBondDistance(HarmonicDistanceConstraints[CurrentSystem][i][CurrentSystem]->Position,HarmonicDistanceConstraints[CurrentSystem][i][1]->Position));
    }
    for(i=0;i<NumberOfHarmonicAngleConstraints[CurrentSystem];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Bond-Angle constraint %d, target: %18.10f final value: %18.10f\n",i,
        HarmonicAngleConstraintParameters[CurrentSystem][i][1],
        ReturnBendAngle(HarmonicAngleConstraints[CurrentSystem][i][CurrentSystem]->Position,HarmonicAngleConstraints[CurrentSystem][i][1]->Position,
        HarmonicAngleConstraints[CurrentSystem][i][2]->Position)*RAD2DEG);
    }
    for(i=0;i<NumberOfHarmonicDihedralConstraints[CurrentSystem];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Dihedral angle constraint %d, target: %18.10f final value: %18.10f\n",i,
        HarmonicDihedralConstraintParameters[CurrentSystem][i][1],
        ReturnDihedralAngle(HarmonicDihedralConstraints[CurrentSystem][i][CurrentSystem]->Position,HarmonicDihedralConstraints[CurrentSystem][i][1]->Position,
             HarmonicDihedralConstraints[CurrentSystem][i][2]->Position,HarmonicDihedralConstraints[CurrentSystem][i][3]->Position)*RAD2DEG);
    }
    fprintf(OutputFilePtr[CurrentSystem],"\n");
  }

  if(NumberOfTwoPointDihedralDefinitions[CurrentSystem]>0)
  {
    fprintf(OutputFilePtr[CurrentSystem],"Status dihedral measurements\n");
    fprintf(OutputFilePtr[CurrentSystem],"============================\n");
    for(i=0;i<NumberOfTwoPointDihedralDefinitions[CurrentSystem];i++)
    {
      dr.x=TwoPointDihedrals[CurrentSystem][i][1]->Position.x-TwoPointDihedrals[CurrentSystem][i][0]->Position.x;
      dr.y=TwoPointDihedrals[CurrentSystem][i][1]->Position.y-TwoPointDihedrals[CurrentSystem][i][0]->Position.y;
      dr.z=TwoPointDihedrals[CurrentSystem][i][1]->Position.z-TwoPointDihedrals[CurrentSystem][i][0]->Position.z;
      dr=ApplyReplicaBoundaryCondition(dr);
      MidPoint1.x=TwoPointDihedrals[CurrentSystem][i][0]->Position.x+0.5*dr.x;
      MidPoint1.y=TwoPointDihedrals[CurrentSystem][i][0]->Position.y+0.5*dr.y;
      MidPoint1.z=TwoPointDihedrals[CurrentSystem][i][0]->Position.z+0.5*dr.z;

      dr.x=TwoPointDihedrals[CurrentSystem][i][5]->Position.x-TwoPointDihedrals[CurrentSystem][i][4]->Position.x;
      dr.y=TwoPointDihedrals[CurrentSystem][i][5]->Position.y-TwoPointDihedrals[CurrentSystem][i][4]->Position.y;
      dr.z=TwoPointDihedrals[CurrentSystem][i][5]->Position.z-TwoPointDihedrals[CurrentSystem][i][4]->Position.z;
      dr=ApplyReplicaBoundaryCondition(dr);
      MidPoint2.x=TwoPointDihedrals[CurrentSystem][i][4]->Position.x+0.5*dr.x;
      MidPoint2.y=TwoPointDihedrals[CurrentSystem][i][4]->Position.y+0.5*dr.y;
      MidPoint2.z=TwoPointDihedrals[CurrentSystem][i][4]->Position.z+0.5*dr.z;

      fprintf(OutputFilePtr[CurrentSystem],"Measured Dihedral angle %d, final value: %18.10f\n",i,
        ReturnDihedralAngle(MidPoint1,TwoPointDihedrals[CurrentSystem][i][2]->Position,
                            TwoPointDihedrals[CurrentSystem][i][3]->Position,MidPoint2)*RAD2DEG);
    }
    fprintf(OutputFilePtr[CurrentSystem],"\n");
  }

  free(x);
  free(Gradient);
  DeleteRealMatrix(Hessian);

  FreeMinimizationLocalMemory();
}


void SetStrainToZero(int np,REAL *x)
{
  // set references boxes
  StoredBox=Box[CurrentSystem];
  StoredReplicaBox=ReplicaBox[CurrentSystem];
  StoredInverseBox=InverseBox[CurrentSystem];

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      x[np]=0.0;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          x[np]=0.0;
          break;
        case ANISOTROPIC:
          switch(Dimension)
          {
            case 3:
              x[np+2]=0.0;
            case 2:
              x[np+1]=0.0;
            case 1:
              x[np]=0.0;
              break;
          }
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          switch(Dimension)
          {
            case 3:
              x[np+5]=0.0;
              x[np+4]=0.0;
              x[np+3]=0.0;
            case 2:
              x[np+2]=0.0;
              x[np+1]=0.0;
            case 1:
              x[np]=0.0;
              break;
          }
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          x[np]=x[np+1]=x[np+2]=x[np+3]=0.0;
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
    default:
      break;
  }
}

void SetWeights(int np,REAL *Weights,REAL *Charges)
{
  int i,m,l,f1,A;
  INT_VECTOR3 index,index2;
  int MolType,AtomType;
  REAL Mass;

  // fill in the generalized coordinates array from the current positions
  index=UNDEFINED_INT_VECTOR3;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        AtomType=Framework[CurrentSystem].Atoms[f1][i].Type;
        Mass=PseudoAtoms[AtomType].Mass;
        index=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;
        if(PseudoAtoms[AtomType].CoreShell==CORE)
        {
          switch(Dimension)
          {
            case 3:
              if(index.z>=0) Weights[index.z]=1.0/sqrt(Mass);
            case 2:
              if(index.y>=0) Weights[index.y]=1.0/sqrt(Mass);
            case 1:
              if(index.x>=0) Weights[index.x]=1.0/sqrt(Mass);
              break;
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    MolType=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[MolType].NumberOfGroups;l++)
    {
      if(Components[MolType].Groups[l].Rigid) // rigid unit
      {
        index=Adsorbates[CurrentSystem][m].Groups[l].HessianIndex;
        index2=Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation;

        Mass=Components[MolType].Groups[l].Mass;
        switch(Dimension)
        {
          case 3:
            if(index.z>=0) Weights[index.z]=1.0/sqrt(Mass);
          case 2:
            if(index.y>=0) Weights[index.y]=1.0/sqrt(Mass);
          case 1:
            if(index.x>=0) Weights[index.x]=1.0/sqrt(Mass);
            break;
        }

        switch(Dimension)
        {
          case 3:
            if(index2.z>=0) Weights[index2.z]=1.0/sqrt(Mass);
          case 2:
            if(index2.y>=0) Weights[index2.y]=1.0/sqrt(Mass);
          case 1:
            if(index2.x>=0) Weights[index2.x]=1.0/sqrt(Mass);
            break;
        }

      }
      else // flexible unit
      {
        for(i=0;i<Components[MolType].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[MolType].Groups[l].Atoms[i];
          AtomType=Adsorbates[CurrentSystem][m].Atoms[A].Type;
          index=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
          Mass=PseudoAtoms[AtomType].Mass;
          if(PseudoAtoms[AtomType].CoreShell==CORE)
          {
            switch(Dimension)
            {
              case 3:
                if(index.z>=0) Weights[index.z]=1.0/sqrt(Mass);
              case 2:
                if(index.y>=0) Weights[index.y]=1.0/sqrt(Mass);
              case 1:
                if(index.x>=0) Weights[index.x]=1.0/sqrt(Mass);
                break;
            }
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    MolType=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[MolType].NumberOfGroups;l++)
    {
      if(Components[MolType].Groups[l].Rigid) // rigid unit
      {
        index=Cations[CurrentSystem][m].Groups[l].HessianIndex;
        index2=Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation;

        Mass=Components[MolType].Groups[l].Mass;
        switch(Dimension)
        {
          case 3:
            if(index.z>=0) Weights[index.z]=1.0/sqrt(Mass);
          case 2:
            if(index.y>=0) Weights[index.y]=1.0/sqrt(Mass);
          case 1:
            if(index.x>=0) Weights[index.x]=1.0/sqrt(Mass);
            break;
        }


        switch(Dimension)
        {
          case 3:
            if(index2.z>=0) Weights[index2.z]=1.0/sqrt(Mass);
          case 2:
            if(index2.y>=0) Weights[index2.y]=1.0/sqrt(Mass);
          case 1:
            if(index2.x>=0) Weights[index2.x]=1.0/sqrt(Mass);
            break;
        }

      }
      else // flexible unit
      {
        for(i=0;i<Components[MolType].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[MolType].Groups[l].Atoms[i];
          AtomType=Cations[CurrentSystem][m].Atoms[A].Type;
          index=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
          Mass=PseudoAtoms[AtomType].Mass;
          if(PseudoAtoms[AtomType].CoreShell==CORE)
          {
            switch(Dimension)
            {
              case 3:
                if(index.z>=0) Weights[index.z]=1.0/sqrt(Mass);
              case 2:
                if(index.y>=0) Weights[index.y]=1.0/sqrt(Mass);
              case 1:
                if(index.x>=0) Weights[index.x]=1.0/sqrt(Mass);
                break;
            }
          }
        }
      }
    }
  }
}

// Create the framework, adsorbate and cation index into gradient and Hessian
int OrderNumberOfMinimiationVariables(void)
{
  int i,m,l,f1,A;
  int index,Type,AtomType;
  int atomic_index;

  NumberOfCoordinatesMinimizationVariables=0;
  NumberOfPositionalMinimizationVariables=0;
  NumberOfOrientationalMinimizationVariables=0;

  // loop over all the framwork-atoms and fix the order and 'index' into the Hessian
  index=0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        AtomType=Framework[CurrentSystem].Atoms[f1][i].Type;
        if(PseudoAtoms[AtomType].CoreShell==CORE)
        {
          switch(Dimension)
          {
            case 3:
              if(Framework[CurrentSystem].Atoms[f1][i].Fixed.x)
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.x=-1;
              else
              {
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.x=index++;
                NumberOfMinimizationVariables++;
                NumberOfPositionalMinimizationVariables++;
                NumberOfCoordinatesMinimizationVariables++;
              }
              if(Framework[CurrentSystem].Atoms[f1][i].Fixed.y)
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.y=-1;
              else
              {
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.y=index++;
                NumberOfMinimizationVariables++;
                NumberOfPositionalMinimizationVariables++;
                NumberOfCoordinatesMinimizationVariables++;
              }
              if(Framework[CurrentSystem].Atoms[f1][i].Fixed.z)
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.z=-1;
              else
              {
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.z=index++;
                NumberOfMinimizationVariables++;
                NumberOfPositionalMinimizationVariables++;
                NumberOfCoordinatesMinimizationVariables++;
              }
              break;
            case 2:
              if(Framework[CurrentSystem].Atoms[f1][i].Fixed.x)
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.x=-1;
              else
              {
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.x=index++;
                NumberOfMinimizationVariables++;
                NumberOfPositionalMinimizationVariables++;
                NumberOfCoordinatesMinimizationVariables++;
              }
              if(Framework[CurrentSystem].Atoms[f1][i].Fixed.y)
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.y=-1;
              else
              {
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.y=index++;
                NumberOfMinimizationVariables++;
                NumberOfPositionalMinimizationVariables++;
                NumberOfCoordinatesMinimizationVariables++;
              }
              break;
            case 1:
              if(Framework[CurrentSystem].Atoms[f1][i].Fixed.x)
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.x=-1;
              else
              {
                Framework[CurrentSystem].Atoms[f1][i].HessianIndex.x=index++;
                NumberOfMinimizationVariables++;
                NumberOfPositionalMinimizationVariables++;
                NumberOfCoordinatesMinimizationVariables++;
              }
              break;
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      // FIX (maybe)
      //Adsorbates[CurrentSystem][m].Groups[l].HessianIndex=index;
      if(Components[Type].Groups[l].Rigid)
      {
        switch(Dimension)
        {
          case 3:
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedCenterOfMass.x)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.x=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedCenterOfMass.y)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.y=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.y=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedCenterOfMass.z)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.z=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.z=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            break;
          case 2:
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedCenterOfMass.x)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.x=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedCenterOfMass.y)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.y=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.y=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            break;
          case 1:
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedCenterOfMass.x)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.x=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndex.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            break;
        }

        switch(Dimension)
        {
          case 3:
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedOrientation.x)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedOrientation.y)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.y=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.y=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedOrientation.z)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.z=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.z=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            break;
          case 2:
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedOrientation.x)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedOrientation.y)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.y=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.y=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            break;
          case 1:
            if(Adsorbates[CurrentSystem][m].Groups[l].FixedOrientation.x)
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=-1;
            else
            {
              Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            break;
        }
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          AtomType=Adsorbates[CurrentSystem][m].Atoms[A].Type;

          if(PseudoAtoms[AtomType].CoreShell==CORE)
          {
            switch(Dimension)
            {
              case 3:
                if(Adsorbates[CurrentSystem][m].Atoms[A].Fixed.x)
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.x=-1;
                else
                {
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.x=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                if(Adsorbates[CurrentSystem][m].Atoms[A].Fixed.y)
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.y=-1;
                else
                {
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.y=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                if(Adsorbates[CurrentSystem][m].Atoms[A].Fixed.z)
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.z=-1;
                else
                {
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.z=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                break;
              case 2:
                if(Adsorbates[CurrentSystem][m].Atoms[A].Fixed.x)
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.x=-1;
                else
                {
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.x=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                if(Adsorbates[CurrentSystem][m].Atoms[A].Fixed.y)
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.y=-1;
                else
                {
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.y=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                break;
              case 1:
                if(Adsorbates[CurrentSystem][m].Atoms[A].Fixed.x)
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.x=-1;
                else
                {
                  Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex.x=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                break;
            }
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      // FIX (maybe)
      //Cations[CurrentSystem][m].Groups[l].HessianIndex=index;
      if(Components[Type].Groups[l].Rigid)
      {
        switch(Dimension)
        {
          case 3:
            if(Cations[CurrentSystem][m].Groups[l].FixedCenterOfMass.x)
              Cations[CurrentSystem][m].Groups[l].HessianIndex.x=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndex.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            if(Cations[CurrentSystem][m].Groups[l].FixedCenterOfMass.y)
              Cations[CurrentSystem][m].Groups[l].HessianIndex.y=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndex.y=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            if(Cations[CurrentSystem][m].Groups[l].FixedCenterOfMass.z)
              Cations[CurrentSystem][m].Groups[l].HessianIndex.z=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndex.z=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            break;
          case 2:
            if(Cations[CurrentSystem][m].Groups[l].FixedCenterOfMass.x)
              Cations[CurrentSystem][m].Groups[l].HessianIndex.x=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndex.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            if(Cations[CurrentSystem][m].Groups[l].FixedCenterOfMass.y)
              Cations[CurrentSystem][m].Groups[l].HessianIndex.y=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndex.y=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            break;
          case 1:
            if(Cations[CurrentSystem][m].Groups[l].FixedCenterOfMass.x)
              Cations[CurrentSystem][m].Groups[l].HessianIndex.x=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndex.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfPositionalMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
            }
            break;
        }

        switch(Dimension)
        {
          case 3:
            if(Cations[CurrentSystem][m].Groups[l].FixedOrientation.x)
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            if(Cations[CurrentSystem][m].Groups[l].FixedOrientation.y)
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.y=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.y=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            if(Cations[CurrentSystem][m].Groups[l].FixedOrientation.z)
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.z=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.z=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            break;
          case 2:
            if(Cations[CurrentSystem][m].Groups[l].FixedOrientation.x)
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            if(Cations[CurrentSystem][m].Groups[l].FixedOrientation.y)
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.y=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.y=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            break;
          case 1:
            if(Cations[CurrentSystem][m].Groups[l].FixedOrientation.x)
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=-1;
            else
            {
              Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation.x=index++;
              NumberOfMinimizationVariables++;
              NumberOfCoordinatesMinimizationVariables++;
              NumberOfOrientationalMinimizationVariables++;

            }
            break;
        }
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          AtomType=Cations[CurrentSystem][m].Atoms[A].Type;

          if(PseudoAtoms[AtomType].CoreShell==CORE)
          {
            switch(Dimension)
            {
              case 3:
                if(Cations[CurrentSystem][m].Atoms[A].Fixed.x)
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.x=-1;
                else
                {
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.x=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                if(Cations[CurrentSystem][m].Atoms[A].Fixed.y)
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.y=-1;
                else
                {
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.y=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                if(Cations[CurrentSystem][m].Atoms[A].Fixed.z)
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.z=-1;
                else
                {
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.z=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                break;
              case 2:
                if(Cations[CurrentSystem][m].Atoms[A].Fixed.x)
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.x=-1;
                else
                {
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.x=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                if(Cations[CurrentSystem][m].Atoms[A].Fixed.y)
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.y=-1;
                else
                {
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.y=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                break;
              case 1:
                if(Cations[CurrentSystem][m].Atoms[A].Fixed.x)
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.x=-1;
                else
                {
                  Cations[CurrentSystem][m].Atoms[A].HessianIndex.x=index++;
                  NumberOfMinimizationVariables++;
                  NumberOfPositionalMinimizationVariables++;
                  NumberOfCoordinatesMinimizationVariables++;
                }
                break;
            }
          }
        }
      }
    }
  }


  ShellIndex=index;
  CoreSize=index;

  // put the shells last
  // ===================
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        AtomType=Framework[CurrentSystem].Atoms[f1][i].Type;
        if(PseudoAtoms[AtomType].CoreShell==SHELL)
        {
          switch(Dimension)
          {
            case 3:
              Framework[CurrentSystem].Atoms[f1][i].HessianIndex.x=index++;
              Framework[CurrentSystem].Atoms[f1][i].HessianIndex.y=index++;
              Framework[CurrentSystem].Atoms[f1][i].HessianIndex.z=index++;
              break;
            case 2:
              Framework[CurrentSystem].Atoms[f1][i].HessianIndex.x=index++;
              Framework[CurrentSystem].Atoms[f1][i].HessianIndex.y=index++;
              break;
            case 1:
              Framework[CurrentSystem].Atoms[f1][i].HessianIndex.x=index++;
              break;
          }
          NumberOfMinimizationVariables+=Dimension;
          NumberOfPositionalMinimizationVariables+=Dimension;
          NumberOfCoordinatesMinimizationVariables+=Dimension;
        }
      }
    }
  }

  ShellSize=index-ShellIndex;

  atomic_index=0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        Framework[CurrentSystem].Atoms[f1][i].HessianAtomIndex=atomic_index++;
    }
  }

  atomic_index=0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    for(i=0;i<Adsorbates[CurrentSystem][m].NumberOfAtoms;i++)
      Adsorbates[CurrentSystem][m].Atoms[i].HessianAtomIndex=atomic_index++;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    for(i=0;i<Cations[CurrentSystem][m].NumberOfAtoms;i++)
      Cations[CurrentSystem][m].Atoms[i].HessianAtomIndex=atomic_index++;
  return index;
}

// Put the positions, orientations and cell-info into the generalized coordinate vector 'x'

void CreateGeneralizedCoordinatesFromPositions(int np,int nb,REAL *x)
{
  int i,m,f1,A,MolType,l;
  INT_VECTOR3 index,index2;
  VECTOR pos,com,EulerAxis,pos1;
  REAL temp;
  REAL_MATRIX3x3 Transform,InverseTransform;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // compute the new Transform matrix based on the new values for the Transform
      switch(Dimension)
      {
        case 2:
          Transform.ax=1.0+x[np]; Transform.bx=0.0;       Transform.cx=0.0;
          Transform.ay=0.0;       Transform.by=1.0+x[np]; Transform.cy=0.0;
          Transform.az=0.0;       Transform.bz=0.0;       Transform.cz=0.0;
          break;
        case 3:
          Transform.ax=1.0+x[np]; Transform.bx=0.0;       Transform.cx=0.0;
          Transform.ay=0.0;       Transform.by=1.0+x[np]; Transform.cy=0.0;
          Transform.az=0.0;       Transform.bz=0.0;       Transform.cz=1.0+x[np];
          break;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          temp=(x[np]+x[np+1]+x[np+2])/3.0;
          Transform.ax=1.0+temp; Transform.bx=0.0;      Transform.cx=0.0;
          Transform.ay=0.0;      Transform.by=1.0+temp; Transform.cy=0.0;
          Transform.az=0.0;      Transform.bz=0.0;      Transform.cz=1.0+temp;
          break;
        case ANISOTROPIC:
          switch(Dimension)
          {
            case 2:
              Transform.ax=1.0+x[np];  Transform.bx=0.0;         Transform.cx=0.0;
              Transform.ay=0.0;        Transform.by=1.0+x[np+1]; Transform.cy=0.0;
              Transform.az=0.0;        Transform.bz=0.0;         Transform.cz=1.0;
              break;
            case 3:
              Transform.ax=1.0+x[np];  Transform.bx=0.0;         Transform.cx=0.0;
              Transform.ay=0.0;        Transform.by=1.0+x[np+1]; Transform.cy=0.0;
              Transform.az=0.0;        Transform.bz=0.0;         Transform.cz=1.0+x[np+2];
              break;
          }
          break;
        case REGULAR:
          switch(Dimension)
          {
            case 2:
              Transform.ax=1.0+x[np];     Transform.bx=0.5*x[np+1];     Transform.cx=0.0;
              Transform.ay=0.5*x[np+1];   Transform.by=1.0+x[np+2];     Transform.cy=0.0;
              Transform.az=0.0;           Transform.bz=0.0;             Transform.cz=0.0;
              break;
            case 3:
              Transform.ax=1.0+x[np];     Transform.bx=0.5*x[np+1];     Transform.cx=0.5*x[np+2];
              Transform.ay=0.5*x[np+1];   Transform.by=1.0+x[np+3];     Transform.cy=0.5*x[np+4];
              Transform.az=0.5*x[np+2];   Transform.bz=0.5*x[np+4];     Transform.cz=1.0+x[np+5];
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          switch(Dimension)
          {
            case 2:
              Transform.ax=1.0+x[np];     Transform.bx=x[np+1];         Transform.cx=0.0;
              Transform.ay=0;             Transform.by=1.0+x[np+2];     Transform.cy=0.0;
              Transform.az=0;             Transform.bz=0;               Transform.cz=0.0;
              break;
            case 3:
              Transform.ax=1.0+x[np];     Transform.bx=x[np+1];         Transform.cx=x[np+2];
              Transform.ay=0;             Transform.by=1.0+x[np+3];     Transform.cy=x[np+4];
              Transform.az=0;             Transform.bz=0;               Transform.cz=1.0+x[np+5];
              break;
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Transform.ax=1.0+x[np];  Transform.bx=0;            Transform.cx=0.0;
              Transform.ay=0.0;        Transform.by=1.0+x[np+1];  Transform.cy=0.5*x[np+2];
              Transform.az=0;          Transform.bz=0.5*x[np+2];  Transform.cz=1.0+x[np+3];
              break;
            case MONOCLINIC_BETA_ANGLE:
              Transform.ax=1.0+x[np];  Transform.bx=0;            Transform.cx=0.5*x[np+1];
              Transform.ay=0.0;        Transform.by=1.0+x[np+2];  Transform.cy=0;
              Transform.az=0.5*x[np+1];Transform.bz=0.0;          Transform.cz=1.0+x[np+3];
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Transform.ax=1.0+x[np];   Transform.bx=0.5*x[np+1];  Transform.cx=0;
              Transform.ay=0.5*x[np+1]; Transform.by=1.0+x[np+2];  Transform.cy=0;
              Transform.az=0;           Transform.bz=0;            Transform.cz=1.0+x[np+3];
              break;
            default:
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Transform.ax=1.0+x[np];  Transform.bx=0;            Transform.cx=0.0;
              Transform.ay=0.0;        Transform.by=1.0+x[np+1];  Transform.cy=x[np+2];
              Transform.az=0;          Transform.bz=0.0;          Transform.cz=1.0+x[np+3];
              break;
            case MONOCLINIC_BETA_ANGLE:
              Transform.ax=1.0+x[np];  Transform.bx=0;            Transform.cx=x[np+1];
              Transform.ay=0.0;        Transform.by=1.0+x[np+2];  Transform.cy=0;
              Transform.az=0.0;        Transform.bz=0.0;          Transform.cz=1.0+x[np+3];
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Transform.ax=1.0+x[np];   Transform.bx=x[np+1];      Transform.cx=0;
              Transform.ay=0;           Transform.by=1.0+x[np+2];  Transform.cy=0;
              Transform.az=0;           Transform.bz=0;            Transform.cz=1.0+x[np+3];
              break;
            default:
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
    default:
      switch(Dimension)
      {
        case 2:
          Transform.ax=1.0;   Transform.bx=0.0;     Transform.cx=0.0;
          Transform.ay=0.0;   Transform.by=1.0;     Transform.cy=0.0;
          Transform.az=0.0;   Transform.bz=0.0;     Transform.cz=0.0;
          break;
        case 3:
          Transform.ax=1.0;   Transform.bx=0.0;     Transform.cx=0.0;
          Transform.ay=0.0;   Transform.by=1.0;     Transform.cy=0.0;
          Transform.az=0.0;   Transform.bz=0.0;     Transform.cz=1.0;
          break;
      }
      break;
  }

  InverseBoxMatrix(&Transform,&InverseTransform);

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        index=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        pos=Framework[CurrentSystem].Atoms[f1][i].Position;

        pos1.x=InverseTransform.ax*pos.x+InverseTransform.bx*pos.y+InverseTransform.cx*pos.z;
        pos1.y=InverseTransform.ay*pos.x+InverseTransform.by*pos.y+InverseTransform.cy*pos.z;
        pos1.z=InverseTransform.az*pos.x+InverseTransform.bz*pos.y+InverseTransform.cz*pos.z;

        if(MinimizationVariables==FRACTIONAL)
          pos=ConvertFromXYZtoABC(pos);
        switch(Dimension)
        {
          case 3:
            if(index.z>=0) x[index.z]=pos1.z;
          case 2:
            if(index.y>=0) x[index.y]=pos1.y;
          case 1:
            if(index.x>=0) x[index.x]=pos1.x;
            break;
        }
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    MolType=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[MolType].NumberOfGroups;l++)
    {
      if(Components[MolType].Groups[l].Rigid) // rigid unit
      {
        // use center-of-mass position as generalized coordinates for translation
        index=Adsorbates[CurrentSystem][m].Groups[l].HessianIndex;
        index2=Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation;

        // convert to fractional positions
        com=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition;

        if(MinimizationVariables==FRACTIONAL)
          com=ConvertFromXYZtoABC(com);

        pos1.x=InverseTransform.ax*com.x+InverseTransform.bx*com.y+InverseTransform.cx*com.z;
        pos1.y=InverseTransform.ay*com.x+InverseTransform.by*com.y+InverseTransform.cy*com.z;
        pos1.z=InverseTransform.az*com.x+InverseTransform.bz*com.y+InverseTransform.cz*com.z;

        switch(Dimension)
        {
          case 3:
            if(index.z>=0) x[index.z]=pos1.z;
          case 2:
            if(index.y>=0) x[index.y]=pos1.y;
          case 1:
            if(index.x>=0) x[index.x]=pos1.x;
            break;
        }


        EulerAxis=Adsorbates[CurrentSystem][m].Groups[l].EulerAxis;
        switch(Dimension)
        {
          case 3:
            if(index2.z>=0) x[index2.z]=EulerAxis.z;
          case 2:
            if(index2.y>=0) x[index2.y]=EulerAxis.y;
          case 1:
            if(index2.x>=0) x[index2.x]=EulerAxis.x;
            break;
        }
      }
      else // flexible unit
      {
        for(i=0;i<Components[MolType].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[MolType].Groups[l].Atoms[i];

          // use atomic position as generalized coordinates for translation
          index=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
          pos=Adsorbates[CurrentSystem][m].Atoms[A].Position;

          pos1.x=InverseTransform.ax*pos.x+InverseTransform.bx*pos.y+InverseTransform.cx*pos.z;
          pos1.y=InverseTransform.ay*pos.x+InverseTransform.by*pos.y+InverseTransform.cy*pos.z;
          pos1.z=InverseTransform.az*pos.x+InverseTransform.bz*pos.y+InverseTransform.cz*pos.z;

          if(MinimizationVariables==FRACTIONAL)
            pos1=ConvertFromXYZtoABC(pos);

          switch(Dimension)
          {
            case 3:
              if(index.z>=0) x[index.z]=pos1.z;
            case 2:
              if(index.y>=0) x[index.y]=pos1.y;
            case 1:
              if(index.x>=0) x[index.x]=pos1.x;
              break;
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    MolType=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[MolType].NumberOfGroups;l++)
    {
      if(Components[MolType].Groups[l].Rigid) // rigid unit
      {
        // use center-of-mass position as generalized coordinates for translation
        index=Cations[CurrentSystem][m].Groups[l].HessianIndex;
        index2=Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation;

        // convert to fractional positions
        com=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition;

        if(MinimizationVariables==FRACTIONAL)
          com=ConvertFromXYZtoABC(com);

        pos1.x=InverseTransform.ax*com.x+InverseTransform.bx*com.y+InverseTransform.cx*com.z;
        pos1.y=InverseTransform.ay*com.x+InverseTransform.by*com.y+InverseTransform.cy*com.z;
        pos1.z=InverseTransform.az*com.x+InverseTransform.bz*com.y+InverseTransform.cz*com.z;

        switch(Dimension)
        {
          case 3:
            if(index.z>=0) x[index.z]=pos1.z;
          case 2:
            if(index.y>=0) x[index.y]=pos1.y;
          case 1:
            if(index.x>=0) x[index.x]=pos1.x;
            break;
        }


        EulerAxis=Cations[CurrentSystem][m].Groups[l].EulerAxis;
        switch(Dimension)
        {
          case 3:
            if(index2.z>=0) x[index2.z]=EulerAxis.z;
          case 2:
            if(index2.y>=0) x[index2.y]=EulerAxis.y;
          case 1:
            if(index2.x>=0) x[index2.x]=EulerAxis.x;
            break;
        }
      }
      else // flexible unit
      {
        for(i=0;i<Components[MolType].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[MolType].Groups[l].Atoms[i];

          // use atomic position as generalized coordinates for translation
          index=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
          pos=Cations[CurrentSystem][m].Atoms[A].Position;

          pos1.x=InverseTransform.ax*pos.x+InverseTransform.bx*pos.y+InverseTransform.cx*pos.z;
          pos1.y=InverseTransform.ay*pos.x+InverseTransform.by*pos.y+InverseTransform.cy*pos.z;
          pos1.z=InverseTransform.az*pos.x+InverseTransform.bz*pos.y+InverseTransform.cz*pos.z;

          if(MinimizationVariables==FRACTIONAL)
            pos1=ConvertFromXYZtoABC(pos);

          switch(Dimension)
          {
            case 3:
              if(index.z>=0) x[index.z]=pos1.z;
            case 2:
              if(index.y>=0) x[index.y]=pos1.y;
            case 1:
              if(index.x>=0) x[index.x]=pos1.x;
              break;
          }
        }
      }
    }
  }

}

// Get the positions, orientation and cell-info from the generalized vector 'x'.
// The vector 'x' is converted into the current state of the system.
// This function is used to compute the energy, gradients, and Hessian matrix at 'x'.
void CreatePositionsFromGeneralizedCoordinates(int np,int nb,REAL *x)
{
  int i,m,l,A,f1,Type;
  INT_VECTOR3 index,index2;
  REAL EulerAngle,det,temp;
  REAL_MATRIX3x3 RotationMatrix,Transform;
  VECTOR p,com,t,pos;
  int ncell,k1,k2,k3;

  Transform.ax=1.0;   Transform.bx=0.0;     Transform.cx=0.0;
  Transform.ay=0.0;   Transform.by=1.0;     Transform.cy=0.0;
  Transform.az=0.0;   Transform.bz=0.0;     Transform.cz=1.0;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // compute the new Transform matrix based on the new values for the Transform
      switch(Dimension)
      {
        case 2:
          Transform.ax=1.0+x[np]; Transform.bx=0.0;       Transform.cx=0.0;
          Transform.ay=0.0;       Transform.by=1.0+x[np]; Transform.cy=0.0;
          Transform.az=0.0;       Transform.bz=0.0;       Transform.cz=0.0;
          break;
        case 3:
          Transform.ax=1.0+x[np]; Transform.bx=0.0;       Transform.cx=0.0;
          Transform.ay=0.0;       Transform.by=1.0+x[np]; Transform.cy=0.0;
          Transform.az=0.0;       Transform.bz=0.0;       Transform.cz=1.0+x[np];
          break;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          temp=(x[np]+x[np+1]+x[np+2])/3.0;
          Transform.ax=1.0+temp; Transform.bx=0.0;      Transform.cx=0.0;
          Transform.ay=0.0;      Transform.by=1.0+temp; Transform.cy=0.0;
          Transform.az=0.0;      Transform.bz=0.0;      Transform.cz=1.0+temp;
          break;
        case ANISOTROPIC:
          switch(Dimension)
          {
            case 2:
              Transform.ax=1.0+x[np];  Transform.bx=0.0;         Transform.cx=0.0;
              Transform.ay=0.0;        Transform.by=1.0+x[np+1]; Transform.cy=0.0;
              Transform.az=0.0;        Transform.bz=0.0;         Transform.cz=1.0;
              break;
            case 3:
              Transform.ax=1.0+x[np];  Transform.bx=0.0;         Transform.cx=0.0;
              Transform.ay=0.0;        Transform.by=1.0+x[np+1]; Transform.cy=0.0;
              Transform.az=0.0;        Transform.bz=0.0;         Transform.cz=1.0+x[np+2];
              break;
          }
          break;
        case REGULAR:
          switch(Dimension)
          {
            case 2:
              Transform.ax=1.0+x[np];     Transform.bx=0.5*x[np+1];     Transform.cx=0.0;
              Transform.ay=0.5*x[np+1];   Transform.by=1.0+x[np+2];     Transform.cy=0.0;
              Transform.az=0.0;           Transform.bz=0.0;             Transform.cz=0.0;
              break;
            case 3:
              Transform.ax=1.0+x[np];     Transform.bx=0.5*x[np+1];     Transform.cx=0.5*x[np+2];
              Transform.ay=0.5*x[np+1];   Transform.by=1.0+x[np+3];     Transform.cy=0.5*x[np+4];
              Transform.az=0.5*x[np+2];   Transform.bz=0.5*x[np+4];     Transform.cz=1.0+x[np+5];
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          switch(Dimension)
          {
            case 2:
              Transform.ax=1.0+x[np];     Transform.bx=x[np+1];         Transform.cx=0.0;
              Transform.ay=0.0;           Transform.by=1.0+x[np+2];     Transform.cy=0.0;
              Transform.az=0.0;           Transform.bz=0.0;             Transform.cz=0.0;
              break;
            case 3:
              Transform.ax=1.0+x[np];     Transform.bx=x[np+1];         Transform.cx=x[np+2];
              Transform.ay=0.0;           Transform.by=1.0+x[np+3];     Transform.cy=x[np+4];
              Transform.az=0.0;           Transform.bz=0.0;             Transform.cz=1.0+x[np+5];
              break;
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Transform.ax=1.0+x[np];  Transform.bx=0;            Transform.cx=0.0;
              Transform.ay=0.0;        Transform.by=1.0+x[np+1];  Transform.cy=0.5*x[np+2];
              Transform.az=0;          Transform.bz=0.5*x[np+2];  Transform.cz=1.0+x[np+3];
              break;
            case MONOCLINIC_BETA_ANGLE:
              Transform.ax=1.0+x[np];  Transform.bx=0;            Transform.cx=0.5*x[np+1];
              Transform.ay=0.0;        Transform.by=1.0+x[np+2];  Transform.cy=0;
              Transform.az=0.5*x[np+1];Transform.bz=0.0;          Transform.cz=1.0+x[np+3];
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Transform.ax=1.0+x[np];   Transform.bx=0.5*x[np+1];  Transform.cx=0;
              Transform.ay=0.5*x[np+1]; Transform.by=1.0+x[np+2];  Transform.cy=0;
              Transform.az=0;           Transform.bz=0;            Transform.cz=1.0+x[np+3];
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Transform.ax=1.0+x[np];  Transform.bx=0;            Transform.cx=0.0;
              Transform.ay=0.0;        Transform.by=1.0+x[np+1];  Transform.cy=x[np+2];
              Transform.az=0;          Transform.bz=0.0;          Transform.cz=1.0+x[np+3];
              break;
            case MONOCLINIC_BETA_ANGLE:
              Transform.ax=1.0+x[np];  Transform.bx=0;            Transform.cx=x[np+1];
              Transform.ay=0.0;        Transform.by=1.0+x[np+2];  Transform.cy=0;
              Transform.az=0.0;        Transform.bz=0.0;          Transform.cz=1.0+x[np+3];
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Transform.ax=1.0+x[np];   Transform.bx=x[np+1];      Transform.cx=0;
              Transform.ay=0;           Transform.by=1.0+x[np+2];  Transform.cy=0;
              Transform.az=0;           Transform.bz=0;            Transform.cz=1.0+x[np+3];
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
    default:
      switch(Dimension)
      {
        case 2:
          Transform.ax=1.0;   Transform.bx=0.0;     Transform.cx=0.0;
          Transform.ay=0.0;   Transform.by=1.0;     Transform.cy=0.0;
          Transform.az=0.0;   Transform.bz=0.0;     Transform.cz=0.0;
          break;
        case 3:
          Transform.ax=1.0;   Transform.bx=0.0;     Transform.cx=0.0;
          Transform.ay=0.0;   Transform.by=1.0;     Transform.cy=0.0;
          Transform.az=0.0;   Transform.bz=0.0;     Transform.cz=1.0;
          break;
      }
      break;
  }

  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(Transform,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);

  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(Transform,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
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

  AlphaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bx);
  BetaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].by);
  GammaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bz);

  SetupKVectors();

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        index=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        pos.x=pos.y=pos.z=0.0;
        switch(Dimension)
        {
          case 3:
            if(index.z>=0) pos.z=x[index.z];
            else pos.z=Framework[CurrentSystem].Atoms[f1][i].Position.z;
          case 2:
            if(index.y>=0) pos.y=x[index.y];
            else pos.y=Framework[CurrentSystem].Atoms[f1][i].Position.y;
          case 1:
            if(index.x>=0) pos.x=x[index.x];
            else pos.x=Framework[CurrentSystem].Atoms[f1][i].Position.x;
            break;
        }


        if(MinimizationVariables==FRACTIONAL)
          Framework[CurrentSystem].Atoms[f1][i].Position=ConvertFromABCtoXYZ(pos);
        else
        {
          Framework[CurrentSystem].Atoms[f1][i].Position.x=Transform.ax*pos.x+Transform.bx*pos.y+Transform.cx*pos.z;
          Framework[CurrentSystem].Atoms[f1][i].Position.y=Transform.ay*pos.x+Transform.by*pos.y+Transform.cy*pos.z;
          Framework[CurrentSystem].Atoms[f1][i].Position.z=Transform.az*pos.x+Transform.bz*pos.y+Transform.cz*pos.z;
        }
      }
    }
  }
  else
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        pos=Framework[CurrentSystem].Atoms[f1][i].Position;
        if(MinimizationVariables==FRACTIONAL)
          Framework[CurrentSystem].Atoms[f1][i].Position=ConvertFromABCtoXYZ(pos);
        else
        {
          Framework[CurrentSystem].Atoms[f1][i].Position.x=Transform.ax*pos.x+Transform.bx*pos.y+Transform.cx*pos.z;
          Framework[CurrentSystem].Atoms[f1][i].Position.y=Transform.ay*pos.x+Transform.by*pos.y+Transform.cy*pos.z;
          Framework[CurrentSystem].Atoms[f1][i].Position.z=Transform.az*pos.x+Transform.bz*pos.y+Transform.cz*pos.z;
        }
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        index=Adsorbates[CurrentSystem][m].Groups[l].HessianIndex;
        index2=Adsorbates[CurrentSystem][m].Groups[l].HessianIndexOrientation;

        pos.x=pos.y=pos.z=0.0;
        switch(Dimension)
        {
          case 3:
            if(index.z>=0) pos.z=x[index.z];
            else pos.z=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z;
          case 2:
            if(index.y>=0) pos.y=x[index.y];
            else pos.y=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y;
          case 1:
            if(index.x>=0) pos.x=x[index.x];
            else pos.x=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x;
            break;
        }

        if(MinimizationVariables==FRACTIONAL)
          Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition=ConvertFromABCtoXYZ(pos);
        else
        {
          Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=Transform.ax*pos.x+Transform.bx*pos.y+Transform.cx*pos.z;
          Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=Transform.ay*pos.x+Transform.by*pos.y+Transform.cy*pos.z;
          Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=Transform.az*pos.x+Transform.bz*pos.y+Transform.cz*pos.z;
        }


        //Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.x=Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.y=Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.z=0.0;
        switch(Dimension)
        {
          case 3:
            if(index2.z>=0) Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.z=x[index2.z];
            //else Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.z=Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.z;
          case 2:
            if(index2.y>=0) Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.y=x[index2.y];
            //else Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.y=Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.y;
          case 1:
            if(index2.x>=0) Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.x=x[index2.x];
            //else Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.x=Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.x;
            break;
        }

        p=Adsorbates[CurrentSystem][m].Groups[l].EulerAxis;
        EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));

        if(EulerAngle<1e-8)
        {
          RotationMatrix.ax=1.0; RotationMatrix.bx=0.0; RotationMatrix.cx=0.0;
          RotationMatrix.ay=0.0; RotationMatrix.by=1.0; RotationMatrix.cy=0.0;
          RotationMatrix.az=0.0; RotationMatrix.bz=0.0; RotationMatrix.cz=1.0;
        }
        else
        {
          RotationMatrix.ax=1.0+(SQR(p.y)+SQR(p.z))*(cos(EulerAngle)-1.0)/SQR(EulerAngle);
          RotationMatrix.ay=(p.x*p.y - p.x*p.y*cos(EulerAngle) + p.z*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.az=(p.x*p.z - p.x*p.z*cos(EulerAngle) - p.y*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);

          RotationMatrix.bx=(p.x*p.y - p.x*p.y*cos(EulerAngle) - p.z*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.by=1.0+((SQR(p.x) + SQR(p.z))*(-1 + cos(EulerAngle)))/SQR(EulerAngle);
          RotationMatrix.bz=(p.y*p.z - p.y*p.z*cos(EulerAngle) + p.x*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);

          RotationMatrix.cx=(p.x*p.z - p.x*p.z*cos(EulerAngle) + p.y*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.cy=(p.y*p.z - p.y*p.z*cos(EulerAngle) - p.x*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.cz=1.0+((SQR(p.x) + SQR(p.y))*(-1.0 + cos(EulerAngle)))/SQR(EulerAngle);
        }

        com=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition;
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          pos=Components[Type].Positions[A];

          t.x=RotationMatrix.ax*pos.x+RotationMatrix.bx*pos.y+RotationMatrix.cx*pos.z;
          t.y=RotationMatrix.ay*pos.x+RotationMatrix.by*pos.y+RotationMatrix.cy*pos.z;
          t.z=RotationMatrix.az*pos.x+RotationMatrix.bz*pos.y+RotationMatrix.cz*pos.z;
          pos.x=com.x+t.x;
          pos.y=com.y+t.y;
          pos.z=com.z+t.z;
          Adsorbates[CurrentSystem][m].Atoms[A].Position=pos;
        }
      }
      else // flexible unit
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          index=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;

          pos.x=pos.y=pos.z=0.0;
          switch(Dimension)
          {
            case 3:
              if(index.z>=0) pos.z=x[index.z];
              else pos.z=Adsorbates[CurrentSystem][m].Atoms[A].Position.z;
            case 2:
              if(index.y>=0) pos.y=x[index.y];
              else pos.y=Adsorbates[CurrentSystem][m].Atoms[A].Position.y;
            case 1:
              if(index.x>=0) pos.x=x[index.x];
              else pos.x=Adsorbates[CurrentSystem][m].Atoms[A].Position.x;
              break;
          }

          if(MinimizationVariables==FRACTIONAL)
            Adsorbates[CurrentSystem][m].Atoms[A].Position=ConvertFromABCtoXYZ(pos);
          else
          {
            Adsorbates[CurrentSystem][m].Atoms[A].Position.x=Transform.ax*pos.x+Transform.bx*pos.y+Transform.cx*pos.z;
            Adsorbates[CurrentSystem][m].Atoms[A].Position.y=Transform.ay*pos.x+Transform.by*pos.y+Transform.cy*pos.z;
            Adsorbates[CurrentSystem][m].Atoms[A].Position.z=Transform.az*pos.x+Transform.bz*pos.y+Transform.cz*pos.z;
          }
        }
      }
    }
    UpdateGroupCenterOfMassAdsorbate(m);
    ComputeQuaternionAdsorbate(m);
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        index=Cations[CurrentSystem][m].Groups[l].HessianIndex;
        index2=Cations[CurrentSystem][m].Groups[l].HessianIndexOrientation;

        pos.x=pos.y=pos.z=0.0;
        switch(Dimension)
        {
          case 3:
            if(index.z>=0) pos.z=x[index.z];
            else pos.z=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z;
          case 2:
            if(index.y>=0) pos.y=x[index.y];
            else pos.y=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y;
          case 1:
            if(index.x>=0) pos.x=x[index.x];
            else pos.x=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x;
            break;
        }

        if(MinimizationVariables==FRACTIONAL)
          Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition=ConvertFromABCtoXYZ(pos);
        else
        {
          Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=Transform.ax*pos.x+Transform.bx*pos.y+Transform.cx*pos.z;
          Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=Transform.ay*pos.x+Transform.by*pos.y+Transform.cy*pos.z;
          Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=Transform.az*pos.x+Transform.bz*pos.y+Transform.cz*pos.z;
        }


        //Cations[CurrentSystem][m].Groups[l].EulerAxis.x=Cations[CurrentSystem][m].Groups[l].EulerAxis.y=Cations[CurrentSystem][m].Groups[l].EulerAxis.z=0.0;
        switch(Dimension)
        {
          case 3:
            if(index2.z>=0) Cations[CurrentSystem][m].Groups[l].EulerAxis.z=x[index2.z];
            //else Cations[CurrentSystem][m].Groups[l].EulerAxis.z=Cations[CurrentSystem][m].Groups[l].EulerAxis.z;
          case 2:
            if(index2.y>=0) Cations[CurrentSystem][m].Groups[l].EulerAxis.y=x[index2.y];
            //else Cations[CurrentSystem][m].Groups[l].EulerAxis.y=Cations[CurrentSystem][m].Groups[l].EulerAxis.y;
          case 1:
            if(index2.x>=0) Cations[CurrentSystem][m].Groups[l].EulerAxis.x=x[index2.x];
            //else Cations[CurrentSystem][m].Groups[l].EulerAxis.x=Cations[CurrentSystem][m].Groups[l].EulerAxis.x;
            break;
        }

        p=Cations[CurrentSystem][m].Groups[l].EulerAxis;
        EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));

        if(EulerAngle<1e-8)
        {
          RotationMatrix.ax=1.0; RotationMatrix.bx=0.0; RotationMatrix.cx=0.0;
          RotationMatrix.ay=0.0; RotationMatrix.by=1.0; RotationMatrix.cy=0.0;
          RotationMatrix.az=0.0; RotationMatrix.bz=0.0; RotationMatrix.cz=1.0;
        }
        else
        {
          RotationMatrix.ax=1.0+(SQR(p.y)+SQR(p.z))*(cos(EulerAngle)-1.0)/SQR(EulerAngle);
          RotationMatrix.ay=(p.x*p.y - p.x*p.y*cos(EulerAngle) + p.z*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.az=(p.x*p.z - p.x*p.z*cos(EulerAngle) - p.y*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);

          RotationMatrix.bx=(p.x*p.y - p.x*p.y*cos(EulerAngle) - p.z*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.by=1.0+((SQR(p.x) + SQR(p.z))*(-1 + cos(EulerAngle)))/SQR(EulerAngle);
          RotationMatrix.bz=(p.y*p.z - p.y*p.z*cos(EulerAngle) + p.x*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);

          RotationMatrix.cx=(p.x*p.z - p.x*p.z*cos(EulerAngle) + p.y*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.cy=(p.y*p.z - p.y*p.z*cos(EulerAngle) - p.x*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.cz=1.0+((SQR(p.x) + SQR(p.y))*(-1.0 + cos(EulerAngle)))/SQR(EulerAngle);
        }

        com=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition;
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          pos=Components[Type].Positions[A];

          t.x=RotationMatrix.ax*pos.x+RotationMatrix.bx*pos.y+RotationMatrix.cx*pos.z;
          t.y=RotationMatrix.ay*pos.x+RotationMatrix.by*pos.y+RotationMatrix.cy*pos.z;
          t.z=RotationMatrix.az*pos.x+RotationMatrix.bz*pos.y+RotationMatrix.cz*pos.z;
          pos.x=com.x+t.x;
          pos.y=com.y+t.y;
          pos.z=com.z+t.z;
          Cations[CurrentSystem][m].Atoms[A].Position=pos;
        }
      }
      else // flexible unit
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          index=Cations[CurrentSystem][m].Atoms[A].HessianIndex;

          pos.x=pos.y=pos.z=0.0;
          switch(Dimension)
          {
            case 3:
              if(index.z>=0) pos.z=x[index.z];
              else pos.z=Cations[CurrentSystem][m].Atoms[A].Position.z;
            case 2:
              if(index.y>=0) pos.y=x[index.y];
              else pos.y=Cations[CurrentSystem][m].Atoms[A].Position.y;
            case 1:
              if(index.x>=0) pos.x=x[index.x];
              else pos.x=Cations[CurrentSystem][m].Atoms[A].Position.x;
              break;
          }

          if(MinimizationVariables==FRACTIONAL)
            Cations[CurrentSystem][m].Atoms[A].Position=ConvertFromABCtoXYZ(pos);
          else
          {
            Cations[CurrentSystem][m].Atoms[A].Position.x=Transform.ax*pos.x+Transform.bx*pos.y+Transform.cx*pos.z;
            Cations[CurrentSystem][m].Atoms[A].Position.y=Transform.ay*pos.x+Transform.by*pos.y+Transform.cy*pos.z;
            Cations[CurrentSystem][m].Atoms[A].Position.z=Transform.az*pos.x+Transform.bz*pos.y+Transform.cz*pos.z;
          }
        }
      }
    }
    UpdateGroupCenterOfMassCation(m);
    ComputeQuaternionCation(m);
  }

}

// Compute the energy, generalized gradient, generalized Hessian matrix
// This routine expects that state 'x' is expanded into the positions and cell-info
void EvaluateDerivatives(int n,REAL *Energy,REAL* Gradient,REAL_MATRIX Hessian,REAL_MATRIX3x3 *StrainFirstDerivative,
                         int ComputeGradient,int ComputeHessian)
{
  REAL ExtPressure;

  *Energy=0.0;
  StrainFirstDerivative->ax=0.0; StrainFirstDerivative->bx=0.0; StrainFirstDerivative->cx=0.0;
  StrainFirstDerivative->ay=0.0; StrainFirstDerivative->by=0.0; StrainFirstDerivative->cy=0.0;
  StrainFirstDerivative->az=0.0; StrainFirstDerivative->bz=0.0; StrainFirstDerivative->cz=0.0;

  ExtPressure=therm_baro_stats.ExternalPressure[CurrentSystem][0];

  // compute the first and second derivatives of the rotation matrix
  PreComputeRotationDerivatives();

  CalculateAnisotropicSites();

  CalculateBondConstraintExclusionHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);

  // compute the intra framework contributions
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    ComputeFrameworkBondHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
    ComputeFrameworkUreyBradleyHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
    ComputeFrameworkBendHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
    ComputeFrameworkInversionBendHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
    ComputeFrameworkTorsionHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
    ComputeFrameworkImproperTorsionHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);

    ComputeFrameworkBendTorsionHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);

    ComputeFrameworkIntraVDWHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
    ComputeFrameworkIntraChargeChargeHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  }

  // compute the intra adsorbate contributions
  CalculateAdsorbateBondHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateAdsorbateBendHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateAdsorbateInversionBendHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateAdsorbateTorsionHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateAdsorbateImproperTorsionHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateAdsorbateIntraVDWHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateAdsorbateIntraCoulombHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateAdsorbateBondBondHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);

  // compute the intra cations contributions
  CalculateCationBondHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateCationBendHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateCationInversionBendHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateCationTorsionHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateCationImproperTorsionHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateCationIntraVDWHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateCationIntraCoulombHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);

  // compute the framework-adsorbate and framework-cations contributions
  ComputeFrameworkAdsorbateVDWHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  ComputeFrameworkAdsorbateChargeChargeHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  ComputeFrameworkCationVDWHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  ComputeFrameworkCationChargeChargeHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);

  // compute the adsorbate-adsorbate, cation-cation and adsorbate-cations intermolecular contributions
  ComputeInterVDWMolecularHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  ComputeInterChargeChargeMolecularHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);

  // compute the Ewald-summation contributions
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    CalculateEwaldFourierDerivatives(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);

  // these are not yet implemented analytically
  if(ComputeHessian)
  {
    AddRemainderOfCrossTermNumerically(Hessian);
    AddRemainderOfBornTermNumerically(Hessian);
  }

  // compute the harmonic constraints contributions
  CalculateHarmonicBondConstraintHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateHarmonicBendConstraintHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  CalculateHarmonicDihedralConstraintHessian(Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);

  if(ComputeGradient)
  {
    switch(Ensemble[CurrentSystem])
    {
      case NVE:
      case NVT:
        break;
      case NPT:
      case NPH:
        Gradient[NumberOfCoordinatesMinimizationVariables]+=ExtPressure*3.0*Volume[CurrentSystem];
        break;
      case NPTPR:
      case NPHPR:
        switch(NPTPRCellType[CurrentSystem])
        {
          case ISOTROPIC:
            break;
          case ANISOTROPIC:
            switch(Dimension)
            {
              case 3:
                Gradient[NumberOfCoordinatesMinimizationVariables+2]+=ExtPressure*Volume[CurrentSystem];
              case 2:
                Gradient[NumberOfCoordinatesMinimizationVariables+1]+=ExtPressure*Volume[CurrentSystem];
              case 1:
                Gradient[NumberOfCoordinatesMinimizationVariables]+=ExtPressure*Volume[CurrentSystem];
                break;
            }
            break;
          case MONOCLINIC:
          case MONOCLINIC_UPPER_TRIANGLE:
            switch(MonoclinicAngleType[CurrentSystem])
            {
              case MONOCLINIC_ALPHA_ANGLE:
                Gradient[NumberOfCoordinatesMinimizationVariables]+=ExtPressure*Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+1]+=ExtPressure*Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                break;
              case MONOCLINIC_BETA_ANGLE:
                Gradient[NumberOfCoordinatesMinimizationVariables]+=ExtPressure*Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+2]+=ExtPressure*Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                break;
              case MONOCLINIC_GAMMA_ANGLE:
                Gradient[NumberOfCoordinatesMinimizationVariables]+=ExtPressure*Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+2]+=ExtPressure*Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                break;
            }
            break;
          case REGULAR:
          case REGULAR_UPPER_TRIANGLE:
            if(therm_baro_stats.UseExternalStress)
            {
              switch(Dimension)
              {
                case 3:
                  Gradient[NumberOfCoordinatesMinimizationVariables]+=therm_baro_stats.ExternalStress[CurrentSystem].ax*Volume[CurrentSystem];
                  Gradient[NumberOfCoordinatesMinimizationVariables+1]+=therm_baro_stats.ExternalStress[CurrentSystem].bx*Volume[CurrentSystem];
                  Gradient[NumberOfCoordinatesMinimizationVariables+2]+=therm_baro_stats.ExternalStress[CurrentSystem].cx*Volume[CurrentSystem];
                  Gradient[NumberOfCoordinatesMinimizationVariables+3]+=therm_baro_stats.ExternalStress[CurrentSystem].by*Volume[CurrentSystem];
                  Gradient[NumberOfCoordinatesMinimizationVariables+4]+=therm_baro_stats.ExternalStress[CurrentSystem].cy*Volume[CurrentSystem];
                  Gradient[NumberOfCoordinatesMinimizationVariables+5]+=therm_baro_stats.ExternalStress[CurrentSystem].cz*Volume[CurrentSystem];
                  break;
                case 2:
                  Gradient[NumberOfCoordinatesMinimizationVariables]+=therm_baro_stats.ExternalStress[CurrentSystem].ax*Volume[CurrentSystem];
                  Gradient[NumberOfCoordinatesMinimizationVariables+1]+=therm_baro_stats.ExternalStress[CurrentSystem].bx*Volume[CurrentSystem];
                  Gradient[NumberOfCoordinatesMinimizationVariables+2]+=therm_baro_stats.ExternalStress[CurrentSystem].by*Volume[CurrentSystem];
                case 1:
                  break;
              }
            }
            else
            {
              switch(Dimension)
              {
                case 3:
                  Gradient[NumberOfCoordinatesMinimizationVariables+5]+=ExtPressure*Volume[CurrentSystem];
                case 2:
                  Gradient[NumberOfCoordinatesMinimizationVariables+Dimension]+=ExtPressure*Volume[CurrentSystem];
                case 1:
                  Gradient[NumberOfCoordinatesMinimizationVariables]+=ExtPressure*Volume[CurrentSystem];
                  break;
              }
            }
            break;
        }
    }
  }
  if(ComputeHessian)
  {
    switch(Ensemble[CurrentSystem])
    {
      case NVE:
      case NVT:
        break;
      case NPT:
      case NPH:
        Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]+=3.0*Dimension*ExtPressure*Volume[CurrentSystem];
        break;
      case NPTPR:
      case NPHPR:
        switch(NPTPRCellType[CurrentSystem])
        {
          case ISOTROPIC:
            break;
          case ANISOTROPIC:
            switch(Dimension)
            {
              case 3:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+2]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+1][NumberOfCoordinatesMinimizationVariables+2]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+2]+=ExtPressure*Volume[CurrentSystem];
              case 2:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+1]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+1][NumberOfCoordinatesMinimizationVariables+1]+=ExtPressure*Volume[CurrentSystem];
              case 1:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]+=ExtPressure*Volume[CurrentSystem];
                break;
            }
            break;
          case MONOCLINIC:
          case MONOCLINIC_UPPER_TRIANGLE:
            switch(MonoclinicAngleType[CurrentSystem])
            {
              case MONOCLINIC_ALPHA_ANGLE:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+1]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+1][NumberOfCoordinatesMinimizationVariables+1]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+1][NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                break;
              case MONOCLINIC_BETA_ANGLE:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+2]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+2]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                break;
              case MONOCLINIC_GAMMA_ANGLE:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+2]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+2]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+3]+=ExtPressure*Volume[CurrentSystem];
                break;
            }
            break;
          case REGULAR:
          case REGULAR_UPPER_TRIANGLE:
            if(therm_baro_stats.UseExternalStress)
            {
              switch(Dimension)
              {
                case 3:
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]+=therm_baro_stats.ExternalStress[CurrentSystem].ax*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+1]+=therm_baro_stats.ExternalStress[CurrentSystem].bx*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+2]+=therm_baro_stats.ExternalStress[CurrentSystem].cx*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+3]+=therm_baro_stats.ExternalStress[CurrentSystem].by*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+4]+=therm_baro_stats.ExternalStress[CurrentSystem].bz*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+5]+=therm_baro_stats.ExternalStress[CurrentSystem].cz*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+3]+=therm_baro_stats.ExternalStress[CurrentSystem].by*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+4]+=therm_baro_stats.ExternalStress[CurrentSystem].bz*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+5]+=therm_baro_stats.ExternalStress[CurrentSystem].cz*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables+5][NumberOfCoordinatesMinimizationVariables+5]+=therm_baro_stats.ExternalStress[CurrentSystem].cz*Volume[CurrentSystem];
                  break;
                case 2:
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]+=therm_baro_stats.ExternalStress[CurrentSystem].ax*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+1]+=therm_baro_stats.ExternalStress[CurrentSystem].bx*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+2]+=therm_baro_stats.ExternalStress[CurrentSystem].by*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+2]+=therm_baro_stats.ExternalStress[CurrentSystem].by*Volume[CurrentSystem];
                  break;
                case 1:
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]+=therm_baro_stats.ExternalStress[CurrentSystem].ax*Volume[CurrentSystem];
                  break;
              }
            }
            else
            {
              switch(Dimension)
              {
                case 3:
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+5]+=ExtPressure*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+5]+=ExtPressure*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables+5][NumberOfCoordinatesMinimizationVariables+5]+=ExtPressure*Volume[CurrentSystem];
                case 2:
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+Dimension]+=ExtPressure*Volume[CurrentSystem];
                  Hessian.element[NumberOfCoordinatesMinimizationVariables+Dimension][NumberOfCoordinatesMinimizationVariables+Dimension]+=ExtPressure*Volume[CurrentSystem];
                case 1:
                  Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]+=ExtPressure*Volume[CurrentSystem];
                  break;
              }
            }
            break;
        }
    }
  }



  int i,j;
  REAL temp;
  REAL energy,pressure;

  energy=pressure=0.0;
  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    for(j=0;j<NumberOfPseudoAtoms;j++)
    {
      if(TailCorrection[i][j])
      {
        temp=2.0*M_PI*(REAL)NumberOfPseudoAtomsType[CurrentSystem][i]*(REAL)NumberOfPseudoAtomsType[CurrentSystem][j]*
                   PotentialCorrection(i,j,CutOffVDW);
        pressure+=(2.0/3.0)*M_PI*NumberOfPseudoAtomsType[CurrentSystem][i]*NumberOfPseudoAtomsType[CurrentSystem][j]*PotentialCorrectionPressure(i,j,CutOffVDW);
        energy+=temp;
      }
    }
  }
  *Energy+=energy/Volume[CurrentSystem];

  StrainFirstDerivative->ax+=pressure/Volume[CurrentSystem];
  StrainFirstDerivative->by+=pressure/Volume[CurrentSystem];
  StrainFirstDerivative->cz+=pressure/Volume[CurrentSystem];

  // correct cell-gradients for tail-correction
  if(ComputeGradient)
  {
    switch(Ensemble[CurrentSystem])
    {
      case NVE:
      case NVT:
        break;
      case NPT:
      case NPH:
        Gradient[NumberOfCoordinatesMinimizationVariables]+=3.0*pressure/Volume[CurrentSystem];
        break;
      case NPTPR:
      case NPHPR:
        switch(NPTPRCellType[CurrentSystem])
        {
          case ANISOTROPIC:
          case ISOTROPIC:
            switch(Dimension)
            {
              case 3:
                Gradient[NumberOfCoordinatesMinimizationVariables+2]+=pressure/Volume[CurrentSystem];
              case 2:
                Gradient[NumberOfCoordinatesMinimizationVariables+1]+=pressure/Volume[CurrentSystem];
              case 1:
                Gradient[NumberOfCoordinatesMinimizationVariables]+=pressure/Volume[CurrentSystem];
                break;
            }
            break;
          case MONOCLINIC:
          case MONOCLINIC_UPPER_TRIANGLE:
            switch(MonoclinicAngleType[CurrentSystem])
            {
              case MONOCLINIC_ALPHA_ANGLE:
                Gradient[NumberOfCoordinatesMinimizationVariables]+=pressure/Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+1]+=pressure/Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+3]+=pressure/Volume[CurrentSystem];
                break;
              case MONOCLINIC_BETA_ANGLE:
                Gradient[NumberOfCoordinatesMinimizationVariables]+=pressure/Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+2]+=pressure/Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+3]+=pressure/Volume[CurrentSystem];
                break;
              case MONOCLINIC_GAMMA_ANGLE:
                Gradient[NumberOfCoordinatesMinimizationVariables]+=pressure/Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+2]+=pressure/Volume[CurrentSystem];
                Gradient[NumberOfCoordinatesMinimizationVariables+3]+=pressure/Volume[CurrentSystem];
                break;
            }
            break;
          case REGULAR:
          case REGULAR_UPPER_TRIANGLE:
            switch(Dimension)
            {
              case 3:
                Gradient[NumberOfCoordinatesMinimizationVariables+5]+=pressure/Volume[CurrentSystem];
              case 2:
                Gradient[NumberOfCoordinatesMinimizationVariables+Dimension]+=pressure/Volume[CurrentSystem];
              case 1:
                Gradient[NumberOfCoordinatesMinimizationVariables]+=pressure/Volume[CurrentSystem];
                break;
            }
            break;
        }
    }
  }

  // correct born-term for tail-correction
  if(ComputeHessian)
  {
    switch(Ensemble[CurrentSystem])
    {
      case NVE:
      case NVT:
        break;
      case NPT:
      case NPH:
        Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]-=3.0*Dimension*pressure/Volume[CurrentSystem];
        break;
      case NPTPR:
      case NPHPR:
        switch(NPTPRCellType[CurrentSystem])
        {
          case ANISOTROPIC:
          case ISOTROPIC:
            switch(Dimension)
            {
              case 3:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+2]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+1][NumberOfCoordinatesMinimizationVariables+2]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+2]-=pressure/Volume[CurrentSystem];
              case 2:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+1]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+1][NumberOfCoordinatesMinimizationVariables+1]-=pressure/Volume[CurrentSystem];
              case 1:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]-=pressure/Volume[CurrentSystem];
                break;
            }
            break;
          case MONOCLINIC:
          case MONOCLINIC_UPPER_TRIANGLE:
            switch(MonoclinicAngleType[CurrentSystem])
            {
              case MONOCLINIC_ALPHA_ANGLE:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+1]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+3]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+1][NumberOfCoordinatesMinimizationVariables+1]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+1][NumberOfCoordinatesMinimizationVariables+3]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+3]-=pressure/Volume[CurrentSystem];
                break;
              case MONOCLINIC_BETA_ANGLE:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+2]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+3]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+2]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+3]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+3]-=pressure/Volume[CurrentSystem];
                break;
              case MONOCLINIC_GAMMA_ANGLE:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+2]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+3]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+2]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+2][NumberOfCoordinatesMinimizationVariables+3]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+3]-=pressure/Volume[CurrentSystem];
                break;
            }
            break;
          case REGULAR:
          case REGULAR_UPPER_TRIANGLE:
            switch(Dimension)
            {
              case 3:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+5]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+3][NumberOfCoordinatesMinimizationVariables+5]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+5][NumberOfCoordinatesMinimizationVariables+5]-=pressure/Volume[CurrentSystem];
              case 2:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables+Dimension]-=pressure/Volume[CurrentSystem];
                Hessian.element[NumberOfCoordinatesMinimizationVariables+Dimension][NumberOfCoordinatesMinimizationVariables+Dimension]-=pressure/Volume[CurrentSystem];
              case 1:
                Hessian.element[NumberOfCoordinatesMinimizationVariables][NumberOfCoordinatesMinimizationVariables]-=pressure/Volume[CurrentSystem];
                break;
            }
            break;
        }
    }
  }
}

void ComputeDerivative(int np,int nb,REAL *x,REAL* Energy,REAL *Gradient,REAL_MATRIX3x3 *StrainFirstDerivative)
{
  int i;
  REAL_MATRIX Hessian;

  *Energy=0.0;
  for(i=0;i<np+nb;i++)
    Gradient[i]=0.0;

  StrainFirstDerivative->ax=0.0; StrainFirstDerivative->bx=0.0; StrainFirstDerivative->cx=0.0;
  StrainFirstDerivative->ay=0.0; StrainFirstDerivative->by=0.0; StrainFirstDerivative->cy=0.0;
  StrainFirstDerivative->az=0.0; StrainFirstDerivative->bz=0.0; StrainFirstDerivative->cz=0.0;

  // generate the positions of all the atoms from array 'x'
  CreatePositionsFromGeneralizedCoordinates(np,nb,x);

  // correct the positions for constraints
  ShakeInMinimization();

  CreateGeneralizedCoordinatesFromPositions(np,nb,x);

  EvaluateDerivatives(np,Energy,Gradient,Hessian,StrainFirstDerivative,TRUE,FALSE);


  // project the constraints from the gradient and Hessian
  ProjectConstraintsFromHessianMatrix(np,nb,Gradient,Hessian,TRUE,FALSE);
}



void ComputeDerivativesMinimization(int np,int nb,REAL *x,REAL* Energy,REAL *Gradient,REAL_MATRIX Hessian,REAL_MATRIX3x3 *StrainFirstDerivative)
{
  int i,j;

  *Energy=0.0;
  for(i=0;i<np+nb;i++)
  {
    Gradient[i]=0.0;
    for(j=0;j<np+nb;j++)
      Hessian.element[i][j]=0.0;
  }

  StrainFirstDerivative->ax=0.0; StrainFirstDerivative->bx=0.0; StrainFirstDerivative->cx=0.0;
  StrainFirstDerivative->ay=0.0; StrainFirstDerivative->by=0.0; StrainFirstDerivative->cy=0.0;
  StrainFirstDerivative->az=0.0; StrainFirstDerivative->bz=0.0; StrainFirstDerivative->cz=0.0;

  // generate the positions of all the atoms from array 'x'
  CreatePositionsFromGeneralizedCoordinates(np,nb,x);

  // correct the positions for constraints
  ShakeInMinimization();

  CreateGeneralizedCoordinatesFromPositions(np,nb,x);

  EvaluateDerivatives(np,Energy,Gradient,Hessian,StrainFirstDerivative,TRUE,TRUE);

  for(i=0;i<np+nb;i++)
    for(j=0;j<np+nb;j++)
      Hessian.element[j][i]=Hessian.element[i][j];
}

void ComputeDerivativesSpectra(int np,int nb,REAL *x,REAL* Energy,REAL *Gradient,REAL_MATRIX Hessian,REAL_MATRIX3x3 *StrainFirstDerivative)
{
  int i,j;

  *Energy=0.0;
  for(i=0;i<np+nb;i++)
  {
    Gradient[i]=0.0;
    for(j=0;j<np+nb;j++)
      Hessian.element[i][j]=0.0;
  }

  StrainFirstDerivative->ax=0.0; StrainFirstDerivative->bx=0.0; StrainFirstDerivative->cx=0.0;
  StrainFirstDerivative->ay=0.0; StrainFirstDerivative->by=0.0; StrainFirstDerivative->cy=0.0;
  StrainFirstDerivative->az=0.0; StrainFirstDerivative->bz=0.0; StrainFirstDerivative->cz=0.0;

  // generate the positions of all the atoms from array 'x'
  CreatePositionsFromGeneralizedCoordinates(np,nb,x);

  // correct the positions for constraints
  ShakeInMinimization();

  CreateGeneralizedCoordinatesFromPositions(np,nb,x);

  EvaluateDerivatives(np,Energy,Gradient,Hessian,StrainFirstDerivative,TRUE,TRUE);

  for(i=0;i<np+nb;i++)
    for(j=0;j<np+nb;j++)
      Hessian.element[j][i]=Hessian.element[i][j];

  if(MinimizationVariables==FRACTIONAL)
    ConvertHessianFromCartesianToFractional(np,nb,Gradient,Hessian);
}



void ConvertGradientFromCartesianToFractional(REAL *Gradient)
{
  int f1,i,l,m,Type,A,ia;
  INT_VECTOR3 index;
  VECTOR grad,fgrad;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        index=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        //if(index>=0)
        {
          grad.x=Gradient[index.x];
          grad.y=Gradient[index.y];
          grad.z=Gradient[index.z];
          fgrad.x=Box[CurrentSystem].ax*grad.x+Box[CurrentSystem].ay*grad.y+Box[CurrentSystem].az*grad.z;
          fgrad.y=Box[CurrentSystem].bx*grad.x+Box[CurrentSystem].by*grad.y+Box[CurrentSystem].bz*grad.z;
          fgrad.z=Box[CurrentSystem].cx*grad.x+Box[CurrentSystem].cy*grad.y+Box[CurrentSystem].cz*grad.z;
          Gradient[index.x]=fgrad.x;
          Gradient[index.y]=fgrad.y;
          Gradient[index.z]=fgrad.z;
        }
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        index=Adsorbates[CurrentSystem][m].Groups[l].HessianIndex;

        //if(index>=0)
        {
          grad.x=Gradient[index.x];
          grad.y=Gradient[index.y];
          grad.z=Gradient[index.z];
          fgrad.x=Box[CurrentSystem].ax*grad.x+Box[CurrentSystem].ay*grad.y+Box[CurrentSystem].az*grad.z;
          fgrad.y=Box[CurrentSystem].bx*grad.x+Box[CurrentSystem].by*grad.y+Box[CurrentSystem].bz*grad.z;
          fgrad.z=Box[CurrentSystem].cx*grad.x+Box[CurrentSystem].cy*grad.y+Box[CurrentSystem].cz*grad.z;
          Gradient[index.x]=fgrad.x;
          Gradient[index.y]=fgrad.y;
          Gradient[index.z]=fgrad.z;
        }
      }
      else
      {
        for(ia=0;ia<Components[Type].Groups[l].NumberOfGroupAtoms;ia++)
        {
          A=Components[Type].Groups[l].Atoms[ia];
          index=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;

          //if(index>=0)
          {
            grad.x=Gradient[index.x];
            grad.y=Gradient[index.y];
            grad.z=Gradient[index.z];
            fgrad.x=Box[CurrentSystem].ax*grad.x+Box[CurrentSystem].ay*grad.y+Box[CurrentSystem].az*grad.z;
            fgrad.y=Box[CurrentSystem].bx*grad.x+Box[CurrentSystem].by*grad.y+Box[CurrentSystem].bz*grad.z;
            fgrad.z=Box[CurrentSystem].cx*grad.x+Box[CurrentSystem].cy*grad.y+Box[CurrentSystem].cz*grad.z;
            Gradient[index.x]=fgrad.x;
            Gradient[index.y]=fgrad.y;
            Gradient[index.z]=fgrad.z;
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        index=Cations[CurrentSystem][m].Groups[l].HessianIndex;

        //if(index>=0)
        {
          grad.x=Gradient[index.x];
          grad.y=Gradient[index.y];
          grad.z=Gradient[index.z];
          fgrad.x=Box[CurrentSystem].ax*grad.x+Box[CurrentSystem].ay*grad.y+Box[CurrentSystem].az*grad.z;
          fgrad.y=Box[CurrentSystem].bx*grad.x+Box[CurrentSystem].by*grad.y+Box[CurrentSystem].bz*grad.z;
          fgrad.z=Box[CurrentSystem].cx*grad.x+Box[CurrentSystem].cy*grad.y+Box[CurrentSystem].cz*grad.z;
          Gradient[index.x]=fgrad.x;
          Gradient[index.y]=fgrad.y;
          Gradient[index.z]=fgrad.z;
        }
      }
      else
      {
        for(ia=0;ia<Components[Type].Groups[l].NumberOfGroupAtoms;ia++)
        {
          A=Components[Type].Groups[l].Atoms[ia];
          index=Cations[CurrentSystem][m].Atoms[A].HessianIndex;

          //if(index>=0)
          {
            grad.x=Gradient[index.x];
            grad.y=Gradient[index.y];
            grad.z=Gradient[index.z];
            fgrad.x=Box[CurrentSystem].ax*grad.x+Box[CurrentSystem].ay*grad.y+Box[CurrentSystem].az*grad.z;
            fgrad.y=Box[CurrentSystem].bx*grad.x+Box[CurrentSystem].by*grad.y+Box[CurrentSystem].bz*grad.z;
            fgrad.z=Box[CurrentSystem].cx*grad.x+Box[CurrentSystem].cy*grad.y+Box[CurrentSystem].cz*grad.z;
            Gradient[index.x]=fgrad.x;
            Gradient[index.y]=fgrad.y;
            Gradient[index.z]=fgrad.z;
          }
        }
      }
    }
  }
}


// correct for using fractional coordinates instead of Cartesian positions
void ConvertHessianFromCartesianToFractional(int np,int nb,REAL *Gradient,REAL_MATRIX Hessian)
{
  int I,J,ig,jg,ia,ja,i,j;
  int f1,f2;
  INT_VECTOR3 index_i,index_j,index;
  int m,l,A,Type;
  VECTOR grad,fgrad;
  REAL_MATRIX3x3 H,fH,tempBox;
  int TypeMolA,TypeMolB;

/*
  // correct the gradient
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        index=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        if(index>=0)
        {
          grad.x=Gradient[index]; grad.y=Gradient[index+1]; grad.z=Gradient[index+2];
          fgrad.x=Box[CurrentSystem].ax*grad.x+Box[CurrentSystem].ay*grad.y+Box[CurrentSystem].az*grad.z;
          fgrad.y=Box[CurrentSystem].bx*grad.x+Box[CurrentSystem].by*grad.y+Box[CurrentSystem].bz*grad.z;
          fgrad.z=Box[CurrentSystem].cx*grad.x+Box[CurrentSystem].cy*grad.y+Box[CurrentSystem].cz*grad.z;
          Gradient[index]=fgrad.x; Gradient[index+1]=fgrad.y; Gradient[index+2]=fgrad.z;
        }
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        index=Adsorbates[CurrentSystem][m].Groups[l].HessianIndex;
        if(index>=0)
        {
          grad.x=Gradient[index]; grad.y=Gradient[index+1]; grad.z=Gradient[index+2];
          fgrad.x=Box[CurrentSystem].ax*grad.x+Box[CurrentSystem].ay*grad.y+Box[CurrentSystem].az*grad.z;
          fgrad.y=Box[CurrentSystem].bx*grad.x+Box[CurrentSystem].by*grad.y+Box[CurrentSystem].bz*grad.z;
          fgrad.z=Box[CurrentSystem].cx*grad.x+Box[CurrentSystem].cy*grad.y+Box[CurrentSystem].cz*grad.z;
          Gradient[index]=fgrad.x; Gradient[index+1]=fgrad.y; Gradient[index+2]=fgrad.z;
        }
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          index=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;

          if(index>=0)
          {
            grad.x=Gradient[index]; grad.y=Gradient[index+1]; grad.z=Gradient[index+2];
            fgrad.x=Box[CurrentSystem].ax*grad.x+Box[CurrentSystem].ay*grad.y+Box[CurrentSystem].az*grad.z;
            fgrad.y=Box[CurrentSystem].bx*grad.x+Box[CurrentSystem].by*grad.y+Box[CurrentSystem].bz*grad.z;
            fgrad.z=Box[CurrentSystem].cx*grad.x+Box[CurrentSystem].cy*grad.y+Box[CurrentSystem].cz*grad.z;
            Gradient[index]=fgrad.x; Gradient[index+1]=fgrad.y; Gradient[index+2]=fgrad.z;
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        index=Cations[CurrentSystem][m].Groups[l].HessianIndex;

        if(index>=0)
        {
          grad.x=Gradient[index]; grad.y=Gradient[index+1]; grad.z=Gradient[index+2];
          fgrad.x=Box[CurrentSystem].ax*grad.x+Box[CurrentSystem].ay*grad.y+Box[CurrentSystem].az*grad.z;
          fgrad.y=Box[CurrentSystem].bx*grad.x+Box[CurrentSystem].by*grad.y+Box[CurrentSystem].bz*grad.z;
          fgrad.z=Box[CurrentSystem].cx*grad.x+Box[CurrentSystem].cy*grad.y+Box[CurrentSystem].cz*grad.z;
          Gradient[index]=fgrad.x; Gradient[index+1]=fgrad.y; Gradient[index+2]=fgrad.z;
        }
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          index=Cations[CurrentSystem][m].Atoms[A].HessianIndex;

          if(index>=0)
          {
            grad.x=Gradient[index]; grad.y=Gradient[index+1]; grad.z=Gradient[index+2];
            fgrad.x=Box[CurrentSystem].ax*grad.x+Box[CurrentSystem].ay*grad.y+Box[CurrentSystem].az*grad.z;
            fgrad.y=Box[CurrentSystem].bx*grad.x+Box[CurrentSystem].by*grad.y+Box[CurrentSystem].bz*grad.z;
            fgrad.z=Box[CurrentSystem].cx*grad.x+Box[CurrentSystem].cy*grad.y+Box[CurrentSystem].cz*grad.z;
            Gradient[index]=fgrad.x; Gradient[index+1]=fgrad.y; Gradient[index+2]=fgrad.z;
          }
        }
      }
    }
  }


  // note: strain derivative needs no correction

  // correct second derivatives
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        for(f2=0;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
        {
          for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
          {
            index_j=Framework[CurrentSystem].Atoms[f2][j].HessianIndex;

            if((index_i>=0)&&(index_j>=0))
            {
              H.ax=Hessian.element[index_i][index_j];   H.bx=Hessian.element[index_i+1][index_j];   H.cx=Hessian.element[index_i+2][index_j];
              H.ay=Hessian.element[index_i][index_j+1]; H.by=Hessian.element[index_i+1][index_j+1]; H.cy=Hessian.element[index_i+2][index_j+1];
              H.az=Hessian.element[index_i][index_j+2]; H.bz=Hessian.element[index_i+1][index_j+2]; H.cz=Hessian.element[index_i+2][index_j+2];

              tempBox=Box[CurrentSystem];
              TransposeMatrix3x3(&tempBox);
              fH=MatrixMatrixMultiplication3x3(tempBox,H);
              TransposeMatrix3x3(&fH);
              fH=MatrixMatrixMultiplication3x3(tempBox,fH);
              TransposeMatrix3x3(&fH);

              Hessian.element[index_i][index_j]=fH.ax;   Hessian.element[index_i+1][index_j]=fH.bx;    Hessian.element[index_i+2][index_j]=fH.cx;
              Hessian.element[index_i][index_j+1]=fH.ay; Hessian.element[index_i+1][index_j+1]=fH.by;  Hessian.element[index_i+2][index_j+1]=fH.cy;
              Hessian.element[index_i][index_j+2]=fH.az; Hessian.element[index_i+1][index_j+2]=fH.bz;  Hessian.element[index_i+2][index_j+2]=fH.cz;
            }
          }
        }

        for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
        {
          TypeMolB=Adsorbates[CurrentSystem][J].Type;
          for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
          {
            for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
            {
              j=Components[TypeMolB].Groups[jg].Atoms[ja];

              if(Components[TypeMolB].Groups[jg].Rigid)
                index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
              else
                index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;

              // adjust position-position
              if((index_i>=0)&&(index_j>=0))
              {
                H.ax=Hessian.element[index_i][index_j];   H.bx=Hessian.element[index_i+1][index_j];   H.cx=Hessian.element[index_i+2][index_j];
                H.ay=Hessian.element[index_i][index_j+1]; H.by=Hessian.element[index_i+1][index_j+1]; H.cy=Hessian.element[index_i+2][index_j+1];
                H.az=Hessian.element[index_i][index_j+2]; H.bz=Hessian.element[index_i+1][index_j+2]; H.cz=Hessian.element[index_i+2][index_j+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);
                TransposeMatrix3x3(&fH);
                fH=MatrixMatrixMultiplication3x3(tempBox,fH);
                TransposeMatrix3x3(&fH);

                Hessian.element[index_i][index_j]=fH.ax;   Hessian.element[index_i+1][index_j]=fH.bx;    Hessian.element[index_i+2][index_j]=fH.cx;
                Hessian.element[index_i][index_j+1]=fH.ay; Hessian.element[index_i+1][index_j+1]=fH.by;  Hessian.element[index_i+2][index_j+1]=fH.cy;
                Hessian.element[index_i][index_j+2]=fH.az; Hessian.element[index_i+1][index_j+2]=fH.bz;  Hessian.element[index_i+2][index_j+2]=fH.cz;
              }

              // adjust position-orientation
              if(Components[TypeMolB].Groups[jg].Rigid)
              {
                H.ax=Hessian.element[index_i+0][index_j+3]; H.bx=Hessian.element[index_i+1][index_j+3]; H.cx=Hessian.element[index_i+2][index_j+3];
                H.ay=Hessian.element[index_i+0][index_j+4]; H.by=Hessian.element[index_i+1][index_j+4]; H.cy=Hessian.element[index_i+2][index_j+4];
                H.az=Hessian.element[index_i+0][index_j+5]; H.bz=Hessian.element[index_i+1][index_j+5]; H.cz=Hessian.element[index_i+2][index_j+5];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                TransposeMatrix3x3(&H);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);
                TransposeMatrix3x3(&fH);

                Hessian.element[index_i+0][index_j+3]=fH.ax; Hessian.element[index_i+1][index_j+3]=fH.bx;  Hessian.element[index_i+2][index_j+3]=fH.cx;
                Hessian.element[index_i+0][index_j+4]=fH.ay; Hessian.element[index_i+1][index_j+4]=fH.by;  Hessian.element[index_i+2][index_j+4]=fH.cy;
                Hessian.element[index_i+0][index_j+5]=fH.az; Hessian.element[index_i+1][index_j+5]=fH.bz;  Hessian.element[index_i+2][index_j+5]=fH.cz;
              }

              // a rigid unit has only one index into the Hessian
              if(Components[TypeMolB].Groups[jg].Rigid) break;
            }
          }
        }
        for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
        {
          TypeMolB=Cations[CurrentSystem][J].Type;
          for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
          {
            for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
            {
              j=Components[TypeMolB].Groups[jg].Atoms[ja];

              if(Components[TypeMolB].Groups[jg].Rigid)
                index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
              else
                index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;

              // adjust position-position
              if((index_i>=0)&&(index_j>=0))
              {
                H.ax=Hessian.element[index_i][index_j];   H.bx=Hessian.element[index_i+1][index_j];   H.cx=Hessian.element[index_i+2][index_j];
                H.ay=Hessian.element[index_i][index_j+1]; H.by=Hessian.element[index_i+1][index_j+1]; H.cy=Hessian.element[index_i+2][index_j+1];
                H.az=Hessian.element[index_i][index_j+2]; H.bz=Hessian.element[index_i+1][index_j+2]; H.cz=Hessian.element[index_i+2][index_j+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);
                TransposeMatrix3x3(&fH);
                fH=MatrixMatrixMultiplication3x3(tempBox,fH);
                TransposeMatrix3x3(&fH);

                Hessian.element[index_i][index_j]=fH.ax;   Hessian.element[index_i+1][index_j]=fH.bx;    Hessian.element[index_i+2][index_j]=fH.cx;
                Hessian.element[index_i][index_j+1]=fH.ay; Hessian.element[index_i+1][index_j+1]=fH.by;  Hessian.element[index_i+2][index_j+1]=fH.cy;
                Hessian.element[index_i][index_j+2]=fH.az; Hessian.element[index_i+1][index_j+2]=fH.bz;  Hessian.element[index_i+2][index_j+2]=fH.cz;
              }

              // adjust position-orientation
              if(Components[TypeMolB].Groups[jg].Rigid)
              {
                H.ax=Hessian.element[index_i+0][index_j+3]; H.bx=Hessian.element[index_i+1][index_j+3]; H.cx=Hessian.element[index_i+2][index_j+3];
                H.ay=Hessian.element[index_i+0][index_j+4]; H.by=Hessian.element[index_i+1][index_j+4]; H.cy=Hessian.element[index_i+2][index_j+4];
                H.az=Hessian.element[index_i+0][index_j+5]; H.bz=Hessian.element[index_i+1][index_j+5]; H.cz=Hessian.element[index_i+2][index_j+5];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                TransposeMatrix3x3(&H);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);
                TransposeMatrix3x3(&fH);

                Hessian.element[index_i+0][index_j+3]=fH.ax; Hessian.element[index_i+1][index_j+3]=fH.bx;  Hessian.element[index_i+2][index_j+3]=fH.cx;
                Hessian.element[index_i+0][index_j+4]=fH.ay; Hessian.element[index_i+1][index_j+4]=fH.by;  Hessian.element[index_i+2][index_j+4]=fH.cy;
                Hessian.element[index_i+0][index_j+5]=fH.az; Hessian.element[index_i+1][index_j+5]=fH.bz;  Hessian.element[index_i+2][index_j+5]=fH.cz;
              }

              // a rigid unit has only one index into the Hessian
              if(Components[TypeMolB].Groups[jg].Rigid) break;
            }
          }
        }
      }
    }
  }


  for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
  {
    TypeMolA=Adsorbates[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(Components[TypeMolA].Groups[ig].Rigid)
          index_i=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndex;
        else
          index_i=Adsorbates[CurrentSystem][I].Atoms[i].HessianIndex;

        for(J=I;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
        {
          TypeMolB=Adsorbates[CurrentSystem][J].Type;
          for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
          {
            for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
            {
              j=Components[TypeMolB].Groups[jg].Atoms[ja];

              if(Components[TypeMolB].Groups[jg].Rigid)
                index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
              else
                index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;

              // adjust position-position
              if((index_i>=0)&&(index_j>=0))
              {
                H.ax=Hessian.element[index_i][index_j];   H.bx=Hessian.element[index_i+1][index_j];   H.cx=Hessian.element[index_i+2][index_j];
                H.ay=Hessian.element[index_i][index_j+1]; H.by=Hessian.element[index_i+1][index_j+1]; H.cy=Hessian.element[index_i+2][index_j+1];
                H.az=Hessian.element[index_i][index_j+2]; H.bz=Hessian.element[index_i+1][index_j+2]; H.cz=Hessian.element[index_i+2][index_j+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);
                TransposeMatrix3x3(&fH);
                fH=MatrixMatrixMultiplication3x3(tempBox,fH);
                TransposeMatrix3x3(&fH);

                Hessian.element[index_i][index_j]=fH.ax;   Hessian.element[index_i+1][index_j]=fH.bx;    Hessian.element[index_i+2][index_j]=fH.cx;
                Hessian.element[index_i][index_j+1]=fH.ay; Hessian.element[index_i+1][index_j+1]=fH.by;  Hessian.element[index_i+2][index_j+1]=fH.cy;
                Hessian.element[index_i][index_j+2]=fH.az; Hessian.element[index_i+1][index_j+2]=fH.bz;  Hessian.element[index_i+2][index_j+2]=fH.cz;
              }

              // adjust orientation-position
              if(Components[TypeMolA].Groups[ig].Rigid)
              {
                H.ax=Hessian.element[index_i+3][index_j+0]; H.bx=Hessian.element[index_i+4][index_j+0]; H.cx=Hessian.element[index_i+5][index_j+0];
                H.ay=Hessian.element[index_i+3][index_j+1]; H.by=Hessian.element[index_i+4][index_j+1]; H.cy=Hessian.element[index_i+5][index_j+1];
                H.az=Hessian.element[index_i+3][index_j+2]; H.bz=Hessian.element[index_i+4][index_j+2]; H.cz=Hessian.element[index_i+5][index_j+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+3][index_j+0]=fH.ax; Hessian.element[index_i+4][index_j+0]=fH.bx;  Hessian.element[index_i+5][index_j+0]=fH.cx;
                Hessian.element[index_i+3][index_j+1]=fH.ay; Hessian.element[index_i+4][index_j+1]=fH.by;  Hessian.element[index_i+5][index_j+1]=fH.cy;
                Hessian.element[index_i+3][index_j+2]=fH.az; Hessian.element[index_i+4][index_j+2]=fH.bz;  Hessian.element[index_i+5][index_j+2]=fH.cz;
              }

              // adjust position-orientation
              if(Components[TypeMolB].Groups[jg].Rigid)
              {
                H.ax=Hessian.element[index_i+0][index_j+3]; H.bx=Hessian.element[index_i+1][index_j+3]; H.cx=Hessian.element[index_i+2][index_j+3];
                H.ay=Hessian.element[index_i+0][index_j+4]; H.by=Hessian.element[index_i+1][index_j+4]; H.cy=Hessian.element[index_i+2][index_j+4];
                H.az=Hessian.element[index_i+0][index_j+5]; H.bz=Hessian.element[index_i+1][index_j+5]; H.cz=Hessian.element[index_i+2][index_j+5];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                TransposeMatrix3x3(&H);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);
                TransposeMatrix3x3(&fH);

                Hessian.element[index_i+0][index_j+3]=fH.ax; Hessian.element[index_i+1][index_j+3]=fH.bx;  Hessian.element[index_i+2][index_j+3]=fH.cx;
                Hessian.element[index_i+0][index_j+4]=fH.ay; Hessian.element[index_i+1][index_j+4]=fH.by;  Hessian.element[index_i+2][index_j+4]=fH.cy;
                Hessian.element[index_i+0][index_j+5]=fH.az; Hessian.element[index_i+1][index_j+5]=fH.bz;  Hessian.element[index_i+2][index_j+5]=fH.cz;
              }

              // note: no need to adjust orientation-orientation (independent of Cartesian or fractional part)

              // a rigid unit has only one index into the Hessian
              if(Components[TypeMolB].Groups[jg].Rigid) break;
            }
          }
        }

        for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
        {
          TypeMolB=Cations[CurrentSystem][J].Type;
          for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
          {
            for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
            {
              j=Components[TypeMolB].Groups[jg].Atoms[ja];

              if(Components[TypeMolB].Groups[jg].Rigid)
                index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
              else
                index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;

              // adjust position-position
              if((index_i>=0)&&(index_j>=0))
              {
                H.ax=Hessian.element[index_i][index_j];   H.bx=Hessian.element[index_i+1][index_j];   H.cx=Hessian.element[index_i+2][index_j];
                H.ay=Hessian.element[index_i][index_j+1]; H.by=Hessian.element[index_i+1][index_j+1]; H.cy=Hessian.element[index_i+2][index_j+1];
                H.az=Hessian.element[index_i][index_j+2]; H.bz=Hessian.element[index_i+1][index_j+2]; H.cz=Hessian.element[index_i+2][index_j+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);
                TransposeMatrix3x3(&fH);
                fH=MatrixMatrixMultiplication3x3(tempBox,fH);
                TransposeMatrix3x3(&fH);

                Hessian.element[index_i][index_j]=fH.ax;   Hessian.element[index_i+1][index_j]=fH.bx;    Hessian.element[index_i+2][index_j]=fH.cx;
                Hessian.element[index_i][index_j+1]=fH.ay; Hessian.element[index_i+1][index_j+1]=fH.by;  Hessian.element[index_i+2][index_j+1]=fH.cy;
                Hessian.element[index_i][index_j+2]=fH.az; Hessian.element[index_i+1][index_j+2]=fH.bz;  Hessian.element[index_i+2][index_j+2]=fH.cz;
              }

              // adjust orientation-position
              if(Components[TypeMolA].Groups[ig].Rigid)
              {
                H.ax=Hessian.element[index_i+3][index_j+0]; H.bx=Hessian.element[index_i+4][index_j+0]; H.cx=Hessian.element[index_i+5][index_j+0];
                H.ay=Hessian.element[index_i+3][index_j+1]; H.by=Hessian.element[index_i+4][index_j+1]; H.cy=Hessian.element[index_i+5][index_j+1];
                H.az=Hessian.element[index_i+3][index_j+2]; H.bz=Hessian.element[index_i+4][index_j+2]; H.cz=Hessian.element[index_i+5][index_j+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+3][index_j+0]=fH.ax; Hessian.element[index_i+4][index_j+0]=fH.bx;  Hessian.element[index_i+5][index_j+0]=fH.cx;
                Hessian.element[index_i+3][index_j+1]=fH.ay; Hessian.element[index_i+4][index_j+1]=fH.by;  Hessian.element[index_i+5][index_j+1]=fH.cy;
                Hessian.element[index_i+3][index_j+2]=fH.az; Hessian.element[index_i+4][index_j+2]=fH.bz;  Hessian.element[index_i+5][index_j+2]=fH.cz;
              }

              // adjust position-orientation
              if(Components[TypeMolB].Groups[jg].Rigid)
              {
                H.ax=Hessian.element[index_i+0][index_j+3]; H.bx=Hessian.element[index_i+1][index_j+3]; H.cx=Hessian.element[index_i+2][index_j+3];
                H.ay=Hessian.element[index_i+0][index_j+4]; H.by=Hessian.element[index_i+1][index_j+4]; H.cy=Hessian.element[index_i+2][index_j+4];
                H.az=Hessian.element[index_i+0][index_j+5]; H.bz=Hessian.element[index_i+1][index_j+5]; H.cz=Hessian.element[index_i+2][index_j+5];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                TransposeMatrix3x3(&H);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);
                TransposeMatrix3x3(&fH);

                Hessian.element[index_i+0][index_j+3]=fH.ax; Hessian.element[index_i+1][index_j+3]=fH.bx;  Hessian.element[index_i+2][index_j+3]=fH.cx;
                Hessian.element[index_i+0][index_j+4]=fH.ay; Hessian.element[index_i+1][index_j+4]=fH.by;  Hessian.element[index_i+2][index_j+4]=fH.cy;
                Hessian.element[index_i+0][index_j+5]=fH.az; Hessian.element[index_i+1][index_j+5]=fH.bz;  Hessian.element[index_i+2][index_j+5]=fH.cz;
              }

              // note: no need to adjust orientation-orientation (independent of Cartesian or fractional part)
              if(Components[TypeMolB].Groups[jg].Rigid) break;
            }
          }
        }

        // a rigid unit has only one index into the Hessian
        if(Components[TypeMolA].Groups[ig].Rigid) break;
      }
    }
  }


  for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
  {
    TypeMolA=Cations[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(Components[TypeMolA].Groups[ig].Rigid)
          index_i=Cations[CurrentSystem][I].Groups[ig].HessianIndex;
        else
          index_i=Cations[CurrentSystem][I].Atoms[i].HessianIndex;

        for(J=I;J<NumberOfCationMolecules[CurrentSystem];J++)
        {
          TypeMolB=Cations[CurrentSystem][J].Type;
          for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
          {
            for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
            {
              j=Components[TypeMolB].Groups[jg].Atoms[ja];

              if(Components[TypeMolB].Groups[jg].Rigid)
                index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
              else
                index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;

              // adjust position-position
              if((index_i>=0)&&(index_j>=0))
              {
                H.ax=Hessian.element[index_i][index_j];   H.bx=Hessian.element[index_i+1][index_j];   H.cx=Hessian.element[index_i+2][index_j];
                H.ay=Hessian.element[index_i][index_j+1]; H.by=Hessian.element[index_i+1][index_j+1]; H.cy=Hessian.element[index_i+2][index_j+1];
                H.az=Hessian.element[index_i][index_j+2]; H.bz=Hessian.element[index_i+1][index_j+2]; H.cz=Hessian.element[index_i+2][index_j+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);
                TransposeMatrix3x3(&fH);
                fH=MatrixMatrixMultiplication3x3(tempBox,fH);
                TransposeMatrix3x3(&fH);

                Hessian.element[index_i][index_j]=fH.ax;   Hessian.element[index_i+1][index_j]=fH.bx;    Hessian.element[index_i+2][index_j]=fH.cx;
                Hessian.element[index_i][index_j+1]=fH.ay; Hessian.element[index_i+1][index_j+1]=fH.by;  Hessian.element[index_i+2][index_j+1]=fH.cy;
                Hessian.element[index_i][index_j+2]=fH.az; Hessian.element[index_i+1][index_j+2]=fH.bz;  Hessian.element[index_i+2][index_j+2]=fH.cz;
              }

              // adjust orientation-position
              if(Components[TypeMolA].Groups[ig].Rigid)
              {
                H.ax=Hessian.element[index_i+3][index_j+0]; H.bx=Hessian.element[index_i+4][index_j+0]; H.cx=Hessian.element[index_i+5][index_j+0];
                H.ay=Hessian.element[index_i+3][index_j+1]; H.by=Hessian.element[index_i+4][index_j+1]; H.cy=Hessian.element[index_i+5][index_j+1];
                H.az=Hessian.element[index_i+3][index_j+2]; H.bz=Hessian.element[index_i+4][index_j+2]; H.cz=Hessian.element[index_i+5][index_j+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+3][index_j+0]=fH.ax; Hessian.element[index_i+4][index_j+0]=fH.bx;  Hessian.element[index_i+5][index_j+0]=fH.cx;
                Hessian.element[index_i+3][index_j+1]=fH.ay; Hessian.element[index_i+4][index_j+1]=fH.by;  Hessian.element[index_i+5][index_j+1]=fH.cy;
                Hessian.element[index_i+3][index_j+2]=fH.az; Hessian.element[index_i+4][index_j+2]=fH.bz;  Hessian.element[index_i+5][index_j+2]=fH.cz;
              }

              // adjust position-orientation
              if(Components[TypeMolB].Groups[jg].Rigid)
              {
                H.ax=Hessian.element[index_i+0][index_j+3]; H.bx=Hessian.element[index_i+1][index_j+3]; H.cx=Hessian.element[index_i+2][index_j+3];
                H.ay=Hessian.element[index_i+0][index_j+4]; H.by=Hessian.element[index_i+1][index_j+4]; H.cy=Hessian.element[index_i+2][index_j+4];
                H.az=Hessian.element[index_i+0][index_j+5]; H.bz=Hessian.element[index_i+1][index_j+5]; H.cz=Hessian.element[index_i+2][index_j+5];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                TransposeMatrix3x3(&H);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);
                TransposeMatrix3x3(&fH);

                Hessian.element[index_i+0][index_j+3]=fH.ax; Hessian.element[index_i+1][index_j+3]=fH.bx;  Hessian.element[index_i+2][index_j+3]=fH.cx;
                Hessian.element[index_i+0][index_j+4]=fH.ay; Hessian.element[index_i+1][index_j+4]=fH.by;  Hessian.element[index_i+2][index_j+4]=fH.cy;
                Hessian.element[index_i+0][index_j+5]=fH.az; Hessian.element[index_i+1][index_j+5]=fH.bz;  Hessian.element[index_i+2][index_j+5]=fH.cz;
              }

              // note: no need to adjust orientation-orientation (independent of Cartesian or fractional part)

              // a rigid unit has only one index into the Hessian
              if(Components[TypeMolB].Groups[jg].Rigid) break;
            }
          }
        }
        // a rigid unit has only one index into the Hessian
        if(Components[TypeMolA].Groups[ig].Rigid) break;
      }
    }
  }




  // correct strain cross-derivative
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
            {
              index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

              if(index_i>=0)
              {
                H.ax=Hessian.element[index_i+0][np+0]; H.bx=Hessian.element[index_i+0][np+1]; H.cx=Hessian.element[index_i+0][np+2];
                H.ay=Hessian.element[index_i+1][np+0]; H.by=Hessian.element[index_i+1][np+1]; H.cy=Hessian.element[index_i+1][np+2];
                H.az=Hessian.element[index_i+2][np+0]; H.bz=Hessian.element[index_i+2][np+1]; H.cz=Hessian.element[index_i+2][np+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+0][np+0]=fH.ax; Hessian.element[index_i+0][np+1]=fH.bx;  Hessian.element[index_i+0][np+2]=fH.cx;
                Hessian.element[index_i+1][np+0]=fH.ay; Hessian.element[index_i+1][np+1]=fH.by;  Hessian.element[index_i+1][np+2]=fH.cy;
                Hessian.element[index_i+2][np+0]=fH.az; Hessian.element[index_i+2][np+1]=fH.bz;  Hessian.element[index_i+2][np+2]=fH.cz;

                H.ax=Hessian.element[index_i+0][np+1]; H.bx=Hessian.element[index_i+0][np+3]; H.cx=Hessian.element[index_i+0][np+4];
                H.ay=Hessian.element[index_i+1][np+1]; H.by=Hessian.element[index_i+1][np+3]; H.cy=Hessian.element[index_i+1][np+4];
                H.az=Hessian.element[index_i+2][np+1]; H.bz=Hessian.element[index_i+2][np+3]; H.cz=Hessian.element[index_i+2][np+4];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+0][np+3]=fH.bx;  Hessian.element[index_i+0][np+4]=fH.cx;
                Hessian.element[index_i+1][np+3]=fH.by;  Hessian.element[index_i+1][np+4]=fH.cy;
                Hessian.element[index_i+2][np+3]=fH.bz;  Hessian.element[index_i+2][np+4]=fH.cz;

                H.ax=Hessian.element[index_i+0][np+2]; H.bx=Hessian.element[index_i+0][np+4]; H.cx=Hessian.element[index_i+0][np+5];
                H.ay=Hessian.element[index_i+1][np+2]; H.by=Hessian.element[index_i+1][np+4]; H.cy=Hessian.element[index_i+1][np+5];
                H.az=Hessian.element[index_i+2][np+2]; H.bz=Hessian.element[index_i+2][np+4]; H.cz=Hessian.element[index_i+2][np+5];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+0][np+5]=fH.cx;
                Hessian.element[index_i+1][np+5]=fH.cy;
                Hessian.element[index_i+2][np+5]=fH.cz;
              }
            }
          }

          for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
          {
            TypeMolA=Adsorbates[CurrentSystem][I].Type;
            for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
            {
              index_i=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndex;

              if(index_i>=0)
              {
                H.ax=Hessian.element[index_i+0][np+0]; H.bx=Hessian.element[index_i+0][np+1]; H.cx=Hessian.element[index_i+0][np+2];
                H.ay=Hessian.element[index_i+1][np+0]; H.by=Hessian.element[index_i+1][np+1]; H.cy=Hessian.element[index_i+1][np+2];
                H.az=Hessian.element[index_i+2][np+0]; H.bz=Hessian.element[index_i+2][np+1]; H.cz=Hessian.element[index_i+2][np+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+0][np+0]=fH.ax; Hessian.element[index_i+0][np+1]=fH.bx;  Hessian.element[index_i+0][np+2]=fH.cx;
                Hessian.element[index_i+1][np+0]=fH.ay; Hessian.element[index_i+1][np+1]=fH.by;  Hessian.element[index_i+1][np+2]=fH.cy;
                Hessian.element[index_i+2][np+0]=fH.az; Hessian.element[index_i+2][np+1]=fH.bz;  Hessian.element[index_i+2][np+2]=fH.cz;

                H.ax=Hessian.element[index_i+0][np+1]; H.bx=Hessian.element[index_i+0][np+3]; H.cx=Hessian.element[index_i+0][np+4];
                H.ay=Hessian.element[index_i+1][np+1]; H.by=Hessian.element[index_i+1][np+3]; H.cy=Hessian.element[index_i+1][np+4];
                H.az=Hessian.element[index_i+2][np+1]; H.bz=Hessian.element[index_i+2][np+3]; H.cz=Hessian.element[index_i+2][np+4];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+0][np+3]=fH.bx;  Hessian.element[index_i+0][np+4]=fH.cx;
                Hessian.element[index_i+1][np+3]=fH.by;  Hessian.element[index_i+1][np+4]=fH.cy;
                Hessian.element[index_i+2][np+3]=fH.bz;  Hessian.element[index_i+2][np+4]=fH.cz;

                H.ax=Hessian.element[index_i+0][np+2]; H.bx=Hessian.element[index_i+0][np+4]; H.cx=Hessian.element[index_i+0][np+5];
                H.ay=Hessian.element[index_i+1][np+2]; H.by=Hessian.element[index_i+1][np+4]; H.cy=Hessian.element[index_i+1][np+5];
                H.az=Hessian.element[index_i+2][np+2]; H.bz=Hessian.element[index_i+2][np+4]; H.cz=Hessian.element[index_i+2][np+5];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+0][np+5]=fH.cx;
                Hessian.element[index_i+1][np+5]=fH.cy;
                Hessian.element[index_i+2][np+5]=fH.cz;
              }
            }
          }

          for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
          {
            TypeMolA=Cations[CurrentSystem][I].Type;
            for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
            {
              index_i=Cations[CurrentSystem][I].Groups[ig].HessianIndex;

              if(index_i>=0)
              {
                H.ax=Hessian.element[index_i+0][np+0]; H.bx=Hessian.element[index_i+0][np+1]; H.cx=Hessian.element[index_i+0][np+2];
                H.ay=Hessian.element[index_i+1][np+0]; H.by=Hessian.element[index_i+1][np+1]; H.cy=Hessian.element[index_i+1][np+2];
                H.az=Hessian.element[index_i+2][np+0]; H.bz=Hessian.element[index_i+2][np+1]; H.cz=Hessian.element[index_i+2][np+2];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+0][np+0]=fH.ax; Hessian.element[index_i+0][np+1]=fH.bx;  Hessian.element[index_i+0][np+2]=fH.cx;
                Hessian.element[index_i+1][np+0]=fH.ay; Hessian.element[index_i+1][np+1]=fH.by;  Hessian.element[index_i+1][np+2]=fH.cy;
                Hessian.element[index_i+2][np+0]=fH.az; Hessian.element[index_i+2][np+1]=fH.bz;  Hessian.element[index_i+2][np+2]=fH.cz;

                H.ax=Hessian.element[index_i+0][np+1]; H.bx=Hessian.element[index_i+0][np+3]; H.cx=Hessian.element[index_i+0][np+4];
                H.ay=Hessian.element[index_i+1][np+1]; H.by=Hessian.element[index_i+1][np+3]; H.cy=Hessian.element[index_i+1][np+4];
                H.az=Hessian.element[index_i+2][np+1]; H.bz=Hessian.element[index_i+2][np+3]; H.cz=Hessian.element[index_i+2][np+4];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+0][np+3]=fH.bx;  Hessian.element[index_i+0][np+4]=fH.cx;
                Hessian.element[index_i+1][np+3]=fH.by;  Hessian.element[index_i+1][np+4]=fH.cy;
                Hessian.element[index_i+2][np+3]=fH.bz;  Hessian.element[index_i+2][np+4]=fH.cz;

                H.ax=Hessian.element[index_i+0][np+2]; H.bx=Hessian.element[index_i+0][np+4]; H.cx=Hessian.element[index_i+0][np+5];
                H.ay=Hessian.element[index_i+1][np+2]; H.by=Hessian.element[index_i+1][np+4]; H.cy=Hessian.element[index_i+1][np+5];
                H.az=Hessian.element[index_i+2][np+2]; H.bz=Hessian.element[index_i+2][np+4]; H.cz=Hessian.element[index_i+2][np+5];

                tempBox=Box[CurrentSystem];
                TransposeMatrix3x3(&tempBox);
                fH=MatrixMatrixMultiplication3x3(tempBox,H);

                Hessian.element[index_i+0][np+5]=fH.cx;
                Hessian.element[index_i+1][np+5]=fH.cy;
                Hessian.element[index_i+2][np+5]=fH.cz;
              }
            }
          }
          break;
        case MONOCLINIC:
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }


  // note: no need to correct strain-strain part

  for(i=0;i<np+nb;i++)
    for(j=0;j<np+nb;j++)
      Hessian.element[j][i]=Hessian.element[i][j];

*/

}


int FOLLOW_MODE;
int FOLLOW_VARIABLE;
int STEP_NUMBER;
int LOCATE_SADDLE;
int USE_NR_STEP;
int QNR;
int CURMOD;
int CURVAR;
int NEGREQ;
int NumberOfNegativeEigenValues;
int NumberOfZeroEigenValues;
REAL *MODVEC;

REAL CalculateBakerSaddlePointStepSize(int n,int NumberOfValidModes,int ISTEP,REAL *dx,REAL *GX,REAL *EigenValues,REAL_MATRIX HessianMatrix)
{
  int NumberOfVariables;
  REAL sum1;
  REAL scale;
  int i,j;
  REAL DXSIZE,Lambda,Lambda1,Lambda2;
  REAL lambdap,lambdan,flam,derivlam,sum2;

  NumberOfVariables=n;

  Lambda=0.0;
  Lambda1=0.0;
  Lambda2=0.0;

  // Initialize the step vector.
  for(i=0;i<NumberOfVariables;i++)
    dx[i]=0.0;

  lambdap=0.5*(EigenValues[0]+sqrt(SQR(EigenValues[0])+4.0*SQR(GX[0])));

  if(EigenValues[1]>=0.01)
    lambdan=0.0;
  else
    lambdan=0.5*(EigenValues[0]+EigenValues[1]);

  do
  {
    sum1=0.0;
    sum2=0.0;

    for(i=0;i<NumberOfVariables;i++)
    {
      sum1=sum1+SQR(GX[i])/(lambdan-EigenValues[i]);
      sum2=sum2+SQR(GX[i])/SQR(lambdan-EigenValues[i]);
    }

    flam=lambdan-sum1;
    derivlam=1.0+sum2;

    lambdan=lambdan-flam/derivlam;
  }while(flam>1e-8);

  Lambda=lambdan;

  scale=GX[0]/(lambdap-EigenValues[0]);
  if(fabs(lambdap-EigenValues[0])<1e-8) scale=0.1;
  for(i=0;i<NumberOfVariables;i++)
    dx[i]+=scale*HessianMatrix.element[0][i];

  // Calculate the step for the minimization.
  for(j=1;j<NumberOfValidModes;j++)
  {
    // skip the modes corresponding to no change in energy (zero eigenvalue)
    scale=GX[j]/(Lambda-EigenValues[j]);
    for(i=0;i<NumberOfVariables;i++)
      dx[i]+=scale*HessianMatrix.element[j][i];
  }


  // Calculate the step size.
  DXSIZE=0.0;
  for(i=0;i<NumberOfVariables;i++)
    DXSIZE+=SQR(dx[i]);
  DXSIZE=sqrt(DXSIZE);

  // Reduce the step size
  if(DXSIZE>MaximumStepLength)
    for(i=0;i<NumberOfVariables;i++)
      dx[i]*=(MaximumStepLength/DXSIZE);

  return Lambda;
}

void BakerSaddlePointSearch(int np,int nb,REAL *x,int run)
{
  int i,j,k;
  int NumberOfVariables;
  REAL Energy;
  REAL *dx,*Gradient,*GX,*EigenValues;
  REAL_MATRIX HessianMatrix;
  REAL RMSGradient,MaxGradient;
  REAL StepSize,Lambda;
  REAL_MATRIX3x3 StrainFirstDerivative;
  int NumberOfValidModes;

  NumberOfVariables=np+nb;

  // Check the input parameters for a minimum search.
  NEGREQ=1;

  CURMOD=0;

  // Check CURVAR.
  if(CURVAR>NumberOfVariables) CURVAR=-1;

  // If a variable is being followed set CURMOD to a value greater than 0.
  if(CURVAR>-1) CURMOD=1;

  // Check CURMOD.
  if((CURMOD<0)||(CURMOD>=NumberOfVariables)) CURMOD=0;

  QNR=FALSE;

  CurrentSystem=0;
  SamplePDBMovies(ALLOCATE,run);
  SamplePDBMovies(INITIALIZE,run);

  fprintf(OutputFilePtr[CurrentSystem],"Beginning saddle point search (Baker's method)\n");
  fprintf(OutputFilePtr[CurrentSystem],"----------------------------------------------\n");
  fflush(OutputFilePtr[CurrentSystem]);

  dx=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Gradient=(REAL*)calloc(NumberOfVariables,sizeof(REAL));  // CHECK
  GX=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  MODVEC=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  EigenValues=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  HessianMatrix=CreateRealMatrix(NumberOfVariables,NumberOfVariables);

  // Initialize the step vector.
  for(i=0;i<NumberOfVariables;i++)
    dx[i]=0.0;

  // Loop over the steps.
  for(k=0;k<MaximumNumberOfMinimizationSteps;k++)
  {
    // Form the new variable vector.
    for(i=0;i<NumberOfVariables;i++)
      x[i]+=dx[i];


    fprintf(OutputFilePtr[CurrentSystem],"Computing generalized Hessian matrix\n");
    fflush(OutputFilePtr[CurrentSystem]);
    ComputeDerivativesMinimization(np,nb,x,&Energy,Gradient,HessianMatrix,&StrainFirstDerivative);

    // project the constraints from the gradient and Hessian
    fprintf(OutputFilePtr[CurrentSystem],"Projecting constraints from generalized Hessian matrix\n");
    fflush(OutputFilePtr[CurrentSystem]);
    ProjectConstraintsFromHessianMatrix(np,nb,Gradient,HessianMatrix,TRUE,TRUE);

    if(MinimizationVariables==FRACTIONAL)
       ConvertHessianFromCartesianToFractional(np,nb,Gradient,HessianMatrix);

    // Diagonalize the Hessian matrix.
    fprintf(OutputFilePtr[CurrentSystem],"Computing eigenvalues and vectors\n");
    fflush(OutputFilePtr[CurrentSystem]);
    SolveEigenValuesAndVectorsHessian(HessianMatrix,EigenValues);

    // Determine the number of negative and zero eigenvalues.
    // remove eigenmodes with zero eigenvalues
    NumberOfValidModes=0;
    NumberOfNegativeEigenValues=0;
    NumberOfZeroEigenValues=0;
    for(i=0;i<NumberOfVariables;i++)
    {
      if(EigenValues[i]<-MINIMUM_EIGEN_VALUE) NumberOfNegativeEigenValues++;
      if(fabs(EigenValues[i])<MINIMUM_EIGEN_VALUE)
        NumberOfZeroEigenValues++;
      else
      {
        if(NumberOfZeroEigenValues>0)
        {
          EigenValues[NumberOfValidModes]=EigenValues[i];
          for(j=0;j<NumberOfVariables;j++)
            HessianMatrix.element[NumberOfValidModes][j]=HessianMatrix.element[i][j];
        }
        NumberOfValidModes++;
      }
    }

    // Transform the gradient vector to the eigenvector coordinate system
    for(j=0;j<NumberOfValidModes;j++)
    {
      GX[j]=0.0;
      for(i=0;i<NumberOfVariables;i++)
        GX[j]+=HessianMatrix.element[j][i]*Gradient[i];
    }

    // Calculate the step vector.
    Lambda=CalculateBakerSaddlePointStepSize(np+nb,NumberOfValidModes,k,dx,GX,EigenValues,HessianMatrix);

    if(k%PrintEvery==0)
      fprintf(OutputFilePtr[CurrentSystem],"Shifting parameter: %18.10f Lowest eigenvalue: %18.10f\n",Lambda,EigenValues[1]);

    StepSize=0.0;
    for(i=0;i<NumberOfVariables;i++)
      StepSize+=SQR(dx[i]);
    StepSize=sqrt(StepSize);

    // Find the RMS gradient.
    RMSGradient=0.0;
    MaxGradient=0.0;
    for(i=0;i<NumberOfVariables;i++)
    {
      if(fabs(Gradient[i])>MaxGradient) MaxGradient=fabs(Gradient[i]);
      RMSGradient+=SQR(Gradient[i]);
    }
    RMSGradient=sqrt(RMSGradient)/(REAL)NumberOfVariables;

    switch(Dimension)
    {
      case 2:
        fprintf(OutputFilePtr[CurrentSystem],"Iteration: %d Energy: %18.10f Area: %18.10f RMS gradient: %g  Max gradient: %g Number of negative eigenvalues: %d Number of zero eigenvalues: %d\n",
          k,(double)(Energy*ENERGY_TO_KELVIN),(double)Volume[CurrentSystem],(double)RMSGradient,(double)MaxGradient,NumberOfNegativeEigenValues,NumberOfZeroEigenValues);
        fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f     Strain derivative: %18.10f %18.10f\n",
            (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,
            (double)StrainFirstDerivative.ax,(double)StrainFirstDerivative.bx);
        fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f                        %18.10f %18.10f\n",
            (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,
            (double)StrainFirstDerivative.ay,(double)StrainFirstDerivative.by);
        fprintf(OutputFilePtr[CurrentSystem],"      Lengths: %18.10f %18.10f, Angle: %18.10f\n",
           (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,GammaAngle[CurrentSystem]*RAD2DEG);
        break;
      case 3:
        fprintf(OutputFilePtr[CurrentSystem],"Iteration: %d Energy: %18.10f Volume: %18.10f RMS gradient: %g  Max gradient: %g Number of negative eigenvalues: %d Number of zero eigenvalues: %d\n",
          k,(double)(Energy*ENERGY_TO_KELVIN),(double)Volume[CurrentSystem],(double)RMSGradient,(double)MaxGradient,NumberOfNegativeEigenValues,NumberOfZeroEigenValues);
        fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f %18.10f     Strain derivative: %18.10f %18.10f %18.10f\n",
            (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx,
            (double)StrainFirstDerivative.ax,(double)StrainFirstDerivative.bx,(double)StrainFirstDerivative.cx);
        fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
            (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy,
            (double)StrainFirstDerivative.ay,(double)StrainFirstDerivative.by,(double)StrainFirstDerivative.cy);
        fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
            (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz,
            (double)StrainFirstDerivative.az,(double)StrainFirstDerivative.bz,(double)StrainFirstDerivative.cz);
        fprintf(OutputFilePtr[CurrentSystem],"      Lengths: %18.10f %18.10f %18.10f, Angles: %18.10f %18.10f %18.10f\n",
           (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,(double)BoxProperties[CurrentSystem].az,
           AlphaAngle[CurrentSystem]*RAD2DEG,BetaAngle[CurrentSystem]*RAD2DEG,GammaAngle[CurrentSystem]*RAD2DEG);
        break;
    }
    fprintf(OutputFilePtr[CurrentSystem],"\n");
    fflush(OutputFilePtr[CurrentSystem]);

    // write the box and positions to a restart-file
    PrintRestartFile();

    // Check for convergence.
    if(((RMSGradient<RMSGradientTolerance)&&(MaxGradient<MaxGradientTolerance))&&(NumberOfNegativeEigenValues==NEGREQ))
    //if(((RMSGradient<RMSGradientTolerance)&&(MaxGradient<MaxGradientTolerance)))
    {
      fprintf(OutputFilePtr[CurrentSystem],"SUCCES: RMS Gradient tolerance %g reached (%g)\n",(double)RMSGradientTolerance,(double)RMSGradient);
      fprintf(OutputFilePtr[CurrentSystem],"        Max Gradient tolerance %g reached (%g)\n",(double)MaxGradientTolerance,(double)MaxGradient);
      break;
    }
    CurrentSystem=0;
    if(k%PrintEvery==0)
      SamplePDBMovies(SAMPLE,run);
  }
  SamplePDBMovies(FINALIZE,run);
}


REAL CalculateBakerStepSize(int n,int NumberOfValidModes,int ISTEP,REAL *dx,REAL *GX,REAL *EigenValues,REAL_MATRIX HessianMatrix)
{
  int NumberOfVariables;
  const int LambdaIterations=100000;
  REAL sum1;
  REAL scale;
  int i,j,k;
  REAL DXSIZE,Lambda,Lambda1,Lambda2;
  REAL fn,d1,d2,diff,rtrm,trm,d11,d22,step,stepmx;
  int LOWER;

  NumberOfVariables=n;

  Lambda=0.0;
  Lambda1=0.0;
  Lambda2=0.0;

  // Initialize the step vector.
  for(i=0;i<NumberOfVariables;i++)
    dx[i]=0.0;

  // Do a N-R step if the number of negative eigenvalues is OK and QNEWTON is on.
  //if(QNR&&(NumberOfNegativeEigenValues==0)||(1==1))
  if(QNR&&(NumberOfNegativeEigenValues==0))
  {
    // Calculate the Newton-Raphson step.
    for(j=0;j<NumberOfValidModes;j++)
    {
      scale=-GX[j]/EigenValues[j];
      for(i=0;i<NumberOfVariables;i++)
        dx[i]+=scale*HessianMatrix.element[j][i];
    }
  }
  else
  {
    // Find the step for the minimization.
    LOWER=0;

    switch(ComputeLambda)
    {
      case LAMBDA_METHOD_1:
        Lambda=EigenValues[LOWER]-1.0;

        for(k=0;k<LambdaIterations;k++)
        {
          fn=Lambda;
          d1=1.0;
          d2=0.0;
          for(j=0;j<NumberOfValidModes;j++)
          {
            trm=SQR(GX[j]);
            diff=Lambda-EigenValues[j];
            if (fabs(diff)<1e-8) diff=SIGN(1e-8,diff);
            rtrm=1.0/diff;
            trm=trm*rtrm;
            fn=fn-trm;
            trm=trm*rtrm;
            d1=d1+trm;
            trm=trm*rtrm;
            d2=d2-2.0*trm;
          }
          d11=2.0*d1*fn;
          d22=2.0*(d1*d1+fn*d2);
          fn=fn*fn;
          step=-d11/fabs(d22);

          if(step<0.0)
            stepmx=100.0;
          else
            stepmx=0.5*(EigenValues[LOWER]-Lambda);

          if(fabs(step)>stepmx)
            step=SIGN(stepmx,step);

          Lambda+=step;
        }
        break;
      case LAMBDA_METHOD_2:
        if(EigenValues[LOWER]<0.0)
        {
          Lambda1=EigenValues[LOWER]-MINIMUM_EIGEN_VALUE;
          Lambda=EigenValues[LOWER]-500.0;
          Lambda2=-1e6;
        }

        for(k=0;k<LambdaIterations;k++)
        {
          sum1=0.0;
          for(i=0;i<NumberOfValidModes;i++)
            sum1+=SQR(GX[i])/(Lambda-EigenValues[i]);

          if(EigenValues[0]>0.0)
            Lambda=sum1;
          else
          {
            if(sum1<Lambda) Lambda1=Lambda;
            if(sum1>Lambda) Lambda2=Lambda;
            if(Lambda2>-1e6)
              Lambda=0.5*(Lambda1+Lambda2);
            else if(Lambda2==-1e6)
              Lambda-=50.0;
          }
        }
        break;
      case LAMBDA_METHOD_3:
        Lambda=-1.0e6;
        do
        {
          sum1=0.0;
          for(i=0;i<NumberOfValidModes;i++)
          {
            diff=Lambda-EigenValues[i];
            if (fabs(diff)<1e-8) diff=SIGN(1e-8,diff);
              sum1+=SQR(GX[i])/diff;
          }
          diff=fabs(Lambda-sum1);
          Lambda=sum1;
        } while(diff>1e-8);
        break;
    }

    sum1=0.0;
    for(i=0;i<NumberOfValidModes;i++)
    {
      diff=Lambda-EigenValues[i];
      if (fabs(diff)<1e-8) diff=SIGN(1e-8,diff);
        sum1+=SQR(GX[i])/diff;
    }

    // check for the value of lambda
    if((Lambda>EigenValues[LOWER])||(EigenValues[LOWER]>0.0&&Lambda>0.0))
      fprintf(stderr, "Baker minimization: Error in lambda\n");

    //if((SimulationType==MINIMIZATION)&&(k%PrintEvery==0))
    //  fprintf(OutputFilePtr[CurrentSystem],"Shifting parameter: %18.10f Lowest eigenvalue: %18.10f\n",Lambda,EigenValues[LOWER]);

    // Modify the CURMOD eigenvalues and eigenvectors.
    if(CURMOD>0)
    {
      EigenValues[CURMOD-1]=Lambda-1.0;
      for(i=0;i<NumberOfVariables;i++)
        HessianMatrix.element[CURMOD-1][i]=0.0;
    }

    // Calculate the step for the minimization.
    for(j=0;j<NumberOfValidModes;j++)
    {
      diff=Lambda-EigenValues[j];
      if (fabs(diff)>MINIMUM_EIGEN_VALUE)
      {
        scale=GX[j]/diff;
        for(i=0;i<NumberOfVariables;i++)
          dx[i]+=scale*HessianMatrix.element[j][i];
      }
    }
  }

  // Calculate the step size.
  DXSIZE=0.0;
  for(i=0;i<NumberOfVariables;i++)
    DXSIZE+=SQR(dx[i]);
  DXSIZE=sqrt(DXSIZE);

  // Reduce the step size
  if(DXSIZE>MaximumStepLength)
    for(i=0;i<NumberOfVariables;i++)
      dx[i]*=(MaximumStepLength/DXSIZE);

  return Lambda;
}

void BakerMinimization(int np,int nb,REAL *x,int run)
{
  int i,j,k,l,f1,A,Type;
  int NumberOfVariables;
  REAL Energy;
  REAL *dx,*Gradient,*GX,*EigenValues;
  REAL_MATRIX HessianMatrix;
  REAL RMSGradient,MaxGradient;
  REAL StepSizePrevious,StepSize;
  REAL_MATRIX3x3 StrainFirstDerivative;
  int NumberOfValidModes;
  VECTOR posA,posB,posC,posD,posE;
  REAL Pressure,Lambda;

  NumberOfVariables=np+nb;

  QNR=FALSE;

  CurrentSystem=0;
  SamplePDBMovies(ALLOCATE,run);
  SamplePDBMovies(INITIALIZE,run);

  fprintf(OutputFilePtr[CurrentSystem],"Beginning Baker minimization:\n");
  fprintf(OutputFilePtr[CurrentSystem],"-----------------------------\n");
  fflush(OutputFilePtr[CurrentSystem]);

  dx=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Gradient=(REAL*)calloc(NumberOfVariables,sizeof(REAL));  // CHECK
  GX=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  MODVEC=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  EigenValues=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  HessianMatrix=CreateRealMatrix(NumberOfVariables,NumberOfVariables);

  StepSize=1.0/pow(Volume[CurrentSystem],1.0/3.0);
  StepSizePrevious=1.0/pow(Volume[CurrentSystem],1.0/3.0);

  // Check the input parameters for a minimum search.
  if (NEGREQ==0)
  {
    // Reset CURMOD and CURVAR.
    CURMOD=-1;
    CURVAR=-1;
  }
  // Check the input parameters for a saddle point search.
  else
  {
    // Check CURVAR.
    if(CURVAR>NumberOfVariables) CURVAR=-1;

    // If a variable is being followed set CURMOD to a value greater than 0.
    if(CURVAR>-1) CURMOD=1;

    // Check CURMOD.
    if((CURMOD<0)||(CURMOD>=NumberOfVariables)) CURMOD=0;
  }

  // Initialize the step vector.
  for(i=0;i<NumberOfVariables;i++)
    dx[i]=0.0;

  // Loop over the steps.
  for(k=0;k<MaximumNumberOfMinimizationSteps;k++)
  {
    // Form the new variable vector.
    for(i=0;i<NumberOfVariables;i++)
      x[i]+=dx[i];

    // Calculate the function and its derivatives.
    //ApplyStrainOnBox(x);
    //ComputeGeneralizedHessian(NumberOfPositionVariables,&Energy,x,Gradient,&StrainDerTensor,HessianMatrix);

    if(k%PrintEvery==0)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Computing generalized Hessian matrix\n");
      fflush(OutputFilePtr[CurrentSystem]);
    }
    ComputeDerivativesMinimization(np,nb,x,&Energy,Gradient,HessianMatrix,&StrainFirstDerivative);

    // project the constraints from the gradient and Hessian
    fprintf(OutputFilePtr[CurrentSystem],"Projecting constraints from generalized Hessian matrix\n");
    fflush(OutputFilePtr[CurrentSystem]);
    ProjectConstraintsFromHessianMatrix(np,nb,Gradient,HessianMatrix,TRUE,TRUE);

    if(MinimizationVariables==FRACTIONAL)
      ConvertHessianFromCartesianToFractional(np,nb,Gradient,HessianMatrix);

    // Diagonalize the Hessian matrix.
    if(k%PrintEvery==0)
    {
    fprintf(OutputFilePtr[CurrentSystem],"Computing eigenvalues and vectors\n");
    fflush(OutputFilePtr[CurrentSystem]);
    }
    SolveEigenValuesAndVectorsHessian(HessianMatrix,EigenValues);

    // Determine the number of negative and zero eigenvalues.
    // remove eigenmodes with zero eigenvalues
    NumberOfValidModes=0;
    NumberOfNegativeEigenValues=0;
    NumberOfZeroEigenValues=0;
    for(i=0;i<NumberOfVariables;i++)
    {
      if(EigenValues[i]<-MINIMUM_EIGEN_VALUE) NumberOfNegativeEigenValues++;
      if(fabs(EigenValues[i])<MINIMUM_EIGEN_VALUE)
        NumberOfZeroEigenValues++;
      else
      {
        if(NumberOfZeroEigenValues>0)
        {
          EigenValues[NumberOfValidModes]=EigenValues[i];
          for(j=0;j<NumberOfVariables;j++)
            HessianMatrix.element[NumberOfValidModes][j]=HessianMatrix.element[i][j];
        }
        NumberOfValidModes++;
      }
    }


    // Transform the gradient vector to the eigenvector coordinate system
    for(j=0;j<NumberOfValidModes;j++)
    {
      GX[j]=0.0;
      for(i=0;i<NumberOfVariables;i++)
        GX[j]+=HessianMatrix.element[j][i]*Gradient[i];
    }

    // Find the RMS gradient.
    RMSGradient=0.0;
    MaxGradient=0.0;
    for(i=0;i<NumberOfVariables;i++)
    {
      if(fabs(Gradient[i])>MaxGradient) MaxGradient=fabs(Gradient[i]);
      RMSGradient+=SQR(Gradient[i]);
    }

    RMSGradient=sqrt(RMSGradient)/(REAL)NumberOfVariables;
    MaximumStepLength=MIN2(MaximumStepLengthInput,pow(RMSGradient,MinimizationConvergenceFactor));

    // Calculate the step vector.
    Lambda=CalculateBakerStepSize(np+nb,NumberOfValidModes,k,dx,GX,EigenValues,HessianMatrix);

    if(k%PrintEvery==0)
      fprintf(OutputFilePtr[CurrentSystem],"Shifting parameter: %18.10f Lowest eigenvalue: %18.10f\n",Lambda,EigenValues[0]);

    StepSize=0.0;
    for(i=0;i<NumberOfVariables;i++)
      StepSize+=SQR(dx[i]);
    StepSize=sqrt(StepSize);

    // Find the RMS gradient.
/*
    RMSGradient=0.0;
    MaxGradient=0.0;
    for(i=0;i<NumberOfVariables;i++)
    {
      if(fabs(Gradient[i])>MaxGradient) MaxGradient=fabs(Gradient[i]);
      RMSGradient+=SQR(Gradient[i]);
    }
    RMSGradient=sqrt(RMSGradient)/(REAL)NumberOfVariables;
*/
    RMSGradient=0.0;
    MaxGradient=0.0;
    for(i=0;i<NumberOfVariables;i++)
    {
      if(fabs(Gradient[i])>MaxGradient) MaxGradient=fabs(Gradient[i]);
      RMSGradient+=SQR(Gradient[i]);
    }
/*
    if(fabs(Gradient[NumberOfCoordinatesMinimizationVariables]+1e4)>MaxGradient) MaxGradient=fabs(Gradient[NumberOfCoordinatesMinimizationVariables]+1e4);
    RMSGradient+=SQR(Gradient[NumberOfCoordinatesMinimizationVariables]+1e4);
    if(fabs(Gradient[NumberOfCoordinatesMinimizationVariables+1]+1e4)>MaxGradient) MaxGradient=fabs(Gradient[NumberOfCoordinatesMinimizationVariables+1]+1e4);
    RMSGradient+=SQR(Gradient[NumberOfCoordinatesMinimizationVariables+1]+1e4);
    if(fabs(Gradient[NumberOfCoordinatesMinimizationVariables+2]+1e4)>MaxGradient) MaxGradient=fabs(Gradient[NumberOfCoordinatesMinimizationVariables+2]+1e4);
    RMSGradient+=SQR(Gradient[NumberOfCoordinatesMinimizationVariables+2]+1e4);
*/

    RMSGradient=sqrt(RMSGradient)/(REAL)NumberOfVariables;


    if(k%PrintEvery==0)
    {
      switch(Dimension)
      {
        case 2:
          fprintf(OutputFilePtr[CurrentSystem],"Iteration: %d Energy: %18.10f Area: %18.10f RMS gradient: %g  Max gradient: %g Number of negative eigenvalues: %d Number of zero eigenvalues: %d\n",
            k,(double)(Energy*ENERGY_TO_KELVIN),(double)Volume[CurrentSystem],(double)RMSGradient,(double)MaxGradient,NumberOfNegativeEigenValues,NumberOfZeroEigenValues);
          fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f     Strain derivative: %18.10f %18.10f\n",
              (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,
              (double)StrainFirstDerivative.ax,(double)StrainFirstDerivative.bx);
          fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f                        %18.10f %18.10f\n",
              (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,
              (double)StrainFirstDerivative.ay,(double)StrainFirstDerivative.by);
          fprintf(OutputFilePtr[CurrentSystem],"      Lengths: %18.10f %18.10f, Angle: %18.10f\n",
             (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,GammaAngle[CurrentSystem]*RAD2DEG);
          break;
        case 3:
          fprintf(OutputFilePtr[CurrentSystem],"Iteration: %d Energy: %18.10f Volume: %18.10f RMS gradient: %g  Max gradient: %g Number of negative eigenvalues: %d Number of zero eigenvalues: %d\n",
            k,(double)(Energy*ENERGY_TO_KELVIN),(double)Volume[CurrentSystem],(double)RMSGradient,(double)MaxGradient,NumberOfNegativeEigenValues,NumberOfZeroEigenValues);
          fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f %18.10f     Strain derivative: %18.10f %18.10f %18.10f\n",
              (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx,
              (double)StrainFirstDerivative.ax,(double)StrainFirstDerivative.bx,(double)StrainFirstDerivative.cx);
          fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
              (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy,
              (double)StrainFirstDerivative.ay,(double)StrainFirstDerivative.by,(double)StrainFirstDerivative.cy);
          fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
              (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz,
              (double)StrainFirstDerivative.az,(double)StrainFirstDerivative.bz,(double)StrainFirstDerivative.cz);
          fprintf(OutputFilePtr[CurrentSystem],"      Lengths: %18.10f %18.10f %18.10f, Angles: %18.10f %18.10f %18.10f\n",
             (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,(double)BoxProperties[CurrentSystem].az,
             AlphaAngle[CurrentSystem]*RAD2DEG,BetaAngle[CurrentSystem]*RAD2DEG,GammaAngle[CurrentSystem]*RAD2DEG);
          break;
      }
    }
    for(i=0;i<NumberOfComponents;i++)
    {
      if(Components[i].RestrictEnantionface)
      {
        posA=Components[i].EnantiofaceAtoms[0]->Position;
        posB=Components[i].EnantiofaceAtoms[1]->Position;
        posC=Components[i].EnantiofaceAtoms[2]->Position;
        posD=Components[i].EnantiofaceAtoms[3]->Position;
        posE=Components[i].EnantiofaceAtoms[4]->Position;
        if(k%PrintEvery==0)
          fprintf(OutputFilePtr[CurrentSystem],"Component: %d, enantionface should be: %s, current enantioface: %s\n",i,
            Components[i].Enantioface==ENANTIOFACE_RE?"Re":"Si",
            CheckEnantioFace(posA,posB,posC,posD,posE)==ENANTIOFACE_RE?"Re":"Si");
      }
    }

    if(k%PrintEvery==0)
    fprintf(OutputFilePtr[CurrentSystem],"\n");

    fflush(OutputFilePtr[CurrentSystem]);

    // write the box and positions to a restart-file
    PrintRestartFile();

    // Check for convergence.
    if(((RMSGradient<RMSGradientTolerance)&&(MaxGradient<MaxGradientTolerance))&&(NumberOfNegativeEigenValues==NEGREQ))
    //if(((RMSGradient<RMSGradientTolerance)&&(MaxGradient<MaxGradientTolerance)))
    {
      fprintf(OutputFilePtr[CurrentSystem],"SUCCES: RMS Gradient tolerance %g reached (%g)\n",(double)RMSGradientTolerance,(double)RMSGradient);
      fprintf(OutputFilePtr[CurrentSystem],"        Max Gradient tolerance %g reached (%g)\n",(double)MaxGradientTolerance,(double)MaxGradient);
      break;
    }
    CurrentSystem=0;
    if(k%PrintEvery==0)
      SamplePDBMovies(SAMPLE,run);
  }
  SamplePDBMovies(FINALIZE,run);

  // if upper-triangle form, convert the cell (and atom positions) to upper-triangle
  //if(Ensemble[CurrentSystem]==NPTPR&&(NPTPRCellType[CurrentSystem]==REGULAR_UPPER_TRIANGLE||NPTPRCellType[CurrentSystem]==MONOCLINIC_UPPER_TRIANGLE))
  //  AdjustToUpperTriangle();

  ComputeDerivativesMinimization(np,nb,x,&Energy,Gradient,HessianMatrix,&StrainFirstDerivative);

  // project the constraints from the gradient and Hessian
  ProjectConstraintsFromHessianMatrix(np,nb,Gradient,HessianMatrix,TRUE,TRUE);

  if(MinimizationVariables==FRACTIONAL)
    ConvertHessianFromCartesianToFractional(np,nb,Gradient,HessianMatrix);

  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
  fprintf(OutputFilePtr[CurrentSystem],"Variable vector:\n");
  fprintf(OutputFilePtr[CurrentSystem],"================\n");
  for(i=0;i<NumberOfVariables;i++)
    fprintf(OutputFilePtr[CurrentSystem],"x: %d = %g\n",i,(double)x[i]);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
  fprintf(OutputFilePtr[CurrentSystem],"Gradient vector:\n");
  fprintf(OutputFilePtr[CurrentSystem],"================\n");
  for(i=0;i<NumberOfVariables;i++)
    fprintf(OutputFilePtr[CurrentSystem],"Gradient: %d = %g\n",i,(double)Gradient[i]);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  SolveEigenValuesAndVectorsHessian(HessianMatrix,EigenValues);

  fprintf(OutputFilePtr[CurrentSystem],"Eigenvalues:\n");
  fprintf(OutputFilePtr[CurrentSystem],"============\n");
  for(i=0;i<NumberOfVariables;i++)
    fprintf(OutputFilePtr[CurrentSystem],"Eigenvalue: %d = % 18.10f\n",i,(double)EigenValues[i]);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  fprintf(OutputFilePtr[CurrentSystem],"Final energy after minimization: %18.10f in %d steps\n",(double)(Energy*ENERGY_TO_KELVIN),k);

  Pressure=-(StrainFirstDerivative.ax+StrainFirstDerivative.by+StrainFirstDerivative.cz)/(Dimension*Volume[CurrentSystem]);
  fprintf(OutputFilePtr[CurrentSystem],"Final pressure after minimization: %18.10f [Pa]\n",(double)(Pressure*PRESSURE_CONVERSION_FACTOR));
  switch(Dimension)
  {
    case 3:
     fprintf(OutputFilePtr[CurrentSystem],"Final stress after minimization: {{%f,%f,%f},{-,%f,%f},{-,-,%f}} %s\n",
          (double)(-StrainFirstDerivative.ax*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.bx*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.cx*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.by*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.bz*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.cz*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          UseReducedUnits?"[-]":"[Pa]");
      break;
    case 2:
     fprintf(OutputFilePtr[CurrentSystem],"Final stress after minimization: {{%f,%f},{-,%f}} %s\n",
          (double)(-StrainFirstDerivative.ax*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.bx*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.by*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          UseReducedUnits?"[-]":"[Pa]");
      break;
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  fprintf(OutputFilePtr[CurrentSystem],"Volume [A^3]: %18.10f\n",Volume[CurrentSystem]);
  fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f %18.10f     Strain derivative: %18.10f %18.10f %18.10f\n",
        (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx,
        (double)StrainFirstDerivative.ax,(double)StrainFirstDerivative.bx,(double)StrainFirstDerivative.cx);
  fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
        (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy,
        (double)StrainFirstDerivative.ay,(double)StrainFirstDerivative.by,(double)StrainFirstDerivative.cy);
  fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
        (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz,
        (double)StrainFirstDerivative.az,(double)StrainFirstDerivative.bz,(double)StrainFirstDerivative.cz);
  fprintf(OutputFilePtr[CurrentSystem],"      Final Lengths: %18.10f %18.10f %18.10f, Angles: %18.10f %18.10f %18.10f\n",
       (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,(double)BoxProperties[CurrentSystem].az,
       AlphaAngle[CurrentSystem]*RAD2DEG,BetaAngle[CurrentSystem]*RAD2DEG,GammaAngle[CurrentSystem]*RAD2DEG);
  fprintf(OutputFilePtr[CurrentSystem],"\n");

  ComputeBornTerm=FALSE;
  CalculateForce();

  ShakeInMinimization();

  fprintf(OutputFilePtr[CurrentSystem],"Strain derivative from force routines (NOTE: different when fixed atoms are present)\n");
  fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %18.10f\n",StrainDerivativeTensor[CurrentSystem].ax,
          StrainDerivativeTensor[CurrentSystem].bx,StrainDerivativeTensor[CurrentSystem].cx);
  fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %18.10f\n",StrainDerivativeTensor[CurrentSystem].ay,
          StrainDerivativeTensor[CurrentSystem].by,StrainDerivativeTensor[CurrentSystem].cy);
  fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %18.10f\n",StrainDerivativeTensor[CurrentSystem].az,
          StrainDerivativeTensor[CurrentSystem].bz,StrainDerivativeTensor[CurrentSystem].cz);
  fprintf(OutputFilePtr[CurrentSystem],"\n");

  fprintf(OutputFilePtr[CurrentSystem],"Forces\n");
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Framework[%d] Atom: %d  %g %g %g\n",f1,i,
        Framework[CurrentSystem].Atoms[f1][i].Fixed.x?0.0:Framework[CurrentSystem].Atoms[f1][i].Force.x,
        Framework[CurrentSystem].Atoms[f1][i].Fixed.y?0.0:Framework[CurrentSystem].Atoms[f1][i].Force.y,
        Framework[CurrentSystem].Atoms[f1][i].Fixed.z?0.0:Framework[CurrentSystem].Atoms[f1][i].Force.z);
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        fprintf(OutputFilePtr[CurrentSystem],"Adsorbate[%d] Group: %d  Com force: %g %g %g\n",i,l,
               Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass.x?0.0:Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.x,
               Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass.y?0.0:Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.y,
               Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass.z?0.0:Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.z);
        fprintf(OutputFilePtr[CurrentSystem],"Adsorbate[%d] Group: %d  torque: %g %g %g\n",i,l,
               Adsorbates[CurrentSystem][i].Groups[l].FixedOrientation.x?0.0:Adsorbates[CurrentSystem][i].Groups[l].Torque.x,
               Adsorbates[CurrentSystem][i].Groups[l].FixedOrientation.y?0.0:Adsorbates[CurrentSystem][i].Groups[l].Torque.y,
               Adsorbates[CurrentSystem][i].Groups[l].FixedOrientation.z?0.0:Adsorbates[CurrentSystem][i].Groups[l].Torque.z);
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          fprintf(OutputFilePtr[CurrentSystem],"Adsorbate[%d] Atom: %d  %g %g %g\n",i,A,
            Adsorbates[CurrentSystem][i].Atoms[A].Fixed.x?0.0:Adsorbates[CurrentSystem][i].Atoms[A].Force.x,
            Adsorbates[CurrentSystem][i].Atoms[A].Fixed.y?0.0:Adsorbates[CurrentSystem][i].Atoms[A].Force.y,
            Adsorbates[CurrentSystem][i].Atoms[A].Fixed.z?0.0:Adsorbates[CurrentSystem][i].Atoms[A].Force.z);
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        fprintf(OutputFilePtr[CurrentSystem],"Cation[%d] Group: %d  Com force: %g %g %g\n",i,l,
               Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass.x?0.0:Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.x,
               Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass.y?0.0:Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.y,
               Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass.z?0.0:Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.z);
        fprintf(OutputFilePtr[CurrentSystem],"Cation[%d] Group: %d  torque: %g %g %g\n",i,l,
               Cations[CurrentSystem][i].Groups[l].FixedOrientation.x?0.0:Cations[CurrentSystem][i].Groups[l].Torque.x,
               Cations[CurrentSystem][i].Groups[l].FixedOrientation.y?0.0:Cations[CurrentSystem][i].Groups[l].Torque.y,
               Cations[CurrentSystem][i].Groups[l].FixedOrientation.z?0.0:Cations[CurrentSystem][i].Groups[l].Torque.z);
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          fprintf(OutputFilePtr[CurrentSystem],"Cation[%d] Atom: %d  %g %g %g\n",i,A,
            Cations[CurrentSystem][i].Atoms[A].Fixed.x?0.0:Cations[CurrentSystem][i].Atoms[A].Force.x,
            Cations[CurrentSystem][i].Atoms[A].Fixed.y?0.0:Cations[CurrentSystem][i].Atoms[A].Force.y,
            Cations[CurrentSystem][i].Atoms[A].Fixed.z?0.0:Cations[CurrentSystem][i].Atoms[A].Force.z);
        }
      }
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  free(dx);
  free(Gradient);
  free(GX);
  free(MODVEC);
  free(EigenValues);
  DeleteRealMatrix(HessianMatrix);

  //ComputeElasticConstantsAfterMinimizationFromEnergyVSUnitCell(np,nb,x);
  ComputeElasticConstantsAfterMinimization(np,nb,x);
}


void NewtonRaphsonMinimization(int np,int nb,REAL *x,int run)
{
  int i,j,k,l,f1,A,Type;
  int NumberOfVariables;
  REAL Energy,scale;
  REAL *dx,*Gradient,*GX,*EigenValues;
  REAL_MATRIX HessianMatrix;
  REAL RMSGradient,MaxGradient;
  REAL StepSizePrevious,StepSize;
  REAL_MATRIX3x3 StrainFirstDerivative;
  int NumberOfValidModes;
  VECTOR posA,posB,posC,posD,posE;

  NumberOfVariables=np+nb;

  QNR=FALSE;

  CurrentSystem=0;
  SamplePDBMovies(ALLOCATE,run);
  SamplePDBMovies(INITIALIZE,run);

  fprintf(OutputFilePtr[CurrentSystem],"Beginning Newton-Raphson minimization:\n");
  fprintf(OutputFilePtr[CurrentSystem],"--------------------------------------\n");
  fflush(OutputFilePtr[CurrentSystem]);

  dx=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Gradient=(REAL*)calloc(NumberOfVariables,sizeof(REAL));  // CHECK
  GX=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  MODVEC=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  EigenValues=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  HessianMatrix=CreateRealMatrix(NumberOfVariables,NumberOfVariables);

  StepSize=1.0/pow(Volume[CurrentSystem],1.0/3.0);
  StepSizePrevious=1.0/pow(Volume[CurrentSystem],1.0/3.0);

  // Check the input parameters for a minimum search.
  if (NEGREQ==0)
  {
    // Reset CURMOD and CURVAR.
    CURMOD=-1;
    CURVAR=-1;
  }
  // Check the input parameters for a saddle point search.
  else
  {
    // Check CURVAR.
    if(CURVAR>NumberOfVariables) CURVAR=-1;

    // If a variable is being followed set CURMOD to a value greater than 0.
    if(CURVAR>-1) CURMOD=1;

    // Check CURMOD.
    if((CURMOD<0)||(CURMOD>=NumberOfVariables)) CURMOD=0;
  }

  // Initialize the step vector.
  for(i=0;i<NumberOfVariables;i++)
    dx[i]=0.0;

  // Loop over the steps.
  for(k=0;k<MaximumNumberOfMinimizationSteps;k++)
  {
    // Form the new variable vector.
    for(i=0;i<NumberOfVariables;i++)
      x[i]+=dx[i];

    // Calculate the function and its derivatives.
    //ApplyStrainOnBox(x);
    //ComputeGeneralizedHessian(NumberOfPositionVariables,&Energy,x,Gradient,&StrainDerTensor,HessianMatrix);


    fprintf(OutputFilePtr[CurrentSystem],"Computing generalized Hessian matrix\n");
    fflush(OutputFilePtr[CurrentSystem]);
    ComputeDerivativesMinimization(np,nb,x,&Energy,Gradient,HessianMatrix,&StrainFirstDerivative);

    // project the constraints from the gradient and Hessian
    fprintf(OutputFilePtr[CurrentSystem],"Projecting constraints from generalized Hessian matrix\n");
    fflush(OutputFilePtr[CurrentSystem]);
    ProjectConstraintsFromHessianMatrix(np,nb,Gradient,HessianMatrix,TRUE,TRUE);

    if(MinimizationVariables==FRACTIONAL)
      ConvertHessianFromCartesianToFractional(np,nb,Gradient,HessianMatrix);

    //for(i=0;i<NumberOfPositionVariables;i++)
    //  fprintf(stderr, "%d gradient %f\n",i,Gradient[i]);

    // Diagonalize the Hessian matrix.
    fprintf(OutputFilePtr[CurrentSystem],"Computing eigenvalues and vectors\n");
    fflush(OutputFilePtr[CurrentSystem]);
    SolveEigenValuesAndVectorsHessian(HessianMatrix,EigenValues);

    // Determine the number of negative and zero eigenvalues.
    // remove eigenmodes with zero eigenvalues
    NumberOfValidModes=0;
    NumberOfNegativeEigenValues=0;
    NumberOfZeroEigenValues=0;
    for(i=0;i<NumberOfVariables;i++)
    {
      if(EigenValues[i]<-MINIMUM_EIGEN_VALUE) NumberOfNegativeEigenValues++;
      if(fabs(EigenValues[i])<MINIMUM_EIGEN_VALUE)
        NumberOfZeroEigenValues++;
      else
      {
        if(NumberOfZeroEigenValues>0)
        {
          EigenValues[NumberOfValidModes]=EigenValues[i];
          for(j=0;j<NumberOfVariables;j++)
            HessianMatrix.element[NumberOfValidModes][j]=HessianMatrix.element[i][j];
        }
        NumberOfValidModes++;
      }
    }


    // Transform the gradient vector to the eigenvector coordinate system
    for(j=0;j<NumberOfValidModes;j++)
    {
      GX[j]=0.0;
      for(i=0;i<NumberOfVariables;i++)
        GX[j]+=HessianMatrix.element[j][i]*Gradient[i];
    }

    // Find the RMS gradient.
    RMSGradient=0.0;
    MaxGradient=0.0;
    for(i=0;i<NumberOfVariables;i++)
    {
      if(fabs(Gradient[i])>MaxGradient) MaxGradient=fabs(Gradient[i]);
      RMSGradient+=SQR(Gradient[i]);
    }
    RMSGradient=sqrt(RMSGradient)/(REAL)NumberOfVariables;
    MaximumStepLength=MIN2(MaximumStepLengthInput,pow(RMSGradient,MinimizationConvergenceFactor));

    // Calculate the Newton-Raphson step.
    for(j=0;j<NumberOfValidModes;j++)
    {
      scale=-GX[j]/EigenValues[j];
      for(i=0;i<NumberOfVariables;i++)
        dx[i]+=scale*HessianMatrix.element[j][i];
    }

    // Calculate the step size.
    scale=0.0;
    for(i=0;i<NumberOfVariables;i++)
      scale+=SQR(dx[i]);
    scale=sqrt(scale);

    // Reduce the step size
    if(scale>MaximumStepLength)
      for(i=0;i<NumberOfVariables;i++)
        dx[i]*=(MaximumStepLength/scale);


    StepSize=0.0;
    for(i=0;i<NumberOfVariables;i++)
      StepSize+=SQR(dx[i]);
    StepSize=sqrt(StepSize);
/*
    if(StepSize>1.5*StepSizePrevious)
    {
      fprintf(stderr, "trimming back: %g %g\n",StepSize,StepSizePrevious);
      for(i=0;i<NumberOfVariables;i++)
        dx[i]*=(1.5*StepSizePrevious/StepSize);
    }
*/
    StepSizePrevious=StepSize;

    // Find the RMS gradient.
    RMSGradient=0.0;
    MaxGradient=0.0;
    for(i=0;i<NumberOfVariables;i++)
    {
      if(fabs(Gradient[i])>MaxGradient) MaxGradient=fabs(Gradient[i]);
      RMSGradient+=SQR(Gradient[i]);
    }
    RMSGradient=sqrt(RMSGradient)/(REAL)NumberOfVariables;

    switch(Dimension)
    {
      case 2:
        fprintf(OutputFilePtr[CurrentSystem],"Iteration: %d Energy: %18.10f Area: %18.10f RMS gradient: %g  Max gradient: %g Number of negative eigenvalues: %d Number of zero eigenvalues: %d\n",
          k,(double)(Energy*ENERGY_TO_KELVIN),(double)Volume[CurrentSystem],(double)RMSGradient,(double)MaxGradient,NumberOfNegativeEigenValues,NumberOfZeroEigenValues);
        fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f     Strain derivative: %18.10f %18.10f\n",
            (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,
            (double)StrainFirstDerivative.ax,(double)StrainFirstDerivative.bx);
        fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f                        %18.10f %18.10f\n",
            (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,
            (double)StrainFirstDerivative.ay,(double)StrainFirstDerivative.by);
        fprintf(OutputFilePtr[CurrentSystem],"      Lengths: %18.10f %18.10f, Angle: %18.10f\n",
           (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,GammaAngle[CurrentSystem]*RAD2DEG);
        break;
      case 3:
        fprintf(OutputFilePtr[CurrentSystem],"Iteration: %d Energy: %18.10f Volume: %18.10f RMS gradient: %g  Max gradient: %g Number of negative eigenvalues: %d Number of zero eigenvalues: %d\n",
          k,(double)(Energy*ENERGY_TO_KELVIN),(double)Volume[CurrentSystem],(double)RMSGradient,(double)MaxGradient,NumberOfNegativeEigenValues,NumberOfZeroEigenValues);
        fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f %18.10f     Strain derivative: %18.10f %18.10f %18.10f\n",
            (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx,
            (double)StrainFirstDerivative.ax,(double)StrainFirstDerivative.bx,(double)StrainFirstDerivative.cx);
        fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
            (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy,
            (double)StrainFirstDerivative.ay,(double)StrainFirstDerivative.by,(double)StrainFirstDerivative.cy);
        fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
            (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz,
            (double)StrainFirstDerivative.az,(double)StrainFirstDerivative.bz,(double)StrainFirstDerivative.cz);
        fprintf(OutputFilePtr[CurrentSystem],"      Lengths: %18.10f %18.10f %18.10f, Angles: %18.10f %18.10f %18.10f\n",
           (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,(double)BoxProperties[CurrentSystem].az,
           AlphaAngle[CurrentSystem]*RAD2DEG,BetaAngle[CurrentSystem]*RAD2DEG,GammaAngle[CurrentSystem]*RAD2DEG);
        break;
    }
    for(i=0;i<NumberOfComponents;i++)
    {
      if(Components[i].RestrictEnantionface)
      {
        posA=Components[i].EnantiofaceAtoms[0]->Position;
        posB=Components[i].EnantiofaceAtoms[1]->Position;
        posC=Components[i].EnantiofaceAtoms[2]->Position;
        posD=Components[i].EnantiofaceAtoms[3]->Position;
        posE=Components[i].EnantiofaceAtoms[4]->Position;
        fprintf(OutputFilePtr[CurrentSystem],"Component: %d, enantionface should be: %s, current enantioface: %s\n",i,
         Components[i].Enantioface==ENANTIOFACE_RE?"Re":"Si",
         CheckEnantioFace(posA,posB,posC,posD,posE)==ENANTIOFACE_RE?"Re":"Si");
      }
    }
    fprintf(OutputFilePtr[CurrentSystem],"\n");
    fflush(OutputFilePtr[CurrentSystem]);

    // write the box and positions to a restart-file
    PrintRestartFile();

    // Check for convergence.
    if((RMSGradient<RMSGradientTolerance)&&(MaxGradient<MaxGradientTolerance))
    //if(((RMSGradient<RMSGradientTolerance)&&(MaxGradient<MaxGradientTolerance)))
    {
      fprintf(OutputFilePtr[CurrentSystem],"SUCCES: RMS Gradient tolerance %g reached (%g)\n",(double)RMSGradientTolerance,(double)RMSGradient);
      fprintf(OutputFilePtr[CurrentSystem],"        Max Gradient tolerance %g reached (%g)\n",(double)MaxGradientTolerance,(double)MaxGradient);
      break;
    }
    CurrentSystem=0;
    if(k%PrintEvery==0)
      SamplePDBMovies(SAMPLE,run);
  }
  SamplePDBMovies(FINALIZE,run);

  // if upper-triangle form, convert the cell (and atom positions) to upper-triangle
  //if(Ensemble[CurrentSystem]==NPTPR&&(NPTPRCellType[CurrentSystem]==REGULAR_UPPER_TRIANGLE||NPTPRCellType[CurrentSystem]==MONOCLINIC_UPPER_TRIANGLE))
  //  AdjustToUpperTriangle();

  ComputeDerivativesMinimization(np,nb,x,&Energy,Gradient,HessianMatrix,&StrainFirstDerivative);

  // project the constraints from the gradient and Hessian
  ProjectConstraintsFromHessianMatrix(np,nb,Gradient,HessianMatrix,TRUE,TRUE);

  if(MinimizationVariables==FRACTIONAL)
    ConvertHessianFromCartesianToFractional(np,nb,Gradient,HessianMatrix);

  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
  fprintf(OutputFilePtr[CurrentSystem],"Variable vector:\n");
  fprintf(OutputFilePtr[CurrentSystem],"================\n");
  for(i=0;i<NumberOfVariables;i++)
    fprintf(OutputFilePtr[CurrentSystem],"x: %d = %g\n",i,(double)x[i]);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
  fprintf(OutputFilePtr[CurrentSystem],"Gradient vector:\n");
  fprintf(OutputFilePtr[CurrentSystem],"================\n");
  for(i=0;i<NumberOfVariables;i++)
    fprintf(OutputFilePtr[CurrentSystem],"Gradient: %d = %g\n",i,(double)Gradient[i]);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  SolveEigenValuesAndVectorsHessian(HessianMatrix,EigenValues);

  fprintf(OutputFilePtr[CurrentSystem],"Eigenvalues:\n");
  fprintf(OutputFilePtr[CurrentSystem],"============\n");
  for(i=0;i<NumberOfVariables;i++)
    fprintf(OutputFilePtr[CurrentSystem],"Eigenvalue: %d = %18.10f\n",i,(double)EigenValues[i]);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  fprintf(OutputFilePtr[CurrentSystem],"Final energy after minimization: %18.10f in %d steps\n",(double)(Energy*ENERGY_TO_KELVIN),k);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  fprintf(OutputFilePtr[CurrentSystem],"Volume [A^3]: %18.10f\n",Volume[CurrentSystem]);
  fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f %18.10f     Strain derivative: %18.10f %18.10f %18.10f\n",
        (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx,
        (double)StrainFirstDerivative.ax,(double)StrainFirstDerivative.bx,(double)StrainFirstDerivative.cx);
  fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
        (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy,
        (double)StrainFirstDerivative.ay,(double)StrainFirstDerivative.by,(double)StrainFirstDerivative.cy);
  fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
        (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz,
        (double)StrainFirstDerivative.az,(double)StrainFirstDerivative.bz,(double)StrainFirstDerivative.cz);
  fprintf(OutputFilePtr[CurrentSystem],"      Final Lengths: %18.10f %18.10f %18.10f, Angles: %18.10f %18.10f %18.10f\n",
       (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,(double)BoxProperties[CurrentSystem].az,
       AlphaAngle[CurrentSystem]*RAD2DEG,BetaAngle[CurrentSystem]*RAD2DEG,GammaAngle[CurrentSystem]*RAD2DEG);
  fprintf(OutputFilePtr[CurrentSystem],"\n");

  ComputeBornTerm=FALSE;
  CalculateForce();

  ShakeInMinimization();

  fprintf(OutputFilePtr[CurrentSystem],"Strain derivative from force routines (NOTE: different when fixed atoms are present)\n");
  fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %18.10f\n",StrainDerivativeTensor[CurrentSystem].ax,
          StrainDerivativeTensor[CurrentSystem].bx,StrainDerivativeTensor[CurrentSystem].cx);
  fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %18.10f\n",StrainDerivativeTensor[CurrentSystem].ay,
          StrainDerivativeTensor[CurrentSystem].by,StrainDerivativeTensor[CurrentSystem].cy);
  fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %18.10f\n",StrainDerivativeTensor[CurrentSystem].az,
          StrainDerivativeTensor[CurrentSystem].bz,StrainDerivativeTensor[CurrentSystem].cz);
  fprintf(OutputFilePtr[CurrentSystem],"\n");

  fprintf(OutputFilePtr[CurrentSystem],"Forces\n");
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Framework[%d] Atom: %d  %g %g %g\n",f1,i,
        Framework[CurrentSystem].Atoms[f1][i].Fixed.x?0.0:Framework[CurrentSystem].Atoms[f1][i].Force.x,
        Framework[CurrentSystem].Atoms[f1][i].Fixed.y?0.0:Framework[CurrentSystem].Atoms[f1][i].Force.y,
        Framework[CurrentSystem].Atoms[f1][i].Fixed.z?0.0:Framework[CurrentSystem].Atoms[f1][i].Force.z);
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        fprintf(OutputFilePtr[CurrentSystem],"Adsorbate[%d] Group: %d  Com force: %g %g %g\n",i,l,
               Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass.x?0.0:Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.x,
               Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass.y?0.0:Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.y,
               Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass.z?0.0:Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.z);
        fprintf(OutputFilePtr[CurrentSystem],"Adsorbate[%d] Group: %d  torque: %g %g %g\n",i,l,
               Adsorbates[CurrentSystem][i].Groups[l].FixedOrientation.x?0.0:Adsorbates[CurrentSystem][i].Groups[l].Torque.x,
               Adsorbates[CurrentSystem][i].Groups[l].FixedOrientation.y?0.0:Adsorbates[CurrentSystem][i].Groups[l].Torque.y,
               Adsorbates[CurrentSystem][i].Groups[l].FixedOrientation.z?0.0:Adsorbates[CurrentSystem][i].Groups[l].Torque.z);
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          fprintf(OutputFilePtr[CurrentSystem],"Adsorbate[%d] Atom: %d  %g %g %g\n",i,A,
            Adsorbates[CurrentSystem][i].Atoms[A].Fixed.x?0.0:Adsorbates[CurrentSystem][i].Atoms[A].Force.x,
            Adsorbates[CurrentSystem][i].Atoms[A].Fixed.y?0.0:Adsorbates[CurrentSystem][i].Atoms[A].Force.y,
            Adsorbates[CurrentSystem][i].Atoms[A].Fixed.z?0.0:Adsorbates[CurrentSystem][i].Atoms[A].Force.z);
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        fprintf(OutputFilePtr[CurrentSystem],"Cation[%d] Group: %d  Com force: %g %g %g\n",i,l,
               Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass.x?0.0:Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.x,
               Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass.y?0.0:Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.y,
               Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass.z?0.0:Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.z);
        fprintf(OutputFilePtr[CurrentSystem],"Cation[%d] Group: %d  torque: %g %g %g\n",i,l,
               Cations[CurrentSystem][i].Groups[l].FixedOrientation.x?0.0:Cations[CurrentSystem][i].Groups[l].Torque.x,
               Cations[CurrentSystem][i].Groups[l].FixedOrientation.y?0.0:Cations[CurrentSystem][i].Groups[l].Torque.y,
               Cations[CurrentSystem][i].Groups[l].FixedOrientation.z?0.0:Cations[CurrentSystem][i].Groups[l].Torque.z);
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          fprintf(OutputFilePtr[CurrentSystem],"Cation[%d] Atom: %d  %g %g %g\n",i,A,
              Cations[CurrentSystem][i].Atoms[A].Fixed.x?0.0:Cations[CurrentSystem][i].Atoms[A].Force.x,
              Cations[CurrentSystem][i].Atoms[A].Fixed.y?0.0:Cations[CurrentSystem][i].Atoms[A].Force.y,
              Cations[CurrentSystem][i].Atoms[A].Fixed.z?0.0:Cations[CurrentSystem][i].Atoms[A].Force.z);
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  ComputeElasticConstantsAfterMinimization(np,nb,x);
}

void ComputeElasticConstantsAfterMinimization(int np,int nb,REAL *x)
{
  int i,j,k,l;
  int NumberOfVariables;
  REAL Energy,fac;
  REAL *Gradient,*EigenValues;
  REAL_MATRIX HessianMatrix,HessianMatrixInverse,CrossTerm;
  REAL_MATRIX9x9 RelaxationTerm,ElasticConstants,PressureCorrection;
  REAL_MATRIX6x6 VoigtMatrix3D,Compliances3D;
  REAL_MATRIX3x3 VoigtMatrix2D,Compliances2D;
  REAL StabilityEigenValues[6];
  REAL KVoigt,GVoigt;
  REAL KReuss,GReuss;
  REAL KHill,GHill;
  REAL_MATRIX3x3 StrainFirstDerivative;
  REAL Factor;

  NumberOfVariables=np+nb;
  Gradient=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  EigenValues=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  HessianMatrix=CreateRealMatrix(NumberOfVariables,NumberOfVariables);

  if(ComputeElasticConstants)
  {
    if((Ensemble[CurrentSystem]==NPTPR)&&(NPTPRCellType[CurrentSystem]==REGULAR_UPPER_TRIANGLE||NPTPRCellType[CurrentSystem]==REGULAR))
    {
      ComputeDerivativesMinimization(np,nb,x,&Energy,Gradient,HessianMatrix,&StrainFirstDerivative);

      // project the constraints from the gradient and Hessian
      ProjectConstraintsFromHessianMatrix(np,nb,Gradient,HessianMatrix,TRUE,TRUE);

      if(MinimizationVariables==FRACTIONAL)
        ConvertHessianFromCartesianToFractional(np,nb,Gradient,HessianMatrix);

      switch(Dimension)
      {
        case 3:
          BornTerm[CurrentSystem].xzxx=HessianMatrix.element[np][np+2];
          BornTerm[CurrentSystem].yzxx=HessianMatrix.element[np][np+4];
          BornTerm[CurrentSystem].zxxx=HessianMatrix.element[np][np+2];
          BornTerm[CurrentSystem].zyxx=HessianMatrix.element[np][np+4];
          BornTerm[CurrentSystem].zzxx=HessianMatrix.element[np][np+5];

          BornTerm[CurrentSystem].xzxy=HessianMatrix.element[np+1][np+2];
          BornTerm[CurrentSystem].yzxy=HessianMatrix.element[np+1][np+4];
          BornTerm[CurrentSystem].zxxy=HessianMatrix.element[np+1][np+2];
          BornTerm[CurrentSystem].zyxy=HessianMatrix.element[np+1][np+4];
          BornTerm[CurrentSystem].zzxy=HessianMatrix.element[np+1][np+5];

          BornTerm[CurrentSystem].xxxz=HessianMatrix.element[np+2][np];
          BornTerm[CurrentSystem].xyxz=HessianMatrix.element[np+2][np+1];
          BornTerm[CurrentSystem].xzxz=HessianMatrix.element[np+2][np+2];
          BornTerm[CurrentSystem].yxxz=HessianMatrix.element[np+2][np+1];
          BornTerm[CurrentSystem].yyxz=HessianMatrix.element[np+2][np+3];
          BornTerm[CurrentSystem].yzxz=HessianMatrix.element[np+2][np+4];
          BornTerm[CurrentSystem].zxxz=HessianMatrix.element[np+2][np+2];
          BornTerm[CurrentSystem].zyxz=HessianMatrix.element[np+2][np+4];
          BornTerm[CurrentSystem].zzxz=HessianMatrix.element[np+2][np+5];

          BornTerm[CurrentSystem].xzyx=HessianMatrix.element[np+1][np+2];
          BornTerm[CurrentSystem].yzyx=HessianMatrix.element[np+1][np+4];
          BornTerm[CurrentSystem].zxyx=HessianMatrix.element[np+1][np+2];
          BornTerm[CurrentSystem].zyyx=HessianMatrix.element[np+1][np+4];
          BornTerm[CurrentSystem].zzyx=HessianMatrix.element[np+1][np+5];

          BornTerm[CurrentSystem].xzyy=HessianMatrix.element[np+3][np+2];
          BornTerm[CurrentSystem].yzyy=HessianMatrix.element[np+3][np+4];
          BornTerm[CurrentSystem].zxyy=HessianMatrix.element[np+3][np+2];
          BornTerm[CurrentSystem].zyyy=HessianMatrix.element[np+3][np+4];
          BornTerm[CurrentSystem].zzyy=HessianMatrix.element[np+3][np+5];

          BornTerm[CurrentSystem].xxyz=HessianMatrix.element[np+4][np];
          BornTerm[CurrentSystem].xyyz=HessianMatrix.element[np+4][np+1];
          BornTerm[CurrentSystem].xzyz=HessianMatrix.element[np+4][np+2];
          BornTerm[CurrentSystem].yxyz=HessianMatrix.element[np+4][np+1];
          BornTerm[CurrentSystem].yyyz=HessianMatrix.element[np+4][np+3];
          BornTerm[CurrentSystem].yzyz=HessianMatrix.element[np+4][np+4];
          BornTerm[CurrentSystem].zxyz=HessianMatrix.element[np+4][np+2];
          BornTerm[CurrentSystem].zyyz=HessianMatrix.element[np+4][np+4];
          BornTerm[CurrentSystem].zzyz=HessianMatrix.element[np+4][np+5];

          BornTerm[CurrentSystem].xxzx=HessianMatrix.element[np+2][np];
          BornTerm[CurrentSystem].xyzx=HessianMatrix.element[np+2][np+1];
          BornTerm[CurrentSystem].xzzx=HessianMatrix.element[np+2][np+2];
          BornTerm[CurrentSystem].yxzx=HessianMatrix.element[np+2][np+1];
          BornTerm[CurrentSystem].yyzx=HessianMatrix.element[np+2][np+3];
          BornTerm[CurrentSystem].yzzx=HessianMatrix.element[np+2][np+4];
          BornTerm[CurrentSystem].zxzx=HessianMatrix.element[np+2][np+2];
          BornTerm[CurrentSystem].zyzx=HessianMatrix.element[np+2][np+4];
          BornTerm[CurrentSystem].zzzx=HessianMatrix.element[np+2][np+5];

          BornTerm[CurrentSystem].xxzy=HessianMatrix.element[np+4][np];
          BornTerm[CurrentSystem].xyzy=HessianMatrix.element[np+4][np+1];
          BornTerm[CurrentSystem].xzzy=HessianMatrix.element[np+4][np+2];
          BornTerm[CurrentSystem].yxzy=HessianMatrix.element[np+4][np+1];
          BornTerm[CurrentSystem].yyzy=HessianMatrix.element[np+4][np+3];
          BornTerm[CurrentSystem].yzzy=HessianMatrix.element[np+4][np+4];
          BornTerm[CurrentSystem].zxzy=HessianMatrix.element[np+4][np+2];
          BornTerm[CurrentSystem].zyzy=HessianMatrix.element[np+4][np+4];
          BornTerm[CurrentSystem].zzzy=HessianMatrix.element[np+4][np+5];

          BornTerm[CurrentSystem].xxzz=HessianMatrix.element[np+5][np];
          BornTerm[CurrentSystem].xyzz=HessianMatrix.element[np+5][np+1];
          BornTerm[CurrentSystem].xzzz=HessianMatrix.element[np+5][np+2];
          BornTerm[CurrentSystem].yxzz=HessianMatrix.element[np+5][np+1];
          BornTerm[CurrentSystem].yyzz=HessianMatrix.element[np+5][np+3];
          BornTerm[CurrentSystem].yzzz=HessianMatrix.element[np+5][np+4];
          BornTerm[CurrentSystem].zxzz=HessianMatrix.element[np+5][np+2];
          BornTerm[CurrentSystem].zyzz=HessianMatrix.element[np+5][np+4];
          BornTerm[CurrentSystem].zzzz=HessianMatrix.element[np+5][np+5];
        case 2:
          BornTerm[CurrentSystem].xyxx=HessianMatrix.element[np][np+1];
          BornTerm[CurrentSystem].yxxx=HessianMatrix.element[np][np+1];
          BornTerm[CurrentSystem].yyxx=HessianMatrix.element[np][np+Dimension];
          BornTerm[CurrentSystem].xxxy=HessianMatrix.element[np+1][np];
          BornTerm[CurrentSystem].xyxy=HessianMatrix.element[np+1][np+1];
          BornTerm[CurrentSystem].yxxy=HessianMatrix.element[np+1][np+1];
          BornTerm[CurrentSystem].yyxy=HessianMatrix.element[np+1][np+Dimension];
          BornTerm[CurrentSystem].xxyx=HessianMatrix.element[np+1][np];
          BornTerm[CurrentSystem].xyyx=HessianMatrix.element[np+1][np+1];
          BornTerm[CurrentSystem].yxyx=HessianMatrix.element[np+1][np+1];
          BornTerm[CurrentSystem].yyyx=HessianMatrix.element[np+1][np+Dimension];
          BornTerm[CurrentSystem].xxyy=HessianMatrix.element[np+Dimension][np];
          BornTerm[CurrentSystem].xyyy=HessianMatrix.element[np+Dimension][np+1];
          BornTerm[CurrentSystem].yxyy=HessianMatrix.element[np+Dimension][np+1];
          BornTerm[CurrentSystem].yyyy=HessianMatrix.element[np+Dimension][np+Dimension];
        case 1:
          BornTerm[CurrentSystem].xxxx=HessianMatrix.element[np][np];
          break;
      }

      // inverse is singular, ommit last particle to remove com-translation
      // the contribution of this particle is still contained in the remainder
      HessianMatrixInverse=CreateRealMatrix(np-Dimension,np-Dimension);
      CrossTerm=CreateRealMatrix(np,SQR(Dimension));

      for(i=0;i<np-Dimension;i++)
        for(j=0;j<np-Dimension;j++)
          HessianMatrixInverse.element[i][j]=HessianMatrix.element[i][j];

      for(i=0;i<np/Dimension;i++)
      {
        switch(Dimension)
        {
          case 3:
            CrossTerm.element[3*i][2]=HessianMatrix.element[3*i][np+2];
            CrossTerm.element[3*i+1][2]=HessianMatrix.element[3*i][np+4];
            CrossTerm.element[3*i+2][0]=HessianMatrix.element[3*i][np+2];
            CrossTerm.element[3*i+2][1]=HessianMatrix.element[3*i][np+4];
            CrossTerm.element[3*i+2][2]=HessianMatrix.element[3*i][np+5];

            CrossTerm.element[3*i][5]=HessianMatrix.element[3*i+1][np+2];
            CrossTerm.element[3*i+1][5]=HessianMatrix.element[3*i+1][np+4];
            CrossTerm.element[3*i+2][3]=HessianMatrix.element[3*i+1][np+2];
            CrossTerm.element[3*i+2][4]=HessianMatrix.element[3*i+1][np+4];
            CrossTerm.element[3*i+2][5]=HessianMatrix.element[3*i+1][np+5];

            CrossTerm.element[3*i][6]=HessianMatrix.element[3*i+2][np];
            CrossTerm.element[3*i][7]=HessianMatrix.element[3*i+2][np+1];
            CrossTerm.element[3*i][8]=HessianMatrix.element[3*i+2][np+2];
            CrossTerm.element[3*i+1][6]=HessianMatrix.element[3*i+2][np+1];
            CrossTerm.element[3*i+1][7]=HessianMatrix.element[3*i+2][np+3];
            CrossTerm.element[3*i+1][8]=HessianMatrix.element[3*i+2][np+4];
            CrossTerm.element[3*i+2][6]=HessianMatrix.element[3*i+2][np+2];
            CrossTerm.element[3*i+2][7]=HessianMatrix.element[3*i+2][np+4];
            CrossTerm.element[3*i+2][8]=HessianMatrix.element[3*i+2][np+5];
          case 2:
            CrossTerm.element[Dimension*i][1]=HessianMatrix.element[Dimension*i][np+1];
            CrossTerm.element[Dimension*i+1][0]=HessianMatrix.element[Dimension*i][np+1];
            CrossTerm.element[Dimension*i+1][1]=HessianMatrix.element[Dimension*i][np+Dimension];
            CrossTerm.element[Dimension*i][Dimension]=HessianMatrix.element[Dimension*i+1][np];
            CrossTerm.element[Dimension*i][Dimension+1]=HessianMatrix.element[Dimension*i+1][np+1];
            CrossTerm.element[Dimension*i+1][Dimension]=HessianMatrix.element[Dimension*i+1][np+1];
            CrossTerm.element[Dimension*i+1][Dimension+1]=HessianMatrix.element[Dimension*i+1][np+Dimension];
          case 1:
            CrossTerm.element[Dimension*i][0]=HessianMatrix.element[Dimension*i][np];
            break;
        }
      }


      //InverseRealMatrix(HessianMatrixInverse);
      SingularValueDecompositionMatrixInversion(HessianMatrixInverse);

      InitializeMatrix9x9(&RelaxationTerm);

      fac=1.0;
      for(i=0;i<(np-Dimension)/Dimension;i++)
      {
        for(j=0;j<(np-Dimension)/Dimension;j++)
        {
          for(k=0;k<Dimension;k++)
          {
            for(l=0;l<Dimension;l++)
            {
              switch(Dimension)
              {
                case 3:
                  RelaxationTerm.zxxx+=fac*CrossTerm.element[3*i][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l];
                  RelaxationTerm.zxyx+=fac*CrossTerm.element[3*i][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+1];
                  RelaxationTerm.xxzx+=fac*CrossTerm.element[3*i][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+2];
                  RelaxationTerm.yxzx+=fac*CrossTerm.element[3*i][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+2];
                  RelaxationTerm.zxzx+=fac*CrossTerm.element[3*i][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+2];
                  RelaxationTerm.zyxx+=fac*CrossTerm.element[3*i+1][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l];
                  RelaxationTerm.zyyx+=fac*CrossTerm.element[3*i+1][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+1];
                  RelaxationTerm.xyzx+=fac*CrossTerm.element[3*i+1][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+2];
                  RelaxationTerm.yyzx+=fac*CrossTerm.element[3*i+1][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+2];
                  RelaxationTerm.zyzx+=fac*CrossTerm.element[3*i+1][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+2];

                  RelaxationTerm.xzxx+=fac*CrossTerm.element[3*i+2][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l];
                  RelaxationTerm.yzxx+=fac*CrossTerm.element[3*i+2][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l];
                  RelaxationTerm.zzxx+=fac*CrossTerm.element[3*i+2][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l];
                  RelaxationTerm.xzyx+=fac*CrossTerm.element[3*i+2][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+1];
                  RelaxationTerm.yzyx+=fac*CrossTerm.element[3*i+2][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+1];
                  RelaxationTerm.zzyx+=fac*CrossTerm.element[3*i+2][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+1];
                  RelaxationTerm.xzzx+=fac*CrossTerm.element[3*i+2][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+2];
                  RelaxationTerm.yzzx+=fac*CrossTerm.element[3*i+2][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+2];
                  RelaxationTerm.zzzx+=fac*CrossTerm.element[3*i+2][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j][3*l+2];

                  RelaxationTerm.zxxy+=fac*CrossTerm.element[3*i][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l];
                  RelaxationTerm.zxyy+=fac*CrossTerm.element[3*i][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+1];
                  RelaxationTerm.xxzy+=fac*CrossTerm.element[3*i][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+2];
                  RelaxationTerm.yxzy+=fac*CrossTerm.element[3*i][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+2];
                  RelaxationTerm.zxzy+=fac*CrossTerm.element[3*i][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+2];

                  RelaxationTerm.zyxy+=fac*CrossTerm.element[3*i+1][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l];
                  RelaxationTerm.zyyy+=fac*CrossTerm.element[3*i+1][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+1];
                  RelaxationTerm.xyzy+=fac*CrossTerm.element[3*i+1][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+2];
                  RelaxationTerm.yyzy+=fac*CrossTerm.element[3*i+1][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+2];
                  RelaxationTerm.zyzy+=fac*CrossTerm.element[3*i+1][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+2];

                  RelaxationTerm.xzxy+=fac*CrossTerm.element[3*i+2][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l];
                  RelaxationTerm.yzxy+=fac*CrossTerm.element[3*i+2][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l];
                  RelaxationTerm.zzxy+=fac*CrossTerm.element[3*i+2][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l];
                  RelaxationTerm.xzyy+=fac*CrossTerm.element[3*i+2][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+1];
                  RelaxationTerm.yzyy+=fac*CrossTerm.element[3*i+2][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+1];
                  RelaxationTerm.zzyy+=fac*CrossTerm.element[3*i+2][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+1];
                  RelaxationTerm.xzzy+=fac*CrossTerm.element[3*i+2][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+2];
                  RelaxationTerm.yzzy+=fac*CrossTerm.element[3*i+2][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+2];
                  RelaxationTerm.zzzy+=fac*CrossTerm.element[3*i+2][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+1][3*l+2];

                  RelaxationTerm.xxxz+=fac*CrossTerm.element[3*i][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l];
                  RelaxationTerm.yxxz+=fac*CrossTerm.element[3*i][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l];
                  RelaxationTerm.zxxz+=fac*CrossTerm.element[3*i][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l];
                  RelaxationTerm.xxyz+=fac*CrossTerm.element[3*i][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+1];
                  RelaxationTerm.yxyz+=fac*CrossTerm.element[3*i][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+1];
                  RelaxationTerm.zxyz+=fac*CrossTerm.element[3*i][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+1];
                  RelaxationTerm.xxzz+=fac*CrossTerm.element[3*i][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+2];
                  RelaxationTerm.yxzz+=fac*CrossTerm.element[3*i][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+2];
                  RelaxationTerm.zxzz+=fac*CrossTerm.element[3*i][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+2];

                  RelaxationTerm.xyxz+=fac*CrossTerm.element[3*i+1][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l];
                  RelaxationTerm.yyxz+=fac*CrossTerm.element[3*i+1][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l];
                  RelaxationTerm.zyxz+=fac*CrossTerm.element[3*i+1][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l];
                  RelaxationTerm.xyyz+=fac*CrossTerm.element[3*i+1][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+1];
                  RelaxationTerm.yyyz+=fac*CrossTerm.element[3*i+1][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+1];
                  RelaxationTerm.zyyz+=fac*CrossTerm.element[3*i+1][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+1];
                  RelaxationTerm.xyzz+=fac*CrossTerm.element[3*i+1][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+2];
                  RelaxationTerm.yyzz+=fac*CrossTerm.element[3*i+1][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+2];
                  RelaxationTerm.zyzz+=fac*CrossTerm.element[3*i+1][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+2];

                  RelaxationTerm.xzxz+=fac*CrossTerm.element[3*i+2][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l];
                  RelaxationTerm.yzxz+=fac*CrossTerm.element[3*i+2][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l];
                  RelaxationTerm.zzxz+=fac*CrossTerm.element[3*i+2][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l];
                  RelaxationTerm.xzyz+=fac*CrossTerm.element[3*i+2][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+1];
                  RelaxationTerm.yzyz+=fac*CrossTerm.element[3*i+2][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+1];
                  RelaxationTerm.zzyz+=fac*CrossTerm.element[3*i+2][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+1];
                  RelaxationTerm.xzzz+=fac*CrossTerm.element[3*i+2][3*k]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+2];
                  RelaxationTerm.yzzz+=fac*CrossTerm.element[3*i+2][3*k+1]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+2];
                  RelaxationTerm.zzzz+=fac*CrossTerm.element[3*i+2][3*k+2]*HessianMatrixInverse.element[3*i+k][3*j+l]*CrossTerm.element[3*j+2][3*l+2];
                case 2:
                  RelaxationTerm.yxxx+=fac*CrossTerm.element[Dimension*i][Dimension*k+1]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j][Dimension*l];
                  RelaxationTerm.xxyx+=fac*CrossTerm.element[Dimension*i][Dimension*k]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j][Dimension*l+1];
                  RelaxationTerm.yxyx+=fac*CrossTerm.element[Dimension*i][Dimension*k+1]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j][Dimension*l+1];
                  RelaxationTerm.xyxx+=fac*CrossTerm.element[Dimension*i+1][Dimension*k]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j][Dimension*l];
                  RelaxationTerm.yyxx+=fac*CrossTerm.element[Dimension*i+1][Dimension*k+1]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j][Dimension*l];
                  RelaxationTerm.xyyx+=fac*CrossTerm.element[Dimension*i+1][Dimension*k]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j][Dimension*l+1];
                  RelaxationTerm.yyyx+=fac*CrossTerm.element[Dimension*i+1][Dimension*k+1]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j][Dimension*l+1];
                  RelaxationTerm.xxxy+=fac*CrossTerm.element[Dimension*i][Dimension*k]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j+1][Dimension*l];
                  RelaxationTerm.yxxy+=fac*CrossTerm.element[Dimension*i][Dimension*k+1]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j+1][Dimension*l];
                  RelaxationTerm.xxyy+=fac*CrossTerm.element[Dimension*i][Dimension*k]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j+1][Dimension*l+1];
                  RelaxationTerm.yxyy+=fac*CrossTerm.element[Dimension*i][Dimension*k+1]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j+1][Dimension*l+1];
                  RelaxationTerm.xyxy+=fac*CrossTerm.element[Dimension*i+1][Dimension*k]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j+1][Dimension*l];
                  RelaxationTerm.yyxy+=fac*CrossTerm.element[Dimension*i+1][Dimension*k+1]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j+1][Dimension*l];
                  RelaxationTerm.xyyy+=fac*CrossTerm.element[Dimension*i+1][Dimension*k]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j+1][Dimension*l+1];
                  RelaxationTerm.yyyy+=fac*CrossTerm.element[Dimension*i+1][Dimension*k+1]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j+1][Dimension*l+1];
                case 1:
                  RelaxationTerm.xxxx+=fac*CrossTerm.element[Dimension*i][Dimension*k]*HessianMatrixInverse.element[Dimension*i+k][Dimension*j+l]*CrossTerm.element[Dimension*j][Dimension*l];
                  break;
              }
            }
          }
        }
      }

      InitializeMatrix9x9(&PressureCorrection);

      PressureCorrection.xxxx-=therm_baro_stats.ExternalPressure[CurrentSystem][0]*Volume[CurrentSystem];
      PressureCorrection.yyyy-=therm_baro_stats.ExternalPressure[CurrentSystem][0]*Volume[CurrentSystem];
      PressureCorrection.zzzz-=therm_baro_stats.ExternalPressure[CurrentSystem][0]*Volume[CurrentSystem];

      PressureCorrection.xyxy-=0.5*therm_baro_stats.ExternalPressure[CurrentSystem][0]*Volume[CurrentSystem];
      PressureCorrection.xzxz-=0.5*therm_baro_stats.ExternalPressure[CurrentSystem][0]*Volume[CurrentSystem];
      PressureCorrection.yzyz-=0.5*therm_baro_stats.ExternalPressure[CurrentSystem][0]*Volume[CurrentSystem];

      if(UseReducedUnits)
        Factor=1.0;
      else
        Factor=1e9/PRESSURE_CONVERSION_FACTOR;

      fprintf(OutputFilePtr[CurrentSystem],"Born-term %s\n",UseReducedUnits?"[-]":"[GPa]");
      fprintf(OutputFilePtr[CurrentSystem],"---------------\n");
      PrintRealMatrix9x9ToFile(&BornTerm[CurrentSystem],OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
      fprintf(OutputFilePtr[CurrentSystem],"\n");


      fprintf(OutputFilePtr[CurrentSystem],"Relaxation-term %s\n",UseReducedUnits?"[-]":"[GPa]");
      fprintf(OutputFilePtr[CurrentSystem],"-----------------------------\n");
      PrintRealMatrix9x9ToFile(&RelaxationTerm,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
      fprintf(OutputFilePtr[CurrentSystem],"\n\n");

      fprintf(OutputFilePtr[CurrentSystem],"Pressure-Corection-term %s\n",UseReducedUnits?"[-]":"[GPa]");
      fprintf(OutputFilePtr[CurrentSystem],"-----------------------------\n");
      PrintRealMatrix9x9ToFile(&PressureCorrection,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
      fprintf(OutputFilePtr[CurrentSystem],"\n\n");


      SubtractRealMatrix9x9(&ElasticConstants,BornTerm[CurrentSystem],RelaxationTerm);
      AddRealMatrix9x9(&ElasticConstants,ElasticConstants,PressureCorrection);

      fprintf(OutputFilePtr[CurrentSystem],"Elastic constant %s\n",UseReducedUnits?"[-]":"[GPa]");
      fprintf(OutputFilePtr[CurrentSystem],"----------------------\n");
      PrintRealMatrix9x9ToFile(&ElasticConstants,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
      fprintf(OutputFilePtr[CurrentSystem],"\n");

      fprintf(OutputFilePtr[CurrentSystem],"Born-term (Voigt notation) %s\n",UseReducedUnits?"[-]":"[GPa]");
      fprintf(OutputFilePtr[CurrentSystem],"-------------------------------\n");
      switch(Dimension)
      {
        case 2:
          VoigtMatrix2D=ConvertToVoigt2D(BornTerm[CurrentSystem]);
          PrintRealMatrix3x3ToFile(&VoigtMatrix2D,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
          break;
        case 3:
          VoigtMatrix3D=ConvertToVoigt3D(BornTerm[CurrentSystem]);
          PrintRealMatrix6x6ToFile(&VoigtMatrix3D,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
          break;
      }
      fprintf(OutputFilePtr[CurrentSystem],"\n");

      fprintf(OutputFilePtr[CurrentSystem],"Relaxation-term (Voigt notation) %s\n",UseReducedUnits?"[-]":"[GPa]");
      fprintf(OutputFilePtr[CurrentSystem],"-------------------------------------\n");
      switch(Dimension)
      {
        case 2:
          VoigtMatrix2D=ConvertToVoigt2D(RelaxationTerm);
          PrintRealMatrix3x3ToFile(&VoigtMatrix2D,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
          break;
        case 3:
          VoigtMatrix3D=ConvertToVoigt3D(RelaxationTerm);
          PrintRealMatrix6x6ToFile(&VoigtMatrix3D,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
          break;
      }
      fprintf(OutputFilePtr[CurrentSystem],"\n");

      fprintf(OutputFilePtr[CurrentSystem],"Pressure-Correction-term (Voigt notation) %s\n",UseReducedUnits?"[-]":"[GPa]");
      fprintf(OutputFilePtr[CurrentSystem],"-------------------------------------\n");
      switch(Dimension)
      {
        case 2:
          VoigtMatrix2D=ConvertToVoigt2D(PressureCorrection);
          PrintRealMatrix3x3ToFile(&VoigtMatrix2D,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
          break;
        case 3:
          VoigtMatrix3D=ConvertToVoigt3D(PressureCorrection);
          PrintRealMatrix6x6ToFile(&VoigtMatrix3D,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
          break;
      }
      fprintf(OutputFilePtr[CurrentSystem],"\n");

      fprintf(OutputFilePtr[CurrentSystem],"Elastic constant (Voigt notation) %s\n",UseReducedUnits?"[-]":"[GPa]");
      fprintf(OutputFilePtr[CurrentSystem],"------------------------------\n");
      switch(Dimension)
      {
        case 2:
          VoigtMatrix2D=ConvertToVoigt2D(ElasticConstants);
          PrintRealMatrix3x3ToFile(&VoigtMatrix2D,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
          break;
        case 3:
          VoigtMatrix3D=ConvertToVoigt3D(ElasticConstants);
          PrintRealMatrix6x6ToFile(&VoigtMatrix3D,OutputFilePtr[CurrentSystem],Factor*Volume[CurrentSystem]);
          break;
      }
      fprintf(OutputFilePtr[CurrentSystem],"\n\n");

      switch(Dimension)
      {
        case 2:
          VoigtMatrix2D=ConvertToVoigt2D(ElasticConstants);
          DivideRealMatrix3x3ByReal(&VoigtMatrix2D,VoigtMatrix2D,Factor*Volume[CurrentSystem]);
          InverseMatrix3x3(VoigtMatrix2D,&Compliances2D);
          break;
        case 3:
          VoigtMatrix3D=ConvertToVoigt3D(ElasticConstants);
          DivideRealMatrix6x6ByReal(&VoigtMatrix3D,VoigtMatrix3D,Factor*Volume[CurrentSystem]);
          InverseMatrix6x6(VoigtMatrix3D,&Compliances3D);
          break;
      }


      fprintf(OutputFilePtr[CurrentSystem],"Elastic compliances (Voigt notation) %s\n",UseReducedUnits?"[-]":"[GPa^-1]");
      fprintf(OutputFilePtr[CurrentSystem],"------------------------------------------\n");
      switch(Dimension)
      {
        case 2:
          fprintf(OutputFilePtr[CurrentSystem],"C11 C12 C16:    % 12.5f % 12.5f % 12.5f\n",Compliances2D.ax,Compliances2D.bx,Compliances2D.cx);
          fprintf(OutputFilePtr[CurrentSystem],"C21 C22 C26:    % 12.5f % 12.5f % 12.5f\n",Compliances2D.ay,Compliances2D.by,Compliances2D.cy);
          fprintf(OutputFilePtr[CurrentSystem],"C61 C62 C66:    % 12.5f % 12.5f % 12.5f\n",Compliances2D.az,Compliances2D.bz,Compliances2D.cz);
          fprintf(OutputFilePtr[CurrentSystem],"\n\n");

          fprintf(OutputFilePtr[CurrentSystem],"Young's modulus x: %18.10f %s\n",1.0/Compliances2D.ax,UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"Young's modulus y: %18.10f %s\n\n",1.0/Compliances2D.by,UseReducedUnits?"[-]":"[GPa]");

          fprintf(OutputFilePtr[CurrentSystem],"Poison's ratios %s\n",UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"---------------------\n");

          fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %s\n",
                  -Compliances2D.ax/Compliances2D.ax,-Compliances2D.bx/Compliances2D.by,UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %s\n\n",
                  -Compliances2D.ay/Compliances2D.ax,-Compliances2D.by/Compliances2D.by,UseReducedUnits?"[-]":"[GPa]");

          VoigtMatrix2D=ConvertToVoigt2D(ElasticConstants);
          EigenSystemVoigt3x3(VoigtMatrix2D,StabilityEigenValues);

          fprintf(OutputFilePtr[CurrentSystem],"Born stability criteria %s\n",UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"-----------------------------\n");
          for(i=0;i<3;i++)
            fprintf(OutputFilePtr[CurrentSystem],"Eigenvalue: %18.10f\n",
                    StabilityEigenValues[i]/(Factor*Volume[CurrentSystem]));

          break;
        case 3:
          PrintRealMatrix6x6ToFile(&Compliances3D,OutputFilePtr[CurrentSystem],1.0);
          fprintf(OutputFilePtr[CurrentSystem],"\n\n");

          fprintf(OutputFilePtr[CurrentSystem],"Young's modulus x: %18.10f %s\n",1.0/Compliances3D.C11,UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"Young's modulus y: %18.10f %s\n",1.0/Compliances3D.C22,UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"Young's modulus z: %18.10f %s\n\n",1.0/Compliances3D.C33,UseReducedUnits?"[-]":"[GPa]");

          fprintf(OutputFilePtr[CurrentSystem],"Poison's ratios %s\n",UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"---------------------\n");
          fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %18.10f %s\n",
                  -Compliances3D.C11/Compliances3D.C11,-Compliances3D.C12/Compliances3D.C22,-Compliances3D.C13/Compliances3D.C33,UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %18.10f %s\n",
                  -Compliances3D.C21/Compliances3D.C11,-Compliances3D.C22/Compliances3D.C22,-Compliances3D.C23/Compliances3D.C33,UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"%18.10f %18.10f %18.10f %s\n\n",
                  -Compliances3D.C31/Compliances3D.C11,-Compliances3D.C32/Compliances3D.C22,-Compliances3D.C33/Compliances3D.C33,UseReducedUnits?"[-]":"[GPa]");


          KVoigt=(VoigtMatrix3D.C11+VoigtMatrix3D.C22+VoigtMatrix3D.C33+2.0*(VoigtMatrix3D.C12+VoigtMatrix3D.C13+VoigtMatrix3D.C23))/9.0;
          fprintf(OutputFilePtr[CurrentSystem],"Bulk modulus K (Voight): %18.10f %s\n",KVoigt,UseReducedUnits?"[-]":"[GPa]");
          GVoigt=(VoigtMatrix3D.C11+VoigtMatrix3D.C22+VoigtMatrix3D.C33+3.0*(VoigtMatrix3D.C44+VoigtMatrix3D.C55+VoigtMatrix3D.C66)
                  -VoigtMatrix3D.C12-VoigtMatrix3D.C13-VoigtMatrix3D.C23)/15.0;
          fprintf(OutputFilePtr[CurrentSystem],"Shear modulus G (Voight): %18.10f %s\n\n",GVoigt,UseReducedUnits?"[-]":"[GPa]");

          KReuss=1.0/(Compliances3D.C11+Compliances3D.C22+Compliances3D.C33+2.0*(Compliances3D.C12+Compliances3D.C13+Compliances3D.C23));
          fprintf(OutputFilePtr[CurrentSystem],"Bulk modulus K (Reuss): %18.10f %s\n",KReuss,UseReducedUnits?"[-]":"[GPa]");
          GReuss=15.0/(4.0*(Compliances3D.C11+Compliances3D.C22+Compliances3D.C33-Compliances3D.C12-Compliances3D.C13-Compliances3D.C23)
                       +3.0*(Compliances3D.C44+Compliances3D.C55+Compliances3D.C66));
          fprintf(OutputFilePtr[CurrentSystem],"Shear modulus G (Reuss): %18.10f %s\n\n",GReuss,UseReducedUnits?"[-]":"[GPa]");

          KHill=0.5*(KVoigt+KReuss);
          GHill=0.5*(GVoigt+GReuss);
          fprintf(OutputFilePtr[CurrentSystem],"Bulk modulus K (Hill): %18.10f %s\n",KHill,UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"Shear modulus G (Hill): %18.10f %s\n\n",GHill,UseReducedUnits?"[-]":"[GPa]");

          VoigtMatrix3D=ConvertToVoigt3D(ElasticConstants);
          EigenSystemVoigt6x6(VoigtMatrix3D,StabilityEigenValues);

          fprintf(OutputFilePtr[CurrentSystem],"Born stability criteria %s\n",UseReducedUnits?"[-]":"[GPa]");
          fprintf(OutputFilePtr[CurrentSystem],"-----------------------------\n");
          for(i=0;i<6;i++)
            fprintf(OutputFilePtr[CurrentSystem],"Eigenvalue: %18.10f\n",
                    StabilityEigenValues[i]/(Factor*Volume[CurrentSystem]));
          break;
      }
    }
  }
}

// Articles:
// 1) J.A. Snyman "A new and dynamic method for unconstained minimization"
//    Appl. Math Modelling, 1982 Vol 6, pp. 449-462
// 2) J.A. Snyman "An improved version of the original leapfrog dynamic method
//                 for unconssrained minimzation: LFOP1 (b)"
//    Appl. Math. Modeliing, Vol 7, 1983, pp. 216-218
// 3) J.A. Snyman "A convergent dynamic method for large minimization problems"
//    Computers Math. Applic., Vol 17, no. 10, pp. 1369-1377

// Snyman's method is inspired by classical mechanics. The method simulates the motion of
// a particle and by monitoring its kinetic energy an interfering strategy is adopted which
// ensures that the potential energy is systematically reduced. The method does not suffer
// from problems resulting from "stiffness" of the system or from excessive oscillations of
// the system. It is a practical, robust, and reliable method and it is particularly usefull
// for large particle systems. However, the initial timestep should not be too large.
// The algorithm terminates when the norm of force-vector is smaller than epsilon.

void SnymanMinimization(int np,int nb,REAL *x,int run)
{
  int i,j,f1,loop;
  REAL *a,*a_new,*v,*x_new,*v_new,*x_bar,*v_bar,*a_bar;
  REAL Energy;
  REAL_MATRIX3x3 StrainFirstDerivative;
  REAL DotProduct,RMSGradient,MaxGradient;
  REAL *EigenValues;
  REAL_MATRIX HessianMatrix;
  REAL Pressure;

  CurrentSystem=0;
  SamplePDBMovies(ALLOCATE,run);
  SamplePDBMovies(INITIALIZE,run);

  fprintf(OutputFilePtr[CurrentSystem],"Beginning Snyman minimization:\n");
  fprintf(OutputFilePtr[CurrentSystem],"------------------------------\n");
  fflush(OutputFilePtr[CurrentSystem]);

  x_new=(REAL*)calloc(np+nb,sizeof(REAL));

  a=(REAL*)calloc(np+nb,sizeof(REAL));
  a_new=(REAL*)calloc(np+nb,sizeof(REAL));

  v=(REAL*)calloc(np+nb,sizeof(REAL));
  v_new=(REAL*)calloc(np+nb,sizeof(REAL));

  x_bar=(REAL*)calloc(np+nb,sizeof(REAL));
  v_bar=(REAL*)calloc(np+nb,sizeof(REAL));
  a_bar=(REAL*)calloc(np+nb,sizeof(REAL));

  ComputeDerivative(np,nb,x,&Energy,a,&StrainFirstDerivative);


  for(i=0;i<np+nb;i++)
    v[i]=-DeltaT*a[i];

  for(loop=0;loop<MaximumNumberOfMinimizationSteps*1000;loop++)
  {
    if(Movies[CurrentSystem]&&(loop%WriteMoviesEvery[CurrentSystem]==0))
      SamplePDBMovies(SAMPLE,run);

    for(i=0;i<np+nb;i++)
      x_new[i]=x[i]+DeltaT*v[i];

    ComputeDerivative(np,nb,x_new,&Energy,a_new,&StrainFirstDerivative);

    // Find the RMS gradient.
    RMSGradient=0.0;
    MaxGradient=0.0;
    for(i=0;i<np+nb;i++)
    {
      if(fabs(a_new[i])>MaxGradient) MaxGradient=fabs(a_new[i]);
      RMSGradient+=SQR(a_new[i]);
    }

    RMSGradient=sqrt(RMSGradient)/(REAL)(np+nb);

    // Check for convergence.
    if(((RMSGradient<RMSGradientTolerance)&&(MaxGradient<MaxGradientTolerance)))
    {
      fprintf(OutputFilePtr[CurrentSystem],"SUCCES: RMS Gradient tolerance %g reached (%g)\n",(double)RMSGradientTolerance,(double)RMSGradient);
      fprintf(OutputFilePtr[CurrentSystem],"        Max Gradient tolerance %g reached (%g)\n",(double)MaxGradientTolerance,(double)MaxGradient);
      break;
    }

    if(loop%PrintEvery==0)
    {
      switch(Dimension)
      {
        case 2:
          fprintf(OutputFilePtr[CurrentSystem],"Iteration: %d Energy: %18.10f Area: %18.10f RMS gradient: %g MaxGradient: %g\n",
            loop,(double)(Energy*ENERGY_TO_KELVIN),(double)Volume[CurrentSystem],(double)RMSGradient,(double)MaxGradient);
          fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f     Strain derivative: %18.10f %18.10f\n",
              (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,
              (double)StrainFirstDerivative.ax,(double)StrainFirstDerivative.bx);
          fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f                        %18.10f %18.10f\n",
              (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,
              (double)StrainFirstDerivative.ay,(double)StrainFirstDerivative.by);
          fprintf(OutputFilePtr[CurrentSystem],"      Lengths: %18.10f %18.10f, Angle: %18.10f\n",
             (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,GammaAngle[CurrentSystem]*RAD2DEG);
          break;
        case 3:
          fprintf(OutputFilePtr[CurrentSystem],"Iteration: %d Energy: %18.10f Volume: %18.10f RMS gradient: %g MaxGradient: %g\n",
            loop,(double)(Energy*ENERGY_TO_KELVIN),(double)Volume[CurrentSystem],(double)RMSGradient,(double)MaxGradient);
          fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f %18.10f     Strain derivative: %18.10f %18.10f %18.10f\n",
              (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx,
              (double)StrainFirstDerivative.ax,(double)StrainFirstDerivative.bx,(double)StrainFirstDerivative.cx);
          fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
              (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy,
              (double)StrainFirstDerivative.ay,(double)StrainFirstDerivative.by,(double)StrainFirstDerivative.cy);
          fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
              (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz,
              (double)StrainFirstDerivative.az,(double)StrainFirstDerivative.bz,(double)StrainFirstDerivative.cz);
          fprintf(OutputFilePtr[CurrentSystem],"      Lengths: %18.10f %18.10f %18.10f, Angles: %18.10f %18.10f %18.10f\n",
             (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,(double)BoxProperties[CurrentSystem].az,
             AlphaAngle[CurrentSystem]*RAD2DEG,BetaAngle[CurrentSystem]*RAD2DEG,GammaAngle[CurrentSystem]*RAD2DEG);
          break;
      }
      fflush(OutputFilePtr[CurrentSystem]);
    }


    DotProduct=0.0;
    for(i=0;i<np+nb;i++)
      DotProduct+=-v[i]*a_new[i];

    if(DotProduct<=0)
    {
      for(i=0;i<np+nb;i++)
        v_bar[i]=-DeltaT*a[i];

      for(i=0;i<np+nb;i++)
        x_bar[i]=x[i]+DeltaT*v_bar[i];

      ComputeDerivative(np,nb,x_bar,&Energy,a_bar,&StrainFirstDerivative);

      DotProduct=0.0;
      for(i=0;i<np+nb;i++)
        DotProduct+=-v_bar[i]*a_bar[i];

      if(DotProduct<=0.0)
      {
        DeltaT=0.5*DeltaT;
        for(i=0;i<np+nb;i++)
          v_new[i]=v[i];
        for(i=0;i<np+nb;i++)
          x_new[i]=x[i];
      }
      else
      {
        for(i=0;i<np+nb;i++)
          x_new[i]=x[i];
        for(i=0;i<np+nb;i++)
          v[i]=0.5*(v[i]-v_bar[i]);
        for(i=0;i<np+nb;i++)
          v_new[i]=v[i]-DeltaT*a[i];

      }
    }
    else
    {
      for(i=0;i<np+nb;i++)
        v_new[i]=v[i]-DeltaT*a_new[i];
    }

    for(i=0;i<np+nb;i++)
    {
      x[i]=x_new[i];
      v[i]=v_new[i];
      a[i]=a_new[i];
    }
  }

  SamplePDBMovies(FINALIZE,run);

  EigenValues=(REAL*)calloc(NumberOfMinimizationVariables,sizeof(REAL));
  HessianMatrix=CreateRealMatrix(NumberOfMinimizationVariables,NumberOfMinimizationVariables);

  ComputeDerivativesMinimization(np,nb,x_new,&Energy,a_new,HessianMatrix,&StrainFirstDerivative);

  // project the constraints from the gradient and Hessian
  ProjectConstraintsFromHessianMatrix(np,nb,a_new,HessianMatrix,TRUE,TRUE);

  if(MinimizationVariables==FRACTIONAL)
    ConvertHessianFromCartesianToFractional(np,nb,a_new,HessianMatrix);

  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
  fprintf(OutputFilePtr[CurrentSystem],"Final energy after minimization: %18.10f\n",(double)(Energy*ENERGY_TO_KELVIN));
  Pressure=-(StrainFirstDerivative.ax+StrainFirstDerivative.by+StrainFirstDerivative.cz)/(Dimension*Volume[CurrentSystem]);
  fprintf(OutputFilePtr[CurrentSystem],"Final pressure after minimization: %18.10f [Pa]\n",(double)(Pressure*PRESSURE_CONVERSION_FACTOR));
  switch(Dimension)
  {
    case 3:
     fprintf(OutputFilePtr[CurrentSystem],"Final stress after minimization: {{%f,%f,%f},{-,%f,%f},{-,-,%f}} %s\n",
          (double)(-StrainFirstDerivative.ax*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.bx*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.cx*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.by*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.bz*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.cz*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          UseReducedUnits?"[-]":"[Pa]");
      break;
    case 2:
     fprintf(OutputFilePtr[CurrentSystem],"Final stress after minimization: {{%f,%f},{-,%f}} %s\n",
          (double)(-StrainFirstDerivative.ax*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.bx*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.by*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          UseReducedUnits?"[-]":"[Pa]");
      break;
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
  fprintf(OutputFilePtr[CurrentSystem],"Variable vector:\n");
  fprintf(OutputFilePtr[CurrentSystem],"================\n");
  for(i=0;i<NumberOfMinimizationVariables;i++)
    fprintf(OutputFilePtr[CurrentSystem],"x: %d = %g\n",i,(double)x[i]);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
  fprintf(OutputFilePtr[CurrentSystem],"Gradient vector:\n");
  fprintf(OutputFilePtr[CurrentSystem],"================\n");
  for(i=0;i<NumberOfMinimizationVariables;i++)
    fprintf(OutputFilePtr[CurrentSystem],"Gradient: %d = %g\n",i,(double)a_new[i]);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  SolveEigenValuesAndVectorsHessian(HessianMatrix,EigenValues);

  fprintf(OutputFilePtr[CurrentSystem],"Eigenvalues:\n");
  fprintf(OutputFilePtr[CurrentSystem],"============\n");
  for(i=0;i<NumberOfMinimizationVariables;i++)
    fprintf(OutputFilePtr[CurrentSystem],"Eigenvalue: %d = %18.10f\n",i,(double)EigenValues[i]);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");


  ComputeBornTerm=FALSE;
  CalculateForce();

  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
  fprintf(OutputFilePtr[CurrentSystem],"Value from force routines:\n");
  fprintf(OutputFilePtr[CurrentSystem],"===========================\n");
  fprintf(OutputFilePtr[CurrentSystem],"Final energy after minimization: %18.10f\n",(double)(UTotal[CurrentSystem]*ENERGY_TO_KELVIN));
  Pressure=-(StrainDerivativeTensor[CurrentSystem].ax+StrainDerivativeTensor[CurrentSystem].by+StrainDerivativeTensor[CurrentSystem].cz)/(Dimension*Volume[CurrentSystem]);
  fprintf(OutputFilePtr[CurrentSystem],"Final pressure after minimization: %18.10f [Pa]\n",(double)(Pressure*PRESSURE_CONVERSION_FACTOR));
  switch(Dimension)
  {
    case 3:
     fprintf(OutputFilePtr[CurrentSystem],"Final stress after minimization: {{%f,%f,%f},{-,%f,%f},{-,-,%f}} %s\n",
          (double)(-StrainFirstDerivative.ax*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.bx*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.cx*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.by*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.bz*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.cz*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          UseReducedUnits?"[-]":"[Pa]");
      break;
    case 2:
     fprintf(OutputFilePtr[CurrentSystem],"Final stress after minimization: {{%f,%f},{-,%f}} %s\n",
          (double)(-StrainFirstDerivative.ax*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.bx*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          (double)(-StrainFirstDerivative.by*PRESSURE_CONVERSION_FACTOR/Volume[CurrentSystem]),
          UseReducedUnits?"[-]":"[Pa]");
      break;
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  fprintf(OutputFilePtr[CurrentSystem],"      Box: %18.10f %18.10f %18.10f     Strain derivative: %18.10f %18.10f %18.10f\n",
        (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx,
        (double)StrainDerivativeTensor[CurrentSystem].ax,(double)StrainDerivativeTensor[CurrentSystem].bx,(double)StrainDerivativeTensor[CurrentSystem].cx);
  fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
        (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy,
        (double)StrainDerivativeTensor[CurrentSystem].ay,(double)StrainDerivativeTensor[CurrentSystem].by,(double)StrainDerivativeTensor[CurrentSystem].cy);
  fprintf(OutputFilePtr[CurrentSystem],"           %18.10f %18.10f %18.10f                        %18.10f %18.10f %18.10f\n",
        (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz,
        (double)StrainDerivativeTensor[CurrentSystem].az,(double)StrainDerivativeTensor[CurrentSystem].bz,(double)StrainDerivativeTensor[CurrentSystem].cz);
  fprintf(OutputFilePtr[CurrentSystem],"      Lengths: %18.10f %18.10f %18.10f, Angles: %18.10f %18.10f %18.10f\n",
       (double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,(double)BoxProperties[CurrentSystem].az,
       AlphaAngle[CurrentSystem]*RAD2DEG,BetaAngle[CurrentSystem]*RAD2DEG,GammaAngle[CurrentSystem]*RAD2DEG);
  fprintf(OutputFilePtr[CurrentSystem],"\n");


  fprintf(OutputFilePtr[CurrentSystem],"Forces\n");
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      fprintf(OutputFilePtr[CurrentSystem],"Framework[%d] Atom: %d  %g %g %g\n",f1,i,
        Framework[CurrentSystem].Atoms[f1][i].Force.x,
        Framework[CurrentSystem].Atoms[f1][i].Force.y,
        Framework[CurrentSystem].Atoms[f1][i].Force.z);
  }
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Adsorbate[%d] Atom: %d  %g %g %g\n",i,j,
        Adsorbates[CurrentSystem][i].Atoms[j].Force.x,
        Adsorbates[CurrentSystem][i].Atoms[j].Force.y,
        Adsorbates[CurrentSystem][i].Atoms[j].Force.z);
    }
  }
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      fprintf(OutputFilePtr[CurrentSystem],"Cations[%d] Atom: %d  %g %g %g\n",i,j,
        Cations[CurrentSystem][i].Atoms[j].Force.x,
        Cations[CurrentSystem][i].Atoms[j].Force.y,
        Cations[CurrentSystem][i].Atoms[j].Force.z);
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  ComputeElasticConstantsAfterMinimization(np,nb,x);

  CurrentSystem=0;
  PrintRestartFile();
}

void AllocateMinimizationLocalMemory(void)
{
  int k,m;
  int total,max;

  max=0;
  for(k=0;k<NumberOfSystems;k++)
  {
    total=0;
    for(m=0;m<NumberOfAdsorbateMolecules[k];m++)
      total+=Components[Adsorbates[k][m].Type].NumberOfRigidAtoms;
    for(m=0;m<NumberOfCationMolecules[k];m++)
      total+=Components[Cations[k][m].Type].NumberOfRigidAtoms;
    if(total>max) max=total;
  }

  DVecX=(VECTOR*)calloc(max,sizeof(VECTOR));
  DVecY=(VECTOR*)calloc(max,sizeof(VECTOR));
  DVecZ=(VECTOR*)calloc(max,sizeof(VECTOR));

  DDVecAX=(VECTOR*)calloc(max,sizeof(VECTOR));
  DDVecBY=(VECTOR*)calloc(max,sizeof(VECTOR));
  DDVecCZ=(VECTOR*)calloc(max,sizeof(VECTOR));
  DDVecAY=(VECTOR*)calloc(max,sizeof(VECTOR));
  DDVecAZ=(VECTOR*)calloc(max,sizeof(VECTOR));
  DDVecBZ=(VECTOR*)calloc(max,sizeof(VECTOR));
}

void FreeMinimizationLocalMemory(void)
{
  free(DVecX);
  free(DVecY);
  free(DVecZ);

  free(DDVecAX);
  free(DDVecBY);
  free(DDVecCZ);
  free(DDVecAY);
  free(DDVecAZ);
  free(DDVecBZ);
}


void AllocateMinimizationMemory(void)
{
  FrameworkFixedInitialization=(int**)calloc(NumberOfSystems,sizeof(int*));
  AdsorbateFixedInitialization=(int*)calloc(NumberOfSystems,sizeof(int));
  CationFixedInitialization=(int*)calloc(NumberOfSystems,sizeof(int));

  NumberOfRattleCyclesStage1=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  MaximumNumberOfRattleCyclesStage1=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfRattleCyclesStage2=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  MaximumNumberOfRattleCyclesStage2=(int*)calloc(NumberOfSystems,sizeof(int));

  NumberOfFixedFrameworkAtoms=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfActiveFrameworkAtoms=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfFixedFrameworkAtomsX=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfActiveFrameworkAtomsX=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfFixedFrameworkAtomsY=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfActiveFrameworkAtomsY=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfFixedFrameworkAtomsZ=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfActiveFrameworkAtomsZ=(int**)calloc(NumberOfSystems,sizeof(int*));

  FixedFrameworkAtoms=(int***)calloc(NumberOfSystems,sizeof(int**));
  ActiveFrameworkAtoms=(int***)calloc(NumberOfSystems,sizeof(int**));
  FixedFrameworkAtomsX=(int***)calloc(NumberOfSystems,sizeof(int**));
  ActiveFrameworkAtomsX=(int***)calloc(NumberOfSystems,sizeof(int**));
  FixedFrameworkAtomsY=(int***)calloc(NumberOfSystems,sizeof(int**));
  ActiveFrameworkAtomsY=(int***)calloc(NumberOfSystems,sizeof(int**));
  FixedFrameworkAtomsZ=(int***)calloc(NumberOfSystems,sizeof(int**));
  ActiveFrameworkAtomsZ=(int***)calloc(NumberOfSystems,sizeof(int**));

  NumberOfFixedAdsorbateMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateMolecules=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfActiveAdsorbateMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateMolecules=(int**)calloc(NumberOfSystems,sizeof(int*));

  NumberOfFixedAdsorbateAtoms=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateAtoms=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedAdsorbateAtomsX=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateAtomsX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedAdsorbateAtomsY=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateAtomsY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedAdsorbateAtomsZ=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateAtomsZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateAtoms=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateAtoms=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateAtomsX=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateAtomsX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateAtomsY=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateAtomsY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateAtomsZ=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateAtomsZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfFixedAdsorbateGroups=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateGroups=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateGroups=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateGroups=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfFixedAdsorbateGroupsCenterOfMass=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateGroupsCenterOfMass=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedAdsorbateGroupsCenterOfMassX=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateGroupsCenterOfMassX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedAdsorbateGroupsCenterOfMassY=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateGroupsCenterOfMassY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedAdsorbateGroupsCenterOfMassZ=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateGroupsCenterOfMassZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfActiveAdsorbateGroupsCenterOfMass=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateGroupsCenterOfMass=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateGroupsCenterOfMassX=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateGroupsCenterOfMassX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateGroupsCenterOfMassY=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateGroupsCenterOfMassY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateGroupsCenterOfMassZ=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateGroupsCenterOfMassZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfFixedAdsorbateGroupsOrientation=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateGroupsOrientation=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedAdsorbateGroupsOrientationX=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateGroupsOrientationX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedAdsorbateGroupsOrientationY=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateGroupsOrientationY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedAdsorbateGroupsOrientationZ=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedAdsorbateGroupsOrientationZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfActiveAdsorbateGroupsOrientation=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateGroupsOrientation=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateGroupsOrientationX=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateGroupsOrientationX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateGroupsOrientationY=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateGroupsOrientationY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveAdsorbateGroupsOrientationZ=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveAdsorbateGroupsOrientationZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfFixedCationMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationMolecules=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfActiveCationMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationMolecules=(int**)calloc(NumberOfSystems,sizeof(int*));

  NumberOfFixedCationAtoms=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationAtoms=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedCationAtomsX=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationAtomsX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedCationAtomsY=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationAtomsY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedCationAtomsZ=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationAtomsZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationAtoms=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationAtoms=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationAtomsX=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationAtomsX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationAtomsY=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationAtomsY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationAtomsZ=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationAtomsZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfFixedCationGroups=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationGroups=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationGroups=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationGroups=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfFixedCationGroupsCenterOfMass=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationGroupsCenterOfMass=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedCationGroupsCenterOfMassX=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationGroupsCenterOfMassX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedCationGroupsCenterOfMassY=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationGroupsCenterOfMassY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedCationGroupsCenterOfMassZ=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationGroupsCenterOfMassZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfActiveCationGroupsCenterOfMass=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationGroupsCenterOfMass=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationGroupsCenterOfMassX=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationGroupsCenterOfMassX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationGroupsCenterOfMassY=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationGroupsCenterOfMassY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationGroupsCenterOfMassZ=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationGroupsCenterOfMassZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfFixedCationGroupsOrientation=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationGroupsOrientation=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedCationGroupsOrientationX=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationGroupsOrientationX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedCationGroupsOrientationY=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationGroupsOrientationY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfFixedCationGroupsOrientationZ=(int*)calloc(NumberOfSystems,sizeof(int));
  FixedCationGroupsOrientationZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfActiveCationGroupsOrientation=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationGroupsOrientation=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationGroupsOrientationX=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationGroupsOrientationX=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationGroupsOrientationY=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationGroupsOrientationY=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));
  NumberOfActiveCationGroupsOrientationZ=(int*)calloc(NumberOfSystems,sizeof(int));
  ActiveCationGroupsOrientationZ=(PAIR**)calloc(NumberOfSystems,sizeof(PAIR*));

  NumberOfDistanceConstraints=(int*)calloc(NumberOfSystems,sizeof(int));
  DistanceConstraints=(ATOM*(**)[2])calloc(NumberOfSystems,sizeof(ATOM*(*)[2]));
  DistanceDefinitions=(int (**)[2][3])calloc(NumberOfSystems,sizeof(int (*)[2][3]));
  DistanceConstraintsDerivatives=(VECTOR (**)[2])calloc(NumberOfSystems,sizeof(VECTOR (*)[2]));

  NumberOfAngleConstraints=(int*)calloc(NumberOfSystems,sizeof(int));
  AngleConstraints=(ATOM*(**)[3])calloc(NumberOfSystems,sizeof(ATOM*(*)[3]));
  AngleDefinitions=(int (**)[3][3])calloc(NumberOfSystems,sizeof(int (*)[3][3]));
  AngleConstraintsDerivatives=(VECTOR (**)[3])calloc(NumberOfSystems,sizeof(VECTOR (*)[3]));

  NumberOfDihedralConstraints=(int*)calloc(NumberOfSystems,sizeof(int));
  DihedralConstraints=(ATOM*(**)[4])calloc(NumberOfSystems,sizeof(ATOM*(*)[4]));
  DihedralDefinitions=(int (**)[4][3])calloc(NumberOfSystems,sizeof(int (*)[4][3]));
  DihedralConstraintsDerivatives=(VECTOR (**)[4])calloc(NumberOfSystems,sizeof(VECTOR (*)[4]));

  NumberOfImproperDihedralConstraints=(int*)calloc(NumberOfSystems,sizeof(int));
  ImproperDihedralConstraints=(ATOM*(**)[4])calloc(NumberOfSystems,sizeof(ATOM*(*)[4]));
  ImproperDihedralDefinitions=(int (**)[4][3])calloc(NumberOfSystems,sizeof(int (*)[4][3]));
  ImproperDihedralConstraintsDerivatives=(VECTOR (**)[4])calloc(NumberOfSystems,sizeof(VECTOR (*)[4]));

  NumberOfInversionBendConstraints=(int*)calloc(NumberOfSystems,sizeof(int));
  InversionBendConstraints=(ATOM*(**)[4])calloc(NumberOfSystems,sizeof(ATOM*(*)[4]));
  InversionBendDefinitions=(int (**)[4][3])calloc(NumberOfSystems,sizeof(int (*)[4][3]));
  InversionBendConstraintsDerivatives=(VECTOR (**)[4])calloc(NumberOfSystems,sizeof(VECTOR (*)[4]));

  NumberOfOutOfPlaneDistanceConstraints=(int*)calloc(NumberOfSystems,sizeof(int));
  OutOfPlaneDistanceConstraints=(ATOM*(**)[4])calloc(NumberOfSystems,sizeof(ATOM*(*)[4]));
  OutOfPlaneDistanceDefinitions=(int (**)[4][3])calloc(NumberOfSystems,sizeof(int (*)[4][3]));
  OutOfPlaneDistanceConstraintsDerivatives=(VECTOR (**)[4])calloc(NumberOfSystems,sizeof(VECTOR (*)[4]));

  NumberOfHarmonicDistanceConstraints=(int*)calloc(NumberOfSystems,sizeof(int));
  HarmonicDistanceConstraints=(ATOM*(**)[2])calloc(NumberOfSystems,sizeof(ATOM*(*)[2]));
  HarmonicDistanceConstraintParameters=(REAL (**)[2])calloc(NumberOfSystems,sizeof(REAL (*)[2]));
  HarmonicDistanceDefinitions=(int (**)[2][3])calloc(NumberOfSystems,sizeof(int (*)[2][3]));

  NumberOfHarmonicAngleConstraints=(int*)calloc(NumberOfSystems,sizeof(int));
  HarmonicAngleConstraints=(ATOM*(**)[3])calloc(NumberOfSystems,sizeof(ATOM*(*)[3]));
  HarmonicAngleConstraintParameters=(REAL (**)[2])calloc(NumberOfSystems,sizeof(REAL (*)[2]));
  HarmonicAngleDefinitions=(int (**)[3][3])calloc(NumberOfSystems,sizeof(int (*)[3][3]));

  NumberOfHarmonicDihedralConstraints=(int*)calloc(NumberOfSystems,sizeof(int));
  HarmonicDihedralConstraints=(ATOM*(**)[4])calloc(NumberOfSystems,sizeof(ATOM*(*)[4]));
  HarmonicDihedralConstraintParameters=(REAL (**)[2])calloc(NumberOfSystems,sizeof(REAL (*)[2]));
  HarmonicDihedralDefinitions=(int (**)[4][3])calloc(NumberOfSystems,sizeof(int (*)[4][3]));

  DistanceConstraintParameter=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  AngleConstraintParameter=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  DihedralConstraintParameter=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ImproperDihedralConstraintParameter=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  InversionBendConstraintParameter=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  OutOfPlaneDistanceConstraintParameter=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  NumberOfTwoPointDihedralDefinitions=(int*)calloc(NumberOfSystems,sizeof(int));
  TwoPointDihedralDefinitions=(int (**)[6][3])calloc(NumberOfSystems,sizeof(int (*)[6][3]));
  TwoPointDihedrals=(ATOM*(**)[6])calloc(NumberOfSystems,sizeof(ATOM*(*)[6]));
}

static int versionNumber=1;

void WriteRestartMinimization(FILE *FilePtr)
{
  REAL Check;

  fwrite(&versionNumber,sizeof(int),1,FilePtr);

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);


}

void ReadRestartMinimization(FILE *FilePtr)
{
  REAL Check;
  int readversionNumber=0;

  AllocateMinimizationMemory();

  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumber)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }


  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartMinimization)\n");
    ContinueAfterCrash=FALSE;
  }
}

