/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'minimization.h' is part of RASPA-2.0

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

#ifndef MINIMIZATION_H
#define MINIMIZATION_H

#define MINIMUM_EIGEN_VALUE 0.001

extern int UseSymmetryInMinimization;

extern int **FrameworkFixedInitialization;
extern int *AdsorbateFixedInitialization;
extern int *CationFixedInitialization;

extern int **NumberOfFixedFrameworkAtoms;
extern int ***FixedFrameworkAtoms;
extern int **NumberOfFixedFrameworkAtomsX;
extern int ***FixedFrameworkAtomsX;
extern int **NumberOfFixedFrameworkAtomsY;
extern int ***FixedFrameworkAtomsY;
extern int **NumberOfFixedFrameworkAtomsZ;
extern int ***FixedFrameworkAtomsZ;

extern int **NumberOfActiveFrameworkAtoms;
extern int ***ActiveFrameworkAtoms;
extern int **NumberOfActiveFrameworkAtomsX;
extern int ***ActiveFrameworkAtomsX;
extern int **NumberOfActiveFrameworkAtomsY;
extern int ***ActiveFrameworkAtomsY;
extern int **NumberOfActiveFrameworkAtomsZ;
extern int ***ActiveFrameworkAtomsZ;

extern int *NumberOfFixedAdsorbateMolecules;
extern int **FixedAdsorbateMolecules;
extern int *NumberOfActiveAdsorbateMolecules;
extern int **ActiveAdsorbateMolecules;

extern int *NumberOfFixedAdsorbateAtoms;
extern PAIR **FixedAdsorbateAtoms;
extern int *NumberOfFixedAdsorbateAtomsX;
extern PAIR **FixedAdsorbateAtomsX;
extern int *NumberOfFixedAdsorbateAtomsY;
extern PAIR **FixedAdsorbateAtomsY;
extern int *NumberOfFixedAdsorbateAtomsZ;
extern PAIR **FixedAdsorbateAtomsZ;
extern int *NumberOfActiveAdsorbateAtoms;
extern PAIR **ActiveAdsorbateAtoms;
extern int *NumberOfActiveAdsorbateAtomsX;
extern PAIR **ActiveAdsorbateAtomsX;
extern int *NumberOfActiveAdsorbateAtomsX;
extern PAIR **ActiveAdsorbateAtomsX;
extern int *NumberOfActiveAdsorbateAtomsY;
extern PAIR **ActiveAdsorbateAtomsY;
extern int *NumberOfActiveAdsorbateAtomsZ;
extern PAIR **ActiveAdsorbateAtomsZ;

extern int *NumberOfFixedAdsorbateGroups;
extern PAIR **FixedAdsorbateGroups;
extern int *NumberOfActiveAdsorbateGroups;
extern PAIR **ActiveAdsorbateGroups;

extern int *NumberOfFixedAdsorbateGroupsCenterOfMass;
extern PAIR **FixedAdsorbateGroupsCenterOfMass;
extern int *NumberOfFixedAdsorbateGroupsCenterOfMassX;
extern PAIR **FixedAdsorbateGroupsCenterOfMassX;
extern int *NumberOfFixedAdsorbateGroupsCenterOfMassY;
extern PAIR **FixedAdsorbateGroupsCenterOfMassY;
extern int *NumberOfFixedAdsorbateGroupsCenterOfMassZ;
extern PAIR **FixedAdsorbateGroupsCenterOfMassZ;
extern int *NumberOfActiveAdsorbateGroupsCenterOfMass;
extern PAIR **ActiveAdsorbateGroupsCenterOfMass;

extern int *NumberOfFixedAdsorbateGroupsOrientation;
extern PAIR **FixedAdsorbateGroupsOrientation;
extern int *NumberOfFixedAdsorbateGroupsOrientationX;
extern PAIR **FixedAdsorbateGroupsOrientationX;
extern int *NumberOfFixedAdsorbateGroupsOrientationY;
extern PAIR **FixedAdsorbateGroupsOrientationY;
extern int *NumberOfFixedAdsorbateGroupsOrientationZ;
extern PAIR **FixedAdsorbateGroupsOrientationZ;
extern int *NumberOfActiveAdsorbateGroupsOrientation;
extern PAIR **ActiveAdsorbateGroupsOrientation;

extern int *NumberOfFixedCationMolecules;
extern int **FixedCationMolecules;
extern int *NumberOfActiveCationMolecules;
extern int **ActiveCationMolecules;

extern int *NumberOfFixedCationAtoms;
extern PAIR **FixedCationAtoms;
extern int *NumberOfFixedCationAtomsX;
extern PAIR **FixedCationAtomsX;
extern int *NumberOfFixedCationAtomsY;
extern PAIR **FixedCationAtomsY;
extern int *NumberOfFixedCationAtomsZ;
extern PAIR **FixedCationAtomsZ;
extern int *NumberOfActiveCationAtoms;
extern PAIR **ActiveCationAtoms;
extern int *NumberOfActiveCationAtomsX;
extern PAIR **ActiveCationAtomsX;
extern int *NumberOfActiveCationAtomsY;
extern PAIR **ActiveCationAtomsY;
extern int *NumberOfActiveCationAtomsZ;
extern PAIR **ActiveCationAtomsZ;

extern int *NumberOfFixedCationGroups;
extern PAIR **FixedCationGroups;
extern int *NumberOfActiveCationGroups;
extern PAIR **ActiveCationGroups;

extern int *NumberOfFixedCationGroupsCenterOfMass;
extern PAIR **FixedCationGroupsCenterOfMass;
extern int *NumberOfFixedCationGroupsCenterOfMassX;
extern PAIR **FixedCationGroupsCenterOfMassX;
extern int *NumberOfFixedCationGroupsCenterOfMassY;
extern PAIR **FixedCationGroupsCenterOfMassY;
extern int *NumberOfFixedCationGroupsCenterOfMassZ;
extern PAIR **FixedCationGroupsCenterOfMassZ;

extern int *NumberOfActiveCationGroupsCenterOfMass;
extern PAIR **ActiveCationGroupsCenterOfMass;
extern int *NumberOfActiveCationGroupsCenterOfMassX;
extern PAIR **ActiveCationGroupsCenterOfMassX;
extern int *NumberOfActiveCationGroupsCenterOfMassY;
extern PAIR **ActiveCationGroupsCenterOfMassY;
extern int *NumberOfActiveCationGroupsCenterOfMassZ;
extern PAIR **ActiveCationGroupsCenterOfMassZ;

extern int *NumberOfFixedCationGroupsOrientation;
extern PAIR **FixedCationGroupsOrientation;
extern int *NumberOfFixedCationGroupsOrientationX;
extern PAIR **FixedCationGroupsOrientationX;
extern int *NumberOfFixedCationGroupsOrientationY;
extern PAIR **FixedCationGroupsOrientationY;
extern int *NumberOfFixedCationGroupsOrientationZ;
extern PAIR **FixedCationGroupsOrientationZ;

extern int *NumberOfActiveCationGroupsOrientation;
extern PAIR **ActiveCationGroupsOrientation;
extern int *NumberOfActiveCationGroupsOrientationX;
extern PAIR **ActiveCationGroupsOrientationX;
extern int *NumberOfActiveCationGroupsOrientationY;
extern PAIR **ActiveCationGroupsOrientationY;
extern int *NumberOfActiveCationGroupsOrientationZ;
extern PAIR **ActiveCationGroupsOrientationZ;


extern int NumberOfFixedAtomTypes;
extern char (*FixedAtomTypes)[256];
extern int NumberOfActiveAtomTypes;
extern char (*ActiveAtomTypes)[256];


// distance constraints between two arbitrary atoms
extern int *NumberOfDistanceConstraints;
extern ATOM* (**DistanceConstraints)[2];
extern REAL **DistanceConstraintParameter;
extern int (**DistanceDefinitions)[2][3];
extern VECTOR (**DistanceConstraintsDerivatives)[2];

// angular constraints between three arbitrary atoms
extern int *NumberOfAngleConstraints;
extern ATOM* (**AngleConstraints)[3];
extern REAL **AngleConstraintParameter;
extern int (**AngleDefinitions)[3][3];
extern VECTOR (**AngleConstraintsDerivatives)[3];

// dihedral constraints between four arbitrary atoms
extern int *NumberOfDihedralConstraints;
extern ATOM* (**DihedralConstraints)[4];
extern REAL **DihedralConstraintParameter;
extern int (**DihedralDefinitions)[4][3];
extern VECTOR (**DihedralConstraintsDerivatives)[4];

// improper dihedral constraints between four arbitrary atoms
extern int *NumberOfImproperDihedralConstraints;
extern ATOM* (**ImproperDihedralConstraints)[4];
extern REAL **ImproperDihedralConstraintParameter;
extern int (**ImproperDihedralDefinitions)[4][3];
extern VECTOR (**ImproperDihedralConstraintsDerivatives)[4];


// Inversion-bend constraints between four arbitrary atoms
extern int *NumberOfInversionBendConstraints;
extern ATOM* (**InversionBendConstraints)[4];
extern REAL **InversionBendConstraintParameter;
extern int (**InversionBendDefinitions)[4][3];
extern VECTOR (**InversionBendConstraintsDerivatives)[4];

// Out-of-plane-distance constraints between four arbitrary atoms
extern int *NumberOfOutOfPlaneDistanceConstraints;
extern ATOM* (**OutOfPlaneDistanceConstraints)[4];
extern REAL **OutOfPlaneDistanceConstraintParameter;
extern int (**OutOfPlaneDistanceDefinitions)[4][3];
extern VECTOR (**OutOfPlaneDistanceConstraintsDerivatives)[4];


// harmonic distance constraints between two arbitrary atoms
extern int *NumberOfHarmonicDistanceConstraints;
extern ATOM* (**HarmonicDistanceConstraints)[2];
extern REAL (**HarmonicDistanceConstraintParameters)[2];
extern int (**HarmonicDistanceDefinitions)[2][3];

// harmonic angular constraints between three arbitrary atoms
extern int *NumberOfHarmonicAngleConstraints;
extern ATOM* (**HarmonicAngleConstraints)[3];
extern REAL (**HarmonicAngleConstraintParameters)[2];
extern int (**HarmonicAngleDefinitions)[3][3];

// harmonic dihedral constraints between four arbitrary atoms
extern int *NumberOfHarmonicDihedralConstraints;
extern ATOM* (**HarmonicDihedralConstraints)[4];
extern REAL (**HarmonicDihedralConstraintParameters)[2];
extern int (**HarmonicDihedralDefinitions)[4][3];

// two point dihedral measuments between four arbitrary points, the first point is the mid-point between atom 1 and 2
extern int *NumberOfTwoPointDihedralDefinitions;
extern int (**TwoPointDihedralDefinitions)[6][3];
extern ATOM* (**TwoPointDihedrals)[6];

enum{STEEPEST_DESCENT_MINIMIZATION,CONJUGATE_GRADIENT_MINIMIZATION,BFGS_MINIMIZATION,SNYMAN_MINIMIZATION,BAKER_MINIMIZATION,NEWTON_RAPHSON_MINIMIZATION,BAKER_SADDLE_POINT};
enum{SEARCH_LOCAL_SADDLE_POINT,SEARCH_LOCAL_MINIMUM};
enum{ANALYTICALLY,NUMERICALLY};
enum{LAMBDA_METHOD_1,LAMBDA_METHOD_2,LAMBDA_METHOD_3};
enum{CARTESIAN,FRACTIONAL};

extern int RemoveTranslationFromHessian;
extern int RemoveRotationFromHessian;

extern int MinimizationMethod;
extern int MinimizationPotentialMethod;
extern int MinimizationType;
extern int MinimizationReferenceParticle;
extern int ComputeLambda;
extern int MinimizationVariables;

extern int ShellIndex;
extern int ShellSize;
extern int CoreSize;

extern int NumberOfMinimizationVariables;                // the total number of variables used in the minimization
extern int NumberOfCellMinimizationVariables;            // the number of cell-variables used in the minimization (0-6)
extern int NumberOfCoordinatesMinimizationVariables;     // the total number of variables used for generalized coordinates in the minimization
extern int NumberOfPositionalMinimizationVariables;      // the total number of variables used for positions in the minimization
extern int NumberOfOrientationalMinimizationVariables;   // the total number of variables used for orientation in the minimization

extern REAL MaximumStepLengthInput;
extern REAL MaximumStepLength;
extern REAL MinimizationConvergenceFactor;

extern int MaximumNumberOfMinimizationSteps;
extern int UseGradientInLineMinimization;

extern REAL RMSGradientTolerance;
extern REAL MaxGradientTolerance;

extern REAL_MATRIX3x3 StoredBox;
extern REAL_MATRIX3x3 StoredReplicaBox;
extern REAL_MATRIX3x3 StoredInverseBox;

extern VECTOR *DVecX;
extern VECTOR *DVecY;
extern VECTOR *DVecZ;
extern VECTOR *DDVecAX;
extern VECTOR *DDVecBY;
extern VECTOR *DDVecCZ;
extern VECTOR *DDVecAY;
extern VECTOR *DDVecAZ;
extern VECTOR *DDVecBZ;

extern REAL UnitCellDeformation;
extern int TransformUnitCell;
extern int ElasticConstantEnergyVolume;
enum {ELASTIC_CONSTANT_CUBIC_C44,ELASTIC_CONSTANT_CUBIC_CS,ELASTIC_CONSTANT_CUBIC_C66,ELASTIC_CONSTANT_HEXAGONAL_C44,ELASTIC_CONSTANT_BULK_MODULUS,
      ELASTIC_CONSTANT_ORTHO_C11,ELASTIC_CONSTANT_ORTHO_C12,ELASTIC_CONSTANT_ORTHO_C13};

void EvaluateDerivatives(int n,REAL *Energy,REAL* Gradient,REAL_MATRIX Hessian,REAL_MATRIX3x3 *StrainFirstDerivative,
                         int ComputeGradient,int ComputeHessian);

void CreateGeneralizedCoordinatesFromPositions(int np,int nb,REAL *x);
void CreatePositionsFromGeneralizedCoordinates(int np,int nb,REAL *x);

void ConvertGradientFromCartesianToFractional(REAL *Gradient);
void ConvertHessianFromCartesianToFractional(int np,int nb,REAL *Gradient,REAL_MATRIX Hessian);

void SetWeights(int np,REAL *Weights,REAL *Charges);
void SetStrainToZero(int np,REAL *x);
int OrderNumberOfMinimiationVariables(void);

void AdjustToUpperTriangle(void);

void ComputeDerivativesMinimization(int np,int nb,REAL *x,REAL* Energy,REAL *Gradient,REAL_MATRIX Hessian,REAL_MATRIX3x3 *StrainFirstDerivative);
void ComputeDerivativesSpectra(int np,int nb,REAL *x,REAL* Energy,REAL *Gradient,REAL_MATRIX Hessian,REAL_MATRIX3x3 *StrainFirstDerivative);

void Minimization(int k);

int BakerMinimizationNoOutput(void);

void AllocateMinimizationLocalMemory(void);
void FreeMinimizationLocalMemory(void);
void AllocateMinimizationMemory(void);

void AllocateMinimizationMemory(void);
void WriteRestartMinimization(FILE *FilePtr);
void ReadRestartMinimization(FILE *FilePtr);

void AllocateMinimizationLocalMemory(void);

#endif

