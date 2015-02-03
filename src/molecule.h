/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2015 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'molecule.h' is part of RASPA-2.0

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

#ifndef MOLECULE_H
#define MOLECULE_H

#include <stdio.h>
#include "utils.h"
#include "simulation.h"
#include "cubic_spline_1d.h"

#define MAX_BOND_POTENTIAL_ARGUMENTS 20
#define MAX_BEND_POTENTIAL_ARGUMENTS 10
#define MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS 10
#define MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS 10
#define MAX_TORSION_POTENTIAL_ARGUMENTS 10
#define MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS 10
#define MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS 10
#define MAX_BOND_BOND_POTENTIAL_ARGUMENTS 10
#define MAX_BOND_BEND_POTENTIAL_ARGUMENTS 10
#define MAX_BEND_BEND_POTENTIAL_ARGUMENTS 10
#define MAX_BOND_TORSION_POTENTIAL_ARGUMENTS 10
#define MAX_BEND_TORSION_POTENTIAL_ARGUMENTS 10
#define MAX_INTRA_VDW_POTENTIAL_ARGUMENTS 10

#define MAX_NUMBER_OF_PRISMS 64
#define MAX_NUMBER_OF_CYLINDERS 64
#define MAX_NUMBER_OF_SPHERES 64

enum {HYBRIDIZATION_UNINITIALIZED,SP3,SP2,SP};

#define MAX_SPLINE_DATA 5000

#define MAX_CONNECTIVITY 8

extern int CurrentComponent;                // the current component

extern int MoleculeId;
extern int AtomId;

extern int WidomParticleInsertionComponent;

enum {CORE,SHELL};

// Bond description A-B
typedef struct pair
{
  int A;
  int B;
} PAIR;

// Bend description A-B-C
typedef struct triple
{
  int A;
  int B;
  int C;
} TRIPLE;

// Torsion/Inversion bend description A-B-C-D
typedef struct quad
{
  int A;
  int B;
  int C;
  int D;
} QUAD;

typedef struct quint
{
  int A;
  int B;
  int C;
  int D;
  int E;
} QUINT;

typedef struct sext
{
  int A;
  int B;
  int C;
  int D;
  int E;
  int F;
} SEXT;


// Pseudoatoms
typedef struct PseudoAtom
{
  char Name[32];                               // the atom type
  char PrintToPDBName[32];                     // the string to print to a pdb-file as name
  int  PrintToPDB;                             // whether to write this atom to the pdf-file or not
  char ChemicalElement[32];                    // the chemical element ("O", "H", etc)
  char OxidationStateString[32];               // the oxydation state ("3+", "2-", etc)
  int OxidationState;
  int Hybridization;                           // Hybridization type (SP3,SP2,SP)
  int ScatteringType;                          // the scattering type (powder diffraction)
  int AnomalousScatteringType;                 // the anmalous scattering type (powder diffraction)
  REAL TemperatureFactor;                      // the temperature factor (powder diffraction)
  REAL_MATRIX3x3 AnisotropicTemperatureFactor; //
  REAL Mass;                                   // the mass of the pseudo-atom
  REAL Charge1;                                // the charge of the pseudo-atom
  REAL_MATRIX3x3 Polarization;                 // the polarization of the atom
  int HasCharges;                              // whether or not the atom-type has charges
  int IsPolarizable;                           // whether or not the atom-type has a induced point dipole
  int Interaction;                             // whether or not the atom-type has interactions
  REAL Radius;                                 // the radius (used for calculating Bonds in the zeolite)
  int Connectivity;                            // the connectivity (used for calculating Bonds/Bends/Torsion in the framework)
  int TinkerType;                              // the tinker atom type for use in txyz-output
  int AnisotropicCorrection;                   // whether the positions of the sites are modified
  REAL AnisotropicDisplacement;                // the anisotropic displacement (in length-units)
  int AnisotropicType;                         // the anisotropic type (relative or absolute)
  int FrameworkAtom;                           // whether or nor the atom-type is a framework atom
  int OxidationNumber;                         // the oxidation number
  REAL Occupancy;
  REAL ScatteringDispersionReal;
  REAL ScatteringDispersionImaginary;
  char ScatteringSource[80];
  int HasVDWInteraction;
  int CoreShell;
  int Swapable;
  int CF;
  int ChargeDefinitionType;
} PSEUDO_ATOM;

extern int NumberOfPseudoAtoms;
extern PSEUDO_ATOM *PseudoAtoms;
extern int **NumberOfPseudoAtomsCount;
extern int **NumberOfPseudoAtomsType;
extern int *NumberOfPseudoAtomsTypeNew;
extern int *NumberOfPseudoAtomsTypeOld;

extern int *MapPseudoAtom;

typedef struct group_definitions
{
  int Rigid;                        // whether or not the group is rigid
  int Type;                         // the type, NONLINEAR_MOLECULE, LINEAR_MOLECULE, POINT_PARTICLE

  REAL Mass;                        // the mass of the group

  int NumberOfGroupAtoms;           // the numer of atoms in the group
  int *Atoms;                       // the atoms in the group

  REAL_MATRIX3x3 InertiaTensor;     // the inertia tensor
  VECTOR InertiaVector;             // the inertia vector
  VECTOR InverseInertiaVector;      // the inverse of inertia vector

  REAL_MATRIX3x3 InverseOriginalCoordinateSystem;  // the rotational matrix
  TRIPLE orientation;               // three atoms to compute quaternions from
  REAL rot_min;

  int HasPermanentDipoles;
  int NumberOfPermanentDipoles;
  VECTOR *PermanentDipoles;
  VECTOR *PermanentDipolePositions;

  int HasPolarizabilities;
  int NumberOfPolarizabilities;
  REAL_MATRIX3x3 *Polarizabilites;
  VECTOR *PolarizabilityPositions;

  int RotationalDegreesOfFreedom;   // the rotational degrees of freedom
} GROUP_DEFINITION;

typedef struct group
{
  REAL Mass;                               // mass of the rigid unit
  QUATERNION Quaternion;                   // orientation of the unit
  QUATERNION QuaternionMomentum;           // quaternion momentum
  QUATERNION QuaternionReferenceMomentum;  // quaternion momentum
  QUATERNION QuaternionForce;              // quaternion force
  VECTOR Torque;                           // torque vector
  VECTOR CenterOfMassPosition;             // the center of mass position
  VECTOR CenterOfMassReferencePosition;    // the reference position for the center of mass
  VECTOR CenterOfMassVelocity;             // the center of mass velocity
  VECTOR CenterOfMassReferenceVelocity;    // the center of mass velocity
  VECTOR CenterOfMassForce;                // the center of mass force
  VECTOR AngularVelocity;                  // the angular velocity of the rigid unit
  VECTOR AngularReferenceVelocity;         // the angular velocity of the rigid unit

  VECTOR EulerAxis;
  VECTOR ReferenceEulerAxis;


  INT_VECTOR3 HessianIndex;
  INT_VECTOR3 HessianIndexOrientation;

  INT_VECTOR3 FixedCenterOfMass;
  INT_VECTOR3 FixedOrientation;
} GROUP;

// Adsorbate
typedef struct atom
{
  int Type;                       // the pseudo-atom type of the atom
  REAL Charge;
  REAL CFVDWScalingParameter;    // the Van der Waals scaling parameter for the Continuous Fraction method
  REAL CFChargeScalingParameter; // the electrostatics scaling parameter for the Continuous Fraction method
  int Modified;                   // not used (yet)
  int OriginalType;               // not used (yet)
  int CreationState;
  int AssymetricType;             // the 'asymmetric' type

  // MC/MD properties
  POINT Position;                 // the position of the atom
  POINT AnisotropicPosition;      // the modified position of the atom (anisotropic model)
  POINT ReferencePosition;        // the 'reference' position of the atom
  POINT ReferenceAnisotropicPosition;        // the 'reference' position of the atom
  POINT RattleReferencePosition;  // the second 'reference' position of the atom

  // MD properties
  VECTOR Velocity;                // the velocity of the atom
  VECTOR ReferenceVelocity;       // the 'reference' velocity of the atom
  VECTOR Force;                   // the force acting on the atom
  VECTOR ReferenceForce;          // the force acting on the atom
  VECTOR RattleGradient;          //

  VECTOR ElectricField;           // the electricfield vector
  VECTOR ReferenceElectricField;  // the 'reference' electricfield vector
  VECTOR InducedElectricField;    // the induced electric field
  VECTOR InducedDipole;           // the induced dipole moment on this atom

  INT_VECTOR3 HessianIndex;       // the index in the Hessian matrix for this atom
  int HessianAtomIndex;           // the index in a hypothetical atomic Hessian matrix

  INT_VECTOR3 Fixed;
} ATOM;

typedef struct adsorbate
{
  int Type;               // the component type of the molecule
  int NumberOfAtoms;      // the number of atoms in the molecule
  REAL Pi;                // the Pi for the fractional molecule, Pi=lambda/no. of atoms
  GROUP *Groups;          // data of the rigid groups
  ATOM *Atoms;            // list of atoms
} ADSORBATE_MOLECULE;

extern int *MaxNumberOfAdsorbateMolecules;
extern int *MaxNumberOfCationMolecules;
extern int *NumberOfFractionalMolecules;
extern int *NumberOfFractionalAdsorbateMolecules;
extern int *NumberOfFractionalCationMolecules;
extern int *NumberOfReactionMolecules;
extern int *NumberOfReactionAdsorbateMolecules;
extern int *NumberOfReactionCationMolecules;

extern int MaxNumberOfCoulombicSites;
extern int MaxNumberOfBondDipoleSites;
extern int LargestNumberOfCoulombicSites;
extern int LargestNumberOfBondDipoleSites;
extern int *NumberOfAtomsPerSystem;
extern int *NumberOfChargesPerSystem;
extern int *NumberOfBondDipolesPerSystem;

extern int *NumberOfAdsorbateMolecules;
extern int CurrentAdsorbateMolecule;
extern ADSORBATE_MOLECULE **Adsorbates;

typedef struct cation
{
  int Type;               // the compount-Type of the molecule
  int NumberOfAtoms;      // the number of atoms in the molecule
  REAL Pi;                // the Pi for the fractional molecule, Pi=lambda/no. of atoms
  GROUP *Groups;          // rigid groups
  ATOM  *Atoms;
} CATION_MOLECULE;

extern int *NumberOfCationMolecules;
extern int CurrentCationMolecule;
extern CATION_MOLECULE **Cations;

typedef struct Component
{
  char Name[256];                     // the name of the component ("methane","C12","propane" etc).
  int NumberOfAtoms;                  // the number of atoms in the component
  int NumberOfCharges;                // the number of atoms in the component
  int NumberOfRigidAtoms;             // the number of atoms in rigid units
  int NumberOfFlexibleAtoms;          // the number of atoms in flexible units
  REAL Mass;                          // the mass of the component
  int *NumberOfMolecules;             // the number of molecules of the component for each system
  int *Fixed;                         // array for the atoms denoting whether to consider the atom fixed in space or free
  int *Type;                          // the pseudo-atom Type of each atom
  REAL* Charge;
  int *Connectivity;                  // the connectivity of each atom
  int **ConnectivityList;             // the neighbors of each atom
  int **ConnectivityMatrix;           // the connectivity matrix
  int HasCharges;                     // whether the molecule contains charges or not
  int IsPolarizable;                  // whether the molecule has point dipoles or not
  int ExtraFrameworkMolecule;         // TRUE: Cation, FALSE: Adsorbate
  int Swapable;                       // whether or not the number of molecules is fluctuating (i.e. GCMC)
  int Widom;                          // whether this component is used for Widom insertions

  int ***ReactantFractionalMolecules; // the index of the fractional molecules for the particular reaction
  int ***ProductFractionalMolecules;  // the index of the fractional molecules for the particular reaction
  REAL PartitionFunction;             // the partition function of the component
  int *RXMCMoleculesPresent;          // whether or not there is a fractional molecule present for each system
  int *NumberOfRXMCMoleculesPresent;  // the number of fractional molecules present of the component for each system

  TRIPLE orientation;              // three atoms to define an orientation for sampling orientations

  int *FractionalMolecule;         // the fractional molecule for each system
  int *CFMoleculePresent;          // whether or not there is a fractional molecule present for each system
  REAL *CFWangLandauScalingFactor; // the CF Wang-Landau scaling factor
  int CFLambdaHistogramSize;       // the size of the histogram and biasing-array
  REAL **CFBiasingFactors;         // the biasing-factors for conventional CF

  REAL *IdealGasRosenbluthWeight;  // the Rosenbluth weight of an ideal-chain per system
  REAL *IdealGasTotalEnergy;       // the total energy of an ideal-chain per system

  REAL *PartialPressure;           // the partial pressure of the component per system
  REAL *FugacityCoefficient;       // the fugacity coefficient of the component per system
  REAL *BulkFluidDensity;          // the bulkfluid-density of the component per system
  REAL *Compressibility;           // the compresibility of the fluid-fase per system
  REAL *MolFraction;               // the mol-fraction of the component per system
  REAL *AmountOfExcessMolecules;   // the amount of excess molecules per syste,

  REAL CriticalTemperature;        // the critical temperature of the component
  REAL CriticalPressure;           // the critical pressure of the component
  REAL AcentricFactor;             // the acentric factor of the component

  int NumberOfGroups;              // the number of groups
  GROUP_DEFINITION *Groups;        // the definition of the groups
  int *group;                      // to which group an atom belongs
  VECTOR *Positions;               // the positions in the body-fixed frame

  int NumberOfHessianIndices;      // the amount of indices of the Hessian matrix

  int AnisotropicType;

  int DegreesOfFreedom;
  int TranslationalDegreesOfFreedom;
  int RotationalDegreesOfFreedom;
  int VibrationalDegreesOfFreedom;
  int ConstraintDegreesOfFreedom;

  char MoleculeDefinition[256];
  REAL *MOLEC_PER_UC_TO_MOL_PER_KG;
  REAL *MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK;
  REAL *MOLEC_PER_UC_TO_CC_STP_G;
  REAL *MOLEC_PER_UC_TO_CC_STP_CC;
  REAL *MOL_PER_KG_TO_CC_STP_G;
  REAL *MOL_PER_KG_TO_CC_STP_CC;

  int ReadBiasingFunction;
  char BiasingFunctionName[256];
  CUBIC_SPLINE BiasingFunction;
  int Biased;
  int BiasingDirection;
  REAL *FreeEnergyXdata;
  REAL QStarA;
  REAL gA,gB;
  REAL LeftBoundary,RightBoundary;
  REAL Normalization;
  REAL TST_estimate;
  REAL RuizMonteroFactor;
  REAL UmbrellaFactor;

  int *BlockPockets;
  int *ComputeFreeEnergyProfile;
  char (*BlockPocketsFilename)[256];
  int *NumberOfBlockCenters;
  REAL **BlockDistance;
  VECTOR **BlockCenters;

  int NumberOfConstraintBonds;
  int NumberOfConstraintBends;
  int NumberOfConstraintInversionBends;
  int NumberOfConstraintTorsions;
  int NumberOfConstraintImproperTorsions;
  int NumberOfConstraintOutOfPlanes;

  int NumberOfBonds;                                    // the number of bonds of the component
  PAIR *Bonds;                                          // the list of bond-pairs
  int *BondType;                                        // the type of the bond for each bond-pair
  REAL (*BondArguments)[MAX_BOND_POTENTIAL_ARGUMENTS];  // the arguments needed for this bond-pair

  int NumberOfBondDipoles;                                       // the number of bond dipoles
  PAIR *BondDipoles;                  // an A-B pair of atoms forming the bond-dipole
  REAL *BondDipoleMagnitude;              // the fixed value of the bond-dipole

  int NumberOfBends;
  QUAD *Bends;
  int *BendType;
  REAL (*BendArguments)[MAX_BEND_POTENTIAL_ARGUMENTS];

  int NumberOfUreyBradleys;
  TRIPLE *UreyBradleys;
  int *UreyBradleyType;
  REAL (*UreyBradleyArguments)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS];

  int NumberOfInversionBends;
  QUAD *InversionBends;
  int *InversionBendType;
  REAL (*InversionBendArguments)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS];

  int NumberOfTorsions;                     // the number of Torsions of the component
  QUAD *Torsions;
  int *TorsionType;
  REAL (*TorsionArguments)[MAX_TORSION_POTENTIAL_ARGUMENTS];

  int NumberOfImproperTorsions;                     // the number of Torsions of the component
  QUAD *ImproperTorsions;
  int *ImproperTorsionType;
  REAL (*ImproperTorsionArguments)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS];

  int NumberOfOutOfPlanes;                     // the number of Torsions of the component
  QUAD *OutOfPlanes;
  int *OutOfPlaneType;
  REAL (*OutOfPlaneArguments)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS];

  int NumberOfBondBonds;                                // the number of Bonds of the component
  TRIPLE *BondBonds;                  // the list of Bond-pairs
  int *BondBondType;                // the Type of the Bond for each Bond-pair
  REAL (*BondBondArguments)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS];

  int NumberOfBondBends;                                // the number of Bonds of the component
  TRIPLE *BondBends;                  // the list of Bond-pairs
  int *BondBendType;                // the Type of the Bond for each Bond-pair
  REAL (*BondBendArguments)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS];

  int NumberOfBendBends;                                // the number of Bonds of the component
  QUAD *BendBends;                  // the list of Bond-pairs
  int *BendBendType;                // the Type of the Bond for each Bond-pair
  REAL (*BendBendArguments)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS];

  int NumberOfBondTorsions;                                // the number of Bonds of the component
  QUAD *BondTorsions;                  // the list of Bond-pairs
  int *BondTorsionType;                // the Type of the Bond for each Bond-pair
  REAL (*BondTorsionArguments)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS];

  int NumberOfBendTorsions;                                // the number of Bonds of the component
  QUAD *BendTorsions;                  // the list of Bond-pairs
  int *BendTorsionType;                // the Type of the Bond for each Bond-pair
  REAL (*BendTorsionArguments)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS];

  int NumberOfIntraVDW;         // the number of intra molecular Lennard-Jones interactions
  PAIR *IntraVDW;               // the list of Lennard-Jones-pairs
  REAL Intra14VDWScalingValue;  // the scaling factor for intra-VDW
  REAL *IntraVDWScaling;        // the scaling factors of the intra-VDW per pairs

  int NumberOfIntraChargeCharge;     // the number of intra molecular Coulombic interactions
  PAIR *IntraChargeCharge;           // the list of Coulombic-pairs
  REAL Intra14ChargeChargeScalingValue;   // the scaling factor for intra charge/charge
  REAL *IntraChargeChargeScaling;    // the scaling factors of the intra-charge/charge pairs

  int NumberOfIntraChargeBondDipole; // the number of intra molecular Coulombic interactions
  PAIR *IntraChargeBondDipole;       // the list of Coulombic-pairs

  int NumberOfIntraBondDipoleBondDipole;   // the number of intra molecular Coulombic interactions
  PAIR *IntraBondDipoleBondDipole;         // the list of Coulombic-pairs

  int NumberOfExcludedIntraChargeCharge;
  PAIR *ExcludedIntraChargeCharge;

  int NumberOfExcludedIntraChargeBondDipole;
  PAIR *ExcludedIntraChargeBondDipole;

  int NumberOfExcludedIntraBondDipoleBondDipole;
  PAIR *ExcludedIntraBondDipoleBondDipole;

  // MC part
  int StartingBead;        // the bead of the molecule used for starting the growing process

  int *CreateNumberOfMolecules;

  int NumberOfGibbsIdentityChanges;
  int *GibbsIdentityChanges;

  int NumberOfConfigMoves;             // the number of config move for the components
  int *NumberOfUnchangedAtomsConfig;   // for every config move the number of unchanged atoms
  int **UnchangedAtomsConfig;   // for every config move a list of the unchanged atoms

  int NumberOfIdentityChanges;
  int *IdentityChanges;
  int NumberOfIdentityConfigMoves;             // the number of config move for the components
  int *NumberOfUnchangedAtomsIdentityConfig;   // for every config move the number of unchanged atoms
  int **UnchangedAtomsIdentityConfig;          // for every config move a list of the unchanged atoms


  int TransmissionCoefficient;

  int TranslationMethod;
  REAL ProbabilityTranslationMove;
  int  TranslationDirection;
  int  SwapEvery;
  REAL_MATRIX3x3 TranslationMatrix;
  REAL ProbabilityRandomTranslationMove;
  REAL ProbabilityRotationMove;
  REAL ProbabilityRandomRotationMove;
  REAL ProbabilityPartialReinsertionMove;
  REAL ProbabilityReinsertionMove;
  REAL ProbabilityReinsertionInPlaceMove;
  REAL ProbabilityReinsertionInPlaneMove;
  REAL ProbabilityIdentityChangeMove;
  REAL ProbabilitySwapMove;
  REAL ProbabilityCFSwapLambdaMove;
  REAL ProbabilityCBCFSwapLambdaMove;
  REAL ProbabilityWidomMove;
  REAL ProbabilitySurfaceAreaMove;
  REAL ProbabilityGibbsChangeMove;
  REAL ProbabilityGibbsIdentityChangeMove;
  REAL ProbabilityCFGibbsChangeMove;
  REAL ProbabilityCBCFGibbsChangeMove;

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
  REAL ProbabilityCFCRXMCLambdaChangeMove;

  REAL *CpuTimeTranslationMove;
  REAL *CpuTimeRandomTranslationMove;
  REAL *CpuTimeRotationMove;
  REAL *CpuTimeRandomRotationMove;
  REAL *CpuTimePartialReinsertionMove;
  REAL *CpuTimeReinsertionMove;
  REAL *CpuTimeReinsertionInPlaceMove;
  REAL *CpuTimeReinsertionInPlaneMove;
  REAL *CpuTimeIdentityChangeMove;
  REAL *CpuTimeSwapMoveInsertion;
  REAL *CpuTimeSwapMoveDeletion;
  REAL *CpuTimeCFSwapLambdaMove;
  REAL *CpuTimeCBCFSwapLambdaMove;
  REAL *CpuTimeWidomMove;
  REAL *CpuTimeSurfaceAreaMove;
  REAL *CpuTimeGibbsChangeMove;
  REAL *CpuTimeGibbsIdentityChangeMove;
  REAL *CpuTimeCFGibbsChangeMove;
  REAL *CpuTimeCBCFGibbsChangeMove;

  int RestrictMovesToBox;
  VECTOR BoxAxisABC_Min,BoxAxisABC_Min2,BoxAxisABC_Min3,BoxAxisABC_Min4;
  VECTOR BoxAxisABC_Max,BoxAxisABC_Max2,BoxAxisABC_Max3,BoxAxisABC_Max4;

  int RestrictMoves;

  int RestrictMovesToPrisms;
  int RestrictMovesToPrism[MAX_NUMBER_OF_PRISMS];
  VECTOR RestrictPrismABC_Min[MAX_NUMBER_OF_PRISMS];
  VECTOR RestrictPrismABC_Max[MAX_NUMBER_OF_PRISMS];

  int RestrictMovesToCylinders;
  int RestrictMovesToCylinder[MAX_NUMBER_OF_CYLINDERS];
  VECTOR RestrictCylinderABC_Min[MAX_NUMBER_OF_CYLINDERS];
  VECTOR RestrictCylinderABC_Max[MAX_NUMBER_OF_CYLINDERS];
  VECTOR RestrictCylinderCenter[MAX_NUMBER_OF_CYLINDERS];
  int RestrictCylinderDirection[MAX_NUMBER_OF_CYLINDERS];
  REAL RestrictCylinderRadius[MAX_NUMBER_OF_CYLINDERS];

  int RestrictMovesToSpheres;
  int RestrictMovesToSphere[MAX_NUMBER_OF_SPHERES];
  VECTOR RestrictSphereCenter[MAX_NUMBER_OF_SPHERES];
  REAL RestrictSphereRadius[MAX_NUMBER_OF_SPHERES];


  REAL FractionOfTranslationMove;
  REAL FractionOfRandomTranslationMove;
  REAL FractionOfRotationMove;
  REAL FractionOfRandomRotationMove;
  REAL FractionOfPartialReinsertionMove;
  REAL FractionOfReinsertionMove;
  REAL FractionOfReinsertionInPlaceMove;
  REAL FractionOfReinsertionInPlaneMove;
  REAL FractionOfIdentityChangeMove;
  REAL FractionOfSwapMove;
  REAL FractionOfCFSwapLambdaMove;
  REAL FractionOfCBCFSwapLambdaMove;
  REAL FractionOfWidomMove;
  REAL FractionOfSurfaceAreaMove;
  REAL FractionOfGibbsChangeMove;
  REAL FractionOfGibbsIdentityChangeMove;
  REAL FractionOfCFGibbsChangeMove;
  REAL FractionOfCBCFGibbsChangeMove;

  REAL FractionOfParallelTemperingMove;
  REAL FractionOfHyperParallelTemperingMove;
  REAL FractionOfParallelMolFractionMove;
  REAL FractionOfChiralInversionMove;
  REAL FractionOfHybridNVEMove;
  REAL FractionOfHybridNPHMove;
  REAL FractionOfHybridNPHPRMove;
  REAL FractionOfVolumeChangeMove;
  REAL FractionOfBoxShapeChangeMove;
  REAL FractionOfGibbsVolumeChangeMove;
  REAL FractionOfFrameworkChangeMove;
  REAL FractionOfFrameworkShiftMove;
  REAL FractionOfCFCRXMCLambdaChangeMove;

  int NumberOfChiralityCenters;

  int *Chirality;             // is this a stereocenter ?
  int *ChiralityType;         // is this the R stereocenter ?

  int *ChiralA;             // atom with highest priority
  int *ChiralB;             // atom with highest priority
  int *ChiralC;             // atom with highest priority
  int *ChiralD;             // atom with highest priority

  int RestrictEnantionface;
  int Enantioface;
  int EnantiofaceAtomDefinitions[5][3];
  ATOM *EnantiofaceAtoms[5];

  int LMCMOL;
  VECTOR *RMCMOL;

  REAL **MaximumCBMCChangeBondLength;       // maximum rotation angle for cone angle
  REAL **MaximumCBMCChangeBendAngle;       // maximum rotation angle for cone angle
  REAL **MaximumCBMCRotationOnCone;       // maximum rotation angle for position on a cone
  REAL **CBMCChangeBondLengthAttempts;     // performance for choise=0
  REAL **CBMCChangeBondLengthAccepted;     // performance for choise=0
  REAL **CBMCChangeBendAngleAttempts;     // performance for choise=1
  REAL **CBMCChangeBendAngleAccepted;     // performance for choise=1
  REAL **CBMCRotationOnConeAttempts;      // performance for choise=2
  REAL **CBMCRotationOnConeAccepted;      // performance for choise=2
  REAL **TotalCBMCChangeBondLengthAttempts;
  REAL **TotalCBMCChangeBendAngleAttempts;
  REAL **TotalCBMCRotationOnConeAttempts;
  REAL **TotalCBMCChangeBondLengthAccepted;
  REAL **TotalCBMCChangeBendAngleAccepted;
  REAL **TotalCBMCRotationOnConeAccepted;
} COMPONENT;

extern int NumberOfComponents;
extern int NumberOfAdsorbateComponents;
extern int NumberOfCationComponents;
extern COMPONENT *Components;

int SelectRandomMoleculeOfTypeExcludingReactionMolecules(int reaction,int **LambdaRetraceMolecules);
int SelectRandomMoleculeOfTypeExcludingProductMolecules(int reaction,int **LambdaRetraceMolecules);

int SelectRandomMoleculeOfType(int comp);
int SelectRandomMoleculeOfTypeExcludingFractionalMolecule(int comp);

int ReturnPseudoAtomNumber(char *buffer);
int ReturnPossiblePseudoAtomNumber(char *buffer);
int CheckPseudoAtomNumber(char *buffer);

int AddPseudoAtom(PSEUDO_ATOM atom);

void ReadPseudoAtomsDefinitions(void);

void AllocateComponentMemory(void);
void ReadComponentDefinition(int comp);

void InsertAdsorbateMolecule(void);
void RemoveAdsorbateMolecule(void);

void InsertCationMolecule(void);
void RemoveCationMolecule(void);

void PrintAdsorbateMolecules(void);
void PrintCationMolecules(void);

void RescaleComponentProbabilities(void);

void UpdateGroupCenterOfMassAdsorbate(int m);
void UpdateGroupCenterOfMassCation(int m);

VECTOR GetCenterOfMassCurrentSystem(void);
VECTOR GetCenterOfCurrentSystem(void);

VECTOR GetCenterOfMassVelocityCurrentSystem(void);
VECTOR GetAdsorbateCenterOfMass(int m);
VECTOR GetAdsorbateCenterOfCharge(int m);
VECTOR GetAdsorbateCenterOfMassVelocity(int m);
VECTOR GetCationCenterOfMassVelocity(int m);
VECTOR GetCationCenterOfMass(int m);
VECTOR GetCationCenterOfCharge(int m);
VECTOR GetAdsorbateCenterOfMassForce(int m);

VECTOR ComputeDipoleMomentAdsorbate(int m);
VECTOR ComputeDipoleMomentCation(int m);

REAL GetAdsorbateMass(int m);
REAL GetCationMass(int m);
REAL GetTotalAdsorbateMass(void);
REAL GetTotalCationMass(void);

VECTOR ComputeTotalDipoleMomentSystemAdsorbates(void);
VECTOR ComputeTotalDipoleMomentSystemCations(void);

void AdjustVelocitiesToTemperature(void);
void InitializeVelocityAdsorbate(int m);
void InitializeAdsorbateVelocities(void);
void InitializeVelocityAdsorbateToZero(int m);
void InitializeAdsorbateVelocitiesToZero(void);
void InitializeVelocityCation(int m);
void InitializeCationVelocities(void);

REAL GetAdsorbateKineticEnergy(void);
REAL GetCationKineticEnergy(void);

VECTOR MeasureVelocityDrift(void);
void RemoveVelocityDrift(void);
void ScaleVelocitesToTemperature(void);
void ScaleAdsorbateVelocitesToTemperature(void);
void ScaleCationVelocitesToTemperature(void);

void SetAdsorbateVelocitesToZero(void);
void SetCationVelocitesToZero(void);

VECTOR MapToBox(VECTOR pos);
void ComputeInertiaTensorGroups(int comp);

VECTOR ComputeDipoleMomentComponent(int,int);
REAL_MATRIX3x3 ComputeQuadrupoleMomentComponent(int,int);

void MoveAdsorbateCenterOfMassBackInBox(int m);
void MoveCationCenterOfMassBackInBox(int m);

void PrintNumericalVersusAnalyticalForces(int m);
void PrintNumericalVersusAnalyticalForcesFramework(void);

void ConstructBondDipolesFromBondsAdsorbates(void);

int ReturnAtomBondedToHydrogen(int Type,int c);

void CalculateAnisotropicSites(void);

void ReadBiasingProfile(int);

int IsFractionalAdsorbateMolecule(int m);
int IsFractionalCationMolecule(int m);

int IsFractionalCFMCAdsorbateMolecule(int m);
int IsFractionalCFMCCationMolecule(int m);

int IsFractionalReactionAdsorbateMolecule(int m);
int IsFractionalReactionCationMolecule(int m);

int ValidFractionalPoint(int i, POINT s);
int ValidCartesianPoint(int i, POINT pos);

void CheckChiralityMolecules(void);

void PrintCPUStatistics(FILE *FilePtr);

void WriteRestartPseudoAtoms(FILE *FilePtr);
void ReadRestartPseudoAtoms(FILE *FilePtr);

void WriteRestartMolecules(FILE *FilePtr);
void ReadRestartMolecules(FILE *FilePtr);

void WriteRestartComponent(FILE *FilePtr);
void ReadRestartComponent(FILE *FilePtr);

#endif
