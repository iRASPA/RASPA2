/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework.h' is part of RASPA-2.0

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

#ifndef FRAMEWORK_H
#define FRAMEWORK_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <molecule.h>
#include "vector.h"

enum {FLEXIBLE_FILE_TYPE_RASPA,FLEXIBLE_FILE_TYPE_DLPOLY};

enum {IMPROPER_TORSION_SCAN_GENERAL,IMPROPER_TORSION_SCAN_UNIQUE};

extern int CorrectNetChargeOnPseudoAtom;
extern REAL CorrectNetChargeOnPseudoAtomValue;

typedef struct citation_information
{
  char CitationId[80];
  char CitationAuthorName[256];
  char CitationCoordinateLinkage[80];
  char CitationTitle[1024];
  char CitationCountry[80];
  char CitationPageFirst[80];
  char CitationPageLast[80];
  char CitationYear[80];
  char CitationJournalAbbrev[80];
  char CitationJournalVolume[80];
  char CitationJournalIssue[80];
  char CitationJournalID_ASTM[80];
  char CitationJournalID_ISSN[80];
  char CitationBookTitle[80];
  char CitationBookPublisher[80];
  char CitationBookID_ISBN[80];
  char CitationSpecialDetails[80];
} CITATION_INFORMATION;

typedef struct crystallographic_statistics
{
  VECTOR *Position;
  VECTOR *PositionSquared;
  VECTOR *Distance;
  REAL *Occupation;
  REAL *AverageDistance;
  REAL Unclassified;
  REAL *RelativeOccupation;
  int *NumberOfCationSites;
  REAL *Count;
  REAL_MATRIX3x3 *TemperatureFactor;
  REAL count;
  REAL count2;
} CRYSTALLOGRAPHIC_STATISTICS;

extern CRYSTALLOGRAPHIC_STATISTICS *crystallographic_stats;

typedef struct framework_asymmetric_atom
{
  char AtomNumberCode[32];
  int Type;                      // the pseudo-atom Type of the atom
  REAL Charge;
  int Replacable;
  int AssymetricType;
  int Identifier;
  POINT Position;
} FRAMEWORK_ASYMMETRIC_ATOM;


extern REAL CutOffIons;

extern int AsymmetricIons;
extern int SpaceGroupIons;

extern int NumberOfSingleSubstitutionRules;
extern char (*SubstitutionSingleFrameworkAtom)[3][256];
extern int (*SubstitutionSingleFrameworkAtomTypes)[3];

extern int NumberOfRandomSubstitutions;
extern int NumberOfSubstitutionRules;
extern char (*SubstitutionFrameworkAtoms)[3][256];
extern int (*SubstitutionFrameworkAtomTypes)[3];

extern int NumberOfSubstitutions;
extern int (*ListOfAtomSubstitutions)[3];

enum {MODIFY_FRAMEWORKATOM_CONNECTED_TO,MODIFY_FRAMEWORKATOM_DIMER,MODIFY_FRAMEWORKATOM_TRIPLE,MODIFY_FRAMEWORKATOM_PLANAR};
extern int NumberOfModificationRules;
extern int *ModificationRuleType;
extern char (*ModifyFrameworkAtoms)[10][256];
extern int (*ModifyFrameworkAtomTypes)[10];

extern int NumberOfForbiddenConnectivityRules;
extern char (*ForbiddenConnectivityAtoms)[3][256];
extern int (*ForbiddenConnectivityTypes)[3];

void CheckFrameworkCharges(void);
void AddAsymmetricAtom(FRAMEWORK_ASYMMETRIC_ATOM atom);
void ReadFrameworkDefinitionCSSR(void);
void ReadFrameworkDefinitionDLPOLY(void);
void CreateAsymetricFrameworkAtoms(void);
void WriteSymmetricFrameworkCssr(void);
void WriteFrameworkCssr(void);

void MakeConnectivityList(void);
void AllocateConnectivityList(void);
void FreeAllocateConnectivityList(void);
void SubstituteAtoms(void);
void ModifyAtomsConnectedToDefinedNeighbours(void);

extern VECTOR *UnitCellSize;
extern INT_VECTOR3 *NumberOfUnitCells;
extern REAL_MATRIX3x3 *UnitCellBox;
extern REAL_MATRIX3x3 *InverseUnitCellBox;
extern int *Lowenstein;

int BlockedPocket(VECTOR pos);
int BlockedPocketGridTest(VECTOR pos);

void ReadCageCenters(void);
void ReadBlockingPockets(void);

extern int RemoveBondNeighboursFromLongRangeInteraction;
extern int RemoveBendNeighboursFromLongRangeInteraction;
extern int RemoveTorsionNeighboursFromLongRangeInteraction;

extern int Remove12NeighboursFromVDWInteraction;
extern int Remove13NeighboursFromVDWInteraction;
extern int Remove14NeighboursFromVDWInteraction;

extern int Remove12NeighboursFromChargeChargeInteraction;
extern int Remove13NeighboursFromChargeChargeInteraction;
extern int Remove14NeighboursFromChargeChargeInteraction;

extern int Remove11NeighboursFromChargeBondDipoleInteraction;
extern int Remove12NeighboursFromChargeBondDipoleInteraction;
extern int Remove13NeighboursFromChargeBondDipoleInteraction;
extern int Remove14NeighboursFromChargeBondDipoleInteraction;

extern int Remove12NeighboursFromBondDipoleBondDipoleInteraction;
extern int Remove13NeighboursFromBondDipoleBondDipoleInteraction;
extern int Remove14NeighboursFromBondDipoleBondDipoleInteraction;

extern int ImproperTorsionScanType;

extern int InternalFrameworkLennardJonesInteractions;

typedef struct FrameworkComponent
{
  char (*Name)[256];                        // the name of the frameworks

  int TotalNumberOfAtoms;                   // the total number of atoms of the frameworks
  int TotalNumberOfCharges;                 // the total number of charges of the frameworks
  int TotalNumberOfUnitCellAtoms;           // the total number of atoms of the unit cell
  REAL FrameworkDensity;                    // the total density of the frameworks
  REAL FrameworkMass;                       // the total mass of the frameworks

  int NumberOfFrameworks;                   // the number of frameworks
  REAL *FrameworkDensityPerComponent;       // the density per framework
  REAL *FrameworkMassPerComponent;          // the mass per framework

  int *NumberOfAtoms;                       // the number of atoms per framework
  int *NumberOfCharges;                     // the number of atoms per framework
  int *NumberOfFixedAtoms;                  // the number of fixed atoms per framework
  int *NumberOfFreeAtoms;                   // the number of free atoms per framework
  int *NumberOfUnitCellAtoms;               // the number of unit cell atoms per framework
  int *MaxNumberOfAtoms;
  ATOM **Atoms;                             // list of framework-atoms per framework

  int **CellListHead;                        // the starting atom per cell
  int **CellList;                            // linked list of framework atoms

  REAL *FrameworkProbability;
  int FrameworkExclusion;
  int RemoveHydrogenDisorder;

  int *NumberOfCitations;
  CITATION_INFORMATION **CitationInformation;

  int FrameworkModel;
  int *FrameworkModels;
  VECTOR *ShiftUnitCell;
  int *Asymmetric;
  int *SpaceGroupIdentifier;
  int *CalculateSpaceGroup;
  int *RemoveAtomNumberCodeFromLabel;
  int *AddAtomNumberCodeToLabel;
  int *InputFileType;
  int ForceSpaceGroupDetection;

  char FrameworkDefinitions[256];
  int FlexibleModelInputType;

  int AnisotropicType;

  VECTOR IntialCenterOfMassPosition;

  REAL SurfaceArea;                        // in units of Angstrom^2
  REAL *SurfaceAreas;
  REAL SurfaceAreaCations;
  char SurfaceAreaProbeAtom[256];
  int SurfaceAreaSamplingPointsPerShere;
  REAL SurfaceAreaProbeDistance;

  REAL PoreSizeDistributionProbeDistance;

  int TranslationDirection;

  int *NumberOfAsymmetricAtoms;
  FRAMEWORK_ASYMMETRIC_ATOM **AtomsAsymmetric;

  char NameIons[256];                        // the Name of the component ("methane","C12","propane" etc).
  int NumberOfIons;
  int MaxNumberOfIons;
  ATOM *Ions;

  int NumberOfAsymmetricIons;
  FRAMEWORK_ASYMMETRIC_ATOM *IonsAsymmetric;

  int RestrictFrameworkAtomsToBox;

  int ReadCIFAsCartesian;

  char ***ExclusionMatrix;

  int **Connectivity;
  int ***Neighbours;

  int *NumberOfCoreShells;
  int **CoreShellConnectivity;

  int NumberOfCoreShellDefinitions;
  PAIR *CoreShellDefinitions;
  int *NumberOfCoreShellsPerType;

  REAL Intra14VDWScalingValue;
  REAL Intra14ChargeChargeScalingValue;

  int NumberOfIntra12Interactions;
  int NumberOfIntra13Interactions;
  int NumberOfIntra14Interactions;
  int NumberOfIntra123Interactions;
  int NumberOfIntra1234Interactions;

  int NumberOfBondsDefinitions;                                  // the number of bond definitions
  int *BondDefinitionType;                                       // the type of the bonds (i.e. HARMONIC_BOND)
  PAIR *BondDefinitions;                                         // the pair of pseudo-atoms in the bond definition
  REAL (*BondArgumentDefinitions)[MAX_BOND_POTENTIAL_ARGUMENTS]; // the parameters corresponding to this bond
  int *NumberOfBondsPerType;                                     // the number of bonds for this type
  int *NumberOfBonds;                                            // the number of Bonds (per framework)
  int *MaxNumberOfBonds;                                         // the currently allocated memory for bonds (per framework)
  PAIR **Bonds;                                                  // the list of bond-pairs (per framework)
  int **BondType;                                                // the type of the bond for each Bond-pair (per framework)
  REAL (**BondArguments)[MAX_BOND_POTENTIAL_ARGUMENTS];          // the bond arguments for this bond-pair (per framework)
  // the bondarguments are sometimes modified, e.g. taking the bond equilibrium distance from the crystal-structure
  // the distances could be different for each bond-pair

  int NumberOfBondDipoleDefinitions;     // the number of bond-dipole definitions
  PAIR *BondDipoleDefinitions;           // the pair of pseudo-atoms in the bond-dipole definition
  REAL *BondDipoleArgumentDefinition;    // the magntitude corresponding to this bond-dipole
  int *NumberOfBondDipolesPerType;       // the number of bond-dipoles for this type

  int *NumberOfBondDipoles;              // the number of bond-dipoles (per framework)
  int *MaxNumberOfBondDipoles;           // the currently allocated memory for bond-dipoles (per framework)
  PAIR **BondDipoles;                    // the list of bond-dipole-pairs (per framework)
  REAL **BondDipoleMagnitude;            // the magntiude of the bond-dipole for this bond-pair (per framework)

  int NumberOfUreyBradleyDefinitions;                                          // the number of Urey-Bradley definitions
  int *UreyBradleyDefinitionType;                                              // the type of the Urey-Bradleys
  TRIPLE *UreyBradleyDefinitions;                                              // the triple of pseudo-atoms in the Urey-Bradley definition
  REAL (*UreyBradleyArgumentDefinitions)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]; // the parameters corresponding to this Urey-Bradley
  int *NumberOfUreyBradleysPerType;                                            // the number of Urey-Bradleys for this type
  int *NumberOfUreyBradleys;                                                   // the number of Urey-Bradleys (per framework)
  int *MaxNumberOfUreyBradleys;                                                // the currently allocated memory for Urey-Bradleys (per framework)
  TRIPLE **UreyBradleys;                                                       // the list of Urey-Bradley-triples (per framework)
  int **UreyBradleyType;                                                       // the type of the Urey-Bradley for each Urey-Bradley-triple (per framework)
  REAL (**UreyBradleyArguments)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS];          // the Urey-Bradley arguments for this triple (per framework)

  int NumberOfBendDefinitions;                                   // the number of bend definitions
  int *BendDefinitionType;                                       // the type of the bend
  QUAD *BendDefinitions;                                         // the triple of pseudo-atoms in the bend definition
  REAL (*BendArgumentDefinitions)[MAX_BEND_POTENTIAL_ARGUMENTS]; // the parameters corresponding to this bend
  int *NumberOfBendsPerType;                                     // the number of bends for this type
  int *NumberOfBends;                                            // the number of bends (per framework)
  int *MaxNumberOfBends;                                         // the currently allocated memory for bends (per framework)
  QUAD **Bends;                                                  // the list of bend-triples (per framework)
  int **BendType;                                                // the type of the bend for each bend-triple (per framework)
  REAL (**BendArguments)[MAX_BEND_POTENTIAL_ARGUMENTS];          // the bend arguments for this triple (per framework)

  int NumberOfInversionBendDefinitions;                                             // the number of inversion bend definitions
  int *InversionBendDefinitionType;                                                 // the type of the inversion bend
  QUAD *InversionBendDefinitions;                                                   // the quad of pseudo-atoms in the inversion bend definition
  REAL (*InversionBendArgumentDefinitions)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]; // the parameters corresponding to this inversion bend
  int *NumberOfInversionBendsPerType;                                               // the number of inversion bends for this type
  int *NumberOfInversionBends;                                                      // the number of inversion bends (per framework)
  int *MaxNumberOfInversionBends;                                                   // the currently allocated memory for inversion bends (per framework)
  QUAD **InversionBends;                                                            // the list of inversion bend-triples (per framework)
  int **InversionBendType;                                                          // the type of the inv bend for each inversion bend-triple (per framework)
  REAL (**InversionBendArguments)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS];          // the inversion bend arguments for this triple (per framework)

  int NumberOfTorsionDefinitions;                                      // the number of torsion definitions
  int *TorsionDefinitionType;                                          // the type of the torsion
  QUAD *TorsionDefinitions;                                            // the quad of pseudo-atoms in the torsion definition
  REAL (*TorsionArgumentDefinitions)[MAX_TORSION_POTENTIAL_ARGUMENTS]; // the parameters corresponding to this torsion
  int *NumberOfTorsionsPerType;                                        // the number of torsions for this type
  int *NumberOfTorsions;                                               // the number of torsions (per framework)
  int *MaxNumberOfTorsions;                                            // the currently allocated memory for torsions (per framework)
  QUAD **Torsions;                                                     // the list of torsion-quads (per framework)
  int **TorsionType;                                                   // the type of the torsions for each torsion-quad (per framework)
  REAL (**TorsionArguments)[MAX_TORSION_POTENTIAL_ARGUMENTS];          // the torsion arguments for this quad (per framework)

  int NumberOfImproperTorsionDefinitions;                                               // the number of improper torsion definitions
  int *ImproperTorsionDefinitionType;                                                   // the type of the improper torsion
  QUAD *ImproperTorsionDefinitions;                                                     // the quad of pseudo-atoms in the improper torsion definition
  REAL (*ImproperTorsionArgumentDefinitions)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]; // the parameters corresponding to this improper torsion
  int *NumberOfImproperTorsionsPerType;                                                 // the number of improper torsions for this type
  int *NumberOfImproperTorsions;                                                        // the number of improper torsions (per framework)
  int *MaxNumberOfImproperTorsions;                                                     // the currently allocated memory for improper torsions (per framework)
  QUAD **ImproperTorsions;                                                              // the list of improper torsion-quads (per framework)
  int **ImproperTorsionType;                                                            // the type of the torsions for each improper torsion-quad (per framework)
  REAL (**ImproperTorsionArguments)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS];          // the improper torsion arguments for this quad (per framework)

  int NumberOfOutOfPlaneDefinitions;                                           // the number of out-of-plane definitions
  int *OutOfPlaneDefinitionType;                                               // the type of the out-of-plane
  QUAD *OutOfPlaneDefinitions;                                                 // the quad of pseudo-atoms in the out-of-plane definition
  REAL (*OutOfPlaneArgumentDefinitions)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]; // the parameters corresponding to this out-of-plane
  int *NumberOfOutOfPlanesPerType;                                             // the number of out-of-planes for this type
  int *NumberOfOutOfPlanes;                                                    // the number of out-of-planes (per framework)
  int *MaxNumberOfOutOfPlanes;                                                 // the currently allocated memory for out-of-plane (per framework)
  QUAD **OutOfPlanes;                                                          // the list of out-of-plane-quads (per framework)
  int **OutOfPlaneType;                                                        // the type of the out-of-plane for each out-of-plane-quad (per framework)
  REAL (**OutOfPlaneArguments)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS];          // the out-of-plane arguments for this quad (per framework)

  int NumberOfBondBondDefinitions;                               // the number of Bonds of the component
  int *BondBondDefinitionType;
  TRIPLE *BondBondDefinitions;
  REAL (*BondBondArgumentDefinitions)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS];
  int *NumberOfBondBondsPerType;
  int *NumberOfBondBonds;
  int *MaxNumberOfBondBonds;
  TRIPLE **BondBonds;
  int **BondBondType;
  REAL (**BondBondArguments)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS];

  int NumberOfBondBendDefinitions;                               // the number of Bonds of the component
  int *BondBendDefinitionType;
  TRIPLE *BondBendDefinitions;
  REAL (*BondBendArgumentDefinitions)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS];
  int *NumberOfBondBendsPerType;
  int *NumberOfBondBends;
  int *MaxNumberOfBondBends;
  TRIPLE **BondBends;
  int **BondBendType;
  REAL (**BondBendArguments)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS];

  int NumberOfBendBendDefinitions;                               // the number of Bonds of the component
  int *BendBendDefinitionType;
  QUAD *BendBendDefinitions;
  REAL (*BendBendArgumentDefinitions)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS];
  int *NumberOfBendBendsPerType;
  int *NumberOfBendBends;
  int *MaxNumberOfBendBends;
  QUAD **BendBends;
  int **BendBendType;
  REAL (**BendBendArguments)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS];

  int NumberOfBondTorsionDefinitions;
  int *BondTorsionDefinitionType;
  QUAD *BondTorsionDefinitions;
  REAL (*BondTorsionArgumentDefinitions)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS];
  int *NumberOfBondTorsionsPerType;
  int *NumberOfBondTorsions;
  int *MaxNumberOfBondTorsions;
  QUAD **BondTorsions;
  int **BondTorsionType;
  REAL (**BondTorsionArguments)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS];

  int NumberOfBendTorsionDefinitions;
  int *BendTorsionDefinitionType;
  QUAD *BendTorsionDefinitions;
  REAL (*BendTorsionArgumentDefinitions)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS];
  int *NumberOfBendTorsionsPerType;
  int *NumberOfBendTorsions;
  int *MaxNumberOfBendTorsions;
  QUAD **BendTorsions;
  int **BendTorsionType;
  REAL (**BendTorsionArguments)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS];

  int *NumberOfIntraVDW;
  int *NumberOfExcludedIntraVDW;

  int *NumberOfIntraCharges;
  int *NumberOfIntraChargeBondDipoles;
  int *NumberOfIntraBondDipoles;

  // lists for the Ewald erf-subtraction
  int *NumberOfExcludedIntraChargeCharge;
  int *MaxNumberOfExcludedIntraChargeCharge;
  PAIR **ExcludedIntraChargeCharge;

  int *NumberOfExcludedIntraChargeBondDipole;
  int *MaxNumberOfExcludedIntraChargeBondDipole;
  PAIR **ExcludedIntraChargeBondDipole;

  int *NumberOfExcludedIntraBondDipoleBondDipole;
  int *MaxNumberOfExcludedIntraBondDipoleBondDipole;
  PAIR **ExcludedIntraBondDipoleBondDipole;

} FRAMEWORK_COMPONENT;

extern FRAMEWORK_COMPONENT *Framework;
extern int CurrentFramework;

void AllocateFrameworkCellList(void);
void MakeFrameworkCellList(void);

void AllocateAnisotropicNeighbors(void);

void ReadFrameworkDefinitionMOL(void);

void ReadFrameworkSpecificDefinition(void);
int ReadFrameworkDefinition(void);
void WriteFrameworkDefinitionDLPOLY(void);

void CellProperties(REAL_MATRIX3x3 *in,REAL_MATRIX3x3 *out,REAL *Volume);

void ReadIonSitingDefinition(void);

void PrintFrameworkStatus(void);
int ReadFrameworkComponentDefinition(void);
void PrintFlexibleFrameworkModel(void);

VECTOR GetFrameworkCenterOfMass(void);
VECTOR GetIndividualFrameworkCenterOfMass(int f1);
VECTOR GetFrameworkCenterOfMassPosition(void);
VECTOR GetFrameworkCenterOfMassVelocity(void);
void InitializeFrameworkVelocities(void);
VECTOR MeasureFrameworkVelocityDrift(void);
void RemoveFrameworkVelocityDrift(void);
REAL GetFrameworkKineticEnergy(void);
void SetFrameworkVelocitesToZero(void);
void ScaleFrameworkVelocitesToTemperature(void);

void UpdateCrystallographics(void);
void GetFAUStructure(void);
void GetFAU2Structure(void);
void GetLTAStructure(void);

void PrintConnectivityList(void);
void PrintReplicaConnectivityList(void);

int BondNeighbours(int A,int B);
int BendNeighbours(int A,int C);
int TorsionNeighbours(int A,int D);

void QuenchCoreSHellVelocities(void);

int ClosestCrystallographicPosition(VECTOR pos);
void ClosestCrystallographicPosition2(VECTOR pos,int *closest,REAL *minimum_distance);

void MakeExclusionMatrix(int system);
void MakeExcludedInteractionLists(int system);

int CheckSurfaceAreaOverlap(int typeA,VECTOR posA,int skipf,int skipfa,int skipcm,int skipca);

void ExpandAsymmetricFrameworkToFullFramework(void);
void ExpandAsymmetricIonsToFullCell(void);

void IndexFrameworkAtomsToAtomIndices(void);

void ConstructBondDipolesFromBondsFramework(void);

int IsDefinedBondType(int,int,int,int);
int IsDefinedBendType(int,int,int,int,int);
int IsDefinedTorsionType(int,int,int,int,int,int);

void WriteRestartFramework(FILE *FilePtr);
void ReadRestartFramework(FILE *FilePtr);

void DetermineSpaceGroup(void);

void GenerateFramework(void);

int ReturnNumberOfAsymmetricAtoms(int sp);

void WriteAsymmetricUnitCell(int sp);

void ReadFrameworkDefinitionCIF(void);
void WriteFrameworkDefinitionCIF(char * string);

void WriteFrameworkDefinitionShell(char * string);

void WriteFrameworkDefinitionCSSR(char *string);
void WriteFrameworkDefinitionGulp(char *string);
void WriteFrameworkDefinitionVASP(char *string);
void WriteFrameworkDefinitionPDB(char *string);
void WriteFrameworkDefinitionTinker(char *string);
void WriteFrameworkDefinitionMOL(char *string);

void PutNoiseOnFrameworkAtomicPositions(void);

#endif
