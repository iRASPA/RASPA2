/*************************************************************************************************************
 The MIT License

 Copyright (C) 2006-2019 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

     D.Dubbeldam@uva.nl            http://www.uva.nl/profiel/d/u/d.dubbeldam/d.dubbeldam.html
     scaldia@upo.es                http://www.upo.es/raspa/
     t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
     don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
     snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************************/

#ifndef CBMC_H
#define CBMC_H

#include "molecule.h"

enum {LJ_BIASING,LJ_AND_REAL_BIASING};

enum {CBMC_INSERTION,CBMC_PARTIAL_INSERTION,CBMC_DELETION,CBMC_RETRACE_REINSERTION};

extern int BiasingMethod;

#define MAX_NUMBER_OF_TRIAL_POSITIONS 2000

// Added by Ambroise de Izarra
//-------------------------------------------------------------------
extern int Itrial;
//-------------------------------------------------------------------

extern int NumberOfTrialPositions;
extern int NumberOfTrialPositionsForTheFirstBead;
extern int NumberOfTrialMovesPerOpenBead;
extern int NumberOfTrialPositionsTorsion;

extern int MaxNumberOfTrialPositions;
extern int NumberOfTrialPositionsReinsertion;
extern int NumberOfTrialPositionsPartialReinsertion;
extern int NumberOfTrialPositionsIdentityChange;
extern int NumberOfTrialPositionsGibbs;
extern int NumberOfTrialPositionsSwap;
extern int NumberOfTrialPositionsWidom;

extern int MaxNumberOfTrialPositionsForTheFirstBead;
extern int NumberOfTrialPositionsForTheFirstBeadReinsertion;
extern int NumberOfTrialPositionsForTheFirstBeadPartialReinsertion;
extern int NumberOfTrialPositionsForTheFirstBeadIdentityChange;
extern int NumberOfTrialPositionsForTheFirstBeadGibbs;
extern int NumberOfTrialPositionsForTheFirstBeadSwap;
extern int NumberOfTrialPositionsForTheFirstBeadWidom;

extern REAL TargetAccRatioSmallMCScheme;
extern REAL EnergyOverlapCriteria;
extern REAL MinimumRosenbluthFactor;

extern int MaxNumberOfBeads;
extern int MaxNumberOfBonds;
extern int MaxNumberOfBondDipoles;
extern int MaxNumberOfBends;
extern int MaxNumberOfBendBends;
extern int MaxNumberOfInversionBends;
extern int MaxNumberOfUreyBradleys;
extern int MaxNumberOfTorsions;
extern int MaxNumberOfImproperTorsions;
extern int MaxNumberOfOutOfPlanes;
extern int MaxNumberOfBondBonds;
extern int MaxNumberOfBondBends;
extern int MaxNumberOfBondTorsions;
extern int MaxNumberOfBendTorsions;
extern int MaxNumberOfIntraVDW;
extern int MaxNumberOfIntraChargeCharge;
extern int MaxNumberOfIntraChargeBondDipole;
extern int MaxNumberOfIntraBondDipoleBondDipole;

extern REAL *UBondNew;
extern REAL *UUreyBradleyNew;
extern REAL *UBendNew;
extern REAL *UBendBendNew;
extern REAL *UInversionBendNew;
extern REAL *UTorsionNew;
extern REAL *UImproperTorsionNew;
extern REAL *UBondBondNew;
extern REAL *UBondBendNew;
extern REAL *UBondTorsionNew;
extern REAL *UBendTorsionNew;
extern REAL *UIntraVDWNew;
extern REAL *UIntraChargeChargeNew;
extern REAL *UIntraChargeBondDipoleNew;
extern REAL *UIntraBondDipoleBondDipoleNew;

extern REAL *UHostVDWNew;
extern REAL *UAdsorbateVDWNew;
extern REAL *UCationVDWNew;
extern REAL *UHostChargeChargeNew;
extern REAL *UAdsorbateChargeChargeNew;
extern REAL *UCationChargeChargeNew;
extern REAL *UHostChargeBondDipoleNew;
extern REAL *UAdsorbateChargeBondDipoleNew;
extern REAL *UCationChargeBondDipoleNew;
extern REAL *UHostBondDipoleBondDipoleNew;
extern REAL *UAdsorbateBondDipoleBondDipoleNew;
extern REAL *UCationBondDipoleBondDipoleNew;

extern REAL *UBondOld;
extern REAL *UUreyBradleyOld;
extern REAL *UBendOld;
extern REAL *UBendBendOld;
extern REAL *UInversionBendOld;
extern REAL *UTorsionOld;
extern REAL *UImproperTorsionOld;
extern REAL *UBondBondOld;
extern REAL *UBondBendOld;
extern REAL *UBondTorsionOld;
extern REAL *UBendTorsionOld;
extern REAL *UIntraVDWOld;
extern REAL *UIntraChargeChargeOld;
extern REAL *UIntraChargeBondDipoleOld;
extern REAL *UIntraBondDipoleBondDipoleOld;

extern REAL *UHostVDWOld;
extern REAL *UAdsorbateVDWOld;
extern REAL *UCationVDWOld;
extern REAL *UHostChargeChargeOld;
extern REAL *UAdsorbateChargeChargeOld;
extern REAL *UCationChargeChargeOld;
extern REAL *UHostChargeBondDipoleOld;
extern REAL *UAdsorbateChargeBondDipoleOld;
extern REAL *UCationChargeBondDipoleOld;
extern REAL *UHostBondDipoleBondDipoleOld;
extern REAL *UAdsorbateBondDipoleBondDipoleOld;
extern REAL *UCationBondDipoleBondDipoleOld;

// trial energies
extern REAL *UBondTrial;
extern REAL *UBendTrial;
extern REAL *UBendBendTrial;
extern REAL *UUreyBradleyTrial;
extern REAL *UInversionBendTrial;
extern REAL *UTorsionTrial;
extern REAL *UImproperTorsionTrial;
extern REAL *UBondBondTrial;
extern REAL *UBondBendTrial;
extern REAL *UBondTorsionTrial;
extern REAL *UBendTorsionTrial;
extern REAL *UIntraVDWTrial;
extern REAL *UIntraChargeChargeTrial;
extern REAL *UIntraChargeBondDipoleTrial;
extern REAL *UIntraBondDipoleBondDipoleTrial;

extern REAL *UHostVDWTrial;
extern REAL *UAdsorbateVDWTrial;
extern REAL *UCationVDWTrial;
extern REAL *UHostChargeChargeTrial;
extern REAL *UAdsorbateChargeChargeTrial;
extern REAL *UCationChargeChargeTrial;
extern REAL *UHostChargeBondDipoleTrial;
extern REAL *UAdsorbateChargeBondDipoleTrial;
extern REAL *UCationChargeBondDipoleTrial;
extern REAL *UHostChargePermanentDipoleTrial;
extern REAL *UAdsorbateChargePermanentDipoleTrial;
extern REAL *UCationChargePermanentDipoleTrial;
extern REAL *UHostBondDipoleBondDipoleTrial;
extern REAL *UAdsorbateBondDipoleBondDipoleTrial;
extern REAL *UCationBondDipoleBondDipoleTrial;
extern REAL *UHostBondDipolePermanentDipoleTrial;
extern REAL *UAdsorbateBondDipolePermanentDipoleTrial;
extern REAL *UCationBondDipolePermanentDipoleTrial;
extern REAL *UHostPermanentDipolePermanentDipoleTrial;
extern REAL *UAdsorbatePermanentDipolePermanentDipoleTrial;
extern REAL *UCationPermanentDipolePermanentDipoleTrial;

void SetConnectivityMatrix(void);
void SetGrowingStatus(void);
void Interactions(void);
int GenerateTrialOrientationsSimpleSphere(int Old);
int GenerateTrialOrientationsMCScheme(int Old);
int ComputeExternalEnergies(void);

// Added by Ambroise de Izarra
//-------------------------------------------------------------------
extern int CurrentBead;       
extern int CurrentGroup;      

extern int NumberOfPreviousBeads;  

extern int PreviousBead;           
extern int PreviousGroup;          

extern int NumberOfBeadsToBePlaced;
extern int *BeadsToBePlaced;

extern int NumberOfBeadsAlreadyPlaced; 
extern int *BeadsAlreadyPlaced; 

extern int NumberOfBranches;           
extern int *NumberOfBranchAtoms;       
extern int **BranchAtoms; 
//-------------------------------------------------------------------

extern VECTOR **NewPosition;
extern VECTOR *OldPosition;
extern VECTOR **NewVelocity;
extern VECTOR **NewForce;
extern VECTOR FirstBeadPosition;
extern VECTOR **TrialPositions;
extern REAL *CFVDWScaling;
extern REAL **CFVDWScalingRXMC;
extern REAL *CFChargeScaling;
extern REAL **CFChargeScalingRXMC;
extern int OVERLAP;

// Added by Ambroise de Izarra
//-------------------------------------------------------------------
int Rosen(void);
//-------------------------------------------------------------------

extern int **MoleculeTodoConnectivity;

void CalculateAnisotropicTrialPositions(int TypeMolA,VECTOR *TrialPosition,VECTOR *TrialAnisotropicPosition);

REAL GrowMolecule(int Iicode);
REAL RetraceMolecule(int Iicode);

void GrowReservoirMolecule(void);
void GrowReservoirMolecules(int reaction);
void GrowReservoirMolecules2(int reaction);

void MakeInitialAdsorbate(void);
void MakeInitialAdsorbates(int n,int type);
void MakeInitialCation(void);
void MakeInitialCations(int n,int type);

void RescaleMaximumRotationAnglesSmallMC(void);
void InitializeSmallMCStatisticsAllSystems(void);
void PrintSmallMCAddStatistics(FILE *FilePtr);

void CheckConfigMoves(void);

void WriteRestartCBMC(FILE *FilePtr);
void AllocateCBMCMemory(void);
void ReadRestartCBMC(FILE *FilePtr);

#endif

