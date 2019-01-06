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

#ifndef EWALD_H
#define EWALD_H

#include "simulation.h"
#include "vector.h"
#include "utils.h"

extern int *NumberOfKVectors;
extern int OmitEwaldFourier;

extern REAL *Alpha;
extern INT_VECTOR3 *kvec;
extern REAL *ReciprocalCutOffSquared;
extern REAL EwaldPrecision;
extern REAL DielectricConstantOfTheMedium;

extern REAL *UAdsorbateAdsorbateChargeChargeFourierDelta;
extern REAL *UCationCationChargeChargeFourierDelta;
extern REAL *UAdsorbateCationChargeChargeFourierDelta;
extern REAL *UHostCationChargeChargeFourierDelta;
extern REAL *UHostAdsorbateChargeChargeFourierDelta;
extern REAL *UHostHostChargeChargeFourierDelta;

extern REAL *UAdsorbateAdsorbateChargeBondDipoleFourierDelta;
extern REAL *UCationCationChargeBondDipoleFourierDelta;
extern REAL *UAdsorbateCationChargeBondDipoleFourierDelta;
extern REAL *UHostCationChargeBondDipoleFourierDelta;
extern REAL *UHostAdsorbateChargeBondDipoleFourierDelta;
extern REAL *UHostHostChargeBondDipoleFourierDelta;

extern REAL *UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta;
extern REAL *UCationCationBondDipoleBondDipoleFourierDelta;
extern REAL *UAdsorbateCationBondDipoleBondDipoleFourierDelta;
extern REAL *UHostCationBondDipoleBondDipoleFourierDelta;
extern REAL *UHostAdsorbateBondDipoleBondDipoleFourierDelta;
extern REAL *UHostHostBondDipoleBondDipoleFourierDelta;

// Maintain these pointers so we can free them on finish
extern COMPLEX * Eikx0,* Eiky0,* Eikz0;
extern COMPLEX * Eikx_bd0,* Eiky_bd0,* Eikz_bd0;
extern COMPLEX * Eikr_bd, *Eikr_xy_bd;

// net charges
extern REAL *NetChargeSystem;
extern REAL *NetChargeFramework;
extern REAL *NetChargeCations;
extern REAL *NetChargeAdsorbates;

extern REAL NetChargeCationDelta;
extern REAL NetChargeAdsorbateDelta;

int AllocateNewEwald(void);
void SetupKVectors(void);
int PrecomputeFixedEwaldContributions(void);
int PrecomputeTotalEwaldContributions(void);
void InitializeEwald(REAL precision,int Automatic);
int EwaldEnergyIon(void);
int EwaldFourierEnergy(void);
int EwaldFourierForce(void);
int EwaldFourierBornTerm(void);
int CalculateEwaldFourierAdsorbate(int New,int old,int mol,int store);
int CalculateEwaldFourierAdsorbateCF(int New,int old,int mol,int store);
int CalculateEwaldFourierAdsorbateTupleLambda(int numberOfEwaldMolecules, int Molecules[],REAL CFChargeScalingNew[], REAL CFChargeScalingOld[],int store);
int CalculateEwaldFourierCation(int New,int old,int mol,int store);
int CalculateEwaldFourierCationCF(int New,int old,int mol,int store);
void AcceptEwaldAdsorbateMove(int A);
void AcceptEwaldCationMove(int A);
void SaveCurrentKVectors(int A,int B);
void RetrieveStoredKVectors(int A,int B);
void SaveCurrentEwaldStructureFactors(int A,int B);
void RetrieveStoredEwaldStructureFactors(int A,int B);
void SwapEwaldSystem(int A,int B);

int AllocateEwaldMemory(void);
int ReallocateEwaldChargeMemory(void);
int ReallocateEwaldBondDipoleMemory(void);

void AcceptEwaldFrameworkDiplacementMove(void);
void AcceptEwaldFrameworkDiplacementMoveRigid(void);

int EwaldFourierStaticElectricField(void);

int CalculateEwaldFourierFrameworkDisplacement(void);
int CalculateEwaldFourierFrameworkAtomTranslate(int index);
void AcceptEwaldFrameworkMove(int A);

// RXMC
int CalculateEwaldFourierAdsorbateRXMC(int reaction,REAL Lambda1,REAL Lambda2,REAL LambdaNew,int **LambdaRetraceMolecules,int direction,int store);
int CalculateEwaldFourierAdsorbateRXMC2(int reaction,REAL LambdaNew,int direction,int store);

int ComputeStaticElectricFieldEwaldAdsorbateDifference(int NewMolecule,int OldMolecule,int mol,int store);

void CalculateEwaldFourierBornTerm(REAL *Energy,REAL* Gradient,REAL_MATRIX3x3 *StrainDerivativeTensor);

void CalculateEwaldFourierCrossTerms(REAL *Energy,REAL* Gradient,REAL_MATRIX3x3 *StrainDerivativeTensor,REAL_MATRIX CrossTerms);

REAL WolfEnergyFull(void);

int ComputeStaticElectricFieldEwald(int New,int excl_ads,int excl_cation);

void WriteRestartEwald(FILE *FilePtr);
void ReadRestartEwald(FILE *FilePtr);

void ComputeElectricFieldFromInducedDipolesEwald(void);
void ComputeElectricFieldFromInducedDipolesEwaldMC(int New,int excl_ads,int excl_cation);

void ComputeInducedDipolesForcesEwald(void);

int CalculateEwaldFourierDerivatives(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainFirstDerivative,int ComputeGradient,int ComputeHessian);

REAL EwaldFourierElectrostaticPotentialAdsorbate(int m,int l);
REAL EwaldFourierElectrostaticPotentialCation(int m,int l);

#endif
