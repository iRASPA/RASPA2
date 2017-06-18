/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'ewald.h' is part of RASPA-2.0

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
