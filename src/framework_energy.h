/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_energy.h' is part of RASPA-2.0

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

#ifndef FRAMEWORK_ENERGY_H
#define FRAMEWORK_ENERGY_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <molecule.h>

REAL CalculateFrameworkVDWEnergyAtPosition(POINT posA,int typeA,REAL scaling);
REAL CalculateFrameworkChargeChargeEnergyAtPosition(POINT pos,int typeA,REAL scaling);

void CalculateFrameworkChargeEnergyAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,REAL scaling);
void CalculateFrameworkBondDipoleEnergyAtPosition(VECTOR posA1,VECTOR posA2,REAL DipoleMagnitudeA,REAL *UChargeBondDipole,REAL *UBondDipoleBondDipole);

REAL CalculateFrameworkBondEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkUreyBradleyEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBendEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkInversionBendEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkTorsionEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkImproperTorsionEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBondBondEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBondBendEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBendBendEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBondTorsionEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBendTorsionEnergy(int flag,int f1,int atom_id);

int CalculateFrameworkIntraVDWEnergy(void);
int CalculateFrameworkIntraReplicaVDWEnergy(void);
int CalculateFrameworkIntraChargeChargeEnergy(void);
int CalculateFrameworkIntraReplicaChargeChargeEnergy(void);
int CalculateFrameworkIntraChargeBondDipoleEnergy(void);
int CalculateFrameworkIntraBondDipoleBondDipoleEnergy(void);

int CalculateFrameworkAdsorbateVDWEnergy(void);
int CalculateFrameworkAdsorbateReplicaVDWEnergy(void);
int CalculateFrameworkCationVDWEnergy(void);
int CalculateFrameworkCationReplicaVDWEnergy(void);

int CalculateFrameworkAdsorbateChargeChargeEnergy(void);
int CalculateFrameworkAdsorbateReplicaChargeChargeEnergy(void);
int CalculateFrameworkCationChargeChargeEnergy(void);
int CalculateFrameworkCationReplicaChargeChargeEnergy(void);

int CalculateFrameworkAdsorbateChargeBondDipoleEnergy(void);
int CalculateFrameworkCationChargeBondDipoleEnergy(void);

int CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergy(void);
int CalculateFrameworkCationBondDipoleBondDipoleEnergy(void);

void CalculateFrameworkShiftEnergyDifferenceAdsorbateVDW(void);
void CalculateFrameworkShiftEnergyDifferenceAdsorbateCharge(void);
void CalculateFrameworkShiftEnergyDifferenceCationVDW(void);
void CalculateFrameworkShiftEnergyDifferenceCationCharge(void);

REAL CalculateEnergyDifferenceFrameworkMoveVDW(int atom_id,VECTOR posA,int typeA);

void CalculateEnergyDifferenceFrameworkMoveCharge(int index);
void CalculateFrameworkEnergyDifferenceShiftedFramework(void);

int CalculateFrameworkAdsorbateVDWEnergyDifference(int m,int comp,int New,int Old,int CanUseGridNew,int CanUseGridOld);
int CalculateFrameworkCationVDWEnergyDifference(int m,int comp,int New,int Old,int CanUseGridNew,int CanUseGridOld);

int CalculateFrameworkAdsorbateChargeChargeEnergyDifference(int m,int comp,int New,int Old,int CanUseGridNew,int CanUseGridOld);
int CalculateFrameworkCationChargeChargeEnergyDifference(int m,int comp,int New,int Old,int CanUseGridNew,int CanUseGridOld);

int CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference(int m,int comp,int New,int Old);
int CalculateFrameworkCationChargeBondDipoleEnergyDifference(int m,int comp,int New,int Old);

int CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference(int m,int comp,int New,int Old);
int CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference(int m,int comp,int New,int Old);


int CalculateFrameworkAdsorbateVDWEnergyDifferenceRXCM(int reaction,REAL Lambda1,REAL Lambda2,REAL LambdaNew,int **LambdaRetraceMolecules,int direction);
int CalculateFrameworkAdsorbateVDWEnergyDifferenceNewRXCM(int reaction,REAL LambdaNew,int direction);
int CalculateFrameworkAdsorbateChargeChargeEnergyDifferenceRXCM(int reaction,REAL Lambda1,REAL Lambda2,REAL LambdaNew,int **LambdaRetraceMolecules,int direction);
int CalculateFrameworkAdsorbateChargeChargeEnergyDifferenceNewRXCM(int reaction,REAL LambdaNew,int direction);

REAL CalculateFrameworkElectrostaticPotential(POINT posA);

REAL CalculateFrameworkVDWEnergyCorrection(VECTOR* Positions,VECTOR *AnisotropicPositions,REAL *scaling);

#endif
