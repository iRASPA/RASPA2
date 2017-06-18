/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'internal_energy.h' is part of RASPA-2.0

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

#ifndef INTERNAL_ENERGY_H
#define INTERNAL_ENERGY_H

REAL GenerateBondlength(int i);
REAL GenerateBendAngle(int i);

REAL CalculateAngle(int A,int B,int C,int Iu);
REAL CalculateAngle2(int A,int B,int group_nr,int Iu);

REAL CalculateBondEnergy(int Itype,int Iu);
REAL CalculateBondEnergyAdsorbate(int m);
void CalculateBondEnergyAdsorbates(void);
REAL CalculateBondEnergyCation(int m);
void CalculateBondEnergyCations(void);

REAL CalculateUreyBradleyEnergy(int Itype,int Iu);
REAL CalculateUreyBradleyEnergyAdsorbate(int m);
void CalculateUreyBradleyEnergyAdsorbates(void);
REAL CalculateUreyBradleyEnergyCation(int m);
void CalculateUreyBradleyEnergyCations(void);


REAL CalculateBendEnergy(int Itype,int Iu);
REAL CalculateBendEnergyAdsorbate(int m);
void CalculateBendEnergyAdsorbates(void);
REAL CalculateBendEnergyCation(int m);
void CalculateBendEnergyCations(void);

REAL CalculateInversionBendEnergy(int Itype,int Iu);
REAL CalculateInversionBendEnergyAdsorbate(int m);
void CalculateInversionBendEnergyAdsorbates(void);
REAL CalculateInversionBendEnergyCation(int m);
void CalculateInversionBendEnergyCations(void);

REAL ComputeTorsionAngle(VECTOR posA,VECTOR posB,VECTOR posC, VECTOR posD);
REAL CalculateTorsionEnergy(int Itype,int Iu);
REAL CalculateTorsionEnergyAdsorbate(int m);
void CalculateTorsionEnergyAdsorbates(void);
REAL CalculateTorsionEnergyCation(int m);
void CalculateTorsionEnergyCations(void);

REAL CalculateImproperTorsionEnergy(int Itype,int Iu);
REAL CalculateImproperTorsionEnergyAdsorbate(int m);
void CalculateImproperTorsionEnergyAdsorbates(void);
REAL CalculateImproperTorsionEnergyCation(int m);
void CalculateImproperTorsionEnergyCations(void);

REAL CalculateBondBondEnergy(int Itype,int Iu);
REAL CalculateBondBondEnergyAdsorbate(int m);
void CalculateBondBondEnergyAdsorbates(void);
REAL CalculateBondBondEnergyCation(int m);
void CalculateBondBondEnergyCations(void);

REAL CalculateBondBendEnergy(int Itype,int Iu);
REAL CalculateBondBendEnergyAdsorbate(int m);
void CalculateBondBendEnergyAdsorbates(void);
REAL CalculateBondBendEnergyCation(int m);
void CalculateBondBendEnergyCations(void);

REAL CalculateBendBendEnergy(int Itype,int Iu);
REAL CalculateBendBendEnergyAdsorbate(int m);
void CalculateBendBendEnergyAdsorbates(void);
REAL CalculateBendBendEnergyCation(int m);
void CalculateBendBendEnergyCations(void);

REAL CalculateBondTorsionEnergy(int Itype,int Iu);
REAL CalculateBondTorsionEnergyAdsorbate(int m);
void CalculateBondTorsionEnergyAdsorbates(void);
REAL CalculateBondTorsionEnergyCation(int m);
void CalculateBondTorsionEnergyCations(void);

REAL CalculateBendTorsionEnergyAdsorbate(int m);
void CalculateBendTorsionEnergyAdsorbates(void);
REAL CalculateBendTorsionEnergyCation(int m);
void CalculateBendTorsionEnergyCations(void);

REAL CalculateBendTorsionEnergy(int Itype,int Iu);
REAL CalculateIntraVDWEnergyAdsorbate(int m);
void CalculateIntraVDWEnergyAdsorbates(void);
REAL CalculateIntraVDWEnergyCation(int m);
void CalculateIntraVDWEnergyCations(void);

REAL CalculateIntraChargeChargeEnergyAdsorbate(int m);
void CalculateIntraChargeChargeEnergyAdsorbates(void);
REAL CalculateIntraChargeChargeEnergyCation(int m);
void CalculateIntraChargeChargeEnergyCations(void);

REAL CalculateIntraChargeBondDipoleEnergyAdsorbate(int m);
void CalculateIntraChargeBondDipoleEnergyAdsorbates(void);
REAL CalculateIntraChargeBondDipoleEnergyCation(int m);
void CalculateIntraChargeBondDipoleEnergyCations(void);

REAL CalculateIntraBondDipoleBondDipoleEnergyAdsorbate(int m);
void CalculateIntraBondDipoleBondDipoleEnergyAdsorbates(void);
REAL CalculateIntraBondDipoleBondDipoleEnergyCation(int m);
void CalculateIntraBondDipoleBondDipoleEnergyCations(void);


void CalculateHarmonicBondConstraintEnergy(void);
void CalculateHarmonicAngleConstraintEnergy(void);
void CalculateHarmonicDihedralConstraintEnergy(void);

#endif
