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
