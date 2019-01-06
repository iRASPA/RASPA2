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

#ifndef INTER_ENERGY_H
#define INTER_ENERGY_H

#include "utils.h"
#include "simulation.h"

int CalculateInterVDWEnergyDifferenceAdsorbateNewRXCM(int reaction,REAL LambdaNew,int direction);
int CalculateInterVDWEnergyDifferenceAdsorbateRXCM(int reaction,REAL Lambda1,REAL Lambda2,REAL LambdaNew,int **LambdaRetraceMolecules,int direction);

int CalculateInterChargeChargeEnergyDifferenceAdsorbateNewRXCM(int reaction,REAL LambdaNew,int direction);
int CalculateInterChargeChargeEnergyDifferenceAdsorbateRXCM(int reaction,REAL Lambda1,REAL Lambda2,REAL LambdaNew,int **LambdaRetraceMolecules,int direction);

REAL CalculateInterVDWSelfEnergyCorrectionNew(void);
REAL CalculateInterVDWSelfEnergyCorrectionAdsorbateOld(int mol);
REAL CalculateInterVDWSelfEnergyCorrectionCationOld(int mol);
REAL CalculateInterChargeChargeSelfEnergyCorrectionNew(void);
REAL CalculateInterChargeChargeSelfEnergyCorrectionAdsorbateOld(int mol);
REAL CalculateInterChargeChargeSelfEnergyCorrectionCationOld(int mol);

REAL CalculateInterVDWEnergyAdsorbateAtPosition(POINT posA,int typeA,int exclude,REAL scaling);
int CalculateInterChargeEnergyAdsorbateAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,int exclude,REAL scaling);
int CalculateInterBondDipoleEnergyAdsorbateAtPosition(POINT posA1,POINT posA2,REAL DipoleMagnitudeA,REAL *UChargeBondDipole,REAL *UBondDipoleBondDipole,int exclude);

REAL CalculateInterVDWEnergyCationAtPosition(POINT posA,int typeA,int exclude,REAL scaling);
int CalculateInterChargeEnergyCationAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,int exclude,REAL scaling);
int CalculateInterBondDipoleEnergyCationAtPosition(POINT posA1,POINT posA2,REAL DipoleMagnitudeA,REAL *UChargeBondDipole,REAL *UBondDipoleBondDipole,int exclude);

int CalculateTotalInterVDWEnergy(void);
int CalculateTotalInterChargeChargeCoulombEnergy(void);
int CalculateTotalInterReplicaChargeChargeCoulombEnergy(void);
int CalculateTotalInterChargeBondDipoleCoulombEnergy(void);
int CalculateTotalInterBondDipoleBondDipoleCoulombEnergy(void);

int CalculateInterVDWEnergyDifferenceAdsorbate(int m,int comp,int New,int Old);
int CalculateInterVDWEnergyDifferenceCation(int m,int comp,int New,int Old);

int CalculateInterChargeChargeEnergyDifferenceAdsorbate(int m,int comp,int New,int Old);
int CalculateInterChargeChargeEnergyDifferenceCation(int m,int comp,int New,int Old);

int CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(int m,int comp,int New,int Old);
int CalculateInterChargeBondDipoleEnergyDifferenceCation(int m,int comp,int New,int Old);

int CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(int m,int comp,int New,int Old);
int CalculateInterBondDipoleBondDipoleEnergyDifferenceCation(int m,int comp,int New,int Old);

REAL CalculateInterChargeElectrostaticPotential(POINT posA);

REAL CalculateInterVDWEnergyCorrectionAdsorbate(VECTOR* Positions,VECTOR *AnisotropicPositions,int exclude);
REAL CalculateInterVDWEnergyCorrectionCation(VECTOR* Positions,VECTOR *AnisotropicPositions,int exclude);

#endif
