/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'internal_force.h' is part of RASPA-2.0

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

#ifndef INTERNAL_FORCE_H
#define INTERNAL_FORCE_H

REAL ReturnBondDistance(VECTOR posA,VECTOR posB);
REAL ReturnBendAngle(VECTOR posA,VECTOR posB,VECTOR posC);
REAL ReturnDihedralAngle(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD);
REAL ReturnInversionBendAngle(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD);
REAL ReturnOutOfPlaneDistance(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD);


void ReturnWilsonVectorsBond(VECTOR posA,VECTOR posB,VECTOR *wa,VECTOR *wb);

void ReturnWilsonVectorsBend(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR *wa,VECTOR *wb,VECTOR *wc);
void ReturnWilsonVectorsTorsion(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR *wa,VECTOR *wb,VECTOR *wc,VECTOR *wd);
void ReturnWilsonVectorsInversionBend(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR *wa,VECTOR *wb,VECTOR *wc,VECTOR *wd);

void CalculateAdsorbateBondForce(int m);
void CalculateAdsorbateUreyBradleyForce(int m);
void CalculateAdsorbateBendForce(int m);
void CalculateAdsorbateInversionBendForce(int m);
void CalculateAdsorbateTorsionForce(int m);
void CalculateAdsorbateImproperTorsionForce(int m);
void CalculateAdsorbateBondBondForce(int m);
void CalculateAdsorbateBondBendForce(int m);
void CalculateAdsorbateBendBendForce(int m);
void CalculateAdsorbateBendTorsionForce(int m);
void CalculateAdsorbateBondTorsionForce(int m);
void CalculateAdsorbateIntraVDWForce(int m);
void CalculateAdsorbateIntraChargeChargeForce(int m);
void CalculateAdsorbateIntraChargeBondDipoleForce(int m);
void CalculateAdsorbateIntraBondDipoleBondDipoleForce(int m);

void CalculateCationBondForce(int m);
void CalculateCationUreyBradleyForce(int m);
void CalculateCationBendForce(int m);
void CalculateCationInversionBendForce(int m);
void CalculateCationTorsionForce(int m);
void CalculateCationImproperTorsionForce(int m);
void CalculateCationBondBondForce(int m);
void CalculateCationBondBendForce(int m);
void CalculateCationBendBendForce(int m);
void CalculateCationBendTorsionForce(int m);
void CalculateCationBondTorsionForce(int m);
void CalculateCationIntraVDWForce(int m);
void CalculateCationIntraChargeChargeForce(int m);
void CalculateCationIntraChargeBondDipoleForce(int m);
void CalculateCationIntraBondDipoleBondDipoleForce(int m);

void CalculateHarmonicBondConstraintForce(void);
void CalculateHarmonicAngleConstraintForce(void);
void CalculateHarmonicDihedralConstraintForce(void);

#endif
