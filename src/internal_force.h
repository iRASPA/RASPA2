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
