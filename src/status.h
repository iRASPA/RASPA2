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

#ifndef STATUS_H
#define STATUS_H

extern int PrintFrameworkBondStatus;
extern int PrintFrameworkUreyBradleyStatus;
extern int PrintFrameworkBendStatus;
extern int PrintFrameworkInversionBendStatus;
extern int PrintFrameworkTorsionStatus;
extern int PrintFrameworkImproperTorsionStatus;
extern int PrintFrameworkBondBondStatus;
extern int PrintFrameworkBondBendStatus;
extern int PrintFrameworkBendBendStatus;
extern int PrintFrameworkBendTorsionStatus;
extern int PrintFrameworkIntraVDWStatus;
extern int PrintFrameworkIntraChargeChargeStatus;
extern int PrintFrameworkIntraChargeBondDipoleStatus;
extern int PrintFrameworkIntraBondDipoleBondDipoleStatus;

extern int PrintAdsorbateBondStatus;
extern int PrintAdsorbateUreyBradleyStatus;
extern int PrintAdsorbateBendStatus;
extern int PrintAdsorbateTorsionStatus;
extern int PrintAdsorbateImproperTorsionStatus;
extern int PrintAdsorbateBondBondStatus;
extern int PrintAdsorbateBondBendStatus;
extern int PrintAdsorbateBendBendStatus;
extern int PrintAdsorbateBondTorsionStatus;
extern int PrintAdsorbateBendTorsionStatus;
extern int PrintAdsorbateIntraVDWStatus;
extern int PrintAdsorbateIntraChargeChargeStatus;
extern int PrintAdsorbateIntraChargeBondDipoleStatus;
extern int PrintAdsorbateIntraBondDipoleBondDipoleStatus;

extern int PrintCationBondStatus;
extern int PrintCationUreyBradleyStatus;
extern int PrintCationBendStatus;
extern int PrintCationTorsionStatus;
extern int PrintCationImproperTorsionStatus;
extern int PrintCationBondBondStatus;
extern int PrintCationBondBendStatus;
extern int PrintCationBendBendStatus;
extern int PrintCationBondTorsionStatus;
extern int PrintCationBendTorsionStatus;
extern int PrintCationIntraVDWStatus;
extern int PrintCationIntraChargeChargeStatus;
extern int PrintCationIntraChargeBondDipoleStatus;
extern int PrintCationIntraBondDipoleBondDipoleStatus;

extern int PrintInterVDWStatus;
extern int PrintInterChargeChargeStatus;
extern int PrintInterChargeBondDipoleStatus;
extern int PrintInterBondDipoleBondDipoleStatus;

extern int PrintFrameworkAdsorbateVDWStatus;
extern int PrintFrameworkAdsorbateChargeChargeStatus;
extern int PrintFrameworkAdsorbateChargeBondDipoleStatus;
extern int PrintFrameworkAdsorbateBondDipoleBondDipoleStatus;

extern int PrintFrameworkCationVDWStatus;
extern int PrintFrameworkCationChargeChargeStatus;
extern int PrintFrameworkCationChargeBondDipoleStatus;
extern int PrintFrameworkCationBondDipoleBondDipoleStatus;


void SetPrintStatusToFalse(void);
void Status(void);

#endif

