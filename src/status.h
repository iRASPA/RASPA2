/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'status.h' is part of RASPA-2.0

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

