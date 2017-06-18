/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_born.h' is part of RASPA-2.0

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

#ifndef FRAMEWORK_BORN_H
#define FRAMEWORK_BORN_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <molecule.h>

void CalculateFrameworkBondBornTerm(void);
void CalculateFrameworkUreyBradleyBornTerm(void);
void CalculateFrameworkBendBornTerm(void);
void CalculateFrameworkTorsionBornTerm(void);
void CalculateFrameworkImproperTorsionBornTerm(void);

int CalculateFrameworkIntraVDWBornTerm(void);
int CalculateFrameworkIntraReplicaVDWBornTerm(void);
int CalculateFrameworkIntraChargeChargeBornTerm(void);
int CalculateFrameworkIntraReplicaChargeChargeBornTerm(void);

void FrameworkAdsorbateVDWBornTerm(void);
void FrameworkAdsorbateChargeChargeBornTerm(void);
void FrameworkCationVDWBornTerm(void);
void FrameworkCationChargeChargeBornTerm(void);

void FrameworkAdsorbateReplicaVDWBornTerm(void);
void FrameworkAdsorbateReplicaChargeChargeBornTerm(void);
void FrameworkCationReplicaVDWBornTerm(void);
void FrameworkCationReplicaChargeChargeBornTerm(void);
#endif
