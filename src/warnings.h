/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'warnings.h' is part of RASPA-2.0

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

#ifndef WARNINGS_H
#define WARNINGS_H

#include "constants.h"
#include "utils.h"
#include "vector.h"

#define MAX_NUMBER_OF_WARNINGS 10
#define MAX_NUMBER_OF_WARNING_ARGUMENTS 100

extern int *NumberOfWarnings;
extern int (*Warnings)[MAX_NUMBER_OF_WARNINGS];
extern int (*NumberOfWarningValues)[MAX_NUMBER_OF_WARNINGS];
extern char (*WarningValues)[MAX_NUMBER_OF_WARNINGS][MAX_NUMBER_OF_WARNING_ARGUMENTS][32];

enum {NET_SYSTEM_CHARGE,LOWENSTEIN_RULE_NOT_OBEYED,UNIT_CELL,OMITTED_HOST_ADSORBATE_VDW_INTERACTIONS,
      OMITTED_HOST_CATION_VDW_INTERACTIONS,OMITTED_ADSORBATE_ADSORBATE_VDW_INTERACTIONS,
      OMITTED_ADSORBATE_CATION_VDW_INTERACTIONS,OMITTED_CATION_CATION_VDW_INTERACTIONS,
      REINSERTION_MOVE_IONS,GRID_ERROR_ENERGY,GRID_ERROR_FORCE,PARTIAL_PRESSURE_NOT_SET,
      MISSING_INTERACTION,ENERGY_DRIFT};

void CreateWarnings(void);
void CheckForErrors(void);
void PrintWarningStatus(void);

void WriteRestartWarnings(FILE *FilePtr);
void AllocateWarningsMemory(void);
void ReadRestartWarnings(FILE *FilePtr);

#endif
