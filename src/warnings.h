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
