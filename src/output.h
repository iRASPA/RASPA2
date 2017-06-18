/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'output.h' is part of RASPA-2.0

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

#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include "simulation.h"

extern FILE **OutputFilePtr;

void OpenOutputFile(void);
void CloseOutputFile(void);

void PrintPreSimulationStatus(void);
void PrintPreSimulationStatusCurrentSystem(int system);

void PrintPostSimulationStatus(void);

void PrintEnergyStatus(FILE *FilePtr,char *string);
void PrintEnergyStatusToStdOut();

void PrintEnergyDriftStatus(FILE *FilePtr);

void PrintRestartFile(void);

void WriteBinaryRestartFiles(void);

void ReadRestartOutput(FILE* FilePtr);
void WriteRestartOutput(FILE* FilePtr);

#endif
