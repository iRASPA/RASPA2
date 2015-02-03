/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2015 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'statistics.h' is part of RASPA-2.0

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

#ifndef STATISTICS_H
#define STATISTICS_H

#include "simulation.h"
#include "molecule.h"

#define NR_BLOCKS 5

extern REAL ***WidomRosenbluthFactorAverage;
extern REAL ***WidomRosenbluthFactorCount;

extern REAL ***WidomEnergyDifferenceAverage;

extern REAL ***WidomEnergyFrameworkAverage;
extern REAL ***WidomEnergyFrameworkCount;

extern REAL **SurfaceAreaFrameworkAverage;
extern REAL ***SurfaceAreaFrameworksAverage;
extern REAL **SurfaceAreaCationsAverage;
extern REAL **SurfaceAreaCount;

REAL GetAverageWPI(int comp);

extern long long *BlockCycle;
extern REAL **BlockCount;
extern int Block;

extern int DensityProfile;

extern VECTOR ***PrincipleMomentsOfInertiaAverage;
extern REAL ***PrincipleMomentsOfInertiaCount;

REAL GetAverageProperty(REAL **Property);
REAL GetAverageComponentProperty(REAL ***Property,int comp);

REAL GetAverageTemperature(void);
REAL GetAverageAdsorbatesTemperature(void);
REAL GetAverageCationsTemperature(void);
REAL GetAverageFrameworkTemperature(void);
REAL_MATRIX3x3 GetAveragePressureTensor(void);

void PrintProperty(FILE *FilePtr,char *string,char *units,REAL conv_factor,REAL **Property);
void PrintEnergies(FILE *FilePtr,char *string,char *units,REAL conv_factor,REAL **Property,
                   REAL **PropertyVDW,REAL **PropertyCoulomb);

void InitializesEnergiesAllSystems(void);
void InitializesEnergiesCurrentSystem(void);
void InitializesEnergyAveragesAllSystems(void);
void UpdateEnergyAveragesCurrentSystem(void);
void PrintIntervalStatusInit(long long CurrentCycle,long long NumberOfCycles, FILE *FilePtr);
void PrintIntervalStatusEquilibration(long long CurrentCycle,long long NumberOfCycles, FILE *FilePtr);
void PrintIntervalStatus(long long CurrentCycle,long long NumberOfCycles, FILE *FilePtr);
void PrintAverageTotalSystemEnergiesMC(FILE *FilePtr);

REAL GetAverageMolecularPressure(void);

void PrintIntervalStatusMD(long long CurrentCycle,FILE *FilePtr);

void PrintPropertyStatus(long long CurrentCycle,long long NumberOfCycles, FILE *FilePtr);
REAL GetAverageIsothermalExpansionCoefficient(void);
REAL GetAverageIsothermalCompressibilityCoefficient(void);
REAL GetAverageHeatCapacity(void);
REAL GetAverageHenryCoefficient(int comp);
REAL GetAverageRosenbluthWeight(int comp);

void WriteRestartStatistics(FILE *FilePtr);
void AllocateStatisticsMemory(void);
void ReadRestartStatistics(FILE *FilePtr);

#endif
