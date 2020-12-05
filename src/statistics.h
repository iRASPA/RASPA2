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

#ifndef STATISTICS_H
#define STATISTICS_H

#include "simulation.h"
#include "molecule.h"

#define NR_BLOCKS 5
#define STUDENTS_T_NUMBER 2.776

#define AVERAGE(sum) (sum/((REAL)NR_BLOCKS))
#define ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared) (2.776*sqrt(fabs((sum_squared/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS))/((REAL)(NR_BLOCKS-1))))

extern REAL ***WidomRosenbluthFactorAccumulated;
extern REAL ***WidomIdealGasAccumulated;
extern REAL ***WidomRosenbluthFactorCount;
extern REAL ***WidomEnergyDifferenceAccumulated;
extern REAL ***WidomEnergyFrameworkAccumulated;
extern REAL ***WidomEnergyFrameworkCount;

extern REAL ***GibbsWidomRosenbluthFactorAccumulated;
extern REAL ***GibbsWidomIdealGasAccumulated;
extern REAL ***GibbsWidomRosenbluthFactorCount;
extern REAL ***GibbsWidomEnergyDifferenceAccumulated;
extern REAL ***GibbsWidomEnergyFrameworkAccumulated;
extern REAL ***GibbsWidomEnergyFrameworkCount;

extern REAL **SurfaceAreaFrameworkAccumulated;
extern REAL ***SurfaceAreaFrameworksAccumulated;
extern REAL **SurfaceAreaCationsAccumulated;
extern REAL **SurfaceAreaCount;

REAL GetAverageWPI(int comp);

extern long long *BlockCycle;
extern REAL **BlockCount;
extern int Block;

extern REAL MeasureLambdaBelow;

extern int DensityProfile;

extern VECTOR ***PrincipleMomentsOfInertiaAccumulated;
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
void PrintIntervalStatusProduction(long long CurrentCycle,long long NumberOfCycles, FILE *FilePtr);
void PrintAverageTotalSystemEnergiesMC(FILE *FilePtr);

REAL GetAverageMolecularPressure(void);

void PrintIntervalStatusMD(long long CurrentCycle,FILE *FilePtr);

void PrintPropertyStatus(long long CurrentCycle,long long NumberOfCycles, FILE *FilePtr);
REAL GetAverageIsothermalExpansionCoefficient(void);
REAL GetAverageIsothermalCompressibilityCoefficient(void);
REAL GetAverageHeatCapacity(void);
REAL GetAverageHenryCoefficient(int comp);
REAL GetAverageRosenbluthWeight(int comp);
REAL GetAverageGibbsWidom(int comp);
REAL GetAverageGibbsInverseDensity(void);
REAL GetAverageWidom(int comp);
REAL GetAverageInverseDensity(void);
REAL GetAverageGibbsInverseDensityForComponent(int comp);

void WriteRestartStatistics(FILE *FilePtr);
void AllocateStatisticsMemory(void);
void ReadRestartStatistics(FILE *FilePtr);

#endif
