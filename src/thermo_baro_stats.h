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

#ifndef THERMO_BARO_STATS_H
#define THERMO_BARO_STATS_H

#define MAX_NUMBER_OF_ISOTHERM_PRESSURES 100
extern int NumberOfIsothermPressures;
extern int CurrentIsothermPressure;

#define MAXIMUM_LENGTH_THERMOSTATS 10
#define MAXIMUM_LENGTH_BAROSTATS 10

typedef struct thermo_barostats
{
  int ThermostatChainLength;
  int BarostatChainLength;
  REAL time_scale_parameter_thermostat;                          // in picoseconds
  REAL time_scale_parameter_barostat;                            // in picoseconds
  REAL *ExternalTemperature;               // in Kelvin
  REAL **ExternalPressure;                 // in Pascal
  int UseExternalStress;
  REAL_MATRIX3x3 *ExternalStress;          // in Pascal
  REAL_MATRIX3x3 *ExternalPressureTensor;  // in Pascal
  VECTOR *ExternalSurfaceTension;          // in Pascal
  int NumberOfYoshidaSuzukiSteps;
  int NumberOfRespaSteps;

  int ThermostatFramework;   // TODO
  int ThermostatAdsorbates;
  int ThermostatCations;
} THERMO_BAROSTATS;

extern THERMO_BAROSTATS therm_baro_stats;

void AdjustCoreShellVelocities(void);

void InitializeNoseHooverCurrentSystem(void);
void InitializeNoseHooverAllSystems(void);
void ComputeNoseHooverEnergySystem(void);

void InitializeBoxVelocities(void);

void NoseHooverNPT(void);
void NoseHooverNPTPR(void);
void UpdatePositions(void);
void UpdateVelocities(void);

void UpdateCellVelocity(void);

REAL_MATRIX3x3 GetKineticStressTensor(void);
REAL GetRotationalKineticEnergy(void);
REAL GetCoreShellTemperature(void);

REAL GetTranslationKineticEnergy(void);
REAL GetTranslationKineticEnergyFramework(void);
REAL GetTranslationKineticEnergyAdsorbates(void);
REAL GetTranslationKineticEnergyAdsorbates2(void);
REAL GetTranslationKineticEnergyCations(void);
REAL GetCellKineticEnergy(void);
REAL GetCellTemperature(void);

void ComputeNoseHooverEnergySystem(void);


void NoseHooverNVTFramework(void);
void NoseHooverNVTAdsorbates(void);
void NoseHooverNVTCations(void);
void NoseHooverNVTRotation(void);

void NoseHooverNPTPR(void);
void UpdatePositionsVelocitiesNPTPR(void);

void NoseHooverNPHPR(void);

void WriteRestartThermoBarostats(FILE *FilePtr);
void AllocateThermoBaroStatMemory(void);
void ReadRestartThermoBarostats(FILE *FilePtr);

#endif
