/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'thermo_baro_stats.h' is part of RASPA-2.0

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
