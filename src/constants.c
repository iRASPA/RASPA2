/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'constants.c' is part of RASPA-2.0

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "utils.h"
#include "ewald.h"

int UseReducedUnits;

REAL LENGTH_UNIT;
REAL TIME_UNIT;
REAL MASS_UNIT;
REAL CHARGE_UNIT;

REAL MASS_CONVERSION_FACTOR;
REAL VOLUME_CONVERSION_FACTOR;
REAL DENSITY_CONVERSION_FACTOR;
REAL ENERGY_CONVERSION_FACTOR;
REAL FORCE_CONVERSION_FACTOR;
REAL PRESSURE_CONVERSION_FACTOR;
REAL POSITION_CONVERSION_FACTOR;
REAL TIME_CONVERSION_FACTOR;
REAL DIFFUSION_CONVERSION_FACTOR;
REAL VELOCITY_CONVERSION_FACTOR;
REAL ACCELERATION_CONVERSION_FACTOR;
REAL DIPOLE_MOMENT_CONVERSION_FACTOR;
REAL DEBYE_CONVERSION_FACTOR;
REAL ELECTRIC_POTENTIAL_CONVERSION_FACTOR;
REAL POLARIZABILITY;
REAL ELECTRIC_FIELD_CONVERSION_FACTOR;
REAL K_B;
REAL COULOMBIC_CONVERSION_FACTOR;
REAL DIELECTRIC_CONSTANT_CONVERSION_FACTOR;
REAL ISOTHERMAL_COMPRESSIBILITY_CONVERSION_FACTOR;
REAL HEAT_CAPACITY_CONVERSION_FACTOR;
REAL VOLUMETRIC_EXPANSION_COEFFICIENT_CONVERSION_FACTOR;
REAL FH_CONVERSION_FACTOR;

REAL ENERGY_TO_KELVIN;
REAL KELVIN_TO_ENERGY;
REAL ENERGY_TO_KJ_PER_MOL;
REAL ENERGY_TO_EV;
REAL ENERGY_TO_KCAL_PER_MOL;
REAL ENERGY_TO_KCAL15_PER_MOL;
REAL KCAL_PER_MOL_TO_ENERGY;
REAL KCAL15_PER_MOL_TO_ENERGY;

REAL TO_WAVENUMBERS;
REAL TO_THZ;


void SetSimulationUnits(void)
{
  // mutual consistent basic set of units
  if(UseReducedUnits)
  {
    LENGTH_UNIT=1.0;
    TIME_UNIT=1.0;
    MASS_UNIT=1.0;
    CHARGE_UNIT=1.0;

    // derived units and their conversion factors
    MASS_CONVERSION_FACTOR=1.0;
    VOLUME_CONVERSION_FACTOR=1.0;
    DENSITY_CONVERSION_FACTOR=1.0;
    ENERGY_CONVERSION_FACTOR=1.0;
    FORCE_CONVERSION_FACTOR=1.0;
    PRESSURE_CONVERSION_FACTOR=1.0;
    POSITION_CONVERSION_FACTOR=1.0;
    TIME_CONVERSION_FACTOR=1.0;
    DIFFUSION_CONVERSION_FACTOR=1.0;
    VELOCITY_CONVERSION_FACTOR=1.0;
    ACCELERATION_CONVERSION_FACTOR=1.0;
    DIPOLE_MOMENT_CONVERSION_FACTOR=1.0;
    DEBYE_CONVERSION_FACTOR=1.0;
    ELECTRIC_POTENTIAL_CONVERSION_FACTOR=1.0;
    POLARIZABILITY=1.0;
    ELECTRIC_FIELD_CONVERSION_FACTOR=1.0;
    K_B=1.0;

    COULOMBIC_CONVERSION_FACTOR=1.0;
    DIELECTRIC_CONSTANT_CONVERSION_FACTOR=1.0;
    ISOTHERMAL_COMPRESSIBILITY_CONVERSION_FACTOR=1.0;
    HEAT_CAPACITY_CONVERSION_FACTOR=1.0;
    VOLUMETRIC_EXPANSION_COEFFICIENT_CONVERSION_FACTOR=1.0;
    FH_CONVERSION_FACTOR=1.0;

    ENERGY_TO_KELVIN=1.0;
    KELVIN_TO_ENERGY=1.0;
    ENERGY_TO_KJ_PER_MOL=1.0;
    ENERGY_TO_EV=1.0;
    ENERGY_TO_KCAL_PER_MOL=1.0;
    ENERGY_TO_KCAL15_PER_MOL=1.0;
    KCAL_PER_MOL_TO_ENERGY=1.0;
    KCAL15_PER_MOL_TO_ENERGY=1.0;

    TO_WAVENUMBERS=1.0;
  }
  else
  {
    LENGTH_UNIT=ANGSTROM;
    TIME_UNIT=PICO_SECOND;
    MASS_UNIT=ATOMIC_MASS_UNIT;
    CHARGE_UNIT=ELECTRONIC_CHARGE_UNIT;

    // derived units and their conversion factors
    MASS_CONVERSION_FACTOR=ATOMIC_MASS_UNIT;                                     // kg
    VOLUME_CONVERSION_FACTOR=(CUBE(LENGTH_UNIT));                                // m^-3
    DENSITY_CONVERSION_FACTOR=(MASS_CONVERSION_FACTOR/VOLUME_CONVERSION_FACTOR); // kg/m^3
    ENERGY_CONVERSION_FACTOR=((MASS_UNIT)*SQR(LENGTH_UNIT)/SQR(TIME_UNIT));      // J
    FORCE_CONVERSION_FACTOR=(ENERGY_CONVERSION_FACTOR/LENGTH_UNIT);              // N (=m/s^2)
    PRESSURE_CONVERSION_FACTOR=(MASS_UNIT/(LENGTH_UNIT*SQR(TIME_UNIT)));         // Pa(=N/m^2)
    POSITION_CONVERSION_FACTOR=(LENGTH_UNIT);                                    // m
    TIME_CONVERSION_FACTOR=(TIME_UNIT);                                          // s
    DIFFUSION_CONVERSION_FACTOR=(SQR(LENGTH_UNIT)/TIME_UNIT);                    // m^2/s
    VELOCITY_CONVERSION_FACTOR=(LENGTH_UNIT/TIME_UNIT);                          // m/s
    ACCELERATION_CONVERSION_FACTOR=(SQR(LENGTH_UNIT)/TIME_UNIT);                 // m^2/s
    DIPOLE_MOMENT_CONVERSION_FACTOR=(CHARGE_UNIT*LENGTH_UNIT);                   // C.m
    DEBYE_CONVERSION_FACTOR=((CHARGE_UNIT*LENGTH_UNIT)/(DEBYE));                 // C.m
    ELECTRIC_POTENTIAL_CONVERSION_FACTOR=(ENERGY_CONVERSION_FACTOR/CHARGE_UNIT); // V
    POLARIZABILITY=(SQR(CHARGE_UNIT)*SQR(LENGTH_UNIT)/ENERGY_CONVERSION_FACTOR);
    ELECTRIC_FIELD_CONVERSION_FACTOR=(ELECTRIC_POTENTIAL_CONVERSION_FACTOR/LENGTH_UNIT);     // V
    K_B=(BOLTZMANN_CONSTANT/ENERGY_CONVERSION_FACTOR);

    COULOMBIC_CONVERSION_FACTOR=(((SQR(CHARGE_UNIT)/(4.0*M_PI*ELECTRIC_CONSTANT*LENGTH_UNIT*ENERGY_CONVERSION_FACTOR)))/(DielectricConstantOfTheMedium));
    DIELECTRIC_CONSTANT_CONVERSION_FACTOR=(SQR(TIME_UNIT)*SQR(CHARGE_UNIT)/(ATOMIC_MASS_UNIT*CUBE(LENGTH_UNIT)));
    ISOTHERMAL_COMPRESSIBILITY_CONVERSION_FACTOR=(1.0/PRESSURE_CONVERSION_FACTOR);
    HEAT_CAPACITY_CONVERSION_FACTOR=(ENERGY_CONVERSION_FACTOR*AVOGADRO_CONSTANT);     // J/mol/K
    VOLUMETRIC_EXPANSION_COEFFICIENT_CONVERSION_FACTOR=(1.0);                         // 1/K
    FH_CONVERSION_FACTOR=(SQR(PLANCK_CONSTANT/(2.0*M_PI)))/(24.0*ATOMIC_MASS_UNIT*BOLTZMANN_CONSTANT*SQR(LENGTH_UNIT));

    // 1 calorie (International Table) = 4.1868 J
    // 1 calorie (thermochemical) = 4.184 J
    // 1 calorie (15°C) = 4.1855 J
    ENERGY_TO_KELVIN=((ENERGY_CONVERSION_FACTOR*AVOGADRO_CONSTANT)/MOLAR_GAS_CONSTANT);
    KELVIN_TO_ENERGY=(MOLAR_GAS_CONSTANT/(ENERGY_CONVERSION_FACTOR*AVOGADRO_CONSTANT));
    ENERGY_TO_KJ_PER_MOL=((ENERGY_CONVERSION_FACTOR*AVOGADRO_CONSTANT)/1000.0);
    ENERGY_TO_EV=(ENERGY_TO_KELVIN/11604.23);
    ENERGY_TO_KCAL_PER_MOL=((ENERGY_CONVERSION_FACTOR*AVOGADRO_CONSTANT)/4184.0);
    ENERGY_TO_KCAL15_PER_MOL=((ENERGY_CONVERSION_FACTOR*AVOGADRO_CONSTANT)/4185.5);
    KCAL_PER_MOL_TO_ENERGY=(4184.0/(ENERGY_CONVERSION_FACTOR*AVOGADRO_CONSTANT));
    KCAL15_PER_MOL_TO_ENERGY=(4185.5/(ENERGY_CONVERSION_FACTOR*AVOGADRO_CONSTANT));

    TO_WAVENUMBERS=(1.0/(100.0*2.0*M_PI*SPEED_OF_LIGHT*TIME_UNIT));
    TO_THZ=(1e-12/(2.0*M_PI*TIME_UNIT));
  }
}

static int versionNumber=1;

void WriteRestartConstants(FILE *FilePtr)
{
  REAL Check;

  fwrite(&versionNumber,sizeof(int),1,FilePtr);
  fwrite(&LENGTH_UNIT,sizeof(REAL),1,FilePtr);
  fwrite(&TIME_UNIT,sizeof(REAL),1,FilePtr);
  fwrite(&MASS_UNIT,sizeof(REAL),1,FilePtr);
  fwrite(&CHARGE_UNIT,sizeof(REAL),1,FilePtr);

  fwrite(&MASS_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&VOLUME_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&DENSITY_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&ENERGY_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&FORCE_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&PRESSURE_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&POSITION_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&TIME_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&DIFFUSION_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&VELOCITY_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&ACCELERATION_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&DIPOLE_MOMENT_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&DEBYE_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&ELECTRIC_POTENTIAL_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&POLARIZABILITY,sizeof(REAL),1,FilePtr);
  fwrite(&ELECTRIC_FIELD_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&K_B,sizeof(REAL),1,FilePtr);

  fwrite(&COULOMBIC_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&DIELECTRIC_CONSTANT_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&ISOTHERMAL_COMPRESSIBILITY_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&HEAT_CAPACITY_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&VOLUMETRIC_EXPANSION_COEFFICIENT_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fwrite(&FH_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);

  fwrite(&ENERGY_TO_KELVIN,sizeof(REAL),1,FilePtr);
  fwrite(&KELVIN_TO_ENERGY,sizeof(REAL),1,FilePtr);
  fwrite(&ENERGY_TO_KJ_PER_MOL,sizeof(REAL),1,FilePtr);
  fwrite(&ENERGY_TO_EV,sizeof(REAL),1,FilePtr);
  fwrite(&ENERGY_TO_KCAL_PER_MOL,sizeof(REAL),1,FilePtr);
  fwrite(&ENERGY_TO_KCAL15_PER_MOL,sizeof(REAL),1,FilePtr);
  fwrite(&KCAL_PER_MOL_TO_ENERGY,sizeof(REAL),1,FilePtr);
  fwrite(&KCAL15_PER_MOL_TO_ENERGY,sizeof(REAL),1,FilePtr);

  fwrite(&TO_WAVENUMBERS,sizeof(REAL),1,FilePtr);
  fwrite(&TO_THZ,sizeof(REAL),1,FilePtr);

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void ReadRestartConstants(FILE *FilePtr)
{
  REAL Check;
  int readversionNumber=0;

  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumber)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  fread(&LENGTH_UNIT,sizeof(REAL),1,FilePtr);
  fread(&TIME_UNIT,sizeof(REAL),1,FilePtr);
  fread(&MASS_UNIT,sizeof(REAL),1,FilePtr);
  fread(&CHARGE_UNIT,sizeof(REAL),1,FilePtr);

  fread(&MASS_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&VOLUME_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&DENSITY_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&ENERGY_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&FORCE_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&PRESSURE_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&POSITION_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&TIME_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&DIFFUSION_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&VELOCITY_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&ACCELERATION_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&DIPOLE_MOMENT_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&DEBYE_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&ELECTRIC_POTENTIAL_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&POLARIZABILITY,sizeof(REAL),1,FilePtr);
  fread(&ELECTRIC_FIELD_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&K_B,sizeof(REAL),1,FilePtr);

  fread(&COULOMBIC_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&DIELECTRIC_CONSTANT_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&ISOTHERMAL_COMPRESSIBILITY_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&HEAT_CAPACITY_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&VOLUMETRIC_EXPANSION_COEFFICIENT_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);
  fread(&FH_CONVERSION_FACTOR,sizeof(REAL),1,FilePtr);

  fread(&ENERGY_TO_KELVIN,sizeof(REAL),1,FilePtr);
  fread(&KELVIN_TO_ENERGY,sizeof(REAL),1,FilePtr);
  fread(&ENERGY_TO_KJ_PER_MOL,sizeof(REAL),1,FilePtr);
  fread(&ENERGY_TO_EV,sizeof(REAL),1,FilePtr);
  fread(&ENERGY_TO_KCAL_PER_MOL,sizeof(REAL),1,FilePtr);
  fread(&ENERGY_TO_KCAL15_PER_MOL,sizeof(REAL),1,FilePtr);
  fread(&KCAL_PER_MOL_TO_ENERGY,sizeof(REAL),1,FilePtr);
  fread(&KCAL15_PER_MOL_TO_ENERGY,sizeof(REAL),1,FilePtr);

  fread(&TO_WAVENUMBERS,sizeof(REAL),1,FilePtr);
  fread(&TO_THZ,sizeof(REAL),1,FilePtr);

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartConstants)\n");
    ContinueAfterCrash=FALSE;
  }
}

