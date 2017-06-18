/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'constants.h' is part of RASPA-2.0

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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <stdio.h>

#define REAL double

#define TWO_PI (2.0*M_PI)
#define FOUR_PI (4.0*M_PI)
#define EIGHT_PI (8.0*M_PI)

// Zero charge is set to a very small value. All electrostatic properties depend on
// a particle having a non-zero charge.
#define ZERO_CHARGE 1e-20

// physical constants
#define ANGSTROM 1e-10
#define NANO_SECOND 1e-9
#define PICO_SECOND 1e-12
#define FEMTO_SECOND 1e-15
#define ATOMIC_MASS_UNIT 1.6605402e-27        // kg
#define ELECTRONIC_CHARGE_UNIT 1.60217733e-19 // C/particle
#define MOLAR_GAS_CONSTANT 8.314464919         // J mol^-1 K^-1
#define BOLTZMANN_CONSTANT 1.380650324e-23    // J K^-1
#define SPEED_OF_LIGHT 299792458.0             // m s^-1
#define AVOGADRO_CONSTANT 6.0221419947e23     // mol^-1
#define ELECTRIC_CONSTANT 8.8541878176e-12     // C^2/(N.m^2)
#define PLANCK_CONSTANT 6.6260687652e-34      // J.s

#define DIELECTRIC_CONSTANT_VACUUM (8.8541878176e-12)   // F/m
#define DEBYE (3.335640952e-30)                // C.m

#define JOULE_TO_CAL 4.184
#define CAL_TO_JOULE (1.0/JOULE_TO_CAL)
#define EV_TO_KJ_MOL (96.48534)
#define EV_TO_KELVIN (EV_TO_KJ_MOL*1000.0/MOLAR_GAS_CONSTANT)

#define MDYNE_PER_ANGSTROM_TO_KCAL_PER_MOL_PER_ANGSTROM_SQUARED (AVOGADRO_CONSTANT/(JOULE_TO_CAL*1e21))
#define MDYNE_PER_ANGSTROM_TO_KJ_PER_MOL_PER_ANGSTROM_SQUARED (AVOGADRO_CONSTANT/(1e21))
#define MDYNE_PER_ANGSTROM_TO_KELVIN_PER_ANGSTROM_SQUARED (AVOGADRO_CONSTANT/(1e18*MOLAR_GAS_CONSTANT))
#define MDYNE_ANGSTROM_PER_RAD_TO_KCAL_PER_MOL_PER_DEGREE_SQUARED (MDYNE_PER_ANGSTROM_TO_KCAL_PER_MOL_PER_ANGSTROM_SQUARED/SQR(RAD2DEG))

#define PA_TO_ATM (1.0/101325.0)
#define PA_TO_BAR (0.00001)
#define PA_TO_TORR 0.0075
#define KPA_TO_ATM (1.0/101.325)
#define KPA_TO_BAR 0.01
#define KPA_TO_TORR 7.5
#define BAR_TO_TORR 750.062
#define BAR_TO_ATM 0.9869
#define BAR_TO_KPA 100.0
#define ATM_TO_TORR 760
#define ATM_TO_BAR (1.0/0.9869)
#define ATM_TO_KPA 101.325
#define ATM_TO_PA 101325.0

#define KELVIN_TO_CELSIUS(x)  ((x)-273.15)
#define CELSIUS_TO_KELVIN(x)  ((x)+273.15)

#define KELVIN_TO_KJ_PER_MOL   (MOLAR_GAS_CONSTANT/1000.0)
#define KELVIN_TO_KCAL_PER_MOL (MOLAR_GAS_CONSTANT/4184.0)
#define KELVIN_TO_KCAL15_PER_MOL (MOLAR_GAS_CONSTANT/4185.5)
#define KELVIN_TO_J_PER_MOL   (MOLAR_GAS_CONSTANT)
#define KELVIN_TO_CAL_PER_MOL (MOLAR_GAS_CONSTANT/4.184)
#define KELVIN_TO_CAL15_PER_MOL (MOLAR_GAS_CONSTANT/4.1855)
#define KELVIN_TO_EV (MOLAR_GAS_CONSTANT/(EV_TO_KJ_MOL*1000.0))
#define KJ_PER_MOL_TO_KELVIN   (1000.0/MOLAR_GAS_CONSTANT)
#define KCAL_PER_MOL_TO_KELVIN (4184.0/MOLAR_GAS_CONSTANT)
#define KCAL15_PER_MOL_TO_KELVIN (4185.5/MOLAR_GAS_CONSTANT)
#define CAL_TO_J (4.184)
#define CAL15_TO_J (4.1855)
#define J_TO_CAL (1.0/4.184)
#define J_TO_CAL15 (1.0/4.1855)

#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

extern int UseReducedUnits;

extern REAL LENGTH_UNIT;
extern REAL TIME_UNIT;
extern REAL MASS_UNIT;
extern REAL CHARGE_UNIT;

extern REAL MASS_CONVERSION_FACTOR;
extern REAL VOLUME_CONVERSION_FACTOR;
extern REAL DENSITY_CONVERSION_FACTOR;
extern REAL ENERGY_CONVERSION_FACTOR;
extern REAL FORCE_CONVERSION_FACTOR;
extern REAL PRESSURE_CONVERSION_FACTOR;
extern REAL POSITION_CONVERSION_FACTOR;
extern REAL TIME_CONVERSION_FACTOR;
extern REAL DIFFUSION_CONVERSION_FACTOR;
extern REAL VELOCITY_CONVERSION_FACTOR;
extern REAL ACCELERATION_CONVERSION_FACTOR;
extern REAL DIPOLE_MOMENT_CONVERSION_FACTOR;
extern REAL DEBYE_CONVERSION_FACTOR;
extern REAL ELECTRIC_POTENTIAL_CONVERSION_FACTOR;
extern REAL POLARIZABILITY;
extern REAL ELECTRIC_FIELD_CONVERSION_FACTOR;
extern REAL K_B;
extern REAL COULOMBIC_CONVERSION_FACTOR;
extern REAL DIELECTRIC_CONSTANT_CONVERSION_FACTOR;
extern REAL ISOTHERMAL_COMPRESSIBILITY_CONVERSION_FACTOR;
extern REAL HEAT_CAPACITY_CONVERSION_FACTOR;
extern REAL VOLUMETRIC_EXPANSION_COEFFICIENT_CONVERSION_FACTOR;
extern REAL FH_CONVERSION_FACTOR;

extern REAL ENERGY_TO_KELVIN;
extern REAL KELVIN_TO_ENERGY;
extern REAL ENERGY_TO_KJ_PER_MOL;
extern REAL ENERGY_TO_EV;
extern REAL ENERGY_TO_KCAL_PER_MOL;
extern REAL ENERGY_TO_KCAL15_PER_MOL;
extern REAL KCAL_PER_MOL_TO_ENERGY;
extern REAL KCAL15_PER_MOL_TO_ENERGY;

extern REAL TO_WAVENUMBERS;
extern REAL TO_THZ;

void SetSimulationUnits(void);

void WriteRestartConstants(FILE *FilePtr);
void ReadRestartConstants(FILE *FilePtr);

#endif
