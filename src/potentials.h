/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'potentials.h' is part of RASPA-2.0

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

#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <stdio.h>
#include "molecule.h"

extern int GeneralMixingRule;
extern int IndividualMixingRules;
extern int IndividualInteractions;

extern int CreateTinkerInput;
extern int CreateDlpolyInput;

enum{NO_MIXING_RULE,LORENTZ_BERTHELOT,JORGENSEN};

enum{NO_BIASING,UMBRELLA,RUIZ_MONTERO};

#define NR_BOND_TYPES 13
enum {HARMONIC_BOND,CORE_SHELL_SPRING,MORSE_BOND,LJ_12_6_BOND,LENNARD_JONES_BOND,BUCKINGHAM_BOND,RESTRAINED_HARMONIC_BOND,
      QUARTIC_BOND,CFF_QUARTIC_BOND,MM3_BOND,RIGID_BOND,FIXED_BOND,MEASURE_BOND};

#define NR_UREYBRADLEY_TYPES 12
enum {HARMONIC_UREYBRADLEY,MORSE_UREYBRADLEY,LJ_12_6_UREYBRADLEY,LENNARD_JONES_UREYBRADLEY,BUCKINGHAM_UREYBRADLEY,RESTRAINED_HARMONIC_UREYBRADLEY,
      QUARTIC_UREYBRADLEY,CFF_QUARTIC_UREYBRADLEY,MM3_UREYBRADLEY,RIGID_UREYBRADLEY,FIXED_UREYBRADLEY,MEASURE_UREYBRADLEY};

#define NR_BEND_TYPES 11
enum {HARMONIC_BEND,CORE_SHELL_BEND,QUARTIC_BEND,CFF_QUARTIC_BEND,HARMONIC_COSINE_BEND,COSINE_BEND,TAFIPOLSKY_BEND,MM3_BEND,MM3_IN_PLANE_BEND,FIXED_BEND,MEASURE_BEND};

#define NR_INVERSION_BEND_TYPES 8
enum {HARMONIC_INVERSION,HARMONIC_COSINE_INVERSION,PLANAR_INVERSION,MM3_INVERSION,HARMONIC_INVERSION2,HARMONIC_COSINE_INVERSION2,PLANAR_INVERSION2,FIXED_INVERSION_BEND};

#define NR_TORSION_TYPES 16
enum {CVFF_DIHEDRAL,HARMONIC_DIHEDRAL,HARMONIC_COSINE_DIHEDRAL,THREE_COSINE_DIHEDRAL,CFF_DIHEDRAL,CFF_DIHEDRAL2,
      FOURIER_SERIES_DIHEDRAL,FOURIER_SERIES_DIHEDRAL2,SIX_COSINE_DIHEDRAL,TRAPPE_DIHEDRAL,TRAPPE_DIHEDRAL_EXTENDED,OPLS_DIHEDRAL,MM3_DIHEDRAL,
      CVFF_BLOCKED_DIHEDRAL,FIXED_DIHEDRAL,MOD_TRAPPE_DIHEDRAL};

#define NR_IMPROPER_TORSION_TYPES 14
enum {CVFF_IMPROPER_DIHEDRAL,HARMONIC_IMPROPER_DIHEDRAL,HARMONIC_COSINE_IMPROPER_DIHEDRAL,THREE_COSINE_IMPROPER_DIHEDRAL,CFF_IMPROPER_DIHEDRAL,
      CFF_IMPROPER_DIHEDRAL2,FOURIER_SERIES_IMPROPER_DIHEDRAL,FOURIER_SERIES_IMPROPER_DIHEDRAL2,SIX_COSINE_IMPROPER_DIHEDRAL,TRAPPE_IMPROPER_DIHEDRAL,
      TRAPPE_IMPROPER_DIHEDRAL_EXTENDED,OPLS_IMPROPER_DIHEDRAL,MM3_IMPROPER_DIHEDRAL,FIXED_IMPROPER_DIHEDRAL};

#define NR_OUT_OF_PLANE_TYPES 1
enum {FIXED_OUT_OF_PLANE_DISTANCE};

#define NR_BOND_BOND_TYPES 3
enum {CVFF_BOND_BOND_CROSS,CFF_BOND_BOND_CROSS,HARMONIC_BOND_BOND_CROSS};

#define NR_BOND_BEND_TYPES 7
enum {CVFF_BOND_BEND_CROSS,CFF_BOND_BEND_CROSS,MM3_BOND_BEND_CROSS,TRUNCATED_HARMONIC,SCREENED_HARMONIC,SCREENED_VESSAL,
      TRUNCATED_VESSAL};

#define NR_BEND_BEND_TYPES 3
enum {CVFF_BEND_BEND_CROSS,CFF_BEND_BEND_CROSS,MM3_BEND_BEND_CROSS};

#define NR_BOND_TORSION_TYPES 1
enum {MM3_BOND_TORSION_CROSS};

#define NR_BEND_TORSION_TYPES 8
enum {SMOOTHED_DIHEDRAL,SMOOTHED_THREE_COSINE_DIHEDRAL,NICHOLAS_DIHEDRAL,SMOOTHED_CFF_DIHEDRAL,SMOOTHED_CFF_DIHEDRAL2,
      CVFF_BEND_TORSION_CROSS,CFF_BEND_TORSION_CROSS,SMOOTHED_CFF_BEND_TORSION_CROSS};

typedef struct potential_Type
{
  int nr_args;
  char Name[50];
} POTENTIAL;

extern POTENTIAL BondTypes[NR_BOND_TYPES];
extern POTENTIAL UreyBradleyTypes[NR_UREYBRADLEY_TYPES];
extern POTENTIAL BendTypes[NR_BEND_TYPES];
extern POTENTIAL InversionBendTypes[NR_INVERSION_BEND_TYPES];
extern POTENTIAL TorsionTypes[NR_TORSION_TYPES];
extern POTENTIAL ImproperTorsionTypes[NR_IMPROPER_TORSION_TYPES];
extern POTENTIAL OutOfPlaneTypes[NR_OUT_OF_PLANE_TYPES];
extern POTENTIAL BondBondTypes[NR_BOND_BOND_TYPES];
extern POTENTIAL BondBendTypes[NR_BOND_BEND_TYPES];
extern POTENTIAL BendBendTypes[NR_BEND_BEND_TYPES];
extern POTENTIAL BondTorsionTypes[NR_BOND_TORSION_TYPES];
extern POTENTIAL BendTorsionTypes[NR_BEND_TORSION_TYPES];

enum {UNDEFINED_POTENTIAL,
      ZERO_POTENTIAL,
      ZERO_POTENTIAL_CONTINUOUS_FRACTIONAL,
      HARD_SPHERE,
      LENNARD_JONES,
      LENNARD_JONES_SMOOTHED3,
      LENNARD_JONES_SMOOTHED5,
      LENNARD_JONES_CONTINUOUS_FRACTIONAL,
      LENNARD_JONES_CONTINUOUS_FRACTIONAL_SMOOTHED3,
      LENNARD_JONES_CONTINUOUS_FRACTIONAL_SMOOTHED5,
      WCA,
      FEYNMAN_HIBBS_LENNARD_JONES,
      FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3,
      FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5,
      FEYNMAN_HIBBS_LENNARD_JONES2,
      FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3,
      FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5,
      LENNARD_JONES_SHIFTED_FORCE,
      LENNARD_JONES_SHIFTED_FORCE2,
      POTENTIAL_12_6,
      POTENTIAL_12_6_SMOOTHED3,
      POTENTIAL_12_6_SMOOTHED5,
      POTENTIAL_12_6_2_0,
      POTENTIAL_12_6_2_0_SMOOTHED3,
      POTENTIAL_12_6_2_0_SMOOTHED5,
      MORSE,
      MORSE_SMOOTHED3,
      MORSE_SMOOTHED5,
      MORSE2,
      MORSE2_SMOOTHED3,
      MORSE2_SMOOTHED5,
      MORSE3,
      MORSE3_SMOOTHED3,
      MORSE3_SMOOTHED5,
      CFF_9_6,
      CFF_9_6_SMOOTHED3,
      CFF_9_6_SMOOTHED5,
      CFF_EPS_SIGMA,
      CFF_EPS_SIGMA_SMOOTHED3,
      CFF_EPS_SIGMA_SMOOTHED5,
      BUCKINGHAM,
      BUCKINGHAM_SMOOTHED3,
      BUCKINGHAM_SMOOTHED5,
      BUCKINGHAM2,
      BUCKINGHAM2_SMOOTHED3,
      BUCKINGHAM2_SMOOTHED5,
      DZUBAK2012,
      MM3_HYDROGEN_VDW,            // obsolete
      MM3_VDW,
      MM3_VDW_SMOOTHED3,
      MM3_VDW_SMOOTHED5,
      MATSUOKA_CLEMENTI_YOSHIMINE,
      MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3,
      MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5,
      MORSE_POTENTIAL,
      MORSE_POTENTIAL_SMOOTHED3,
      MORSE_POTENTIAL_SMOOTHED5,
      GENERIC,
      GENERIC_SMOOTHED3,
      GENERIC_SMOOTHED5,
      PELLENQ_NICHOLSON,
      PELLENQ_NICHOLSON_SMOOTHED3,
      PELLENQ_NICHOLSON_SMOOTHED5,
      HYDRATED_ION_WATER,
      HYDRATED_ION_WATER_SMOOTHED3,
      HYDRATED_ION_WATER_SMOOTHED5,
      MIE,
      MIE_SMOOTHED3,
      MIE_SMOOTHED5,
      MIE_CUTOFF,
      MIE_SMOOTHED3_CUTOFF,
      MIE_SMOOTHED5_CUTOFF,
      BORN_HUGGINS_MEYER,
      BORN_HUGGINS_MEYER_SMOOTHED3,
      BORN_HUGGINS_MEYER_SMOOTHED5,
      HYDROGEN,
      HYDROGEN_SMOOTHED3,
      HYDROGEN_SMOOTHED5};

#define MAX_NUMBER_OF_POTENTIAL_ARGUMENTS 10
extern REAL (**PotentialParms)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS];
extern int **PotentialType;

extern int TailCorrections;
extern int ShiftPotentials;
extern int **TailCorrection;
extern int **ShiftPotential;

extern REAL SwitchingVDWFactors3[4];
extern REAL SwitchingVDWFactors5[6];
extern REAL SwitchingVDWFactors7[8];
extern REAL SwitchingChargeChargeFactors3[4];
extern REAL SwitchingChargeChargeFactors5[6];
extern REAL SwitchingChargeChargeFactors7[8];
extern REAL SwitchingChargeBondDipoleFactors3[4];
extern REAL SwitchingChargeBondDipoleFactors5[6];
extern REAL SwitchingChargeBondDipoleFactors7[8];
extern REAL SwitchingBondDipoleBondDipoleFactors3[4];
extern REAL SwitchingBondDipoleBondDipoleFactors5[6];
extern REAL SwitchingBondDipoleBondDipoleFactors7[8];

extern int NumberOfPseudoAtomsWithoutVDWInteraction;
extern char (*PseudoAtomsWithoutVDWInteraction)[32];

void ChangeVDWtoCFVDW(void);

void ReadForceFieldDefinitionsMixingRules(void);
void ReadForceFieldDefinitions(void);
void ComputePotentialShifts(void);

VECTOR ConvertFromABCtoXYZ(VECTOR t);
VECTOR ConvertFromXYZtoABC(VECTOR t);

VECTOR ConvertFromABCtoXYZUnitCell(VECTOR t);
VECTOR ConvertFromXYZtoABCUnitCell(VECTOR t);

//extern VECTOR ApplyBoundaryCondition(VECTOR dr);

static inline VECTOR ApplyBoundaryCondition(VECTOR dr)
{
  VECTOR s,t;

  switch(BoundaryCondition[CurrentSystem])
  {
    case FINITE:
      break;
    case RECTANGULAR:
    case CUBIC:
      dr.x-=Box[CurrentSystem].ax*(REAL)NINT(dr.x*InverseBox[CurrentSystem].ax);
      dr.y-=Box[CurrentSystem].by*(REAL)NINT(dr.y*InverseBox[CurrentSystem].by);
      dr.z-=Box[CurrentSystem].cz*(REAL)NINT(dr.z*InverseBox[CurrentSystem].cz);
      break;
    case TRICLINIC:
      // convert from xyz to abc
      s.x=InverseBox[CurrentSystem].ax*dr.x+InverseBox[CurrentSystem].bx*dr.y+InverseBox[CurrentSystem].cx*dr.z;
      s.y=InverseBox[CurrentSystem].ay*dr.x+InverseBox[CurrentSystem].by*dr.y+InverseBox[CurrentSystem].cy*dr.z;
      s.z=InverseBox[CurrentSystem].az*dr.x+InverseBox[CurrentSystem].bz*dr.y+InverseBox[CurrentSystem].cz*dr.z;

      // apply boundary condition
      t.x=s.x-(REAL)NINT(s.x);
      t.y=s.y-(REAL)NINT(s.y);
      t.z=s.z-(REAL)NINT(s.z);

      // convert from abc to xyz
      dr.x=Box[CurrentSystem].ax*t.x+Box[CurrentSystem].bx*t.y+Box[CurrentSystem].cx*t.z;
      dr.y=Box[CurrentSystem].ay*t.x+Box[CurrentSystem].by*t.y+Box[CurrentSystem].cy*t.z;
      dr.z=Box[CurrentSystem].az*t.x+Box[CurrentSystem].bz*t.y+Box[CurrentSystem].cz*t.z;
      break;
    default:
      fprintf(stderr,"Error: Unkown boundary condition....\n");
      exit(0);
      break;
  }
  return dr;
}



extern VECTOR ApplyReplicaBoundaryCondition(VECTOR dr);
extern VECTOR ApplyBoundaryConditionUnitCell(VECTOR dr);

extern REAL PotentialValue(int a,int b,REAL rr,REAL scaling);
void PotentialGradient(int a,int b,REAL rr,REAL *energy,REAL *force_factor,REAL scaling);
void PotentialSecondDerivative(int typeA,int typeB,REAL rr,REAL *energy,REAL *factor1,REAL *factor2,REAL scaling);
void PotentialThirdDerivative(int typeA,int typeB,REAL rr,REAL *energy,REAL *factor1,REAL *factor2,REAL *factor3);

REAL PotentialCorrection(int a,int b,REAL r);
REAL PotentialCorrectionPressure(int typeA,int typeB,REAL r);
void CalculateTailCorrection(void);
REAL TailMolecularEnergyDifferenceRXMX(int reaction,int direction);
REAL TailMolecularEnergyDifference(int ComponentToAdd,int ComponentToRemove,int Add,int Remove);
REAL TailMolecularEnergyDifferenceAdd(void);
REAL TailMolecularEnergyDifferenceRemove(void);
REAL TailMolecularEnergyDifferenceAddRemove(int Add,int Remove);

REAL BiasingPotential(int i,VECTOR pos);
REAL BiasingPotentialDerivatives(int i,REAL q,REAL *Derivatives);
REAL BiasingPotentialQ(int i,REAL q);
REAL BiasingPotentialUmbrellaQ(int i,REAL q);
REAL BiasingPotentialUmbrella(int i,VECTOR pos);
REAL BiasingPotentialRuizMontero(int i,VECTOR pos);

void ComputeDummyInteractions(void);

REAL SwitchingFunctionBend(REAL theta);
REAL SwitchingFunctionBendDerivative(REAL theta);

void ComputeDampingCoefficients(REAL r, REAL b,REAL *f6,REAL *f8,REAL *f10);
void ComputeSwitchingFactors(void);

REAL PotentialValueCoulombic(REAL chargeA,REAL chargeB,REAL r);
void PotentialGradientCoulombic(REAL chargeA,REAL chargeB,REAL rr,REAL *energy,REAL *force_factor);
void PotentialSecondDerivativeCoulombic(REAL chargeA,REAL chargeB,REAL rr,REAL *energy,REAL *force_factor1,REAL *force_factor2);

REAL PotentialValueChargeBondDipole(REAL chargeA,VECTOR dipoleB,VECTOR dr,REAL r);
void PotentialGradientChargeBondDipole(REAL chargeB,REAL DipoleMagnitudeA,REAL length,VECTOR dipoleA,VECTOR dr,REAL rr,
                                       REAL *energy,VECTOR *fa1,VECTOR *fa2,VECTOR *fb1,VECTOR *term);

REAL PotentialValueBondDipoleBondDipole(VECTOR dipoleA,VECTOR dipoleB,VECTOR dr,REAL r);
void PotentialGradientBondDipoleBondDipole(REAL DipoleMagnitudeA,REAL ri2,VECTOR dipoleA,REAL DipoleMagnitudeB,REAL rk2,VECTOR dipoleB,VECTOR dr,REAL rr,
                                           REAL *energy,VECTOR *fa1,VECTOR *fa2,VECTOR *fb1,VECTOR *fb2,VECTOR *term);


void PotentialElectricFieldBondDipoleBondDipole(VECTOR dipoleA,VECTOR dipoleB,VECTOR dr,REAL rr,VECTOR *termA,VECTOR *termB);
void PotentialGradientInducedDipoleInducedDipole(VECTOR dipoleA,VECTOR dipoleB,VECTOR dr,REAL rr,REAL *energy,VECTOR *term);

void WriteTinkerParameterFile(void);
void WriteTinkerKeyFile(void);

#endif
