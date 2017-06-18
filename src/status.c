/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'status.c' is part of RASPA-2.0

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
#include <math.h>
#include <string.h>
#include "constants.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "simulation.h"
#include "ewald.h"
#include "potentials.h"
#include "utils.h"
#include "output.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "spacegroup.h"


int PrintFrameworkBondStatus;
int PrintFrameworkUreyBradleyStatus;
int PrintFrameworkBendStatus;
int PrintFrameworkInversionBendStatus;
int PrintFrameworkTorsionStatus;
int PrintFrameworkImproperTorsionStatus;
int PrintFrameworkBondBondStatus;
int PrintFrameworkBondBendStatus;
int PrintFrameworkBendBendStatus;
int PrintFrameworkBendTorsionStatus;
int PrintFrameworkIntraVDWStatus;
int PrintFrameworkIntraChargeChargeStatus;
int PrintFrameworkIntraChargeBondDipoleStatus;
int PrintFrameworkIntraBondDipoleBondDipoleStatus;

int PrintAdsorbateBondStatus;
int PrintAdsorbateUreyBradleyStatus;
int PrintAdsorbateBendStatus;
int PrintAdsorbateTorsionStatus;
int PrintAdsorbateImproperTorsionStatus;
int PrintAdsorbateBondBondStatus;
int PrintAdsorbateBondBendStatus;
int PrintAdsorbateBendBendStatus;
int PrintAdsorbateBondTorsionStatus;
int PrintAdsorbateBendTorsionStatus;
int PrintAdsorbateIntraVDWStatus;
int PrintAdsorbateIntraChargeChargeStatus;
int PrintAdsorbateIntraChargeBondDipoleStatus;
int PrintAdsorbateIntraBondDipoleBondDipoleStatus;

int PrintCationBondStatus;
int PrintCationUreyBradleyStatus;
int PrintCationBendStatus;
int PrintCationTorsionStatus;
int PrintCationImproperTorsionStatus;
int PrintCationBondBondStatus;
int PrintCationBondBendStatus;
int PrintCationBendBendStatus;
int PrintCationBondTorsionStatus;
int PrintCationBendTorsionStatus;
int PrintCationIntraVDWStatus;
int PrintCationIntraChargeChargeStatus;
int PrintCationIntraChargeBondDipoleStatus;
int PrintCationIntraBondDipoleBondDipoleStatus;

int PrintInterVDWStatus;
int PrintInterChargeChargeStatus;
int PrintInterChargeBondDipoleStatus;
int PrintInterBondDipoleBondDipoleStatus;

int PrintFrameworkAdsorbateVDWStatus;
int PrintFrameworkAdsorbateChargeChargeStatus;
int PrintFrameworkAdsorbateChargeBondDipoleStatus;
int PrintFrameworkAdsorbateBondDipoleBondDipoleStatus;

int PrintFrameworkCationVDWStatus;
int PrintFrameworkCationChargeChargeStatus;
int PrintFrameworkCationChargeBondDipoleStatus;
int PrintFrameworkCationBondDipoleBondDipoleStatus;


void SetPrintStatusToFalse(void)
{
  PrintFrameworkBondStatus=FALSE;
  PrintFrameworkUreyBradleyStatus=FALSE;
  PrintFrameworkBendStatus=FALSE;
  PrintFrameworkInversionBendStatus=FALSE;
  PrintFrameworkTorsionStatus=FALSE;
  PrintFrameworkImproperTorsionStatus=FALSE;
  PrintFrameworkBondBondStatus=FALSE;
  PrintFrameworkBondBendStatus=FALSE;
  PrintFrameworkBendBendStatus=FALSE;
  PrintFrameworkBendTorsionStatus=FALSE;
  PrintFrameworkIntraVDWStatus=FALSE;
  PrintFrameworkIntraChargeChargeStatus=FALSE;
  PrintFrameworkIntraChargeBondDipoleStatus=FALSE;
  PrintFrameworkIntraBondDipoleBondDipoleStatus=FALSE;

  PrintAdsorbateBondStatus=FALSE;
  PrintAdsorbateUreyBradleyStatus=FALSE;
  PrintAdsorbateBendStatus=FALSE;
  PrintAdsorbateTorsionStatus=FALSE;
  PrintAdsorbateImproperTorsionStatus=FALSE;
  PrintAdsorbateBondBondStatus=FALSE;
  PrintAdsorbateBondBendStatus=FALSE;
  PrintAdsorbateBendBendStatus=FALSE;
  PrintAdsorbateBondTorsionStatus=FALSE;
  PrintAdsorbateBendTorsionStatus=FALSE;
  PrintAdsorbateIntraVDWStatus=FALSE;
  PrintAdsorbateIntraChargeChargeStatus=FALSE;
  PrintAdsorbateIntraChargeBondDipoleStatus=FALSE;
  PrintAdsorbateIntraBondDipoleBondDipoleStatus=FALSE;

  PrintCationBondStatus=FALSE;
  PrintCationUreyBradleyStatus=FALSE;
  PrintCationBendStatus=FALSE;
  PrintCationTorsionStatus=FALSE;
  PrintCationImproperTorsionStatus=FALSE;
  PrintCationBondBondStatus=FALSE;
  PrintCationBondBendStatus=FALSE;
  PrintCationBendBendStatus=FALSE;
  PrintCationBondTorsionStatus=FALSE;
  PrintCationBendTorsionStatus=FALSE;
  PrintCationIntraVDWStatus=FALSE;
  PrintCationIntraChargeChargeStatus=FALSE;
  PrintCationIntraChargeBondDipoleStatus=FALSE;
  PrintCationIntraBondDipoleBondDipoleStatus=FALSE;

  PrintInterVDWStatus=FALSE;
  PrintInterChargeChargeStatus=FALSE;
  PrintInterChargeBondDipoleStatus=FALSE;
  PrintInterBondDipoleBondDipoleStatus=FALSE;

  PrintFrameworkAdsorbateVDWStatus=FALSE;
  PrintFrameworkAdsorbateChargeChargeStatus=FALSE;
  PrintFrameworkAdsorbateChargeBondDipoleStatus=FALSE;
  PrintFrameworkAdsorbateBondDipoleBondDipoleStatus=FALSE;

  PrintFrameworkCationVDWStatus=FALSE;
  PrintFrameworkCationChargeChargeStatus=FALSE;
  PrintFrameworkCationChargeBondDipoleStatus=FALSE;
  PrintFrameworkCationBondDipoleBondDipoleStatus=FALSE;
}


// =============================================================================================================================
// auxilary routines for printing
// =============================================================================================================================

REAL PrintBondEnergyStatus(int nr,char *string,int BondType,REAL *parms,REAL r)
{
  REAL rr,r1,U;
  REAL temp,temp2,exp_term;

  rr=SQR(r);

  switch(BondType)
  {
    case HARMONIC_BOND:
      // 0.5*p0*SQR(r-p1);
      // ===============================================
      // p_0/k_B [K/A^2]   force constant
      // p_1     [A]       reference bond distance
      U=0.5*parms[0]*SQR(r-parms[1]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_BOND %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CORE_SHELL_SPRING:
      U=0.5*parms[0]*SQR(r);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CORE_SHELL_SPRING %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MORSE_BOND:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference bond distance
      temp=exp(parms[1]*(parms[2]-r));
      U=parms[0]*(SQR(1.0-temp)-1.0);
      fprintf(OutputFilePtr[CurrentSystem],"%4d MORSE_BOND %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1]  p_2=%8.5f [A], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case LJ_12_6_BOND:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      temp=CUBE(1.0/rr);
      U=parms[0]*SQR(temp)-parms[1]*temp;
      fprintf(OutputFilePtr[CurrentSystem],"%4d LJ_12_6_BOND %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [K A^6], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case LENNARD_JONES_BOND:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A]
      temp=CUBE(parms[1]/rr);
      U=4.0*parms[0]*(temp*(temp-1.0));
      fprintf(OutputFilePtr[CurrentSystem],"%4d LENNARD_JONES_BOND %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case BUCKINGHAM_BOND:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      temp=parms[2]*CUBE(1.0/rr);
      exp_term=parms[0]*exp(-parms[1]*r);
      U=-temp+exp_term;
      fprintf(OutputFilePtr[CurrentSystem],"%4d BUCKINGHAM_BOND %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^6], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*ENERGY_TO_KELVIN,
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case RESTRAINED_HARMONIC_BOND:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      r1=r-parms[1];
      U=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
            +parms[0]*parms[2]*MAX2(fabs(r1)-parms[2],(REAL)0.0);
      fprintf(OutputFilePtr[CurrentSystem],"%4d RESTRAINED_HARMONIC_BOND %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], p_2=%8.5f [K], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case QUARTIC_BOND:
      // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
      // ===========================================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2/k_B [K/A^3]
      // p_3/k_B [K/A^4]
      temp=r-parms[1];
      temp2=SQR(r-parms[1]);
      U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d QUARTIC_BOND %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], p_2/k_B=%8.5f [K/A^3], p_3/k_B=%8.5f [K/A^4], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_QUARTIC_BOND:
      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2/k_B [K/A^3]
      // p_3/k_B [K/A^4]
      temp=r-parms[1];
      temp2=SQR(r-parms[1]);
      U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_QUARTIC_BOND %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], p_2/k_B=%8.5f [K/A^3], p_3/k_B=%8.5f [K/A^4], Distance %10.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_BOND:
      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
      // =================================================================
      // p_0     [mdyne/A molecule]
      // p_1     [A]
      temp=r-parms[1];
      temp2=SQR(temp);
      U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_BOND %s, p_0=%8.5f [mdyne/A molecule], p_1=%8.5f [A], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]/(71.94*KCAL_PER_MOL_TO_ENERGY),
                parms[1],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MEASURE_BOND:
      U=0.0;
      break;
    case RIGID_BOND:
      U=0.0;
      break;
    case FIXED_BOND:
      U=0.0;
      break;
    default:
      fprintf(stderr, "Undefined Bond potential in routine 'PrintBondEnergyStatus' ('status.c')\n");
      exit(0);
      break;
  }
  return U;
}

REAL PrintUreyBradleyEnergyStatus(int nr,char *string,int UreyBradleyType,REAL *parms,REAL r)
{
  REAL rr,r1,U;
  REAL temp,temp2,exp_term;

  rr=SQR(r);

  switch(UreyBradleyType)
  {
    case HARMONIC_UREYBRADLEY:
      // 0.5*p0*SQR(r-p1);
      // ===============================================
      // p_0/k_B [K/A^2]   force constant
      // p_1     [A]       reference bond distance
      U=0.5*parms[0]*SQR(r-parms[1]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_UREYBRADLEY %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MORSE_UREYBRADLEY:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference bond distance
      temp=exp(parms[1]*(parms[2]-r));
      U=parms[0]*(SQR(1.0-temp)-1.0);
      fprintf(OutputFilePtr[CurrentSystem],"%4d MORSE_UREYBRADLEY %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1]  p_2=%8.5f [A], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case LJ_12_6_UREYBRADLEY:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      temp=CUBE(1.0/rr);
      U=parms[0]*SQR(temp)-parms[1]*temp;
      fprintf(OutputFilePtr[CurrentSystem],"%4d LJ_12_6_UREYBRADLEY %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [K A^6], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case LENNARD_JONES_UREYBRADLEY:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A]
      temp=CUBE(parms[1]/rr);
      U=4.0*parms[0]*(temp*(temp-1.0));
      fprintf(OutputFilePtr[CurrentSystem],"%4d LENNARD_JONES_UREYBRADLEY %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case BUCKINGHAM_UREYBRADLEY:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      temp=parms[2]*CUBE(1.0/rr);
      exp_term=parms[0]*exp(-parms[1]*r);
      U=-temp+exp_term;
      fprintf(OutputFilePtr[CurrentSystem],"%4d BUCKINGHAM_UREYBRADLEY %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^6], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*ENERGY_TO_KELVIN,
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case RESTRAINED_HARMONIC_UREYBRADLEY:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      r1=r-parms[1];
      U=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
            +parms[0]*parms[2]*MAX2(fabs(r1)-parms[2],(REAL)0.0);
      fprintf(OutputFilePtr[CurrentSystem],"%4d RESTRAINED_HARMONIC_UREYBRADLEY %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], p_2=%8.5f [K], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case QUARTIC_UREYBRADLEY:
      // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
      // ===========================================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2/k_B [K/A^3]
      // p_3/k_B [K/A^4]
      temp=r-parms[1];
      temp2=SQR(r-parms[1]);
      U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d QUARTIC_UREYBRADLEY %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], p_2/k_B=%8.5f [K/A^3], p_3/k_B=%8.5f [K/A^4], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_QUARTIC_UREYBRADLEY:
      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2/k_B [K/A^3]
      // p_3/k_B [K/A^4]
      temp=r-parms[1];
      temp2=SQR(r-parms[1]);
      U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_QUARTIC_UREYBRADLEY %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], p_2/k_B=%8.5f [K/A^3], p_3/k_B=%8.5f [K/A^4], Distance %10.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_UREYBRADLEY:
      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
      // =================================================================
      // p_0     [mdyne/A molecule]
      // p_1     [A]
      temp=r-parms[1];
      temp2=SQR(temp);
      U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_UREYBRADLEY %s, p_0=%8.5f [mdyne/A molecule], p_1=%8.5f [A], Distance %8.5f [K], Energy: %12.5f [K] %9.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]/(71.94*KCAL_PER_MOL_TO_ENERGY),
                parms[1],
                r,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MEASURE_UREYBRADLEY:
      U=0.0;
      break;
    case RIGID_UREYBRADLEY:
      U=0.0;
      break;
    case FIXED_UREYBRADLEY:
      U=0.0;
      break;
    default:
      fprintf(stderr, "Undefined UreyBradley potential in routine 'PrintUreyBradleyEnergyStatus' ('status.c')\n");
      exit(0);
      break;
  }
  return U;
}

REAL PrintBendEnergyStatus(int nr,char *string,int BendType,REAL *parms,REAL Theta)
{
  REAL U,temp,temp2;
  REAL CosTheta;

  CosTheta=cos(Theta);
  switch(BendType)
  {
    case HARMONIC_BEND:
      // (1/2)p_0*(Theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      U=0.5*parms[0]*SQR(Theta-parms[1]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_BEND %s, p_0/k_B=%8.5f [K/rad^2], p_1=%8.5f [degrees], Theta: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                Theta*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CORE_SHELL_BEND:
      // (1/2)p_0*(Theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      U=0.5*parms[0]*SQR(Theta-parms[1]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CORE_SHELL_BEND %s, p_0/k_B=%8.5f [K/rad^2], p_1=%8.5f [degrees], Theta: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                Theta*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case QUARTIC_BEND:
      // (1/2)p_0*(Theta-p_1)^2+(1/3)*p_2*(Theta-p_1)^3+(1/4)*p_2*(Theta-p_1)^4
      // ======================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      temp=(Theta-parms[1]);
      temp2=SQR(temp);
      U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_QUARTIC_BEND %s, p_0/k_B=%8.5f [K/rad^2], p_1=%8.5f [degrees], p_2/k_B=%8.5f [K/rad^3], p_3/k_B=%8.5f [K/rad^4], Theta: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                Theta*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_QUARTIC_BEND:
      // p_0*(Theta-p_1)^2+p_2*(Theta-p_1)^3+p_3*(Theta-p_1)^4
      // =====================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      temp=(Theta-parms[1]);
      temp2=SQR(temp);
      U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_QUARTIC_BEND %s, p_0/k_B=%8.5f [K/rad^2], p_1=%8.5f [degrees], p_2/k_B=%8.5f [K/rad^3], p_3/k_B=%8.5f [K/rad^4], Theta: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                Theta*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HARMONIC_COSINE_BEND:
      // (1/2)*p_0*(cos(Theta)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      U=0.5*parms[0]*SQR(CosTheta-parms[1]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_COSINE_BEND %s, p_0/k_B=%8.5f [K], p_1=%8.5f [degrees], Theta: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                acos(parms[1])*RAD2DEG,
                Theta*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case COSINE_BEND:
      // p_0*(1+cos(p_1*Theta-p_2))
      // ===============================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      temp=parms[1]*Theta-parms[2];
      U=parms[0]*(1.0+cos(temp));
      fprintf(OutputFilePtr[CurrentSystem],"%4d COSINE_BEND %s, p_0/k_B=%8.5f [K], p_1=%8.5f [-], p_2=%8.5f [degrees], Theta: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*RAD2DEG,
                Theta*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case TAFIPOLSKY_BEND:
      // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
      // ===============================================
      // p_0/k_B [K]
      U=0.5*parms[0]*(1+cos(Theta))*(1+2.0*cos(Theta));
      fprintf(OutputFilePtr[CurrentSystem],"%4d TAFIPOLSKY_BEND %s, p_0/k_B=%8.5f [K/rad^2], Theta: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                Theta*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_BEND:
      // p_0*(Theta-p_1)^2(1-0.014*(Theta-p_1)+5.6e-5*(Theta-p_1)^2-7e-7*(Theta-p_1)^3+2.2e-8(Theta-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      temp=RAD2DEG*(Theta-parms[1]);
      temp2=SQR(temp);
      U=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_BEND %s, p_0=%8.5f [mdyne A/rad^2], p_1=%8.5f [degrees], Theta: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]/(0.02191418*KCAL_PER_MOL_TO_ENERGY),
                parms[1]*RAD2DEG,
                Theta*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_IN_PLANE_BEND:
      // p_0*(Theta-p_1)^2(1-0.014*(Theta-p_1)+5.6e-5*(Theta-p_1)^2-7e-7*(Theta-p_1)^3+2.2e-8(Theta-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      temp=RAD2DEG*(Theta-parms[1]);
      temp2=SQR(temp);
      U=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_IN_PLANE_BEND %s, p_0=%8.5f [mdyne A/rad^2], p_1=%8.5f [degrees], Theta: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]/(0.02191418*KCAL_PER_MOL_TO_ENERGY),
                parms[1]*RAD2DEG,
                Theta*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case FIXED_BEND:
      U=0.0;
      break;
    case MEASURE_BEND:
      U=0.0;
      break;
    default:
      fprintf(stderr, "Undefined Bend potential in routine 'CalculateFrameworkBendEnergy' ('framework_energy.c')\n");
      exit(0);
      break;
  }
  return U;
}

REAL PrintInversionBendEnergyStatus(int nr,char *string,int InversionBendType,REAL *parms,REAL Chi)
{
  REAL energy,CosChi;
  REAL temp,temp2;

  CosChi=cos(Chi);
  switch(InversionBendType)
  {
    case HARMONIC_INVERSION:
      // (1/2)*p_0*(chi-p_1)^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      energy=0.5*parms[0]*SQR(Chi-parms[1]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_INVERSION %s,, p_0/k_B=%8.5f [K], p_1=%8.5f [degrees], Chi (Wilson): %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
              nr,
              string,
              parms[0]*ENERGY_TO_KELVIN,
              parms[1]*RAD2DEG,
              Chi*RAD2DEG,
              energy*ENERGY_TO_KELVIN,
              energy*ENERGY_TO_KJ_PER_MOL,
              energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HARMONIC_INVERSION2:
      // (1/2)*p_0*(chi-p_1)^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      energy=0.5*parms[0]*SQR(Chi-parms[1]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_INVERSION2 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [degrees], Chi (Allinger): %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
              nr,
              string,
              parms[0]*ENERGY_TO_KELVIN,
              parms[1]*RAD2DEG,
              Chi*RAD2DEG,
              energy*ENERGY_TO_KELVIN,
              energy*ENERGY_TO_KJ_PER_MOL,
              energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HARMONIC_COSINE_INVERSION:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      energy=0.5*parms[0]*SQR(CosChi-parms[1]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_COSINE_INVERSION %s, p_0/k_B=%8.5f [K], p_1=%8.5f [degrees], Theta (Wilson): %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
              nr,
              string,
              parms[0]*ENERGY_TO_KELVIN,
              acos(parms[1])*RAD2DEG,
              Chi*RAD2DEG,
              energy*ENERGY_TO_KELVIN,
              energy*ENERGY_TO_KJ_PER_MOL,
              energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HARMONIC_COSINE_INVERSION2:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      energy=0.5*parms[0]*SQR(CosChi-parms[1]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_COSINE_INVERSION2 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [degrees], Theta (Allinger): %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
              nr,
              string,
              parms[0]*ENERGY_TO_KELVIN,
              acos(parms[1])*RAD2DEG,
              Chi*RAD2DEG,
              energy*ENERGY_TO_KELVIN,
              energy*ENERGY_TO_KJ_PER_MOL,
              energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case PLANAR_INVERSION:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      energy=parms[0]*(1.0-CosChi);
      fprintf(OutputFilePtr[CurrentSystem],"%4d PLANAR_INVERSION %s, p_0/k_B=%8.5f [K], Theta (Wilson): %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
              nr,
              string,
              parms[0]*ENERGY_TO_KELVIN,
              Chi*RAD2DEG,
              energy*ENERGY_TO_KELVIN,
              energy*ENERGY_TO_KJ_PER_MOL,
              energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case PLANAR_INVERSION2:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      energy=parms[0]*(1.0-CosChi);
      fprintf(OutputFilePtr[CurrentSystem],"%4d PLANAR_INVERSION2 %s, p_0/k_B=%8.5f [K], Chi (Allinger): %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
              nr,
              string,
              parms[0]*ENERGY_TO_KELVIN,
              Chi*RAD2DEG,
              energy*ENERGY_TO_KELVIN,
              energy*ENERGY_TO_KJ_PER_MOL,
              energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_INVERSION:
      // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      temp=RAD2DEG*(Chi-parms[1]);
      temp2=SQR(temp);
      energy=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_INVERSION %s, p_0/k_B=%8.5f [mdyne A/rad^2], p_1=%8.5f [degrees], Chi (Allinger): %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
              nr,
              string,
              parms[0]/(0.02191418*KCAL_PER_MOL_TO_ENERGY),
              parms[1]*RAD2DEG,
              Chi*RAD2DEG,
              energy*ENERGY_TO_KELVIN,
              energy*ENERGY_TO_KJ_PER_MOL,
              energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
      exit(0);
      break;
  }
 return energy;
}

REAL PrintTorsionEnergyStatus(int nr,char *string,int TorsionType,REAL *parms,REAL Phi)
{
  REAL U,CosPhi,CosPhi2;

  CosPhi=cos(Phi);
  CosPhi2=SQR(CosPhi);
  switch(TorsionType)
  {
    case HARMONIC_DIHEDRAL:
      // (1/2)*p_0*(phi-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // potential defined in terms of 'phi' and therefore contains a singularity
      // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
      // same direction as Rbc, and negative otherwise
      Phi-=parms[1];
      Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
      U=0.5*parms[0]*SQR(Phi);
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_DIHEDRAL %s, p_0/k_B=%8.5f [K/rad^2],  p_1=%8.5f [A], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HARMONIC_COSINE_DIHEDRAL:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      U=0.5*parms[0]*SQR(CosPhi-cos(parms[1]));
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_COSINE_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1=%8.5f [degrees], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                acos(parms[1])*RAD2DEG,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case THREE_COSINE_DIHEDRAL:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d THREE_COSINE_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_DIHEDRAL:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0     [kcal/mol]
      // p_1     [kcal/mol]
      // p_2     [kcal/mol]
      U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_DIHEDRAL %s, p_0=%8.5f [kcal/mol], p_1=%8.5f [kcal/mol], p_2=%8.5f [kcal/mol], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]/KCAL_PER_MOL_TO_ENERGY,
                parms[1]/KCAL_PER_MOL_TO_ENERGY,
                parms[2]/KCAL_PER_MOL_TO_ENERGY,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CVFF_BLOCKED_DIHEDRAL:
      //
      // ========================================================================
      // p_0     [rad]
      // p_1     [K]
      // p_2     [-]
      // p_3     [rad]
      // p_4     [rad]
      break;
    case CFF_DIHEDRAL:
      // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_DIHEDRAL2:
      // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_DIHEDRAL2 %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SIX_COSINE_DIHEDRAL:
      // Prod_i=0^5 p_i*cos(phi)^i
      // =========================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
      // between the first and last atoms of the dihedral, and Phi'=Phi-pi is defined accoording to the
      // polymer convention Phi'(trans)=0.
      U=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
             parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
      fprintf(OutputFilePtr[CurrentSystem],"%4d SIX_COSINE_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "p_3/k_B=%8.5f [K], p_4/k_B=%8.5f [K], p_5/k_B=%8.5f [K], Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                parms[4]*ENERGY_TO_KELVIN,
                parms[5]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case TRAPPE_DIHEDRAL:
      // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
      // =============================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
      fprintf(OutputFilePtr[CurrentSystem],"%4d TRAPPE_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "p_3/k_B=%8.5f [K], Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CVFF_DIHEDRAL:
      // p_0*(1+cos(p_1*phi-p_2))
      // ========================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
      fprintf(OutputFilePtr[CurrentSystem],"%4d CVFF_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1=%8.5f [-], p_2=%8.5f [degrees], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*RAD2DEG,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case OPLS_DIHEDRAL:
      // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
      // =================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
      fprintf(OutputFilePtr[CurrentSystem],"%4d OPLS_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "p_3/k_B=%8.5f [K], Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case FOURIER_SERIES_DIHEDRAL:
      // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
      // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
      // =======================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      U=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
             2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi2+4.0*(parms[2]-5.0*parms[4])*CosPhi2*CosPhi-
             8.0*(parms[3]-6.0*parms[5])*SQR(CosPhi2)+16.0*parms[4]*SQR(CosPhi2)*CosPhi-32.0*parms[5]*CUBE(CosPhi2));
      fprintf(OutputFilePtr[CurrentSystem],"%4d FOURIER_SERIES_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "p_3/k_B=%8.5f [K], p_4/k_B=%8.5f [K], p_5/k_B=%8.5f [K], Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                parms[4]*ENERGY_TO_KELVIN,
                parms[5]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case FOURIER_SERIES_DIHEDRAL2:
      // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
      // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
      // =======================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      U=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
        2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
        2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
      fprintf(OutputFilePtr[CurrentSystem],"%4d FOURIER_SERIES_DIHEDRAL2 %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "p_3/k_B=%8.5f [K], p_4/k_B=%8.5f [K], p_5/k_B=%8.5f [K], Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                parms[4]*ENERGY_TO_KELVIN,
                parms[5]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Undefined Torsion potential in routine 'CalculateFrameworkTorsionEnergy' ('framework_energy.c')\n");
      exit(0);
      break;
  }
  return U;
}

REAL PrintImproperTorsionEnergyStatus(int nr,char *string,int TorsionType,REAL *parms,REAL Phi)
{
  REAL U,CosPhi,CosPhi2;

  CosPhi=cos(Phi);
  CosPhi2=SQR(CosPhi);
  switch(TorsionType)
  {
    case HARMONIC_IMPROPER_DIHEDRAL:
      // (1/2)*p_0*(phi-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // potential defined in terms of 'phi' and therefore contains a singularity
      // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
      // same direction as Rbc, and negative otherwise
      Phi-=parms[1];
      Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
      U=0.5*parms[0]*SQR(Phi);
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_IMPROPER_DIHEDRAL %s, p_0/k_B=%8.5f [K/rad^2],  p_1=%8.5f [A], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      U=0.5*parms[0]*SQR(CosPhi-cos(parms[1]));
      fprintf(OutputFilePtr[CurrentSystem],"%4d HARMONIC_COSINE_IMPROPER_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1=%8.5f [degrees], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                acos(parms[1])*RAD2DEG,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case THREE_COSINE_IMPROPER_DIHEDRAL:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d THREE_COSINE_IMPROPER_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_IMPROPER_DIHEDRAL:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0     [kcal/mol]
      // p_1     [kcal/mol]
      // p_2     [kcal/mol]
      U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_IMPROPER_DIHEDRAL %s, p_0=%8.5f [kcal/mol], p_1=%8.5f [kcal/mol], p_2=%8.5f [kcal/mol], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]/KCAL_PER_MOL_TO_ENERGY,
                parms[1]/KCAL_PER_MOL_TO_ENERGY,
                parms[2]/KCAL_PER_MOL_TO_ENERGY,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_IMPROPER_DIHEDRAL:
      // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_IMPROPER_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_IMPROPER_DIHEDRAL2:
      // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_IMPROPER_DIHEDRAL2 %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SIX_COSINE_IMPROPER_DIHEDRAL:
      // Prod_i=0^5 p_i*cos(phi)^i
      // =========================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
      // between the first and last atoms of the dihedral, and Phi'=Phi-pi is defined accoording to the
      // polymer convention Phi'(trans)=0.
      U=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
             parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
      fprintf(OutputFilePtr[CurrentSystem],"%4d SIX_COSINE_IMPROPER_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "p_3/k_B=%8.5f [K], p_4/k_B=%8.5f [K], p_5/k_B=%8.5f [K], Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                parms[4]*ENERGY_TO_KELVIN,
                parms[5]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case TRAPPE_IMPROPER_DIHEDRAL:
      // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
      // =============================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
      fprintf(OutputFilePtr[CurrentSystem],"%4d TRAPPE_IMPROPER_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "p_3/k_B=%8.5f [K], Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CVFF_IMPROPER_DIHEDRAL:
      // p_0*(1+cos(p_1*phi-p_2))
      // ========================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
      fprintf(OutputFilePtr[CurrentSystem],"%4d CVFF_IMPROPER_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1=%8.5f [-], p_2=%8.5f [degrees], "
                                           "Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*RAD2DEG,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case OPLS_IMPROPER_DIHEDRAL:
      // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
      // =================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
      fprintf(OutputFilePtr[CurrentSystem],"%4d OPLS_IMPROPER_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "p_3/k_B=%8.5f [K], Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case FOURIER_SERIES_IMPROPER_DIHEDRAL:
      // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
      // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
      // =======================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      U=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
             2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi2+4.0*(parms[2]-5.0*parms[4])*CosPhi2*CosPhi-
             8.0*(parms[3]-6.0*parms[5])*SQR(CosPhi2)+16.0*parms[4]*SQR(CosPhi2)*CosPhi-32.0*parms[5]*CUBE(CosPhi2));
      fprintf(OutputFilePtr[CurrentSystem],"%4d FOURIER_SERIES_IMPROPER_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "p_3/k_B=%8.5f [K], p_4/k_B=%8.5f [K], p_5/k_B=%8.5f [K], Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                parms[4]*ENERGY_TO_KELVIN,
                parms[5]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case FOURIER_SERIES_IMPROPER_DIHEDRAL2:
      // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
      // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
      // =======================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      U=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
        2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
        2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
      fprintf(OutputFilePtr[CurrentSystem],"%4d FOURIER_SERIES_IMPROPER_DIHEDRAL2 %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "p_3/k_B=%8.5f [K], p_4/k_B=%8.5f [K], p_5/k_B=%8.5f [K], Phi: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                parms[3]*ENERGY_TO_KELVIN,
                parms[4]*ENERGY_TO_KELVIN,
                parms[5]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Undefined Torsion potential in routine 'CalculateFrameworkTorsionEnergy' ('framework_energy.c')\n");
      exit(0);
      break;
  }
  return U;
}

REAL PrintBondBondEnergyStatus(int nr,char *string, int BondBondType,REAL *parms,REAL rab,REAL rbc)
{
  REAL energy;

  switch(BondBondType)
  {
    case CVFF_BOND_BOND_CROSS:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      energy=parms[0]*(rab-parms[1])*(rbc-parms[2]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CVFF_BOND_BOND_CROSS %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], p_2=%8.5f [A], "
                                           "Distance %8.5f [A], Distance %8.5f [A], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr++,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2],
                rab,
                rbc,
                energy*ENERGY_TO_KELVIN,
                energy*ENERGY_TO_KJ_PER_MOL,
                energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_BOND_BOND_CROSS:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      energy=parms[0]*(rab-parms[1])*(rbc-parms[2]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_BOND_BOND_CROSS %s, p_0/k_B=%8.5f [K/A^2], p_1=%8.5f [A], p_2=%8.5f [A], "
                                           "Distance %8.5f [A], Distance %8.5f [A], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr++,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2],
                rab,
                rbc,
                energy*ENERGY_TO_KELVIN,
                energy*ENERGY_TO_KJ_PER_MOL,
                energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Undefined Bond-Bond potential in routine 'CalculateFrameworkBondBondEnergy' ('status.c')\n");
      exit(0);
      break;
  }
  return energy;
}

REAL PrintBondBendEnergyStatus(int nr,char *string,int BondBendType,REAL *parms,REAL rab,REAL rbc,REAL Theta)
{
  REAL U;

  switch(BondBendType)
  {
    case CVFF_BOND_BEND_CROSS:
      // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
      // =========================================
      // p_0     [degrees]
      // p_1/k_B [K/A/rad]
      // p_2     [A]
      // p_3/k_B [K/A/rad]
      // p_4     [A]
      U=(Theta-parms[0])*(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]));
      fprintf(OutputFilePtr[CurrentSystem],"%4d CVFF_BOND_BEND_CROSS %s, p_0=%8.5f [degrees], p_1/k_B=%8.5f [K/A/rad], p_2=%8.5f [A], "
                                           "p_3/k_B=%8.5f [K/A/rad], p_4=%8.5f [A], Theta: %8.5f [degrees], Distance %8.5f [A], Distance %8.5f [A], "
                                           "Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*RAD2DEG,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2],
                parms[3]*ENERGY_TO_KELVIN,
                parms[4],
                Theta,
                rab,
                rbc,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_BOND_BEND_CROSS:
      // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
      // =========================================
      // p_0     [degrees]
      // p_1/k_B [K/A/rad]
      // p_2     [A]
      // p_3/k_B [K/A/rad]
      // p_4     [A]
      U=(Theta-parms[0])*(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]));
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_BOND_BEND_CROSS %s, p_0=%8.5f [degrees], p_1/k_B=%8.5f [K/A/rad], p_2=%8.5f [A], "
                                           "p_3/k_B=%8.5f [K/A/rad], p_4=%8.5f [A], Theta: %8.5f [degrees], Distance %8.5f [A], Distance %8.5f [A], "
                                           "Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*RAD2DEG,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2],
                parms[3]*ENERGY_TO_KELVIN,
                parms[4],
                Theta,
                rab,
                rbc,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_BOND_BEND_CROSS:
      // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
      // =====================================
      // p_0     [mdyne/rad]
      // p_1     [A]
      // p_2     [A]
      // p_3     [degrees]
      U=parms[0]*((rab-parms[1])+(rbc-parms[2]))*RAD2DEG*(Theta-parms[3]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_BOND_BEND_CROSS %s, p_0/k_B=%8.5f [mdyne/rad], p_1=%8.5f [A], p_2=%8.5f [A], "
                                           "p_3=%8.5f [degrees], Theta: %8.5f [degrees], Distance %8.5f [A], Distance %8.5f [A], "
                                           "Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]/(2.51118*KCAL_PER_MOL_TO_ENERGY),
                parms[1],
                parms[2],
                parms[3]*RAD2DEG,
                Theta,
                rab,
                rbc,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case TRUNCATED_HARMONIC:
      // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
      // ================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      U=0.5*parms[0]*SQR(Theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8));
      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_HARMONIC %s, p_0/k_B=%8.5f [K/rad^2], p_1=%8.5f [degrees], p_2=%8.5f [A], "
                                           "Theta: %8.5f [degrees], Distance %8.5f [A], Distance %8.5f [A], "
                                           "Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2],
                Theta,
                rab,
                rbc,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SCREENED_HARMONIC:
     // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
     // ===============================================
     // p_0/k_B [K/rad^2]
     // p_1     [degrees]
     // p_2     [A]
     // p_3     [A]
     U=0.5*parms[0]*SQR(Theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3]));
     fprintf(OutputFilePtr[CurrentSystem],"%4d SCREENED_HARMONIC %s, p_0/k_B=%8.5f [K/rad^2], p_1=%8.5f [degrees], p_2=%8.5f [A], "
                                           "p_3=%8.5f [A], Theta: %8.5f [degrees], Distance %8.5f [A], Distance %8.5f [A], "
                                           "Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2],
                parms[3],
                Theta,
                rab,
                rbc,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SCREENED_VESSAL:
      // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
      // ============================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      // p_3     [A]
      U=(parms[0]/(8.0*SQR(Theta-M_PI)))*SQR(SQR(parms[1]-M_PI)-SQR(Theta-M_PI))
            *exp(-(rab/parms[2]+rbc/parms[3]));
      fprintf(OutputFilePtr[CurrentSystem],"%4d SCREENED_VESSAL %s, p_0/k_B=%8.5f [K/rad^2], p_1=%8.5f [degrees], p_2=%8.5f [A], "
                                           "p_3=%8.5f [A], Theta: %8.5f [degrees], Distance %8.5f [A], Distance %8.5f [A], "
                                           "Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2],
                parms[3],
                Theta,
                rab,
                rbc,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case TRUNCATED_VESSAL:
      // p_0*[pow(Theta,p_2)*(Theta-p_1)^2*(Theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(Theta-p_1)^2*pow(PI-p_1,3)]
      //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
      // ============================================================================
      // p_0/k_B [K/rad^(4+p_2)]
      // p_1     [degrees]
      // p_2     [-]
      // p_3     [A]
      U=parms[0]*(pow(Theta,parms[2])*SQR(Theta-parms[1])*SQR(Theta+parms[1]-2.0*M_PI)
            -0.5*parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*SQR(Theta-parms[1])*pow(M_PI-parms[1],(REAL)3.0))
            *exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8));
      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_VESSAL %s, p_0/k_B=%8.5f [K/rad^(4+p_2)], p_1=%8.5f [degrees], p_2=%8.5f [-], "
                                           "p_3=%8.5f [A], Theta: %8.5f [degrees], Distance %8.5f [A], Distance %8.5f [A], "
                                           "Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2],
                parms[3],
                Theta,
                rab,
                rbc,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Undefined Bond-Bend potential in routine 'CalculateFrameworkBondBendEnergy' ('status.c')\n");
      exit(0);
      break;
  }
  return U;
}

REAL PrintBendBendEnergyStatus(int nr,char *string,int BendBendType,REAL *parms,REAL Theta1,REAL Theta2)
{
  REAL U;

  switch(BendBendType)
  {
    case CVFF_BEND_BEND_CROSS:
      // p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0/k_B [K/rad^2)]
      // p_1     [degrees]
      // p_2     [degrees]
      U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CVFF_BEND_BEND_CROSS %s, p_0/k_B=%8.5f [K], p_1=%8.5f [degrees], p_2=%8.5f [degrees], "
                                           "Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2]*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_BEND_BEND_CROSS:
      // p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0/k_B [K/rad^2)]
      // p_1     [degrees]
      // p_2     [degrees]
      U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_BEND_BEND_CROSS %s, p_0/k_B=%8.5f [K/rad^2], p_1=%8.5f [degrees], p_2=%8.5f [degrees], "
                                           "Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2]*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_BEND_BEND_CROSS:
      // -p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0     [mdyne A/rad^2]
      // p_1     [degrees]
      // p_2     [degrees]
      U=-parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])*(Theta2-parms[2]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_BEND_BEND_CROSS %s, p_0=%8.5f [mdyne A/rad^2], p_1=%8.5f [degrees], "
                                           "p_2=%8.5f [degrees], Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]/(0.02191418*KCAL_PER_MOL_TO_ENERGY),
                parms[1]*RAD2DEG,
                parms[2]*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Undefined Bend-Bend potential in routine 'CalculateFrameworkBendBendEnergy' ('framework_energy.c')\n");
      exit(0);
      break;
  }
  return U;
}

REAL PrintBondTorsionEnergyStatus(int nr,char *string,int BondTorsionType,REAL *parms,REAL rbc,REAL Phi)
{
  REAL U,temp;
  REAL CosPhi,CosPhi2;

  CosPhi=cos(Phi);
  CosPhi2=SQR(CosPhi);

  switch(BondTorsionType)
  {
    case MM3_BOND_TORSION_CROSS:
      // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
      // =====================================================================================
      // p_0     [kcal/A mole]
      // p_1     [kcal/A mole]
      // p_2     [kcal/A mole]
      // p_3     [A]
      temp=(rbc-parms[3]);
      U=parms[0]*temp*CosPhi+parms[1]*temp*(2.0*CosPhi2-1.0)+parms[2]*temp*(4.0*CosPhi2*CosPhi-3.0*CosPhi);
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_BOND_TORSION_CROSS %s, p_0=%8.5f [kcal/A mole], p_1=%8.5f [kcal/A mole], p_2=%8.5f [kcal/A mole], "
                                           "p_3=%8.5f [kcal/A mole], Phi: %8.5f [degrees], Distance %8.5f [A], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]/(KCAL_PER_MOL_TO_ENERGY),
                parms[1]/(KCAL_PER_MOL_TO_ENERGY),
                parms[2]/(KCAL_PER_MOL_TO_ENERGY),
                parms[3],
                Phi*RAD2DEG,
                rbc,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Undefined Bond-Torsion potential in routine 'CalculateFrameworkBondTorsionEnergy' ('status.c')\n");
      exit(0);
      break;
  }
  return U;
}

REAL PrintBendTorsionEnergyStatus(int nr,char *string,int BendTorsionType,REAL *parms,REAL Theta1,REAL Theta2,REAL Phi)
{
  REAL U;
  REAL CosPhi,CosPhi2;

  CosPhi=cos(Phi);
  CosPhi2=SQR(CosPhi);

  switch(BendTorsionType)
  {
    case CVFF_BEND_TORSION_CROSS:
      // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
      // =====================================================================================
      // p_0/k_B [K/rad^3]
      // p_1     [degrees]
      // p_2     [degrees]
      U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi;
      fprintf(OutputFilePtr[CurrentSystem],"%4d CVFF_BEND_TORSION_CROSS %s, p_0/k_B=%8.5f [K/rad^3], p_1=%8.5f [degrees], p_2=%8.5f [degrees], "
                                           "Phi: %8.5f [degrees], Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2]*RAD2DEG,
                Phi*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_BEND_TORSION_CROSS:
      // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
      // =====================================================================================
      // p_0/k_B [K/rad^3]
      // p_1     [degrees]
      // p_2     [degrees]
      U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi;
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_BEND_TORSION_CROSS %s, p_0/k_B=%8.5f [K/rad^2], p_1=%8.5f [degrees], p_2=%8.5f [degrees], "
                                           "Phi: %8.5f [degrees], Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2]*RAD2DEG,
                Phi*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SMOOTHED_DIHEDRAL:
      // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2)]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [-]
      // p_2     [degrees]
      U=parms[0]*(1.0+cos(parms[2]*Phi-parms[1]))*Smoothing(Theta1)*Smoothing(Theta2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_DIHEDRAL %s, p_0/k_B=%8.5f [K/rad^2], p_1=%8.5f [-], p_2=%8.5f [degrees], "
                                           "Phi: %8.5f [degrees], Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*RAD2DEG,
                Phi*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SMOOTHED_THREE_COSINE_DIHEDRAL:
      // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      U=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
             Smoothing(Theta1)*Smoothing(Theta2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_THREE_COSINE_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "Phi: %8.5f [degrees], Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1],
                parms[2]*RAD2DEG,
                Phi*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case NICHOLAS_DIHEDRAL:
      // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      U=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
             Smoothing(Theta1);
      fprintf(OutputFilePtr[CurrentSystem],"%4d NICHOLAS_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "Phi: %8.5f [degrees], Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SMOOTHED_CFF_DIHEDRAL:
      // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      U=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_CFF_DIHEDRAL %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "Phi: %8.5f [degrees], Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SMOOTHED_CFF_DIHEDRAL2:
      // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      U=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*Smoothing(Theta1)*Smoothing(Theta2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_CFF_DIHEDRAL2 %s, p_0/k_B=%8.5f [K], p_1/k_B=%8.5f [K], p_2/k_B=%8.5f [K], "
                                           "Phi: %8.5f [degrees], Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*ENERGY_TO_KELVIN,
                parms[2]*ENERGY_TO_KELVIN,
                Phi*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SMOOTHED_CFF_BEND_TORSION_CROSS:
      // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K/rad^3]
      // p_1     [degrees]
      // p_2     [degrees]
      U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi*Smoothing(Theta1)*Smoothing(Theta2);
      fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_CFF_BEND_TORSION_CROSS %s, p_0/k_B=%8.5f [K], p_1=%8.5f [degrees], p_2=%8.5f [degrees], "
                                           "Phi: %8.5f [degrees], Theta1: %8.5f [degrees], Theta2: %8.5f [degrees], Energy: %8.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
                nr,
                string,
                parms[0]*ENERGY_TO_KELVIN,
                parms[1]*RAD2DEG,
                parms[2]*RAD2DEG,
                Phi*RAD2DEG,
                Theta1*RAD2DEG,
                Theta2*RAD2DEG,
                U*ENERGY_TO_KELVIN,
                U*ENERGY_TO_KJ_PER_MOL,
                U*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Undefined Bend-Torsion potential in routine 'CalculateFrameworkBendTorsionEnergy' ('status.c')\n");
      exit(0);
      break;
  }
  return U;
}

void PrintVDWEnergyStatus(int nr,char *string,int typeA,int typeB,REAL r,REAL energy)
{
  REAL arg1,arg2,arg3,arg4,arg5,arg6,arg7;
  REAL f6,f8,f10;

  switch(PotentialType[typeA][typeB])
  {
    case UNDEFINED_POTENTIAL:
    case ZERO_POTENTIAL:
      fprintf(OutputFilePtr[CurrentSystem],"%4d ZERO_POTENTIAL %s\n",
          nr++,string);
      break;
    case LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]
      // p_2     [A]
      // p_3/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      fprintf(OutputFilePtr[CurrentSystem],"%4d LENNARD_JONES %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A], shift/k_B=%8.5f [K], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case LENNARD_JONES_SMOOTHED3:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d LENNARD_JONES_SMOOTHED3 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case LENNARD_JONES_SMOOTHED5:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d LENNARD_JONES_SMOOTHED5 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+4 [p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2]*(p_2/T)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      // p_2     [A^2]
      // p_3/k_B [K]  (non-zero for a shifted potential)
      // T is the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      fprintf(OutputFilePtr[CurrentSystem],"%4d FEYNMAN_HIBBS_LENNARD_JONES %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A], p_2=%8.5f [A^2], shift/k_B=%8.5f [K]"
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3,
          arg4*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+4 [p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2]*(p_2/T)*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      // p_2     [A^2]
      // T is the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      fprintf(OutputFilePtr[CurrentSystem],"%4d FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A], p_2=%8.5f [A^2]"
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+4 [p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2]*(p_2/T)*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      // p_2     [A^2]
      // T is the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      fprintf(OutputFilePtr[CurrentSystem],"%4d FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A], p_2=%8.5f [A^2]"
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case LENNARD_JONES_SHIFTED_FORCE:
      // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]-[(p_1/rc)^12-(p_1/rc)^6]}+[12*(p_1/rc)^12-6*(p_1/rc)^6]*(r-rc)/rc
      // ===============================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d LENNARD_JONES_SHIFTED_FORCE %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case LENNARD_JONES_SHIFTED_FORCE2:
      // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]+[6*(p_1/rc)^12-3*(p_1/rc)^6]}*r^2/rc^2+[7*(p_1/rc)^12+4*(p_1/rc)^6]
      // =================================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d LENNARD_JONES_SHIFTED_FORCE2 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case POTENTIAL_12_6:
      // p_0/r^12-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      fprintf(OutputFilePtr[CurrentSystem],"%4d POTENTIAL_12_6 %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [K A^6], shift/k_B=%8.5f [K], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case POTENTIAL_12_6_SMOOTHED3:
      // {p_0/r^12-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d POTENTIAL_12_6_SMOOTHED3 %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [K A^6], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case POTENTIAL_12_6_SMOOTHED5:
      // {p_0/r^12-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d POTENTIAL_12_6_SMOOTHED5 %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [K A^6], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case POTENTIAL_12_6_2_0:
      // p_0/r^12+p_1/r^6+p_2/r^2+p_3
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      fprintf(OutputFilePtr[CurrentSystem],"%4d POTENTIAL_12_6_2_0 %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [K A^6], p_2/k_B=%8.5f [K A^2], p_3/k_B=%8.5f [K], "
                                           "shift/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case POTENTIAL_12_6_2_0_SMOOTHED3:
      // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      fprintf(OutputFilePtr[CurrentSystem],"%4d POTENTIAL_12_6_2_0_SMOOTHED3 %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [K A^6], p_2/k_B=%8.5f [K A^2], p_3/k_B=%8.5f [K], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case POTENTIAL_12_6_2_0_SMOOTHED5:
      // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      fprintf(OutputFilePtr[CurrentSystem],"%4d POTENTIAL_12_6_2_0_SMOOTHED5 %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [K A^6], p_2/k_B=%8.5f [K A^2], p_3/k_B=%8.5f [K], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_9_6:
      // p_0/r^9-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_9_6 %s, p_0/k_B=%8.5f [K A^9], p_1/k_B=%8.5f [K A^6], shift/k_B=%8.5f [K], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_9_6_SMOOTHED3:
      // {p_0/r^9-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_9_6_SMOOTHED3 %s, p_0/k_B=%8.5f [K A^9], p_1/k_B=%8.5f [K A^6], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_9_6_SMOOTHED5:
      // {p_0/r^9-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_9_6_SMOOTHED5 %s, p_0/k_B=%8.5f [K A^9], p_1/k_B=%8.5f [K A^6], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_EPS_SIGMA:
      // p_0*[2*(p_1/r)^9-3(p_1/r)^6]
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_EPS_SIGMA %s, p_0/k_B=%8.5f [K A^9], p_1=%8.5f [A], shift/k_B=%8.5f [K], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_EPS_SIGMA_SMOOTHED3:
      // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_EPS_SIGMA_SMOOTHED3 %s, p_0/k_B=%8.5f [K A^9], p_1=%8.5f [A], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case CFF_EPS_SIGMA_SMOOTHED5:
      // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d CFF_EPS_SIGMA_SMOOTHED5 %s, p_0/k_B=%8.5f [K A^9], p_1=%8.5f [A], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case BUCKINGHAM:
      // p_0*exp(-p_1*r)-p_2/r^6
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      fprintf(OutputFilePtr[CurrentSystem],"%4d BUCKINGHAM %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^6], shift/k_B=%8.5f [K], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case BUCKINGHAM_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      fprintf(OutputFilePtr[CurrentSystem],"%4d BUCKINGHAM_SMOOTHED3 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^6], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case BUCKINGHAM_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      fprintf(OutputFilePtr[CurrentSystem],"%4d BUCKINGHAM_SMOOTHED5 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^6], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_VDW:
      // sqrt(p_0^i*p_0^j)*[1.84e-5*exp(-12/P)-2.25*P^6]  if P>=3.02
      // sqrt(p_0^i*p_0^j)*192.27*P^2                     if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      // p_3     [kcal/mol]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_VDW %s, p_0=%8.5f [kcal/mol], p_1=%8.5f [A], shift=%8.5f [kcal/mol], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1/KCAL_PER_MOL_TO_ENERGY,
          arg2,
          arg3/KCAL_PER_MOL_TO_ENERGY,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_VDW_SMOOTHED3:
      // {sqrt(p_0^i*p_0^j)*[1.84e-5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
      // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                     if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      if(r<CutOffVDWSwitch)
        fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_VDW_SMOOTHED3 %s, p_0=%8.5f [kcal/mol], p_1=%8.5f [A], "
                                             "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
            nr++,
            string,
            arg1/KCAL_PER_MOL_TO_ENERGY,
            arg2,
            r,
            energy*ENERGY_TO_KELVIN,
            energy*ENERGY_TO_KJ_PER_MOL,
            energy*ENERGY_TO_KCAL_PER_MOL);
         else
        fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_VDW_SMOOTHED3 (Switching) %s, p_0=%8.5f [kcal/mol], p_1=%8.5f [A], "
                                             "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
            nr++,
            string,
            arg1/KCAL_PER_MOL_TO_ENERGY,
            arg2,
            r,
            energy*ENERGY_TO_KELVIN,
            energy*ENERGY_TO_KJ_PER_MOL,
            energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MM3_VDW_SMOOTHED5:
      // {sqrt(p_0^i*p_0^j)*[1.84e-5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
      // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                     if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      if(r<CutOffVDWSwitch)
        fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_VDW_SMOOTHED5 %s, p_0=%8.5f [kcal/mol], p_1=%8.5f [A], "
                                             "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
            nr++,
            string,
            arg1/KCAL_PER_MOL_TO_ENERGY,
            arg2,
            r,
            energy*ENERGY_TO_KELVIN,
            energy*ENERGY_TO_KJ_PER_MOL,
            energy*ENERGY_TO_KCAL_PER_MOL);
         else
        fprintf(OutputFilePtr[CurrentSystem],"%4d MM3_VDW_SMOOTHED5 (Switching) %s, p_0=%8.5f [kcal/mol], p_1=%8.5f [A], "
                                             "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
            nr++,
            string,
            arg1/KCAL_PER_MOL_TO_ENERGY,
            arg2,
            r,
            energy*ENERGY_TO_KELVIN,
            energy*ENERGY_TO_KJ_PER_MOL,
            energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MATSUOKA_CLEMENTI_YOSHIMINE:
      // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      fprintf(OutputFilePtr[CurrentSystem],"%4d MATSUOKA_CLEMENTI_YOSHIMINE %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K], p_3=%8.5f [A^-1], "
                                           "shift/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4,
          arg5*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3:
      // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      fprintf(OutputFilePtr[CurrentSystem],"%4d MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K], p_3=%8.5f [A^-1], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5:
      // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      fprintf(OutputFilePtr[CurrentSystem],"%4d MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K], p_3=%8.5f [A^-1], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case GENERIC:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      // p_6/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      arg7=PotentialParms[typeA][typeB][6];
      fprintf(OutputFilePtr[CurrentSystem],"%4d GENERIC %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^4], p_3/k_B=%8.5f [K A^6], p_4/k_B=%8.5f [K A^8],"
                                           "p_5/k_B=%8.5f [K A^10], shift/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          arg6*ENERGY_TO_KELVIN,
          arg7*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case GENERIC_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      fprintf(OutputFilePtr[CurrentSystem],"%4d GENERIC_SMOOTHED3 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^4], p_3/k_B=%8.5f [K A^6], p_4/k_B=%8.5f [K A^8],"
                                           "p_5/k_B=%8.5f [K A^10], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          arg6*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case GENERIC_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      fprintf(OutputFilePtr[CurrentSystem],"%4d GENERIC_SMOOTHED5 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^4], p_3/k_B=%8.5f [K A^6], p_4/k_B=%8.5f [K A^8],"
                                           "p_5/k_B=%8.5f [K A^10], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          arg6*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case PELLENQ_NICHOLSON:
      // p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      fprintf(OutputFilePtr[CurrentSystem],"%4d PELLENQ_NICHOLSON %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^6], p_3/k_B=%8.5f [K A^8], p_4/k_B=%8.5f [K A^10], "
                                           "f_6/k_B=%8.5f [K], f_8/k_B=%8.5f [K], f_10/k_B=%8.5f [K], shift/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          f6*ENERGY_TO_KELVIN,
          f8*ENERGY_TO_KELVIN,
          f10*ENERGY_TO_KELVIN,
          arg6*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case PELLENQ_NICHOLSON_SMOOTHED3:
      // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      fprintf(OutputFilePtr[CurrentSystem],"%4d PELLENQ_NICHOLSON_SMOOTHED3 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^6], p_3/k_B=%8.5f [K A^8], p_4/k_B=%8.5f [K A^10], "
                                           "f_6/k_B=%8.5f [K], f_8/k_B=%8.5f [K], f_10/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          f6*ENERGY_TO_KELVIN,
          f8*ENERGY_TO_KELVIN,
          f10*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case PELLENQ_NICHOLSON_SMOOTHED5:
      // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      fprintf(OutputFilePtr[CurrentSystem],"%4d PELLENQ_NICHOLSON_SMOOTHED5 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^6], p_3/k_B=%8.5f [K A^8], p_4/k_B=%8.5f [K A^10], "
                                           "f_6/k_B=%8.5f [K], f_8/k_B=%8.5f [K], f_10/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          f6*ENERGY_TO_KELVIN,
          f8*ENERGY_TO_KELVIN,
          f10*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HYDRATED_ION_WATER:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      fprintf(OutputFilePtr[CurrentSystem],"%4d HYDRATED_ION_WATER %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^4], p_3/k_B=%8.5f [K A^6], "
                                           "p_4/k_B=%8.5f [K A^12], shift/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          arg6*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HYDRATED_ION_WATER_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      fprintf(OutputFilePtr[CurrentSystem],"%4d HYDRATED_ION_WATER_SMOOTHED3 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^4], p_3/k_B=%8.5f [K A^6], "
                                           "p_4/k_B=%8.5f [K A^12], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HYDRATED_ION_WATER_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      fprintf(OutputFilePtr[CurrentSystem],"%4d HYDRATED_ION_WATER_SMOOTHED5 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^-1], p_2/k_B=%8.5f [K A^4], p_3/k_B=%8.5f [K A^6], "
                                           "p_4/k_B=%8.5f [K A^12], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3*ENERGY_TO_KELVIN,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MIE:
      // p_0*[p_1/r^p_2-p_1/r^p_3]
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K A^p_2]
      // p_2     [-]
      // p_3     [-]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      fprintf(OutputFilePtr[CurrentSystem],"%4d MIE %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^p_2], p_2=%8.5f [-], p_3=%8.5f [-], "
                                           "shift/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3,
          arg4,
          arg5*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MIE_SMOOTHED3:
      // {p_0*[p_1/r^p_2-p_1/r^p_3]}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K A^p_2]
      // p_2     [-]
      // p_3     [-]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      fprintf(OutputFilePtr[CurrentSystem],"%4d MIE_SMOOTHED3 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^p_2], p_2=%8.5f [-], p_3=%8.5f [-], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3,
          arg4,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MIE_SMOOTHED5:
      // {p_0*[p_1/r^p_2-p_1/r^p_3]}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K A^p_2]
      // p_2     [-]
      // p_3     [-]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      fprintf(OutputFilePtr[CurrentSystem],"%4d MIE_SMOOTHED5 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^p_2], p_2=%8.5f [-], p_3=%8.5f [-], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3,
          arg4,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MIE_CUTOFF:
      // p_0*[p_1/r^p_2-p_1/r^p_3] if r < p_4, otherwise 0
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K A^p_2]
      // p_2     [-]
      // p_3     [-]
      // p_4     [-]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      fprintf(OutputFilePtr[CurrentSystem],"%4d MIE_CUTOFF %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^p_2], p_2=%8.5f [-], p_3=%8.5f [-], p_4=%8.5f [A],"
                                           "shift/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3,
          arg4,
          arg5,
          arg6*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MIE_SMOOTHED3_CUTOFF:
      // {p_0*[p_1/r^p_2-p_1/r^p_3]}*S(r) if r < p_4, otherwise 0
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K A^p_2]
      // p_2     [-]
      // p_3     [-]
      // p_4     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      fprintf(OutputFilePtr[CurrentSystem],"%4d MIE_SMOOTHED3_CUTOFF %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^p_2], p_2=%8.5f [-], p_3=%8.5f [-], p_4=%8.5f [A], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3,
          arg4,
          arg5,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case MIE_SMOOTHED5_CUTOFF:
      // {p_0*[p_1/r^p_2-p_1/r^p_3]}*S(r) if r < p_4, otherwise 0
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K A^p_2]
      // p_2     [-]
      // p_3     [-]
      // p_4     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      fprintf(OutputFilePtr[CurrentSystem],"%4d MIE_SMOOTHED5_CUTOFF %s, p_0/k_B=%8.5f [K], p_1=%8.5f [A^p_2], p_2=%8.5f [-], p_3=%8.5f [-], p_4=%8.5f [A], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3,
          arg4,
          arg5,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case BORN_HUGGINS_MEYER:
      // p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      fprintf(OutputFilePtr[CurrentSystem],"%4d BORN_HUGGINS_MEYER %s, p_0/k_B=%8.5f [K], p_1=%8.5f [-], p_2=%8.5f [A], p_3/k_B=%8.5f [K A^6], "
                                           "p_4/k_B=%8.5f [K A^8], shift/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          arg6*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case BORN_HUGGINS_MEYER_SMOOTHED3:
      // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      fprintf(OutputFilePtr[CurrentSystem],"%4d BORN_HUGGINS_MEYER_SMOOTHED3 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [-], p_2=%8.5f [A], p_3/k_B=%8.5f [K A^6], "
                                           "p_4/k_B=%8.5f [K A^8], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case BORN_HUGGINS_MEYER_SMOOTHED5:
      // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      fprintf(OutputFilePtr[CurrentSystem],"%4d BORN_HUGGINS_MEYER_SMOOTHED5 %s, p_0/k_B=%8.5f [K], p_1=%8.5f [-], p_2=%8.5f [A], p_3/k_B=%8.5f [K A^6], "
                                           "p_4/k_B=%8.5f [K A^8], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2,
          arg3,
          arg4*ENERGY_TO_KELVIN,
          arg5*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HYDROGEN:
      // p_0/r^12-p_1/r^10
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      fprintf(OutputFilePtr[CurrentSystem],"%4d HYDROGEN %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [A A^10], "
                                           "shift/k_B=%8.5f [K], Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          arg3*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HYDROGEN_SMOOTHED3:
      // {p_0/r^12-p_1/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d HYDROGEN_SMOOTHED3 %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [A A^10], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case HYDROGEN_SMOOTHED5:
      // {p_0/r^12-p_1/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      fprintf(OutputFilePtr[CurrentSystem],"%4d HYDROGEN_SMOOTHED5 %s, p_0/k_B=%8.5f [K A^12], p_1/k_B=%8.5f [A A^10], "
                                           "Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr++,
          string,
          arg1*ENERGY_TO_KELVIN,
          arg2*ENERGY_TO_KELVIN,
          r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Undefined VDW potential in routine 'PrintVDWEnergyStatus' ('status.c')\n");
      exit(0);
      break;
  }
}


REAL PrintChargeChargeEnergyStatus(int nr,char *string,REAL chargeA,REAL chargeB,REAL r)
{
  REAL energy,rr;
  REAL SwitchingValue;
  REAL TranslationValue;

  rr=SQR(r);
  switch(ChargeMethod)
  {
    case NONE:
      energy=0.0;
      break;
    case TRUNCATED_COULOMB:
      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_COULOMB %s, Charge: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,chargeA,chargeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);

      break;
    case SHIFTED_COULOMB:
      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge[CurrentSystem]);
      fprintf(OutputFilePtr[CurrentSystem],"%4d SHIFTED_COULOMB %s, Charge: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,chargeA,chargeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SMOOTHED_COULOMB:
      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch[CurrentSystem]+CutOffChargeCharge[CurrentSystem]));
      if(rr>CutOffChargeChargeSwitchSquared[CurrentSystem])
      {
        SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                       SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

        TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                         (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                          SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                          SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
        energy=energy*SwitchingValue+TranslationValue;
        fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_COULOMB (Switching) %s, Charge: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
            nr,string,chargeA,chargeB,r,
            energy*ENERGY_TO_KELVIN,
            energy*ENERGY_TO_KJ_PER_MOL,
            energy*ENERGY_TO_KCAL_PER_MOL);
      }
      else
      {
        fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_COULOMB %s, Charge: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
            nr,string,chargeA,chargeB,r,
            energy*ENERGY_TO_KELVIN,
            energy*ENERGY_TO_KJ_PER_MOL,
            energy*ENERGY_TO_KCAL_PER_MOL);
      }
      break;
    case WOLFS_METHOD_DAMPED_FG:
      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
               -erfc(Alpha[CurrentSystem]*CutOffChargeCharge[CurrentSystem])*InverseCutOffChargeCharge[CurrentSystem]+
                (r-CutOffChargeCharge[CurrentSystem])*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge[CurrentSystem])*SQR(InverseCutOffChargeCharge[CurrentSystem])+
                (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared[CurrentSystem])*M_1_SQRTPI*InverseCutOffChargeCharge[CurrentSystem])));
      fprintf(OutputFilePtr[CurrentSystem],"%4d WOLFS_METHOD_DAMPED_FG %s, Charge: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,chargeA,chargeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case EWALD:
      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                       (erfc(Alpha[CurrentSystem]*r)/r);
      fprintf(OutputFilePtr[CurrentSystem],"%4d EWALD %s, Charge: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,chargeA,chargeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Unknown charge-charge method in 'CalculateFrameworkIntraChargeChargeEnergy'\n");
      exit(0);
      break;
  }

  return energy;
}

REAL PrintChargeBondDipoleStatus(int nr,char *string,REAL cosA,REAL magntitude, REAL ChargeB,REAL r)
{
  REAL energy,rr;
  REAL SwitchingValue;
  REAL Bt0,Bt1,Bt2;

  rr=SQR(r);

  switch(ChargeMethod)
  {
    case NONE:
      Bt0=Bt1=Bt2=0.0;
      energy=0.0;
      fprintf(OutputFilePtr[CurrentSystem],"%4d NONE %s, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case TRUNCATED_COULOMB:
      Bt0=1.0/(r);
      Bt1=1.0/(r*rr);
      Bt2=3.0/(r*rr*rr);
      energy=-COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_COULOMB %s, Magnitude: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,magntitude,ChargeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);

      break;
    case SHIFTED_COULOMB:
      Bt0=1.0/(r);
      Bt1=1.0/(r*rr);
      Bt2=3.0/(r*rr*rr);
      energy=-COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
      fprintf(OutputFilePtr[CurrentSystem],"%4d SHIFTED_COULOMB %s, Magnitude: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,magntitude,ChargeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SMOOTHED_COULOMB:
      Bt0=1.0/(r);
      Bt1=1.0/(r*rr);
      Bt2=3.0/(r*rr*rr);
      energy=-COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
      if(rr>CutOffChargeBondDipoleSwitchSquared)
      {
        SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                       SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
        energy*=SwitchingValue;
        fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_COULOMB (Switching) %s, Magnitude: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
            nr,string,magntitude,ChargeB,r,
            energy*ENERGY_TO_KELVIN,
            energy*ENERGY_TO_KJ_PER_MOL,
            energy*ENERGY_TO_KCAL_PER_MOL);
      }
      else
      {
        fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_COULOMB %s, Magnitude: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
            nr,string,magntitude,ChargeB,r,
            energy*ENERGY_TO_KELVIN,
            energy*ENERGY_TO_KJ_PER_MOL,
            energy*ENERGY_TO_KCAL_PER_MOL);
      }
      break;
    case EWALD:
      Bt0=erfc(Alpha[CurrentSystem]*r)/r;
      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          erfc(Alpha[CurrentSystem]*r)/(rr*r);
      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
      energy=-COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
      fprintf(OutputFilePtr[CurrentSystem],"%4d EWALD %s, Magnitude: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,magntitude,ChargeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    default:
      fprintf(stderr, "Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombEnergy'\n");
      exit(0);
      break;
  }
  return energy;
}

REAL PrintBondDipoleBondDipoleStatus(int nr,char *string,REAL cosA,REAL cosB,REAL cosAB,REAL magntitudeA,REAL magnitudeB,REAL r)
{
  REAL Bt1,Bt2,rr;
  REAL SwitchingValue;
  REAL energy;

  rr=SQR(r);

  switch(ChargeMethod)
  {
    case NONE:
      energy=0.0;
      Bt1=Bt2=0.0;
      fprintf(OutputFilePtr[CurrentSystem],"%4d NONE %s, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case TRUNCATED_COULOMB:
      Bt1=1.0/(rr*r);
      Bt2=3.0/(SQR(rr)*r);
      energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_COULOMB %s, Magnitude: %8.5f, Magnitude: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,magntitudeA,magnitudeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SHIFTED_COULOMB:
      Bt1=1.0/(rr*r);
      Bt2=3.0/(SQR(rr)*r);
      energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_COULOMB %s, Magnitude: %8.5f, Magnitude: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,magntitudeA,magnitudeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);
      break;
    case SMOOTHED_COULOMB:
      Bt1=1.0/(r*rr);
      Bt2=3.0/(r*rr*rr);
      energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
      if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
      {
        SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                       SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
        energy*=SwitchingValue;
        fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_COULOMB (Switching) %s, Magnitude: %8.5f, Magnitude: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
            nr,string,magntitudeA,magnitudeB,r,
            energy*ENERGY_TO_KELVIN,
            energy*ENERGY_TO_KJ_PER_MOL,
            energy*ENERGY_TO_KCAL_PER_MOL);

      }
      else
      {
        fprintf(OutputFilePtr[CurrentSystem],"%4d SMOOTHED_COULOMB %s, Magnitude: %8.5f, Magnitude: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
            nr,string,magntitudeA,magnitudeB,r,
            energy*ENERGY_TO_KELVIN,
            energy*ENERGY_TO_KJ_PER_MOL,
            energy*ENERGY_TO_KCAL_PER_MOL);
      }
      break;
    case EWALD:
      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          erfc(Alpha[CurrentSystem]*r)/(rr*r);
      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
      energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
      break;
    default:
      fprintf(stderr, "Unknown bonddipole-bonddipole method in 'CalculateFrameworkIntraBondDipoleBondDipoleEnergy'\n");
      exit(0);
      break;
  }
  return energy;
}


// =============================================================================================================================
// framework intra-molecular energies
// =============================================================================================================================

void PrintFrameworkBondEnergyStatus(void)
{
  int i,f1;
  REAL r,rr,U;
  REAL *parms,UHostBond;
  VECTOR posA,posB,dr;
  char string[256];
  int A,B,nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostBond:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostBond=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBonds[f1];i++)
    {
      A=Framework[CurrentSystem].Bonds[f1][i].A;
      B=Framework[CurrentSystem].Bonds[f1][i].B;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      parms=(REAL*)&Framework[CurrentSystem].BondArguments[f1][i];

      sprintf(string,"(F%1d,[%-4d,%-4d]) (%-6s,%-6s)",f1,A,B,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name);

      U=PrintBondEnergyStatus(nr++,string,Framework[CurrentSystem].BondType[f1][i],parms,r);

      // add contribution to the Adsorbate stretch energy
      UHostBond+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework bonds: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostBond Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostBond*ENERGY_TO_KELVIN,
          UHostBond*ENERGY_TO_KJ_PER_MOL,
          UHostBond*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}

void PrintFrameworkUreyBradleyEnergyStatus(void)
{
  int i,f1;
  REAL r,rr,U;
  VECTOR dr,posA,posC;
  char string[256];
  REAL *parms,UHostUreyBradley;
  int A,C,nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostUreyBradley:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostUreyBradley=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfUreyBradleys[f1];i++)
    {
      A=Framework[CurrentSystem].UreyBradleys[f1][i].A;
      C=Framework[CurrentSystem].UreyBradleys[f1][i].C;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;

      dr.x=posA.x-posC.x;
      dr.y=posA.y-posC.y;
      dr.z=posA.z-posC.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      parms=(REAL*)&Framework[CurrentSystem].UreyBradleyArguments[f1][i];

      sprintf(string,"(F%1d,[%-4d,%-4d]) (%-6s,%-6s)",f1,A,C,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name);
      U=PrintUreyBradleyEnergyStatus(nr++,string,Framework[CurrentSystem].UreyBradleyType[f1][i],parms,r);

      // add contribution to the Adsorbate Urey-Bradley energy
      UHostUreyBradley+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework Urey-Bradleys: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostUreyBradley Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostUreyBradley*ENERGY_TO_KELVIN,
          UHostUreyBradley*ENERGY_TO_KJ_PER_MOL,
          UHostUreyBradley*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}



void PrintFrameworkBendEnergyStatus(void)
{
  int i,A,B,C,D,f1,nr;
  REAL *parms,U;
  REAL CosTheta,Theta;
  REAL rab,rbc,rac,UHostBend;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL delta,rt2,rap2,rcp2;
  char string[256];

  fprintf(OutputFilePtr[CurrentSystem],"UHostBend:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostBend=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBends[f1];i++)
    {
      A=Framework[CurrentSystem].Bends[f1][i].A;
      B=Framework[CurrentSystem].Bends[f1][i].B;
      C=Framework[CurrentSystem].Bends[f1][i].C;
      D=-1;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;

      switch(Framework[CurrentSystem].BendType[f1][i])
      {
        case MM3_IN_PLANE_BEND:
          D=Framework[CurrentSystem].Bends[f1][i].D;
          posD=Framework[CurrentSystem].Atoms[f1][D].Position;
          Rad.x=posA.x-posD.x;
          Rad.y=posA.y-posD.y;
          Rad.z=posA.z-posD.z;
          Rad=ApplyBoundaryCondition(Rad);

          Rbd.x=posB.x-posD.x;
          Rbd.y=posB.y-posD.y;
          Rbd.z=posB.z-posD.z;
          Rbd=ApplyBoundaryCondition(Rbd);

          Rcd.x=posC.x-posD.x;
          Rcd.y=posC.y-posD.y;
          Rcd.z=posC.z-posD.z;
          Rcd=ApplyBoundaryCondition(Rcd);

          t.x=Rad.y*Rcd.z-Rad.z*Rcd.y;
          t.y=Rad.z*Rcd.x-Rad.x*Rcd.z;
          t.z=Rad.x*Rcd.y-Rad.y*Rcd.x;
          rt2=t.x*t.x+t.y*t.y+t.z*t.z;
          delta=-(t.x*Rbd.x+t.y*Rbd.y+t.z*Rbd.z)/rt2;

          ip.x=posB.x+t.x*delta;
          ip.y=posB.y+t.y*delta;
          ip.z=posB.z+t.z*delta;
          ap.x=posA.x-ip.x;
          ap.y=posA.y-ip.y;
          ap.z=posA.z-ip.z;
          cp.x=posC.x-ip.x;
          cp.y=posC.y-ip.y;
          cp.z=posC.z-ip.z;
          ap=ApplyBoundaryCondition(ap);
          cp=ApplyBoundaryCondition(cp);

          rap2=ap.x*ap.x+ap.y*ap.y+ap.z*ap.z;
          rcp2=cp.x*cp.x+cp.y*cp.y+cp.z*cp.z;

          CosTheta=(ap.x*cp.x+ap.y*cp.y+ap.z*cp.z)/sqrt(rap2*rcp2);
          break;
        default:
          Rab.x=posA.x-posB.x;
          Rab.y=posA.y-posB.y;
          Rab.z=posA.z-posB.z;
          Rab=ApplyBoundaryCondition(Rab);
          rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
          Rab.x/=rab;
          Rab.y/=rab;
          Rab.z/=rab;

          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;
          Rbc=ApplyBoundaryCondition(Rbc);
          rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
          Rbc.x/=rbc;
          Rbc.y/=rbc;
          Rbc.z/=rbc;

          Rac.x=posC.x-posA.x;
          Rac.y=posC.y-posA.y;
          Rac.z=posC.z-posA.z;
          Rac=ApplyBoundaryCondition(Rac);
          rac=sqrt(SQR(Rac.x)+SQR(Rac.y)+SQR(Rac.z));
          Rac.x/=rac;
          Rac.y/=rac;
          Rac.z/=rac;

          CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
          break;
      }

      CosTheta=MIN2(1.0,MAX2(-1.0,CosTheta));
      Theta=acos(CosTheta);

      parms=Framework[CurrentSystem].BendArguments[f1][i];

      switch(Framework[CurrentSystem].BendType[f1][i])
      {
         case MM3_IN_PLANE_BEND:
           sprintf(string,"(F%1d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",f1,A,B,C,D,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][D].Type].Name);
           break;
         default:
           sprintf(string,"(F%1d,[%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s)",f1,A,B,C,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name,
                     PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name);
           break;
      }

      U=PrintBendEnergyStatus(nr++,string,Framework[CurrentSystem].BendType[f1][i],parms,Theta);

      // add contribution to the energy
      UHostBend+=U;
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework bends: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostBend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostBend*ENERGY_TO_KELVIN,
          UHostBend*ENERGY_TO_KJ_PER_MOL,
          UHostBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

// Out-of-plane angle  A-B<(C,D)
// ============================
// R.E. Tuzun, D.W. Noid, B.G. Sumpter
// Journal of computational chemistry, vol 18(14), pages 1804 -1811, 1997
// Efficient Treatment of Out-of-Plane Bend and Improper Torsion Interactions in MM2, MM3, and MM4 Molecular Mechanics Calculations
//
// MM3: angle between the C-D-A plane and the A-B vector
// others: angle between the C-B-D plane and the B-A vector

void PrintFrameworkInversionBendEnergyStatus(void)
{
  int i,A,B,C,D,f1;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2;
  REAL CosChi,Chi,energy;
  REAL UHostInversionBend;
  VECTOR Rab,Rbc,Rbd,Rcd,Rad;
  POINT posA,posB,posC,posD;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostInversionBend\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostInversionBend=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfInversionBends[f1];i++)
    {
      A=Framework[CurrentSystem].InversionBends[f1][i].A;
      B=Framework[CurrentSystem].InversionBends[f1][i].B;
      C=Framework[CurrentSystem].InversionBends[f1][i].C;
      D=Framework[CurrentSystem].InversionBends[f1][i].D;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD=Framework[CurrentSystem].Atoms[f1][D].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      Rab=ApplyBoundaryCondition(Rab);
      rab2=Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z;
      rrab=sqrt(rab2);

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      Rbc=ApplyBoundaryCondition(Rbc);

      Rbd.x=posD.x-posB.x;
      Rbd.y=posD.y-posB.y;
      Rbd.z=posD.z-posB.z;
      Rbd=ApplyBoundaryCondition(Rbd);

      Rcd.x=posD.x-posC.x;
      Rcd.y=posD.y-posC.y;
      Rcd.z=posD.z-posC.z;
      Rcd=ApplyBoundaryCondition(Rcd);

      Rad.x=posD.x-posA.x;
      Rad.y=posD.y-posA.y;
      Rad.z=posD.z-posA.z;
      Rad=ApplyBoundaryCondition(Rad);

      switch(Framework[CurrentSystem].InversionBendType[f1][i])
      {
        case HARMONIC_INVERSION:
        case HARMONIC_COSINE_INVERSION:
        case PLANAR_INVERSION:
          // w is a vector perpendicular to the B-C-D plane
          // c=w.w=(Rbc x Rbd).(Rbc x Rbd)= r_bc^2 r_bd^2 - (r_cb . r_bd)^2
          c=(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z)*(Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z)-SQR(Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z);
          break;
        case HARMONIC_INVERSION2:
        case HARMONIC_COSINE_INVERSION2:
        case PLANAR_INVERSION2:
        case MM3_INVERSION:
          // w is a vector perpendicular to the A-C-D plane
          // c=w.w=(Rcd x Rad).(Rcd x Rad)=r_cd^2 r_ad^2 - (r_da . r_cd)^2
          c=(Rcd.x*Rcd.x+Rcd.y*Rcd.y+Rcd.z*Rcd.z)*(Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z)-SQR(Rad.x*Rcd.x+Rad.y*Rcd.y+Rad.z*Rcd.z);
          break;
        default:
          fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
          exit(0);
          break;
      }

      e=Rab.x*(Rbd.y*Rbc.z-Rbd.z*Rbc.y)+Rab.y*(Rbd.z*Rbc.x-Rbd.x*Rbc.z)+Rab.z*(Rbd.x*Rbc.y-Rbd.y*Rbc.x);
      CosChi=sqrt((Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z)-SQR(e)/c)/rrab;

      // Ensure CosChi is between -1 and 1.
      CosChi=SIGN(MIN2(fabs(CosChi),(REAL)1.0),CosChi);

      Chi=acos(CosChi);

      parms=Framework[CurrentSystem].InversionBendArguments[f1][i];

      sprintf(string,"(F%1d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",f1,A,B,C,D,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][D].Type].Name);

      energy=PrintInversionBendEnergyStatus(nr++,string,Framework[CurrentSystem].InversionBendType[f1][i],parms,Chi);

      // energy
      UHostInversionBend+=energy;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework inversion-bends: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostInversionBend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostInversionBend*ENERGY_TO_KELVIN,
          UHostInversionBend*ENERGY_TO_KJ_PER_MOL,
          UHostInversionBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintFrameworkTorsionEnergyStatus(void)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL rbc,U;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,UHostTorsion;
  VECTOR Pb,Pc;
  char string[256];
  REAL *parms;
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostTorsion\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  UHostTorsion=0.0;
  nr=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[f1];i++)
    {
      A=Framework[CurrentSystem].Torsions[f1][i].A;
      B=Framework[CurrentSystem].Torsions[f1][i].B;
      C=Framework[CurrentSystem].Torsions[f1][i].C;
      D=Framework[CurrentSystem].Torsions[f1][i].D;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD=Framework[CurrentSystem].Atoms[f1][D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      Dab=ApplyBoundaryCondition(Dab);

      Dcb.x=posC.x-posB.x;
      Dcb.y=posC.y-posB.y;
      Dcb.z=posC.z-posB.z;
      Dcb=ApplyBoundaryCondition(Dcb);
      rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
      Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

      Ddc.x=posD.x-posC.x;
      Ddc.y=posD.y-posC.y;
      Ddc.z=posD.z-posC.z;
      Ddc=ApplyBoundaryCondition(Ddc);

      dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
      dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

      dr.x=Dab.x-dot_ab*Dcb.x;
      dr.y=Dab.y-dot_ab*Dcb.y;
      dr.z=Dab.z-dot_ab*Dcb.z;
      r=MAX2((REAL)1.0e-8,sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z)));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Ddc.x-dot_cd*Dcb.x;
      ds.y=Ddc.y-dot_cd*Dcb.y;
      ds.z=Ddc.z-dot_cd*Dcb.z;
      s=MAX2((REAL)1.0e-8,sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z)));
      ds.x/=s; ds.y/=s; ds.z/=s;

      // compute Cos(Phi)
      // Phi is defined in protein convention Phi(trans)=Pi
      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      // Ensure CosPhi is between -1 and 1.
      CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
      CosPhi2=SQR(CosPhi);

      Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
      Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
      Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
      Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
      Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
      Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
      sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);

      parms=Framework[CurrentSystem].TorsionArguments[f1][i];

      sprintf(string,"(F%1d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",f1,A,B,C,D,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][D].Type].Name);

      U=PrintTorsionEnergyStatus(nr++,string,Framework[CurrentSystem].TorsionType[f1][i],parms,Phi);


      // energy
      UHostTorsion+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework torsions: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostTorsion*ENERGY_TO_KELVIN,
          UHostTorsion*ENERGY_TO_KJ_PER_MOL,
          UHostTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintFrameworkImproperTorsionEnergyStatus(void)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL rbc,U;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,UHostImproperTorsion;
  char string[256];
  VECTOR Pb,Pc;
  REAL *parms;
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostImproperTorsion\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  UHostImproperTorsion=0.0;
  nr=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfImproperTorsions[f1];i++)
    {
      A=Framework[CurrentSystem].ImproperTorsions[f1][i].A;
      B=Framework[CurrentSystem].ImproperTorsions[f1][i].B;
      C=Framework[CurrentSystem].ImproperTorsions[f1][i].C;
      D=Framework[CurrentSystem].ImproperTorsions[f1][i].D;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD=Framework[CurrentSystem].Atoms[f1][D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      Dab=ApplyBoundaryCondition(Dab);

      Dcb.x=posC.x-posB.x;
      Dcb.y=posC.y-posB.y;
      Dcb.z=posC.z-posB.z;
      Dcb=ApplyBoundaryCondition(Dcb);
      rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
      Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

      Ddc.x=posD.x-posC.x;
      Ddc.y=posD.y-posC.y;
      Ddc.z=posD.z-posC.z;
      Ddc=ApplyBoundaryCondition(Ddc);

      dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
      dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

      dr.x=Dab.x-dot_ab*Dcb.x;
      dr.y=Dab.y-dot_ab*Dcb.y;
      dr.z=Dab.z-dot_ab*Dcb.z;
      r=MAX2((REAL)1.0e-8,sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z)));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Ddc.x-dot_cd*Dcb.x;
      ds.y=Ddc.y-dot_cd*Dcb.y;
      ds.z=Ddc.z-dot_cd*Dcb.z;
      s=MAX2((REAL)1.0e-8,sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z)));
      ds.x/=s; ds.y/=s; ds.z/=s;

      // compute Cos(Phi)
      // Phi is defined in protein convention Phi(trans)=Pi
      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      // Ensure CosPhi is between -1 and 1.
      CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
      CosPhi2=SQR(CosPhi);

      Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
      Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
      Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
      Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
      Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
      Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
      sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);

      parms=Framework[CurrentSystem].ImproperTorsionArguments[f1][i];

      sprintf(string,"(F%1d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",f1,A,B,C,D,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][D].Type].Name);

      U=PrintImproperTorsionEnergyStatus(nr++,string,Framework[CurrentSystem].ImproperTorsionType[f1][i],parms,Phi);

      // energy
      UHostImproperTorsion+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework improper torsions: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostImproperTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostImproperTorsion*ENERGY_TO_KELVIN,
          UHostImproperTorsion*ENERGY_TO_KJ_PER_MOL,
          UHostImproperTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintFrameworkBondBondEnergyStatus(void)
{
  int i,A,B,C,f1;
  REAL *parms;
  REAL energy;
  REAL rab,rbc,UHostBondBond;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostBondBond\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostBondBond=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondBonds[f1];i++)
    {
      A=Framework[CurrentSystem].BondBonds[f1][i].A;
      B=Framework[CurrentSystem].BondBonds[f1][i].B;
      C=Framework[CurrentSystem].BondBonds[f1][i].C;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      Rab=ApplyBoundaryCondition(Rab);
      rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      Rbc=ApplyBoundaryCondition(Rbc);
      rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));

      parms=Framework[CurrentSystem].BondBondArguments[f1][i];

      sprintf(string,"(F%1d,[%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s)",f1,A,B,C,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name);

      energy=PrintBondBondEnergyStatus(nr++,string,Framework[CurrentSystem].BondBondType[f1][i],parms,rab,rbc);


      // add contribution to the energy
      UHostBondBond+=energy;
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework bond-bond cross terms: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostBondBond Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostBondBond*ENERGY_TO_KELVIN,
          UHostBondBond*ENERGY_TO_KJ_PER_MOL,
          UHostBondBond*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintFrameworkBondBendEnergyStatus(void)
{
  int i,A,B,C,f1;
  REAL *parms,U;
  REAL cost,Theta,sint;
  REAL rab,rbc,rac,UHostBondBend;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc,Rac;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostBondBend\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostBondBend=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondBends[f1];i++)
    {
      A=Framework[CurrentSystem].BondBends[f1][i].A;
      B=Framework[CurrentSystem].BondBends[f1][i].B;
      C=Framework[CurrentSystem].BondBends[f1][i].C;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      Rab=ApplyBoundaryCondition(Rab);
      rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
      Rab.x/=rab;
      Rab.y/=rab;
      Rab.z/=rab;

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      Rbc=ApplyBoundaryCondition(Rbc);
      rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
      Rbc.x/=rbc;
      Rbc.y/=rbc;
      Rbc.z/=rbc;

      Rac.x=posC.x-posA.x;
      Rac.y=posC.y-posA.y;
      Rac.z=posC.z-posA.z;
      Rac=ApplyBoundaryCondition(Rac);
      rac=sqrt(SQR(Rac.x)+SQR(Rac.y)+SQR(Rac.z));
      Rac.x/=rac;
      Rac.y/=rac;
      Rac.z/=rac;

      cost=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
      cost=SIGN(MIN2(fabs(cost),(REAL)1.0),cost);
      Theta=acos(cost);
      sint=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(cost)));

      parms=Framework[CurrentSystem].BondBendArguments[f1][i];

      sprintf(string,"(F%1d,[%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s)",f1,A,B,C,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name);

      U=PrintBondBendEnergyStatus(nr++,string,Framework[CurrentSystem].BondBendType[f1][i],parms,rab,rbc,Theta);

      // energy
      UHostBondBend+=U;
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework bond-bend cross terms: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostBondBend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostBondBend*ENERGY_TO_KELVIN,
          UHostBondBend*ENERGY_TO_KJ_PER_MOL,
          UHostBondBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

// Bend/Bend cross term for a centered atom B
// the first angle is A-B-C
// the second angle is A-B-D
void PrintFrameworkBendBendEnergyStatus()
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL U,rab,rbc,rbd;
  VECTOR Dab,Dbc,Dbd;
  REAL dot_abc,dot_abd,UHostBendBend;
  REAL CosTheta1,CosTheta2,Theta1,Theta2;
  char string[256];
  REAL *parms;
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostBendBend\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostBendBend=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBendBends[f1];i++)
    {
      A=Framework[CurrentSystem].BendBends[f1][i].A;
      B=Framework[CurrentSystem].BendBends[f1][i].B;
      C=Framework[CurrentSystem].BendBends[f1][i].C;
      D=Framework[CurrentSystem].BendBends[f1][i].D;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD=Framework[CurrentSystem].Atoms[f1][D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      Dab=ApplyBoundaryCondition(Dab);
      rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));
      Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      Dbc=ApplyBoundaryCondition(Dbc);
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dbd.x=posD.x-posB.x;
      Dbd.y=posD.y-posB.y;
      Dbd.z=posD.z-posB.z;
      Dbd=ApplyBoundaryCondition(Dbd);
      rbd=sqrt(SQR(Dbd.x)+SQR(Dbd.y)+SQR(Dbd.z));
      Dbd.x/=rbd; Dbd.y/=rbd; Dbd.z/=rbd;

      dot_abc=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      CosTheta1=dot_abc;
      CosTheta1=SIGN(MIN2(fabs(CosTheta1),(REAL)1.0),CosTheta1);
      Theta1=acos(CosTheta1);

      dot_abd=Dab.x*Dbd.x+Dab.y*Dbd.y+Dab.z*Dbd.z;
      CosTheta2=dot_abd;
      CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
      Theta2=acos(CosTheta2);

      parms=Framework[CurrentSystem].BendBendArguments[f1][i];

      sprintf(string,"(F%1d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",f1,A,B,C,D,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][D].Type].Name);

      U=PrintBendBendEnergyStatus(nr++,string,Framework[CurrentSystem].BendBendType[f1][i],parms,Theta1,Theta2);


      UHostBendBend+=U;
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework bend-bend cross terms: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostBendBend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostBendBend*ENERGY_TO_KELVIN,
          UHostBendBend*ENERGY_TO_KJ_PER_MOL,
          UHostBendBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintFrameworkBondTorsionEnergyStatus(void)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL U,rab,rbc,rcd;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL Phi,CosPhi,CosPhi2,sign;
  REAL *parms,UHostBondTorsion;
  char string[256];
  VECTOR Pb,Pc;
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostBondTorsion\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostBondTorsion=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondTorsions[f1];i++)
    {
      A=Framework[CurrentSystem].BondTorsions[f1][i].A;
      B=Framework[CurrentSystem].BondTorsions[f1][i].B;
      C=Framework[CurrentSystem].BondTorsions[f1][i].C;
      D=Framework[CurrentSystem].BondTorsions[f1][i].D;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD=Framework[CurrentSystem].Atoms[f1][D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      Dab=ApplyBoundaryCondition(Dab);
      rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

      Dcb.x=posC.x-posB.x;
      Dcb.y=posC.y-posB.y;
      Dcb.z=posC.z-posB.z;
      Dcb=ApplyBoundaryCondition(Dcb);
      rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
      Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

      Ddc.x=posD.x-posC.x;
      Ddc.y=posD.y-posC.y;
      Ddc.z=posD.z-posC.z;
      Ddc=ApplyBoundaryCondition(Ddc);
      rcd=sqrt(SQR(Ddc.x)+SQR(Ddc.y)+SQR(Ddc.z));

      dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
      dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

      dr.x=Dab.x-dot_ab*Dcb.x;
      dr.y=Dab.y-dot_ab*Dcb.y;
      dr.z=Dab.z-dot_ab*Dcb.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Ddc.x-dot_cd*Dcb.x;
      ds.y=Ddc.y-dot_cd*Dcb.y;
      ds.z=Ddc.z-dot_cd*Dcb.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      // compute Cos(Phi)
      // Phi is defined in protein convention Phi(trans)=Pi
      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      // Ensure CosPhi is between -1 and 1.
      CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
      CosPhi2=SQR(CosPhi);

      Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
      Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
      Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
      Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
      Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
      Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
      sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);

      parms=(REAL*)&Framework[CurrentSystem].BondTorsionArguments[f1][i];

      sprintf(string,"(F%1d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",f1,A,B,C,D,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][D].Type].Name);

      U=PrintBondTorsionEnergyStatus(nr++,string,Framework[CurrentSystem].BondTorsionType[f1][i],parms,rbc,Phi);

      // energy
      UHostBondTorsion+=U;
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework bond-torsion cross terms: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostBondTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostBondTorsion*ENERGY_TO_KELVIN,
          UHostBondTorsion*ENERGY_TO_KJ_PER_MOL,
          UHostBondTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintFrameworkBendTorsionEnergyStatus(void)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL U,rab,rbc,rcd;
  VECTOR Dab,Dbc,Dcd,dr,ds,Pb,Pc;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2;
  REAL CosTheta1,CosTheta2,Theta1,Theta2;
  REAL sign,Phi;
  REAL *parms,UHostBendTorsion;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostBendTorsion\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostBendTorsion=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBendTorsions[f1];i++)
    {
      A=Framework[CurrentSystem].BendTorsions[f1][i].A;
      B=Framework[CurrentSystem].BendTorsions[f1][i].B;
      C=Framework[CurrentSystem].BendTorsions[f1][i].C;
      D=Framework[CurrentSystem].BendTorsions[f1][i].D;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD=Framework[CurrentSystem].Atoms[f1][D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      Dab=ApplyBoundaryCondition(Dab);
      rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      Dbc=ApplyBoundaryCondition(Dbc);
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dcd.x=posD.x-posC.x;
      Dcd.y=posD.y-posC.y;
      Dcd.z=posD.z-posC.z;
      Dcd=ApplyBoundaryCondition(Dcd);
      rcd=sqrt(SQR(Dcd.x)+SQR(Dcd.y)+SQR(Dcd.z));

      dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      CosTheta1=dot_ab/rab;
      CosTheta1=MAX2(MIN2(CosTheta1,(REAL)1.0),-1.0);
      Theta1=acos(CosTheta1);

      dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;
      CosTheta2=-dot_cd/rcd;
      CosTheta2=MAX2(MIN2(CosTheta2,(REAL)1.0),-1.0);
      Theta2=acos(CosTheta2);

      dr.x=Dab.x-dot_ab*Dbc.x;
      dr.y=Dab.y-dot_ab*Dbc.y;
      dr.z=Dab.z-dot_ab*Dbc.z;
      r=MAX2((REAL)1.0e-8,sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z)));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Dcd.x-dot_cd*Dbc.x;
      ds.y=Dcd.y-dot_cd*Dbc.y;
      ds.z=Dcd.z-dot_cd*Dbc.z;
      s=MAX2((REAL)1.0e-8,sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z)));
      ds.x/=s; ds.y/=s; ds.z/=s;

      // compute Cos(Phi)
      // Phi is defined in protein convention Phi(trans)=Pi
      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      // Ensure CosPhi is between -1 and 1.
      CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
      CosPhi2=SQR(CosPhi);

      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);

      parms=Framework[CurrentSystem].BendTorsionArguments[f1][i];

      sprintf(string,"(F%1d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",f1,A,B,C,D,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][C].Type].Name,
                 PseudoAtoms[Framework[CurrentSystem].Atoms[f1][D].Type].Name);

      U=PrintBendTorsionEnergyStatus(nr++,string,Framework[CurrentSystem].BendTorsionType[f1][i],parms,Theta1,Theta2,Phi);

      // energy
      UHostBendTorsion+=U;
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of framework bend-torsion cross terms: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostBendTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostBendTorsion*ENERGY_TO_KELVIN,
          UHostBendTorsion*ENERGY_TO_KJ_PER_MOL,
          UHostBendTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


int PrintFrameworkIntraVDWEnergyStatus(void)
{
  int i,j,typeA,typeB,start;
  REAL rr,energy,UHostHostVDW;
  VECTOR posA,posB,dr;
  char string[256];
  int f1,f2;
  int nr,nr_excluded;

  // Framework-Framework energy
  UHostHostVDW=0.0;

  if(!InternalFrameworkLennardJonesInteractions) return 0;

  fprintf(OutputFilePtr[CurrentSystem],"UHostHostVDW\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  nr_excluded=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

        if(f1==f2) start=i+1;
        else start=0;
        for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
        {
          if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],0))
          {
            typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
            posB=Framework[CurrentSystem].Atoms[f2][j].AnisotropicPosition;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              energy=PotentialValue(typeA,typeB,rr,1.0);
              UHostHostVDW+=energy;

              sprintf(string,"(F%1d,%-4d)-(F%1d,%-4d) (%-6s,%-6s)",f1,i,f2,j,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);

              PrintVDWEnergyStatus(nr++,string,typeA,typeB,sqrt(rr),energy);
            }
          }
          else
            nr_excluded++;
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d (excluded: %d)\n",nr,nr_excluded);
  fprintf(OutputFilePtr[CurrentSystem],"Total Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostHostVDW*ENERGY_TO_KELVIN,
          UHostHostVDW*ENERGY_TO_KJ_PER_MOL,
          UHostHostVDW*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
  return 0;
}



void PrintFrameworkIntraChargeChargeEnergyStatus(void)
{
  int i,j,typeA,typeB,start;
  REAL chargeA,chargeB;
  REAL r,rr,energy;
  REAL UHostHostChargeChargeReal;
  VECTOR posA,posB,dr;
  char string[256];
  int f1,f2,nr,nr_excluded;

  fprintf(OutputFilePtr[CurrentSystem],"UHostHostChargeChargeReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  // Framework-Framework energy
  UHostHostChargeChargeReal=0.0;

  nr=0;
  nr_excluded=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        chargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;

        if(f1==f2) start=i+1;
        else start=0;
        for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
        {
          if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],1))
          {
            typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
            posB=Framework[CurrentSystem].Atoms[f2][j].Position;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              r=sqrt(rr);
              chargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;

              sprintf(string,"(F%1d,%-4d)-(F%1d,%-4d) (%-6s,%-6s)",f1,i,f2,j,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
              energy=PrintChargeChargeEnergyStatus(nr++,string,chargeA,chargeB,r);

              UHostHostChargeChargeReal+=energy;
            }
          }
          else
            nr_excluded++;
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d (excluded: %d)\n",nr,nr_excluded);
  fprintf(OutputFilePtr[CurrentSystem],"Total Energy: %18.10f [K] %18.10f [kJ/mol] %18.10f [kcal/mol]\n",
          UHostHostChargeChargeReal*ENERGY_TO_KELVIN,
          UHostHostChargeChargeReal*ENERGY_TO_KJ_PER_MOL,
          UHostHostChargeChargeReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}


void PrintFrameworkIntraChargeBondDipoleEnergyStatus(void)
{
  int i,j,f1,f2;
  int A1,A2;
  int typeA1,typeA2,typeB,nr;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL r,rr,ri2,cosA,energy,temp,length,chargeB;
  VECTOR dipoleA;
  REAL Bt1,UHostHostChargeBondDipoleReal;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;
  char string[256];
  int nr_excluded;

  Bt1=0.0;
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  fprintf(OutputFilePtr[CurrentSystem],"UHostHostChargeBondDipoleReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  nr_excluded=0;
  UHostHostChargeBondDipoleReal=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=0;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
      {
        DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
        A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
        A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
        typeA1=Framework[CurrentSystem].Atoms[f1][A1].Type;
        typeA2=Framework[CurrentSystem].Atoms[f1][A2].Type;
        posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
        posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        dipoleA=ApplyBoundaryCondition(dipoleA);
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
        length=sqrt(ri2);
        temp=DipoleMagnitudeA/length;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
        {
          // if framework are different, or if they are the same but not excluded within the framework
          if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][j][i],2))
          {
            typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
            if(PseudoAtoms[typeB].HasCharges)
            {
              posB=Framework[CurrentSystem].Atoms[f2][j].Position;
              chargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeBondDipoleSquared)
              {
                r=sqrt(rr);
                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

                sprintf(string,"(F%1d,[%-4d,%-4d])-(F%1d,%-4d) ([%-6s,%-6s],%-6s)",f1,A1,A2,f2,j,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,PseudoAtoms[typeB].Name);
                energy=PrintChargeBondDipoleStatus(nr++,string,cosA,DipoleMagnitudeA,chargeB,r);

                UHostHostChargeBondDipoleReal+=energy;
              }
            }
          }
          else
            nr_excluded++;
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d (excluded: %d)\n",nr,nr_excluded);
  fprintf(OutputFilePtr[CurrentSystem],"Total Energy: %18.10f [K] %18.10f [kJ/mol] %18.10f [kcal/mol]\n",
          UHostHostChargeBondDipoleReal*ENERGY_TO_KELVIN,
          UHostHostChargeBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UHostHostChargeBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}

int PrintFrameworkIntraBondDipoleBondDipoleEnergyStatus(void)
{
  int i,j,f1,f2;
  int A1,A2,B1,B2;
  int start;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,rr,ri2,rk2,cosAB,cosA,cosB,energy,temp,length;
  VECTOR dipoleA,dipoleB;
  int typeA1,typeA2,typeB1,typeB2;
  REAL Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL UHostHostBondDipoleBondDipoleReal;
  char string[256];
  REAL_MATRIX3x3 v;
  int nr_excluded,nr;

  Bt1=Bt2=0.0;
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  UHostHostBondDipoleBondDipoleReal=0.0;
  if(ChargeMethod==NONE) return 0;

  fprintf(OutputFilePtr[CurrentSystem],"UHostHostBondDipoleBondDipoleReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  nr_excluded=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
      {
        DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
        A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
        A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
        posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
        posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
        typeA1=Framework[CurrentSystem].Atoms[f1][A1].Type;
        typeA2=Framework[CurrentSystem].Atoms[f1][A2].Type;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        dipoleA=ApplyBoundaryCondition(dipoleA);
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
        length=sqrt(ri2);
        temp=DipoleMagnitudeA/length;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        if(f1==f2) start=i+1;
        else start=0;
        for(j=start;j<Framework[CurrentSystem].NumberOfBondDipoles[f2];j++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f2][j];
          B1=Framework[CurrentSystem].BondDipoles[f2][j].A;
          B2=Framework[CurrentSystem].BondDipoles[f2][j].B;

          // if framework are different, or if they are the same but not excluded within the framework
          if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f2][i][j],3))
          {
            posB1=Framework[CurrentSystem].Atoms[f2][B1].Position;
            posB2=Framework[CurrentSystem].Atoms[f2][B2].Position;
            typeB1=Framework[CurrentSystem].Atoms[f2][B1].Type;
            typeB2=Framework[CurrentSystem].Atoms[f2][B2].Type;

            dipoleB.x=posB2.x-posB1.x;
            dipoleB.y=posB2.y-posB1.y;
            dipoleB.z=posB2.z-posB1.z;
            dipoleB=ApplyBoundaryCondition(dipoleB);
            posB.x=posB1.x+0.5*dipoleB.x;
            posB.y=posB1.y+0.5*dipoleB.y;
            posB.z=posB1.z+0.5*dipoleB.z;
            rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
            length=sqrt(rk2);
            temp=DipoleMagnitudeB/length;
            dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              r=sqrt(rr);
              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

              sprintf(string,"(F%1d,[%-4d,%-4d])-(F%1d,[%-4d,-%4d]) ([%-6s,%-6s],[%-6s,%-6s])",f1,A1,A2,f2,B1,B2,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,
                               PseudoAtoms[typeB1].Name,PseudoAtoms[typeB2].Name);
              energy=PrintBondDipoleBondDipoleStatus(nr++,string,cosA,cosB,cosAB,DipoleMagnitudeA,DipoleMagnitudeB,r);

              UHostHostBondDipoleBondDipoleReal+=energy;
            }
          }
          else
            nr_excluded++;
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d (excluded: %d)\n",nr,nr_excluded);
  fprintf(OutputFilePtr[CurrentSystem],"Total Energy: %18.10f [K] %18.10f [kJ/mol] %18.10f [kcal/mol]\n",
          UHostHostBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,
          UHostHostBondDipoleBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UHostHostBondDipoleBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

  return 0;
}

// =============================================================================================================================
// adsorbate intra-molecular energies
// =============================================================================================================================


void PrintAdsorbateBondEnergyStatus(void)
{
  int m,i,nr,Type,A,B;
  REAL U,rr,UAdsorbateBond;
  REAL *parms;
  char string[256];
  POINT posA,posB;
  VECTOR dr;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateBond:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateBond=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBonds;i++)
    {
      A=Components[Type].Bonds[i].A;
      B=Components[Type].Bonds[i].B;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      parms=(REAL*)&Components[Type].BondArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d]) (%-6s,%-6s)",m,A,B,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name);
      U=PrintBondEnergyStatus(nr++,string,Components[Type].BondType[i],parms,sqrt(rr));

      UAdsorbateBond+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate bonds: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateBond Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateBond*ENERGY_TO_KELVIN,
          UAdsorbateBond*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateBond*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintAdsorbateUreyBradleyEnergyStatus(void)
{
  int m,i,nr,Type,A,C;
  REAL U,rr,UAdsorbateUreyBradley;
  REAL *parms;
  char string[256];
  POINT posA,posC;
  VECTOR dr;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateUreyBradley:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateUreyBradley=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBonds;i++)
    {
      A=Components[Type].UreyBradleys[i].A;
      C=Components[Type].UreyBradleys[i].C;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

      dr.x=posA.x-posC.x;
      dr.y=posA.y-posC.y;
      dr.z=posA.z-posC.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      parms=Components[Type].UreyBradleyArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d]) (%-6s,%-6s)",m,A,C,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name);
      U=PrintUreyBradleyEnergyStatus(nr++,string,Components[Type].UreyBradleyType[i],parms,sqrt(rr));

      UAdsorbateUreyBradley+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate UreyBradleys: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateUreyBradley Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateUreyBradley*ENERGY_TO_KELVIN,
          UAdsorbateUreyBradley*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateUreyBradley*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintAdsorbateBendEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UAdsorbateBend;
  char string[256];
  REAL *parms,U;
  REAL CosTheta,Theta;
  REAL rab,rbc,rac;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL delta,rt2,rap2,rcp2;


  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateBend:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateBend=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBends;i++)
    {
      A=Components[Type].Bends[i].A;
      B=Components[Type].Bends[i].B;
      C=Components[Type].Bends[i].C;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

      parms=Components[Type].BendArguments[i];

      switch(Components[Type].BendType[i])
      {
        case MM3_IN_PLANE_BEND:
          D=Components[Type].Bends[i].D;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;
          Rad.x=posA.x-posD.x;
          Rad.y=posA.y-posD.y;
          Rad.z=posA.z-posD.z;

          Rbd.x=posB.x-posD.x;
          Rbd.y=posB.y-posD.y;
          Rbd.z=posB.z-posD.z;

          Rcd.x=posC.x-posD.x;
          Rcd.y=posC.y-posD.y;
          Rcd.z=posC.z-posD.z;

          t.x=Rad.y*Rcd.z-Rad.z*Rcd.y;
          t.y=Rad.z*Rcd.x-Rad.x*Rcd.z;
          t.z=Rad.x*Rcd.y-Rad.y*Rcd.x;
          rt2=t.x*t.x+t.y*t.y+t.z*t.z;
          delta=-(t.x*Rbd.x+t.y*Rbd.y+t.z*Rbd.z)/rt2;

          ip.x=posB.x+t.x*delta;
          ip.y=posB.y+t.y*delta;
          ip.z=posB.z+t.z*delta;
          ap.x=posA.x-ip.x;
          ap.y=posA.y-ip.y;
          ap.z=posA.z-ip.z;
          cp.x=posC.x-ip.x;
          cp.y=posC.y-ip.y;
          cp.z=posC.z-ip.z;

          rap2=ap.x*ap.x+ap.y*ap.y+ap.z*ap.z;
          rcp2=cp.x*cp.x+cp.y*cp.y+cp.z*cp.z;

          CosTheta=(ap.x*cp.x+ap.y*cp.y+ap.z*cp.z)/sqrt(rap2*rcp2);
          break;
        default:
          Rab.x=posA.x-posB.x;
          Rab.y=posA.y-posB.y;
          Rab.z=posA.z-posB.z;
          rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
          Rab.x/=rab;
          Rab.y/=rab;
          Rab.z/=rab;

          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;
          rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
          Rbc.x/=rbc;
          Rbc.y/=rbc;
          Rbc.z/=rbc;

          Rac.x=posC.x-posA.x;
          Rac.y=posC.y-posA.y;
          Rac.z=posC.z-posA.z;
          rac=sqrt(SQR(Rac.x)+SQR(Rac.y)+SQR(Rac.z));
          Rac.x/=rac;
          Rac.y/=rac;
          Rac.z/=rac;

          CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
          break;
      }

      CosTheta=MIN2(1.0,MAX2(-1.0,CosTheta));
      Theta=acos(CosTheta);

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s)",m,A,B,C,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name);

      U=PrintBendEnergyStatus(nr++,string,Components[Type].BendType[i],parms,Theta);

      UAdsorbateBend+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate bends: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateBend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateBend*ENERGY_TO_KELVIN,
          UAdsorbateBend*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintAdsorbateInversionBendEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UAdsorbateInversionBend;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC,posD;
  REAL rrab,rab2;
  REAL CosChi,Chi;
  VECTOR Rab,Rbc,Rbd,Rcd,Rad;
  REAL c,e;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateInversionBend:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateInversionBend=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfInversionBends;i++)
    {
      A=Components[Type].InversionBends[i].A;
      B=Components[Type].InversionBends[i].B;
      C=Components[Type].InversionBends[i].C;
      D=Components[Type].InversionBends[i].D;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
      posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      rab2=Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z;
      rrab=sqrt(rab2);

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;

      Rbd.x=posD.x-posB.x;
      Rbd.y=posD.y-posB.y;
      Rbd.z=posD.z-posB.z;

      Rcd.x=posD.x-posC.x;
      Rcd.y=posD.y-posC.y;
      Rcd.z=posD.z-posC.z;

      Rad.x=posD.x-posA.x;
      Rad.y=posD.y-posA.y;
      Rad.z=posD.z-posA.z;

      switch(Components[Type].InversionBendType[i])
      {
        case HARMONIC_INVERSION:
        case HARMONIC_COSINE_INVERSION:
        case PLANAR_INVERSION:
          // w is a vector perpendicular to the B-C-D plane
          // c=w.w=(Rbc x Rbd).(Rbc x Rbd)= r_bc^2 r_bd^2 - (r_cb . r_bd)^2
          c=(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z)*(Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z)-SQR(Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z);
          break;
        case HARMONIC_INVERSION2:
        case HARMONIC_COSINE_INVERSION2:
        case PLANAR_INVERSION2:
        case MM3_INVERSION:
          // w is a vector perpendicular to the A-C-D plane
          // c=w.w=(Rcd x Rad).(Rcd x Rad)=r_cd^2 r_ad^2 - (r_da . r_cd)^2
          c=(Rcd.x*Rcd.x+Rcd.y*Rcd.y+Rcd.z*Rcd.z)*(Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z)-SQR(Rad.x*Rcd.x+Rad.y*Rcd.y+Rad.z*Rcd.z);
          break;
        default:
          fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateInversionBendEnergyAdsorbate' ('internal_energy.c')\n");
          exit(0);
          break;
      }

      e=Rab.x*(Rbd.y*Rbc.z-Rbd.z*Rbc.y)+Rab.y*(Rbd.z*Rbc.x-Rbd.x*Rbc.z)+Rab.z*(Rbd.x*Rbc.y-Rbd.y*Rbc.x);
      CosChi=sqrt((Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z)-SQR(e)/c)/rrab;

      // Ensure CosChi is between -1 and 1.
      CosChi=SIGN(MIN2(fabs(CosChi),(REAL)1.0),CosChi);
      Chi=acos(CosChi);

      parms=Components[Type].InversionBendArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintInversionBendEnergyStatus(nr++,string,Components[Type].InversionBendType[i],parms,Chi);

      UAdsorbateInversionBend+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate inversion-bends: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateBInversionend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateInversionBend*ENERGY_TO_KELVIN,
          UAdsorbateInversionBend*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateInversionBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintAdsorbateTorsionEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UAdsorbateTorsion;
  char string[256];
  REAL *parms,U;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,rbc;
  VECTOR Pb,Pc;
  VECTOR posA,posB,posC,posD;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateTorsion:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateTorsion=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfTorsions;i++)
    {
      A=Components[Type].Torsions[i].A;
      B=Components[Type].Torsions[i].B;
      C=Components[Type].Torsions[i].C;
      D=Components[Type].Torsions[i].D;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
      posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dcd.x=posD.x-posC.x;
      Dcd.y=posD.y-posC.y;
      Dcd.z=posD.z-posC.z;

      dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;

      dr.x=Dab.x-dot_ab*Dbc.x;
      dr.y=Dab.y-dot_ab*Dbc.y;
      dr.z=Dab.z-dot_ab*Dbc.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Dcd.x-dot_cd*Dbc.x;
      ds.y=Dcd.y-dot_cd*Dbc.y;
      ds.z=Dcd.z-dot_cd*Dbc.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      // compute Cos(Phi)
      // Phi is defined in protein convention Phi(trans)=Pi
      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);

      parms=Components[Type].TorsionArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintTorsionEnergyStatus(nr++,string,Components[Type].TorsionType[i],parms,Phi);

      UAdsorbateTorsion+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate torsions: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateBTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateTorsion*ENERGY_TO_KELVIN,
          UAdsorbateTorsion*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}



void PrintAdsorbateImproperTorsionEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UAdsorbateImproperTorsion;
  char string[256];
  REAL *parms,U;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,rbc;
  VECTOR Pb,Pc;
  VECTOR posA,posB,posC,posD;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateImproperTorsion:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateImproperTorsion=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfImproperTorsions;i++)
    {
      A=Components[Type].ImproperTorsions[i].A;
      B=Components[Type].ImproperTorsions[i].B;
      C=Components[Type].ImproperTorsions[i].C;
      D=Components[Type].ImproperTorsions[i].D;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
      posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dcd.x=posD.x-posC.x;
      Dcd.y=posD.y-posC.y;
      Dcd.z=posD.z-posC.z;

      dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;

      dr.x=Dab.x-dot_ab*Dbc.x;
      dr.y=Dab.y-dot_ab*Dbc.y;
      dr.z=Dab.z-dot_ab*Dbc.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Dcd.x-dot_cd*Dbc.x;
      ds.y=Dcd.y-dot_cd*Dbc.y;
      ds.z=Dcd.z-dot_cd*Dbc.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      // compute Cos(Phi)
      // Phi is defined in protein convention Phi(trans)=Pi
      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);

      parms=Components[Type].ImproperTorsionArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintImproperTorsionEnergyStatus(nr++,string,Components[Type].ImproperTorsionType[i],parms,Phi);

      UAdsorbateImproperTorsion+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate improper torsions: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateBImproperTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateImproperTorsion*ENERGY_TO_KELVIN,
          UAdsorbateImproperTorsion*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateImproperTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintAdsorbateBondBondEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C;
  REAL UAdsorbateBondBond;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc;
  REAL rab,rbc;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateBondBond:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateBondBond=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBondBonds;i++)
    {
      A=Components[Type].BondBonds[i].A;
      B=Components[Type].BondBonds[i].B;
      C=Components[Type].BondBonds[i].C;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));

      parms=Components[Type].BondBondArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s)",m,A,B,C,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name);

      U=PrintBondBondEnergyStatus(nr++,string,Components[Type].BondBondType[i],parms,rab,rbc);

      UAdsorbateBondBond+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate bond-bonds: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateBondBond Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateBondBond*ENERGY_TO_KELVIN,
          UAdsorbateBondBond*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateBondBond*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintAdsorbateBondBendEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C;
  REAL UAdsorbateBondBend;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc;
  REAL rab,rbc,CosTheta,Theta;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateBondBend:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateBondBend=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBondBends;i++)
    {
      A=Components[Type].BondBends[i].A;
      B=Components[Type].BondBends[i].B;
      C=Components[Type].BondBends[i].C;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
      Rab.x/=rab;
      Rab.y/=rab;
      Rab.z/=rab;

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
      Rbc.x/=rbc;
      Rbc.y/=rbc;
      Rbc.z/=rbc;

      CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
      Theta=acos(CosTheta);

      parms=Components[Type].BondBendArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s)",m,A,B,C,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name);

      U=PrintBondBendEnergyStatus(nr++,string,Components[Type].BondBendType[i],parms,rab,rbc,Theta);

      UAdsorbateBondBend+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate bond-bends: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateBondBend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateBondBend*ENERGY_TO_KELVIN,
          UAdsorbateBondBend*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateBondBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintAdsorbateBendBendEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UAdsorbateBendBend;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC,posD;
  VECTOR Dab,Dbc,Dbd;
  REAL rab,rbc,rbd;
  REAL dot_abc,dot_abd;
  REAL Theta1,Theta2;
  REAL CosTheta1,CosTheta2;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateBendBend:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateBendBend=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBendBends;i++)
    {
      A=Components[Type].BendBends[i].A;
      B=Components[Type].BendBends[i].B;
      C=Components[Type].BendBends[i].C;
      D=Components[Type].BendBends[i].D;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
      posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));
      Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dbd.x=posD.x-posB.x;
      Dbd.y=posD.y-posB.y;
      Dbd.z=posD.z-posB.z;
      rbd=sqrt(SQR(Dbd.x)+SQR(Dbd.y)+SQR(Dbd.z));
      Dbd.x/=rbd; Dbd.y/=rbd; Dbd.z/=rbd;

      dot_abc=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      CosTheta1=dot_abc;
      CosTheta1=SIGN(MIN2(fabs(CosTheta1),(REAL)1.0),CosTheta1);
      Theta1=acos(CosTheta1);

      dot_abd=Dab.x*Dbd.x+Dab.y*Dbd.y+Dab.z*Dbd.z;
      CosTheta2=dot_abd;
      CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
      Theta2=acos(CosTheta2);

      parms=Components[Type].BendBendArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintBendBendEnergyStatus(nr++,string,Components[Type].BendBendType[i],parms,Theta1,Theta2);

      UAdsorbateBendBend+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate bend-bends: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateBendBend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateBendBend*ENERGY_TO_KELVIN,
          UAdsorbateBendBend*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateBendBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintAdsorbateBondTorsionEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UAdsorbateBondTorsion;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC,posD;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL rbc,s,r,sign;
  REAL dot_ab,dot_cd,CosPhi,Phi;
  VECTOR Pb,Pc;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateBondTorsion:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateBondTorsion=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBondTorsions;i++)
    {
      A=Components[Type].BondTorsions[i].A;
      B=Components[Type].BondTorsions[i].B;
      C=Components[Type].BondTorsions[i].C;
      D=Components[Type].BondTorsions[i].D;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
      posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dcd.x=posD.x-posC.x;
      Dcd.y=posD.y-posC.y;
      Dcd.z=posD.z-posC.z;

      dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;

      dr.x=Dab.x-dot_ab*Dbc.x;
      dr.y=Dab.y-dot_ab*Dbc.y;
      dr.z=Dab.z-dot_ab*Dbc.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Dcd.x-dot_cd*Dbc.x;
      ds.y=Dcd.y-dot_cd*Dbc.y;
      ds.z=Dcd.z-dot_cd*Dbc.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      // Ensure CosPhi is between -1 and 1.
      CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);

      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);


      parms=Components[Type].BondTorsionArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintBondTorsionEnergyStatus(nr++,string,Components[Type].BondTorsionType[i],parms,rbc,Phi);

      UAdsorbateBondTorsion+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate bond-torsions: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateBondTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateBondTorsion*ENERGY_TO_KELVIN,
          UAdsorbateBondTorsion*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateBondTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintAdsorbateBendTorsionEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UAdsorbateBendTorsion;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC,posD;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL rab,rbc,rcd,s,r,sign;
  REAL dot_ab,dot_cd,CosPhi,Phi;
  REAL CosTheta1,CosTheta2;
  REAL Theta1,Theta2;
  VECTOR Pb,Pc;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateBendTorsion:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateBendTorsion=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBendTorsions;i++)
    {
      A=Components[Type].BendTorsions[i].A;
      B=Components[Type].BendTorsions[i].B;
      C=Components[Type].BendTorsions[i].C;
      D=Components[Type].BendTorsions[i].D;

      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
      posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dcd.x=posD.x-posC.x;
      Dcd.y=posD.y-posC.y;
      Dcd.z=posD.z-posC.z;
      rcd=sqrt(SQR(Dcd.x)+SQR(Dcd.y)+SQR(Dcd.z));

      dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      CosTheta1=dot_ab/rab;
      CosTheta1=SIGN(MIN2(fabs(CosTheta1),(REAL)1.0),CosTheta1);
      Theta1=acos(CosTheta1);

      dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;
      CosTheta2=-dot_cd/rcd;
      CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
      Theta2=acos(CosTheta2);

      dr.x=Dab.x-dot_ab*Dbc.x;
      dr.y=Dab.y-dot_ab*Dbc.y;
      dr.z=Dab.z-dot_ab*Dbc.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Dcd.x-dot_cd*Dbc.x;
      ds.y=Dcd.y-dot_cd*Dbc.y;
      ds.z=Dcd.z-dot_cd*Dbc.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      // Ensure CosPhi is between -1 and 1.
      CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);

      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);

      parms=Components[Type].BendTorsionArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintBendTorsionEnergyStatus(nr++,string,Components[Type].BendTorsionType[i],parms,Theta1,Theta2,Phi);

      UAdsorbateBendTorsion+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate bend-torsions: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateBendTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateBendTorsion*ENERGY_TO_KELVIN,
          UAdsorbateBendTorsion*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateBendTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintAdsorbateIntraVDWEnergyStatus(void)
{
  int i,NumberOfIntraLJ;
  int Type,TypeA,TypeB,A,B;
  int nr,m;
  REAL rr,Scaling;
  VECTOR dr;
  POINT posA,posB;
  REAL UAdsorbateIntraVDW,energy;
  char string[256];

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateIntraVDW:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateIntraVDW=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfIntraLJ=Components[Type].NumberOfIntraVDW;
    for(i=0;i<NumberOfIntraLJ;i++)
    {
      A=Components[Type].IntraVDW[i].A;
      B=Components[Type].IntraVDW[i].B;
      Scaling=Components[Type].IntraVDWScaling[i];
      TypeA=Adsorbates[CurrentSystem][m].Atoms[A].Type;
      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      if(rr<CutOffVDWSquared)
      {
        TypeB=Adsorbates[CurrentSystem][m].Atoms[B].Type;
        energy=Scaling*PotentialValue(TypeA,TypeB,rr,1.0);

        sprintf(string,"(Scaling: %8.5f) (A%-4d,[%-4d,%-4d]) (%-6s,%-6s)",Scaling,m,A,B,PseudoAtoms[TypeA].Name,PseudoAtoms[TypeB].Name);
        PrintVDWEnergyStatus(nr,string,TypeA,TypeB,sqrt(rr),energy);

        UAdsorbateIntraVDW+=energy;
      }
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate intra-VDW: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateIntraVDW Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateIntraVDW*ENERGY_TO_KELVIN,
          UAdsorbateIntraVDW*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateIntraVDW*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}

void PrintAdsorbateIntraChargeChargeEnergyStatus(void)
{
  int m,i,NumberOfIntraChargeChargeCoulomb,Type,TypeA,TypeB,A,B;
  REAL r,rr,chargeA,chargeB,Scaling;
  REAL UAdsorbateIntraChargeCharge,energy;
  VECTOR dr;
  POINT posA,posB;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateIntraChargeCharge:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateIntraChargeCharge=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfIntraChargeChargeCoulomb=Components[Type].NumberOfIntraChargeCharge;
    for(i=0;i<NumberOfIntraChargeChargeCoulomb;i++)
    {
      A=Components[Type].IntraChargeCharge[i].A;
      B=Components[Type].IntraChargeCharge[i].B;
      Scaling=Components[Type].IntraChargeChargeScaling[i];
      TypeA=Adsorbates[CurrentSystem][m].Atoms[A].Type;
      chargeA=Adsorbates[CurrentSystem][m].Atoms[A].Charge;
      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      TypeB=Adsorbates[CurrentSystem][m].Atoms[B].Type;
      chargeB=Adsorbates[CurrentSystem][m].Atoms[B].Charge;
      energy=Scaling*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;

      sprintf(string,"(A%-4d,[%-4d,%-4d)) (%-6s,%-6s)",m,A,B,PseudoAtoms[TypeA].Name,PseudoAtoms[TypeB].Name);
      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_COULOMBIC %s, Charge: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,chargeA,chargeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);

      UAdsorbateIntraChargeCharge+=energy;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate intra-charge-charge: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateIntraChargeCharge Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateIntraChargeCharge*ENERGY_TO_KELVIN,
          UAdsorbateIntraChargeCharge*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateIntraChargeCharge*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}


void PrintAdsorbateIntraChargeBondDipoleEnergyStatus(void)
{
  int m,i,NumberOfIntraChargeBondDipoleCoulomb,Type,A,B,B1,B2;
  REAL r,rr,ri2,length;
  REAL Bt1,cosB;
  VECTOR dr,dipoleB;
  POINT posA,posB,posB1,posB2;
  REAL ChargeA,temp,UAdsorbateIntraChargeBondDipole;
  REAL DipoleMagnitudeB,energy;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateIntraChargeBondDipole:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateIntraChargeBondDipole=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfIntraChargeBondDipoleCoulomb=Components[Type].NumberOfIntraChargeBondDipole;
    for(i=0;i<NumberOfIntraChargeBondDipoleCoulomb;i++)
    {
      A=Components[Type].IntraChargeBondDipole[i].A;
      B=Components[Type].IntraChargeBondDipole[i].B;

      ChargeA=Adsorbates[CurrentSystem][m].Atoms[A].Charge;
      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;

      B1=Components[Type].BondDipoles[B].A;
      B2=Components[Type].BondDipoles[B].B;
      DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
      posB1=Adsorbates[CurrentSystem][m].Atoms[B1].Position;
      posB2=Adsorbates[CurrentSystem][m].Atoms[B2].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      ri2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeB/length;
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt1=1.0/(r*rr);
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy=-COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);

      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_COULOMB %s, Magnitude: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,DipoleMagnitudeB,ChargeA,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);

      UAdsorbateIntraChargeBondDipole+=energy;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate intra-charge-bonddipole: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateIntraChargeBondDipole Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateIntraChargeBondDipole*ENERGY_TO_KELVIN,
          UAdsorbateIntraChargeBondDipole*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateIntraChargeBondDipole*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintAdsorbateIntraBondDipoleBondDipoleEnergyStatus(void)
{
  int i,m,NumberOfIntraBondDipoleBondDipoleCoulomb,Type,A,B,A1,A2,B1,B2;
  REAL r,rr,ri2,rk2,length;
  REAL Bt0,Bt1,Bt2,Bt3,cosA,cosB,cosAB;
  VECTOR dr,dipoleA,dipoleB;
  POINT posA,posB,posA1,posA2,posB1,posB2;
  REAL temp,UAdsorbateIntraBondDipoleBondDipole;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,energy;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UAdsorbateIntraBondDipoleBondDipole:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;


  UAdsorbateIntraBondDipoleBondDipole=0.0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfIntraBondDipoleBondDipoleCoulomb=Components[Type].NumberOfIntraBondDipoleBondDipole;
    for(i=0;i<NumberOfIntraBondDipoleBondDipoleCoulomb;i++)
    {
      A=Components[Type].IntraBondDipoleBondDipole[i].A;
      B=Components[Type].IntraBondDipoleBondDipole[i].B;

      A1=Components[Type].BondDipoles[A].A;
      A2=Components[Type].BondDipoles[A].B;
      posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
      DipoleMagnitudeA=Components[Type].BondDipoleMagnitude[A];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      B1=Components[Type].BondDipoles[B].A;
      B2=Components[Type].BondDipoles[B].B;
      posB1=Adsorbates[CurrentSystem][m].Atoms[B1].Position;
      posB2=Adsorbates[CurrentSystem][m].Atoms[B2].Position;
      DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
      length=sqrt(rk2);
      temp=DipoleMagnitudeB/length;
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=1.0/(r);
      Bt1=1.0/(r*rr);
      Bt2=3.0/(r*rr*rr);
      Bt3=15.0/(r*rr*rr*rr);

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);

      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_COULOMB %s, Magnitude: %8.5f, Magnitude: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,DipoleMagnitudeA,DipoleMagnitudeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);

      UAdsorbateIntraBondDipoleBondDipole+=energy;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of adsorbate intra-bonddipole-bonddipole: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateIntraBondDipoleBondDipole Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateIntraBondDipoleBondDipole*ENERGY_TO_KELVIN,
          UAdsorbateIntraBondDipoleBondDipole*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateIntraBondDipoleBondDipole*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


// =============================================================================================================================
// cation intra-molecular energies
// =============================================================================================================================


void PrintCationBondEnergyStatus(void)
{
  int m,i,nr,Type,A,B;
  REAL U,rr,UCationBond;
  REAL *parms;
  char string[256];
  POINT posA,posB;
  VECTOR dr;

  fprintf(OutputFilePtr[CurrentSystem],"UCationBond:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationBond=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBonds;i++)
    {
      A=Components[Type].Bonds[i].A;
      B=Components[Type].Bonds[i].B;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      parms=(REAL*)&Components[Type].BondArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d]) (%-6s,%-6s)",m,A,B,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name);
      U=PrintBondEnergyStatus(nr++,string,Components[Type].BondType[i],parms,sqrt(rr));

      UCationBond+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation bonds: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationBond Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationBond*ENERGY_TO_KELVIN,
          UCationBond*ENERGY_TO_KJ_PER_MOL,
          UCationBond*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintCationUreyBradleyEnergyStatus(void)
{
  int m,i,nr,Type,A,C;
  REAL U,rr,UCationUreyBradley;
  REAL *parms;
  char string[256];
  POINT posA,posC;
  VECTOR dr;

  fprintf(OutputFilePtr[CurrentSystem],"UCationUreyBradley:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationUreyBradley=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBonds;i++)
    {
      A=Components[Type].UreyBradleys[i].A;
      C=Components[Type].UreyBradleys[i].C;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posC=Cations[CurrentSystem][m].Atoms[C].Position;

      dr.x=posA.x-posC.x;
      dr.y=posA.y-posC.y;
      dr.z=posA.z-posC.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      parms=Components[Type].UreyBradleyArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d]) (%-6s,%-6s)",m,A,C,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name);
      U=PrintUreyBradleyEnergyStatus(nr++,string,Components[Type].UreyBradleyType[i],parms,sqrt(rr));

      UCationUreyBradley+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation UreyBradleys: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationUreyBradley Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationUreyBradley*ENERGY_TO_KELVIN,
          UCationUreyBradley*ENERGY_TO_KJ_PER_MOL,
          UCationUreyBradley*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintCationBendEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UCationBend;
  char string[256];
  REAL *parms,U;
  REAL CosTheta,Theta;
  REAL rab,rbc,rac;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL delta,rt2,rap2,rcp2;


  fprintf(OutputFilePtr[CurrentSystem],"UCationBend:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationBend=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBends;i++)
    {
      A=Components[Type].Bends[i].A;
      B=Components[Type].Bends[i].B;
      C=Components[Type].Bends[i].C;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=Cations[CurrentSystem][m].Atoms[C].Position;

      parms=Components[Type].BendArguments[i];

      switch(Components[Type].BendType[i])
      {
        case MM3_IN_PLANE_BEND:
          D=Components[Type].Bends[i].D;
          posD=Cations[CurrentSystem][m].Atoms[D].Position;
          Rad.x=posA.x-posD.x;
          Rad.y=posA.y-posD.y;
          Rad.z=posA.z-posD.z;

          Rbd.x=posB.x-posD.x;
          Rbd.y=posB.y-posD.y;
          Rbd.z=posB.z-posD.z;

          Rcd.x=posC.x-posD.x;
          Rcd.y=posC.y-posD.y;
          Rcd.z=posC.z-posD.z;

          t.x=Rad.y*Rcd.z-Rad.z*Rcd.y;
          t.y=Rad.z*Rcd.x-Rad.x*Rcd.z;
          t.z=Rad.x*Rcd.y-Rad.y*Rcd.x;
          rt2=t.x*t.x+t.y*t.y+t.z*t.z;
          delta=-(t.x*Rbd.x+t.y*Rbd.y+t.z*Rbd.z)/rt2;

          ip.x=posB.x+t.x*delta;
          ip.y=posB.y+t.y*delta;
          ip.z=posB.z+t.z*delta;
          ap.x=posA.x-ip.x;
          ap.y=posA.y-ip.y;
          ap.z=posA.z-ip.z;
          cp.x=posC.x-ip.x;
          cp.y=posC.y-ip.y;
          cp.z=posC.z-ip.z;

          rap2=ap.x*ap.x+ap.y*ap.y+ap.z*ap.z;
          rcp2=cp.x*cp.x+cp.y*cp.y+cp.z*cp.z;

          CosTheta=(ap.x*cp.x+ap.y*cp.y+ap.z*cp.z)/sqrt(rap2*rcp2);
          break;
        default:
          Rab.x=posA.x-posB.x;
          Rab.y=posA.y-posB.y;
          Rab.z=posA.z-posB.z;
          rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
          Rab.x/=rab;
          Rab.y/=rab;
          Rab.z/=rab;

          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;
          rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
          Rbc.x/=rbc;
          Rbc.y/=rbc;
          Rbc.z/=rbc;

          Rac.x=posC.x-posA.x;
          Rac.y=posC.y-posA.y;
          Rac.z=posC.z-posA.z;
          rac=sqrt(SQR(Rac.x)+SQR(Rac.y)+SQR(Rac.z));
          Rac.x/=rac;
          Rac.y/=rac;
          Rac.z/=rac;

          CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
          break;
      }

      CosTheta=MIN2(1.0,MAX2(-1.0,CosTheta));
      Theta=acos(CosTheta);

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s)",m,A,B,C,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name);

      U=PrintBendEnergyStatus(nr++,string,Components[Type].BendType[i],parms,Theta);

      UCationBend+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation bends: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationBend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationBend*ENERGY_TO_KELVIN,
          UCationBend*ENERGY_TO_KJ_PER_MOL,
          UCationBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintCationInversionBendEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UCationInversionBend;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC,posD;
  REAL rrab,rab2;
  REAL CosChi,Chi;
  VECTOR Rab,Rbc,Rbd,Rcd,Rad;
  REAL c,e;

  fprintf(OutputFilePtr[CurrentSystem],"UCationInversionBend:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationInversionBend=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfInversionBends;i++)
    {
      A=Components[Type].InversionBends[i].A;
      B=Components[Type].InversionBends[i].B;
      C=Components[Type].InversionBends[i].C;
      D=Components[Type].InversionBends[i].D;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=Cations[CurrentSystem][m].Atoms[C].Position;
      posD=Cations[CurrentSystem][m].Atoms[D].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      rab2=Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z;
      rrab=sqrt(rab2);

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;

      Rbd.x=posD.x-posB.x;
      Rbd.y=posD.y-posB.y;
      Rbd.z=posD.z-posB.z;

      Rcd.x=posD.x-posC.x;
      Rcd.y=posD.y-posC.y;
      Rcd.z=posD.z-posC.z;

      Rad.x=posD.x-posA.x;
      Rad.y=posD.y-posA.y;
      Rad.z=posD.z-posA.z;

      switch(Components[Type].InversionBendType[i])
      {
        case HARMONIC_INVERSION:
        case HARMONIC_COSINE_INVERSION:
        case PLANAR_INVERSION:
          // w is a vector perpendicular to the B-C-D plane
          // c=w.w=(Rbc x Rbd).(Rbc x Rbd)= r_bc^2 r_bd^2 - (r_cb . r_bd)^2
          c=(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z)*(Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z)-SQR(Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z);
          break;
        case HARMONIC_INVERSION2:
        case HARMONIC_COSINE_INVERSION2:
        case PLANAR_INVERSION2:
        case MM3_INVERSION:
          // w is a vector perpendicular to the A-C-D plane
          // c=w.w=(Rcd x Rad).(Rcd x Rad)=r_cd^2 r_ad^2 - (r_da . r_cd)^2
          c=(Rcd.x*Rcd.x+Rcd.y*Rcd.y+Rcd.z*Rcd.z)*(Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z)-SQR(Rad.x*Rcd.x+Rad.y*Rcd.y+Rad.z*Rcd.z);
          break;
        default:
          fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateInversionBendEnergyCation' ('internal_energy.c')\n");
          exit(0);
          break;
      }

      e=Rab.x*(Rbd.y*Rbc.z-Rbd.z*Rbc.y)+Rab.y*(Rbd.z*Rbc.x-Rbd.x*Rbc.z)+Rab.z*(Rbd.x*Rbc.y-Rbd.y*Rbc.x);
      CosChi=sqrt((Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z)-SQR(e)/c)/rrab;

      // Ensure CosChi is between -1 and 1.
      CosChi=SIGN(MIN2(fabs(CosChi),(REAL)1.0),CosChi);
      Chi=acos(CosChi);

      parms=Components[Type].InversionBendArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintInversionBendEnergyStatus(nr++,string,Components[Type].InversionBendType[i],parms,Chi);

      UCationInversionBend+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation inversion-bends: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationBInversionend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationInversionBend*ENERGY_TO_KELVIN,
          UCationInversionBend*ENERGY_TO_KJ_PER_MOL,
          UCationInversionBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintCationTorsionEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UCationTorsion;
  char string[256];
  REAL *parms,U;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,rbc;
  VECTOR Pb,Pc;
  VECTOR posA,posB,posC,posD;

  fprintf(OutputFilePtr[CurrentSystem],"UCationTorsion:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationTorsion=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfTorsions;i++)
    {
      A=Components[Type].Torsions[i].A;
      B=Components[Type].Torsions[i].B;
      C=Components[Type].Torsions[i].C;
      D=Components[Type].Torsions[i].D;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=Cations[CurrentSystem][m].Atoms[C].Position;
      posD=Cations[CurrentSystem][m].Atoms[D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dcd.x=posD.x-posC.x;
      Dcd.y=posD.y-posC.y;
      Dcd.z=posD.z-posC.z;

      dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;

      dr.x=Dab.x-dot_ab*Dbc.x;
      dr.y=Dab.y-dot_ab*Dbc.y;
      dr.z=Dab.z-dot_ab*Dbc.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Dcd.x-dot_cd*Dbc.x;
      ds.y=Dcd.y-dot_cd*Dbc.y;
      ds.z=Dcd.z-dot_cd*Dbc.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      // compute Cos(Phi)
      // Phi is defined in protein convention Phi(trans)=Pi
      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);

      parms=Components[Type].TorsionArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintTorsionEnergyStatus(nr++,string,Components[Type].TorsionType[i],parms,Phi);

      UCationTorsion+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation torsions: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationBTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationTorsion*ENERGY_TO_KELVIN,
          UCationTorsion*ENERGY_TO_KJ_PER_MOL,
          UCationTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}



void PrintCationImproperTorsionEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UCationImproperTorsion;
  char string[256];
  REAL *parms,U;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,rbc;
  VECTOR Pb,Pc;
  VECTOR posA,posB,posC,posD;

  fprintf(OutputFilePtr[CurrentSystem],"UCationImproperTorsion:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationImproperTorsion=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfImproperTorsions;i++)
    {
      A=Components[Type].ImproperTorsions[i].A;
      B=Components[Type].ImproperTorsions[i].B;
      C=Components[Type].ImproperTorsions[i].C;
      D=Components[Type].ImproperTorsions[i].D;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=Cations[CurrentSystem][m].Atoms[C].Position;
      posD=Cations[CurrentSystem][m].Atoms[D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dcd.x=posD.x-posC.x;
      Dcd.y=posD.y-posC.y;
      Dcd.z=posD.z-posC.z;

      dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;

      dr.x=Dab.x-dot_ab*Dbc.x;
      dr.y=Dab.y-dot_ab*Dbc.y;
      dr.z=Dab.z-dot_ab*Dbc.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Dcd.x-dot_cd*Dbc.x;
      ds.y=Dcd.y-dot_cd*Dbc.y;
      ds.z=Dcd.z-dot_cd*Dbc.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      // compute Cos(Phi)
      // Phi is defined in protein convention Phi(trans)=Pi
      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);

      parms=Components[Type].ImproperTorsionArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintImproperTorsionEnergyStatus(nr++,string,Components[Type].ImproperTorsionType[i],parms,Phi);

      UCationImproperTorsion+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation improper torsions: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationBImproperTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationImproperTorsion*ENERGY_TO_KELVIN,
          UCationImproperTorsion*ENERGY_TO_KJ_PER_MOL,
          UCationImproperTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintCationBondBondEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C;
  REAL UCationBondBond;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc;
  REAL rab,rbc;

  fprintf(OutputFilePtr[CurrentSystem],"UCationBondBond:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationBondBond=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBondBonds;i++)
    {
      A=Components[Type].BondBonds[i].A;
      B=Components[Type].BondBonds[i].B;
      C=Components[Type].BondBonds[i].C;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=Cations[CurrentSystem][m].Atoms[C].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));

      parms=Components[Type].BondBondArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s)",m,A,B,C,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name);

      U=PrintBondBondEnergyStatus(nr++,string,Components[Type].BondBondType[i],parms,rab,rbc);

      UCationBondBond+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation bond-bonds: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationBondBond Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationBondBond*ENERGY_TO_KELVIN,
          UCationBondBond*ENERGY_TO_KJ_PER_MOL,
          UCationBondBond*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintCationBondBendEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C;
  REAL UCationBondBend;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc;
  REAL rab,rbc,CosTheta,Theta;

  fprintf(OutputFilePtr[CurrentSystem],"UCationBondBend:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationBondBend=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBondBends;i++)
    {
      A=Components[Type].BondBends[i].A;
      B=Components[Type].BondBends[i].B;
      C=Components[Type].BondBends[i].C;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=Cations[CurrentSystem][m].Atoms[C].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
      Rab.x/=rab;
      Rab.y/=rab;
      Rab.z/=rab;

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
      Rbc.x/=rbc;
      Rbc.y/=rbc;
      Rbc.z/=rbc;

      CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
      Theta=acos(CosTheta);

      parms=Components[Type].BondBendArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s)",m,A,B,C,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name);

      U=PrintBondBendEnergyStatus(nr++,string,Components[Type].BondBendType[i],parms,rab,rbc,Theta);

      UCationBondBend+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation bond-bends: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationBondBend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationBondBend*ENERGY_TO_KELVIN,
          UCationBondBend*ENERGY_TO_KJ_PER_MOL,
          UCationBondBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintCationBendBendEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UCationBendBend;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC,posD;
  VECTOR Dab,Dbc,Dbd;
  REAL rab,rbc,rbd;
  REAL dot_abc,dot_abd;
  REAL Theta1,Theta2;
  REAL CosTheta1,CosTheta2;

  fprintf(OutputFilePtr[CurrentSystem],"UCationBendBend:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationBendBend=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBendBends;i++)
    {
      A=Components[Type].BendBends[i].A;
      B=Components[Type].BendBends[i].B;
      C=Components[Type].BendBends[i].C;
      D=Components[Type].BendBends[i].D;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=Cations[CurrentSystem][m].Atoms[C].Position;
      posD=Cations[CurrentSystem][m].Atoms[D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));
      Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dbd.x=posD.x-posB.x;
      Dbd.y=posD.y-posB.y;
      Dbd.z=posD.z-posB.z;
      rbd=sqrt(SQR(Dbd.x)+SQR(Dbd.y)+SQR(Dbd.z));
      Dbd.x/=rbd; Dbd.y/=rbd; Dbd.z/=rbd;

      dot_abc=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      CosTheta1=dot_abc;
      CosTheta1=SIGN(MIN2(fabs(CosTheta1),(REAL)1.0),CosTheta1);
      Theta1=acos(CosTheta1);

      dot_abd=Dab.x*Dbd.x+Dab.y*Dbd.y+Dab.z*Dbd.z;
      CosTheta2=dot_abd;
      CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
      Theta2=acos(CosTheta2);

      parms=Components[Type].BendBendArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintBendBendEnergyStatus(nr++,string,Components[Type].BendBendType[i],parms,Theta1,Theta2);

      UCationBendBend+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation bend-bends: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationBendBend Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationBendBend*ENERGY_TO_KELVIN,
          UCationBendBend*ENERGY_TO_KJ_PER_MOL,
          UCationBendBend*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintCationBondTorsionEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UCationBondTorsion;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC,posD;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL rbc,s,r,sign;
  REAL dot_ab,dot_cd,CosPhi,Phi;
  VECTOR Pb,Pc;

  fprintf(OutputFilePtr[CurrentSystem],"UCationBondTorsion:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationBondTorsion=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBondTorsions;i++)
    {
      A=Components[Type].BondTorsions[i].A;
      B=Components[Type].BondTorsions[i].B;
      C=Components[Type].BondTorsions[i].C;
      D=Components[Type].BondTorsions[i].D;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=Cations[CurrentSystem][m].Atoms[C].Position;
      posD=Cations[CurrentSystem][m].Atoms[D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dcd.x=posD.x-posC.x;
      Dcd.y=posD.y-posC.y;
      Dcd.z=posD.z-posC.z;

      dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;

      dr.x=Dab.x-dot_ab*Dbc.x;
      dr.y=Dab.y-dot_ab*Dbc.y;
      dr.z=Dab.z-dot_ab*Dbc.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Dcd.x-dot_cd*Dbc.x;
      ds.y=Dcd.y-dot_cd*Dbc.y;
      ds.z=Dcd.z-dot_cd*Dbc.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      // Ensure CosPhi is between -1 and 1.
      CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);

      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);


      parms=Components[Type].BondTorsionArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintBondTorsionEnergyStatus(nr++,string,Components[Type].BondTorsionType[i],parms,rbc,Phi);

      UCationBondTorsion+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation bond-torsions: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationBondTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationBondTorsion*ENERGY_TO_KELVIN,
          UCationBondTorsion*ENERGY_TO_KJ_PER_MOL,
          UCationBondTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintCationBendTorsionEnergyStatus(void)
{
  int m,i,nr,Type,A,B,C,D;
  REAL UCationBendTorsion;
  char string[256];
  REAL *parms,U;
  VECTOR posA,posB,posC,posD;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL rab,rbc,rcd,s,r,sign;
  REAL dot_ab,dot_cd,CosPhi,Phi;
  REAL CosTheta1,CosTheta2;
  REAL Theta1,Theta2;
  VECTOR Pb,Pc;

  fprintf(OutputFilePtr[CurrentSystem],"UCationBendTorsion:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationBendTorsion=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfBendTorsions;i++)
    {
      A=Components[Type].BendTorsions[i].A;
      B=Components[Type].BendTorsions[i].B;
      C=Components[Type].BendTorsions[i].C;
      D=Components[Type].BendTorsions[i].D;

      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=Cations[CurrentSystem][m].Atoms[C].Position;
      posD=Cations[CurrentSystem][m].Atoms[D].Position;

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

      Dbc.x=posC.x-posB.x;
      Dbc.y=posC.y-posB.y;
      Dbc.z=posC.z-posB.z;
      rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
      Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

      Dcd.x=posD.x-posC.x;
      Dcd.y=posD.y-posC.y;
      Dcd.z=posD.z-posC.z;
      rcd=sqrt(SQR(Dcd.x)+SQR(Dcd.y)+SQR(Dcd.z));

      dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
      CosTheta1=dot_ab/rab;
      CosTheta1=SIGN(MIN2(fabs(CosTheta1),(REAL)1.0),CosTheta1);
      Theta1=acos(CosTheta1);

      dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;
      CosTheta2=-dot_cd/rcd;
      CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
      Theta2=acos(CosTheta2);

      dr.x=Dab.x-dot_ab*Dbc.x;
      dr.y=Dab.y-dot_ab*Dbc.y;
      dr.z=Dab.z-dot_ab*Dbc.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Dcd.x-dot_cd*Dbc.x;
      ds.y=Dcd.y-dot_cd*Dbc.y;
      ds.z=Dcd.z-dot_cd*Dbc.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      // Ensure CosPhi is between -1 and 1.
      CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);

      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);

      parms=Components[Type].BendTorsionArguments[i];

      sprintf(string,"(A%-4d,[%-4d,%-4d,%-4d,%-4d]) (%-6s,%-6s,%-6s,%-6s)",m,A,B,C,D,
                     PseudoAtoms[Components[Type].Type[A]].Name,
                     PseudoAtoms[Components[Type].Type[B]].Name,
                     PseudoAtoms[Components[Type].Type[C]].Name,
                     PseudoAtoms[Components[Type].Type[D]].Name);

      U=PrintBendTorsionEnergyStatus(nr++,string,Components[Type].BendTorsionType[i],parms,Theta1,Theta2,Phi);

      UCationBendTorsion+=U;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation bend-torsions: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationBendTorsion Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationBendTorsion*ENERGY_TO_KELVIN,
          UCationBendTorsion*ENERGY_TO_KJ_PER_MOL,
          UCationBendTorsion*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintCationIntraVDWEnergyStatus(void)
{
  int i,NumberOfIntraLJ;
  int Type,TypeA,TypeB,A,B;
  int nr,m;
  REAL rr,Scaling;
  VECTOR dr;
  POINT posA,posB;
  REAL UCationIntraVDW,energy;
  char string[256];

  fprintf(OutputFilePtr[CurrentSystem],"UCationIntraVDW:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationIntraVDW=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    NumberOfIntraLJ=Components[Type].NumberOfIntraVDW;
    for(i=0;i<NumberOfIntraLJ;i++)
    {
      A=Components[Type].IntraVDW[i].A;
      B=Components[Type].IntraVDW[i].B;
      Scaling=Components[Type].IntraVDWScaling[i];
      TypeA=Cations[CurrentSystem][m].Atoms[A].Type;
      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      if(rr<CutOffVDWSquared)
      {
        TypeB=Cations[CurrentSystem][m].Atoms[B].Type;
        energy=Scaling*PotentialValue(TypeA,TypeB,rr,1.0);

        sprintf(string,"(Scaling: %8.5f) (A%-4d,[%-4d,%-4d]) (%-6s,%-6s)",Scaling,m,A,B,PseudoAtoms[TypeA].Name,PseudoAtoms[TypeB].Name);
        PrintVDWEnergyStatus(nr,string,TypeA,TypeB,sqrt(rr),energy);

        UCationIntraVDW+=energy;
      }
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation intra-VDW: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationIntraVDW Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationIntraVDW*ENERGY_TO_KELVIN,
          UCationIntraVDW*ENERGY_TO_KJ_PER_MOL,
          UCationIntraVDW*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}



void PrintCationIntraChargeChargeEnergyStatus(void)
{
  int m,i,NumberOfIntraChargeChargeCoulomb,Type,TypeA,TypeB,A,B;
  REAL r,rr,chargeA,chargeB,Scaling;
  REAL UCationIntraChargeCharge,energy;
  VECTOR dr;
  POINT posA,posB;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UCationIntraChargeCharge:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationIntraChargeCharge=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    NumberOfIntraChargeChargeCoulomb=Components[Type].NumberOfIntraChargeCharge;
    for(i=0;i<NumberOfIntraChargeChargeCoulomb;i++)
    {
      A=Components[Type].IntraChargeCharge[i].A;
      B=Components[Type].IntraChargeCharge[i].B;
      Scaling=Components[Type].IntraChargeChargeScaling[i];
      TypeA=Cations[CurrentSystem][m].Atoms[A].Type;
      chargeA=Cations[CurrentSystem][m].Atoms[A].Charge;
      posA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=Cations[CurrentSystem][m].Atoms[B].Position;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      TypeB=Cations[CurrentSystem][m].Atoms[B].Type;
      chargeB=Cations[CurrentSystem][m].Atoms[B].Charge;
      energy=Scaling*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;

      sprintf(string,"(A%-4d,[%-4d,%-4d)) (%-6s,%-6s)",m,A,B,PseudoAtoms[TypeA].Name,PseudoAtoms[TypeB].Name);
      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_COULOMBIC %s, Charge: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,chargeA,chargeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);

      UCationIntraChargeCharge+=energy;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation intra-charge-charge: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationIntraChargeCharge Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationIntraChargeCharge*ENERGY_TO_KELVIN,
          UCationIntraChargeCharge*ENERGY_TO_KJ_PER_MOL,
          UCationIntraChargeCharge*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}


void PrintCationIntraChargeBondDipoleEnergyStatus(void)
{
  int m,i,NumberOfIntraChargeBondDipoleCoulomb,Type,A,B,B1,B2;
  REAL r,rr,ri2,length;
  REAL Bt1,cosB;
  VECTOR dr,dipoleB;
  POINT posA,posB,posB1,posB2;
  REAL ChargeA,temp,UCationIntraChargeBondDipole;
  REAL DipoleMagnitudeB,energy;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UCationIntraChargeBondDipole:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationIntraChargeBondDipole=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    NumberOfIntraChargeBondDipoleCoulomb=Components[Type].NumberOfIntraChargeBondDipole;
    for(i=0;i<NumberOfIntraChargeBondDipoleCoulomb;i++)
    {
      A=Components[Type].IntraChargeBondDipole[i].A;
      B=Components[Type].IntraChargeBondDipole[i].B;

      ChargeA=Cations[CurrentSystem][m].Atoms[A].Charge;
      posA=Cations[CurrentSystem][m].Atoms[A].Position;

      B1=Components[Type].BondDipoles[B].A;
      B2=Components[Type].BondDipoles[B].B;
      DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
      posB1=Cations[CurrentSystem][m].Atoms[B1].Position;
      posB2=Cations[CurrentSystem][m].Atoms[B2].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      ri2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeB/length;
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt1=1.0/(r*rr);
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy=-COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);

      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_COULOMB %s, Magnitude: %8.5f, Charge: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,DipoleMagnitudeB,ChargeA,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);

      UCationIntraChargeBondDipole+=energy;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cation intra-charge-bonddipole: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationIntraChargeBondDipole Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationIntraChargeBondDipole*ENERGY_TO_KELVIN,
          UCationIntraChargeBondDipole*ENERGY_TO_KJ_PER_MOL,
          UCationIntraChargeBondDipole*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintCationIntraBondDipoleBondDipoleEnergyStatus(void)
{
  int i,m,NumberOfIntraBondDipoleBondDipoleCoulomb,Type,A,B,A1,A2,B1,B2;
  REAL r,rr,ri2,rk2,length;
  REAL Bt0,Bt1,Bt2,Bt3,cosA,cosB,cosAB;
  VECTOR dr,dipoleA,dipoleB;
  POINT posA,posB,posA1,posA2,posB1,posB2;
  REAL temp,UCationIntraBondDipoleBondDipole;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,energy;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UCationIntraBondDipoleBondDipole:\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UCationIntraBondDipoleBondDipole=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    NumberOfIntraBondDipoleBondDipoleCoulomb=Components[Type].NumberOfIntraBondDipoleBondDipole;
    for(i=0;i<NumberOfIntraBondDipoleBondDipoleCoulomb;i++)
    {
      A=Components[Type].IntraBondDipoleBondDipole[i].A;
      B=Components[Type].IntraBondDipoleBondDipole[i].B;

      A1=Components[Type].BondDipoles[A].A;
      A2=Components[Type].BondDipoles[A].B;
      posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
      DipoleMagnitudeA=Components[Type].BondDipoleMagnitude[A];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      B1=Components[Type].BondDipoles[B].A;
      B2=Components[Type].BondDipoles[B].B;
      posB1=Cations[CurrentSystem][m].Atoms[B1].Position;
      posB2=Cations[CurrentSystem][m].Atoms[B2].Position;
      DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
      length=sqrt(rk2);
      temp=DipoleMagnitudeB/length;
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=1.0/(r);
      Bt1=1.0/(r*rr);
      Bt2=3.0/(r*rr*rr);
      Bt3=15.0/(r*rr*rr*rr);

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);

      fprintf(OutputFilePtr[CurrentSystem],"%4d TRUNCATED_COULOMB %s, Magnitude: %8.5f, Magnitude: %8.5f, Distance %8.5f [A], Energy: %10.5f [K] %8.5f [kJ/mol] %8.5f [kcal/mol]\n",
          nr,string,DipoleMagnitudeA,DipoleMagnitudeB,r,
          energy*ENERGY_TO_KELVIN,
          energy*ENERGY_TO_KJ_PER_MOL,
          energy*ENERGY_TO_KCAL_PER_MOL);

      UCationIntraBondDipoleBondDipole+=energy;
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\nNumber of cations intra-bonddipole-bonddipole: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationIntraBondDipoleBondDipole Energy: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationIntraBondDipoleBondDipole*ENERGY_TO_KELVIN,
          UCationIntraBondDipoleBondDipole*ENERGY_TO_KJ_PER_MOL,
          UCationIntraBondDipoleBondDipole*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

// =============================================================================================================================
// inter-molecular energies
// =============================================================================================================================


void PrintInterVDWEnergyStatus(void)
{
  int i,j,k,l;
  int typeA,typeB,TypeMolA,TypeMolB;
  REAL rr,energy;
  REAL UAdsorbateAdsorbateVDW;
  REAL UAdsorbateCationVDW;
  REAL UCationCationVDW;
  VECTOR posA,posB,dr;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UInterVDW\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateAdsorbateVDW=0.0;
  UAdsorbateCationVDW=0.0;
  UCationCationVDW=0.0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeMolA=Adsorbates[CurrentSystem][i].Type;
    for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[k].AnisotropicPosition;

      // loop over adsorbant molecules
      if(!OmitAdsorbateAdsorbateVDWInteractions)
      {
        for(j=i+1;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
        {
          TypeMolB=Adsorbates[CurrentSystem][j].Type;
          for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
          {
            typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
            posB=Adsorbates[CurrentSystem][j].Atoms[l].AnisotropicPosition;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              energy=PotentialValue(typeA,typeB,rr,1.0);
              UAdsorbateAdsorbateVDW+=energy;

              sprintf(string,"(A%-4d,%-4d)-(C%-4d,%-4d) (%-6s,%-6s)",i,k,j,l,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
              PrintVDWEnergyStatus(nr++,string,typeA,typeB,sqrt(rr),energy);
            }
          }
        }
      }
      for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
      {
        TypeMolB=Cations[CurrentSystem][j].Type;
        for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
        {
          posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
          typeB=Cations[CurrentSystem][j].Atoms[l].Type;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
          {
            energy=PotentialValue(typeA,typeB,rr,1.0);
            UAdsorbateCationVDW+=energy;
            sprintf(string,"(A%-4d,%-4d)-(A%-4d,%-4d) (%-6s,%-6s)",i,k,j,l,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
            PrintVDWEnergyStatus(nr++,string,typeA,typeB,sqrt(rr),energy);
          }
        }
      }
    }
  }

  if(!OmitCationCationVDWInteractions)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      TypeMolA=Cations[CurrentSystem][i].Type;
      for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Cations[CurrentSystem][i].Atoms[k].Type;
        posA=Cations[CurrentSystem][i].Atoms[k].AnisotropicPosition;

        // loop over cation molecules
        for(j=i+1;j<NumberOfCationMolecules[CurrentSystem];j++)
        {
          TypeMolB=Cations[CurrentSystem][j].Type;
          for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
          {
            posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
            typeB=Cations[CurrentSystem][j].Atoms[l].Type;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
            {
              energy=PotentialValue(typeA,typeB,rr,1.0);
              UCationCationVDW+=energy;
              sprintf(string,"(C%-4d,%-4d)-(C%-4d,%-4d) (%-6s,%-6s)",i,k,j,l,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
              PrintVDWEnergyStatus(nr++,string,typeA,typeB,sqrt(rr),energy);
            }
          }
        }
      }
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateAdsorbateVDW: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateAdsorbateVDW*ENERGY_TO_KELVIN,
          UAdsorbateAdsorbateVDW*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateAdsorbateVDW*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateCationVDW: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateCationVDW*ENERGY_TO_KELVIN,
          UAdsorbateCationVDW*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateCationVDW*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationCationVDW: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationCationVDW*ENERGY_TO_KELVIN,
          UCationCationVDW*ENERGY_TO_KJ_PER_MOL,
          UCationCationVDW*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintInterChargeChargeEnergyStatus(void)
{
  int i,j,k,l;
  int typeA,typeB;
  REAL r,r2;
  REAL chargeA,chargeB;
  REAL energy;
  VECTOR posA,posB,dr;
  REAL NetChargeB,UWolfCorrection;
  REAL UAdsorbateAdsorbateChargeChargeReal;
  REAL UAdsorbateCationChargeChargeReal;
  REAL UCationCationChargeChargeReal;
  char string[256];
  int nr;

  nr=0;
  UAdsorbateAdsorbateChargeChargeReal=0.0;
  UAdsorbateCationChargeChargeReal=0.0;
  UCationCationChargeChargeReal=0.0;

  fprintf(OutputFilePtr[CurrentSystem],"UInterChargeChargeReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[k].Position;
      chargeA=Adsorbates[CurrentSystem][i].Atoms[k].Charge;

      // loop over adsorbant molecules
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(j=i+1;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
        {
          for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
          {
            posB=Adsorbates[CurrentSystem][j].Atoms[l].Position;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(r2<CutOffChargeChargeSquared[CurrentSystem])
            {
              r=sqrt(r2);
              typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
              chargeB=Adsorbates[CurrentSystem][j].Atoms[l].Charge;

              sprintf(string,"(A%-4d,%-4d)-(A%-4d,%-4d) (%-6s,%-6s)",i,k,j,l,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
              energy=PrintChargeChargeEnergyStatus(nr++,string,chargeA,chargeB,r);

              UAdsorbateAdsorbateChargeChargeReal+=energy;
            }
          }
        }
      }

      for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
      {
        for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
        {
          posB=Cations[CurrentSystem][j].Atoms[l].Position;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(r2<CutOffChargeChargeSquared[CurrentSystem])
          {
            r=sqrt(r2);
            typeB=Cations[CurrentSystem][j].Atoms[l].Type;
            chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;

            sprintf(string,"(A%-4d,%-4d)-(C%-4d,%-4d) (%-6s,%-6s)",i,k,j,l,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
            energy=PrintChargeChargeEnergyStatus(nr++,string,chargeA,chargeB,r);

            UAdsorbateCationChargeChargeReal+=energy;
          }
        }
      }
    }
  }

  if(!OmitCationCationCoulombInteractions)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Cations[CurrentSystem][i].Atoms[k].Type;
        posA=Cations[CurrentSystem][i].Atoms[k].Position;
        chargeA=Cations[CurrentSystem][i].Atoms[k].Charge;

        // loop over cation molecules
        for(j=i+1;j<NumberOfCationMolecules[CurrentSystem];j++)
        {
          for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
          {
            posB=Cations[CurrentSystem][j].Atoms[l].Position;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(r2<CutOffChargeChargeSquared[CurrentSystem])
            {
              r=sqrt(r2);
              typeB=Cations[CurrentSystem][j].Atoms[l].Type;
              chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;

              sprintf(string,"(C%-4d,%-4d)-(C%-4d,%-4d) (%-6s,%-6s)",i,k,j,l,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
              energy=PrintChargeChargeEnergyStatus(nr++,string,chargeA,chargeB,r);

              UCationCationChargeChargeReal+=energy;
            }
          }
        }
      }
    }
  }

  switch(ChargeMethod)
  {
    case WOLFS_METHOD:
      UWolfCorrection=0.0;
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=0.5*COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*InverseCutOffChargeCharge[CurrentSystem];
      }
      UAdsorbateAdsorbateChargeChargeReal+=UWolfCorrection;

      UWolfCorrection=0.0;
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Cations[CurrentSystem][i].Atoms[j].Type;
          chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=0.5*COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*InverseCutOffChargeCharge[CurrentSystem];
      }
      UCationCationChargeChargeReal+=UWolfCorrection;
      break;
    case WOLFS_METHOD_DAMPED:
    case WOLFS_METHOD_DAMPED_FG:
      UWolfCorrection=0.0;
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*
             (0.5*erfc(Alpha[CurrentSystem]*CutOffChargeCharge[CurrentSystem])*InverseCutOffChargeCharge[CurrentSystem]+
              Alpha[CurrentSystem]*M_1_SQRTPI);
      }
      UAdsorbateAdsorbateChargeChargeReal+=UWolfCorrection;

      UWolfCorrection=0.0;
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Cations[CurrentSystem][i].Atoms[j].Type;
          chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*
             (0.5*erfc(Alpha[CurrentSystem]*CutOffChargeCharge[CurrentSystem])*InverseCutOffChargeCharge[CurrentSystem]+
              Alpha[CurrentSystem]*M_1_SQRTPI);
      }
      UCationCationChargeChargeReal+=UWolfCorrection;
      break;
    default:
      break;
  }

  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateAdsorbateChargeChargeReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateAdsorbateChargeChargeReal*ENERGY_TO_KELVIN,
          UAdsorbateAdsorbateChargeChargeReal*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateAdsorbateChargeChargeReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateCationChargeChargeReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateCationChargeChargeReal*ENERGY_TO_KELVIN,
          UAdsorbateCationChargeChargeReal*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateCationChargeChargeReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationCationChargeChargeReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationCationChargeChargeReal*ENERGY_TO_KELVIN,
          UCationCationChargeChargeReal*ENERGY_TO_KJ_PER_MOL,
          UCationCationChargeChargeReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}

void PrintInterChargeBondDipoleEnergyStatus(void)
{
  int i,j,k,l;
  int A1,A2;
  int TypeA,typeA1,typeA2,typeB;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL r,rr,ri2,cosA,energy,temp,length,ChargeB;
  VECTOR dipoleA;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA;
  REAL UAdsorbateAdsorbateChargeBondDipoleReal;
  REAL UAdsorbateCationChargeBondDipoleReal;
  REAL UCationCationChargeBondDipoleReal;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UInterChargeBondDipoleReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UAdsorbateAdsorbateChargeBondDipoleReal=0.0;
  UAdsorbateCationChargeBondDipoleReal=0.0;
  UCationCationChargeBondDipoleReal=0.0;
  Bt0=Bt1=Bt2=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeA=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      posA1=Adsorbates[CurrentSystem][i].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][i].Atoms[A2].Position;
      typeA1=Adsorbates[CurrentSystem][i].Atoms[A1].Type;
      typeA2=Adsorbates[CurrentSystem][i].Atoms[A2].Type;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(i!=k)
          {
            for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
            {
              typeB=Adsorbates[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[typeB].HasCharges)
              {
                posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeBondDipoleSquared)
                {
                  r=sqrt(rr);
                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

                  sprintf(string,"(A%-4d,[%-4d,%-4d])-(A%-4d,%-4d) ([%-6s,%-6s],%-6s)",i,A1,A2,k,l,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,PseudoAtoms[typeB].Name);
                  energy=PrintChargeBondDipoleStatus(nr++,string,cosA,DipoleMagnitudeA,ChargeB,r);

                  UAdsorbateAdsorbateChargeBondDipoleReal+=energy;
                }
              }
            }
          }
        }
      }

      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
        {
          typeB=Cations[CurrentSystem][k].Atoms[l].Type;
          if(PseudoAtoms[typeB].HasCharges)
          {
            posB=Cations[CurrentSystem][k].Atoms[l].Position;
            ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              sprintf(string,"(A%-4d,[%-4d,%-4d])-(C%-4d,%-4d) ([%-6s,%-6s],%-6s)",i,A1,A2,k,l,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,PseudoAtoms[typeB].Name);
              energy=PrintChargeBondDipoleStatus(nr++,string,cosA,DipoleMagnitudeA,ChargeB,r);
              UAdsorbateCationChargeBondDipoleReal+=energy;
            }
          }
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    TypeA=Cations[CurrentSystem][i].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      posA1=Cations[CurrentSystem][i].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][i].Atoms[A2].Position;
      typeA1=Cations[CurrentSystem][i].Atoms[A1].Type;
      typeA2=Cations[CurrentSystem][i].Atoms[A2].Type;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        if(i!=k)
        {
          if(!OmitCationCationCoulombInteractions)
          {
            for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
            {
              typeB=Cations[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[typeB].HasCharges)
              {
                posB=Cations[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeBondDipoleSquared)
                {
                  r=sqrt(rr);
                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                  sprintf(string,"(C%-4d,[%-4d,%-4d])-(C%-4d,%-4d) ([%-6s,%-6s],%-6s)",i,A1,A2,k,l,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,PseudoAtoms[typeB].Name);
                  energy=PrintChargeBondDipoleStatus(nr++,string,cosA,DipoleMagnitudeA,ChargeB,r);
                  UCationCationChargeBondDipoleReal+=energy;
                }
              }
            }
          }
        }
      }

      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
        {
          typeB=Adsorbates[CurrentSystem][k].Atoms[l].Type;
          if(PseudoAtoms[typeB].HasCharges)
          {
            posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
            ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              sprintf(string,"(C%-4d,[%-4d,%-4d])-(A%-4d,%-4d) ([%-6s,%-6s],%-6s)",i,A1,A2,k,l,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,PseudoAtoms[typeB].Name);
              energy=PrintChargeBondDipoleStatus(nr++,string,cosA,DipoleMagnitudeA,ChargeB,r);
              UAdsorbateCationChargeBondDipoleReal+=energy;
            }
          }
        }
      }
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateAdsorbateChargeBondDipoleReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateAdsorbateChargeBondDipoleReal*ENERGY_TO_KELVIN,
          UAdsorbateAdsorbateChargeBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateAdsorbateChargeBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateCationChargeBondDipoleReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateCationChargeBondDipoleReal*ENERGY_TO_KELVIN,
          UAdsorbateCationChargeBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateCationChargeBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationCationChargeBondDipoleReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationCationChargeBondDipoleReal*ENERGY_TO_KELVIN,
          UCationCationChargeBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UCationCationChargeBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}

void PrintInterBondDipoleBondDipoleEnergyStatus(void)
{
  int i,j,k,l;
  int A1,A2,B1,B2;
  int TypeA,TypeB;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,rr,ri2,rk2,cosAB,cosA,cosB,energy,temp,length;
  VECTOR dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL UAdsorbateAdsorbateBondDipoleBondDipoleReal;
  REAL UAdsorbateCationBondDipoleBondDipoleReal;
  REAL UCationCationBondDipoleBondDipoleReal;
  int typeA1,typeA2,typeB1,typeB2;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UInterBondDipoleBondDipoleReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  Bt0=Bt1=Bt2=0.0;
  UAdsorbateAdsorbateBondDipoleBondDipoleReal=0.0;
  UAdsorbateCationBondDipoleBondDipoleReal=0.0;
  UCationCationBondDipoleBondDipoleReal=0.0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeA=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      posA1=Adsorbates[CurrentSystem][i].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][i].Atoms[A2].Position;
      typeA1=Adsorbates[CurrentSystem][i].Atoms[A1].Type;
      typeA2=Adsorbates[CurrentSystem][i].Atoms[A2].Type;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=i+1;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          TypeB=Adsorbates[CurrentSystem][k].Type;
          for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
          {
            B1=Components[TypeB].BondDipoles[l].A;
            B2=Components[TypeB].BondDipoles[l].B;
            posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
            posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
            typeB1=Adsorbates[CurrentSystem][k].Atoms[B1].Type;
            typeB2=Adsorbates[CurrentSystem][k].Atoms[B2].Type;
            DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
            dipoleB.x=posB2.x-posB1.x;
            dipoleB.y=posB2.y-posB1.y;
            dipoleB.z=posB2.z-posB1.z;
            posB.x=posB1.x+0.5*dipoleB.x;
            posB.y=posB1.y+0.5*dipoleB.y;
            posB.z=posB1.z+0.5*dipoleB.z;
            rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
            length=sqrt(rk2);
            temp=DipoleMagnitudeB/length;
            dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              r=sqrt(rr);

              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

              sprintf(string,"(F%1d,[%-4d,%-4d])-(F%1d,[%-4d,-%4d]) ([%-6s,%-6s],[%-6s,%-6s])",i,A1,A2,k,B1,B2,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,
                               PseudoAtoms[typeB1].Name,PseudoAtoms[typeB2].Name);
              energy=PrintBondDipoleBondDipoleStatus(nr++,string,cosA,cosB,cosAB,DipoleMagnitudeA,DipoleMagnitudeB,r);
              UAdsorbateAdsorbateBondDipoleBondDipoleReal+=energy;
            }
          }
        }
      }

      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        TypeB=Cations[CurrentSystem][k].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
          posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
          typeB1=Cations[CurrentSystem][k].Atoms[B1].Type;
          typeB2=Cations[CurrentSystem][k].Atoms[B2].Type;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(rr);

            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            sprintf(string,"(F%1d,[%-4d,%-4d])-(F%1d,[%-4d,-%4d]) ([%-6s,%-6s],[%-6s,%-6s])",i,A1,A2,k,B1,B2,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,
                             PseudoAtoms[typeB1].Name,PseudoAtoms[typeB2].Name);
            energy=PrintBondDipoleBondDipoleStatus(nr++,string,cosA,cosB,cosAB,DipoleMagnitudeA,DipoleMagnitudeB,r);

            UAdsorbateCationBondDipoleBondDipoleReal+=energy;
          }
        }
      }
    }
  }

  if(!OmitCationCationCoulombInteractions)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      TypeA=Cations[CurrentSystem][i].Type;
      for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
      {
        A1=Components[TypeA].BondDipoles[j].A;
        A2=Components[TypeA].BondDipoles[j].B;
        posA1=Cations[CurrentSystem][i].Atoms[A1].Position;
        posA2=Cations[CurrentSystem][i].Atoms[A2].Position;
        typeA1=Cations[CurrentSystem][i].Atoms[A1].Type;
        typeA2=Cations[CurrentSystem][i].Atoms[A2].Type;
        DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
        length=sqrt(ri2);
        temp=DipoleMagnitudeA/length;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        for(k=i+1;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          TypeB=Cations[CurrentSystem][k].Type;
          for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
          {
            B1=Components[TypeB].BondDipoles[l].A;
            B2=Components[TypeB].BondDipoles[l].B;
            posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
            posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
            typeB1=Cations[CurrentSystem][k].Atoms[B1].Type;
            typeB2=Cations[CurrentSystem][k].Atoms[B2].Type;
            DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
            dipoleB.x=posB2.x-posB1.x;
            dipoleB.y=posB2.y-posB1.y;
            dipoleB.z=posB2.z-posB1.z;
            posB.x=posB1.x+0.5*dipoleB.x;
            posB.y=posB1.y+0.5*dipoleB.y;
            posB.z=posB1.z+0.5*dipoleB.z;
            rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
            length=sqrt(rk2);
            temp=DipoleMagnitudeB/length;
            dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              r=sqrt(rr);
              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              sprintf(string,"(F%1d,[%-4d,%-4d])-(F%1d,[%-4d,-%4d]) ([%-6s,%-6s],[%-6s,%-6s])",i,A1,A2,k,B1,B2,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,
                               PseudoAtoms[typeB1].Name,PseudoAtoms[typeB2].Name);
              energy=PrintBondDipoleBondDipoleStatus(nr++,string,cosA,cosB,cosAB,DipoleMagnitudeA,DipoleMagnitudeB,r);

              UCationCationBondDipoleBondDipoleReal+=energy;
            }
          }
        }
      }
    }
  }

  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateAdsorbateBondDipoleBondDipoleReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateAdsorbateBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,
          UAdsorbateAdsorbateBondDipoleBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateAdsorbateBondDipoleBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"Total UAdsorbateCationBondDipoleBondDipoleReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UAdsorbateCationBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,
          UAdsorbateCationBondDipoleBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UAdsorbateCationBondDipoleBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"Total UCationCationBondDipoleBondDipoleReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UCationCationBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,
          UCationCationBondDipoleBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UCationCationBondDipoleBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintFrameworkAdsorbateVDWEnergyStatus(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,UHostAdsorbateVDW,energy;
  VECTOR posA,posB,dr;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostAdsorbateVDW\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostAdsorbateVDW=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
         for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
          typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffVDWSquared)
          {
            energy=PotentialValue(typeA,typeB,rr,1.0);
            UHostAdsorbateVDW+=energy;
            sprintf(string,"(A%-4d,%-4d)-(F%1d,%-4d) (%-6s,%-6s)",i,j,f1,k,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
            PrintVDWEnergyStatus(nr++,string,typeA,typeB,sqrt(rr),energy);
          }
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostAdsorbateVDW: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostAdsorbateVDW*ENERGY_TO_KELVIN,
          UHostAdsorbateVDW*ENERGY_TO_KJ_PER_MOL,
          UHostAdsorbateVDW*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintFrameworkAdsorbateChargeChargeEnergyStatus(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  REAL r,rr,energy;
  REAL chargeA,chargeB;
  VECTOR posA,posB,dr;
  REAL UHostAdsorbateChargeChargeReal;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostAdsorbateChargeChargeReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0.0;
  UHostAdsorbateChargeChargeReal=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            r=sqrt(rr);
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

            sprintf(string,"(A%-4d,%-4d)-(F%1d,%-4d) (%-6s,%-6s)",i,j,f1,k,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
            energy=PrintChargeChargeEnergyStatus(nr++,string,chargeA,chargeB,r);

            UHostAdsorbateChargeChargeReal+=energy;
          }
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostAdsorbateChargeChargeReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostAdsorbateChargeChargeReal*ENERGY_TO_KELVIN,
          UHostAdsorbateChargeChargeReal*ENERGY_TO_KJ_PER_MOL,
          UHostAdsorbateChargeChargeReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintFrameworkAdsorbateChargeBondDipoleEnergyStatus(void)
{
  int i,j,k,f1;
  int A1,A2;
  int typeA1,typeA2,typeB,TypeA;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL r,rr,ri2,cosA,energy,temp,length,chargeB;
  VECTOR dipoleA;
  REAL Bt0,Bt1;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;
  REAL UHostAdsorbateChargeBondDipoleReal;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostAdsorbateChargeBondDipoleReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostAdsorbateChargeBondDipoleReal=0.0;

  Bt0=Bt1=0.0;
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
      A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
      A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
      posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
      typeA1=Framework[CurrentSystem].Atoms[f1][A1].Type;
      typeA2=Framework[CurrentSystem].Atoms[f1][A2].Type;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
      {
        for(k=0;k<Adsorbates[CurrentSystem][j].NumberOfAtoms;k++)
        {
          typeB=Adsorbates[CurrentSystem][j].Atoms[k].Type;
          if(PseudoAtoms[typeB].HasCharges)
          {
            posB=Adsorbates[CurrentSystem][j].Atoms[k].Position;
            chargeB=Adsorbates[CurrentSystem][j].Atoms[k].Charge;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

              sprintf(string,"(F%1d,[%-4d,%-4d])-(A%-4d,%-4d) ([%-6s,%-6s],%-6s)",f1,A1,A2,j,k,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,PseudoAtoms[typeB].Name);
              energy=PrintChargeBondDipoleStatus(nr++,string,cosA,DipoleMagnitudeA,chargeB,r);

              UHostAdsorbateChargeBondDipoleReal+=energy;
            }
          }
        }
      }
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeA=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      posA1=Adsorbates[CurrentSystem][i].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][i].Atoms[A2].Position;
      typeA1=Adsorbates[CurrentSystem][i].Atoms[A1].Type;
      typeA2=Adsorbates[CurrentSystem][i].Atoms[A2].Type;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

            sprintf(string,"(C%-4d,[%-4d,%-4d])-(A%-4d,%-4d) ([%-6s,%-6s],%-6s)",f1,A1,A2,i,j,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,PseudoAtoms[typeB].Name);
            energy=PrintChargeBondDipoleStatus(nr++,string,cosA,DipoleMagnitudeA,chargeB,r);
            UHostAdsorbateChargeBondDipoleReal+=energy;
          }
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostAdsorbateChargeBondDipoleReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostAdsorbateChargeBondDipoleReal*ENERGY_TO_KELVIN,
          UHostAdsorbateChargeBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UHostAdsorbateChargeBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}

void PrintFrameworkAdsorbateBondDipoleBondDipoleEnergyStatus(void)
{
  int i,k,l,f1;
  int A1,A2,B1,B2;
  int TypeB,typeA1,typeA2,typeB1,typeB2;;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,r2,ri2,rk2,cosA,cosB,cosAB,energy,temp,length;
  VECTOR dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL UHostAdsorbateBondDipoleBondDipoleReal;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostAdsorbateBondDipoleBondDipoleReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostAdsorbateBondDipoleBondDipoleReal=0.0;

  Bt0=Bt1=Bt2=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
      A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
      A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
      posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
      typeA1=Framework[CurrentSystem].Atoms[f1][A1].Type;
      typeA2=Framework[CurrentSystem].Atoms[f1][A2].Type;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        TypeB=Adsorbates[CurrentSystem][k].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
          posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
          typeB1=Adsorbates[CurrentSystem][k].Atoms[B1].Type;
          typeB2=Adsorbates[CurrentSystem][k].Atoms[B2].Type;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(r2<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(r2);
            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

            sprintf(string,"(F%1d,[%-4d,%-4d])-(F%1d,[%-4d,-%4d]) ([%-6s,%-6s],[%-6s,%-6s])",i,A1,A2,k,B1,B2,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,
                             PseudoAtoms[typeB1].Name,PseudoAtoms[typeB2].Name);
            energy=PrintBondDipoleBondDipoleStatus(nr++,string,cosA,cosB,cosAB,DipoleMagnitudeA,DipoleMagnitudeB,r);
            UHostAdsorbateBondDipoleBondDipoleReal+=energy;
          }
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostAdsorbateBondDipoleBondDipoleReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostAdsorbateBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,
          UHostAdsorbateBondDipoleBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UHostAdsorbateBondDipoleBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}


void PrintFrameworkCationVDWEnergyStatus(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,UHostCationVDW,energy;
  VECTOR posA,posB,dr;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostCationVDW\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostCationVDW=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
         for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
          typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffVDWSquared)
          {
            energy=PotentialValue(typeA,typeB,rr,1.0);
            UHostCationVDW+=energy;
            sprintf(string,"(A%-4d,%-4d)-(F%1d,%-4d) (%-6s,%-6s)",i,j,f1,k,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
            PrintVDWEnergyStatus(nr++,string,typeA,typeB,sqrt(rr),energy);
          }
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostCationVDW: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostCationVDW*ENERGY_TO_KELVIN,
          UHostCationVDW*ENERGY_TO_KJ_PER_MOL,
          UHostCationVDW*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}


void PrintFrameworkCationChargeChargeEnergyStatus(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  REAL r,rr,energy;
  REAL chargeA,chargeB;
  VECTOR posA,posB,dr;
  REAL UHostCationChargeChargeReal;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostCationChargeChargeReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0.0;
  UHostCationChargeChargeReal=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;
      posA=Cations[CurrentSystem][i].Atoms[j].Position;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            r=sqrt(rr);
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

            sprintf(string,"(A%-4d,%-4d)-(F%1d,%-4d) (%-6s,%-6s)",i,j,f1,k,PseudoAtoms[typeA].Name,PseudoAtoms[typeB].Name);
            energy=PrintChargeChargeEnergyStatus(nr++,string,chargeA,chargeB,r);

            UHostCationChargeChargeReal+=energy;
          }
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostCationChargeChargeReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostCationChargeChargeReal*ENERGY_TO_KELVIN,
          UHostCationChargeChargeReal*ENERGY_TO_KJ_PER_MOL,
          UHostCationChargeChargeReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void PrintFrameworkCationChargeBondDipoleEnergyStatus(void)
{
  int i,j,k,f1;
  int A1,A2;
  int typeA1,typeA2,typeB,TypeA;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL r,rr,ri2,cosA,energy,temp,length,chargeB;
  VECTOR dipoleA;
  REAL Bt0,Bt1;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;
  REAL UHostCationChargeBondDipoleReal;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostCationChargeBondDipoleReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostCationChargeBondDipoleReal=0.0;

  Bt0=Bt1=0.0;
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
      A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
      A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
      posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
      typeA1=Framework[CurrentSystem].Atoms[f1][A1].Type;
      typeA2=Framework[CurrentSystem].Atoms[f1][A2].Type;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
      {
        for(k=0;k<Cations[CurrentSystem][j].NumberOfAtoms;k++)
        {
          typeB=Cations[CurrentSystem][j].Atoms[k].Type;
          if(PseudoAtoms[typeB].HasCharges)
          {
            posB=Cations[CurrentSystem][j].Atoms[k].Position;
            chargeB=Cations[CurrentSystem][j].Atoms[k].Charge;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

              sprintf(string,"(F%1d,[%-4d,%-4d])-(A%-4d,%-4d) ([%-6s,%-6s],%-6s)",f1,A1,A2,j,k,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,PseudoAtoms[typeB].Name);
              energy=PrintChargeBondDipoleStatus(nr++,string,cosA,DipoleMagnitudeA,chargeB,r);

              UHostCationChargeBondDipoleReal+=energy;
            }
          }
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    TypeA=Cations[CurrentSystem][i].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      posA1=Cations[CurrentSystem][i].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][i].Atoms[A2].Position;
      typeA1=Cations[CurrentSystem][i].Atoms[A1].Type;
      typeA2=Cations[CurrentSystem][i].Atoms[A2].Type;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

            sprintf(string,"(C%-4d,[%-4d,%-4d])-(A%-4d,%-4d) ([%-6s,%-6s],%-6s)",f1,A1,A2,i,j,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,PseudoAtoms[typeB].Name);
            energy=PrintChargeBondDipoleStatus(nr++,string,cosA,DipoleMagnitudeA,chargeB,r);
            UHostCationChargeBondDipoleReal+=energy;
          }
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostCationChargeBondDipoleReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostCationChargeBondDipoleReal*ENERGY_TO_KELVIN,
          UHostCationChargeBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UHostCationChargeBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");

}

void PrintFrameworkCationBondDipoleBondDipoleEnergyStatus(void)
{
  int i,k,l,f1;
  int A1,A2,B1,B2;
  int TypeB,typeA1,typeA2,typeB1,typeB2;;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,r2,ri2,rk2,cosA,cosB,cosAB,energy,temp,length;
  VECTOR dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL UHostCationBondDipoleBondDipoleReal;
  char string[256];
  int nr;

  fprintf(OutputFilePtr[CurrentSystem],"UHostCationBondDipoleBondDipoleReal\n");
  fprintf(OutputFilePtr[CurrentSystem],"==========================================================================\n");

  nr=0;
  UHostCationBondDipoleBondDipoleReal=0.0;

  Bt0=Bt1=Bt2=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
      A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
      A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
      posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
      typeA1=Framework[CurrentSystem].Atoms[f1][A1].Type;
      typeA2=Framework[CurrentSystem].Atoms[f1][A2].Type;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        TypeB=Cations[CurrentSystem][k].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
          posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
          typeB1=Cations[CurrentSystem][k].Atoms[B1].Type;
          typeB2=Cations[CurrentSystem][k].Atoms[B2].Type;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(r2<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(r2);
            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

            sprintf(string,"(F%1d,[%-4d,%-4d])-(F%1d,[%-4d,-%4d]) ([%-6s,%-6s],[%-6s,%-6s])",i,A1,A2,k,B1,B2,PseudoAtoms[typeA1].Name,PseudoAtoms[typeA2].Name,
                             PseudoAtoms[typeB1].Name,PseudoAtoms[typeB2].Name);
            energy=PrintBondDipoleBondDipoleStatus(nr++,string,cosA,cosB,cosAB,DipoleMagnitudeA,DipoleMagnitudeB,r);
            UHostCationBondDipoleBondDipoleReal+=energy;
          }
        }
      }
    }
  }
  fprintf(OutputFilePtr[CurrentSystem],"\n");
  fprintf(OutputFilePtr[CurrentSystem],"Total interactions within cutoff-range: %d\n",nr);
  fprintf(OutputFilePtr[CurrentSystem],"Total UHostCationBondDipoleBondDipoleReal: %10.6f [K] %10.6f [kJ/mol] %10.6f [kcal/mol]\n",
          UHostCationBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,
          UHostCationBondDipoleBondDipoleReal*ENERGY_TO_KJ_PER_MOL,
          UHostCationBondDipoleBondDipoleReal*ENERGY_TO_KCAL_PER_MOL);
  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
}

void Status(void)
{

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    if(PrintFrameworkBondStatus)
      PrintFrameworkBondEnergyStatus();
    if(PrintFrameworkUreyBradleyStatus)
      PrintFrameworkUreyBradleyEnergyStatus();
    if(PrintFrameworkBendStatus)
      PrintFrameworkBendEnergyStatus();
    if(PrintFrameworkInversionBendStatus)
      PrintFrameworkInversionBendEnergyStatus();
    if(PrintFrameworkTorsionStatus)
      PrintFrameworkTorsionEnergyStatus();
    if(PrintFrameworkImproperTorsionStatus)
      PrintFrameworkImproperTorsionEnergyStatus();
    if(PrintFrameworkBondBondStatus)
      PrintFrameworkBondBondEnergyStatus();
    if(PrintFrameworkBondBendStatus)
      PrintFrameworkBondBendEnergyStatus();
    if(PrintFrameworkBendBendStatus)
      PrintFrameworkBendBendEnergyStatus();
    if(PrintFrameworkBendTorsionStatus)
      PrintFrameworkBendTorsionEnergyStatus();
    if(PrintFrameworkIntraVDWStatus)
      PrintFrameworkIntraVDWEnergyStatus();
    if(PrintFrameworkIntraChargeChargeStatus)
      PrintFrameworkIntraChargeChargeEnergyStatus();
    if(PrintFrameworkIntraChargeBondDipoleStatus)
      PrintFrameworkIntraChargeBondDipoleEnergyStatus();
    if(PrintFrameworkIntraBondDipoleBondDipoleStatus)
      PrintFrameworkIntraBondDipoleBondDipoleEnergyStatus();

    if(PrintAdsorbateBondStatus)
      PrintAdsorbateBondEnergyStatus();
    if(PrintAdsorbateUreyBradleyStatus)
      PrintAdsorbateUreyBradleyEnergyStatus();
    if(PrintAdsorbateBendStatus)
      PrintAdsorbateBendEnergyStatus();
    if(PrintAdsorbateTorsionStatus)
      PrintAdsorbateTorsionEnergyStatus();
    if(PrintAdsorbateImproperTorsionStatus)
      PrintAdsorbateImproperTorsionEnergyStatus();
    if(PrintAdsorbateBondBondStatus)
      PrintAdsorbateBondBondEnergyStatus();
    if(PrintAdsorbateBondBendStatus)
      PrintAdsorbateBondBendEnergyStatus();
    if(PrintAdsorbateBendBendStatus)
      PrintAdsorbateBendBendEnergyStatus();
    if(PrintAdsorbateBondTorsionStatus)
      PrintAdsorbateBondTorsionEnergyStatus();
    if(PrintAdsorbateBendTorsionStatus)
      PrintAdsorbateBendTorsionEnergyStatus();
    if(PrintAdsorbateIntraVDWStatus)
      PrintAdsorbateIntraVDWEnergyStatus();
    if(PrintAdsorbateIntraChargeChargeStatus)
      PrintAdsorbateIntraChargeChargeEnergyStatus();
    if(PrintAdsorbateIntraChargeBondDipoleStatus)
      PrintAdsorbateIntraChargeBondDipoleEnergyStatus();
    if(PrintAdsorbateIntraBondDipoleBondDipoleStatus)
      PrintAdsorbateIntraBondDipoleBondDipoleEnergyStatus();

    if(PrintCationBondStatus)
      PrintCationBondEnergyStatus();
    if(PrintCationUreyBradleyStatus)
      PrintCationUreyBradleyEnergyStatus();
    if(PrintCationBendStatus)
      PrintCationBendEnergyStatus();
    if(PrintCationTorsionStatus)
      PrintCationTorsionEnergyStatus();
    if(PrintCationImproperTorsionStatus)
      PrintCationImproperTorsionEnergyStatus();
    if(PrintCationBondBondStatus)
      PrintCationBondBondEnergyStatus();
    if(PrintCationBondBendStatus)
      PrintCationBondBendEnergyStatus();
    if(PrintCationBendBendStatus)
      PrintCationBendBendEnergyStatus();
    if(PrintCationBondTorsionStatus)
      PrintCationBondTorsionEnergyStatus();
    if(PrintCationBendTorsionStatus)
      PrintCationBendTorsionEnergyStatus();
    if(PrintCationIntraVDWStatus)
      PrintCationIntraVDWEnergyStatus();
    if(PrintCationIntraChargeChargeStatus)
      PrintCationIntraChargeChargeEnergyStatus();
    if(PrintCationIntraChargeBondDipoleStatus)
      PrintCationIntraChargeBondDipoleEnergyStatus();
    if(PrintCationIntraBondDipoleBondDipoleStatus)
      PrintCationIntraBondDipoleBondDipoleEnergyStatus();

    if(PrintInterVDWStatus)
      PrintInterVDWEnergyStatus();
    if(PrintInterChargeChargeStatus)
      PrintInterChargeChargeEnergyStatus();
    if(PrintInterChargeBondDipoleStatus)
      PrintInterChargeBondDipoleEnergyStatus();
    if(PrintInterBondDipoleBondDipoleStatus)
      PrintInterBondDipoleBondDipoleEnergyStatus();

    if(PrintFrameworkAdsorbateVDWStatus)
      PrintFrameworkAdsorbateVDWEnergyStatus();
    if(PrintFrameworkAdsorbateChargeChargeStatus)
      PrintFrameworkAdsorbateChargeChargeEnergyStatus();
    if(PrintFrameworkAdsorbateChargeBondDipoleStatus)
      PrintFrameworkAdsorbateChargeBondDipoleEnergyStatus();
    if(PrintFrameworkAdsorbateBondDipoleBondDipoleStatus)
      PrintFrameworkAdsorbateBondDipoleBondDipoleEnergyStatus();

    if(PrintFrameworkCationVDWStatus)
      PrintFrameworkCationVDWEnergyStatus();
    if(PrintFrameworkCationChargeChargeStatus)
      PrintFrameworkCationChargeChargeEnergyStatus();
    if(PrintFrameworkCationChargeBondDipoleStatus)
      PrintFrameworkCationChargeBondDipoleEnergyStatus();
    if(PrintFrameworkCationBondDipoleBondDipoleStatus)
      PrintFrameworkCationBondDipoleBondDipoleEnergyStatus();
  }
}


