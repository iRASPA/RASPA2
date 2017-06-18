/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_born.c' is part of RASPA-2.0

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
#include "framework_force.h"
#include "simulation.h"
#include "ewald.h"
#include "potentials.h"
#include "utils.h"
#include "output.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "spacegroup.h"
#include "spectra.h"

void CalculateFrameworkBondBornTerm(void)
{
  int i;   // loop variable
  REAL r;  // distance
  REAL rr; // distance squared
  REAL temp,temp2; // temporary
  REAL exp_term;   // temporary
  REAL U;  // energy of a specific interaction
  REAL DF; // first derivative
  REAL DDF;  // second derivative
  VECTOR dr,f; // atoms separation vector
  int A,B;   // atom-indices
  REAL *parms;  // pointer to potential parameter

  for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBonds[CurrentFramework];i++)
    {
      A=Framework[CurrentSystem].Bonds[CurrentFramework][i].A;
      B=Framework[CurrentSystem].Bonds[CurrentFramework][i].B;

      dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
           Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
      dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
           Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
      dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
           Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;

      // apply boundary condition
      dr=ApplyBoundaryCondition(dr);

      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      parms=(REAL*)&Framework[CurrentSystem].BondArguments[CurrentFramework][i];

      switch(Framework[CurrentSystem].BondType[CurrentFramework][i])
      {
        case HARMONIC_BOND:
          // 0.5*p0*SQR(r-p1);
          // ===============================================
          // p_0/k_B [K/A^2]   force constant
          // p_1     [A]       reference bond distance
          U=0.5*parms[0]*SQR(r-parms[1]);
          DF=parms[0]*(r-parms[1])/r;
          DDF=(parms[0]*parms[1])/(r*rr);
          break;
       case CORE_SHELL_SPRING:
          U=0.5*parms[0]*SQR(r);
          DF=parms[0];
          DDF=0.0;
          break;
        case MORSE_BOND:
          // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
          // ===============================================
          // p_0/k_B [K]       force constant
          // p_1     [A^-1]    parameter
          // p_2     [A]       reference bond distance
          temp=exp(parms[1]*(parms[2]-r));
          U=parms[0]*(SQR(1.0-temp)-1.0);
          DF=2.0*parms[0]*parms[1]*(1.0-temp)*temp/r;
          DDF=2.0*parms[0]*parms[1]*temp*((1.0+2.0*parms[1]*r)*temp-parms[1]*r-1.0)/(r*rr);
          break;
        case LJ_12_6_BOND:
          // A/r_ij^12-B/r_ij^6
          // ===============================================
          // p_0/k_B [K A^12]
          // p_1/k_B [K A^6]
          temp=CUBE(1.0/rr);
          U=parms[0]*SQR(temp)-parms[1]*temp;
          DF=6.0*(parms[1]*temp-2.0*parms[0]*SQR(temp))/rr;
          DDF=24.0*(7.0*parms[0]*SQR(temp)-2.0*parms[1]*temp)/SQR(rr);
          break;
        case LENNARD_JONES_BOND:
          // 4*p_0*((p_1/r)^12-(p_1/r)^6)
          // ===============================================
          // p_0/k_B [K]
          // p_1     [A]
          temp=CUBE(parms[1]/rr);
          U=4.0*parms[0]*(temp*(temp-1.0));
          DF=24.0*parms[0]*(temp*(1.0-2.0*temp))/rr;
          DDF=96.0*parms[0]*(temp*(7.0*temp-2.0))/SQR(rr);
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
          DF=(6.0/rr)*temp-parms[1]*exp_term/r;
          DDF=(-48.0*temp/rr+parms[1]*(1.0+parms[1]*r)*exp_term/r)/rr;
          break;
        case RESTRAINED_HARMONIC_BOND:
          // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
          // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
          // ===============================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2     [A]
          temp=r-parms[1];
          U=0.5*parms[0]*SQR(MIN2(fabs(temp),parms[2]))
                +parms[0]*parms[2]*MAX2(fabs(temp)-parms[2],(REAL)0.0);
          DF=-parms[0]*(SIGN(MIN2(fabs(temp),parms[2]),temp))/r;
          DDF=fabs(temp)<parms[2]?-parms[0]*parms[1]/(r*rr):parms[0]*SIGN(parms[2],temp)/(r*rr);
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
          DF=temp*(parms[0]+parms[2]*temp+parms[3]*temp2)/r;
          DDF=2.0*parms[3]+(parms[2]-3.0*parms[1]*parms[3])/r+(parms[1]*(parms[0]+parms[1]*(parms[1]*parms[3]-parms[2])))/(r*rr);
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
          DF=temp*(2.0*parms[0]+3.0*parms[2]*temp+4.0*parms[3]*temp2)/r;
          DDF=8.0*parms[3]+(3.0*parms[2]-3.0*parms[1]*4.0*parms[3])/r+(parms[1]*(2.0*parms[0]+parms[1]*(parms[1]*4.0*parms[3]-3.0*parms[2])))/(r*rr);
          break;
        case MM3_BOND:
          // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
          // =================================================================
          // p_0     [mdyne/A molecule]
          // p_1     [A]
          temp=r-parms[1];
          temp2=SQR(temp);
          U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
          DF=parms[0]*(2.0+2.55*(4.0*2.55*(7.0/12.0)*temp-3.0)*temp)*temp/r;
          DDF=(parms[0]*(SQR(2.55)*4.0*7.0*temp2*(parms[1]+2.0*r)+12.0*(2.0*parms[1]+2.55*3.0*(SQR(parms[1])-SQR(r)))))/(12.0*SQR(r)*r);
          break;
        case RIGID_BOND:
          U=DF=DDF=0.0;
          break;
        case FIXED_BOND:
          U=DF=DDF=0.0;
          break;
        case MEASURE_BOND:
          U=DF=DDF=0.0;
          break;
        default:
          U=DF=DDF=0.0;
          fprintf(stderr, "Undefined Bond potential (Bond born framework)\n");
          exit(0);
          break;
      }

      // add contribution to the Adsorbate stretch energy
      UHostBond[CurrentSystem]+=U;

      // forces are oppositely directed to the gradient
      f.x=-DF*dr.x;
      f.y=-DF*dr.y;
      f.z=-DF*dr.z;

      // add contribution to the forces
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.x+=f.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.y+=f.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.z+=f.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.x-=f.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.y-=f.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.z-=f.z;

      // add contribution to the strain derivative tensor
      StrainDerivativeTensor[CurrentSystem].ax-=f.x*dr.x;
      StrainDerivativeTensor[CurrentSystem].bx-=f.y*dr.x;
      StrainDerivativeTensor[CurrentSystem].cx-=f.z*dr.x;

      StrainDerivativeTensor[CurrentSystem].ay-=f.x*dr.y;
      StrainDerivativeTensor[CurrentSystem].by-=f.y*dr.y;
      StrainDerivativeTensor[CurrentSystem].cy-=f.z*dr.y;

      StrainDerivativeTensor[CurrentSystem].az-=f.x*dr.z;
      StrainDerivativeTensor[CurrentSystem].bz-=f.y*dr.z;
      StrainDerivativeTensor[CurrentSystem].cz-=f.z*dr.z;

      // add contribution to the born term
      AddContributionToBornTerm(DDF,DF,dr);
    }
  }
}

void CalculateFrameworkUreyBradleyBornTerm(void)
{
  int i;   // loop variable
  REAL r;  // distance
  REAL rr;  // distance squared
  REAL temp,temp2; // temporary
  REAL exp_term; // temporary
  REAL U;  // energy of a specific interaction
  REAL DF; // first derivative
  REAL DDF;  // second derivative
  VECTOR dr,f; // atoms separation vector
  int A,C;   // atom-indices
  REAL *parms;  // pointer to potential parameter

  UHostUreyBradley[CurrentSystem]=0.0;
  for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfUreyBradleys[CurrentFramework];i++)
    {
      A=Framework[CurrentSystem].UreyBradleys[CurrentFramework][i].A;
      C=Framework[CurrentSystem].UreyBradleys[CurrentFramework][i].C;

      dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.x;
      dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.y;
      dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      parms=(REAL*)&Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][i];

      switch(Framework[CurrentSystem].UreyBradleyType[CurrentFramework][i])
      {
        case HARMONIC_UREYBRADLEY:
          // 0.5*p0*SQR(r-p1);
          // ===============================================
          // p_0/k_B [K/A^2]   force constant
          // p_1     [A]       reference bond distance
          U=0.5*parms[0]*SQR(r-parms[1]);
          DF=parms[0]*(r-parms[1])/r;
          DDF=(parms[0]*parms[1])/(r*rr);
          break;
        case MORSE_UREYBRADLEY:
          // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
          // ===============================================
          // p_0/k_B [K]       force constant
          // p_1     [A^-1]    parameter
          // p_2     [A]       reference bond distance
          temp=exp(parms[1]*(parms[2]-r));
          U=parms[0]*(SQR(1.0-temp)-1.0);
          DF=2.0*parms[0]*parms[1]*(1.0-temp)*temp/r;
          DDF=2.0*parms[0]*parms[1]*temp*((1.0+2.0*parms[1]*r)*temp-parms[1]*r-1.0)/(r*rr);
          break;
        case LJ_12_6_UREYBRADLEY:
          // A/r_ij^12-B/r_ij^6
          // ===============================================
          // p_0/k_B [K A^12]
          // p_1/k_B [K A^6]
          temp=CUBE(1.0/rr);
          U=parms[0]*SQR(temp)-parms[1]*temp;
          DF=6.0*(parms[1]*temp-2.0*parms[0]*SQR(temp))/rr;
          DDF=24.0*(7.0*parms[0]*SQR(temp)-2.0*parms[1]*temp)/SQR(rr);
          break;
        case LENNARD_JONES_UREYBRADLEY:
          // 4*p_0*((p_1/r)^12-(p_1/r)^6)
          // ===============================================
          // p_0/k_B [K]
          // p_1     [A]
          temp=CUBE(parms[1]/rr);
          U=4.0*parms[0]*(temp*(temp-1.0));
          DF=24.0*parms[0]*(temp*(1.0-2.0*temp))/rr;
          DDF=96.0*parms[0]*(temp*(7.0*temp-2.0))/SQR(rr);
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
          DF=(6.0/rr)*temp-parms[1]*exp_term/r;
          DDF=(-48.0*temp/rr+parms[1]*(1.0+parms[1]*r)*exp_term/r)/rr;
          break;
        case RESTRAINED_HARMONIC_UREYBRADLEY:
          // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
          // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
          // ===============================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2     [A]
          temp=r-parms[1];
          U=0.5*parms[0]*SQR(MIN2(fabs(temp),parms[2]))
                +parms[0]*parms[2]*MAX2(fabs(temp)-parms[2],(REAL)0.0);
          DF=-parms[0]*(SIGN(MIN2(fabs(temp),parms[2]),temp))/r;
          DDF=fabs(temp)<parms[2]?-parms[0]*parms[1]/(r*rr):parms[0]*SIGN(parms[2],temp)/(r*rr);
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
          DF=temp*(parms[0]+parms[2]*temp+parms[3]*temp2)/r;
          DDF=2.0*parms[3]+(parms[2]-3.0*parms[1]*parms[3])/r+(parms[1]*(parms[0]+parms[1]*(parms[1]*parms[3]-parms[2])))/(r*rr);
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
          DF=temp*(2.0*parms[0]+3.0*parms[2]*temp+4.0*parms[3]*temp2)/r;
          DDF=8.0*parms[3]+(3.0*parms[2]-3.0*parms[1]*4.0*parms[3])/r+(parms[1]*(2.0*parms[0]+parms[1]*(parms[1]*4.0*parms[3]-3.0*parms[2])))/(r*rr);
          break;
        case MM3_UREYBRADLEY:
          // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
          // =================================================================
          // p_0     [mdyne/A molecule]
          // p_1     [A]
          temp=r-parms[1];
          temp2=SQR(temp);
          U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
          DF=parms[0]*(2.0+2.55*(4.0*2.55*(7.0/12.0)*temp-3.0)*temp)*temp/r;
          DDF=(parms[0]*(SQR(2.55)*4.0*7.0*temp2*(parms[1]+2.0*r)+12.0*(2.0*parms[1]+2.55*3.0*(SQR(parms[1])-SQR(r)))))/(12.0*SQR(r)*r);
          break;
        case RIGID_UREYBRADLEY:
          U=DF=DDF=0.0;
          break;
        case FIXED_UREYBRADLEY:
          U=DF=DDF=0.0;
          break;
        case MEASURE_UREYBRADLEY:
          U=DF=DDF=0.0;
          break;
        default:
          U=DF=DDF=0.0;
          fprintf(stderr, "Undefined Urey-Bradley potential\n");
          exit(0);
          break;
      }

      // add contribution to the energy
      UHostUreyBradley[CurrentSystem]+=U;

      // forces are oppositely directed to the gradient
      f.x=-DF*dr.x;
      f.y=-DF*dr.y;
      f.z=-DF*dr.z;

      // add contribution to the forces
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.x+=f.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.y+=f.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.z+=f.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.x-=f.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.y-=f.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.z-=f.z;

      // add contribution to the strain derivative tensor
      StrainDerivativeTensor[CurrentSystem].ax-=f.x*dr.x;
      StrainDerivativeTensor[CurrentSystem].bx-=f.y*dr.x;
      StrainDerivativeTensor[CurrentSystem].cx-=f.z*dr.x;

      StrainDerivativeTensor[CurrentSystem].ay-=f.x*dr.y;
      StrainDerivativeTensor[CurrentSystem].by-=f.y*dr.y;
      StrainDerivativeTensor[CurrentSystem].cy-=f.z*dr.y;

      StrainDerivativeTensor[CurrentSystem].az-=f.x*dr.z;
      StrainDerivativeTensor[CurrentSystem].bz-=f.y*dr.z;
      StrainDerivativeTensor[CurrentSystem].cz-=f.z*dr.z;

      // add contribution to the born term
      AddContributionToBornTerm(DDF,DF,dr);
    }
  }
}

void CalculateFrameworkBendBornTerm(void)
{
  int i,A,B,C;
  REAL *parms,U;
  REAL CosTheta,Theta,SinTheta,temp,temp2;
  REAL rab,rbc,rac;
  POINT posA,posB,posC;
  VECTOR Rab,Rbc,Rac,fa,fb,fc;
  REAL DTDX,DF,DDF,DF2;
  VECTOR dtA,dtB,dtC;
  REAL_MATRIX3x3 S,T;
  VECTOR vec_u,vec_v;
  REAL u,v;

  for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBends[CurrentFramework];i++)
    {
      A=Framework[CurrentSystem].Bends[CurrentFramework][i].A;
      B=Framework[CurrentSystem].Bends[CurrentFramework][i].B;
      C=Framework[CurrentSystem].Bends[CurrentFramework][i].C;
      parms=Framework[CurrentSystem].BendArguments[CurrentFramework][i];

      posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
      posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;
      posC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      Rab=ApplyBoundaryCondition(Rab);
      vec_u=Rab;
      rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
      u=rab;
      Rab.x/=rab;
      Rab.y/=rab;
      Rab.z/=rab;

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      Rbc=ApplyBoundaryCondition(Rbc);
      vec_v=Rbc;
      rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
      v=rbc;
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
      Theta=acos(CosTheta);
      SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));
      DTDX=-1.0/SinTheta;

      switch(Framework[CurrentSystem].BendType[CurrentFramework][i])
      {
        case HARMONIC_BEND:
          // (1/2)p_0*(theta-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          temp=Theta-parms[1];
          U=0.5*parms[0]*SQR(temp);
          DF=parms[0]*temp*DTDX;
          DF2=-parms[0]*temp*DTDX;
          DDF=parms[0]*SQR(DTDX)+parms[0]*temp*CosTheta*CUBE(DTDX);
          break;
        case CORE_SHELL_BEND:
          // (1/2)p_0*(theta-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          temp=Theta-parms[1];
          U=0.5*parms[0]*SQR(temp);
          DF=parms[0]*temp*DTDX;
          DF2=-parms[0]*temp*DTDX;
          DDF=parms[0]*SQR(DTDX)+parms[0]*temp*CosTheta*CUBE(DTDX);
          break;
        case QUARTIC_BEND:
          // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
          // ======================================================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // p_2/k_B [K/rad^3]
          // p_3/k_B [K/rad^4]
          temp=(Theta-parms[1]);
          temp2=SQR(temp);
          U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
          DF=(parms[0]*temp+parms[2]*temp2+parms[3]*temp*temp2)*DTDX;
          DDF=parms[0]*SQR(DTDX)+parms[0]*temp*CosTheta*CUBE(DTDX)
              +2.0*parms[2]*temp*SQR(DTDX)+parms[2]*temp2*CosTheta*CUBE(DTDX)
              +3.0*parms[3]*temp2*SQR(DTDX)+parms[3]*temp*temp2*CosTheta*CUBE(DTDX);
          break;
        case CFF_QUARTIC_BEND:
          // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
          // =====================================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // p_2/k_B [K/rad^3]
          // p_3/k_B [K/rad^4]
          temp=(Theta-parms[1]);
          temp2=SQR(temp);
          U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
          DF=(2.0*parms[0]*temp+3.0*parms[2]*temp2+4.0*parms[3]*temp*temp2)*DTDX;
          DDF=2.0*parms[0]*SQR(DTDX)+2.0*parms[0]*temp*CosTheta*CUBE(DTDX)
              +6.0*parms[2]*temp*SQR(DTDX)+3.0*parms[2]*temp2*CosTheta*CUBE(DTDX)
              +12.0*parms[3]*temp2*SQR(DTDX)+4.0*parms[3]*temp*temp2*CosTheta*CUBE(DTDX);
          break;
        case HARMONIC_COSINE_BEND:
          // (1/2)*p_0*(cos(theta)-cos(p_1))^2
          // ===============================================
          // p_0/k_B [K]
          // p_1     [degrees]
          U=0.5*parms[0]*SQR(CosTheta-parms[1]);
          DF=parms[0]*(CosTheta-parms[1]);
          DF2=-parms[0]*(CosTheta-parms[1]);
          DDF=parms[0];
          break;
        case COSINE_BEND:
          // p_0*(1+cos(p_1*theta-p_2))
          // ===============================================
          // p_0/k_B [K]
          // p_1     [-]
          // p_2     [degrees]
          temp=parms[1]*Theta-parms[2];
          U=parms[0]*(1.0+cos(temp));
          DF=-(parms[0]*parms[1]*sin(temp))*DTDX;
          DDF=-parms[0]*parms[1]*(parms[1]*cos(temp)+sin(temp)*CosTheta*DTDX)*SQR(DTDX);
          break;
        case TAFIPOLSKY_BEND:
          // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
          // ===============================================
          // p_0/k_B [K]
          U=0.5*parms[0]*(1.0+cos(Theta))*(1.0+cos(2.0*Theta));
          DF=parms[0]*CosTheta*(2.0+3.0*CosTheta);
          DDF=2.0*parms[0]*(1.0+3.0*CosTheta);
          break;
        case MM3_BEND:
        case MM3_IN_PLANE_BEND:
          // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
          // =================================================================================================
          // p_0/k_B [mdyne A/rad^2]
          // p_1     [degrees]
          temp=RAD2DEG*(Theta-parms[1]);
          temp2=SQR(temp);
          U=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
          DF=parms[0]*RAD2DEG*(2.0-(3.0*0.014-(4.0*5.6e-5-(5.0*7.0e-7-6.0*2.2e-8*temp)*temp)*temp)*temp)*temp*DTDX;
          fprintf(stderr, "TO BE DONE!\n");
          exit(0);
          break;
        case FIXED_BEND:
          U=DF=DDF=0.0;
          break;
        case MEASURE_BEND:
          U=DF=DDF=0.0;
          break;
        default:
          U=DF=DDF=0.0;
          fprintf(stderr, "Undefined Bend potential in routine 'CalculateFrameworkBendBornTerm'\n");
          exit(0);
          break;
      }

      // add contribution to the energy
      UHostBend[CurrentSystem]+=U;

      // Calculate the components of the derivatives.
      dtA.x=(Rbc.x-CosTheta*Rab.x)/rab;
      dtA.y=(Rbc.y-CosTheta*Rab.y)/rab;
      dtA.z=(Rbc.z-CosTheta*Rab.z)/rab;

      dtC.x=(Rab.x-CosTheta*Rbc.x)/rbc;
      dtC.y=(Rab.y-CosTheta*Rbc.y)/rbc;
      dtC.z=(Rab.z-CosTheta*Rbc.z)/rbc;

      dtB.x=-(dtA.x+dtC.x);
      dtB.y=-(dtA.y+dtC.y);
      dtB.z=-(dtA.z+dtC.z);

      // forces are oppositely directed to the gradient
      fa.x=-DF*dtA.x;
      fa.y=-DF*dtA.y;
      fa.z=-DF*dtA.z;

      fb.x=-DF*dtB.x;
      fb.y=-DF*dtB.y;
      fb.z=-DF*dtB.z;

      fc.x=-DF*dtC.x;
      fc.y=-DF*dtC.y;
      fc.z=-DF*dtC.z;

      // add contribution to the forces
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.x+=fa.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.y+=fa.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.z+=fa.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.x-=(fa.x+fc.x);
      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.y-=(fa.y+fc.y);
      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.z-=(fa.z+fc.z);

      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.x+=fc.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.y+=fc.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.z+=fc.z;

      // add contribution to the strain derivative tensor
      // Note: rab and rbc are here because the vectors were normalized before
      S.ax=rab*Rab.x*DF*dtA.x+rbc*Rbc.x*DF*dtC.x;
      S.bx=rab*Rab.y*DF*dtA.x+rbc*Rbc.y*DF*dtC.x;
      S.cx=rab*Rab.z*DF*dtA.x+rbc*Rbc.z*DF*dtC.x;

      S.ay=rab*Rab.x*DF*dtA.y+rbc*Rbc.x*DF*dtC.y;
      S.by=rab*Rab.y*DF*dtA.y+rbc*Rbc.y*DF*dtC.y;
      S.cy=rab*Rab.z*DF*dtA.y+rbc*Rbc.z*DF*dtC.y;

      S.az=rab*Rab.x*DF*dtA.z+rbc*Rbc.x*DF*dtC.z;
      S.bz=rab*Rab.y*DF*dtA.z+rbc*Rbc.y*DF*dtC.z;
      S.cz=rab*Rab.z*DF*dtA.z+rbc*Rbc.z*DF*dtC.z;

      StrainDerivativeTensor[CurrentSystem].ax+=S.ax;
      StrainDerivativeTensor[CurrentSystem].bx+=S.bx;
      StrainDerivativeTensor[CurrentSystem].cx+=S.cx;

      StrainDerivativeTensor[CurrentSystem].ay+=S.ay;
      StrainDerivativeTensor[CurrentSystem].by+=S.by;
      StrainDerivativeTensor[CurrentSystem].cy+=S.cy;

      StrainDerivativeTensor[CurrentSystem].az+=S.az;
      StrainDerivativeTensor[CurrentSystem].bz+=S.bz;
      StrainDerivativeTensor[CurrentSystem].cz+=S.cz;

      T.ax=(vec_u.x*vec_v.x+vec_v.x*vec_u.x)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.x*vec_u.x)/SQR(u)+(vec_v.x*vec_v.x)/SQR(v));
      T.ay=(vec_u.x*vec_v.y+vec_v.x*vec_u.y)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.x*vec_u.y)/SQR(u)+(vec_v.x*vec_v.y)/SQR(v));
      T.az=(vec_u.x*vec_v.z+vec_v.x*vec_u.z)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.x*vec_u.z)/SQR(u)+(vec_v.x*vec_v.z)/SQR(v));

      T.bx=(vec_u.y*vec_v.x+vec_v.y*vec_u.x)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.y*vec_u.x)/SQR(u)+(vec_v.y*vec_v.x)/SQR(v));
      T.by=(vec_u.y*vec_v.y+vec_v.y*vec_u.y)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.y*vec_u.y)/SQR(u)+(vec_v.y*vec_v.y)/SQR(v));
      T.bz=(vec_u.y*vec_v.z+vec_v.y*vec_u.z)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.y*vec_u.z)/SQR(u)+(vec_v.y*vec_v.z)/SQR(v));

      T.cx=(vec_u.z*vec_v.x+vec_v.z*vec_u.x)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.z*vec_u.x)/SQR(u)+(vec_v.z*vec_v.x)/SQR(v));
      T.cy=(vec_u.z*vec_v.y+vec_v.z*vec_u.y)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.z*vec_u.y)/SQR(u)+(vec_v.z*vec_v.y)/SQR(v));
      T.cz=(vec_u.z*vec_v.z+vec_v.z*vec_u.z)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.z*vec_u.z)/SQR(u)+(vec_v.z*vec_v.z)/SQR(v));

      BornTerm[CurrentSystem].xxxx+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                                    +T.ax*T.ax)+2.0*S.ax;
      BornTerm[CurrentSystem].xxyy+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                                    +T.ax*T.by);
      BornTerm[CurrentSystem].xxzz+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                                    +T.ax*T.cz);
      BornTerm[CurrentSystem].xxyz+=DDF*CosTheta*T.ax*CosTheta*T.bz+DF*CosTheta*(
                                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.z)/SQR(SQR(v)))
                                    +T.ax*T.bz);
      BornTerm[CurrentSystem].xxzx+=DDF*CosTheta*T.ax*CosTheta*T.cx+DF*CosTheta*(
                                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.x)/SQR(SQR(v)))
                                    +T.ax*T.cx)+S.az;
      BornTerm[CurrentSystem].xxxy+=DDF*CosTheta*T.ax*CosTheta*T.ay+DF*CosTheta*(
                                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.y)/SQR(SQR(v)))
                                    +T.ax*T.ay)+S.ay;

      BornTerm[CurrentSystem].yyyy+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                                    -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                                    +T.by*T.by)+2.0*S.by;
      BornTerm[CurrentSystem].yyzz+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                                    -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                                    +T.by*T.cz);
      BornTerm[CurrentSystem].yyyz+=DDF*CosTheta*T.by*CosTheta*T.bz+DF*CosTheta*(
                                    -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                                    +T.by*T.bz)+S.bz;
      BornTerm[CurrentSystem].yyzx+=DDF*CosTheta*T.by*CosTheta*T.cx+DF*CosTheta*(
                                    -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.x)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.x)/SQR(SQR(v)))
                                    +T.by*T.cx);
      BornTerm[CurrentSystem].yyxy+=DDF*CosTheta*T.by*CosTheta*T.ay+DF*CosTheta*(
                                    -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.y*vec_u.y*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.x*vec_v.y)/SQR(SQR(v)))
                                    +T.by*T.ay)+S.bx;

      BornTerm[CurrentSystem].zzzz+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                                    -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                                    +T.cz*T.cz)+2.0*S.cz;
      BornTerm[CurrentSystem].zzyz+=DDF*CosTheta*T.cz*CosTheta*T.bz+DF*CosTheta*(
                                    -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.z*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                                    +T.cz*T.bz)+S.cy;
      BornTerm[CurrentSystem].zzzx+=DDF*CosTheta*T.cz*CosTheta*T.cx+DF*CosTheta*(
                                    -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.x)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.x)/SQR(SQR(v)))
                                    +T.cz*T.cx)+S.cx;
      BornTerm[CurrentSystem].zzxy+=DDF*CosTheta*T.cz*CosTheta*T.ay+DF*CosTheta*(
                                    -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.z*vec_u.z*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.x*vec_v.y)/SQR(SQR(v)))
                                    +T.cz*T.ay);

      BornTerm[CurrentSystem].yzyz+=DDF*CosTheta*T.bz*CosTheta*T.bz+DF*CosTheta*(
                                    -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.y*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                                    +T.bz*T.bz)+0.5*(S.by+S.cz);
      BornTerm[CurrentSystem].yzzx+=DDF*CosTheta*T.bz*CosTheta*T.cx+DF*CosTheta*(
                                    -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.y*vec_u.z*vec_u.z*vec_u.x)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.z*vec_v.x)/SQR(SQR(v)))
                                    +T.bz*T.cx)+0.5*S.bx;
      BornTerm[CurrentSystem].yzxy+=DDF*CosTheta*T.bz*CosTheta*T.ay+DF*CosTheta*(
                                    -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.y*vec_u.z*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.x*vec_v.y)/SQR(SQR(v)))
                                    +T.bz*T.ay)+0.5*S.az;

      BornTerm[CurrentSystem].zxzx+=DDF*CosTheta*T.cx*CosTheta*T.cx+DF*CosTheta*(
                                    -(vec_v.z*vec_u.x+vec_u.z*vec_v.x)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.z*vec_u.x*vec_u.z*vec_u.x)/SQR(SQR(u))+(vec_v.z*vec_v.x*vec_v.z*vec_v.x)/SQR(SQR(v)))
                                    +T.cx*T.cx)+0.5*(S.cz+S.ax);
      BornTerm[CurrentSystem].zxxy+=DDF*CosTheta*T.cx*CosTheta*T.ay+DF*CosTheta*(
                                    -(vec_v.z*vec_u.x+vec_u.z*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.z*vec_u.x*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.z*vec_v.x*vec_v.x*vec_v.y)/SQR(SQR(v)))
                                    +T.cx*T.ay)+0.5*S.cy;

      BornTerm[CurrentSystem].xyxy+=DDF*CosTheta*T.ay*CosTheta*T.ay+DF*CosTheta*(
                                    -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                                    2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.y)/SQR(SQR(v)))
                                    +T.ay*T.ay)+0.5*(S.ax+S.by);
    }
  }
}

void CalculateFrameworkTorsionBornTerm(void)
{
  int i,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds,fa,fb,fc,fd;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  REAL ShiftedCosPhi,ShiftedCosPhi2,ShiftedSinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL *parms;
  REAL DF,DDF;
  VECTOR vec_u,vec_v,vec_w,vec_n,vec_m;
  REAL u,w,v,m2,n2,dot_mn;
  REAL_MATRIX3x3 dme,dne,dmn,S,T;

  UHostTorsion[CurrentSystem]=0.0;
  for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[CurrentFramework];i++)
    {
      A=Framework[CurrentSystem].Torsions[CurrentFramework][i].A;
      B=Framework[CurrentSystem].Torsions[CurrentFramework][i].B;
      C=Framework[CurrentSystem].Torsions[CurrentFramework][i].C;
      D=Framework[CurrentSystem].Torsions[CurrentFramework][i].D;
      parms=Framework[CurrentSystem].TorsionArguments[CurrentFramework][i];

      posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
      posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;
      posC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position;
      posD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Position;

      vec_u.x=posB.x-posC.x;
      vec_u.y=posB.y-posC.y;
      vec_u.z=posB.z-posC.z;
      vec_u=ApplyBoundaryCondition(vec_u);
      u=sqrt(SQR(vec_u.x)+SQR(vec_u.y)+SQR(vec_u.z));

      vec_v.x=posD.x-posC.x;
      vec_v.y=posD.y-posC.y;
      vec_v.z=posD.z-posC.z;
      vec_v=ApplyBoundaryCondition(vec_v);
      v=sqrt(SQR(vec_v.x)+SQR(vec_v.y)+SQR(vec_v.z));

      vec_w.x=posA.x-posB.x;
      vec_w.y=posA.y-posB.y;
      vec_w.z=posA.z-posB.z;
      vec_w=ApplyBoundaryCondition(vec_w);
      w=sqrt(SQR(vec_w.x)+SQR(vec_w.y)+SQR(vec_w.z));

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

      switch(Framework[CurrentSystem].TorsionType[CurrentFramework][i])
      {
        case HARMONIC_DIHEDRAL:
          // (1/2)*p_0*(phi-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // potential defined in terms of 'phi' and therefore contains a singularity
          // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
          // same direction as Rbc, and negative otherwise
          Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
          Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
          Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
          Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
          Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
          Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
          sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
          Phi=SIGN(acos(CosPhi),sign);
          SinPhi=sin(Phi);
          Phi-=parms[1];
          Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
          SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
          U=0.5*parms[0]*SQR(Phi);
          DF=-parms[0]*Phi/SinPhi;
          DDF=-parms[0]*(Phi*CosPhi-SinPhi)/CUBE(SinPhi);
          break;
        case HARMONIC_COSINE_DIHEDRAL:
          // (1/2)*p_0*(cos(phi)-cos(p_1))^2
          // ===============================================
          // p_0/k_B [K]
          // p_1     [degrees]
          U=0.5*parms[0]*SQR(CosPhi-parms[1]);
          DF=parms[0]*(CosPhi-parms[1]);
          DDF=parms[0];
          break;
        case THREE_COSINE_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
          DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[1]-6.0*parms[2]*CosPhi);
          break;
        case MM3_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0     [kcal/mol]
          // p_1     [kcal/mol]
          // p_2     [kcal/mol]
          U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
          DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[1]-6.0*parms[2]*CosPhi);
          break;
        case CVFF_BLOCKED_DIHEDRAL:
          // 
          // ========================================================================
          // p_0     [rad]
          // p_1     [K]
          // p_2     [-]
          // p_3     [rad]
          // p_4     [rad]
          U=0.0;
          DF=0.0;
          DDF=0.0;
          break;
        case CFF_DIHEDRAL:
          // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
          DF=-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2);
          DDF=-4.0*(parms[1]+6.0*parms[2]*CosPhi);
          break;
        case CFF_DIHEDRAL2:
          // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
          DF=parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi);
          DDF=4.0*(parms[1]+6.0*parms[2]*CosPhi);
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
          DF=-parms[1]+2.0*parms[2]*CosPhi-3.0*parms[3]*CosPhi2+4.0*parms[4]*CosPhi2*CosPhi-5.0*parms[5]*SQR(CosPhi2);
          DDF=2.0*parms[2]-6.0*parms[3]*CosPhi+12.0*parms[4]*CosPhi2-20.0*parms[5]*CosPhi2*CosPhi;
          break;
        case TRAPPE_DIHEDRAL:
          // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
          // ==========================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
          DF=parms[1]-4.0*parms[2]*CosPhi+3.0*parms[3]*(4.0*CosPhi2-1.0);
          DDF=-4.0*(parms[2]-6.0*parms[3]*CosPhi);
          break;
        case TRAPPE_DIHEDRAL_EXTENDED:
          // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
          // =============================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          U=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
          DF=parms[1]-3.0*parms[3]+4.0*(parms[2]-4.0*parms[4])*CosPhi+12.0*parms[3]*CosPhi2+32.0*parms[4]*CosPhi2*CosPhi;
          DDF=4.0*parms[2]-16.0*parms[4]+24.0*parms[3]*CosPhi+96.0*parms[4]*CosPhi2;
          break;
        case MOD_TRAPPE_DIHEDRAL:
          /* Salvador modification: 16/08/2016
           add phase in cos function:
           p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
          */
          Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
          Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
          Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
          Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
          Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
          Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
          sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
          Phi=SIGN(acos(CosPhi),sign);
          SinPhi=sin(Phi);
          SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
          Phi-=parms[4];           // shift Phi as Phi+parms[4]
          Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
          ShiftedCosPhi=cos(Phi);
          ShiftedSinPhi=sin(Phi);
          ShiftedCosPhi2=SQR(ShiftedCosPhi);
          U=parms[0]+parms[1]+parms[3]+(parms[1]-3.0*parms[3])*ShiftedCosPhi -2.0*parms[2]*ShiftedCosPhi2 +4.0*parms[3]*ShiftedCosPhi*ShiftedCosPhi2;
          DF=((parms[1]-3.0*parms[3])*ShiftedSinPhi -4.0*parms[2]*ShiftedCosPhi*ShiftedSinPhi + 12.0*parms[3]*ShiftedCosPhi2*ShiftedSinPhi)/SinPhi;
          DDF=-(parms[1]-3.0*parms[3])*sin(parms[4])/CUBE(SinPhi)+
              4.0*parms[2]*(ShiftedCosPhi2-ShiftedCosPhi*ShiftedSinPhi*CosPhi/SinPhi-SQR(ShiftedSinPhi))/SQR(SinPhi)+
              12.0*parms[3]*ShiftedCosPhi*(-ShiftedCosPhi2+ShiftedCosPhi*ShiftedSinPhi*CosPhi/SinPhi+2.0*SQR(ShiftedSinPhi))/SQR(SinPhi);
          break;
        case CVFF_DIHEDRAL:
          // p_0*(1+cos(p_1*phi-p_2))
          // ========================
          // p_0/k_B [K]
          // p_1     [-]
          // p_2     [degrees]
          // potential defined in terms of 'phi' and therefore contains a singularity
          // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
          // same direction as Dbc, and negative otherwise
          Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
          Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
          Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
          Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
          Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
          Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
          sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
          Phi=SIGN(acos(CosPhi),sign);
          SinPhi=sin(Phi);
          U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
          if(fabs(SinPhi)<1.0e-8)
          {
            SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
            DF=parms[0]*parms[1]*cos(parms[2])*SIGN(parms[1],sin(parms[1]*Phi)*Phi)+parms[0]*parms[1]*cos(parms[1]*Phi)*sin(parms[2])/SinPhi;
          }
          else
            DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
            DDF=-parms[0]*parms[1]*(parms[1]*cos(-parms[1]*Phi+parms[2])+sin(-parms[1]*Phi+parms[2])*CosPhi/SinPhi)/SQR(SinPhi);
          break;
        case OPLS_DIHEDRAL:
          // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
          // =================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
          DF=0.5*parms[1]-2.0*parms[2]*CosPhi+1.5*parms[3]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[2]-6.0*parms[3]*CosPhi);
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
          DF=0.5*(parms[0]-3.0*parms[2]+5.0*parms[4])-2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi+
             6.0*(parms[2]-5.0*parms[4])*CosPhi2-16.0*(parms[3]-6.0*parms[5])*CosPhi2*CosPhi+40.0*parms[4]*SQR(CosPhi2)-96.0*parms[5]*CosPhi2*CUBE(CosPhi);
          DDF=-2.0*parms[1]+8.0*parms[3]-18.0*parms[5]+12.0*(parms[2]-5.0*parms[4])*CosPhi
              -48.0*(parms[3]-6.0*parms[5])*CosPhi2+160.0*parms[4]*CosPhi2*CosPhi-480.0*parms[5]*SQR(CosPhi2);
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
          DF=0.5*parms[0]+parms[2]*(6.0*CosPhi2-1.5)+parms[4]*(2.5-30.0*CosPhi2+40.0*SQR(CosPhi2))+
             CosPhi*(-2.0*parms[1]+parms[3]*(16.0*CosPhi2-8.0)+parms[5]*(18.0-96.0*CosPhi2+96.0*SQR(CosPhi2)));
          DDF=-2.0*parms[1]+12.0*parms[2]*CosPhi+parms[3]*(48.0*CosPhi2-8.0)+parms[4]*CosPhi*(160.0*CosPhi2-60.0)
              +parms[5]*(18.0-288.0*CosPhi2+480.0*SQR(CosPhi2));
          break;
        default:
          fprintf(stderr, "Undefined Torsion potential in 'framework_born.c'\n");
          U=DF=DDF=0.0;
          exit(0);
          break;
      }

      // energy
      UHostTorsion[CurrentSystem]+=U;

      // Calculate the first derivative vectors.
      d=dot_ab/rbc;
      e=dot_cd/rbc;

      dtA.x=(ds.x-CosPhi*dr.x)/r;
      dtA.y=(ds.y-CosPhi*dr.y)/r;
      dtA.z=(ds.z-CosPhi*dr.z)/r;

      dtD.x=(dr.x-CosPhi*ds.x)/s;
      dtD.y=(dr.y-CosPhi*ds.y)/s;
      dtD.z=(dr.z-CosPhi*ds.z)/s;

      dtB.x=dtA.x*(d-1.0)+e*dtD.x;
      dtB.y=dtA.y*(d-1.0)+e*dtD.y;
      dtB.z=dtA.z*(d-1.0)+e*dtD.z;

      dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
      dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
      dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

      // forces are oppositely directed to the gradient
      fa.x=-DF*dtA.x;
      fa.y=-DF*dtA.y;
      fa.z=-DF*dtA.z;

      fb.x=-DF*dtB.x;
      fb.y=-DF*dtB.y;
      fb.z=-DF*dtB.z;

      fc.x=-DF*dtC.x;
      fc.y=-DF*dtC.y;
      fc.z=-DF*dtC.z;

      fd.x=-DF*dtD.x;
      fd.y=-DF*dtD.y;
      fd.z=-DF*dtD.z;

      // add contribution to the forces
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.x+=fa.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.y+=fa.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.z+=fa.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.x+=fb.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.y+=fb.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.z+=fb.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.x+=fc.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.y+=fc.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.z+=fc.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][D].Force.x+=fd.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][D].Force.y+=fd.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][D].Force.z+=fd.z;

      // add contribution to the stress tensor
      // Note: rbc is here because the vector was normalized before
      S.ax=-(Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x);
      S.bx=-(Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x);
      S.cx=-(Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x);

      S.ay=-(Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y);
      S.by=-(Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y);
      S.cy=-(Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y);

      S.az=-(Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z);
      S.bz=-(Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z);
      S.cz=-(Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z);

      StrainDerivativeTensor[CurrentSystem].ax+=S.ax;
      StrainDerivativeTensor[CurrentSystem].bx+=S.bx;
      StrainDerivativeTensor[CurrentSystem].cx+=S.cx;

      StrainDerivativeTensor[CurrentSystem].ay+=S.bx;
      StrainDerivativeTensor[CurrentSystem].by+=S.by;
      StrainDerivativeTensor[CurrentSystem].cy+=S.cy;

      StrainDerivativeTensor[CurrentSystem].az+=S.az;
      StrainDerivativeTensor[CurrentSystem].bz+=S.bz;
      StrainDerivativeTensor[CurrentSystem].cz+=S.cz;

      dme.ax=2.0*((SQR(u))*vec_v.x*vec_v.x+(SQR(v))*vec_u.x*vec_u.x-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x));
      dme.ay=2.0*((SQR(u))*vec_v.x*vec_v.y+(SQR(v))*vec_u.x*vec_u.y-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y));
      dme.az=2.0*((SQR(u))*vec_v.x*vec_v.z+(SQR(v))*vec_u.x*vec_u.z-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z));
      dme.bx=2.0*((SQR(u))*vec_v.y*vec_v.x+(SQR(v))*vec_u.y*vec_u.x-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.y*vec_u.x+vec_u.y*vec_v.x));
      dme.by=2.0*((SQR(u))*vec_v.y*vec_v.y+(SQR(v))*vec_u.y*vec_u.y-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y));
      dme.bz=2.0*((SQR(u))*vec_v.y*vec_v.z+(SQR(v))*vec_u.y*vec_u.z-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z));
      dme.cx=2.0*((SQR(u))*vec_v.z*vec_v.x+(SQR(v))*vec_u.z*vec_u.x-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x));
      dme.cy=2.0*((SQR(u))*vec_v.z*vec_v.y+(SQR(v))*vec_u.z*vec_u.y-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.z*vec_u.y+vec_u.z*vec_v.y));
      dme.cz=2.0*((SQR(u))*vec_v.z*vec_v.z+(SQR(v))*vec_u.z*vec_u.z-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z));

      dne.ax=2.0*((SQR(u))*vec_w.x*vec_w.x+(SQR(w))*vec_u.x*vec_u.x-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x));
      dne.ay=2.0*((SQR(u))*vec_w.x*vec_w.y+(SQR(w))*vec_u.x*vec_u.y-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y));
      dne.az=2.0*((SQR(u))*vec_w.x*vec_w.z+(SQR(w))*vec_u.x*vec_u.z-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z));
      dne.bx=2.0*((SQR(u))*vec_w.y*vec_w.x+(SQR(w))*vec_u.y*vec_u.x-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.y*vec_u.x+vec_u.y*vec_w.x));
      dne.by=2.0*((SQR(u))*vec_w.y*vec_w.y+(SQR(w))*vec_u.y*vec_u.y-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y));
      dne.bz=2.0*((SQR(u))*vec_w.y*vec_w.z+(SQR(w))*vec_u.y*vec_u.z-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z));
      dne.cx=2.0*((SQR(u))*vec_w.z*vec_w.x+(SQR(w))*vec_u.z*vec_u.x-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x));
      dne.cy=2.0*((SQR(u))*vec_w.z*vec_w.y+(SQR(w))*vec_u.z*vec_u.y-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.z*vec_u.y+vec_u.z*vec_w.y));
      dne.cz=2.0*((SQR(u))*vec_w.z*vec_w.z+(SQR(w))*vec_u.z*vec_u.z-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z));

      dmn.ax=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.x*vec_w.x+vec_w.x*vec_u.x)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x);
      dmn.ay=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.x*vec_w.y+vec_w.x*vec_u.y)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y);
      dmn.az=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.x*vec_w.z+vec_w.x*vec_u.z)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z);

      dmn.bx=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.y*vec_w.x+vec_w.y*vec_u.x)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.y*vec_u.x+vec_u.y*vec_v.x)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.y*vec_u.x+vec_u.y*vec_u.x)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.y*vec_w.x+vec_w.y*vec_v.x);
      dmn.by=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.y*vec_w.y+vec_w.y*vec_u.y)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y);
      dmn.bz=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.y*vec_w.z+vec_w.y*vec_u.z)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z);

      dmn.cx=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.z*vec_w.x+vec_w.z*vec_u.x)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x);
      dmn.cy=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.z*vec_w.y+vec_w.z*vec_u.y)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.z*vec_u.y+vec_u.z*vec_v.y)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.z*vec_u.y+vec_u.z*vec_u.y)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.z*vec_w.y+vec_w.z*vec_v.y);
      dmn.cz=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.z*vec_w.z+vec_w.z*vec_u.z)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z);

      vec_m.x=vec_v.y*vec_u.z-vec_v.z*vec_u.y;
      vec_m.y=vec_v.z*vec_u.x-vec_v.x*vec_u.z;
      vec_m.z=vec_v.x*vec_u.y-vec_v.y*vec_u.x;
      m2=SQR(vec_m.x)+SQR(vec_m.y)+SQR(vec_m.z);

      vec_n.x=vec_u.y*vec_w.z-vec_u.z*vec_w.y;
      vec_n.y=vec_u.z*vec_w.x-vec_u.x*vec_w.z;
      vec_n.z=vec_u.x*vec_w.y-vec_u.y*vec_w.x;
      n2=SQR(vec_n.x)+SQR(vec_n.y)+SQR(vec_n.z);

      dot_mn=(vec_m.x*vec_n.x+vec_m.y*vec_n.y+vec_m.z*vec_n.z);

      T.ax=0.5*CosPhi*(2.0*dmn.ax/dot_mn-dme.ax/m2-dne.ax/n2);
      T.ay=0.5*CosPhi*(2.0*dmn.ay/dot_mn-dme.ay/m2-dne.ay/n2);
      T.az=0.5*CosPhi*(2.0*dmn.az/dot_mn-dme.az/m2-dne.az/n2);

      T.bx=0.5*CosPhi*(2.0*dmn.bx/dot_mn-dme.bx/m2-dne.bx/n2);
      T.by=0.5*CosPhi*(2.0*dmn.by/dot_mn-dme.by/m2-dne.by/n2);
      T.bz=0.5*CosPhi*(2.0*dmn.bz/dot_mn-dme.bz/m2-dne.bz/n2);

      T.cx=0.5*CosPhi*(2.0*dmn.cx/dot_mn-dme.cx/m2-dne.cx/n2);
      T.cy=0.5*CosPhi*(2.0*dmn.cy/dot_mn-dme.cy/m2-dne.cy/n2);
      T.cz=0.5*CosPhi*(2.0*dmn.cz/dot_mn-dme.cz/m2-dne.cz/n2);

      BornTerm[CurrentSystem].xxxx+=
        DDF*T.ax*T.ax+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
        +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
        +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
      BornTerm[CurrentSystem].xxyy+=
        DDF*T.ax*T.by+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.by/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
        +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
        +DF*(1.0/CosPhi)*T.ax*T.by;
      BornTerm[CurrentSystem].xxzz+=
        DDF*T.ax*T.cz+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
        +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
        +DF*(1.0/CosPhi)*T.ax*T.cz;
      BornTerm[CurrentSystem].xxyz+=
        DDF*T.ax*T.bz+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.bz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.bz)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
        +(dme.ax*dme.bz/SQR(m2)+dne.ax*dne.bz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
        +DF*(1.0/CosPhi)*T.ax*T.bz;
      BornTerm[CurrentSystem].xxzx+=
        DDF*T.ax*T.cx+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.cx/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cx)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x)))
        +(dme.ax*dme.cx/SQR(m2)+dne.ax*dne.cx/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.x+vec_u.x*vec_u.x*vec_v.z*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.x+vec_u.x*vec_u.x*vec_w.z*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)))
        +DF*(1.0/CosPhi)*T.ax*T.cx+S.az;
      BornTerm[CurrentSystem].xxxy+=
        DDF*T.ax*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ay)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.ax*dme.ay/SQR(m2)+dne.ax*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.ax*T.ay+S.ay;

      BornTerm[CurrentSystem].yyyy+=
        DDF*T.by*T.by+
        DF*0.5*CosPhi*(
        -4.0*dmn.by*dmn.by/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
        ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
        -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
        +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
        +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
      BornTerm[CurrentSystem].yyzz+=
        DDF*T.by*T.cz+
        DF*0.5*CosPhi*(
        -4.0*dmn.by*dmn.cz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
        ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
        -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
        +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
        +DF*(1.0/CosPhi)*T.by*T.cz;
      BornTerm[CurrentSystem].yyyz+=
        DDF*T.by*T.bz+
        DF*0.5*CosPhi*(
        -4.0*dmn.by*dmn.bz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.bz)+2.0*dot_mn*
        ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
        -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
        +(dme.by*dme.bz/SQR(m2)+dne.by*dne.bz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
        +DF*(1.0/CosPhi)*T.by*T.bz+S.bz;
      BornTerm[CurrentSystem].yyzx+=
        DDF*T.by*T.cx+
        DF*0.5*CosPhi*(
        -4.0*dmn.by*dmn.cx/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cx)+2.0*dot_mn*
        ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
        -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x)))
        +(dme.by*dme.cx/SQR(m2)+dne.by*dne.cx/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.x+vec_u.y*vec_u.y*vec_v.z*vec_v.x)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.x+vec_u.y*vec_u.y*vec_w.z*vec_w.x)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)))
        +DF*(1.0/CosPhi)*T.by*T.cx;
      BornTerm[CurrentSystem].yyxy+=
        DDF*T.by*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.by*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.ay)+2.0*dot_mn*
        ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.by*dme.ay/SQR(m2)+dne.by*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.x*vec_u.y+vec_u.y*vec_u.y*vec_v.x*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.x*vec_u.y+vec_u.y*vec_u.y*vec_w.x*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.by*T.ay+S.bx;

      BornTerm[CurrentSystem].zzzz+=
        DDF*T.cz*T.cz+
        DF*0.5*CosPhi*(
        -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
        ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
        -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
        +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
        +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
      BornTerm[CurrentSystem].zzyz+=
        DDF*T.cz*T.bz+
        DF*0.5*CosPhi*(
        -4.0*dmn.cz*dmn.bz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.bz)+2.0*dot_mn*
        ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
        -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
        +(dme.cz*dme.bz/SQR(m2)+dne.cz*dne.bz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.y*vec_u.z+vec_u.z*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.y*vec_u.z+vec_u.z*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
        +DF*(1.0/CosPhi)*T.cz*T.bz+S.cy;
      BornTerm[CurrentSystem].zzzx+=
        DDF*T.cz*T.cx+
        DF*0.5*CosPhi*(
        -4.0*dmn.cz*dmn.cx/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cx)+2.0*dot_mn*
        ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
        -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x)))
        +(dme.cz*dme.cx/SQR(m2)+dne.cz*dne.cx/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.x+vec_u.z*vec_u.z*vec_v.z*vec_v.x)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.x+vec_u.z*vec_u.z*vec_w.z*vec_w.x)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)))
        +DF*(1.0/CosPhi)*T.cz*T.cx+S.cx;
      BornTerm[CurrentSystem].zzxy+=
        DDF*T.cz*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.cz*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.ay)+2.0*dot_mn*
        ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.cz*dme.ay/SQR(m2)+dne.cz*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.x*vec_u.y+vec_u.z*vec_u.z*vec_v.x*vec_v.y)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.x*vec_u.y+vec_u.z*vec_u.z*vec_w.x*vec_w.y)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.cz*T.ay;

      BornTerm[CurrentSystem].yzyz+=
        DDF*T.bz*T.bz+
        DF*0.5*CosPhi*(
        -4.0*dmn.bz*dmn.bz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.bz)+2.0*dot_mn*
        ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
        -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
        +(dme.bz*dme.bz/SQR(m2)+dne.bz*dne.bz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
        +DF*(1.0/CosPhi)*T.bz*T.bz+0.5*(S.by+S.cz);
      BornTerm[CurrentSystem].yzzx+=
        DDF*T.bz*T.cx+
        DF*0.5*CosPhi*(
        -4.0*dmn.bz*dmn.cx/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.cx)+2.0*dot_mn*
        ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
        -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x)))
        +(dme.bz*dme.cx/SQR(m2)+dne.bz*dne.cx/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.z*vec_u.x+vec_u.y*vec_u.z*vec_v.z*vec_v.x)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.z*vec_u.x+vec_u.y*vec_u.z*vec_w.z*vec_w.x)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)))
        +DF*(1.0/CosPhi)*T.bz*T.cx+0.5*S.bx;
      BornTerm[CurrentSystem].yzxy+=
        DDF*T.bz*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.bz*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.ay)+2.0*dot_mn*
        ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.bz*dme.ay/SQR(m2)+dne.bz*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.x*vec_u.y+vec_u.y*vec_u.z*vec_v.x*vec_v.y)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.x*vec_u.y+vec_u.y*vec_u.z*vec_w.x*vec_w.y)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.bz*T.ay+0.5*S.az;

      BornTerm[CurrentSystem].zxzx+=
        DDF*T.cx*T.cx+
        DF*0.5*CosPhi*(
        -4.0*dmn.cx*dmn.cx/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cx*dmn.cx)+2.0*dot_mn*
        ((vec_v.z*vec_u.x+vec_u.z*vec_v.x)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)+(vec_w.z*vec_u.x+vec_u.z*vec_w.x)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
        -(vec_v.z*vec_w.x+vec_w.z*vec_v.x)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)-(vec_u.z*vec_u.x+vec_u.z*vec_u.x)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x)))
        +(dme.cx*dme.cx/SQR(m2)+dne.cx*dne.cx/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.x*vec_u.z*vec_u.x+vec_u.z*vec_u.x*vec_v.z*vec_v.x)-2.0*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.x*vec_u.z*vec_u.x+vec_u.z*vec_u.x*vec_w.z*vec_w.x)-2.0*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)))
        +DF*(1.0/CosPhi)*T.cx*T.cx+0.5*(S.cz+S.ax);
      BornTerm[CurrentSystem].zxxy+=
        DDF*T.cx*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.cx*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cx*dmn.ay)+2.0*dot_mn*
        ((vec_v.z*vec_u.x+vec_u.z*vec_v.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.z*vec_u.x+vec_u.z*vec_w.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.z*vec_w.x+vec_w.z*vec_v.x)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.z*vec_u.x+vec_u.z*vec_u.x)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.cx*dme.ay/SQR(m2)+dne.cx*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.x*vec_u.x*vec_u.y+vec_u.z*vec_u.x*vec_v.x*vec_v.y)-2.0*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.x*vec_u.x*vec_u.y+vec_u.z*vec_u.x*vec_w.x*vec_w.y)-2.0*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.cx*T.ay+0.5*S.cy;

      BornTerm[CurrentSystem].xyxy+=
        DDF*T.ay*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.ay*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.ay)+2.0*dot_mn*
        ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.ay*dme.ay/SQR(m2)+dne.ay*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.ay*T.ay+0.5*(S.ax+S.by);
    }
  }
}

void CalculateFrameworkImproperTorsionBornTerm(void)
{
  int i,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds,fa,fb,fc,fd;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL *parms;
  REAL DF,DDF;
  VECTOR vec_u,vec_v,vec_w,vec_n,vec_m;
  REAL u,w,v,m2,n2,dot_mn;
  REAL_MATRIX3x3 dme,dne,dmn,S,T;


  UHostImproperTorsion[CurrentSystem]=0.0;
  for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfImproperTorsions[CurrentFramework];i++)
    {
      A=Framework[CurrentSystem].ImproperTorsions[CurrentFramework][i].A;
      B=Framework[CurrentSystem].ImproperTorsions[CurrentFramework][i].B;
      C=Framework[CurrentSystem].ImproperTorsions[CurrentFramework][i].C;
      D=Framework[CurrentSystem].ImproperTorsions[CurrentFramework][i].D;
      parms=Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][i];

      posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
      posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;
      posC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position;
      posD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Position;

      vec_u.x=posB.x-posC.x;
      vec_u.y=posB.y-posC.y;
      vec_u.z=posB.z-posC.z;
      vec_u=ApplyBoundaryCondition(vec_u);
      u=sqrt(SQR(vec_u.x)+SQR(vec_u.y)+SQR(vec_u.z));

      vec_v.x=posD.x-posC.x;
      vec_v.y=posD.y-posC.y;
      vec_v.z=posD.z-posC.z;
      vec_v=ApplyBoundaryCondition(vec_v);
      v=sqrt(SQR(vec_v.x)+SQR(vec_v.y)+SQR(vec_v.z));

      vec_w.x=posA.x-posB.x;
      vec_w.y=posA.y-posB.y;
      vec_w.z=posA.z-posB.z;
      vec_w=ApplyBoundaryCondition(vec_w);
      w=sqrt(SQR(vec_w.x)+SQR(vec_w.y)+SQR(vec_w.z));

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

      switch(Framework[CurrentSystem].ImproperTorsionType[CurrentFramework][i])
      {
        case HARMONIC_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(phi-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // potential defined in terms of 'phi' and therefore contains a singularity
          // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
          // same direction as Rbc, and negative otherwise
          Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
          Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
          Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
          Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
          Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
          Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
          sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
          Phi=SIGN(acos(CosPhi),sign);
          SinPhi=sin(Phi);
          Phi-=parms[1];
          Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
          SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
          U=0.5*parms[0]*SQR(Phi);
          DF=-parms[0]*Phi/SinPhi;
          DDF=-parms[0]*(Phi*CosPhi-SinPhi)/CUBE(SinPhi);
          break;
        case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(cos(phi)-cos(p_1))^2
          // ===============================================
          // p_0/k_B [K]
          // p_1     [degrees]
          U=0.5*parms[0]*SQR(CosPhi-parms[1]);
          DF=parms[0]*(CosPhi-parms[1]);
          DDF=parms[0];
          break;
        case THREE_COSINE_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
          DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[1]-6.0*parms[2]*CosPhi);
          break;
        case MM3_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0     [kcal/mol]
          // p_1     [kcal/mol]
          // p_2     [kcal/mol]
          U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
          DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[1]-6.0*parms[2]*CosPhi);
          break;
        case CFF_IMPROPER_DIHEDRAL:
          // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
          DF=-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2);
          DDF=-4.0*(parms[1]+6.0*parms[2]*CosPhi);
          break;
        case CFF_IMPROPER_DIHEDRAL2:
          // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
          DF=parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi);
          DDF=4.0*(parms[1]+6.0*parms[2]*CosPhi);
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
          DF=-parms[1]+2.0*parms[2]*CosPhi-3.0*parms[3]*CosPhi2+4.0*parms[4]*CosPhi2*CosPhi-5.0*parms[5]*SQR(CosPhi2);
          DDF=2.0*parms[2]-6.0*parms[3]*CosPhi+12.0*parms[4]*CosPhi2-20.0*parms[5]*CosPhi2*CosPhi;
          break;
        case TRAPPE_IMPROPER_DIHEDRAL:
          // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
          // ==========================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
          DF=parms[1]-4.0*parms[2]*CosPhi+3.0*parms[3]*(4.0*CosPhi2-1.0);
          DDF=-4.0*(parms[2]-6.0*parms[3]*CosPhi);
          break;
        case TRAPPE_IMPROPER_DIHEDRAL_EXTENDED:
          // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
          // =============================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          U=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
          DF=parms[1]-3.0*parms[3]+4.0*(parms[2]-4.0*parms[4])*CosPhi+12.0*parms[3]*CosPhi2+32.0*parms[4]*CosPhi2*CosPhi;
          DDF=4.0*parms[2]-16.0*parms[4]+24.0*parms[3]*CosPhi+96.0*parms[4]*CosPhi2;
          break;
        case CVFF_IMPROPER_DIHEDRAL:
          // p_0*(1+cos(p_1*phi-p_2))
          // ========================
          // p_0/k_B [K]
          // p_1     [-]
          // p_2     [degrees]
          // potential defined in terms of 'phi' and therefore contains a singularity
          // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
          // same direction as Dbc, and negative otherwise
          Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
          Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
          Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
          Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
          Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
          Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
          sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
          Phi=SIGN(acos(CosPhi),sign);
          SinPhi=sin(Phi);
          SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
          U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
          DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
          DDF=-(SinPhi*parms[0]*SQR(parms[1])*cos(parms[1]*Phi-parms[2])-parms[0]*parms[1]*sin(parms[1]*Phi-parms[2])*CosPhi)/CUBE(SinPhi);
          break;
        case OPLS_IMPROPER_DIHEDRAL:
          // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
          // =================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
          DF=0.5*parms[1]-2.0*parms[2]*CosPhi+1.5*parms[3]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[2]-6.0*parms[3]*CosPhi);
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
          DF=0.5*(parms[0]-3.0*parms[2]+5.0*parms[4])-2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi+
             6.0*(parms[2]-5.0*parms[4])*CosPhi2-16.0*(parms[3]-6.0*parms[5])*CosPhi2*CosPhi+40.0*parms[4]*SQR(CosPhi2)-96.0*parms[5]*CosPhi2*CUBE(CosPhi);
          DDF=-2.0*parms[1]+8.0*parms[3]-18.0*parms[5]+12.0*(parms[2]-5.0*parms[4])*CosPhi
              -48.0*(parms[3]-6.0*parms[5])*CosPhi2+160.0*parms[4]*CosPhi2*CosPhi-480.0*parms[5]*SQR(CosPhi2);
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
          DF=0.5*parms[0]+parms[2]*(6.0*CosPhi2-1.5)+parms[4]*(2.5-30.0*CosPhi2+40.0*SQR(CosPhi2))+
             CosPhi*(-2.0*parms[1]+parms[3]*(16.0*CosPhi2-8.0)+parms[5]*(18.0-96.0*CosPhi2+96.0*SQR(CosPhi2)));
          DDF=-2.0*parms[1]+12.0*parms[2]*CosPhi+parms[3]*(48.0*CosPhi2-8.0)+parms[4]*CosPhi*(160.0*CosPhi2-60.0)
              +parms[5]*(18.0-288.0*CosPhi2+480.0*SQR(CosPhi2));
          break;
        default:
          fprintf(stderr, "Undefined Improper Torsion potential\n");
          U=DF=DDF=0.0;
          exit(0);
          break;
      }

      // energy
      UHostImproperTorsion[CurrentSystem]+=U;


      // Calculate the first derivative vectors.
      d=dot_ab/rbc;
      e=dot_cd/rbc;

      dtA.x=(ds.x-CosPhi*dr.x)/r;
      dtA.y=(ds.y-CosPhi*dr.y)/r;
      dtA.z=(ds.z-CosPhi*dr.z)/r;

      dtD.x=(dr.x-CosPhi*ds.x)/s;
      dtD.y=(dr.y-CosPhi*ds.y)/s;
      dtD.z=(dr.z-CosPhi*ds.z)/s;

      dtB.x=dtA.x*(d-1.0)+e*dtD.x;
      dtB.y=dtA.y*(d-1.0)+e*dtD.y;
      dtB.z=dtA.z*(d-1.0)+e*dtD.z;

      dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
      dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
      dtC.z=-dtD.z*(e+1.0)-d*dtA.z;


      // forces are oppositely directed to the gradient
      fa.x=-DF*dtA.x;
      fa.y=-DF*dtA.y;
      fa.z=-DF*dtA.z;

      fb.x=-DF*dtB.x;
      fb.y=-DF*dtB.y;
      fb.z=-DF*dtB.z;

      fc.x=-DF*dtC.x;
      fc.y=-DF*dtC.y;
      fc.z=-DF*dtC.z;

      fd.x=-DF*dtD.x;
      fd.y=-DF*dtD.y;
      fd.z=-DF*dtD.z;

      // add contribution to the forces
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.x+=fa.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.y+=fa.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.z+=fa.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.x+=fb.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.y+=fb.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.z+=fb.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.x+=fc.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.y+=fc.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][C].Force.z+=fc.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][D].Force.x+=fd.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][D].Force.y+=fd.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][D].Force.z+=fd.z;

      // add contribution to the stress tensor
      // Note: rbc is here because the vector was normalized before
      S.ax=-(Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x);
      S.bx=-(Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x);
      S.cx=-(Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x);

      S.ay=-(Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y);
      S.by=-(Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y);
      S.cy=-(Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y);

      S.az=-(Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z);
      S.bz=-(Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z);
      S.cz=-(Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z);

      StrainDerivativeTensor[CurrentSystem].ax+=S.ax;
      StrainDerivativeTensor[CurrentSystem].bx+=S.bx;
      StrainDerivativeTensor[CurrentSystem].cx+=S.cx;

      StrainDerivativeTensor[CurrentSystem].ay+=S.bx;
      StrainDerivativeTensor[CurrentSystem].by+=S.by;
      StrainDerivativeTensor[CurrentSystem].cy+=S.cy;

      StrainDerivativeTensor[CurrentSystem].az+=S.az;
      StrainDerivativeTensor[CurrentSystem].bz+=S.bz;
      StrainDerivativeTensor[CurrentSystem].cz+=S.cz;

      dme.ax=2.0*((SQR(u))*vec_v.x*vec_v.x+(SQR(v))*vec_u.x*vec_u.x-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x));
      dme.ay=2.0*((SQR(u))*vec_v.x*vec_v.y+(SQR(v))*vec_u.x*vec_u.y-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y));
      dme.az=2.0*((SQR(u))*vec_v.x*vec_v.z+(SQR(v))*vec_u.x*vec_u.z-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z));
      dme.bx=2.0*((SQR(u))*vec_v.y*vec_v.x+(SQR(v))*vec_u.y*vec_u.x-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.y*vec_u.x+vec_u.y*vec_v.x));
      dme.by=2.0*((SQR(u))*vec_v.y*vec_v.y+(SQR(v))*vec_u.y*vec_u.y-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y));
      dme.bz=2.0*((SQR(u))*vec_v.y*vec_v.z+(SQR(v))*vec_u.y*vec_u.z-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z));
      dme.cx=2.0*((SQR(u))*vec_v.z*vec_v.x+(SQR(v))*vec_u.z*vec_u.x-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x));
      dme.cy=2.0*((SQR(u))*vec_v.z*vec_v.y+(SQR(v))*vec_u.z*vec_u.y-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.z*vec_u.y+vec_u.z*vec_v.y));
      dme.cz=2.0*((SQR(u))*vec_v.z*vec_v.z+(SQR(v))*vec_u.z*vec_u.z-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z));

      dne.ax=2.0*((SQR(u))*vec_w.x*vec_w.x+(SQR(w))*vec_u.x*vec_u.x-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x));
      dne.ay=2.0*((SQR(u))*vec_w.x*vec_w.y+(SQR(w))*vec_u.x*vec_u.y-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y));
      dne.az=2.0*((SQR(u))*vec_w.x*vec_w.z+(SQR(w))*vec_u.x*vec_u.z-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z));
      dne.bx=2.0*((SQR(u))*vec_w.y*vec_w.x+(SQR(w))*vec_u.y*vec_u.x-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.y*vec_u.x+vec_u.y*vec_w.x));
      dne.by=2.0*((SQR(u))*vec_w.y*vec_w.y+(SQR(w))*vec_u.y*vec_u.y-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y));
      dne.bz=2.0*((SQR(u))*vec_w.y*vec_w.z+(SQR(w))*vec_u.y*vec_u.z-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z));
      dne.cx=2.0*((SQR(u))*vec_w.z*vec_w.x+(SQR(w))*vec_u.z*vec_u.x-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x));
      dne.cy=2.0*((SQR(u))*vec_w.z*vec_w.y+(SQR(w))*vec_u.z*vec_u.y-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.z*vec_u.y+vec_u.z*vec_w.y));
      dne.cz=2.0*((SQR(u))*vec_w.z*vec_w.z+(SQR(w))*vec_u.z*vec_u.z-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z));

      dmn.ax=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.x*vec_w.x+vec_w.x*vec_u.x)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x);
      dmn.ay=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.x*vec_w.y+vec_w.x*vec_u.y)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y);
      dmn.az=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.x*vec_w.z+vec_w.x*vec_u.z)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z);

      dmn.bx=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.y*vec_w.x+vec_w.y*vec_u.x)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.y*vec_u.x+vec_u.y*vec_v.x)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.y*vec_u.x+vec_u.y*vec_u.x)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.y*vec_w.x+vec_w.y*vec_v.x);
      dmn.by=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.y*vec_w.y+vec_w.y*vec_u.y)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y);
      dmn.bz=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.y*vec_w.z+vec_w.y*vec_u.z)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z);

      dmn.cx=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.z*vec_w.x+vec_w.z*vec_u.x)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x);
      dmn.cy=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.z*vec_w.y+vec_w.z*vec_u.y)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.z*vec_u.y+vec_u.z*vec_v.y)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.z*vec_u.y+vec_u.z*vec_u.y)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.z*vec_w.y+vec_w.z*vec_v.y);
      dmn.cz=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.z*vec_w.z+vec_w.z*vec_u.z)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z);

      vec_m.x=vec_v.y*vec_u.z-vec_v.z*vec_u.y;
      vec_m.y=vec_v.z*vec_u.x-vec_v.x*vec_u.z;
      vec_m.z=vec_v.x*vec_u.y-vec_v.y*vec_u.x;
      m2=SQR(vec_m.x)+SQR(vec_m.y)+SQR(vec_m.z);

      vec_n.x=vec_u.y*vec_w.z-vec_u.z*vec_w.y;
      vec_n.y=vec_u.z*vec_w.x-vec_u.x*vec_w.z;
      vec_n.z=vec_u.x*vec_w.y-vec_u.y*vec_w.x;
      n2=SQR(vec_n.x)+SQR(vec_n.y)+SQR(vec_n.z);

      dot_mn=(vec_m.x*vec_n.x+vec_m.y*vec_n.y+vec_m.z*vec_n.z);

      T.ax=0.5*CosPhi*(2.0*dmn.ax/dot_mn-dme.ax/m2-dne.ax/n2);
      T.ay=0.5*CosPhi*(2.0*dmn.ay/dot_mn-dme.ay/m2-dne.ay/n2);
      T.az=0.5*CosPhi*(2.0*dmn.az/dot_mn-dme.az/m2-dne.az/n2);

      T.bx=0.5*CosPhi*(2.0*dmn.bx/dot_mn-dme.bx/m2-dne.bx/n2);
      T.by=0.5*CosPhi*(2.0*dmn.by/dot_mn-dme.by/m2-dne.by/n2);
      T.bz=0.5*CosPhi*(2.0*dmn.bz/dot_mn-dme.bz/m2-dne.bz/n2);

      T.cx=0.5*CosPhi*(2.0*dmn.cx/dot_mn-dme.cx/m2-dne.cx/n2);
      T.cy=0.5*CosPhi*(2.0*dmn.cy/dot_mn-dme.cy/m2-dne.cy/n2);
      T.cz=0.5*CosPhi*(2.0*dmn.cz/dot_mn-dme.cz/m2-dne.cz/n2);

      BornTerm[CurrentSystem].xxxx+=
        DDF*T.ax*T.ax+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
        +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
        +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
      BornTerm[CurrentSystem].xxyy+=
        DDF*T.ax*T.by+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.by/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
        +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
        +DF*(1.0/CosPhi)*T.ax*T.by;
      BornTerm[CurrentSystem].xxzz+=
        DDF*T.ax*T.cz+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
        +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
        +DF*(1.0/CosPhi)*T.ax*T.cz;
      BornTerm[CurrentSystem].xxyz+=
        DDF*T.ax*T.bz+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.bz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.bz)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
        +(dme.ax*dme.bz/SQR(m2)+dne.ax*dne.bz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
        +DF*(1.0/CosPhi)*T.ax*T.bz;
      BornTerm[CurrentSystem].xxzx+=
        DDF*T.ax*T.cx+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.cx/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cx)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x)))
        +(dme.ax*dme.cx/SQR(m2)+dne.ax*dne.cx/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.x+vec_u.x*vec_u.x*vec_v.z*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.x+vec_u.x*vec_u.x*vec_w.z*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)))
        +DF*(1.0/CosPhi)*T.ax*T.cx+S.az;
      BornTerm[CurrentSystem].xxxy+=
        DDF*T.ax*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.ax*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ay)+2.0*dot_mn*
        ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.ax*dme.ay/SQR(m2)+dne.ax*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.ax*T.ay+S.ay;

      BornTerm[CurrentSystem].yyyy+=
        DDF*T.by*T.by+
        DF*0.5*CosPhi*(
        -4.0*dmn.by*dmn.by/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
        ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
        -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
        +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
        +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
      BornTerm[CurrentSystem].yyzz+=
        DDF*T.by*T.cz+
        DF*0.5*CosPhi*(
        -4.0*dmn.by*dmn.cz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
        ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
        -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
        +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
        +DF*(1.0/CosPhi)*T.by*T.cz;
      BornTerm[CurrentSystem].yyyz+=
        DDF*T.by*T.bz+
        DF*0.5*CosPhi*(
        -4.0*dmn.by*dmn.bz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.bz)+2.0*dot_mn*
        ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
        -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
        +(dme.by*dme.bz/SQR(m2)+dne.by*dne.bz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
        +DF*(1.0/CosPhi)*T.by*T.bz+S.bz;
      BornTerm[CurrentSystem].yyzx+=
        DDF*T.by*T.cx+
        DF*0.5*CosPhi*(
        -4.0*dmn.by*dmn.cx/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cx)+2.0*dot_mn*
        ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
        -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x)))
        +(dme.by*dme.cx/SQR(m2)+dne.by*dne.cx/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.x+vec_u.y*vec_u.y*vec_v.z*vec_v.x)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.x+vec_u.y*vec_u.y*vec_w.z*vec_w.x)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)))
        +DF*(1.0/CosPhi)*T.by*T.cx;
      BornTerm[CurrentSystem].yyxy+=
        DDF*T.by*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.by*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.ay)+2.0*dot_mn*
        ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.by*dme.ay/SQR(m2)+dne.by*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.x*vec_u.y+vec_u.y*vec_u.y*vec_v.x*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.x*vec_u.y+vec_u.y*vec_u.y*vec_w.x*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.by*T.ay+S.bx;

      BornTerm[CurrentSystem].zzzz+=
        DDF*T.cz*T.cz+
        DF*0.5*CosPhi*(
        -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
        ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
        -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
        +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
        +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
      BornTerm[CurrentSystem].zzyz+=
        DDF*T.cz*T.bz+
        DF*0.5*CosPhi*(
        -4.0*dmn.cz*dmn.bz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.bz)+2.0*dot_mn*
        ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
        -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
        +(dme.cz*dme.bz/SQR(m2)+dne.cz*dne.bz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.y*vec_u.z+vec_u.z*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.y*vec_u.z+vec_u.z*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
        +DF*(1.0/CosPhi)*T.cz*T.bz+S.cy;
      BornTerm[CurrentSystem].zzzx+=
        DDF*T.cz*T.cx+
        DF*0.5*CosPhi*(
        -4.0*dmn.cz*dmn.cx/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cx)+2.0*dot_mn*
        ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
        -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x)))
        +(dme.cz*dme.cx/SQR(m2)+dne.cz*dne.cx/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.x+vec_u.z*vec_u.z*vec_v.z*vec_v.x)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.x+vec_u.z*vec_u.z*vec_w.z*vec_w.x)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)))
        +DF*(1.0/CosPhi)*T.cz*T.cx+S.cx;
      BornTerm[CurrentSystem].zzxy+=
        DDF*T.cz*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.cz*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.ay)+2.0*dot_mn*
        ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.cz*dme.ay/SQR(m2)+dne.cz*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.x*vec_u.y+vec_u.z*vec_u.z*vec_v.x*vec_v.y)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.x*vec_u.y+vec_u.z*vec_u.z*vec_w.x*vec_w.y)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.cz*T.ay;

      BornTerm[CurrentSystem].yzyz+=
        DDF*T.bz*T.bz+
        DF*0.5*CosPhi*(
        -4.0*dmn.bz*dmn.bz/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.bz)+2.0*dot_mn*
        ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
        -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
        +(dme.bz*dme.bz/SQR(m2)+dne.bz*dne.bz/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
        +DF*(1.0/CosPhi)*T.bz*T.bz+0.5*(S.by+S.cz);
      BornTerm[CurrentSystem].yzzx+=
        DDF*T.bz*T.cx+
        DF*0.5*CosPhi*(
        -4.0*dmn.bz*dmn.cx/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.cx)+2.0*dot_mn*
        ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
        -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x)))
        +(dme.bz*dme.cx/SQR(m2)+dne.bz*dne.cx/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.z*vec_u.x+vec_u.y*vec_u.z*vec_v.z*vec_v.x)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.z*vec_u.x+vec_u.y*vec_u.z*vec_w.z*vec_w.x)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)))
        +DF*(1.0/CosPhi)*T.bz*T.cx+0.5*S.bx;
      BornTerm[CurrentSystem].yzxy+=
        DDF*T.bz*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.bz*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.ay)+2.0*dot_mn*
        ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.bz*dme.ay/SQR(m2)+dne.bz*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.x*vec_u.y+vec_u.y*vec_u.z*vec_v.x*vec_v.y)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.x*vec_u.y+vec_u.y*vec_u.z*vec_w.x*vec_w.y)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.bz*T.ay+0.5*S.az;

      BornTerm[CurrentSystem].zxzx+=
        DDF*T.cx*T.cx+
        DF*0.5*CosPhi*(
        -4.0*dmn.cx*dmn.cx/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cx*dmn.cx)+2.0*dot_mn*
        ((vec_v.z*vec_u.x+vec_u.z*vec_v.x)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)+(vec_w.z*vec_u.x+vec_u.z*vec_w.x)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
        -(vec_v.z*vec_w.x+vec_w.z*vec_v.x)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)-(vec_u.z*vec_u.x+vec_u.z*vec_u.x)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x)))
        +(dme.cx*dme.cx/SQR(m2)+dne.cx*dne.cx/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.x*vec_u.z*vec_u.x+vec_u.z*vec_u.x*vec_v.z*vec_v.x)-2.0*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.x*vec_u.z*vec_u.x+vec_u.z*vec_u.x*vec_w.z*vec_w.x)-2.0*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)))
        +DF*(1.0/CosPhi)*T.cx*T.cx+0.5*(S.cz+S.ax);
      BornTerm[CurrentSystem].zxxy+=
        DDF*T.cx*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.cx*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.cx*dmn.ay)+2.0*dot_mn*
        ((vec_v.z*vec_u.x+vec_u.z*vec_v.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.z*vec_u.x+vec_u.z*vec_w.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.z*vec_w.x+vec_w.z*vec_v.x)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.z*vec_u.x+vec_u.z*vec_u.x)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.cx*dme.ay/SQR(m2)+dne.cx*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.z*vec_v.x*vec_u.x*vec_u.y+vec_u.z*vec_u.x*vec_v.x*vec_v.y)-2.0*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.z*vec_w.x*vec_u.x*vec_u.y+vec_u.z*vec_u.x*vec_w.x*vec_w.y)-2.0*(vec_w.z*vec_u.x+vec_u.z*vec_w.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.cx*T.ay+0.5*S.cy;

      BornTerm[CurrentSystem].xyxy+=
        DDF*T.ay*T.ay+
        DF*0.5*CosPhi*(
        -4.0*dmn.ay*dmn.ay/SQR(dot_mn)
        +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.ay)+2.0*dot_mn*
        ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
        -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
        +(dme.ay*dme.ay/SQR(m2)+dne.ay*dne.ay/SQR(n2))
        -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
        -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
        +DF*(1.0/CosPhi)*T.ay*T.ay+0.5*(S.ax+S.by);
    }
  }
}


int CalculateFrameworkIntraVDWBornTerm(void)
{
  int i,j,typeA,typeB,start;
  REAL chargeA,chargeB,energy,DF,DDF;
  REAL rr;
  VECTOR posA,posB,dr,f;
  int f1,f2;

  // Framework-Framework energy
  UHostHostVDW[CurrentSystem]=0.0;

  if(!InternalFrameworkLennardJonesInteractions) return 0;

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
          if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],0))
          {
            typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
            posB=Framework[CurrentSystem].Atoms[f2][j].Position;
            chargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,1.0);

              UHostHostVDW[CurrentSystem]+=energy;


              f.x=DF*dr.x;
              f.y=DF*dr.y;
              f.z=DF*dr.z;

              Framework[CurrentSystem].Atoms[f1][i].Force.x-=f.x;
              Framework[CurrentSystem].Atoms[f1][i].Force.y-=f.y;
              Framework[CurrentSystem].Atoms[f1][i].Force.z-=f.z;

              Framework[CurrentSystem].Atoms[f2][j].Force.x+=f.x;
              Framework[CurrentSystem].Atoms[f2][j].Force.y+=f.y;
              Framework[CurrentSystem].Atoms[f2][j].Force.z+=f.z;

              // add contribution to the strain derivative tensor
              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

              // add contribution to the born term
              AddContributionToBornTerm(DDF,DF,dr);
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraChargeChargeBornTerm(void)
{
  int i,j,typeA,typeB,start;
  REAL chargeA,chargeB,DF,DDF;
  REAL rr,energy;
  VECTOR posA,posB,dr,f;
  int f1,f2;

  UHostHostChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

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
            chargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              PotentialSecondDerivativeCoulombic(chargeA,chargeB,rr,&energy,&DF,&DDF);

              UHostHostChargeChargeReal[CurrentSystem]+=energy;

              // add contribution to the energy
              f.x=DF*dr.x;
              f.y=DF*dr.y;
              f.z=DF*dr.z;

              // forces
              Framework[CurrentSystem].Atoms[f1][i].Force.x-=f.x;
              Framework[CurrentSystem].Atoms[f1][i].Force.y-=f.y;
              Framework[CurrentSystem].Atoms[f1][i].Force.z-=f.z;

              Framework[CurrentSystem].Atoms[f2][j].Force.x+=f.x;
              Framework[CurrentSystem].Atoms[f2][j].Force.y+=f.y;
              Framework[CurrentSystem].Atoms[f2][j].Force.z+=f.z;

              // add contribution to the strain derivative tensor
              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

              // add contribution to the born term
              AddContributionToBornTerm(DDF,DF,dr);
            }
          }
        }
      }
    }
  }
  return 0;
}

void FrameworkAdsorbateVDWBornTerm(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  REAL rr;
  REAL energy,DF,DDF;
  VECTOR posA,posB,dr,f;
  REAL scalingA;

  UHostAdsorbateVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
      scalingA=Adsorbates[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      // energy
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);

          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,scalingA);

            // energy
            UHostAdsorbateVDW[CurrentSystem]+=energy;

            f.x=-DF*dr.x;
            f.y=-DF*dr.y;
            f.z=-DF*dr.z;

            Adsorbates[CurrentSystem][i].Atoms[j].Force.x+=f.x;
            Adsorbates[CurrentSystem][i].Atoms[j].Force.y+=f.y;
            Adsorbates[CurrentSystem][i].Atoms[j].Force.z+=f.z;

            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              Framework[CurrentSystem].Atoms[f1][k].Force.x-=f.x;
              Framework[CurrentSystem].Atoms[f1][k].Force.y-=f.y;
              Framework[CurrentSystem].Atoms[f1][k].Force.z-=f.z;
            }


            // add contribution to the strain derivative tensor
            StrainDerivativeTensor[CurrentSystem].ax-=f.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx-=f.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx-=f.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay-=f.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by-=f.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy-=f.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az-=f.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz-=f.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz-=f.z*dr.z;

            // add contribution to the born term
            AddContributionToBornTerm(DDF,DF,dr);
          }
        }
      }
    }
  }
}

void FrameworkCationVDWBornTerm(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  REAL rr;
  REAL energy,DF,DDF;
  VECTOR posA,posB,dr,f;
  REAL scalingA;

  UHostCationVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].Position;
      scalingA=Cations[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      // energy
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);

          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,scalingA);

            // energy
            UHostCationVDW[CurrentSystem]+=energy;

            f.x=-DF*dr.x;
            f.y=-DF*dr.y;
            f.z=-DF*dr.z;

            Cations[CurrentSystem][i].Atoms[j].Force.x+=f.x;
            Cations[CurrentSystem][i].Atoms[j].Force.y+=f.y;
            Cations[CurrentSystem][i].Atoms[j].Force.z+=f.z;

            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              Framework[CurrentSystem].Atoms[f1][k].Force.x-=f.x;
              Framework[CurrentSystem].Atoms[f1][k].Force.y-=f.y;
              Framework[CurrentSystem].Atoms[f1][k].Force.z-=f.z;
            }


            // add contribution to the strain derivative tensor
            StrainDerivativeTensor[CurrentSystem].ax-=f.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx-=f.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx-=f.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay-=f.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by-=f.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy-=f.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az-=f.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz-=f.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz-=f.z*dr.z;

            // add contribution to the born term
            AddContributionToBornTerm(DDF,DF,dr);
          }
        }
      }
    }
  }
}

void FrameworkAdsorbateChargeChargeBornTerm(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  REAL rr;
  REAL chargeA,chargeB;
  REAL energy,DF,DDF;
  VECTOR posA,posB,dr,f;

  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;

      // energy
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);

          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

            PotentialSecondDerivativeCoulombic(chargeA,chargeB,rr,&energy,&DF,&DDF);

            UHostAdsorbateChargeChargeReal[CurrentSystem]+=energy;

            f.x=DF*dr.x;
            f.y=DF*dr.y;
            f.z=DF*dr.z;

            Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=f.x;
            Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=f.y;
            Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=f.z;

            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              Framework[CurrentSystem].Atoms[f1][k].Force.x+=f.x;
              Framework[CurrentSystem].Atoms[f1][k].Force.y+=f.y;
              Framework[CurrentSystem].Atoms[f1][k].Force.z+=f.z;
            }

            // add contribution to the strain derivative tensor
            StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

            // add contribution to the born term
            AddContributionToBornTerm(DDF,DF,dr);
          }
        }
      }
    }
  }
}

void FrameworkCationChargeChargeBornTerm(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  REAL rr;
  REAL chargeA,chargeB;
  REAL energy,DF,DDF;
  VECTOR posA,posB,dr,f;

  UHostCationChargeChargeReal[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;
      posA=Cations[CurrentSystem][i].Atoms[j].Position;

      // energy
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);

          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

            PotentialSecondDerivativeCoulombic(chargeA,chargeB,rr,&energy,&DF,&DDF);

            UHostCationChargeChargeReal[CurrentSystem]+=energy;

            f.x=DF*dr.x;
            f.y=DF*dr.y;
            f.z=DF*dr.z;

            Cations[CurrentSystem][i].Atoms[j].Force.x-=f.x;
            Cations[CurrentSystem][i].Atoms[j].Force.y-=f.y;
            Cations[CurrentSystem][i].Atoms[j].Force.z-=f.z;

            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              Framework[CurrentSystem].Atoms[f1][k].Force.x+=f.x;
              Framework[CurrentSystem].Atoms[f1][k].Force.y+=f.y;
              Framework[CurrentSystem].Atoms[f1][k].Force.z+=f.z;
            }

            // add contribution to the strain derivative tensor
            StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

            // add contribution to the born term
            AddContributionToBornTerm(DDF,DF,dr);
          }
        }
      }
    }
  }
}

int CalculateFrameworkIntraReplicaVDWBornTerm(void)
{
  int i,j,typeA,typeB,start;
  REAL chargeA,chargeB,energy,DF,DDF;
  REAL rr;
  VECTOR posA,posB,dr,f;
  int ncell,k1,k2,k3,index_j;
  int f1,f2;

  // Framework-Framework energy
  UHostHostVDW[CurrentSystem]=0.0;

  if(!InternalFrameworkLennardJonesInteractions) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        chargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;

        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if((f1==f2)&&(ncell==0)) start=i+1;
              else start=0;
              for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
              {
                index_j=j+ncell*Framework[CurrentSystem].NumberOfAtoms[f2];
                if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][index_j],0))
                {
                  typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
                  posB=Framework[CurrentSystem].Atoms[f2][j].Position;
                  chargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffVDWSquared)
                  {
                    PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,1.0);

                    if(ncell==0)
                      UHostHostVDW[CurrentSystem]+=energy;
                    else
                      UHostHostVDW[CurrentSystem]+=0.5*energy;

                    f.x=DF*dr.x;
                    f.y=DF*dr.y;
                    f.z=DF*dr.z;

                    Framework[CurrentSystem].Atoms[f1][i].Force.x-=f.x;
                    Framework[CurrentSystem].Atoms[f1][i].Force.y-=f.y;
                    Framework[CurrentSystem].Atoms[f1][i].Force.z-=f.z;

                    if(ncell==0)
                    {
                      Framework[CurrentSystem].Atoms[f2][j].Force.x+=f.x;
                      Framework[CurrentSystem].Atoms[f2][j].Force.y+=f.y;
                      Framework[CurrentSystem].Atoms[f2][j].Force.z+=f.z;

                      // add contribution to the strain derivative tensor
                      StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                      // add contribution to the born term
                      AddContributionToBornTerm(DDF,DF,dr);
                    }
                    else
                    {
                      // add contribution to the strain derivative tensor
                      StrainDerivativeTensor[CurrentSystem].ax+=0.5*f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx+=0.5*f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx+=0.5*f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay+=0.5*f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by+=0.5*f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy+=0.5*f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az+=0.5*f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz+=0.5*f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz+=0.5*f.z*dr.z;

                      // add contribution to the born term
                      AddContributionToBornTerm(0.5*DDF,0.5*DF,dr);
                    }
                  }
                }
              }
              ncell++;
            }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraReplicaChargeChargeBornTerm(void)
{
  int i,j,typeA,typeB,start;
  REAL chargeA,chargeB,energy,DF,DDF;
  REAL rr;
  VECTOR posA,posB,dr,f;
  int f1,f2,ncell,k1,k2,k3,index_j;

  // Framework-Adsorbate energy
  UHostHostChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        chargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;

        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if((f1==f2)&&(ncell==0)) start=i+1;
              else start=0;
              for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
              {
                index_j=j+ncell*Framework[CurrentSystem].NumberOfAtoms[f2];
                if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][index_j],1))
                {
                  typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
                  posB=Framework[CurrentSystem].Atoms[f2][j].Position;
                  chargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  {
                    PotentialSecondDerivativeCoulombic(chargeA,chargeB,rr,&energy,&DF,&DDF);

                    if(ncell==0)
                      UHostHostChargeChargeReal[CurrentSystem]+=energy;
                    else
                      UHostHostChargeChargeReal[CurrentSystem]+=0.5*energy;

                    // add contribution to the energy
                    f.x=DF*dr.x;
                    f.y=DF*dr.y;
                    f.z=DF*dr.z;

                    // forces
                    Framework[CurrentSystem].Atoms[f1][i].Force.x-=f.x;
                    Framework[CurrentSystem].Atoms[f1][i].Force.y-=f.y;
                    Framework[CurrentSystem].Atoms[f1][i].Force.z-=f.z;

                    if(ncell==0)
                    {
                      Framework[CurrentSystem].Atoms[f2][j].Force.x+=f.x;
                      Framework[CurrentSystem].Atoms[f2][j].Force.y+=f.y;
                      Framework[CurrentSystem].Atoms[f2][j].Force.z+=f.z;

                      // add contribution to the strain derivative tensor
                      StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                      // add contribution to the born term
                      AddContributionToBornTerm(DDF,DF,dr);
                    }
                    else
                    {
                      // add contribution to the strain derivative tensor
                      StrainDerivativeTensor[CurrentSystem].ax+=0.5*f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx+=0.5*f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx+=0.5*f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay+=0.5*f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by+=0.5*f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy+=0.5*f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az+=0.5*f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz+=0.5*f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz+=0.5*f.z*dr.z;

                      // add contribution to the born term
                      AddContributionToBornTerm(0.5*DDF,0.5*DF,dr);
                    }
                  }
                }
              }
              ncell++;
            }
      }
    }
  }
  return 0;
}


void FrameworkAdsorbateReplicaVDWBornTerm(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  int ncell,k1,k2,k3;
  REAL rr;
  REAL energy,DF,DDF;
  VECTOR posA,posB,dr,f;
  REAL scalingA;

  UHostAdsorbateVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
      scalingA=Adsorbates[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      // energy
      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
            {
              for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
              {
                posB=Framework[CurrentSystem].Atoms[f1][k].Position;
                posB.x+=ReplicaShift[ncell].x;
                posB.y+=ReplicaShift[ncell].y;
                posB.z+=ReplicaShift[ncell].z;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyReplicaBoundaryCondition(dr);

                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                {
                  typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
                  PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,scalingA);

                  // energy
                  UHostAdsorbateVDW[CurrentSystem]+=energy;

                  f.x=DF*dr.x;
                  f.y=DF*dr.y;
                  f.z=DF*dr.z;

                  Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=f.x;
                  Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=f.y;
                  Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=f.z;

                  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                  {
                    Framework[CurrentSystem].Atoms[f1][k].Force.x+=f.x;
                    Framework[CurrentSystem].Atoms[f1][k].Force.y+=f.y;
                    Framework[CurrentSystem].Atoms[f1][k].Force.z+=f.z;
                  }


                  // add contribution to the strain derivative tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                  // add contribution to the born term
                  AddContributionToBornTerm(DDF,DF,dr);
                }
              }
            }
            ncell++;
         }
    }
  }
}

void FrameworkCationReplicaVDWBornTerm(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  int ncell,k1,k2,k3;
  REAL rr;
  REAL energy,DF,DDF;
  VECTOR posA,posB,dr,f;
  REAL scalingA;

  UHostCationVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].Position;
      scalingA=Cations[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      // energy
      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
            {
              for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
              {
                posB=Framework[CurrentSystem].Atoms[f1][k].Position;
                posB.x+=ReplicaShift[ncell].x;
                posB.y+=ReplicaShift[ncell].y;
                posB.z+=ReplicaShift[ncell].z;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyReplicaBoundaryCondition(dr);

                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                {
                  typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
                  PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,scalingA);

                  // energy
                  UHostCationVDW[CurrentSystem]+=energy;

                  f.x=DF*dr.x;
                  f.y=DF*dr.y;
                  f.z=DF*dr.z;

                  Cations[CurrentSystem][i].Atoms[j].Force.x-=f.x;
                  Cations[CurrentSystem][i].Atoms[j].Force.y-=f.y;
                  Cations[CurrentSystem][i].Atoms[j].Force.z-=f.z;

                  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                  {
                    Framework[CurrentSystem].Atoms[f1][k].Force.x+=f.x;
                    Framework[CurrentSystem].Atoms[f1][k].Force.y+=f.y;
                    Framework[CurrentSystem].Atoms[f1][k].Force.z+=f.z;
                  }


                  // add contribution to the strain derivative tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                  // add contribution to the born term
                  AddContributionToBornTerm(DDF,DF,dr);
                }
              }
            }
            ncell++;
         }
    }
  }
}

void FrameworkAdsorbateReplicaChargeChargeBornTerm(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  int ncell,k1,k2,k3;
  REAL r,rr;
  REAL chargeA,chargeB;
  REAL energy,DF,DDF;
  VECTOR posA,posB,dr,f;

  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;

      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
            {
              for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
              {
                posB=Framework[CurrentSystem].Atoms[f1][k].Position;

                posB.x+=ReplicaShift[ncell].x;
                posB.y+=ReplicaShift[ncell].y;
                posB.z+=ReplicaShift[ncell].z;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyReplicaBoundaryCondition(dr);

                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffChargeChargeSquared[CurrentSystem])
                {
                  r=sqrt(rr);
                  typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
                  chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

                  PotentialSecondDerivativeCoulombic(chargeA,chargeB,rr,&energy,&DF,&DDF);

                  UHostAdsorbateChargeChargeReal[CurrentSystem]+=energy;

                  f.x=DF*dr.x;
                  f.y=DF*dr.y;
                  f.z=DF*dr.z;

                  Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=f.x;
                  Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=f.y;
                  Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=f.z;

                  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                  {
                    Framework[CurrentSystem].Atoms[f1][k].Force.x+=f.x;
                    Framework[CurrentSystem].Atoms[f1][k].Force.y+=f.y;
                    Framework[CurrentSystem].Atoms[f1][k].Force.z+=f.z;
                  }

                  // add contribution to the strain derivative tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                  // add contribution to the born term
                  AddContributionToBornTerm(DDF,DF,dr);
                }
              }
            }
            ncell++;
          }
    }
  }
}

void FrameworkCationReplicaChargeChargeBornTerm(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  int ncell,k1,k2,k3;
  REAL rr;
  REAL chargeA,chargeB;
  REAL energy,DF,DDF;
  VECTOR posA,posB,dr,f;

  UHostCationChargeChargeReal[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;
      posA=Cations[CurrentSystem][i].Atoms[j].Position;

      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
            {
              for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
              {
                posB=Framework[CurrentSystem].Atoms[f1][k].Position;

                posB.x+=ReplicaShift[ncell].x;
                posB.y+=ReplicaShift[ncell].y;
                posB.z+=ReplicaShift[ncell].z;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyReplicaBoundaryCondition(dr);

                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffChargeChargeSquared[CurrentSystem])
                {
                  typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
                  chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

                  PotentialSecondDerivativeCoulombic(chargeA,chargeB,rr,&energy,&DF,&DDF);

                  UHostCationChargeChargeReal[CurrentSystem]+=energy;

                  f.x=DF*dr.x;
                  f.y=DF*dr.y;
                  f.z=DF*dr.z;

                  Cations[CurrentSystem][i].Atoms[j].Force.x-=f.x;
                  Cations[CurrentSystem][i].Atoms[j].Force.y-=f.y;
                  Cations[CurrentSystem][i].Atoms[j].Force.z-=f.z;

                  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                  {
                    Framework[CurrentSystem].Atoms[f1][k].Force.x+=f.x;
                    Framework[CurrentSystem].Atoms[f1][k].Force.y+=f.y;
                    Framework[CurrentSystem].Atoms[f1][k].Force.z+=f.z;
                  }

                  // add contribution to the strain derivative tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                  // add contribution to the born term
                  AddContributionToBornTerm(DDF,DF,dr);
                }
              }
            }
            ncell++;
          }
    }
  }
}
