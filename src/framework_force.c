/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2015 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_force.c' is part of RASPA.

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
 *****************************************************************************************************/

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
#include "molecule.h"

// not used
void CalculateFrameworkFullForce(int m)
{
  POINT pos1,pos2;
  REAL energy,force_factor,rr;
  int i,j,f1,typeA,typeB;
  VECTOR dr,f;

  UAdsorbateAdsorbateVDW[CurrentSystem]=0.0;
  for(i=0;i<Adsorbates[CurrentSystem][m].NumberOfAtoms;i++)
  {
    pos1=Adsorbates[CurrentSystem][m].Atoms[i].Position;
    typeA=Adsorbates[CurrentSystem][m].Atoms[i].Type;

    Adsorbates[CurrentSystem][m].Atoms[i].Force.x=0.0;
    Adsorbates[CurrentSystem][m].Atoms[i].Force.y=0.0;
    Adsorbates[CurrentSystem][m].Atoms[i].Force.z=0.0;

    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
      {
        pos2=Framework[CurrentSystem].Atoms[f1][j].Position;
        dr.x=pos1.x-pos2.x;
        dr.y=pos1.y-pos2.y;
        dr.z=pos1.z-pos2.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffVDWSquared)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
          PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);

          // energy
          UAdsorbateAdsorbateVDW[CurrentSystem]+=energy;

          // forces
          f.x=force_factor*dr.x;
          f.y=force_factor*dr.y;
          f.z=force_factor*dr.z;

          Adsorbates[CurrentSystem][m].Atoms[i].Force.x-=f.x;
          Adsorbates[CurrentSystem][m].Atoms[i].Force.y-=f.y;
          Adsorbates[CurrentSystem][m].Atoms[i].Force.z-=f.z;
        }
      }
    }
  }
}

// Calculate the Energy of a bead of type "typeA" at a position "posA".
// posA is defined in xyz (also for MONOCLINIC).
void CalculateFrameworkForceAtPosition(POINT pos,int typeA,REAL *UVDW,VECTOR *Force,REAL *UCoulomb,VECTOR *CoulombForce)
{
  int i,typeB,f1;
  REAL rr,chargeA,chargeB;
  REAL energy,force_factor,temp;
  VECTOR dr;

  *UVDW=0.0;
  *UCoulomb=0.0;
  Force->x=0.0;
  Force->y=0.0;
  Force->z=0.0;
  CoulombForce->x=0.0;
  CoulombForce->y=0.0;
  CoulombForce->z=0.0;

  if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
  {
    *UVDW=InterpolateVDWForceGrid(typeA,pos,Force);
    if((ChargeMethod!=NONE)&&CoulombGrid)
      *UCoulomb=InterpolateCoulombForceGrid(typeA,pos,CoulombForce);
  }
  else
  {
    chargeA=PseudoAtoms[typeA].Charge1;
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeB=Framework[CurrentSystem].Atoms[f1][i].Type;
        dr.x=pos.x-Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition.x;
        dr.y=pos.y-Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition.y;
        dr.z=pos.z-Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<CutOffVDWSquared)
        {
          PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);
          (*UVDW)+=energy;
          Force->x-=force_factor*dr.x;
          Force->y-=force_factor*dr.y;
          Force->z-=force_factor*dr.z;
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            Framework[CurrentSystem].Atoms[f1][i].Force.x+=force_factor*dr.x;
            Framework[CurrentSystem].Atoms[f1][i].Force.y+=force_factor*dr.y;
            Framework[CurrentSystem].Atoms[f1][i].Force.z+=force_factor*dr.z;
          }
        }
        if(rr<CutOffChargeChargeSquared[CurrentSystem])
        {
          chargeB=Framework[CurrentSystem].Atoms[f1][i].Charge;
          PotentialGradientCoulombic(chargeA,chargeB,rr,&energy,&temp);

          (*UCoulomb)+=energy;

          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            Framework[CurrentSystem].Atoms[f1][i].Force.x+=temp*dr.x;
            Framework[CurrentSystem].Atoms[f1][i].Force.y+=temp*dr.y;
            Framework[CurrentSystem].Atoms[f1][i].Force.z+=temp*dr.z;
          }
          CoulombForce->x-=temp*dr.x;
          CoulombForce->y-=temp*dr.y;
          CoulombForce->z-=temp*dr.z;
        }
      }
    }
  }
}

void CalculateFrameworkBondForce(void)
{
  int i,f1;
  REAL r,rr,r1,DF,U;
  REAL temp,temp2,exp_term;
  REAL *parms;
  VECTOR posA,posB,dr,f;
  int A,B;

  UHostBond[CurrentSystem]=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
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

        switch(Framework[CurrentSystem].BondType[f1][i])
        {
          case HARMONIC_BOND:
            // 0.5*p0*SQR(r-p1);
            // ===============================================
            // p_0/k_B [K/A^2]   force constant
            // p_1     [A]       reference bond distance
            U=0.5*parms[0]*SQR(r-parms[1]);
            DF=parms[0]*(r-parms[1])/r;
            break;
          case CORE_SHELL_SPRING:
            U=0.5*parms[0]*SQR(r);
            DF=parms[0];
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
            break;
          case LJ_12_6_BOND:
            // A/r_ij^12-B/r_ij^6
            // ===============================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^6]
            temp=CUBE(1.0/rr);
            U=parms[0]*SQR(temp)-parms[1]*temp;
            DF=6.0*(parms[1]*temp-2.0*parms[0]*SQR(temp))/rr;
            break;
          case LENNARD_JONES_BOND:
            // 4*p_0*((p_1/r)^12-(p_1/r)^6)
            // ===============================================
            // p_0/k_B [K]
            // p_1     [A]
            temp=CUBE(parms[1]/rr);
            U=4.0*parms[0]*(temp*(temp-1.0));
            DF=24.0*parms[0]*(temp*(1.0-2.0*temp))/rr;
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
            DF=-parms[0]*(SIGN(MIN2(fabs(r1),parms[2]),r1))/r;
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
            break;
          case MEASURE_BOND:
            U=DF=0.0;
            break;
          case RIGID_BOND:
            U=DF=0.0;
            break;
          case FIXED_BOND:
            U=DF=0.0;
            break;
          default:
            fprintf(stderr, "Undefined Bond potential in routine 'CalculateFrameworkBondForce' ('framework_force.c')\n");
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
        Framework[CurrentSystem].Atoms[f1][A].Force.x+=f.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=f.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=f.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x-=f.x;
        Framework[CurrentSystem].Atoms[f1][B].Force.y-=f.y;
        Framework[CurrentSystem].Atoms[f1][B].Force.z-=f.z;

        // add contribution to the stress tensor
        StrainDerivativeTensor[CurrentSystem].ax-=dr.x*f.x;
        StrainDerivativeTensor[CurrentSystem].bx-=dr.y*f.x;
        StrainDerivativeTensor[CurrentSystem].cx-=dr.z*f.x;

        StrainDerivativeTensor[CurrentSystem].ay-=dr.x*f.y;
        StrainDerivativeTensor[CurrentSystem].by-=dr.y*f.y;
        StrainDerivativeTensor[CurrentSystem].cy-=dr.z*f.y;

        StrainDerivativeTensor[CurrentSystem].az-=dr.x*f.z;
        StrainDerivativeTensor[CurrentSystem].bz-=dr.y*f.z;
        StrainDerivativeTensor[CurrentSystem].cz-=dr.z*f.z;
      }
    }
  }
}

void CalculateFrameworkUreyBradleyForce(void)
{
  int i,f1;
  REAL r,rr,r1,DF,U;
  REAL temp,temp2,exp_term;
  VECTOR dr,f,posA,posC;
  REAL *parms;
  int A,C;

  UHostUreyBradley[CurrentSystem]=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
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

        switch(Framework[CurrentSystem].UreyBradleyType[f1][i])
        {
          case HARMONIC_UREYBRADLEY:
            // 0.5*p0*SQR(r-p1);
            // ===============================================
            // p_0/k_B [K/A^2]   force constant
            // p_1     [A]       reference bond distance
            U=0.5*parms[0]*SQR(r-parms[1]);
            DF=parms[0]*(r-parms[1])/r;
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
            break;
          case LJ_12_6_UREYBRADLEY:
            // A/r_ij^12-B/r_ij^6
            // ===============================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^6]
            temp=CUBE(1.0/rr);
            U=parms[0]*SQR(temp)-parms[1]*temp;
            DF=6.0*(parms[1]*temp-2.0*parms[0]*SQR(temp))/rr;
            break;
          case LENNARD_JONES_UREYBRADLEY:
            // 4*p_0*((p_1/r)^12-(p_1/r)^6)
            // ===============================================
            // p_0/k_B [K]
            // p_1     [A]
            temp=CUBE(parms[1]/rr);
            U=4.0*parms[0]*(temp*(temp-1.0));
            DF=24.0*parms[0]*(temp*(1.0-2.0*temp))/rr;
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
            DF=-parms[0]*(SIGN(MIN2(fabs(r1),parms[2]),r1))/r;
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
            break;
          case RIGID_UREYBRADLEY:
            U=DF=0.0;
            break;
          case FIXED_UREYBRADLEY:
            U=DF=0.0;
            break;
          default:
            fprintf(stderr, "Undefined Urey-Bradley potential in routine 'CalculateFrameworkUreyBradleyForce' ('framework_force.c')\n");
            exit(0);
            break;
        }

        // add contribution to the Adsorbate Urey-Bradley energy
        UHostUreyBradley[CurrentSystem]+=U;

        // forces are oppositely directed to the gradient
        f.x=-DF*dr.x;
        f.y=-DF*dr.y;
        f.z=-DF*dr.z;

        // add contribution to the forces
        Framework[CurrentSystem].Atoms[f1][A].Force.x+=f.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=f.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=f.z;

        Framework[CurrentSystem].Atoms[f1][C].Force.x-=f.x;
        Framework[CurrentSystem].Atoms[f1][C].Force.y-=f.y;
        Framework[CurrentSystem].Atoms[f1][C].Force.z-=f.z;

        // add contribution to the stress tensor
        StrainDerivativeTensor[CurrentSystem].ax-=dr.x*f.x;
        StrainDerivativeTensor[CurrentSystem].bx-=dr.y*f.x;
        StrainDerivativeTensor[CurrentSystem].cx-=dr.z*f.x;

        StrainDerivativeTensor[CurrentSystem].ay-=dr.x*f.y;
        StrainDerivativeTensor[CurrentSystem].by-=dr.y*f.y;
        StrainDerivativeTensor[CurrentSystem].cy-=dr.z*f.y;

        StrainDerivativeTensor[CurrentSystem].az-=dr.x*f.z;
        StrainDerivativeTensor[CurrentSystem].bz-=dr.y*f.z;
        StrainDerivativeTensor[CurrentSystem].cz-=dr.z*f.z;
      }
    }
  }
}

void CalculateFrameworkBendForce(void)
{
  int i,A,B,C,D,f1;
  REAL *parms,DF,U,temp,temp2;
  REAL CosTheta,Theta,SinTheta;
  REAL rab,rbc,rac,DTDX;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac,fa,fb,fc,dtA,dtB,dtC;;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL rt2,rap2,rcp2;
  REAL terma,termc,delta,delta2,rm,ptrt2,term;
  VECTOR dpdia,dpdic,m,dedip,fd;

  Rab.x=Rab.y=Rab.z=0.0;
  Rad.x=Rad.y=Rad.z=0.0;
  Rcd.x=Rcd.y=Rcd.z=0.0;
  Rbd.x=Rbd.y=Rbd.z=0.0;
  Rbc.x=Rbc.y=Rbc.z=0.0;
  ap.x=ap.y=ap.z=0.0;
  cp.x=cp.y=cp.z=0.0;
  t.x=t.y=t.z=0.0;
  delta=rcp2=rap2=rt2=0.0;
  rbc=rab=temp=0.0;

  UHostBend[CurrentSystem]=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBends[f1];i++)
      {
        A=Framework[CurrentSystem].Bends[f1][i].A;
        B=Framework[CurrentSystem].Bends[f1][i].B;
        C=Framework[CurrentSystem].Bends[f1][i].C;
        D=Framework[CurrentSystem].Bends[f1][i].D;
        parms=Framework[CurrentSystem].BendArguments[f1][i];

        posA=Framework[CurrentSystem].Atoms[f1][A].Position;
        posB=Framework[CurrentSystem].Atoms[f1][B].Position;
        posC=Framework[CurrentSystem].Atoms[f1][C].Position;
        posD=Framework[CurrentSystem].Atoms[f1][D].Position;

        switch(Framework[CurrentSystem].BendType[f1][i])
        {
          case MM3_IN_PLANE_BEND:
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
        SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));
        DTDX=-1.0/sqrt(1.0-SQR(CosTheta));

        parms=Framework[CurrentSystem].BendArguments[f1][i];

        switch(Framework[CurrentSystem].BendType[f1][i])
        {
          case HARMONIC_BEND:
            // (1/2)p_0*(theta-p_1)^2
            // ===============================================
            // p_0/k_B [K/rad^2]
            // p_1     [degrees]
            U=0.5*parms[0]*SQR(Theta-parms[1]);
            DF=parms[0]*(Theta-parms[1])*DTDX;
            break;
          case CORE_SHELL_BEND:
            // (1/2)p_0*(theta-p_1)^2
            // ===============================================
            // p_0/k_B [K/rad^2]
            // p_1     [degrees]
            U=0.5*parms[0]*SQR(Theta-parms[1]);
            DF=parms[0]*(Theta-parms[1])*DTDX;
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
            break;
          case HARMONIC_COSINE_BEND:
            // (1/2)*p_0*(cos(theta)-cos(p_1))^2
            // ===============================================
            // p_0/k_B [K]
            // p_1     [degrees]
            U=0.5*parms[0]*SQR(CosTheta-parms[1]);
            DF=parms[0]*(CosTheta-parms[1]);
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
            break;
          case TAFIPOLSKY_BEND:
            // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
            // ===============================================
            // p_0/k_B [K]
            U=0.5*parms[0]*(1.0+cos(Theta))*(1.0+cos(2.0*Theta));
            DF=parms[0]*CosTheta*(2.0+3.0*CosTheta);
            break;
          case MM3_BEND:
          case MM3_IN_PLANE_BEND:
            // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
            // =================================================================================================
            // p_0/k_B [mdyne A/rad^2]
            // p_1     [degrees]
            U=DF=0;
            temp=RAD2DEG*(Theta-parms[1]);
            temp2=SQR(temp);
            U=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
            DF=parms[0]*RAD2DEG*(2.0-(3.0*0.014-(4.0*5.6e-5-(5.0*7.0e-7-6.0*2.2e-8*temp)*temp)*temp)*temp)*temp*DTDX;
            break;
          case FIXED_BEND:
            U=DF=0.0;
            break;
          case MEASURE_BEND:
            U=DF=0.0;
            break;
          default:
            fprintf(stderr, "Undefined Bend potential in routine 'CalculateFrameworkBendForce' ('framework_force.c')\n");
            exit(0);
            break;
        }

        // add contribution to the energy
        UHostBend[CurrentSystem]+=U;

        switch(Framework[CurrentSystem].BendType[f1][i])
        {
          case MM3_IN_PLANE_BEND:
            DF=parms[0]*RAD2DEG*(2.0-(3.0*0.014-(4.0*5.6e-5-(5.0*7.0e-7-6.0*2.2e-8*temp)*temp)*temp)*temp)*temp;
            m.x=cp.y*ap.z-cp.z*ap.y;
            m.y=cp.z*ap.x-cp.x*ap.z;
            m.z=cp.x*ap.y-cp.y*ap.x;
            rm=sqrt(m.x*m.x+m.y*m.y+m.z*m.z);
            rm=MAX2(rm,0.000001);

            terma=-DF/(rap2*rm);
            termc=DF/(rcp2*rm);
            fa.x=terma*(ap.y*m.z-ap.z*m.y);
            fa.y=terma*(ap.z*m.x-ap.x*m.z);
            fa.z=terma*(ap.x*m.y-ap.y*m.x);
            fc.x=termc*(cp.y*m.z-cp.z*m.y);
            fc.y=termc*(cp.z*m.x-cp.x*m.z);
            fc.z=termc*(cp.x*m.y-cp.y*m.x);
            dedip.x=-(fa.x+fc.x);
            dedip.y=-(fa.y+fc.y);
            dedip.z=-(fa.z+fc.z);

            delta2=2.0*delta;
            ptrt2=(dedip.x*t.x+dedip.y*t.y+dedip.z*t.z)/rt2;

            term=(Rcd.z*Rbd.y-Rcd.y*Rbd.z)+delta2*(t.y*Rcd.z-t.z*Rcd.y);
            dpdia.x=delta*(Rcd.y*dedip.z-Rcd.z*dedip.y)+term*ptrt2;

            term=(Rcd.x*Rbd.z-Rcd.z*Rbd.x)+delta2*(t.z*Rcd.x-t.x*Rcd.z);
            dpdia.y=delta*(Rcd.z*dedip.x-Rcd.x*dedip.z)+term*ptrt2;

            term=(Rcd.y*Rbd.x-Rcd.x*Rbd.y)+delta2*(t.x*Rcd.y-t.y*Rcd.x);
            dpdia.z=delta*(Rcd.x*dedip.y-Rcd.y*dedip.x)+term*ptrt2;

            term=(Rad.y*Rbd.z-Rad.z*Rbd.y)+delta2*(t.z*Rad.y-t.y*Rad.z);
            dpdic.x=delta*(Rad.z*dedip.y-Rad.y*dedip.z)+term*ptrt2;

            term=(Rad.z*Rbd.x-Rad.x*Rbd.z)+delta2*(t.x*Rad.z-t.z*Rad.x);
            dpdic.y=delta*(Rad.x*dedip.z-Rad.z*dedip.x)+term*ptrt2;

            term=(Rad.x*Rbd.y-Rad.y*Rbd.x)+delta2*(t.y*Rad.x-t.x*Rad.y);
            dpdic.z=delta*(Rad.y*dedip.x-Rad.x*dedip.y)+term*ptrt2;

            fa.x=fa.x+dpdia.x;
            fa.y=fa.y+dpdia.y;
            fa.z=fa.z+dpdia.z;
            fb.x=dedip.x;
            fb.y=dedip.y;
            fb.z=dedip.z;
            fc.x=fc.x+dpdic.x;
            fc.y=fc.y+dpdic.y;
            fc.z=fc.z+dpdic.z;
            fd.x=-(fa.x+fb.x+fc.x);
            fd.y=-(fa.y+fb.y+fc.y);
            fd.z=-(fa.z+fb.z+fc.z);

            // add contribution to the forces
            Framework[CurrentSystem].Atoms[f1][A].Force.x-=fa.x;
            Framework[CurrentSystem].Atoms[f1][A].Force.y-=fa.y;
            Framework[CurrentSystem].Atoms[f1][A].Force.z-=fa.z;

            Framework[CurrentSystem].Atoms[f1][B].Force.x-=fb.x;
            Framework[CurrentSystem].Atoms[f1][B].Force.y-=fb.y;
            Framework[CurrentSystem].Atoms[f1][B].Force.z-=fb.z;

            Framework[CurrentSystem].Atoms[f1][C].Force.x-=fc.x;
            Framework[CurrentSystem].Atoms[f1][C].Force.y-=fc.y;
            Framework[CurrentSystem].Atoms[f1][C].Force.z-=fc.z;

            Framework[CurrentSystem].Atoms[f1][D].Force.x-=fd.x;
            Framework[CurrentSystem].Atoms[f1][D].Force.y-=fd.y;
            Framework[CurrentSystem].Atoms[f1][D].Force.z-=fd.z;

            // add contribution to the stress tensor
            StrainDerivativeTensor[CurrentSystem].ax+=Rad.x*fa.x+Rbd.x*fb.x+Rcd.x*fc.x;
            StrainDerivativeTensor[CurrentSystem].bx+=Rad.y*fa.x+Rbd.y*fb.x+Rcd.y*fc.x;
            StrainDerivativeTensor[CurrentSystem].cx+=Rad.z*fa.x+Rbd.z*fb.x+Rcd.z*fc.x;

            StrainDerivativeTensor[CurrentSystem].ay+=Rad.x*fa.y+Rbd.x*fb.y+Rcd.x*fc.y;
            StrainDerivativeTensor[CurrentSystem].by+=Rad.y*fa.y+Rbd.y*fb.y+Rcd.y*fc.y;
            StrainDerivativeTensor[CurrentSystem].cy+=Rad.z*fa.y+Rbd.z*fb.y+Rcd.z*fc.y;

            StrainDerivativeTensor[CurrentSystem].az+=Rad.x*fa.z+Rbd.x*fb.z+Rcd.x*fc.z;
            StrainDerivativeTensor[CurrentSystem].bz+=Rad.y*fa.z+Rbd.y*fb.z+Rcd.y*fc.z;
            StrainDerivativeTensor[CurrentSystem].cz+=Rad.z*fa.z+Rbd.z*fb.z+Rcd.z*fc.z;
            break;
          default:
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
            Framework[CurrentSystem].Atoms[f1][A].Force.x+=fa.x;
            Framework[CurrentSystem].Atoms[f1][A].Force.y+=fa.y;
            Framework[CurrentSystem].Atoms[f1][A].Force.z+=fa.z;

            Framework[CurrentSystem].Atoms[f1][B].Force.x-=(fa.x+fc.x);
            Framework[CurrentSystem].Atoms[f1][B].Force.y-=(fa.y+fc.y);
            Framework[CurrentSystem].Atoms[f1][B].Force.z-=(fa.z+fc.z);

            Framework[CurrentSystem].Atoms[f1][C].Force.x+=fc.x;
            Framework[CurrentSystem].Atoms[f1][C].Force.y+=fc.y;
            Framework[CurrentSystem].Atoms[f1][C].Force.z+=fc.z;

            // add contribution to the stress tensor
            // Note: rab and rbc are here because the vectors were normalized before
            StrainDerivativeTensor[CurrentSystem].ax-=rab*Rab.x*fa.x+rbc*Rbc.x*fc.x;
            StrainDerivativeTensor[CurrentSystem].bx-=rab*Rab.y*fa.x+rbc*Rbc.y*fc.x;
            StrainDerivativeTensor[CurrentSystem].cx-=rab*Rab.z*fa.x+rbc*Rbc.z*fc.x;

            StrainDerivativeTensor[CurrentSystem].ay-=rab*Rab.x*fa.y+rbc*Rbc.x*fc.y;
            StrainDerivativeTensor[CurrentSystem].by-=rab*Rab.y*fa.y+rbc*Rbc.y*fc.y;
            StrainDerivativeTensor[CurrentSystem].cy-=rab*Rab.z*fa.y+rbc*Rbc.z*fc.y;

            StrainDerivativeTensor[CurrentSystem].az-=rab*Rab.x*fa.z+rbc*Rbc.x*fc.z;
            StrainDerivativeTensor[CurrentSystem].bz-=rab*Rab.y*fa.z+rbc*Rbc.y*fc.z;
            StrainDerivativeTensor[CurrentSystem].cz-=rab*Rab.z*fa.z+rbc*Rbc.z*fc.z;
            break;
        }
      }
    }
  }
}

void CalculateFrameworkInversionBendForce(void)
{
  int i,A,B,C,D,f1;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2,rrbc,rbc2,rbd2,rad2,rac2,dot;
  REAL CosChi,Chi,energy;
  REAL temp,temp2;
  VECTOR Rab,Rbc,Rbd,Rad,Rac;
  POINT posA,posB,posC,posD;
  VECTOR fa,fb,fc,fd;
  REAL term,dedcos;
  VECTOR dccdia,dccdic,dccdid;
  VECTOR deedia,deedic,deedid;


  UHostInversionBend[CurrentSystem]=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
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
        rbc2=Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z;
        rrbc=sqrt(rbc2);

        Rbd.x=posD.x-posB.x;
        Rbd.y=posD.y-posB.y;
        Rbd.z=posD.z-posB.z;
        Rbd=ApplyBoundaryCondition(Rbd);
        rbd2=Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z;

        Rac.x=posC.x-posA.x;
        Rac.y=posC.y-posA.y;
        Rac.z=posC.z-posA.z;
        Rac=ApplyBoundaryCondition(Rac);
        rac2=Rac.x*Rac.x+Rac.y*Rac.y+Rac.z*Rac.z;

        Rad.x=posD.x-posA.x;
        Rad.y=posD.y-posA.y;
        Rad.z=posD.z-posA.z;
        Rad=ApplyBoundaryCondition(Rad);
        rad2=Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z;

        switch(Framework[CurrentSystem].InversionBendType[f1][i])
        {
          case HARMONIC_INVERSION:
          case HARMONIC_COSINE_INVERSION:
          case PLANAR_INVERSION:
            // w is a vector perpendicular to the B-C-D plane
            // c=w.w=(Rbc x Rbd).(Rbc x Rbd)= r_bc^2 r_bd^2 - (r_cb . r_bd)^2
            dot=Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z;
            c=rbc2*rbd2-SQR(dot);
            //c=(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z)*(Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z)-SQR(Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z);
            break;
          case MM3_INVERSION:
          case HARMONIC_INVERSION2:
          case HARMONIC_COSINE_INVERSION2:
          case PLANAR_INVERSION2:
            // w is a vector perpendicular to the A-C-D plane
            // c=w.w=(Rcd x Rad).(Rcd x Rad)=r_cd^2 r_ad^2 - (r_da . r_cd)^2
            dot=Rac.x*Rad.x+Rac.y*Rad.y+Rac.z*Rad.z;
            c=rac2*rad2-SQR(dot);
            //c=(Rcd.x*Rcd.x+Rcd.y*Rcd.y+Rcd.z*Rcd.z)*(Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z)-SQR(Rad.x*Rcd.x+Rad.y*Rcd.y+Rad.z*Rcd.z);
            break;
          default:
            fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
            exit(0);
            break;
        }

        e=Rab.x*(Rbc.y*Rbd.z-Rbc.z*Rbd.y)+Rab.y*(Rbc.z*Rbd.x-Rbc.x*Rbd.z)+Rab.z*(Rbc.x*Rbd.y-Rbc.y*Rbd.x);
        CosChi=sqrt((Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z)-SQR(e)/c)/rrab;

        // Ensure CosChi is between -1 and 1.
        CosChi=SIGN(MIN2(fabs(CosChi),(REAL)1.0),CosChi);

        parms=Framework[CurrentSystem].InversionBendArguments[f1][i];

        switch(Framework[CurrentSystem].InversionBendType[f1][i])
        {
          case HARMONIC_INVERSION:
          case HARMONIC_INVERSION2:
            // (1/2)*p_0*(chi-p_1)^2
            // ===============================================
            // p_0/k_B [K]
            // p_1     [degrees]
            Chi=acos(CosChi);
            energy=0.5*parms[0]*SQR(Chi-parms[1]);
            dedcos=-SIGN(1.0,e)*(parms[0]*(Chi-parms[1])/sqrt(c*(rab2-e*e/c)));
            break;
          case HARMONIC_COSINE_INVERSION:
          case HARMONIC_COSINE_INVERSION2:
            // (1/2)*p_0*(cos(phi)-cos(p_1))^2
            // ===============================================
            // p_0/k_B [K]
            // p_1     [degrees]
            Chi=acos(CosChi);
            energy=0.5*parms[0]*SQR(CosChi-parms[1]);
            dedcos=SIGN(1.0,e)*parms[0]*(CosChi-parms[1])*sin(Chi)/sqrt(c*(rab2-e*e/c));
            break;
          case PLANAR_INVERSION:
          case PLANAR_INVERSION2:
            // (1/2)*p_0*(1-cos(phi))
            // ===============================================
            // p_0/k_B [K]
            Chi=acos(CosChi);
            energy=parms[0]*(1.0-CosChi);
            dedcos=-SIGN(1.0,e)*parms[0]*sin(Chi)/sqrt(c*(rab2-e*e/c));
            break;
          case MM3_INVERSION:
            // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
            // =================================================================================================
            // p_0/k_B [mdyne A/rad^2]
            // p_1     [degrees]
            Chi=acos(CosChi);
            temp=RAD2DEG*(Chi-parms[1]);
            temp2=SQR(temp);
            energy=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
            dedcos=-SIGN(1.0,e)*parms[0]*temp*RAD2DEG*(2.0-3.0*0.014*temp+4.0*5.6e-5*temp2-5.0*7.0e-7*temp*temp2+6.0*2.2e-8*SQR(temp2))/sqrt(c*(rab2-e*e/c));
            break;
          default:
            fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
            exit(0);
            break;
        }

        // energy
        UHostInversionBend[CurrentSystem]+=energy;

        switch(Framework[CurrentSystem].InversionBendType[f1][i])
        {
          case HARMONIC_COSINE_INVERSION:
          case PLANAR_INVERSION:
          case HARMONIC_INVERSION:
            term=e/c;
            dccdia.x=0.0;
            dccdia.y=0.0;
            dccdia.z=0.0;
            dccdic.x=(Rbc.x*rbd2-Rbd.x*dot)*term;
            dccdic.y=(Rbc.y*rbd2-Rbd.y*dot)*term;
            dccdic.z=(Rbc.z*rbd2-Rbd.z*dot)*term;
            dccdid.x=(Rbd.x*rbc2-Rbc.x*dot)*term;
            dccdid.y=(Rbd.y*rbc2-Rbc.y*dot)*term;
            dccdid.z=(Rbd.z*rbc2-Rbc.z*dot)*term;
            break;
          case HARMONIC_INVERSION2:
          case HARMONIC_COSINE_INVERSION2:
          case PLANAR_INVERSION2:
          case MM3_INVERSION:
            term=e/c;
            dccdic.x=(Rac.x*rad2-Rad.x*dot)*term;
            dccdic.y=(Rac.y*rad2-Rad.y*dot)*term;
            dccdic.z=(Rac.z*rad2-Rad.z*dot)*term;
            dccdid.x=(Rad.x*rac2-Rac.x*dot)*term;
            dccdid.y=(Rad.y*rac2-Rac.y*dot)*term;
            dccdid.z=(Rad.z*rac2-Rac.z*dot)*term;
            dccdia.x=-(dccdic.x+dccdid.x);
            dccdia.y=-(dccdic.y+dccdid.y);
            dccdia.z=-(dccdic.z+dccdid.z);
            break;
          default:
            fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
            exit(0);
            break;
        }


        term=e/rab2;
        deedia.x=Rbd.y*Rbc.z-Rbd.z*Rbc.y+Rab.x*term;
        deedia.y=Rbd.z*Rbc.x-Rbd.x*Rbc.z+Rab.y*term;
        deedia.z=Rbd.x*Rbc.y-Rbd.y*Rbc.x+Rab.z*term;
        deedic.x=Rab.y*Rbd.z-Rab.z*Rbd.y;
        deedic.y=Rab.z*Rbd.x-Rab.x*Rbd.z;
        deedic.z=Rab.x*Rbd.y-Rab.y*Rbd.x;
        deedid.x=Rbc.y*Rab.z-Rbc.z*Rab.y;
        deedid.y=Rbc.z*Rab.x-Rbc.x*Rab.z;
        deedid.z=Rbc.x*Rab.y-Rbc.y*Rab.x;

        fa.x=dedcos*(dccdia.x+deedia.x);
        fa.y=dedcos*(dccdia.y+deedia.y);
        fa.z=dedcos*(dccdia.z+deedia.z);
        fc.x=dedcos*(dccdic.x+deedic.x);
        fc.y=dedcos*(dccdic.y+deedic.y);
        fc.z=dedcos*(dccdic.z+deedic.z);
        fd.x=dedcos*(dccdid.x+deedid.x);
        fd.y=dedcos*(dccdid.y+deedid.y);
        fd.z=dedcos*(dccdid.z+deedid.z);

        fb.x=-(fa.x+fc.x+fd.x);
        fb.y=-(fa.y+fc.y+fd.y);
        fb.z=-(fa.z+fc.z+fd.z);

        // add contribution to the forces
        Framework[CurrentSystem].Atoms[f1][A].Force.x-=fa.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y-=fa.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z-=fa.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x-=fb.x;
        Framework[CurrentSystem].Atoms[f1][B].Force.y-=fb.y;
        Framework[CurrentSystem].Atoms[f1][B].Force.z-=fb.z;

        Framework[CurrentSystem].Atoms[f1][C].Force.x-=fc.x;
        Framework[CurrentSystem].Atoms[f1][C].Force.y-=fc.y;
        Framework[CurrentSystem].Atoms[f1][C].Force.z-=fc.z;

        Framework[CurrentSystem].Atoms[f1][D].Force.x-=fd.x;
        Framework[CurrentSystem].Atoms[f1][D].Force.y-=fd.y;
        Framework[CurrentSystem].Atoms[f1][D].Force.z-=fd.z;

        // add contribution to the stress tensor
        StrainDerivativeTensor[CurrentSystem].ax+=Rab.x*fa.x+Rbc.x*fc.x+Rbd.x*fd.x;
        StrainDerivativeTensor[CurrentSystem].ay+=Rab.x*fa.y+Rbc.x*fc.y+Rbd.x*fd.y;
        StrainDerivativeTensor[CurrentSystem].az+=Rab.x*fa.z+Rbc.x*fc.z+Rbd.x*fd.z;

        StrainDerivativeTensor[CurrentSystem].bx+=Rab.y*fa.x+Rbc.y*fc.x+Rbd.y*fd.x;
        StrainDerivativeTensor[CurrentSystem].by+=Rab.y*fa.y+Rbc.y*fc.y+Rbd.y*fd.y;
        StrainDerivativeTensor[CurrentSystem].bz+=Rab.y*fa.z+Rbc.y*fc.z+Rbd.y*fd.z;

        StrainDerivativeTensor[CurrentSystem].cx+=Rab.z*fa.x+Rbc.z*fc.x+Rbd.z*fd.x;
        StrainDerivativeTensor[CurrentSystem].cy+=Rab.z*fa.y+Rbc.z*fc.y+Rbd.z*fd.y;
        StrainDerivativeTensor[CurrentSystem].cz+=Rab.z*fa.z+Rbc.z*fc.z+Rbd.z*fd.z;
      }
    }
  }
}

void CalculateFrameworkTorsionForce(void)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc,U,DF;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  VECTOR fa,fb,fc,fd;
  REAL *parms;

  UHostTorsion[CurrentSystem]=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[f1];i++)
      {
        A=Framework[CurrentSystem].Torsions[f1][i].A;
        B=Framework[CurrentSystem].Torsions[f1][i].B;
        C=Framework[CurrentSystem].Torsions[f1][i].C;
        D=Framework[CurrentSystem].Torsions[f1][i].D;
        parms=Framework[CurrentSystem].TorsionArguments[f1][i];

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

        switch(Framework[CurrentSystem].TorsionType[f1][i])
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
            break;
          case HARMONIC_COSINE_DIHEDRAL:
            // (1/2)*p_0*(cos(phi)-cos(p_1))^2
            // ===============================================
            // p_0/k_B [K]
            // p_1     [degrees]
            U=0.5*parms[0]*SQR(CosPhi-parms[1]);
            DF=parms[0]*(CosPhi-parms[1]);
            break;
          case THREE_COSINE_DIHEDRAL:
            // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
            // ========================================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
            DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
            break;
          case MM3_DIHEDRAL:
            // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
            // ========================================================================
            // p_0     [kcal/mol]
            // p_1     [kcal/mol]
            // p_2     [kcal/mol]
            U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
            DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
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
            break;
          case CFF_DIHEDRAL:
            // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
            // ======================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
            DF=-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2);
            break;
          case CFF_DIHEDRAL2:
            // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
            // ======================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
            DF=parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi);
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
            break;
          case TRAPPE_DIHEDRAL:
            // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
            // =============================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            // p_3/k_B [K]
            U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
            DF=parms[1]-4.0*parms[2]*CosPhi+3.0*parms[3]*(4.0*CosPhi2-1.0);
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
            SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
            U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
            DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
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
            break;
          default:
            fprintf(stderr, "Undefined Torsion potential in routine 'CalculateFrameworkTorsionForce' ('framework_force.c')\n");
            exit(0);
            break;
        }

        // energy
        UHostTorsion[CurrentSystem]+=U;

        // virial
        // pure torsion has *no* contributions to the virial

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
        Framework[CurrentSystem].Atoms[f1][A].Force.x+=fa.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=fa.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=fa.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x+=fb.x;
        Framework[CurrentSystem].Atoms[f1][B].Force.y+=fb.y;
        Framework[CurrentSystem].Atoms[f1][B].Force.z+=fb.z;

        Framework[CurrentSystem].Atoms[f1][C].Force.x+=fc.x;
        Framework[CurrentSystem].Atoms[f1][C].Force.y+=fc.y;
        Framework[CurrentSystem].Atoms[f1][C].Force.z+=fc.z;

        Framework[CurrentSystem].Atoms[f1][D].Force.x+=fd.x;
        Framework[CurrentSystem].Atoms[f1][D].Force.y+=fd.y;
        Framework[CurrentSystem].Atoms[f1][D].Force.z+=fd.z;

        // add contribution to the stress tensor
        // Note: rbc is here because the vector was normalized before
        StrainDerivativeTensor[CurrentSystem].ax-=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
        StrainDerivativeTensor[CurrentSystem].bx-=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
        StrainDerivativeTensor[CurrentSystem].cx-=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

        StrainDerivativeTensor[CurrentSystem].ay-=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
        StrainDerivativeTensor[CurrentSystem].by-=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
        StrainDerivativeTensor[CurrentSystem].cy-=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

        StrainDerivativeTensor[CurrentSystem].az-=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
        StrainDerivativeTensor[CurrentSystem].bz-=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
        StrainDerivativeTensor[CurrentSystem].cz-=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;
      }
    }
  }
}

void CalculateFrameworkImproperTorsionForce(void)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc,U,DF;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  VECTOR fa,fb,fc,fd;
  REAL *parms;

  UHostImproperTorsion[CurrentSystem]=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfImproperTorsions[f1];i++)
      {
        A=Framework[CurrentSystem].ImproperTorsions[f1][i].A;
        B=Framework[CurrentSystem].ImproperTorsions[f1][i].B;
        C=Framework[CurrentSystem].ImproperTorsions[f1][i].C;
        D=Framework[CurrentSystem].ImproperTorsions[f1][i].D;
        parms=Framework[CurrentSystem].ImproperTorsionArguments[f1][i];

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

        switch(Framework[CurrentSystem].ImproperTorsionType[f1][i])
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
            break;
          case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
            // (1/2)*p_0*(cos(phi)-cos(p_1))^2
            // ===============================================
            // p_0/k_B [K]
            // p_1     [degrees]
            U=0.5*parms[0]*SQR(CosPhi-parms[1]);
            DF=parms[0]*(CosPhi-parms[1]);
            break;
          case THREE_COSINE_IMPROPER_DIHEDRAL:
            // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
            // ========================================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
            DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
            break;
          case MM3_IMPROPER_DIHEDRAL:
            // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
            // ========================================================================
            // p_0     [kcal/mol]
            // p_1     [kcal/mol]
            // p_2     [kcal/mol]
            U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
            DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
            break;
          case CFF_IMPROPER_DIHEDRAL:
            // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
            // ======================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
            DF=-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2);
            break;
          case CFF_IMPROPER_DIHEDRAL2:
            // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
            // ======================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
            DF=parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi);
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
            break;
          case TRAPPE_IMPROPER_DIHEDRAL:
            // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
            // =============================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            // p_3/k_B [K]
            U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
            DF=parms[1]-4.0*parms[2]*CosPhi+3.0*parms[3]*(4.0*CosPhi2-1.0);
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
            break;
          case FIXED_IMPROPER_DIHEDRAL:
            U=DF=0.0;
            break;
          default:
            fprintf(stderr, "Undefined Improper-Torsion potential in routine 'CalculateFrameworkImproperTorsionForce' ('framework_force.c')\n");
            exit(0);
            break;
        }

        // energy
        UHostImproperTorsion[CurrentSystem]+=U;

        // virial
        // pure torsion has *no* contributions to the virial

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
        Framework[CurrentSystem].Atoms[f1][A].Force.x+=fa.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=fa.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=fa.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x+=fb.x;
        Framework[CurrentSystem].Atoms[f1][B].Force.y+=fb.y;
        Framework[CurrentSystem].Atoms[f1][B].Force.z+=fb.z;

        Framework[CurrentSystem].Atoms[f1][C].Force.x+=fc.x;
        Framework[CurrentSystem].Atoms[f1][C].Force.y+=fc.y;
        Framework[CurrentSystem].Atoms[f1][C].Force.z+=fc.z;

        Framework[CurrentSystem].Atoms[f1][D].Force.x+=fd.x;
        Framework[CurrentSystem].Atoms[f1][D].Force.y+=fd.y;
        Framework[CurrentSystem].Atoms[f1][D].Force.z+=fd.z;

        // add contribution to the stress tensor
        // Note: rbc is here because the vector was normalized before
        StrainDerivativeTensor[CurrentSystem].ax-=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
        StrainDerivativeTensor[CurrentSystem].bx-=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
        StrainDerivativeTensor[CurrentSystem].cx-=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

        StrainDerivativeTensor[CurrentSystem].ay-=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
        StrainDerivativeTensor[CurrentSystem].by-=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
        StrainDerivativeTensor[CurrentSystem].cy-=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

        StrainDerivativeTensor[CurrentSystem].az-=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
        StrainDerivativeTensor[CurrentSystem].bz-=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
        StrainDerivativeTensor[CurrentSystem].cz-=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;
      }
    }
  }
}


void CalculateFrameworkBondBondForce(void)
{
  int i,A,B,C,f1;
  REAL *parms,gamma,gammc;
  REAL energy;
  REAL rab,rbc;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc,fa,fc;

  UHostBondBond[CurrentSystem]=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondBonds[f1];i++)
      {
        A=Framework[CurrentSystem].BondBonds[f1][i].A;
        B=Framework[CurrentSystem].BondBonds[f1][i].B;
        C=Framework[CurrentSystem].BondBonds[f1][i].C;
        parms=Framework[CurrentSystem].BondBondArguments[f1][i];

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

        switch(Framework[CurrentSystem].BondBondType[f1][i])
        {
          case CVFF_BOND_BOND_CROSS:
          case CFF_BOND_BOND_CROSS:
            // p_0*(rab-p_1)*(rbc-p_2)
            // =======================
            // p_0/k_B [K/A^2]
            // p_1     [A]
            // p_2     [A]
            energy=parms[0]*(rab-parms[1])*(rbc-parms[2]);
            gamma=-parms[0]*(rbc-parms[2])/rab;
            gammc=-parms[0]*(rab-parms[1])/rbc;
          break;
          default:
            fprintf(stderr, "Undefined Bond-Bond potential in routine 'CalculateFrameworkBondBondForce' ('framework_force.c')\n");
            exit(0);
            break;
        }

        // add contribution to the energy
        UHostBondBond[CurrentSystem]+=energy;

        // forces
        fa.x=gamma*Rab.x;
        fa.y=gamma*Rab.y;
        fa.z=gamma*Rab.z;

        fc.x=gammc*Rbc.x;
        fc.y=gammc*Rbc.y;
        fc.z=gammc*Rbc.z;

        // add contribution to the forces
        Framework[CurrentSystem].Atoms[f1][A].Force.x+=fa.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=fa.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=fa.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x-=(fa.x+fc.x);
        Framework[CurrentSystem].Atoms[f1][B].Force.y-=(fa.y+fc.y);
        Framework[CurrentSystem].Atoms[f1][B].Force.z-=(fa.z+fc.z);

        Framework[CurrentSystem].Atoms[f1][C].Force.x+=fc.x;
        Framework[CurrentSystem].Atoms[f1][C].Force.y+=fc.y;
        Framework[CurrentSystem].Atoms[f1][C].Force.z+=fc.z;

        // add contribution to the stress tensor
        StrainDerivativeTensor[CurrentSystem].ax-=Rab.x*fa.x+Rbc.x*fc.x;
        StrainDerivativeTensor[CurrentSystem].bx-=Rab.y*fa.x+Rbc.y*fc.x;
        StrainDerivativeTensor[CurrentSystem].cx-=Rab.z*fa.x+Rbc.z*fc.x;

        StrainDerivativeTensor[CurrentSystem].ay-=Rab.x*fa.y+Rbc.x*fc.y;
        StrainDerivativeTensor[CurrentSystem].by-=Rab.y*fa.y+Rbc.y*fc.y;
        StrainDerivativeTensor[CurrentSystem].cy-=Rab.z*fa.y+Rbc.z*fc.y;

        StrainDerivativeTensor[CurrentSystem].az-=Rab.x*fa.z+Rbc.x*fc.z;
        StrainDerivativeTensor[CurrentSystem].bz-=Rab.y*fa.z+Rbc.y*fc.z;
        StrainDerivativeTensor[CurrentSystem].cz-=Rab.z*fa.z+Rbc.z*fc.z;
      }
    }
  }
}

void CalculateFrameworkBondBendForce(void)
{
  int i,A,B,C,f1;
  REAL *parms,gamma,gamsa,gamsc,pterm,vterm;
  REAL cost,theta,sint;
  REAL rab,rbc,rac;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc,Rac,fa,fc;

  UHostBondBend[CurrentSystem]=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondBends[f1];i++)
      {
        A=Framework[CurrentSystem].BondBends[f1][i].A;
        B=Framework[CurrentSystem].BondBends[f1][i].B;
        C=Framework[CurrentSystem].BondBends[f1][i].C;
        parms=Framework[CurrentSystem].BondBendArguments[f1][i];

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
        cost=MIN2(1.0,MAX2(-1.0,cost));
        theta=acos(cost);
        sint=sqrt(1.0-SQR(cost));
        //sint=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(cost)));

        switch(Framework[CurrentSystem].BondBendType[f1][i])
        {
          case CVFF_BOND_BEND_CROSS:
          case CFF_BOND_BEND_CROSS:
            // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
            // =========================================
            // p_0     [degrees]
            // p_1/k_B [K/A/rad]
            // p_2     [A]
            // p_3/k_B [K/A/rad]
            // p_4     [A]
            pterm=(theta-parms[0])*(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]));
            vterm=0.0;
            gamma=(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]))/sint;
            gamsa=-parms[1]*(theta-parms[0]);
            gamsc=-parms[3]*(theta-parms[0]);
            break;
          case MM3_BOND_BEND_CROSS:
            // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
            // =====================================
            // p_0     [mdyne/rad]
            // p_1     [A]
            // p_2     [A]
            // p_3     [degrees]
            pterm=parms[0]*((rab-parms[1])+(rbc-parms[2]))*RAD2DEG*(theta-parms[3]);
            vterm=0.0;
            gamma=parms[0]*RAD2DEG*((rab-parms[1])+(rbc-parms[2]))/sint;
            gamsa=-parms[0]*RAD2DEG*(theta-parms[3]);
            gamsc=-parms[0]*RAD2DEG*(theta-parms[3]);
            break;
          case TRUNCATED_HARMONIC:
            // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
            // ================================================================
            // p_0/k_B [K/rad^2]
            // p_1     [degrees]
            // p_2     [A]
            pterm=0.5*parms[0]*SQR(theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8));
            vterm=-8.0*pterm*(pow(rab,8)+pow(rbc,8))/pow(parms[2],8);
            gamma=(parms[0]*(theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8)))/sint;
            gamsa=(8.0*pterm/pow(parms[2],8))*pow(rab,7);
            gamsc=(8.0*pterm/pow(parms[2],8))*pow(rbc,7);
            break;
          case SCREENED_HARMONIC:
            // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
            // ===============================================
            // p_0/k_B [K/rad^2]
            // p_1     [degrees]
            // p_2     [A]
            // p_3     [A]
            pterm=0.5*parms[0]*SQR(theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3]));
            vterm=-pterm*(rab/parms[2]+rbc/parms[3]);
            gamma=(parms[0]*(theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3])))/sint;
            gamsa=(pterm/parms[2]);
            gamsc=(pterm/parms[3]);
            break;
          case SCREENED_VESSAL:
            // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
            // ============================================================================
            // p_0/k_B [K/rad^2]
            // p_1     [degrees]
            // p_2     [A]
            // p_3     [A]
            pterm=(parms[0]/(8.0*SQR(theta-M_PI)))*SQR(SQR(parms[1]-M_PI)-SQR(theta-M_PI))
                  *exp(-(rab/parms[2]+rbc/parms[3]));
            vterm=-pterm*(rab/parms[2]+rbc/parms[3]);
            gamma=(parms[0]/(2.0*SQR(theta-M_PI))*(SQR(parms[1]-M_PI)-SQR(theta-M_PI))*(theta-M_PI)
                  *exp(-(rab/parms[2]+rbc/parms[3])))/sint;
            gamsa=(pterm/parms[2]);
            gamsc=(pterm/parms[3]);
            break;
          case TRUNCATED_VESSAL:
            // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
            //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
            // ============================================================================
            // p_0/k_B [K/rad^(4+p_2)]
            // p_1     [degrees]
            // p_2     [-]
            // p_3     [A]
            pterm=parms[0]*(pow(theta,parms[2])*SQR(theta-parms[1])*SQR(theta+parms[1]-2.0*M_PI)
                  -0.5*parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*SQR(theta-parms[1])*pow(M_PI-parms[1],(REAL)3.0))
                  *exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8));
            vterm=-8.0*pterm*(pow(rab,8)+pow(rbc,8))/pow(parms[3],8);
            gamma=(parms[0]*(pow(theta,(parms[2]-1.0))*(theta-parms[1])*(theta+parms[1]-2.0*M_PI)
                   *((parms[2]+4.0)*SQR(theta)-2.0*M_PI*(parms[2]+2.0)*theta+parms[2]*parms[1]*
                   (2.0*M_PI-parms[1]))-parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*(theta-parms[1])*
                   pow(M_PI-parms[1],3))*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8)))/sint;
            gamsa=(8.0*pterm/pow(parms[3],8))*pow(rab,7);
            gamsc=(8.0*pterm/pow(parms[3],8))*pow(rbc,7);
            break;
          default:
            fprintf(stderr, "Undefined Bond-Bend potential in routine 'CalculateFrameworkBondBendForce' ('framework_force.c')\n");
            exit(0);
            break;
        }

        // energy
        UHostBondBend[CurrentSystem]+=pterm;

        // forces
        fa.x=gamma*(Rbc.x-Rab.x*cost)/rab+gamsa*Rab.x;
        fa.y=gamma*(Rbc.y-Rab.y*cost)/rab+gamsa*Rab.y;
        fa.z=gamma*(Rbc.z-Rab.z*cost)/rab+gamsa*Rab.z;

        fc.x=gamma*(Rab.x-Rbc.x*cost)/rbc+gamsc*Rbc.x;
        fc.y=gamma*(Rab.y-Rbc.y*cost)/rbc+gamsc*Rbc.y;
        fc.z=gamma*(Rab.z-Rbc.z*cost)/rbc+gamsc*Rbc.z;

        // add contribution to the forces
        Framework[CurrentSystem].Atoms[f1][A].Force.x+=fa.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=fa.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=fa.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x-=(fa.x+fc.x);
        Framework[CurrentSystem].Atoms[f1][B].Force.y-=(fa.y+fc.y);
        Framework[CurrentSystem].Atoms[f1][B].Force.z-=(fa.z+fc.z);

        Framework[CurrentSystem].Atoms[f1][C].Force.x+=fc.x;
        Framework[CurrentSystem].Atoms[f1][C].Force.y+=fc.y;
        Framework[CurrentSystem].Atoms[f1][C].Force.z+=fc.z;

        // add contribution to the stress tensor
        StrainDerivativeTensor[CurrentSystem].ax-=(rab*Rab.x*fa.x+rbc*Rbc.x*fc.x);
        StrainDerivativeTensor[CurrentSystem].ay-=(rab*Rab.x*fa.y+rbc*Rbc.x*fc.y);
        StrainDerivativeTensor[CurrentSystem].az-=(rab*Rab.x*fa.z+rbc*Rbc.x*fc.z);

        StrainDerivativeTensor[CurrentSystem].bx-=(rab*Rab.y*fa.x+rbc*Rbc.y*fc.x);
        StrainDerivativeTensor[CurrentSystem].by-=(rab*Rab.y*fa.y+rbc*Rbc.y*fc.y);
        StrainDerivativeTensor[CurrentSystem].bz-=(rab*Rab.y*fa.z+rbc*Rbc.y*fc.z);

        StrainDerivativeTensor[CurrentSystem].cx-=(rab*Rab.z*fa.x+rbc*Rbc.z*fc.x);
        StrainDerivativeTensor[CurrentSystem].cy-=(rab*Rab.z*fa.y+rbc*Rbc.z*fc.y);
        StrainDerivativeTensor[CurrentSystem].cz-=(rab*Rab.z*fa.z+rbc*Rbc.z*fc.z);
      }
    }
  }
}

// Bend/Bend cross term for a centered atom B
// the first angle is A-B-C
// the second angle is A-B-D
void CalculateFrameworkBendBendForce(void)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL energy,rab,rbc,rbd;
  VECTOR Dab,Dbc,Dbd;
  REAL dot_abc,dot_abd;
  VECTOR fa,fb,fc,fd;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL DTheta1,DTheta2;
  REAL *parms;

  UHostBendBend[CurrentSystem]=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBendBends[f1];i++)
      {
        A=Framework[CurrentSystem].BendBends[f1][i].A;
        B=Framework[CurrentSystem].BendBends[f1][i].B;
        C=Framework[CurrentSystem].BendBends[f1][i].C;
        D=Framework[CurrentSystem].BendBends[f1][i].D;
        parms=(REAL*)&Framework[CurrentSystem].BendBendArguments[f1][i];

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
        SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

        dot_abd=Dab.x*Dbd.x+Dab.y*Dbd.y+Dab.z*Dbd.z;
        CosTheta2=dot_abd;
        CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
        Theta2=acos(CosTheta2);
        SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

        switch(Framework[CurrentSystem].BendBendType[f1][i])
        {
          case CVFF_BEND_BEND_CROSS:
          case CFF_BEND_BEND_CROSS:
            // p_0*(Theta1-p_1)*(Theta2-p_2)
            // ===================================
            // p_0/k_B [K/rad^2)]
            // p_1     [degrees]
            // p_2     [degrees]
            energy=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
            DTheta1=parms[0]*(Theta2-parms[2])/SinTheta1;
            DTheta2=parms[0]*(Theta1-parms[1])/SinTheta2;
            break;
          case MM3_BEND_BEND_CROSS:
            // -p_0*(Theta1-p_1)*(Theta2-p_2)
            // ===================================
            // p_0     [mdyne A/rad^2]
            // p_1     [degrees]
            // p_2     [degrees]
            energy=-parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])*(Theta2-parms[2]);
            DTheta1=-parms[0]*SQR(RAD2DEG)*(Theta2-parms[2])/SinTheta1;
            DTheta2=-parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])/SinTheta2;
            break;
          default:
            fprintf(stderr, "Undefined Bend-Bend potential in routine 'CalculateFrameworkBendBendForce' ('framework_force.c')\n");
            exit(0);
            break;
        }
        // energy
        UHostBendBend[CurrentSystem]+=energy;

        // forces bend
        fa.x=DTheta1*(Dbc.x-CosTheta1*Dab.x)/rab;
        fa.y=DTheta1*(Dbc.y-CosTheta1*Dab.y)/rab;
        fa.z=DTheta1*(Dbc.z-CosTheta1*Dab.z)/rab;

        fc.x=DTheta1*(Dab.x-CosTheta1*Dbc.x)/rbc;
        fc.y=DTheta1*(Dab.y-CosTheta1*Dbc.y)/rbc;
        fc.z=DTheta1*(Dab.z-CosTheta1*Dbc.z)/rbc;

        fb.x=-DTheta1*(fa.x+fc.x);
        fb.y=-DTheta1*(fa.y+fc.y);
        fb.z=-DTheta1*(fa.z+fc.z);

        Framework[CurrentSystem].Atoms[f1][A].Force.x+=fa.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=fa.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=fa.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x-=fa.x+fc.x;
        Framework[CurrentSystem].Atoms[f1][B].Force.y-=fa.y+fc.y;
        Framework[CurrentSystem].Atoms[f1][B].Force.z-=fa.z+fc.z;

        Framework[CurrentSystem].Atoms[f1][C].Force.x+=fc.x;
        Framework[CurrentSystem].Atoms[f1][C].Force.y+=fc.y;
        Framework[CurrentSystem].Atoms[f1][C].Force.z+=fc.z;

        StrainDerivativeTensor[CurrentSystem].ax-=rab*Dab.x*fa.x+rbc*Dbc.x*fc.x;
        StrainDerivativeTensor[CurrentSystem].bx-=rab*Dab.y*fa.x+rbc*Dbc.y*fc.x;
        StrainDerivativeTensor[CurrentSystem].cx-=rab*Dab.z*fa.x+rbc*Dbc.z*fc.x;

        StrainDerivativeTensor[CurrentSystem].ay-=rab*Dab.x*fa.y+rbc*Dbc.x*fc.y;
        StrainDerivativeTensor[CurrentSystem].by-=rab*Dab.y*fa.y+rbc*Dbc.y*fc.y;
        StrainDerivativeTensor[CurrentSystem].cy-=rab*Dab.z*fa.y+rbc*Dbc.z*fc.y;

        StrainDerivativeTensor[CurrentSystem].az-=rab*Dab.x*fa.z+rbc*Dbc.x*fc.z;
        StrainDerivativeTensor[CurrentSystem].bz-=rab*Dab.y*fa.z+rbc*Dbc.y*fc.z;
        StrainDerivativeTensor[CurrentSystem].cz-=rab*Dab.z*fa.z+rbc*Dbc.z*fc.z;

        fa.x=DTheta2*(Dbd.x-CosTheta2*Dab.x)/rab;
        fa.y=DTheta2*(Dbd.y-CosTheta2*Dab.y)/rab;
        fa.z=DTheta2*(Dbd.z-CosTheta2*Dab.z)/rab;

        fd.x=DTheta2*(Dab.x-CosTheta2*Dbd.x)/rbd;
        fd.y=DTheta2*(Dab.y-CosTheta2*Dbd.y)/rbd;
        fd.z=DTheta2*(Dab.z-CosTheta2*Dbd.z)/rbd;

        fb.x=-DTheta2*(fa.x+fd.x);
        fb.y=-DTheta2*(fa.y+fd.y);
        fb.z=-DTheta2*(fa.z+fd.z);

        Framework[CurrentSystem].Atoms[f1][A].Force.x+=fa.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=fa.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=fa.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x-=fa.x+fd.x;
        Framework[CurrentSystem].Atoms[f1][B].Force.y-=fa.y+fd.y;
        Framework[CurrentSystem].Atoms[f1][B].Force.z-=fa.z+fd.z;

        Framework[CurrentSystem].Atoms[f1][D].Force.x+=fd.x;
        Framework[CurrentSystem].Atoms[f1][D].Force.y+=fd.y;
        Framework[CurrentSystem].Atoms[f1][D].Force.z+=fd.z;

        StrainDerivativeTensor[CurrentSystem].ax-=rab*Dab.x*fa.x+rbd*Dbd.x*fd.x;
        StrainDerivativeTensor[CurrentSystem].bx-=rab*Dab.y*fa.x+rbd*Dbd.y*fd.x;
        StrainDerivativeTensor[CurrentSystem].cx-=rab*Dab.z*fa.x+rbd*Dbd.z*fd.x;

        StrainDerivativeTensor[CurrentSystem].ay-=rab*Dab.x*fa.y+rbd*Dbd.x*fd.y;
        StrainDerivativeTensor[CurrentSystem].by-=rab*Dab.y*fa.y+rbd*Dbd.y*fd.y;
        StrainDerivativeTensor[CurrentSystem].cy-=rab*Dab.z*fa.y+rbd*Dbd.z*fd.y;

        StrainDerivativeTensor[CurrentSystem].az-=rab*Dab.x*fa.z+rbd*Dbd.x*fd.z;
        StrainDerivativeTensor[CurrentSystem].bz-=rab*Dab.y*fa.z+rbd*Dbd.y*fd.z;
        StrainDerivativeTensor[CurrentSystem].cz-=rab*Dab.z*fa.z+rbd*Dbd.z*fd.z;
      }
    }
  }
}

// the Bend/Torsion cross term makes no contribution to the stress
void CalculateFrameworkBendTorsionForce(void)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL d,e, energy,rab,rbc,rcd;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2,DCos;
  VECTOR dtA,dtB,dtC,dtD,fa,fb,fc,fd,Pb,Pc;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL DTheta1,DTheta2,sign,Phi,SinPhi;
  REAL *parms;

  UHostBendTorsion[CurrentSystem]=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBendTorsions[f1];i++)
      {
        A=Framework[CurrentSystem].BendTorsions[f1][i].A;
        B=Framework[CurrentSystem].BendTorsions[f1][i].B;
        C=Framework[CurrentSystem].BendTorsions[f1][i].C;
        D=Framework[CurrentSystem].BendTorsions[f1][i].D;
        parms=(REAL*)&Framework[CurrentSystem].BendTorsionArguments[f1][i];

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
        SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

        dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;
        CosTheta2=-dot_cd/rcd;
        CosTheta2=MAX2(MIN2(CosTheta2,(REAL)1.0),-1.0);
        Theta2=acos(CosTheta2);
        SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

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

        switch(Framework[CurrentSystem].BendTorsionType[f1][i])
        {
          case CVFF_BEND_TORSION_CROSS:
          case CFF_BEND_TORSION_CROSS:
            // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
            // =====================================================================================
            // p_0/k_B [K/rad^3]
            // p_1     [degrees]
            // p_2     [degrees]
            energy=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi;
            DCos=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
            DTheta1=parms[0]*(Theta2-parms[2])*CosPhi/SinTheta1;
            DTheta2=parms[0]*(Theta1-parms[1])*CosPhi/SinTheta2;
            break;
          case SMOOTHED_DIHEDRAL:
            // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2)]*S(Theta2)
            // ======================================================================================
            // p_0/k_B [K/rad^2]
            // p_1     [-]
            // p_2     [degrees]
            Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
            Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
            Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
            Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
            Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
            Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
            sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                  +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
            Phi=SIGN(acos(CosPhi),sign);
            SinPhi=sin(Phi);
            SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
            energy=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2);
            DCos=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2)/SinPhi;
            DTheta1=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
            DTheta2=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*SmoothingDerivative(Theta2)/SinTheta2;
            break;
          case SMOOTHED_THREE_COSINE_DIHEDRAL:
            // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            energy=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                   Smoothing(Theta1)*Smoothing(Theta2);
            DCos=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0)*Smoothing(Theta1)*Smoothing(Theta2);
            DTheta1=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                     SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
            DTheta2=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                     SmoothingDerivative(Theta2)*Smoothing(Theta1)/SinTheta2;
            break;
          case NICHOLAS_DIHEDRAL:
            // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            energy=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                   Smoothing(Theta1);
            DCos=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0)*Smoothing(Theta1);
            DTheta1=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                    SmoothingDerivative(Theta1)/SinTheta1;
            DTheta2=0.0;
            break;
          case SMOOTHED_CFF_DIHEDRAL:
            // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            energy=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
            DCos=(-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
            DTheta1=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*
                    SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
            DTheta2=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*
                    Smoothing(Theta1)*SmoothingDerivative(Theta2)/SinTheta2;
            break;
          case SMOOTHED_CFF_DIHEDRAL2:
            // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            energy=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*Smoothing(Theta1)*Smoothing(Theta2);
            DCos=(parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi))*Smoothing(Theta1)*Smoothing(Theta2);
            DTheta1=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*
                    SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
            DTheta2=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*
                    Smoothing(Theta1)*SmoothingDerivative(Theta2)/SinTheta2;
            break;
          case SMOOTHED_CFF_BEND_TORSION_CROSS:
            // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
            // ======================================================================================
            // p_0/k_B [K/rad^3]
            // p_1     [degrees]
            // p_2     [degrees]
            energy=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi*Smoothing(Theta1)*Smoothing(Theta2);
            DCos=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*Smoothing(Theta1)*Smoothing(Theta2);
            DTheta1=CosPhi*parms[0]*(Theta2-parms[2])*Smoothing(Theta2)*(Smoothing(Theta1)+(Theta1-parms[1])*SmoothingDerivative(Theta1))/SinTheta1;
            DTheta2=CosPhi*parms[0]*(Theta1-parms[1])*Smoothing(Theta1)*(Smoothing(Theta2)+(Theta2-parms[2])*SmoothingDerivative(Theta2))/SinTheta2;
            break;
          default:
            fprintf(stderr, "Undefined Bend-Torsion potential in routine 'CalculateFrameworkBendTorsionForce' ('framework_force.c')\n");
            exit(0);
            break;
        }

        // energy
        UHostBendTorsion[CurrentSystem]+=energy;

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
        fa.x=-DCos*dtA.x;
        fa.y=-DCos*dtA.y;
        fa.z=-DCos*dtA.z;

        fb.x=-DCos*dtB.x;
        fb.y=-DCos*dtB.y;
        fb.z=-DCos*dtB.z;

        fc.x=-DCos*dtC.x;
        fc.y=-DCos*dtC.y;
        fc.z=-DCos*dtC.z;

        fd.x=-DCos*dtD.x;
        fd.y=-DCos*dtD.y;
        fd.z=-DCos*dtD.z;

        // forces torsion
        Framework[CurrentSystem].Atoms[f1][A].Force.x+=fa.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=fa.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=fa.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x+=fb.x;
        Framework[CurrentSystem].Atoms[f1][B].Force.y+=fb.y;
        Framework[CurrentSystem].Atoms[f1][B].Force.z+=fb.z;

        Framework[CurrentSystem].Atoms[f1][C].Force.x+=fc.x;
        Framework[CurrentSystem].Atoms[f1][C].Force.y+=fc.y;
        Framework[CurrentSystem].Atoms[f1][C].Force.z+=fc.z;

        Framework[CurrentSystem].Atoms[f1][D].Force.x+=fd.x;
        Framework[CurrentSystem].Atoms[f1][D].Force.y+=fd.y;
        Framework[CurrentSystem].Atoms[f1][D].Force.z+=fd.z;

        // add contribution to the stress tensor
        // Note: rbc is here because the vector was normalized before
        StrainDerivativeTensor[CurrentSystem].ax-=Dab.x*fa.x+Dbc.x*rbc*(fc.x+fd.x)+Dcd.x*fd.x;
        StrainDerivativeTensor[CurrentSystem].bx-=Dab.y*fa.x+Dbc.y*rbc*(fc.x+fd.x)+Dcd.y*fd.x;
        StrainDerivativeTensor[CurrentSystem].cx-=Dab.z*fa.x+Dbc.z*rbc*(fc.x+fd.x)+Dcd.z*fd.x;

        StrainDerivativeTensor[CurrentSystem].ay-=Dab.x*fa.y+Dbc.x*rbc*(fc.y+fd.y)+Dcd.x*fd.y;
        StrainDerivativeTensor[CurrentSystem].by-=Dab.y*fa.y+Dbc.y*rbc*(fc.y+fd.y)+Dcd.y*fd.y;
        StrainDerivativeTensor[CurrentSystem].cy-=Dab.z*fa.y+Dbc.z*rbc*(fc.y+fd.y)+Dcd.z*fd.y;

        StrainDerivativeTensor[CurrentSystem].az-=Dab.x*fa.z+Dbc.x*rbc*(fc.z+fd.z)+Dcd.x*fd.z;
        StrainDerivativeTensor[CurrentSystem].bz-=Dab.y*fa.z+Dbc.y*rbc*(fc.z+fd.z)+Dcd.y*fd.z;
        StrainDerivativeTensor[CurrentSystem].cz-=Dab.z*fa.z+Dbc.z*rbc*(fc.z+fd.z)+Dcd.z*fd.z;

        Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;
        Dcd.x/=rcd; Dcd.y/=rcd; Dcd.z/=rcd;

        // forces bends
        fa.x=DTheta1*(Dbc.x-Dab.x*CosTheta1)/rab;
        fa.y=DTheta1*(Dbc.y-Dab.y*CosTheta1)/rab;
        fa.z=DTheta1*(Dbc.z-Dab.z*CosTheta1)/rab;

        fb.x=DTheta2*(Dcd.x+Dbc.x*CosTheta2)/rbc;
        fb.y=DTheta2*(Dcd.y+Dbc.y*CosTheta2)/rbc;
        fb.z=DTheta2*(Dcd.z+Dbc.z*CosTheta2)/rbc;

        fc.x=DTheta1*(Dab.x-Dbc.x*CosTheta1)/rbc;
        fc.y=DTheta1*(Dab.y-Dbc.y*CosTheta1)/rbc;
        fc.z=DTheta1*(Dab.z-Dbc.z*CosTheta1)/rbc;

        fd.x=DTheta2*(-Dbc.x-Dcd.x*CosTheta2)/rcd;
        fd.y=DTheta2*(-Dbc.y-Dcd.y*CosTheta2)/rcd;
        fd.z=DTheta2*(-Dbc.z-Dcd.z*CosTheta2)/rcd;

        Framework[CurrentSystem].Atoms[f1][A].Force.x+=fa.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=fa.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=fa.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x+=fb.x-fa.x-fc.x;
        Framework[CurrentSystem].Atoms[f1][B].Force.y+=fb.y-fa.y-fc.y;
        Framework[CurrentSystem].Atoms[f1][B].Force.z+=fb.z-fa.z-fc.z;

        Framework[CurrentSystem].Atoms[f1][C].Force.x+=fc.x-fb.x-fd.x;
        Framework[CurrentSystem].Atoms[f1][C].Force.y+=fc.y-fb.y-fd.y;
        Framework[CurrentSystem].Atoms[f1][C].Force.z+=fc.z-fb.z-fd.z;

        Framework[CurrentSystem].Atoms[f1][D].Force.x+=fd.x;
        Framework[CurrentSystem].Atoms[f1][D].Force.y+=fd.y;
        Framework[CurrentSystem].Atoms[f1][D].Force.z+=fd.z;

        StrainDerivativeTensor[CurrentSystem].ax-=rab*Dab.x*fa.x+rbc*Dbc.x*fc.x+rbc*Dbc.x*fb.x+rcd*Dcd.x*fd.x;
        StrainDerivativeTensor[CurrentSystem].bx-=rab*Dab.y*fa.x+rbc*Dbc.y*fc.x+rbc*Dbc.y*fb.x+rcd*Dcd.y*fd.x;
        StrainDerivativeTensor[CurrentSystem].cx-=rab*Dab.z*fa.x+rbc*Dbc.z*fc.x+rbc*Dbc.z*fb.x+rcd*Dcd.z*fd.x;

        StrainDerivativeTensor[CurrentSystem].ay-=rab*Dab.x*fa.y+rbc*Dbc.x*fc.y+rbc*Dbc.x*fb.y+rcd*Dcd.x*fd.y;
        StrainDerivativeTensor[CurrentSystem].by-=rab*Dab.y*fa.y+rbc*Dbc.y*fc.y+rbc*Dbc.y*fb.y+rcd*Dcd.y*fd.y;
        StrainDerivativeTensor[CurrentSystem].cy-=rab*Dab.z*fa.y+rbc*Dbc.z*fc.y+rbc*Dbc.z*fb.y+rcd*Dcd.z*fd.y;

        StrainDerivativeTensor[CurrentSystem].az-=rab*Dab.x*fa.z+rbc*Dbc.x*fc.z+rbc*Dbc.x*fb.z+rcd*Dcd.x*fd.z;
        StrainDerivativeTensor[CurrentSystem].bz-=rab*Dab.y*fa.z+rbc*Dbc.y*fc.z+rbc*Dbc.y*fb.z+rcd*Dcd.y*fd.z;
        StrainDerivativeTensor[CurrentSystem].cz-=rab*Dab.z*fa.z+rbc*Dbc.z*fc.z+rbc*Dbc.z*fb.z+rcd*Dcd.z*fd.z;
      }
    }
  }
}

void CalculateFrameworkBondTorsionForce(void)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL d,e, energy,rab,rbc,rcd;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2,DCos,temp;
  VECTOR dtA,dtB,dtC,dtD,fa,fb,fc,fd;
  REAL *parms,gamsa,gamsb,gamsc;

  UHostBondTorsion[CurrentSystem]=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondTorsions[f1];i++)
      {
        A=Framework[CurrentSystem].BondTorsions[f1][i].A;
        B=Framework[CurrentSystem].BondTorsions[f1][i].B;
        C=Framework[CurrentSystem].BondTorsions[f1][i].C;
        D=Framework[CurrentSystem].BondTorsions[f1][i].D;
        parms=(REAL*)&Framework[CurrentSystem].BondTorsionArguments[f1][i];

        posA=Framework[CurrentSystem].Atoms[f1][A].Position;
        posB=Framework[CurrentSystem].Atoms[f1][B].Position;
        posC=Framework[CurrentSystem].Atoms[f1][C].Position;
        posD=Framework[CurrentSystem].Atoms[f1][C].Position;

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

        switch(Framework[CurrentSystem].BondTorsionType[f1][i])
        {
          case MM3_BOND_TORSION_CROSS:
            // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
            // =====================================================================================
            // p_0     [kcal/A mole]
            // p_1     [kcal/A mole]
            // p_2     [kcal/A mole]
            // p_3     [A]
            temp=(rbc-parms[3]);
            energy=parms[0]*temp*CosPhi+parms[1]*temp*(2.0*CosPhi2-1.0)+parms[2]*temp*(4.0*CosPhi2*CosPhi-3.0*CosPhi);
            DCos=parms[0]*temp+4.0*parms[1]*temp*CosPhi+parms[2]*temp*(12.0*CosPhi2-3.0);
            gamsa=0.0;
            gamsb=-(parms[0]*CosPhi+parms[1]*(2.0*CosPhi2-1.0)+parms[2]*(4.0*CosPhi2*CosPhi-3.0*CosPhi));
            gamsc=0.0;
            break;
          default:
            fprintf(stderr, "Undefined Bond-Torsion potential in routine 'CalculateFrameworkBondTorsionForce' ('framework_force.c')\n");
            exit(0);
            break;
        }

        // energy
        UHostBondTorsion[CurrentSystem]+=energy;

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

        Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;
        Ddc.x/=rcd; Ddc.y/=rcd; Ddc.z/=rcd;

        fa.x=DCos*dtA.x-gamsa*Dab.x;
        fa.y=DCos*dtA.y-gamsa*Dab.y;
        fa.z=DCos*dtA.z-gamsa*Dab.z;

        fb.x=DCos*dtB.x+gamsb*Dcb.x+gamsa*Dab.x;
        fb.y=DCos*dtB.y+gamsb*Dcb.y+gamsa*Dab.y;
        fb.z=DCos*dtB.z+gamsb*Dcb.z+gamsa*Dab.z;

        fc.x=DCos*dtC.x-gamsb*Dcb.x+gamsc*Ddc.x;
        fc.y=DCos*dtC.y-gamsb*Dcb.y+gamsc*Ddc.y;
        fc.z=DCos*dtC.z-gamsb*Dcb.z+gamsc*Ddc.z;

        fd.x=DCos*dtD.x-gamsc*Ddc.x;
        fd.y=DCos*dtD.y-gamsc*Ddc.y;
        fd.z=DCos*dtD.z-gamsc*Ddc.z;

        // forces torsion
        Framework[CurrentSystem].Atoms[f1][A].Force.x-=fa.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y-=fa.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z-=fa.z;

        Framework[CurrentSystem].Atoms[f1][B].Force.x-=fb.x;
        Framework[CurrentSystem].Atoms[f1][B].Force.y-=fb.y;
        Framework[CurrentSystem].Atoms[f1][B].Force.z-=fb.z;

        Framework[CurrentSystem].Atoms[f1][C].Force.x-=fc.x;
        Framework[CurrentSystem].Atoms[f1][C].Force.y-=fc.y;
        Framework[CurrentSystem].Atoms[f1][C].Force.z-=fc.z;

        Framework[CurrentSystem].Atoms[f1][D].Force.x-=fd.x;
        Framework[CurrentSystem].Atoms[f1][D].Force.y-=fd.y;
        Framework[CurrentSystem].Atoms[f1][D].Force.z-=fd.z;

        // add contribution to the stress tensor
        StrainDerivativeTensor[CurrentSystem].ax+=rbc*Dcb.x*(fc.x+fd.x)+rab*Dab.x*fa.x+rcd*Ddc.x*fd.x;
        StrainDerivativeTensor[CurrentSystem].ay+=rbc*Dcb.x*(fc.y+fd.y)+rab*Dab.x*fa.y+rcd*Ddc.x*fd.y;
        StrainDerivativeTensor[CurrentSystem].az+=rbc*Dcb.x*(fc.z+fd.z)+rab*Dab.x*fa.z+rcd*Ddc.x*fd.z;

        StrainDerivativeTensor[CurrentSystem].bx+=rbc*Dcb.y*(fc.x+fd.x)+rab*Dab.y*fa.x+rcd*Ddc.y*fd.x;
        StrainDerivativeTensor[CurrentSystem].by+=rbc*Dcb.y*(fc.y+fd.y)+rab*Dab.y*fa.y+rcd*Ddc.y*fd.y;
        StrainDerivativeTensor[CurrentSystem].bz+=rbc*Dcb.y*(fc.z+fd.z)+rab*Dab.y*fa.z+rcd*Ddc.y*fd.z;

        StrainDerivativeTensor[CurrentSystem].cx+=rbc*Dcb.z*(fc.x+fd.x)+rab*Dab.z*fa.x+rcd*Ddc.z*fd.x;
        StrainDerivativeTensor[CurrentSystem].cy+=rbc*Dcb.z*(fc.y+fd.y)+rab*Dab.z*fa.y+rcd*Ddc.z*fd.y;
        StrainDerivativeTensor[CurrentSystem].cz+=rbc*Dcb.z*(fc.z+fd.z)+rab*Dab.z*fa.z+rcd*Ddc.z*fd.z;
      }
    }
  }
}


int CalculateFrameworkIntraVDWForce(void)
{
  int i,j,typeA,typeB;
  REAL energy,force_factor;
  REAL rr,ReductionA,ReductionB;
  VECTOR pos,posA,posB,dr,f;
  VECTOR posA1,posA2,posB1,posB2;
  VECTOR drA,drB,fa,fb,v,w;
  int ConnectedAtomA1,ConnectedAtomA2;
  int ConnectedAtomB1,ConnectedAtomB2;
  REAL ra,rb,length_v,length_w,dot_product;
  int f1,f2,A,B;
  REAL *parms;

  // Framework-Framework energy
  UHostHostVDW[CurrentSystem]=0.0;

  if(!InternalFrameworkLennardJonesInteractions) return 0;

  // contributions from intra-framework
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

        ReductionA=0.0;
        ConnectedAtomA1=ConnectedAtomA2=-1;
        ra=drA.x=drA.y=drA.z=0.0;
        v.x=v.y=v.z=length_v=0.0;
        if(PseudoAtoms[typeA].AnisotropicCorrection)
        {
          switch(Framework[CurrentSystem].Connectivity[f1][i])
          {
            case 0:
              break;
            case 1:
              ConnectedAtomA1=Framework[CurrentSystem].Neighbours[f1][i][0];
              if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
                ReductionA=1.0+PseudoAtoms[typeA].AnisotropicDisplacement;
              else
              {
                posA1=Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Position;
                pos=Framework[CurrentSystem].Atoms[f1][i].Position;
                drA.x=posA1.x-pos.x;
                drA.y=posA1.y-pos.y;
                drA.z=posA1.z-pos.z;
                drA=ApplyBoundaryCondition(drA);
                ra=sqrt(SQR(drA.x)+SQR(drA.y)+SQR(drA.z));
                ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
              }
              break;
            case 2:
              switch(Framework[CurrentSystem].AnisotropicType)
              {
                case ANISOTROPIC_BISECTION:
                  fprintf(stderr, "ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                  exit(0);
                  break;
                case ANISOTROPIC_MID_POINT:
                  ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
                  ConnectedAtomA1=Framework[CurrentSystem].Neighbours[f1][i][0];
                  ConnectedAtomA2=Framework[CurrentSystem].Neighbours[f1][i][1];
                  pos=Framework[CurrentSystem].Atoms[f1][i].Position;
                  posA1=Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Position;
                  posA2=Framework[CurrentSystem].Atoms[f1][ConnectedAtomA2].Position;
                  v.x=pos.x-0.5*(posA1.x+posA2.x);
                  v.y=pos.y-0.5*(posA1.y+posA2.y);
                  v.z=pos.z-0.5*(posA1.z+posA2.z);
                  v=ApplyBoundaryCondition(v);
                  length_v=sqrt(SQR(v.x)+SQR(v.y)+SQR(v.z));
                  break;
                default:
                  fprintf(stderr, "ERROR: unknown Anisotropy type in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                  exit(0);
                  break;
              }
              break;
            default:
              break;
          }
        }

        for(j=i+1;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],0))
          {
            typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
            posB=Framework[CurrentSystem].Atoms[f1][j].AnisotropicPosition;

            ReductionB=0.0;
            ConnectedAtomB1=ConnectedAtomB2=-1;
            rb=drB.x=drB.y=drB.z=0.0;
            w.x=w.y=w.z=length_w=0.0;
            if(PseudoAtoms[typeB].AnisotropicCorrection)
            {
              switch(Framework[CurrentSystem].Connectivity[f1][j])
              {
                case 0:
                  break;
                case 1:
                  ConnectedAtomB1=Framework[CurrentSystem].Neighbours[f1][j][0];
                  if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                    ReductionB=1.0+PseudoAtoms[typeB].AnisotropicDisplacement;
                  else
                  {
                    posB1=Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Position;
                    pos=Framework[CurrentSystem].Atoms[f1][j].Position;
                    drB.x=posB1.x-pos.x;
                    drB.y=posB1.y-pos.y;
                    drB.z=posB1.z-pos.z;
                    drB=ApplyBoundaryCondition(drB);
                    rb=sqrt(SQR(drB.x)+SQR(drB.y)+SQR(drB.z));
                    ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                  }
                  break;
                case 2:
                  switch(Framework[CurrentSystem].AnisotropicType)
                  {
                    case ANISOTROPIC_BISECTION:
                      fprintf(stderr, "ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                      exit(0);
                      break;
                    case ANISOTROPIC_MID_POINT:
                      ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                      ConnectedAtomB1=Framework[CurrentSystem].Neighbours[f1][j][0];
                      ConnectedAtomB2=Framework[CurrentSystem].Neighbours[f1][j][1];
                      pos=Framework[CurrentSystem].Atoms[f1][j].Position;
                      posB1=Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Position;
                      posB2=Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Position;
                      w.x=pos.x-0.5*(posB1.x+posB2.x);
                      w.y=pos.y-0.5*(posB1.y+posB2.y);
                      w.z=pos.z-0.5*(posB1.z+posB2.z);
                      w=ApplyBoundaryCondition(w);
                      length_w=sqrt(SQR(w.x)+SQR(w.y)+SQR(w.z));
                      break;
                    default:
                      fprintf(stderr, "ERROR: unknown Anisotropy type in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                      exit(0);
                      break;
                  }
                  break;
                default:
                  break;
              }
            }

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);

              UHostHostVDW[CurrentSystem]+=energy;

              f.x=force_factor*dr.x;
              f.y=force_factor*dr.y;
              f.z=force_factor*dr.z;

              if(PseudoAtoms[typeA].AnisotropicCorrection)
              {
                switch(Framework[CurrentSystem].Connectivity[f1][i])
                {
                  case 0:
                    break;
                  case 1:
                    if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
                    {
                      Framework[CurrentSystem].Atoms[f1][i].Force.x-=ReductionA*f.x;
                      Framework[CurrentSystem].Atoms[f1][i].Force.y-=ReductionA*f.y;
                      Framework[CurrentSystem].Atoms[f1][i].Force.z-=ReductionA*f.z;

                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.x-=(1.0-ReductionA)*f.x;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.y-=(1.0-ReductionA)*f.y;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.z-=(1.0-ReductionA)*f.z;
                    }
                    else
                    {
                      dot_product=dr.x*drA.x+dr.y*drA.y+dr.z*drA.z;

                      fa.x=(-(ReductionA*drA.x*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.x))*force_factor;
                      fa.y=(-(ReductionA*drA.y*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.y))*force_factor;
                      fa.z=(-(ReductionA*drA.z*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.z))*force_factor;

                      fb.x=(ReductionA*drA.x*dot_product/CUBE(ra)-(ReductionA/ra)*dr.x)*force_factor;
                      fb.y=(ReductionA*drA.y*dot_product/CUBE(ra)-(ReductionA/ra)*dr.y)*force_factor;
                      fb.z=(ReductionA*drA.z*dot_product/CUBE(ra)-(ReductionA/ra)*dr.z)*force_factor;

                      Framework[CurrentSystem].Atoms[f1][i].Force.x-=fa.x;
                      Framework[CurrentSystem].Atoms[f1][i].Force.y-=fa.y;
                      Framework[CurrentSystem].Atoms[f1][i].Force.z-=fa.z;

                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.x-=fb.x;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.y-=fb.y;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.z-=fb.z;
                    }
                    break;
                  case 2:
                    dot_product=v.x*dr.x+v.y*dr.y+v.z*dr.z;

                    fa.x=0.5*(ReductionA*v.x*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.x)*force_factor;
                    fa.y=0.5*(ReductionA*v.y*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.y)*force_factor;
                    fa.z=0.5*(ReductionA*v.z*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.z)*force_factor;

                    fb.x=(-ReductionA*v.x*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.x)*force_factor;
                    fb.y=(-ReductionA*v.y*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.y)*force_factor;
                    fb.z=(-ReductionA*v.z*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.z)*force_factor;

                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.x-=fa.x;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.y-=fa.y;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.z-=fa.z;

                    Framework[CurrentSystem].Atoms[f1][i].Force.x-=fb.x;
                    Framework[CurrentSystem].Atoms[f1][i].Force.y-=fb.y;
                    Framework[CurrentSystem].Atoms[f1][i].Force.z-=fb.z;

                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA2].Force.x-=fa.x;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA2].Force.y-=fa.y;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA2].Force.z-=fa.z;
                    break;
                  default:
                    fprintf(stderr, "ERROR: Not yet implemented in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                    exit(0);
                    break;
                }
              }
              else
              {
                Framework[CurrentSystem].Atoms[f1][i].Force.x-=f.x;
                Framework[CurrentSystem].Atoms[f1][i].Force.y-=f.y;
                Framework[CurrentSystem].Atoms[f1][i].Force.z-=f.z;
              }

              if(PseudoAtoms[typeB].AnisotropicCorrection)
              {
                switch(Framework[CurrentSystem].Connectivity[f1][j])
                {
                  case 0:
                    break;
                  case 1:
                    if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                    {
                      Framework[CurrentSystem].Atoms[f1][j].Force.x+=ReductionB*f.x;
                      Framework[CurrentSystem].Atoms[f1][j].Force.y+=ReductionB*f.y;
                      Framework[CurrentSystem].Atoms[f1][j].Force.z+=ReductionB*f.z;

                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.x+=(1.0-ReductionB)*f.x;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.y+=(1.0-ReductionB)*f.y;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.z+=(1.0-ReductionB)*f.z;
                    }
                    else
                    {
                      dot_product=dr.x*drB.x+dr.y*drB.y+dr.z*drB.z;

                      fa.x=(-(ReductionB*drA.x*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.x))*force_factor;
                      fa.y=(-(ReductionB*drA.y*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.y))*force_factor;
                      fa.z=(-(ReductionB*drA.z*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.z))*force_factor;

                      fb.x=(ReductionB*drB.x*dot_product/CUBE(rb)-(ReductionB/rb)*dr.x)*force_factor;
                      fb.y=(ReductionB*drB.y*dot_product/CUBE(rb)-(ReductionB/rb)*dr.y)*force_factor;
                      fb.z=(ReductionB*drB.z*dot_product/CUBE(rb)-(ReductionB/rb)*dr.z)*force_factor;

                      Framework[CurrentSystem].Atoms[f1][j].Force.x+=fa.x;
                      Framework[CurrentSystem].Atoms[f1][j].Force.y+=fa.y;
                      Framework[CurrentSystem].Atoms[f1][j].Force.z+=fa.z;

                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.x+=fb.x;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.y+=fb.y;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.z+=fb.z;
                    }
                    break;
                  case 2:
                    dot_product=w.x*dr.x+w.y*dr.y+w.z*dr.z;

                    fa.x=0.5*(ReductionB*w.x*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.x)*force_factor;
                    fa.y=0.5*(ReductionB*w.y*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.y)*force_factor;
                    fa.z=0.5*(ReductionB*w.z*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.z)*force_factor;

                    fb.x=(-ReductionB*w.x*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.x)*force_factor;
                    fb.y=(-ReductionB*w.y*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.y)*force_factor;
                    fb.z=(-ReductionB*w.z*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.z)*force_factor;

                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.x+=fa.x;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.y+=fa.y;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.z+=fa.z;

                    Framework[CurrentSystem].Atoms[f1][j].Force.x+=fb.x;
                    Framework[CurrentSystem].Atoms[f1][j].Force.y+=fb.y;
                    Framework[CurrentSystem].Atoms[f1][j].Force.z+=fb.z;

                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Force.x+=fa.x;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Force.y+=fa.y;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Force.z+=fa.z;
                    break;
                  default:
                    fprintf(stderr, "ERROR: Not yet implemented in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                    exit(0);
                    break;
                }
              }
              else
              {
                Framework[CurrentSystem].Atoms[f1][j].Force.x+=f.x;
                Framework[CurrentSystem].Atoms[f1][j].Force.y+=f.y;
                Framework[CurrentSystem].Atoms[f1][j].Force.z+=f.z;
              }

              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
            }
          }
        }
      }
    }
  }

  // contributions from interactions between the frameworks
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1+1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

        ReductionA=0.0;
        ConnectedAtomA1=ConnectedAtomA2=-1;
        ra=drA.x=drA.y=drA.z=0.0;
        v.x=v.y=v.z=length_v=0.0;
        if(PseudoAtoms[typeA].AnisotropicCorrection)
        {
          switch(Framework[CurrentSystem].Connectivity[f1][i])
          {
            case 0:
              break;
            case 1:
              ConnectedAtomA1=Framework[CurrentSystem].Neighbours[f1][i][0];
              if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
                ReductionA=1.0+PseudoAtoms[typeA].AnisotropicDisplacement;
              else
              {
                posA1=Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Position;
                pos=Framework[CurrentSystem].Atoms[f1][i].Position;
                drA.x=posA1.x-pos.x;
                drA.y=posA1.y-pos.y;
                drA.z=posA1.z-pos.z;
                drA=ApplyBoundaryCondition(drA);
                ra=sqrt(SQR(drA.x)+SQR(drA.y)+SQR(drA.z));
                ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
              }
              break;
            case 2:
              switch(Framework[CurrentSystem].AnisotropicType)
              {
                case ANISOTROPIC_BISECTION:
                  fprintf(stderr, "ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                  exit(0);
                  break;
                case ANISOTROPIC_MID_POINT:
                  ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
                  ConnectedAtomA1=Framework[CurrentSystem].Neighbours[f1][i][0];
                  ConnectedAtomA2=Framework[CurrentSystem].Neighbours[f1][i][1];
                  pos=Framework[CurrentSystem].Atoms[f1][i].Position;
                  posA1=Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Position;
                  posA2=Framework[CurrentSystem].Atoms[f1][ConnectedAtomA2].Position;
                  v.x=pos.x-0.5*(posA1.x+posA2.x);
                  v.y=pos.y-0.5*(posA1.y+posA2.y);
                  v.z=pos.z-0.5*(posA1.z+posA2.z);
                  v=ApplyBoundaryCondition(v);
                  length_v=sqrt(SQR(v.x)+SQR(v.y)+SQR(v.z));
                  break;
                default:
                  fprintf(stderr, "ERROR: unknown Anisotropy type in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                  exit(0);
                  break;
              }
              break;
            default:
              break;
          }
        }


        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
        {
          typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
          posB=Framework[CurrentSystem].Atoms[f2][j].AnisotropicPosition;

          ReductionB=0.0;
          ConnectedAtomB1=ConnectedAtomB2=-1;
          rb=drB.x=drB.y=drB.z=0.0;
          w.x=w.y=w.z=length_w=0.0;
          if(PseudoAtoms[typeB].AnisotropicCorrection)
          {
            switch(Framework[CurrentSystem].Connectivity[f2][j])
            {
              case 0:
                break;
              case 1:
                ConnectedAtomB1=Framework[CurrentSystem].Neighbours[f2][j][0];
                if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                  ReductionB=1.0+PseudoAtoms[typeB].AnisotropicDisplacement;
                else
                {
                  posB1=Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Position;
                  pos=Framework[CurrentSystem].Atoms[f2][j].Position;
                  drB.x=posB1.x-pos.x;
                  drB.y=posB1.y-pos.y;
                  drB.z=posB1.z-pos.z;
                  drB=ApplyBoundaryCondition(drB);
                  rb=sqrt(SQR(drB.x)+SQR(drB.y)+SQR(drB.z));
                  ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                }
                break;
              case 2:
                switch(Framework[CurrentSystem].AnisotropicType)
                {
                  case ANISOTROPIC_BISECTION:
                    fprintf(stderr, "ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                    exit(0);
                    break;
                  case ANISOTROPIC_MID_POINT:
                    ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                    ConnectedAtomB1=Framework[CurrentSystem].Neighbours[f2][j][0];
                    ConnectedAtomB2=Framework[CurrentSystem].Neighbours[f2][j][1];
                    pos=Framework[CurrentSystem].Atoms[f2][j].Position;
                    posB1=Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Position;
                    posB2=Framework[CurrentSystem].Atoms[f2][ConnectedAtomB2].Position;
                    w.x=pos.x-0.5*(posB1.x+posB2.x);
                    w.y=pos.y-0.5*(posB1.y+posB2.y);
                    w.z=pos.z-0.5*(posB1.z+posB2.z);
                    w=ApplyBoundaryCondition(w);
                    length_w=sqrt(SQR(w.x)+SQR(w.y)+SQR(w.z));
                    break;
                  default:
                    fprintf(stderr, "ERROR: unknown Anisotropy type in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                    exit(0);
                    break;
                }
                break;
              default:
                break;
            }
          }

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffVDWSquared)
          {
            PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);

            UHostHostVDW[CurrentSystem]+=energy;

            f.x=force_factor*dr.x;
            f.y=force_factor*dr.y;
            f.z=force_factor*dr.z;

            if(PseudoAtoms[typeA].AnisotropicCorrection)
            {
              switch(Framework[CurrentSystem].Connectivity[f1][i])
              {
                case 0:
                  break;
                case 1:
                  if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
                  {
                    Framework[CurrentSystem].Atoms[f1][i].Force.x-=ReductionA*f.x;
                    Framework[CurrentSystem].Atoms[f1][i].Force.y-=ReductionA*f.y;
                    Framework[CurrentSystem].Atoms[f1][i].Force.z-=ReductionA*f.z;

                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.x-=(1.0-ReductionA)*f.x;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.y-=(1.0-ReductionA)*f.y;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.z-=(1.0-ReductionA)*f.z;
                  }
                  else
                  {
                    dot_product=dr.x*drA.x+dr.y*drA.y+dr.z*drA.z;

                    fa.x=(-(ReductionA*drA.x*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.x))*force_factor;
                    fa.y=(-(ReductionA*drA.y*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.y))*force_factor;
                    fa.z=(-(ReductionA*drA.z*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.z))*force_factor;

                    fb.x=(ReductionA*drA.x*dot_product/CUBE(ra)-(ReductionA/ra)*dr.x)*force_factor;
                    fb.y=(ReductionA*drA.y*dot_product/CUBE(ra)-(ReductionA/ra)*dr.y)*force_factor;
                    fb.z=(ReductionA*drA.z*dot_product/CUBE(ra)-(ReductionA/ra)*dr.z)*force_factor;

                    Framework[CurrentSystem].Atoms[f1][i].Force.x-=fa.x;
                    Framework[CurrentSystem].Atoms[f1][i].Force.y-=fa.y;
                    Framework[CurrentSystem].Atoms[f1][i].Force.z-=fa.z;

                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.x-=fb.x;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.y-=fb.y;
                    Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.z-=fb.z;
                  }
                  break;
                case 2:
                  dot_product=v.x*dr.x+v.y*dr.y+v.z*dr.z;

                  fa.x=0.5*(ReductionA*v.x*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.x)*force_factor;
                  fa.y=0.5*(ReductionA*v.y*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.y)*force_factor;
                  fa.z=0.5*(ReductionA*v.z*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.z)*force_factor;

                  fb.x=(-ReductionA*v.x*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.x)*force_factor;
                  fb.y=(-ReductionA*v.y*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.y)*force_factor;
                  fb.z=(-ReductionA*v.z*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.z)*force_factor;

                  Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.x-=fa.x;
                  Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.y-=fa.y;
                  Framework[CurrentSystem].Atoms[f1][ConnectedAtomA1].Force.z-=fa.z;

                  Framework[CurrentSystem].Atoms[f1][i].Force.x-=fb.x;
                  Framework[CurrentSystem].Atoms[f1][i].Force.y-=fb.y;
                  Framework[CurrentSystem].Atoms[f1][i].Force.z-=fb.z;

                  Framework[CurrentSystem].Atoms[f1][ConnectedAtomA2].Force.x-=fa.x;
                  Framework[CurrentSystem].Atoms[f1][ConnectedAtomA2].Force.y-=fa.y;
                  Framework[CurrentSystem].Atoms[f1][ConnectedAtomA2].Force.z-=fa.z;
                  break;
                default:
                  fprintf(stderr, "ERROR: Not yet implemented in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                  exit(0);
                  break;
              }
            }
            else
            {
              Framework[CurrentSystem].Atoms[f1][i].Force.x-=f.x;
              Framework[CurrentSystem].Atoms[f1][i].Force.y-=f.y;
              Framework[CurrentSystem].Atoms[f1][i].Force.z-=f.z;
            }

            if(PseudoAtoms[typeB].AnisotropicCorrection)
            {
              switch(Framework[CurrentSystem].Connectivity[f2][j])
              {
                case 0:
                  break;
                case 1:
                  if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                  {
                    Framework[CurrentSystem].Atoms[f2][j].Force.x+=ReductionB*f.x;
                    Framework[CurrentSystem].Atoms[f2][j].Force.y+=ReductionB*f.y;
                    Framework[CurrentSystem].Atoms[f2][j].Force.z+=ReductionB*f.z;

                    Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Force.x+=(1.0-ReductionB)*f.x;
                    Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Force.y+=(1.0-ReductionB)*f.y;
                    Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Force.z+=(1.0-ReductionB)*f.z;
                  }
                  else
                  {
                    dot_product=dr.x*drB.x+dr.y*drB.y+dr.z*drB.z;

                    fa.x=(-(ReductionB*drA.x*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.x))*force_factor;
                    fa.y=(-(ReductionB*drA.y*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.y))*force_factor;
                    fa.z=(-(ReductionB*drA.z*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.z))*force_factor;

                    fb.x=(ReductionB*drB.x*dot_product/CUBE(rb)-(ReductionB/rb)*dr.x)*force_factor;
                    fb.y=(ReductionB*drB.y*dot_product/CUBE(rb)-(ReductionB/rb)*dr.y)*force_factor;
                    fb.z=(ReductionB*drB.z*dot_product/CUBE(rb)-(ReductionB/rb)*dr.z)*force_factor;

                    Framework[CurrentSystem].Atoms[f2][j].Force.x+=fa.x;
                    Framework[CurrentSystem].Atoms[f2][j].Force.y+=fa.y;
                    Framework[CurrentSystem].Atoms[f2][j].Force.z+=fa.z;

                    Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Force.x+=fb.x;
                    Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Force.y+=fb.y;
                    Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Force.z+=fb.z;
                  }
                  break;
                case 2:
                  dot_product=w.x*dr.x+w.y*dr.y+w.z*dr.z;

                  fa.x=0.5*(ReductionB*w.x*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.x)*force_factor;
                  fa.y=0.5*(ReductionB*w.y*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.y)*force_factor;
                  fa.z=0.5*(ReductionB*w.z*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.z)*force_factor;

                  fb.x=(-ReductionB*w.x*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.x)*force_factor;
                  fb.y=(-ReductionB*w.y*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.y)*force_factor;
                  fb.z=(-ReductionB*w.z*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.z)*force_factor;

                  Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Force.x+=fa.x;
                  Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Force.y+=fa.y;
                  Framework[CurrentSystem].Atoms[f2][ConnectedAtomB1].Force.z+=fa.z;

                  Framework[CurrentSystem].Atoms[f2][j].Force.x+=fb.x;
                  Framework[CurrentSystem].Atoms[f2][j].Force.y+=fb.y;
                  Framework[CurrentSystem].Atoms[f2][j].Force.z+=fb.z;

                  Framework[CurrentSystem].Atoms[f2][ConnectedAtomB2].Force.x+=fa.x;
                  Framework[CurrentSystem].Atoms[f2][ConnectedAtomB2].Force.y+=fa.y;
                  Framework[CurrentSystem].Atoms[f2][ConnectedAtomB2].Force.z+=fa.z;
                  break;
                default:
                  fprintf(stderr, "ERROR: Not yet implemented in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                  exit(0);
                  break;
              }
            }
            else
            {
              Framework[CurrentSystem].Atoms[f2][j].Force.x+=f.x;
              Framework[CurrentSystem].Atoms[f2][j].Force.y+=f.y;
              Framework[CurrentSystem].Atoms[f2][j].Force.z+=f.z;
            }

            StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
          }
        }
      }
    }
  }

  // contributions from interactions 1-4 torsions
  // TODO: fix to handle anisotropic sites
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[f1];i++)
      {
        parms=Framework[CurrentSystem].TorsionArguments[f1][i];
        if(fabs(parms[6])>1e-8)
        {
          A=Framework[CurrentSystem].Torsions[f1][i].A;
          typeA=Framework[CurrentSystem].Atoms[f1][A].Type;
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;

          B=Framework[CurrentSystem].Torsions[f1][i].D;
          typeB=Framework[CurrentSystem].Atoms[f1][B].Type;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;


          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);

          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);

          UHostHostVDW[CurrentSystem]+=parms[6]*energy;

          force_factor=parms[6]*force_factor;

          StrainDerivativeTensor[CurrentSystem].ax+=force_factor*dr.x*dr.x;
          StrainDerivativeTensor[CurrentSystem].bx+=force_factor*dr.y*dr.x;
          StrainDerivativeTensor[CurrentSystem].cx+=force_factor*dr.z*dr.x;

          StrainDerivativeTensor[CurrentSystem].ay+=force_factor*dr.x*dr.y;
          StrainDerivativeTensor[CurrentSystem].by+=force_factor*dr.y*dr.y;
          StrainDerivativeTensor[CurrentSystem].cy+=force_factor*dr.z*dr.y;

          StrainDerivativeTensor[CurrentSystem].az+=force_factor*dr.x*dr.z;
          StrainDerivativeTensor[CurrentSystem].bz+=force_factor*dr.y*dr.z;
          StrainDerivativeTensor[CurrentSystem].cz+=force_factor*dr.z*dr.z;

          // forces
          f.x=force_factor*dr.x;
          f.y=force_factor*dr.y;
          f.z=force_factor*dr.z;

          Framework[CurrentSystem].Atoms[f1][A].Force.x-=f.x;
          Framework[CurrentSystem].Atoms[f1][A].Force.y-=f.y;
          Framework[CurrentSystem].Atoms[f1][A].Force.z-=f.z;

          Framework[CurrentSystem].Atoms[f1][B].Force.x+=f.x;
          Framework[CurrentSystem].Atoms[f1][B].Force.y+=f.y;
          Framework[CurrentSystem].Atoms[f1][B].Force.z+=f.z;
        }
      }
    }
  }

  return 0;
}

int CalculateFrameworkIntraChargeChargeForce(void)
{
  int i,j,typeA,typeB;
  REAL chargeA,chargeB;
  REAL r,rr,energy,DF;
  VECTOR posA,posB,dr,f;
  int f1,f2,A,B;
  REAL *parms;
  REAL alpha;

  // Framework-Framework energy
  UHostHostChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  // contributions from intra-framework
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        chargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;

        for(j=i+1;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],1))
          {
            typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
            posB=Framework[CurrentSystem].Atoms[f1][j].Position;
            chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              PotentialGradientCoulombic(chargeA,chargeB,rr,&energy,&DF);

              UHostHostChargeChargeReal[CurrentSystem]+=energy;

              f.x=-DF*dr.x;
              f.y=-DF*dr.y;
              f.z=-DF*dr.z;

              Framework[CurrentSystem].Atoms[f1][i].Force.x+=f.x;
              Framework[CurrentSystem].Atoms[f1][i].Force.y+=f.y;
              Framework[CurrentSystem].Atoms[f1][i].Force.z+=f.z;

              Framework[CurrentSystem].Atoms[f1][j].Force.x-=f.x;
              Framework[CurrentSystem].Atoms[f1][j].Force.y-=f.y;
              Framework[CurrentSystem].Atoms[f1][j].Force.z-=f.z;

              StrainDerivativeTensor[CurrentSystem].ax-=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx-=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx-=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay-=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by-=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy-=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az-=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz-=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz-=f.z*dr.z;
            }
          }
        }
      }
    }
  }

  // contributions from interactions between the frameworks
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1+1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        chargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;

        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
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
            PotentialGradientCoulombic(chargeA,chargeB,rr,&energy,&DF);

            UHostHostChargeChargeReal[CurrentSystem]+=energy;

            f.x=-DF*dr.x;
            f.y=-DF*dr.y;
            f.z=-DF*dr.z;

            Framework[CurrentSystem].Atoms[f1][i].Force.x+=f.x;
            Framework[CurrentSystem].Atoms[f1][i].Force.y+=f.y;
            Framework[CurrentSystem].Atoms[f1][i].Force.z+=f.z;

            Framework[CurrentSystem].Atoms[f2][j].Force.x-=f.x;
            Framework[CurrentSystem].Atoms[f2][j].Force.y-=f.y;
            Framework[CurrentSystem].Atoms[f2][j].Force.z-=f.z;

            StrainDerivativeTensor[CurrentSystem].ax-=f.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx-=f.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx-=f.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay-=f.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by-=f.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy-=f.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az-=f.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz-=f.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz-=f.z*dr.z;
          }
        }
      }
    }
  }

  // contributions from interactions 1-4 torsions
  // TODO: fix to handle anisotropic sites
  alpha=Alpha[CurrentSystem];
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[f1];i++)
      {
        parms=Framework[CurrentSystem].TorsionArguments[f1][i];

        if(fabs(parms[7])>1e-8)
        {
          A=Framework[CurrentSystem].Torsions[f1][i].A;
          chargeA=Framework[CurrentSystem].Atoms[f1][A].Charge;
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;

          B=Framework[CurrentSystem].Torsions[f1][i].D;
          chargeB=Framework[CurrentSystem].Atoms[f1][B].Charge;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;


          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          r=sqrt(rr);

          // note: no cutoff used here
          switch(ChargeMethod)
          {
            case NONE:
              DF=0.0;
              break;
            case SHIFTED_COULOMB:
            case TRUNCATED_COULOMB:
            case EWALD:
            default:
              UHostHostChargeChargeReal[CurrentSystem]+=parms[7]*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
              DF=parms[7]*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
              break;
          }

          StrainDerivativeTensor[CurrentSystem].ax-=DF*dr.x*dr.x;
          StrainDerivativeTensor[CurrentSystem].bx-=DF*dr.y*dr.x;
          StrainDerivativeTensor[CurrentSystem].cx-=DF*dr.z*dr.x;

          StrainDerivativeTensor[CurrentSystem].ay-=DF*dr.x*dr.y;
          StrainDerivativeTensor[CurrentSystem].by-=DF*dr.y*dr.y;
          StrainDerivativeTensor[CurrentSystem].cy-=DF*dr.z*dr.y;

          StrainDerivativeTensor[CurrentSystem].az-=DF*dr.x*dr.z;
          StrainDerivativeTensor[CurrentSystem].bz-=DF*dr.y*dr.z;
          StrainDerivativeTensor[CurrentSystem].cz-=DF*dr.z*dr.z;

          Framework[CurrentSystem].Atoms[f1][A].Force.x+=DF*dr.x;
          Framework[CurrentSystem].Atoms[f1][A].Force.y+=DF*dr.y;
          Framework[CurrentSystem].Atoms[f1][A].Force.z+=DF*dr.z;

          Framework[CurrentSystem].Atoms[f1][B].Force.x-=DF*dr.x;
          Framework[CurrentSystem].Atoms[f1][B].Force.y-=DF*dr.y;
          Framework[CurrentSystem].Atoms[f1][B].Force.z-=DF*dr.z;
        }
      }
    }
  }

  return 0;
}

int CalculateFrameworkIntraChargeBondDipoleForce(void)
{
  int i,j,f1,f2;
  int A1,A2;
  int Type;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL rr,ri2,energy,temp,length,chargeB;
  VECTOR dipoleA,fb1,fa1,fa2,term;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;

  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  UHostHostChargeBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  // contributions from intra-framework
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
      {
        DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
        A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
        A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
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

        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          // if framework are different, or if they are the same but not excluded within the framework
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][j][i],2))
          {
            Type=Framework[CurrentSystem].Atoms[f1][j].Type;
            if(PseudoAtoms[Type].HasCharges)
            {
              posB=Framework[CurrentSystem].Atoms[f1][j].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeBondDipoleSquared)
              {
                PotentialGradientChargeBondDipole(chargeB,DipoleMagnitudeA,length,dipoleA,dr,rr,&energy,&fa1,&fa2,&fb1,&term);

                UHostHostChargeBondDipoleReal[CurrentSystem]-=energy;

                Framework[CurrentSystem].Atoms[f1][j].Force.x+=fb1.x;
                Framework[CurrentSystem].Atoms[f1][j].Force.y+=fb1.y;
                Framework[CurrentSystem].Atoms[f1][j].Force.z+=fb1.z;

                Framework[CurrentSystem].Atoms[f1][A1].Force.x+=fa1.x;
                Framework[CurrentSystem].Atoms[f1][A1].Force.y+=fa1.y;
                Framework[CurrentSystem].Atoms[f1][A1].Force.z+=fa1.z;

                Framework[CurrentSystem].Atoms[f1][A2].Force.x+=fa2.x;
                Framework[CurrentSystem].Atoms[f1][A2].Force.y+=fa2.y;
                Framework[CurrentSystem].Atoms[f1][A2].Force.z+=fa2.z;

                // convert forces on atoms to molecular virial
                v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x);
                v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x);
                v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x);

                v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y);
                v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y);
                v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y);

                v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z);
                v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z);
                v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z);

                // the strain derivative
                StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
                StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
                StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

                StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
                StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
                StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

                StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
                StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
                StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
              }
            }
          }
        }
      }
    }
  }

  // contributions from interactions between the frameworks
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=0;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      if(f1!=f2)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
        {
          DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
          A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
          A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
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
            Type=Framework[CurrentSystem].Atoms[f2][j].Type;
            if(PseudoAtoms[Type].HasCharges)
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
                PotentialGradientChargeBondDipole(chargeB,DipoleMagnitudeA,length,dipoleA,dr,rr,&energy,&fa1,&fa2,&fb1,&term);

                UHostHostChargeBondDipoleReal[CurrentSystem]-=energy;

                Framework[CurrentSystem].Atoms[f1][j].Force.x+=fb1.x;
                Framework[CurrentSystem].Atoms[f1][j].Force.y+=fb1.y;
                Framework[CurrentSystem].Atoms[f1][j].Force.z+=fb1.z;

                Framework[CurrentSystem].Atoms[f2][A1].Force.x+=fa1.x;
                Framework[CurrentSystem].Atoms[f2][A1].Force.y+=fa1.y;
                Framework[CurrentSystem].Atoms[f2][A1].Force.z+=fa1.z;

                Framework[CurrentSystem].Atoms[f2][A2].Force.x+=fa2.x;
                Framework[CurrentSystem].Atoms[f2][A2].Force.y+=fa2.y;
                Framework[CurrentSystem].Atoms[f2][A2].Force.z+=fa2.z;

                // convert forces on atoms to molecular virial
                v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x);
                v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x);
                v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x);

                v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y);
                v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y);
                v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y);

                v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z);
                v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z);
                v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z);

                // the strain derivative
                StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
                StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
                StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

                StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
                StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
                StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

                StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
                StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
                StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraBondDipoleBondDipoleForce(void)
{
  int i,j,f1,f2;
  int A1,A2,B1,B2;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL ri2,rk2,energy,temp,length;
  VECTOR dipoleA,dipoleB,fb1,fb2,fa1,fa2,term;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL_MATRIX3x3 v;
  REAL r2;

  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;


  UHostHostBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  // contributions from intra-framework
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
      {
        DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
        A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
        A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
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

        for(j=i+1;j<Framework[CurrentSystem].NumberOfBondDipoles[f1];j++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][j];
          B1=Framework[CurrentSystem].BondDipoles[f1][j].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][j].B;

          // if framework are different, or if they are the same but not excluded within the framework
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],3))
          {
            posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
            posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
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
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(r2<CutOffBondDipoleBondDipoleSquared)
            {
              PotentialGradientBondDipoleBondDipole(DipoleMagnitudeA,ri2,dipoleA,DipoleMagnitudeB,rk2,dipoleB,dr,r2,
                                                    &energy,&fa1,&fa2,&fb1,&fb2,&term);

              UHostHostBondDipoleBondDipoleReal[CurrentSystem]+=energy;

              Framework[CurrentSystem].Atoms[f1][A1].Force.x+=fa1.x;
              Framework[CurrentSystem].Atoms[f1][A1].Force.y+=fa1.y;
              Framework[CurrentSystem].Atoms[f1][A1].Force.z+=fa1.z;

              Framework[CurrentSystem].Atoms[f1][A2].Force.x+=fa2.x;
              Framework[CurrentSystem].Atoms[f1][A2].Force.y+=fa2.y;
              Framework[CurrentSystem].Atoms[f1][A2].Force.z+=fa2.z;

              Framework[CurrentSystem].Atoms[f1][B1].Force.x+=fb1.x;
              Framework[CurrentSystem].Atoms[f1][B1].Force.y+=fb1.y;
              Framework[CurrentSystem].Atoms[f1][B1].Force.z+=fb1.z;

              Framework[CurrentSystem].Atoms[f1][B2].Force.x+=fb2.x;
              Framework[CurrentSystem].Atoms[f1][B2].Force.y+=fb2.y;
              Framework[CurrentSystem].Atoms[f1][B2].Force.z+=fb2.z;

              // convert forces on atoms to molecular virial
              // usually this conversion produces a torque on the center of mass and an asymmetric stress
              // by making it symmetric we regain the 'atomic' strain derivative
              v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
              v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
              v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

              v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
              v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
              v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

              v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
              v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
              v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

              // the strain derivative
              StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
              StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
              StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

              StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
              StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
              StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

              StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
              StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
              StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
            }
          }
        }
      }
    }
  }


  // contributions from interactions between the frameworks
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1+1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
      {
        DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
        A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
        A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
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

        for(j=0;j<Framework[CurrentSystem].NumberOfBondDipoles[f2];j++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f2][j];
          B1=Framework[CurrentSystem].BondDipoles[f2][j].A;
          B2=Framework[CurrentSystem].BondDipoles[f2][j].B;

          posB1=Framework[CurrentSystem].Atoms[f2][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f2][B2].Position;
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
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(r2<CutOffBondDipoleBondDipoleSquared)
          {
            PotentialGradientBondDipoleBondDipole(DipoleMagnitudeA,ri2,dipoleA,DipoleMagnitudeB,rk2,dipoleB,dr,r2,
                                                    &energy,&fa1,&fa2,&fb1,&fb2,&term);

            UHostHostBondDipoleBondDipoleReal[CurrentSystem]+=energy;

            Framework[CurrentSystem].Atoms[f1][A1].Force.x+=fa1.x;
            Framework[CurrentSystem].Atoms[f1][A1].Force.y+=fa1.y;
            Framework[CurrentSystem].Atoms[f1][A1].Force.z+=fa1.z;

            Framework[CurrentSystem].Atoms[f1][A2].Force.x+=fa2.x;
            Framework[CurrentSystem].Atoms[f1][A2].Force.y+=fa2.y;
            Framework[CurrentSystem].Atoms[f1][A2].Force.z+=fa2.z;

            Framework[CurrentSystem].Atoms[f2][B1].Force.x+=fb1.x;
            Framework[CurrentSystem].Atoms[f2][B1].Force.y+=fb1.y;
            Framework[CurrentSystem].Atoms[f2][B1].Force.z+=fb1.z;

            Framework[CurrentSystem].Atoms[f2][B2].Force.x+=fb2.x;
            Framework[CurrentSystem].Atoms[f2][B2].Force.y+=fb2.y;
            Framework[CurrentSystem].Atoms[f2][B2].Force.z+=fb2.z;

            // convert forces on atoms to molecular virial
            // usually this conversion produces a torque on the center of mass and an asymmetric stress
            // by making it symmetric we regain the 'atomic' strain derivative
            v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
            v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
            v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

            v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
            v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
            v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

            v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
            v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
            v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

            // the strain derivative
            StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
            StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
            StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

            StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
            StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
            StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

            StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
            StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
            StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
          }
        }
      }
    }
  }
  return 0;
}

void CalculateFrameworkAdsorbateVDWForce(void)
{
  int i,j,k,f1;
  int typeA,typeB,TypeMolA;
  REAL rr,scalingA;
  REAL energy,force_factor;
  VECTOR posA,posB,dr,f;
  REAL ReductionA,ReductionB;
  VECTOR drA,drB,fa,fb;
  VECTOR pos,posA1,posA2,posB1,posB2;
  REAL ra,rb;
  int ConnectedAtomA1,ConnectedAtomA2;
  int ConnectedAtomB1,ConnectedAtomB2;
  VECTOR v,w;
  REAL length_v,length_w;
  REAL dot_product;

  UHostAdsorbateVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeMolA=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
      scalingA=Adsorbates[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      ReductionA=0.0;
      ConnectedAtomA1=ConnectedAtomA2=-1;
      ra=drA.x=drA.y=drA.z=0.0;
      v.x=v.y=v.z=length_v=0.0;
      if(PseudoAtoms[typeA].AnisotropicCorrection)
      {
        switch(Components[TypeMolA].Connectivity[j])
        {
          case 0:
            break;
          case 1:
            ConnectedAtomA1=Components[TypeMolA].ConnectivityList[j][0];
            if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
              ReductionA=1.0+PseudoAtoms[typeA].AnisotropicDisplacement;
            else
            {
              posA1=Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Position;
              pos=Adsorbates[CurrentSystem][i].Atoms[j].Position;
              drA.x=posA1.x-pos.x;
              drA.y=posA1.y-pos.y;
              drA.z=posA1.z-pos.z;
              ra=sqrt(SQR(drA.x)+SQR(drA.y)+SQR(drA.z));
              ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
            }
            break;
          case 2:
            switch(Components[TypeMolA].AnisotropicType)
            {
              case ANISOTROPIC_BISECTION:
                fprintf(stderr, "ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateTotalInterVDWForce' (inter_force.c)\n");
                exit(0);
                break;
              case ANISOTROPIC_MID_POINT:
                ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
                ConnectedAtomA1=Components[TypeMolA].ConnectivityList[j][0];
                ConnectedAtomA2=Components[TypeMolA].ConnectivityList[j][1];
                pos=Adsorbates[CurrentSystem][i].Atoms[j].Position;
                posA1=Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Position;
                posA2=Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Position;
                v.x=pos.x-0.5*(posA1.x+posA2.x);
                v.y=pos.y-0.5*(posA1.y+posA2.y);
                v.z=pos.z-0.5*(posA1.z+posA2.z);
                length_v=sqrt(SQR(v.x)+SQR(v.y)+SQR(v.z));
                break;
              default:
                fprintf(stderr, "ERROR!\n");
                exit(0);
                break;
            }
            break;
          default:
            break;
        }
      }

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA])&&(!IsFractionalAdsorbateMolecule(i)))
      {
        UHostAdsorbateVDW[CurrentSystem]+=InterpolateVDWForceGrid(typeA,posA,&f);
        Adsorbates[CurrentSystem][i].Atoms[j].Force.x+=f.x;
        Adsorbates[CurrentSystem][i].Atoms[j].Force.y+=f.y;
        Adsorbates[CurrentSystem][i].Atoms[j].Force.z+=f.z;
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

            ReductionB=0.0;
            ConnectedAtomB1=ConnectedAtomB2=-1;
            rb=drB.x=drB.y=drB.z=0.0;
            w.x=w.y=w.z=length_w=0.0;
            if(PseudoAtoms[typeB].AnisotropicCorrection)
            {
              switch(Framework[CurrentSystem].Connectivity[f1][k])
              {
                case 0:
                  break;
                case 1:
                  ConnectedAtomB1=Framework[CurrentSystem].Neighbours[f1][k][0];
                  if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                    ReductionB=1.0+PseudoAtoms[typeB].AnisotropicDisplacement;
                  else
                  {
                    posB1=Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Position;
                    pos=Framework[CurrentSystem].Atoms[f1][k].Position;
                    drB.x=posB1.x-pos.x;
                    drB.y=posB1.y-pos.y;
                    drB.z=posB1.z-pos.z;
                    drB=ApplyBoundaryCondition(drB);
                    rb=sqrt(SQR(drB.x)+SQR(drB.y)+SQR(drB.z));
                    ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                  }
                  break;
                case 2:
                  switch(Framework[CurrentSystem].AnisotropicType)
                  {
                    case ANISOTROPIC_BISECTION:
                      fprintf(stderr, "ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                      exit(0);
                      break;
                    case ANISOTROPIC_MID_POINT:
                      ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                      ConnectedAtomB1=Framework[CurrentSystem].Neighbours[f1][k][0];
                      ConnectedAtomB2=Framework[CurrentSystem].Neighbours[f1][k][1];
                      pos=Framework[CurrentSystem].Atoms[f1][k].Position;
                      posB1=Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Position;
                      posB2=Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Position;
                      w.x=pos.x-0.5*(posB1.x+posB2.x);
                      w.y=pos.y-0.5*(posB1.y+posB2.y);
                      w.z=pos.z-0.5*(posB1.z+posB2.z);
                      w=ApplyBoundaryCondition(w);
                      length_w=sqrt(SQR(w.x)+SQR(w.y)+SQR(w.z));
                      break;
                    default:
                      fprintf(stderr, "ERROR: unknown Anisotropy type in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                      exit(0);
                      break;
                  }
                  break;
                default:
                  break;
              }
            }

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

              // compute the energy and force-factor
              PotentialGradient(typeA,typeB,rr,&energy,&force_factor,scalingA);

              // energy
              UHostAdsorbateVDW[CurrentSystem]+=energy;

              // forces
              f.x=force_factor*dr.x;
              f.y=force_factor*dr.y;
              f.z=force_factor*dr.z;

              if(PseudoAtoms[typeA].AnisotropicCorrection)
              {
                switch(Components[TypeMolA].Connectivity[j])
                {
                  case 0:
                    break;
                  case 1:
                    if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
                    {
                      Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=ReductionA*f.x;
                      Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=ReductionA*f.y;
                      Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=ReductionA*f.z;

                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=(1.0-ReductionA)*f.x;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=(1.0-ReductionA)*f.y;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=(1.0-ReductionA)*f.z;
                    }
                    else
                    {
                      dot_product=dr.x*drA.x+dr.y*drA.y+dr.z*drA.z;

                      fa.x=(-ReductionA*drA.x*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.x)*force_factor;
                      fa.y=(-ReductionA*drA.y*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.y)*force_factor;
                      fa.z=(-ReductionA*drA.z*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.z)*force_factor;

                      fb.x=(ReductionA*drA.x*dot_product/CUBE(ra)-(ReductionA/ra)*dr.x)*force_factor;
                      fb.y=(ReductionA*drA.y*dot_product/CUBE(ra)-(ReductionA/ra)*dr.y)*force_factor;
                      fb.z=(ReductionA*drA.z*dot_product/CUBE(ra)-(ReductionA/ra)*dr.z)*force_factor;

                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=fb.x;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=fb.y;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=fb.z;

                      Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=fa.x;
                      Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=fa.y;
                      Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=fa.z;
                    }
                    break;
                  case 2:
                    dot_product=v.x*dr.x+v.y*dr.y+v.z*dr.z;

                    fa.x=0.5*(ReductionA*v.x*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.x)*force_factor;
                    fa.y=0.5*(ReductionA*v.y*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.y)*force_factor;
                    fa.z=0.5*(ReductionA*v.z*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.z)*force_factor;

                    fb.x=(-ReductionA*v.x*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.x)*force_factor;
                    fb.y=(-ReductionA*v.y*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.y)*force_factor;
                    fb.z=(-ReductionA*v.z*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.z)*force_factor;

                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=fa.x;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=fa.y;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=fa.z;

                    Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=fb.x;
                    Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=fb.y;
                    Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=fb.z;

                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.x-=fa.x;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.y-=fa.y;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.z-=fa.z;
                    break;
                  default:
                    fprintf(stderr, "Not yet implemented in routine 'CalculateTotalInterVDWForce'\n");
                    exit(0);
                    break;
                }
              }
              else
              {
                Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=f.x;
                Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=f.y;
                Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=f.z;
              }

              // only if the framework is flexible do we need to compute the forces on the framework atoms
              if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
              {
                if(PseudoAtoms[typeB].AnisotropicCorrection)
                {
                  switch(Framework[CurrentSystem].Connectivity[f1][k])
                  {
                    case 0:
                      break;
                    case 1:
                      if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                      {
                        Framework[CurrentSystem].Atoms[f1][k].Force.x+=ReductionB*f.x;
                        Framework[CurrentSystem].Atoms[f1][k].Force.y+=ReductionB*f.y;
                        Framework[CurrentSystem].Atoms[f1][k].Force.z+=ReductionB*f.z;

                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.x+=(1.0-ReductionB)*f.x;
                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.y+=(1.0-ReductionB)*f.y;
                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.z+=(1.0-ReductionB)*f.z;
                      }
                      else
                      {
                        dot_product=dr.x*drB.x+dr.y*drB.y+dr.z*drB.z;

                        fa.x=(-(ReductionB*drA.x*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.x))*force_factor;
                        fa.y=(-(ReductionB*drA.y*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.y))*force_factor;
                        fa.z=(-(ReductionB*drA.z*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.z))*force_factor;

                        fb.x=(ReductionB*drB.x*dot_product/CUBE(rb)-(ReductionB/rb)*dr.x)*force_factor;
                        fb.y=(ReductionB*drB.y*dot_product/CUBE(rb)-(ReductionB/rb)*dr.y)*force_factor;
                        fb.z=(ReductionB*drB.z*dot_product/CUBE(rb)-(ReductionB/rb)*dr.z)*force_factor;

                        Framework[CurrentSystem].Atoms[f1][k].Force.x+=fa.x;
                        Framework[CurrentSystem].Atoms[f1][k].Force.y+=fa.y;
                        Framework[CurrentSystem].Atoms[f1][k].Force.z+=fa.z;

                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.x+=fb.x;
                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.y+=fb.y;
                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.z+=fb.z;
                      }
                      break;
                    case 2:
                      dot_product=w.x*dr.x+w.y*dr.y+w.z*dr.z;

                      fa.x=0.5*(ReductionB*w.x*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.x)*force_factor;
                      fa.y=0.5*(ReductionB*w.y*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.y)*force_factor;
                      fa.z=0.5*(ReductionB*w.z*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.z)*force_factor;

                      fb.x=(-ReductionB*w.x*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.x)*force_factor;
                      fb.y=(-ReductionB*w.y*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.y)*force_factor;
                      fb.z=(-ReductionB*w.z*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.z)*force_factor;

                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.x+=fa.x;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.y+=fa.y;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.z+=fa.z;

                      Framework[CurrentSystem].Atoms[f1][k].Force.x+=fb.x;
                      Framework[CurrentSystem].Atoms[f1][k].Force.y+=fb.y;
                      Framework[CurrentSystem].Atoms[f1][k].Force.z+=fb.z;

                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Force.x+=fa.x;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Force.y+=fa.y;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Force.z+=fa.z;
                      break;
                    default:
                      fprintf(stderr, "ERROR: Not yet implemented in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                      exit(0);
                      break;
                  }
                }
                else
                {
                  Framework[CurrentSystem].Atoms[f1][k].Force.x+=f.x;
                  Framework[CurrentSystem].Atoms[f1][k].Force.y+=f.y;
                  Framework[CurrentSystem].Atoms[f1][k].Force.z+=f.z;
                }
              }

              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
            }
          }
        }
      }
    }
  }
}

void CalculateFrameworkCationVDWForce(void)
{
  int i,j,k,f1;
  int typeA,typeB,TypeMolA;
  REAL rr,scalingA;
  REAL energy,force_factor;
  VECTOR posA,posB,dr,f;
  REAL ReductionA,ReductionB;
  VECTOR drA,drB,fa,fb;
  VECTOR pos,posA1,posA2,posB1,posB2;
  REAL ra,rb;
  int ConnectedAtomA1,ConnectedAtomA2;
  int ConnectedAtomB1,ConnectedAtomB2;
  VECTOR v,w;
  REAL length_v,length_w;
  REAL dot_product;

  UHostCationVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    TypeMolA=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
      scalingA=Cations[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      ReductionA=0.0;
      ConnectedAtomA1=ConnectedAtomA2=-1;
      ra=drA.x=drA.y=drA.z=0.0;
      v.x=v.y=v.z=length_v=0.0;
      if(PseudoAtoms[typeA].AnisotropicCorrection)
      {
        switch(Components[TypeMolA].Connectivity[j])
        {
          case 0:
            break;
          case 1:
            ConnectedAtomA1=Components[TypeMolA].ConnectivityList[j][0];
            if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
              ReductionA=1.0+PseudoAtoms[typeA].AnisotropicDisplacement;
            else
            {
              posA1=Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Position;
              pos=Cations[CurrentSystem][i].Atoms[j].Position;
              drA.x=posA1.x-pos.x;
              drA.y=posA1.y-pos.y;
              drA.z=posA1.z-pos.z;
              ra=sqrt(SQR(drA.x)+SQR(drA.y)+SQR(drA.z));
              ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
            }
            break;
          case 2:
            switch(Components[TypeMolA].AnisotropicType)
            {
              case ANISOTROPIC_BISECTION:
                fprintf(stderr, "ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateTotalInterVDWForce' (inter_force.c)\n");
                exit(0);
                break;
              case ANISOTROPIC_MID_POINT:
                ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
                ConnectedAtomA1=Components[TypeMolA].ConnectivityList[j][0];
                ConnectedAtomA2=Components[TypeMolA].ConnectivityList[j][1];
                pos=Cations[CurrentSystem][i].Atoms[j].Position;
                posA1=Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Position;
                posA2=Cations[CurrentSystem][i].Atoms[ConnectedAtomA2].Position;
                v.x=pos.x-0.5*(posA1.x+posA2.x);
                v.y=pos.y-0.5*(posA1.y+posA2.y);
                v.z=pos.z-0.5*(posA1.z+posA2.z);
                length_v=sqrt(SQR(v.x)+SQR(v.y)+SQR(v.z));
                break;
              default:
                fprintf(stderr, "ERROR!\n");
                exit(0);
                break;
            }
            break;
          default:
            break;
        }
      }

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA])&&(!IsFractionalCationMolecule(i)))
      {
        UHostCationVDW[CurrentSystem]+=InterpolateVDWForceGrid(typeA,posA,&f);
        Cations[CurrentSystem][i].Atoms[j].Force.x+=f.x;
        Cations[CurrentSystem][i].Atoms[j].Force.y+=f.y;
        Cations[CurrentSystem][i].Atoms[j].Force.z+=f.z;
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

            ReductionB=0.0;
            ConnectedAtomB1=ConnectedAtomB2=-1;
            rb=drB.x=drB.y=drB.z=0.0;
            w.x=w.y=w.z=length_w=0.0;
            if(PseudoAtoms[typeB].AnisotropicCorrection)
            {
              switch(Framework[CurrentSystem].Connectivity[f1][k])
              {
                case 0:
                  break;
                case 1:
                  ConnectedAtomB1=Framework[CurrentSystem].Neighbours[f1][k][0];
                  if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                    ReductionB=1.0+PseudoAtoms[typeB].AnisotropicDisplacement;
                  else
                  {
                    posB1=Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Position;
                    pos=Framework[CurrentSystem].Atoms[f1][k].Position;
                    drB.x=posB1.x-pos.x;
                    drB.y=posB1.y-pos.y;
                    drB.z=posB1.z-pos.z;
                    drB=ApplyBoundaryCondition(drB);
                    rb=sqrt(SQR(drB.x)+SQR(drB.y)+SQR(drB.z));
                    ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                  }
                  break;
                case 2:
                  switch(Framework[CurrentSystem].AnisotropicType)
                  {
                    case ANISOTROPIC_BISECTION:
                      fprintf(stderr, "ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                      exit(0);
                      break;
                    case ANISOTROPIC_MID_POINT:
                      ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                      ConnectedAtomB1=Framework[CurrentSystem].Neighbours[f1][k][0];
                      ConnectedAtomB2=Framework[CurrentSystem].Neighbours[f1][k][1];
                      pos=Framework[CurrentSystem].Atoms[f1][k].Position;
                      posB1=Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Position;
                      posB2=Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Position;
                      w.x=pos.x-0.5*(posB1.x+posB2.x);
                      w.y=pos.y-0.5*(posB1.y+posB2.y);
                      w.z=pos.z-0.5*(posB1.z+posB2.z);
                      w=ApplyBoundaryCondition(w);
                      length_w=sqrt(SQR(w.x)+SQR(w.y)+SQR(w.z));
                      break;
                    default:
                      fprintf(stderr, "ERROR: unknown Anisotropy type in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                      exit(0);
                      break;
                  }
                  break;
                default:
                  break;
              }
            }

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

              // compute the energy and force-factor
              PotentialGradient(typeA,typeB,rr,&energy,&force_factor,scalingA);

              // energy
              UHostCationVDW[CurrentSystem]+=energy;

              // forces
              f.x=force_factor*dr.x;
              f.y=force_factor*dr.y;
              f.z=force_factor*dr.z;

              if(PseudoAtoms[typeA].AnisotropicCorrection)
              {
                switch(Components[TypeMolA].Connectivity[j])
                {
                  case 0:
                    break;
                  case 1:
                    if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
                    {
                      Cations[CurrentSystem][i].Atoms[j].Force.x-=ReductionA*f.x;
                      Cations[CurrentSystem][i].Atoms[j].Force.y-=ReductionA*f.y;
                      Cations[CurrentSystem][i].Atoms[j].Force.z-=ReductionA*f.z;

                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=(1.0-ReductionA)*f.x;
                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=(1.0-ReductionA)*f.y;
                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=(1.0-ReductionA)*f.z;
                    }
                    else
                    {
                      dot_product=dr.x*drA.x+dr.y*drA.y+dr.z*drA.z;

                      fa.x=(-ReductionA*drA.x*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.x)*force_factor;
                      fa.y=(-ReductionA*drA.y*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.y)*force_factor;
                      fa.z=(-ReductionA*drA.z*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.z)*force_factor;

                      fb.x=(ReductionA*drA.x*dot_product/CUBE(ra)-(ReductionA/ra)*dr.x)*force_factor;
                      fb.y=(ReductionA*drA.y*dot_product/CUBE(ra)-(ReductionA/ra)*dr.y)*force_factor;
                      fb.z=(ReductionA*drA.z*dot_product/CUBE(ra)-(ReductionA/ra)*dr.z)*force_factor;

                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=fb.x;
                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=fb.y;
                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=fb.z;

                      Cations[CurrentSystem][i].Atoms[j].Force.x-=fa.x;
                      Cations[CurrentSystem][i].Atoms[j].Force.y-=fa.y;
                      Cations[CurrentSystem][i].Atoms[j].Force.z-=fa.z;
                    }
                    break;
                  case 2:
                    dot_product=v.x*dr.x+v.y*dr.y+v.z*dr.z;

                    fa.x=0.5*(ReductionA*v.x*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.x)*force_factor;
                    fa.y=0.5*(ReductionA*v.y*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.y)*force_factor;
                    fa.z=0.5*(ReductionA*v.z*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.z)*force_factor;

                    fb.x=(-ReductionA*v.x*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.x)*force_factor;
                    fb.y=(-ReductionA*v.y*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.y)*force_factor;
                    fb.z=(-ReductionA*v.z*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.z)*force_factor;

                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=fa.x;
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=fa.y;
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=fa.z;

                    Cations[CurrentSystem][i].Atoms[j].Force.x-=fb.x;
                    Cations[CurrentSystem][i].Atoms[j].Force.y-=fb.y;
                    Cations[CurrentSystem][i].Atoms[j].Force.z-=fb.z;

                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.x-=fa.x;
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.y-=fa.y;
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.z-=fa.z;
                    break;
                  default:
                    fprintf(stderr, "Not yet implemented in routine 'CalculateTotalInterVDWForce'\n");
                    exit(0);
                    break;
                }
              }
              else
              {
                Cations[CurrentSystem][i].Atoms[j].Force.x-=f.x;
                Cations[CurrentSystem][i].Atoms[j].Force.y-=f.y;
                Cations[CurrentSystem][i].Atoms[j].Force.z-=f.z;
              }

              // only if the framework is flexible do we need to compute the forces on the framework atoms
              if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
              {
                if(PseudoAtoms[typeB].AnisotropicCorrection)
                {
                  switch(Framework[CurrentSystem].Connectivity[f1][k])
                  {
                    case 0:
                      break;
                    case 1:
                      if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                      {
                        Framework[CurrentSystem].Atoms[f1][k].Force.x+=ReductionB*f.x;
                        Framework[CurrentSystem].Atoms[f1][k].Force.y+=ReductionB*f.y;
                        Framework[CurrentSystem].Atoms[f1][k].Force.z+=ReductionB*f.z;

                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.x+=(1.0-ReductionB)*f.x;
                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.y+=(1.0-ReductionB)*f.y;
                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.z+=(1.0-ReductionB)*f.z;
                      }
                      else
                      {
                        dot_product=dr.x*drB.x+dr.y*drB.y+dr.z*drB.z;

                        fa.x=(-(ReductionB*drA.x*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.x))*force_factor;
                        fa.y=(-(ReductionB*drA.y*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.y))*force_factor;
                        fa.z=(-(ReductionB*drA.z*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.z))*force_factor;

                        fb.x=(ReductionB*drB.x*dot_product/CUBE(rb)-(ReductionB/rb)*dr.x)*force_factor;
                        fb.y=(ReductionB*drB.y*dot_product/CUBE(rb)-(ReductionB/rb)*dr.y)*force_factor;
                        fb.z=(ReductionB*drB.z*dot_product/CUBE(rb)-(ReductionB/rb)*dr.z)*force_factor;

                        Framework[CurrentSystem].Atoms[f1][k].Force.x+=fa.x;
                        Framework[CurrentSystem].Atoms[f1][k].Force.y+=fa.y;
                        Framework[CurrentSystem].Atoms[f1][k].Force.z+=fa.z;

                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.x+=fb.x;
                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.y+=fb.y;
                        Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.z+=fb.z;
                      }
                      break;
                    case 2:
                      dot_product=w.x*dr.x+w.y*dr.y+w.z*dr.z;

                      fa.x=0.5*(ReductionB*w.x*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.x)*force_factor;
                      fa.y=0.5*(ReductionB*w.y*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.y)*force_factor;
                      fa.z=0.5*(ReductionB*w.z*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.z)*force_factor;

                      fb.x=(-ReductionB*w.x*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.x)*force_factor;
                      fb.y=(-ReductionB*w.y*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.y)*force_factor;
                      fb.z=(-ReductionB*w.z*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.z)*force_factor;

                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.x+=fa.x;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.y+=fa.y;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB1].Force.z+=fa.z;

                      Framework[CurrentSystem].Atoms[f1][k].Force.x+=fb.x;
                      Framework[CurrentSystem].Atoms[f1][k].Force.y+=fb.y;
                      Framework[CurrentSystem].Atoms[f1][k].Force.z+=fb.z;

                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Force.x+=fa.x;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Force.y+=fa.y;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomB2].Force.z+=fa.z;
                      break;
                    default:
                      fprintf(stderr, "ERROR: Not yet implemented in routine 'CalculateFrameworkIntraVDWForce' (framework_force.c)\n");
                      exit(0);
                      break;
                  }
                }
                else
                {
                  Framework[CurrentSystem].Atoms[f1][k].Force.x+=f.x;
                  Framework[CurrentSystem].Atoms[f1][k].Force.y+=f.y;
                  Framework[CurrentSystem].Atoms[f1][k].Force.z+=f.z;
                }
              }

              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
            }
          }
        }
      }
    }
  }
}


int CalculateFrameworkAdsorbateChargeChargeForce(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  REAL rr;
  REAL energy,chargeA,chargeB,temp;
  VECTOR posA,posB,dr,f;
  REAL scalingA;
  int TypeMolA;

  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeMolA=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      scalingA=Adsorbates[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
      chargeA=scalingA*Adsorbates[CurrentSystem][i].Atoms[j].Charge;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid)&&(!IsFractionalAdsorbateMolecule(i)))
      {
        UHostAdsorbateChargeChargeReal[CurrentSystem]+=InterpolateCoulombForceGrid(typeA,posA,&f);
        Adsorbates[CurrentSystem][i].Atoms[j].Force.x+=f.x;
        Adsorbates[CurrentSystem][i].Atoms[j].Force.y+=f.y;
        Adsorbates[CurrentSystem][i].Atoms[j].Force.z+=f.z;
      }
      else
      {
        // energy
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
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              PotentialGradientCoulombic(chargeA,chargeB,rr,&energy,&temp);

              UHostAdsorbateChargeChargeReal[CurrentSystem]+=energy;

              // forces
              f.x=temp*dr.x;
              f.y=temp*dr.y;
              f.z=temp*dr.z;

              Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=f.x;
              Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=f.y;
              Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=f.z;

              if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
              {
                Framework[CurrentSystem].Atoms[f1][k].Force.x+=f.x;
                Framework[CurrentSystem].Atoms[f1][k].Force.y+=f.y;
                Framework[CurrentSystem].Atoms[f1][k].Force.z+=f.z;
              }

              // stress-tensor
              StrainDerivativeTensor[CurrentSystem].ax+=temp*dr.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=temp*dr.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=temp*dr.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=temp*dr.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=temp*dr.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=temp*dr.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=temp*dr.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=temp*dr.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=temp*dr.z*dr.z;
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationChargeChargeForce(void)
{
  int i,j,k,f1;
  int typeA,typeB;
  REAL energy,rr;
  REAL chargeA,chargeB,temp;
  VECTOR posA,posB,dr,f;
  REAL scalingA;
  int TypeMolA;

  UHostCationChargeChargeReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    TypeMolA=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      scalingA=Cations[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
      chargeA=scalingA*Cations[CurrentSystem][i].Atoms[j].Charge;
      posA=Cations[CurrentSystem][i].Atoms[j].Position;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid)&&(!IsFractionalCationMolecule(i)))
      {
        UHostCationChargeChargeReal[CurrentSystem]+=InterpolateCoulombForceGrid(typeA,posA,&f);
        Cations[CurrentSystem][i].Atoms[j].Force.x+=f.x;
        Cations[CurrentSystem][i].Atoms[j].Force.y+=f.y;
        Cations[CurrentSystem][i].Atoms[j].Force.z+=f.z;
      }
      else
      {
        // energy
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
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              PotentialGradientCoulombic(chargeA,chargeB,rr,&energy,&temp);

              UHostCationChargeChargeReal[CurrentSystem]+=energy;

              // forces
              f.x=temp*dr.x;
              f.y=temp*dr.y;
              f.z=temp*dr.z;

              Cations[CurrentSystem][i].Atoms[j].Force.x-=f.x;
              Cations[CurrentSystem][i].Atoms[j].Force.y-=f.y;
              Cations[CurrentSystem][i].Atoms[j].Force.z-=f.z;

              if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
              {
                Framework[CurrentSystem].Atoms[f1][k].Force.x+=f.x;
                Framework[CurrentSystem].Atoms[f1][k].Force.y+=f.y;
                Framework[CurrentSystem].Atoms[f1][k].Force.z+=f.z;
              }

              // stress-tensor
              StrainDerivativeTensor[CurrentSystem].ax+=temp*dr.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=temp*dr.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=temp*dr.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=temp*dr.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=temp*dr.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=temp*dr.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=temp*dr.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=temp*dr.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=temp*dr.z*dr.z;
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateChargeBondDipoleForce(void)
{
  int i,j,k,f1;
  int A1,A2;
  int Type,TypeA;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL rr,ri2,energy,temp,length,chargeB;
  VECTOR dipoleA,fb1,fa1,fa2,term;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;

  UHostAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

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

      for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
      {
        for(k=0;k<Adsorbates[CurrentSystem][j].NumberOfAtoms;k++)
        {
          Type=Adsorbates[CurrentSystem][j].Atoms[k].Type;
          if(PseudoAtoms[Type].HasCharges)
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
              PotentialGradientChargeBondDipole(chargeB,DipoleMagnitudeA,length,dipoleA,dr,rr,&energy,&fa1,&fa2,&fb1,&term);

              UHostAdsorbateChargeBondDipoleReal[CurrentSystem]-=energy;

              Adsorbates[CurrentSystem][j].Atoms[k].Force.x+=fb1.x;
              Adsorbates[CurrentSystem][j].Atoms[k].Force.y+=fb1.y;
              Adsorbates[CurrentSystem][j].Atoms[k].Force.z+=fb1.z;

              if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
              {
                Framework[CurrentSystem].Atoms[f1][A1].Force.x+=fa1.x;
                Framework[CurrentSystem].Atoms[f1][A1].Force.y+=fa1.y;
                Framework[CurrentSystem].Atoms[f1][A1].Force.z+=fa1.z;

                Framework[CurrentSystem].Atoms[f1][A2].Force.x+=fa2.x;
                Framework[CurrentSystem].Atoms[f1][A2].Force.y+=fa2.y;
                Framework[CurrentSystem].Atoms[f1][A2].Force.z+=fa2.z;

                // convert forces on atoms to molecular virial
                v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x);
                v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x);
                v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x);

                v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y);
                v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y);
                v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y);

                v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z);
                v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z);
                v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z);
              }


              // the strain derivative
              StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
              StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
              StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

              StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
              StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
              StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

              StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
              StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
              StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
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
          Type=Framework[CurrentSystem].Atoms[f1][k].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            PotentialGradientChargeBondDipole(chargeB,DipoleMagnitudeA,length,dipoleA,dr,rr,&energy,&fa1,&fa2,&fb1,&term);

            UHostAdsorbateChargeBondDipoleReal[CurrentSystem]-=energy;

            if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
            {
              Framework[CurrentSystem].Atoms[f1][k].Force.x+=fb1.x;
              Framework[CurrentSystem].Atoms[f1][k].Force.y+=fb1.y;
              Framework[CurrentSystem].Atoms[f1][k].Force.z+=fb1.z;
            }

            Adsorbates[CurrentSystem][i].Atoms[A1].Force.x+=fa1.x;
            Adsorbates[CurrentSystem][i].Atoms[A1].Force.y+=fa1.y;
            Adsorbates[CurrentSystem][i].Atoms[A1].Force.z+=fa1.z;

            Adsorbates[CurrentSystem][i].Atoms[A2].Force.x+=fa2.x;
            Adsorbates[CurrentSystem][i].Atoms[A2].Force.y+=fa2.y;
            Adsorbates[CurrentSystem][i].Atoms[A2].Force.z+=fa2.z;

            // convert forces on atoms to molecular virial
            v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x);
            v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x);
            v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x);

            v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y);
            v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y);
            v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y);

            v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z);
            v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z);
            v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z);

            // the strain derivative
            StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
            StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
            StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

            StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
            StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
            StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

            StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
            StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
            StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationChargeBondDipoleForce(void)
{
  int i,j,k,f1;
  int A1,A2;
  int Type,TypeA;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL rr,ri2,energy,temp,length,chargeB;
  VECTOR dipoleA,fb1,fa1,fa2,term;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;

  UHostCationChargeBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

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

      for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
      {
        for(k=0;k<Cations[CurrentSystem][j].NumberOfAtoms;k++)
        {
          Type=Cations[CurrentSystem][j].Atoms[k].Type;
          if(PseudoAtoms[Type].HasCharges)
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
              PotentialGradientChargeBondDipole(chargeB,DipoleMagnitudeA,length,dipoleA,dr,rr,&energy,&fa1,&fa2,&fb1,&term);

              Cations[CurrentSystem][j].Atoms[k].Force.x+=fb1.x;
              Cations[CurrentSystem][j].Atoms[k].Force.y+=fb1.y;
              Cations[CurrentSystem][j].Atoms[k].Force.z+=fb1.z;

              UHostCationChargeBondDipoleReal[CurrentSystem]-=energy;

              if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
              {
                Framework[CurrentSystem].Atoms[f1][A1].Force.x+=fa1.x;
                Framework[CurrentSystem].Atoms[f1][A1].Force.y+=fa1.y;
                Framework[CurrentSystem].Atoms[f1][A1].Force.z+=fa1.z;

                Framework[CurrentSystem].Atoms[f1][A2].Force.x+=fa2.x;
                Framework[CurrentSystem].Atoms[f1][A2].Force.y+=fa2.y;
                Framework[CurrentSystem].Atoms[f1][A2].Force.z+=fa2.z;
              }

              // convert forces on atoms to molecular virial
              v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x);
              v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x);
              v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x);

              v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y);
              v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y);
              v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y);

              v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z);
              v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z);
              v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z);

              // the strain derivative
              StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
              StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
              StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

              StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
              StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
              StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

              StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
              StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
              StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
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
          Type=Framework[CurrentSystem].Atoms[f1][k].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            PotentialGradientChargeBondDipole(chargeB,DipoleMagnitudeA,length,dipoleA,dr,rr,&energy,&fa1,&fa2,&fb1,&term);

            UHostCationChargeBondDipoleReal[CurrentSystem]-=energy;

            if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
            {
              Framework[CurrentSystem].Atoms[f1][k].Force.x+=fb1.x;
              Framework[CurrentSystem].Atoms[f1][k].Force.y+=fb1.y;
              Framework[CurrentSystem].Atoms[f1][k].Force.z+=fb1.z;
            }

            Cations[CurrentSystem][i].Atoms[A1].Force.x+=fa1.x;
            Cations[CurrentSystem][i].Atoms[A1].Force.y+=fa1.y;
            Cations[CurrentSystem][i].Atoms[A1].Force.z+=fa1.z;

            Cations[CurrentSystem][i].Atoms[A2].Force.x+=fa2.x;
            Cations[CurrentSystem][i].Atoms[A2].Force.y+=fa2.y;
            Cations[CurrentSystem][i].Atoms[A2].Force.z+=fa2.z;

            // convert forces on atoms to molecular virial
            v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x);
            v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x);
            v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x);

            v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y);
            v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y);
            v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y);

            v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z);
            v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z);
            v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z);

            // the strain derivative
            StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
            StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
            StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

            StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
            StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
            StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

            StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
            StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
            StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateBondDipoleBondDipoleForce(void)
{
  int i,k,l,f1;
  int A1,A2,B1,B2;
  int TypeB;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL rr,ri2,rk2,energy,temp,length;
  VECTOR dipoleA,dipoleB,fb1,fb2,fa1,fa2,term;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL_MATRIX3x3 v;

  UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

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
            PotentialGradientBondDipoleBondDipole(DipoleMagnitudeA,ri2,dipoleA,DipoleMagnitudeB,rk2,dipoleB,dr,rr,
                                                    &energy,&fa1,&fa2,&fb1,&fb2,&term);

            UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=energy;

            if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
            {
              Framework[CurrentSystem].Atoms[f1][A1].Force.x+=fa1.x;
              Framework[CurrentSystem].Atoms[f1][A1].Force.y+=fa1.y;
              Framework[CurrentSystem].Atoms[f1][A1].Force.z+=fa1.z;

              Framework[CurrentSystem].Atoms[f1][A2].Force.x+=fa2.x;
              Framework[CurrentSystem].Atoms[f1][A2].Force.y+=fa2.y;
              Framework[CurrentSystem].Atoms[f1][A2].Force.z+=fa2.z;
            }

            Adsorbates[CurrentSystem][k].Atoms[B1].Force.x+=fb1.x;
            Adsorbates[CurrentSystem][k].Atoms[B1].Force.y+=fb1.y;
            Adsorbates[CurrentSystem][k].Atoms[B1].Force.z+=fb1.z;

            Adsorbates[CurrentSystem][k].Atoms[B2].Force.x+=fb2.x;
            Adsorbates[CurrentSystem][k].Atoms[B2].Force.y+=fb2.y;
            Adsorbates[CurrentSystem][k].Atoms[B2].Force.z+=fb2.z;

            // convert forces on atoms to molecular virial
            // usually this conversion produces a torque on the center of mass and an asymmetric stress
            // by making it symmetric we regain the 'atomic' strain derivative
            v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
            v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
            v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

            v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
            v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
            v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

            v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
            v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
            v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

            // the strain derivative
            StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
            StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
            StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

            StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
            StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
            StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

            StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
            StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
            StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationBondDipoleBondDipoleForce(void)
{
  int i,k,l,f1;
  int A1,A2,B1,B2;
  int TypeB;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL rr,ri2,rk2,energy,temp,length;
  VECTOR dipoleA,dipoleB,fb1,fb2,fa1,fa2,term;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL_MATRIX3x3 v;

  UHostCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

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
            PotentialGradientBondDipoleBondDipole(DipoleMagnitudeA,ri2,dipoleA,DipoleMagnitudeB,rk2,dipoleB,dr,rr,
                                                    &energy,&fa1,&fa2,&fb1,&fb2,&term);

            UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=energy;

            if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
            {
              Framework[CurrentSystem].Atoms[f1][A1].Force.x+=fa1.x;
              Framework[CurrentSystem].Atoms[f1][A1].Force.y+=fa1.y;
              Framework[CurrentSystem].Atoms[f1][A1].Force.z+=fa1.z;

              Framework[CurrentSystem].Atoms[f1][A2].Force.x+=fa2.x;
              Framework[CurrentSystem].Atoms[f1][A2].Force.y+=fa2.y;
              Framework[CurrentSystem].Atoms[f1][A2].Force.z+=fa2.z;
            }

            Cations[CurrentSystem][k].Atoms[B1].Force.x+=fb1.x;
            Cations[CurrentSystem][k].Atoms[B1].Force.y+=fb1.y;
            Cations[CurrentSystem][k].Atoms[B1].Force.z+=fb1.z;

            Cations[CurrentSystem][k].Atoms[B2].Force.x+=fb2.x;
            Cations[CurrentSystem][k].Atoms[B2].Force.y+=fb2.y;
            Cations[CurrentSystem][k].Atoms[B2].Force.z+=fb2.z;

            // convert forces on atoms to molecular virial
            // usually this conversion produces a torque on the center of mass and an asymmetric stress
            // by making it symmetric we regain the 'atomic' strain derivative
            v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
            v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
            v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

            v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
            v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
            v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

            v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
            v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
            v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

            // the strain derivative
            StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
            StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
            StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

            StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
            StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
            StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

            StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
            StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
            StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
          }
        }
      }
    }
  }
  return 0;
}


int CalculateFrameworkChargeChargeElectricFieldMC(int New,int excl_ads,int excl_cation)
{
  int i,j,k,f1;
  int typeA,typeB;
  REAL rr;
  REAL chargeA,chargeB,temp;
  REAL  force_factor_A,force_factor_B;
  VECTOR posA,posB,dr;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    if(i!=excl_ads)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
        chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;

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

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              temp=PotentialValueCoulombic(1.0,1.0,sqrt(rr));

              // forces
              force_factor_A=temp*chargeB;
              force_factor_B=temp*chargeA;

              Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x-=force_factor_A*dr.x;
              Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y-=force_factor_A*dr.y;
              Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z-=force_factor_A*dr.z;

              if((BackPolarization)||(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE))
              {
                Framework[CurrentSystem].Atoms[f1][k].ElectricField.x+=force_factor_B*dr.x;
                Framework[CurrentSystem].Atoms[f1][k].ElectricField.y+=force_factor_B*dr.y;
                Framework[CurrentSystem].Atoms[f1][k].ElectricField.z+=force_factor_B*dr.z;
              }
            }
          }
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    if(i!=excl_cation)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        typeA=Cations[CurrentSystem][i].Atoms[j].Type;
        posA=Cations[CurrentSystem][i].Atoms[j].Position;
        chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;

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

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              temp=PotentialValueCoulombic(1.0,1.0,sqrt(rr));

              // forces
              force_factor_A=temp*chargeB;
              force_factor_B=temp*chargeA;

              Cations[CurrentSystem][i].Atoms[j].ElectricField.x-=force_factor_A*dr.x;
              Cations[CurrentSystem][i].Atoms[j].ElectricField.y-=force_factor_A*dr.y;
              Cations[CurrentSystem][i].Atoms[j].ElectricField.z-=force_factor_A*dr.z;

              if((BackPolarization)||(Framework[CurrentSystem].FrameworkModel==FLEXIBLE))
              {
                Framework[CurrentSystem].Atoms[f1][k].ElectricField.x+=force_factor_B*dr.x;
                Framework[CurrentSystem].Atoms[f1][k].ElectricField.y+=force_factor_B*dr.y;
                Framework[CurrentSystem].Atoms[f1][k].ElectricField.z+=force_factor_B*dr.z;
              }
            }
          }
        }
      }
    }
  }

  if(New)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      posA=TrialPosition[CurrentSystem][i];
      typeA=Components[CurrentComponent].Type[i];
      chargeA=PseudoAtoms[typeA].Charge1;

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

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            temp=PotentialValueCoulombic(1.0,1.0,sqrt(rr));

            // forces
            force_factor_A=temp*chargeB;
            force_factor_B=temp*chargeA;

            ElectricFieldAtTrialPosition[CurrentSystem][i].x-=force_factor_A*dr.x;
            ElectricFieldAtTrialPosition[CurrentSystem][i].y-=force_factor_A*dr.y;
            ElectricFieldAtTrialPosition[CurrentSystem][i].z-=force_factor_A*dr.z;

            if((BackPolarization)||(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE))
            {
              Framework[CurrentSystem].Atoms[f1][k].ElectricField.x+=force_factor_B*dr.x;
              Framework[CurrentSystem].Atoms[f1][k].ElectricField.y+=force_factor_B*dr.y;
              Framework[CurrentSystem].Atoms[f1][k].ElectricField.z+=force_factor_B*dr.z;
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateElectricFieldFromInducedDipole(void)
{
  int i,k,l,f1;
  int TypeB;
  VECTOR posA,posB,dr;
  REAL rr;
  VECTOR dipoleA,dipoleB;
  VECTOR termA,termB;

  if(ChargeMethod==NONE) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      posA=Framework[CurrentSystem].Atoms[f1][i].Position;
      dipoleA=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;

      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        TypeB=Adsorbates[CurrentSystem][k].Type;
        for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
        {
          posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
          dipoleB=Adsorbates[CurrentSystem][k].Atoms[l].InducedDipole;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            PotentialElectricFieldBondDipoleBondDipole(dipoleA,dipoleB,dr,rr,&termA,&termB);

            if((BackPolarization)||(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE))
            {
              Framework[CurrentSystem].Atoms[f1][i].ElectricField.x+=termA.x;
              Framework[CurrentSystem].Atoms[f1][i].ElectricField.y+=termA.y;
              Framework[CurrentSystem].Atoms[f1][i].ElectricField.z+=termA.z;
            }

            Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.x+=termB.x;
            Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.y+=termB.y;
            Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.z+=termB.z;
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationElectricFieldFromInducedDipole(void)
{
  int i,k,l,f1;
  int TypeB;
  VECTOR posA,posB,dr;
  REAL rr;
  VECTOR dipoleA,dipoleB;
  VECTOR termA,termB;

  if(ChargeMethod==NONE) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      posA=Framework[CurrentSystem].Atoms[f1][i].Position;
      dipoleA=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;

      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        TypeB=Cations[CurrentSystem][k].Type;
        for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
        {
          posB=Cations[CurrentSystem][k].Atoms[l].Position;
          dipoleB=Cations[CurrentSystem][k].Atoms[l].InducedDipole;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            PotentialElectricFieldBondDipoleBondDipole(dipoleA,dipoleB,dr,rr,&termA,&termB);

            if((BackPolarization)||(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE))
            {
              Framework[CurrentSystem].Atoms[f1][i].ElectricField.x+=termA.x;
              Framework[CurrentSystem].Atoms[f1][i].ElectricField.y+=termA.y;
              Framework[CurrentSystem].Atoms[f1][i].ElectricField.z+=termA.z;
            }

            Cations[CurrentSystem][k].Atoms[l].ElectricField.x+=termB.x;
            Cations[CurrentSystem][k].Atoms[l].ElectricField.y+=termB.y;
            Cations[CurrentSystem][k].Atoms[l].ElectricField.z+=termB.z;
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkMoleculeElectricFieldFromInducedDipoleMC(int New,int excl_ads,int excl_cation)
{
  int i,k,l,f1;
  int TypeB;
  VECTOR posA,posB,dr;
  REAL rr;
  VECTOR dipoleA,dipoleB;
  VECTOR termA,termB;

  if(ChargeMethod==NONE) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      posA=Framework[CurrentSystem].Atoms[f1][i].Position;
      dipoleA=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;

      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        if(k!=excl_ads)
        {
          TypeB=Adsorbates[CurrentSystem][k].Type;
          for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
          {
            posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
            dipoleB=Adsorbates[CurrentSystem][k].Atoms[l].InducedDipole;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              PotentialElectricFieldBondDipoleBondDipole(dipoleA,dipoleB,dr,rr,&termA,&termB);

              if((BackPolarization)||(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE))
              {
                Framework[CurrentSystem].Atoms[f1][i].ElectricField.x+=termB.x;
                Framework[CurrentSystem].Atoms[f1][i].ElectricField.y+=termB.y;
                Framework[CurrentSystem].Atoms[f1][i].ElectricField.z+=termB.z;
              }

              Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.x+=termA.x;
              Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.y+=termA.y;
              Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.z+=termA.z;
            }
          }
        }
      }

      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        if(k!=excl_cation)
        {
          TypeB=Cations[CurrentSystem][k].Type;
          for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
          {
            posB=Cations[CurrentSystem][k].Atoms[l].Position;
            dipoleB=Cations[CurrentSystem][k].Atoms[l].InducedDipole;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              PotentialElectricFieldBondDipoleBondDipole(dipoleA,dipoleB,dr,rr,&termA,&termB);

              if((BackPolarization)||(Framework[CurrentSystem].FrameworkModel==FLEXIBLE))
              {
                Framework[CurrentSystem].Atoms[f1][i].ElectricField.x+=termB.x;
                Framework[CurrentSystem].Atoms[f1][i].ElectricField.y+=termB.y;
                Framework[CurrentSystem].Atoms[f1][i].ElectricField.z+=termB.z;
              }

              Cations[CurrentSystem][k].Atoms[l].ElectricField.x+=termA.x;
              Cations[CurrentSystem][k].Atoms[l].ElectricField.y+=termA.y;
              Cations[CurrentSystem][k].Atoms[l].ElectricField.z+=termA.z;
            }
          }
        }
      }

    }
  }

  if(New)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      posA=TrialPosition[CurrentSystem][i];
      dipoleA=InducedDipoleAtTrialPosition[CurrentSystem][i];

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          dipoleB=Framework[CurrentSystem].Atoms[f1][k].InducedDipole;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            PotentialElectricFieldBondDipoleBondDipole(dipoleA,dipoleB,dr,rr,&termA,&termB);

            ElectricFieldAtTrialPosition[CurrentSystem][i].x+=termB.x;
            ElectricFieldAtTrialPosition[CurrentSystem][i].y+=termB.y;
            ElectricFieldAtTrialPosition[CurrentSystem][i].z+=termB.z;

            if((BackPolarization)||(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE))
            {
              Framework[CurrentSystem].Atoms[f1][k].ElectricField.x+=termA.x;
              Framework[CurrentSystem].Atoms[f1][k].ElectricField.y+=termA.y;
              Framework[CurrentSystem].Atoms[f1][k].ElectricField.z+=termA.z;
            }


          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateChargeInducedDipoleForce(void)
{
  int i,j,k,f1;
  int Type,TypeA;
  VECTOR posA,posB,dr;
  REAL r,rr,chargeB;
  VECTOR dipoleA,termA,termB;

  if(ChargeMethod==NONE) return 0;

  if(BackPolarization)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        dipoleA=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;

        for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
        {
          for(k=0;k<Adsorbates[CurrentSystem][j].NumberOfAtoms;k++)
          {
            Type=Adsorbates[CurrentSystem][j].Atoms[k].Type;
            if(PseudoAtoms[Type].HasCharges)
            {
              posB=Adsorbates[CurrentSystem][j].Atoms[k].Position;
              chargeB=Adsorbates[CurrentSystem][j].Atoms[k].Charge;

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared[CurrentSystem])
              {
                PotentialElectricFieldBondDipoleBondDipole(dipoleA,dipoleA,dr,rr,&termA,&termB);

                if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
                {
                  Framework[CurrentSystem].Atoms[f1][i].Force.x-=chargeB*termA.x;
                  Framework[CurrentSystem].Atoms[f1][i].Force.y-=chargeB*termA.y;
                  Framework[CurrentSystem].Atoms[f1][i].Force.z-=chargeB*termA.z;
                }

                Adsorbates[CurrentSystem][j].Atoms[k].Force.x+=chargeB*termA.x;
                Adsorbates[CurrentSystem][j].Atoms[k].Force.y+=chargeB*termA.y;
                Adsorbates[CurrentSystem][j].Atoms[k].Force.z+=chargeB*termA.z;
              }
            }
          }
        }
      }
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeA=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
      dipoleA=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          Type=Framework[CurrentSystem].Atoms[f1][k].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          r=sqrt(rr);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            PotentialElectricFieldBondDipoleBondDipole(dipoleA,dipoleA,dr,rr,&termA,&termB);

            Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=chargeB*termA.x;
            Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=chargeB*termA.y;
            Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=chargeB*termA.z;

            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              Framework[CurrentSystem].Atoms[f1][k].Force.x+=chargeB*termA.x;
              Framework[CurrentSystem].Atoms[f1][k].Force.y+=chargeB*termA.y;
              Framework[CurrentSystem].Atoms[f1][k].Force.z+=chargeB*termA.z;
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateInducedDipoleInducedDipoleForce(void)
{
  int i,k,l,f1;
  int TypeB;
  VECTOR posA,posB,dr;
  REAL energy,rr;
  VECTOR dipoleA,dipoleB,term;

  if(ChargeMethod==NONE) return 0;
  if(!BackPolarization) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      posA=Framework[CurrentSystem].Atoms[f1][i].Position;
      dipoleA=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;

      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        TypeB=Adsorbates[CurrentSystem][k].Type;
        for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
        {
          posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
          dipoleB=Adsorbates[CurrentSystem][k].Atoms[l].InducedDipole;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            PotentialGradientInducedDipoleInducedDipole(dipoleA,dipoleB,dr,rr,&energy,&term);

            if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
            {
              Framework[CurrentSystem].Atoms[f1][i].Force.x-=term.x;
              Framework[CurrentSystem].Atoms[f1][i].Force.y-=term.y;
              Framework[CurrentSystem].Atoms[f1][i].Force.z-=term.z;
            }

            Adsorbates[CurrentSystem][k].Atoms[l].Force.x+=term.x;
            Adsorbates[CurrentSystem][k].Atoms[l].Force.y+=term.y;
            Adsorbates[CurrentSystem][k].Atoms[l].Force.z+=term.z;

          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraReplicaVDWForce(void)
{
  int i,j,typeA,typeB,start;
  REAL energy,force_factor;
  REAL rr,ReductionA,ReductionB;
  VECTOR posA,posB,dr,f;
  int ConnectedAtomA,ConnectedAtomB;
  int f1,f2,k1,k2,k3,ncell,index_j;
  int A,B;
  REAL *parms;

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
        posA=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

        ReductionA=1.0;
        ConnectedAtomA=-1;

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
                  posB=Framework[CurrentSystem].Atoms[f2][j].AnisotropicPosition;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  ReductionB=1.0;
                  ConnectedAtomB=-1;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffVDWSquared)
                  {
                    PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);

                    if(ncell==0)
                      UHostHostVDW[CurrentSystem]+=energy;
                    else
                      UHostHostVDW[CurrentSystem]+=0.5*energy;

                    f.x=force_factor*dr.x;
                    f.y=force_factor*dr.y;
                    f.z=force_factor*dr.z;

                    Framework[CurrentSystem].Atoms[f1][i].Force.x-=ReductionA*f.x;
                    Framework[CurrentSystem].Atoms[f1][i].Force.y-=ReductionA*f.y;
                    Framework[CurrentSystem].Atoms[f1][i].Force.z-=ReductionA*f.z;

                    if(ConnectedAtomA>=0)
                    {
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomA].Force.x-=(1.0-ReductionA)*f.x;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomA].Force.y-=(1.0-ReductionA)*f.y;
                      Framework[CurrentSystem].Atoms[f1][ConnectedAtomA].Force.z-=(1.0-ReductionA)*f.z;
                    }

                    if(ncell==0)
                    {
                      Framework[CurrentSystem].Atoms[f2][j].Force.x+=ReductionB*f.x;
                      Framework[CurrentSystem].Atoms[f2][j].Force.y+=ReductionB*f.y;
                      Framework[CurrentSystem].Atoms[f2][j].Force.z+=ReductionB*f.z;

                      if(ConnectedAtomB>=0)
                      {
                        Framework[CurrentSystem].Atoms[f2][ConnectedAtomB].Force.x+=(1.0-ReductionB)*f.x;
                        Framework[CurrentSystem].Atoms[f2][ConnectedAtomB].Force.y+=(1.0-ReductionB)*f.y;
                        Framework[CurrentSystem].Atoms[f2][ConnectedAtomB].Force.z+=(1.0-ReductionB)*f.z;
                      }
                      StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
                    }
                    else
                    {
                      // no force on particle outside the main unit cell
                      // only the interaction on 'i' in the unit cell counts
                      StrainDerivativeTensor[CurrentSystem].ax+=0.5*f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx+=0.5*f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx+=0.5*f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay+=0.5*f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by+=0.5*f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy+=0.5*f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az+=0.5*f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz+=0.5*f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz+=0.5*f.z*dr.z;
                    }
                  }
                }
              }
              ncell++;
            }
      }
    }
  }

  // contributions from interactions 1-4 torsions
  // TODO: fix to handle anisotropic sites
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[f1];i++)
      {
        parms=Framework[CurrentSystem].TorsionArguments[f1][i];
        if(fabs(parms[6])>1e-8)
        {
          A=Framework[CurrentSystem].Torsions[f1][i].A;
          typeA=Framework[CurrentSystem].Atoms[f1][A].Type;
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;

          B=Framework[CurrentSystem].Torsions[f1][i].D;
          typeB=Framework[CurrentSystem].Atoms[f1][B].Type;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;


          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);

          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);

          UHostHostVDW[CurrentSystem]+=parms[6]*energy;

          force_factor=parms[6]*force_factor;
          StrainDerivativeTensor[CurrentSystem].ax+=force_factor*dr.x*dr.x;
          StrainDerivativeTensor[CurrentSystem].bx+=force_factor*dr.y*dr.x;
          StrainDerivativeTensor[CurrentSystem].cx+=force_factor*dr.z*dr.x;

          StrainDerivativeTensor[CurrentSystem].ay+=force_factor*dr.x*dr.y;
          StrainDerivativeTensor[CurrentSystem].by+=force_factor*dr.y*dr.y;
          StrainDerivativeTensor[CurrentSystem].cy+=force_factor*dr.z*dr.y;

          StrainDerivativeTensor[CurrentSystem].az+=force_factor*dr.x*dr.z;
          StrainDerivativeTensor[CurrentSystem].bz+=force_factor*dr.y*dr.z;
          StrainDerivativeTensor[CurrentSystem].cz+=force_factor*dr.z*dr.z;

          // forces
          f.x=force_factor*dr.x;
          f.y=force_factor*dr.y;
          f.z=force_factor*dr.z;

          Framework[CurrentSystem].Atoms[f1][A].Force.x-=f.x;
          Framework[CurrentSystem].Atoms[f1][A].Force.y-=f.y;
          Framework[CurrentSystem].Atoms[f1][A].Force.z-=f.z;

          Framework[CurrentSystem].Atoms[f1][B].Force.x+=f.x;
          Framework[CurrentSystem].Atoms[f1][B].Force.y+=f.y;
          Framework[CurrentSystem].Atoms[f1][B].Force.z+=f.z;
        }
      }
    }
  }


  return 0;
}

int CalculateFrameworkIntraReplicaChargeChargeForce(void)
{
  int i,j,typeA,typeB,start;
  REAL chargeA,chargeB;
  REAL rr,energy,DF;
  VECTOR posA,posB,dr,f;
  int f1,f2,ncell,k1,k2,k3,index_j;
  int A,B;
  REAL *parms;
  REAL alpha,r;

  // Framework-Framework energy
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

        ncell=-1;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              ncell++;
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
                    PotentialGradientCoulombic(chargeA,chargeB,rr,&energy,&DF);

                    if(ncell==0)
                      UHostHostChargeChargeReal[CurrentSystem]+=energy;
                    else
                      UHostHostChargeChargeReal[CurrentSystem]+=0.5*energy;

                    f.x=-DF*dr.x;
                    f.y=-DF*dr.y;
                    f.z=-DF*dr.z;

                    Framework[CurrentSystem].Atoms[f1][i].Force.x+=f.x;
                    Framework[CurrentSystem].Atoms[f1][i].Force.y+=f.y;
                    Framework[CurrentSystem].Atoms[f1][i].Force.z+=f.z;

                    if(ncell==0)
                    {
                      Framework[CurrentSystem].Atoms[f2][j].Force.x-=f.x;
                      Framework[CurrentSystem].Atoms[f2][j].Force.y-=f.y;
                      Framework[CurrentSystem].Atoms[f2][j].Force.z-=f.z;

                      StrainDerivativeTensor[CurrentSystem].ax-=f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx-=f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx-=f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay-=f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by-=f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy-=f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az-=f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz-=f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz-=f.z*dr.z;
                    }
                    else
                    {
                      StrainDerivativeTensor[CurrentSystem].ax-=0.5*f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx-=0.5*f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx-=0.5*f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay-=0.5*f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by-=0.5*f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy-=0.5*f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az-=0.5*f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz-=0.5*f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz-=0.5*f.z*dr.z;
                    }
                  }
                }
              }
            }
      }
    }
  }


  // contributions from interactions 1-4 torsions
  // TODO: fix to handle anisotropic sites
  alpha=Alpha[CurrentSystem];
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[f1];i++)
      {
        parms=Framework[CurrentSystem].TorsionArguments[f1][i];

        if(fabs(parms[7])>1e-8)
        {
          A=Framework[CurrentSystem].Torsions[f1][i].A;
          chargeA=Framework[CurrentSystem].Atoms[f1][A].Charge;
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;

          B=Framework[CurrentSystem].Torsions[f1][i].D;
          chargeB=Framework[CurrentSystem].Atoms[f1][B].Charge;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;


          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          r=sqrt(rr);

          // note: no cutoff used here
          switch(ChargeMethod)
          {
            case NONE:
              DF=0.0;
              break;
            case SHIFTED_COULOMB:
            case TRUNCATED_COULOMB:
            case EWALD:
            default:
              UHostHostChargeChargeReal[CurrentSystem]+=parms[7]*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
              DF=parms[7]*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
              break;
          }
          StrainDerivativeTensor[CurrentSystem].ax-=DF*dr.x*dr.x;
          StrainDerivativeTensor[CurrentSystem].bx-=DF*dr.y*dr.x;
          StrainDerivativeTensor[CurrentSystem].cx-=DF*dr.z*dr.x;

          StrainDerivativeTensor[CurrentSystem].ay-=DF*dr.x*dr.y;
          StrainDerivativeTensor[CurrentSystem].by-=DF*dr.y*dr.y;
          StrainDerivativeTensor[CurrentSystem].cy-=DF*dr.z*dr.y;

          StrainDerivativeTensor[CurrentSystem].az-=DF*dr.x*dr.z;
          StrainDerivativeTensor[CurrentSystem].bz-=DF*dr.y*dr.z;
          StrainDerivativeTensor[CurrentSystem].cz-=DF*dr.z*dr.z;

          Framework[CurrentSystem].Atoms[f1][A].Force.x+=DF*dr.x;
          Framework[CurrentSystem].Atoms[f1][A].Force.y+=DF*dr.y;
          Framework[CurrentSystem].Atoms[f1][A].Force.z+=DF*dr.z;

          Framework[CurrentSystem].Atoms[f1][B].Force.x-=DF*dr.x;
          Framework[CurrentSystem].Atoms[f1][B].Force.y-=DF*dr.y;
          Framework[CurrentSystem].Atoms[f1][B].Force.z-=DF*dr.z;
        }
      }
    }
  }

  return 0;
}

void CalculateFrameworkAdsorbateReplicaVDWForce(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,scaling;
  REAL energy,force_factor;
  VECTOR posA,posB,dr,f;
  int ConnectedAtomA,ConnectedAtomB;
  REAL ReductionA,ReductionB;
  int ncell,k1,k2,k3;

  UHostAdsorbateVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
      scaling=Adsorbates[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      ReductionA=1.0;
      ConnectedAtomA=-1;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
      {
        UHostAdsorbateVDW[CurrentSystem]+=InterpolateVDWForceGrid(typeA,posA,&f);
        Adsorbates[CurrentSystem][i].Atoms[j].Force.x+=f.x;
        Adsorbates[CurrentSystem][i].Atoms[j].Force.y+=f.y;
        Adsorbates[CurrentSystem][i].Atoms[j].Force.z+=f.z;
      }
      else
      {
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
                {
                  posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
                  typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  ReductionB=1.0;
                  ConnectedAtomB=-1;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffVDWSquared)
                  {
                    typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

                    // compute the energy and force-factor
                    PotentialGradient(typeA,typeB,rr,&energy,&force_factor,scaling);

                    // energy
                    UHostAdsorbateVDW[CurrentSystem]+=energy;

                    // forces
                    f.x=force_factor*dr.x;
                    f.y=force_factor*dr.y;
                    f.z=force_factor*dr.z;

                    Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=ReductionA*f.x;
                    Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=ReductionA*f.y;
                    Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=ReductionA*f.z;

                   if(ConnectedAtomA>=0)
                   {
                     Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.x-=(1.0-ReductionA)*f.x;
                     Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.y-=(1.0-ReductionA)*f.y;
                     Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.z-=(1.0-ReductionA)*f.z;
                   }

                   if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                   {
                     Framework[CurrentSystem].Atoms[f1][k].Force.x+=ReductionB*f.x;
                     Framework[CurrentSystem].Atoms[f1][k].Force.y+=ReductionB*f.y;
                     Framework[CurrentSystem].Atoms[f1][k].Force.z+=ReductionB*f.z;

                     if(ConnectedAtomB>=0)
                     {
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.x+=(1.0-ReductionB)*f.x;
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.y+=(1.0-ReductionB)*f.y;
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.z+=(1.0-ReductionB)*f.z;
                     }
                   }

                   StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                   StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                   StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                   StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                   StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                   StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                   StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                   StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                   StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
                 }
                }
              }
              ncell++;
            }
      }
    }
  }
}

void CalculateFrameworkCationReplicaVDWForce(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,scaling;
  REAL energy,force_factor;
  VECTOR posA,posB,dr,f;
  int ConnectedAtomA,ConnectedAtomB;
  REAL ReductionA,ReductionB;
  int ncell,k1,k2,k3;

  UHostCationVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
      scaling=Cations[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      ReductionA=1.0;
      ConnectedAtomA=-1;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
      {
        UHostCationVDW[CurrentSystem]+=InterpolateVDWForceGrid(typeA,posA,&f);
        Cations[CurrentSystem][i].Atoms[j].Force.x+=f.x;
        Cations[CurrentSystem][i].Atoms[j].Force.y+=f.y;
        Cations[CurrentSystem][i].Atoms[j].Force.z+=f.z;
      }
      else
      {
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
                {
                  posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
                  typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  ReductionB=1.0;
                  ConnectedAtomB=-1;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffVDWSquared)
                  {
                    typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

                    // compute the energy and force-factor
                    PotentialGradient(typeA,typeB,rr,&energy,&force_factor,scaling);

                    // energy
                    UHostCationVDW[CurrentSystem]+=energy;

                    // forces
                    f.x=force_factor*dr.x;
                    f.y=force_factor*dr.y;
                    f.z=force_factor*dr.z;

                    Cations[CurrentSystem][i].Atoms[j].Force.x-=ReductionA*f.x;
                    Cations[CurrentSystem][i].Atoms[j].Force.y-=ReductionA*f.y;
                    Cations[CurrentSystem][i].Atoms[j].Force.z-=ReductionA*f.z;

                   if(ConnectedAtomA>=0)
                   {
                     Cations[CurrentSystem][i].Atoms[ConnectedAtomA].Force.x-=(1.0-ReductionA)*f.x;
                     Cations[CurrentSystem][i].Atoms[ConnectedAtomA].Force.y-=(1.0-ReductionA)*f.y;
                     Cations[CurrentSystem][i].Atoms[ConnectedAtomA].Force.z-=(1.0-ReductionA)*f.z;
                   }

                   if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                   {
                     Framework[CurrentSystem].Atoms[f1][k].Force.x+=ReductionB*f.x;
                     Framework[CurrentSystem].Atoms[f1][k].Force.y+=ReductionB*f.y;
                     Framework[CurrentSystem].Atoms[f1][k].Force.z+=ReductionB*f.z;

                     if(ConnectedAtomB>=0)
                     {
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.x+=(1.0-ReductionB)*f.x;
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.y+=(1.0-ReductionB)*f.y;
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.z+=(1.0-ReductionB)*f.z;
                     }
                   }

                   StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                   StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                   StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                   StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                   StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                   StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                   StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                   StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                   StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
                 }
                }
              }
              ncell++;
            }
      }
    }
  }
}

int CalculateFrameworkAdsorbateReplicaChargeChargeForce(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,chargeA,chargeB;
  REAL energy,DF,scaling;
  VECTOR posA,posB,dr,f;
  int ConnectedAtomA,ConnectedAtomB;
  REAL ReductionA,ReductionB;
  int ncell,k1,k2,k3;

  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
      scaling=Adsorbates[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
      chargeA=scaling*Adsorbates[CurrentSystem][i].Atoms[j].Charge;

      ReductionA=1.0;
      ConnectedAtomA=-1;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
      {
        UHostAdsorbateChargeChargeReal[CurrentSystem]+=InterpolateCoulombForceGrid(typeA,posA,&f);
        Adsorbates[CurrentSystem][i].Atoms[j].Force.x+=f.x;
        Adsorbates[CurrentSystem][i].Atoms[j].Force.y+=f.y;
        Adsorbates[CurrentSystem][i].Atoms[j].Force.z+=f.z;
      }
      else
      {
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
                {
                  posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
                  typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
                  chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  ReductionB=1.0;
                  ConnectedAtomB=-1;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  {
                    PotentialGradientCoulombic(chargeA,chargeB,rr,&energy,&DF);

                    // energy
                    UHostAdsorbateChargeChargeReal[CurrentSystem]+=energy;

                    // forces
                    f.x=DF*dr.x;
                    f.y=DF*dr.y;
                    f.z=DF*dr.z;

                    Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=ReductionA*f.x;
                    Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=ReductionA*f.y;
                    Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=ReductionA*f.z;

                   if(ConnectedAtomA>=0)
                   {
                     Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.x-=(1.0-ReductionA)*f.x;
                     Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.y-=(1.0-ReductionA)*f.y;
                     Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.z-=(1.0-ReductionA)*f.z;
                   }

                   if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                   {
                     Framework[CurrentSystem].Atoms[f1][k].Force.x+=ReductionB*f.x;
                     Framework[CurrentSystem].Atoms[f1][k].Force.y+=ReductionB*f.y;
                     Framework[CurrentSystem].Atoms[f1][k].Force.z+=ReductionB*f.z;

                     if(ConnectedAtomB>=0)
                     {
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.x+=(1.0-ReductionB)*f.x;
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.y+=(1.0-ReductionB)*f.y;
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.z+=(1.0-ReductionB)*f.z;
                     }
                   }

                   StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                   StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                   StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                   StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                   StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                   StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                   StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                   StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                   StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
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

int CalculateFrameworkCationReplicaChargeChargeForce(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,chargeA,chargeB;
  REAL energy,DF,scaling;
  VECTOR posA,posB,dr,f;
  int ConnectedAtomA,ConnectedAtomB;
  REAL ReductionA,ReductionB;
  int ncell,k1,k2,k3;

  UHostCationChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
      scaling=Cations[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
      chargeA=scaling*Cations[CurrentSystem][i].Atoms[j].Charge;

      ReductionA=1.0;
      ConnectedAtomA=-1;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
      {
        UHostCationChargeChargeReal[CurrentSystem]+=InterpolateCoulombForceGrid(typeA,posA,&f);
        Cations[CurrentSystem][i].Atoms[j].Force.x+=f.x;
        Cations[CurrentSystem][i].Atoms[j].Force.y+=f.y;
        Cations[CurrentSystem][i].Atoms[j].Force.z+=f.z;
      }
      else
      {
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
                {
                  posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
                  typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
                  chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  ReductionB=1.0;
                  ConnectedAtomB=-1;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  {
                    PotentialGradientCoulombic(chargeA,chargeB,rr,&energy,&DF);

                    // energy
                    UHostCationChargeChargeReal[CurrentSystem]+=energy;

                    // forces
                    f.x=DF*dr.x;
                    f.y=DF*dr.y;
                    f.z=DF*dr.z;

                    Cations[CurrentSystem][i].Atoms[j].Force.x-=ReductionA*f.x;
                    Cations[CurrentSystem][i].Atoms[j].Force.y-=ReductionA*f.y;
                    Cations[CurrentSystem][i].Atoms[j].Force.z-=ReductionA*f.z;

                   if(ConnectedAtomA>=0)
                   {
                     Cations[CurrentSystem][i].Atoms[ConnectedAtomA].Force.x-=(1.0-ReductionA)*f.x;
                     Cations[CurrentSystem][i].Atoms[ConnectedAtomA].Force.y-=(1.0-ReductionA)*f.y;
                     Cations[CurrentSystem][i].Atoms[ConnectedAtomA].Force.z-=(1.0-ReductionA)*f.z;
                   }

                   if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                   {
                     Framework[CurrentSystem].Atoms[f1][k].Force.x+=ReductionB*f.x;
                     Framework[CurrentSystem].Atoms[f1][k].Force.y+=ReductionB*f.y;
                     Framework[CurrentSystem].Atoms[f1][k].Force.z+=ReductionB*f.z;

                     if(ConnectedAtomB>=0)
                     {
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.x+=(1.0-ReductionB)*f.x;
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.y+=(1.0-ReductionB)*f.y;
                       Framework[CurrentSystem].Atoms[f1][ConnectedAtomB].Force.z+=(1.0-ReductionB)*f.z;
                     }
                   }

                   StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                   StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                   StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                   StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                   StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                   StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                   StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                   StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                   StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
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
