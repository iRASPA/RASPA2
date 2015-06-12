/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2015 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_energy.c' is part of RASPA.

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
#include "simulation.h"
#include "ewald.h"
#include "potentials.h"
#include "utils.h"
#include "output.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "spacegroup.h"

// the routines for VDW and charge-charge interaction are splitted. The VDW interactions can potentially shift the
// actual interaction site across a bond (see e.g. MM3).

// Calculate the Energy of a bead of type "typeA" at a position "posA".
// posA is defined in xyz (also for triclinic).
REAL CalculateFrameworkVDWEnergyAtPosition(POINT posA,int typeA,REAL scaling)
{
  int i,j,typeB,f1;
  REAL rr;
  REAL UVDW;
  VECTOR posB,dr,s;
  int icell,icell0;
  int ncell;

  UVDW=0.0;

  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))   // grid-interpolating for rigid frameworks
    {
      UVDW=InterpolateVDWGrid(typeA,posA);
      return UVDW;
    }
    else if(UseCellLists[CurrentSystem]) // energy using cell-lists
    {
      // convert from xyz to abc
      s.x=InverseBox[CurrentSystem].ax*posA.x+InverseBox[CurrentSystem].bx*posA.y+InverseBox[CurrentSystem].cx*posA.z;
      s.y=InverseBox[CurrentSystem].ay*posA.x+InverseBox[CurrentSystem].by*posA.y+InverseBox[CurrentSystem].cy*posA.z;
      s.z=InverseBox[CurrentSystem].az*posA.x+InverseBox[CurrentSystem].bz*posA.y+InverseBox[CurrentSystem].cz*posA.z;

      // apply boundary condition
      s.x-=(REAL)NINT(s.x);
      s.y-=(REAL)NINT(s.y);
      s.z-=(REAL)NINT(s.z);

      // s between 0 and 1
      s.x+=0.5;
      s.y+=0.5;
      s.z+=0.5;

      // compute the corresponding cell-id
      icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
             ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
             ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        // loop over cells
        for(i=0;i<27;i++)
        {
          icell=CellListMap[CurrentSystem][icell0][i];

          j=Framework[CurrentSystem].CellListHead[f1][icell];

          while(j>=0)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
            posB=Framework[CurrentSystem].Atoms[f1][j].AnisotropicPosition;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
              UVDW+=PotentialValue(typeA,typeB,rr,scaling);

            j=Framework[CurrentSystem].CellList[f1][j];
          }
        }
      }
    }
    else if(UseReplicas[CurrentSystem]) // energy using replicas
    {
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][i].Type;
          posB=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
              UVDW+=PotentialValue(typeA,typeB,rr,scaling);
          }
        }
      }
    }
    else // regular energy calculation
    {
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][i].Type;
          posB=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffVDWSquared)
            UVDW+=PotentialValue(typeA,typeB,rr,scaling);
        }
      }
    }
  }
  return UVDW;
}


void CalculateFrameworkChargeEnergyAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UchargeBondDipole,REAL scaling)
{
  int i,typeB,f1,k;
  int B1,B2;
  REAL r,rr,chargeA,chargeB;
  REAL DipoleMagnitudeB,temp;
  VECTOR posB,posB1,posB2,dr,dipoleB;
  int ncell;

  (*UChargeCharge)=0.0;
  (*UchargeBondDipole)=0.0;

  if(ChargeMethod==NONE) return;

  if(!PseudoAtoms[typeA].HasCharges) return;

  chargeA=scaling*PseudoAtoms[typeA].Charge1;
  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
    {
      (*UChargeCharge)=InterpolateCoulombGrid(typeA,posA);
    }
    else if(UseReplicas[CurrentSystem])
    {
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][i].Type;
          posB=Framework[CurrentSystem].Atoms[f1][i].Position;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              r=sqrt(rr);
              chargeB=Framework[CurrentSystem].Atoms[f1][i].Charge;

              (*UChargeCharge)+=PotentialValueCoulombic(chargeA,chargeB,r);
            }
          }
        }

        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyReplicaBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
              (*UchargeBondDipole)+=PotentialValueChargeBondDipole(chargeA,dipoleB,dr,sqrt(rr));
          }
        }
      }
    }
    else
    {
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][i].Type;
          posB=Framework[CurrentSystem].Atoms[f1][i].Position;
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            r=sqrt(rr);
            chargeB=Framework[CurrentSystem].Atoms[f1][i].Charge;

            (*UChargeCharge)+=PotentialValueCoulombic(chargeA,chargeB,r);
          }
        }

        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
            (*UchargeBondDipole)+=PotentialValueChargeBondDipole(chargeA,dipoleB,dr,sqrt(rr));
        }
      }
    }
  }
}

void CalculateFrameworkBondDipoleEnergyAtPosition(VECTOR posA1,VECTOR posA2,REAL DipoleMagnitudeA,REAL *UchargeBondDipole,REAL *UBondDipoleBondDipole)
{
  int k,f1;
  int B1,B2;
  int TypeB;
  VECTOR posA,posB,posB1,posB2,dr;
  REAL DipoleMagnitudeB,chargeB;
  REAL ri2,rk2,rr,length,temp;
  VECTOR dipoleA,dipoleB;
  int ncell;

  (*UchargeBondDipole)=0.0;
  (*UBondDipoleBondDipole)=0.0;

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

  if(UseReplicas[CurrentSystem])
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
      {
        TypeB=Framework[CurrentSystem].Atoms[f1][k].Type;

        if(PseudoAtoms[TypeB].HasCharges)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].x);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
              (*UchargeBondDipole)-=PotentialValueChargeBondDipole(chargeB,dipoleA,dr,sqrt(rr));
          }
        }
      }

      for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
      {
        DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
        B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
        B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
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

        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
          dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
          dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
          dr=ApplyReplicaBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffBondDipoleBondDipoleSquared)
            (*UBondDipoleBondDipole)+=PotentialValueBondDipoleBondDipole(dipoleA,dipoleB,dr,sqrt(rr));
        }
      }
    }
  }
  else
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
      {
        TypeB=Framework[CurrentSystem].Atoms[f1][k].Type;

        if(PseudoAtoms[TypeB].HasCharges)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
            (*UchargeBondDipole)-=PotentialValueChargeBondDipole(chargeB,dipoleA,dr,sqrt(rr));
        }
      }

      for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
      {
        DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
        B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
        B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
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
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffBondDipoleBondDipoleSquared)
          (*UBondDipoleBondDipole)+=PotentialValueBondDipoleBondDipole(dipoleA,dipoleB,dr,sqrt(rr));
      }
    }
  }
}


// Monte-Carlo full energy routines
// ================================

REAL CalculateFrameworkBondEnergy(int flag,int f2,int atom_id)
{
  int i,f1;
  REAL r,rr,r1,U;
  REAL temp,temp2,exp_term;
  REAL *parms,UHostBond;
  VECTOR posA,posB,dr;
  int A,B;

  UHostBond=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBonds[f1];i++)
      {
        A=Framework[CurrentSystem].Bonds[f1][i].A;
        B=Framework[CurrentSystem].Bonds[f1][i].B;

        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id)))
        {
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
              break;
            case CORE_SHELL_SPRING:
              U=0.5*parms[0]*SQR(r);
              break;
            case MORSE_BOND:
              // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
              // ===============================================
              // p_0/k_B [K]       force constant
              // p_1     [A^-1]    parameter
              // p_2     [A]       reference bond distance
              temp=exp(parms[1]*(parms[2]-r));
              U=parms[0]*(SQR(1.0-temp)-1.0);
              break;
            case LJ_12_6_BOND:
              // A/r_ij^12-B/r_ij^6
              // ===============================================
              // p_0/k_B [K A^12]
              // p_1/k_B [K A^6]
              temp=CUBE(1.0/rr);
              U=parms[0]*SQR(temp)-parms[1]*temp;
              break;
            case LENNARD_JONES_BOND:
              // 4*p_0*((p_1/r)^12-(p_1/r)^6)
              // ===============================================
              // p_0/k_B [K]
              // p_1     [A]
              temp=CUBE(parms[1]/rr);
              U=4.0*parms[0]*(temp*(temp-1.0));
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
              break;
            case MM3_BOND:
              // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
              // =================================================================
              // p_0     [mdyne/A molecule]
              // p_1     [A]
              temp=r-parms[1];
              temp2=SQR(temp);
              U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
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
              fprintf(stderr, "Undefined Bond potential in routine 'CalculateFrameworkBondEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }

          // add contribution to the Adsorbate stretch energy
          UHostBond+=U;
        }
      }
    }
  }
  return UHostBond;
}

REAL CalculateFrameworkUreyBradleyEnergy(int flag,int f2,int atom_id)
{
  int i,f1;
  REAL r,rr,r1,U;
  REAL temp,temp2,exp_term;
  VECTOR dr,posA,posC;
  REAL *parms,UHostUreyBradley;
  int A,C;

  UHostUreyBradley=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfUreyBradleys[f1];i++)
      {
        A=Framework[CurrentSystem].UreyBradleys[f1][i].A;
        C=Framework[CurrentSystem].UreyBradleys[f1][i].C;

        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||C==atom_id)))
        {
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
              break;
            case MORSE_UREYBRADLEY:
              // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
              // ===============================================
              // p_0/k_B [K]       force constant
              // p_1     [A^-1]    parameter
              // p_2     [A]       reference bond distance
              temp=exp(parms[1]*(parms[2]-r));
              U=parms[0]*(SQR(1.0-temp)-1.0);
              break;
            case LJ_12_6_UREYBRADLEY:
              // A/r_ij^12-B/r_ij^6
              // ===============================================
              // p_0/k_B [K A^12]
              // p_1/k_B [K A^6]
              temp=CUBE(1.0/rr);
              U=parms[0]*SQR(temp)-parms[1]*temp;
              break;
            case LENNARD_JONES_UREYBRADLEY:
              // 4*p_0*((p_1/r)^12-(p_1/r)^6)
              // ===============================================
              // p_0/k_B [K]
              // p_1     [A]
              temp=CUBE(parms[1]/rr);
              U=4.0*parms[0]*(temp*(temp-1.0));
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
              break;
            case MM3_UREYBRADLEY:
              // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
              // =================================================================
              // p_0     [mdyne/A molecule]
              // p_1     [A]
              temp=r-parms[1];
              temp2=SQR(temp);
              U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
              break;
            case RIGID_UREYBRADLEY:
              U=0.0;
              break;
            case FIXED_UREYBRADLEY:
              U=0.0;
              break;
            default:
              fprintf(stderr, "Undefined Urey-Bradley potential in routine 'CalculateFrameworkUreyBradleyEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }

          // add contribution to the Adsorbate Urey-Bradley energy
          UHostUreyBradley+=U;
        }
      }
    }
  }
  return UHostUreyBradley;
}

REAL CalculateFrameworkBendEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  REAL *parms,U,temp,temp2;
  REAL CosTheta,Theta,SinTheta;
  REAL rab,rbc,rac,DTDX,UHostBend;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL delta,rt2,rap2,rcp2;

  UHostBend=0.0;
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

        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&((A==atom_id)||(B==atom_id)||(C==atom_id)||(D==atom_id))))
        {
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
              break;
            case CORE_SHELL_BEND:
              // (1/2)p_0*(theta-p_1)^2
              // ===============================================
              // p_0/k_B [K/rad^2]
              // p_1     [degrees]
              U=0.5*parms[0]*SQR(Theta-parms[1]);
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
              break;
            case HARMONIC_COSINE_BEND:
              // (1/2)*p_0*(cos(theta)-cos(p_1))^2
              // ===============================================
              // p_0/k_B [K]
              // p_1     [degrees]
              U=0.5*parms[0]*SQR(CosTheta-parms[1]);
              break;
            case COSINE_BEND:
              // p_0*(1+cos(p_1*theta-p_2))
              // ===============================================
              // p_0/k_B [K]
              // p_1     [-]
              // p_2     [degrees]
              temp=parms[1]*Theta-parms[2];
              U=parms[0]*(1.0+cos(temp));
              break;
            case TAFIPOLSKY_BEND:
              // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
              // ===============================================
              // p_0/k_B [K]
              U=0.5*parms[0]*(1.0+cos(Theta))*(1.0+cos(2.0*Theta));
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

          // add contribution to the energy
          UHostBend+=U;
        }
      }
    }
  }
  return UHostBend;
}

REAL CalculateFrameworkInversionBendEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2;
  REAL CosChi,Chi,energy;
  REAL temp,temp2,UHostInversionBend;
  VECTOR Rab,Rbc,Rbd,Rcd,Rad;
  POINT posA,posB,posC,posD;

  UHostInversionBend=0.0;
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

        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
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
              break;
            case HARMONIC_COSINE_INVERSION:
            case HARMONIC_COSINE_INVERSION2:
              // (1/2)*p_0*(cos(phi)-cos(p_1))^2
              // ===============================================
              // p_0/k_B [K]
              // p_1     [degrees]
              energy=0.5*parms[0]*SQR(CosChi-parms[1]);
              break;
            case PLANAR_INVERSION:
            case PLANAR_INVERSION2:
              // (1/2)*p_0*(1-cos(phi))
              // ===============================================
              // p_0/k_B [K]
              energy=parms[0]*(1.0-CosChi);
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
              break;
            default:
              fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }

          // energy
          UHostInversionBend+=energy;
        }
      }
    }
  }
  return UHostInversionBend;
}

REAL CalculateFrameworkTorsionEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL rbc,U;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi,UHostTorsion;
  VECTOR Pb,Pc;
  REAL *parms;

  UHostTorsion=0.0;
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

        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
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

          parms=Framework[CurrentSystem].TorsionArguments[f1][i];

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
              break;
            case HARMONIC_COSINE_DIHEDRAL:
              // (1/2)*p_0*(cos(phi)-cos(p_1))^2
              // ===============================================
              // p_0/k_B [K]
              // p_1     [degrees]
              U=0.5*parms[0]*SQR(CosPhi-cos(parms[1]));
              break;
            case THREE_COSINE_DIHEDRAL:
              // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
              // ========================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
              break;
            case MM3_DIHEDRAL:
              // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
              // ========================================================================
              // p_0     [kcal/mol]
              // p_1     [kcal/mol]
              // p_2     [kcal/mol]
              U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
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
              break;
            case CFF_DIHEDRAL:
              // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
              // ======================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
              break;
            case CFF_DIHEDRAL2:
              // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
              // ======================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
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
              break;
            case TRAPPE_DIHEDRAL:
              // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
              // ==========================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              // p_3/k_B [K]
              U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
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
              break;
            case OPLS_DIHEDRAL:
              // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
              // =================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              // p_3/k_B [K]
              U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
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
              break;
            default:
              fprintf(stderr, "Undefined Torsion potential in routine 'CalculateFrameworkTorsionEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }

          // energy
          UHostTorsion+=U;
        }
      }
    }
  }
  return UHostTorsion;
}

REAL CalculateFrameworkImproperTorsionEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL rbc,U;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi,UHostImproperTorsion;
  VECTOR Pb,Pc;
  REAL *parms;

  UHostImproperTorsion=0.0;
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

        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
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

          parms=Framework[CurrentSystem].ImproperTorsionArguments[f1][i];

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
              break;
            case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
              // (1/2)*p_0*(cos(phi)-cos(p_1))^2
              // ===============================================
              // p_0/k_B [K]
              // p_1     [degrees]
              U=0.5*parms[0]*SQR(CosPhi-cos(parms[1]));
              break;
            case THREE_COSINE_IMPROPER_DIHEDRAL:
              // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
              // ========================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
              break;
            case MM3_IMPROPER_DIHEDRAL:
              // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
              // ========================================================================
              // p_0     [kcal/mol]
              // p_1     [kcal/mol]
              // p_2     [kcal/mol]
              U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
              break;
            case CFF_IMPROPER_DIHEDRAL:
              // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
              // ======================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
              break;
            case CFF_IMPROPER_DIHEDRAL2:
              // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
              // ======================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
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
              break;
            case TRAPPE_IMPROPER_DIHEDRAL:
              // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
              // ==========================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              // p_3/k_B [K]
              U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
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
              break;
            case OPLS_IMPROPER_DIHEDRAL:
              // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
              // =================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              // p_3/k_B [K]
              U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
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
              break;
            case FIXED_IMPROPER_DIHEDRAL:
              U=0.0;
              break;
            default:
              fprintf(stderr, "Undefined Improper-Torsion potential in routine 'CalculateFrameworkImproperTorsionEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }

          // energy
          UHostImproperTorsion+=U;
        }
      }
    }
  }
  return UHostImproperTorsion;
}

REAL CalculateFrameworkBondBondEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,f1;
  REAL *parms;
  REAL energy;
  REAL rab,rbc,UHostBondBond;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc;

  UHostBondBond=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondBonds[f1];i++)
      {
        A=Framework[CurrentSystem].BondBonds[f1][i].A;
        B=Framework[CurrentSystem].BondBonds[f1][i].B;
        C=Framework[CurrentSystem].BondBonds[f1][i].C;

        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id)))
        {
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
            break;
            default:
              fprintf(stderr, "Undefined Bond-Bond potential in routine 'CalculateFrameworkBondBondEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }

          // add contribution to the energy
          UHostBondBond+=energy;
        }
      }
    }
  }
  return UHostBondBond;
}

REAL CalculateFrameworkBondBendEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,f1;
  REAL *parms,U;
  REAL cost,theta,sint;
  REAL rab,rbc,rac,UHostBondBend;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc,Rac;

  UHostBondBend=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondBends[f1];i++)
      {
        A=Framework[CurrentSystem].BondBends[f1][i].A;
        B=Framework[CurrentSystem].BondBends[f1][i].B;
        C=Framework[CurrentSystem].BondBends[f1][i].C;

        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id)))
        {
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
          theta=acos(cost);
          sint=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(cost)));

          parms=Framework[CurrentSystem].BondBendArguments[f1][i];

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
              U=(theta-parms[0])*(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]));
              break;
            case MM3_BOND_BEND_CROSS:
              // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
              // =====================================
              // p_0     [mdyne/rad]
              // p_1     [A]
              // p_2     [A]
              // p_3     [degrees]
              U=parms[0]*((rab-parms[1])+(rbc-parms[2]))*RAD2DEG*(theta-parms[3]);
              break;
            case TRUNCATED_HARMONIC:
              // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
              // ================================================================
              // p_0/k_B [K/rad^2]
              // p_1     [degrees]
              // p_2     [A]
              U=0.5*parms[0]*SQR(theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8));
              break;
            case SCREENED_HARMONIC:
              // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
              // ===============================================
              // p_0/k_B [K/rad^2]
              // p_1     [degrees]
              // p_2     [A]
              // p_3     [A]
              U=0.5*parms[0]*SQR(theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3]));
              break;
            case SCREENED_VESSAL:
              // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
              // ============================================================================
              // p_0/k_B [K/rad^2]
              // p_1     [degrees]
              // p_2     [A]
              // p_3     [A]
              U=(parms[0]/(8.0*SQR(theta-M_PI)))*SQR(SQR(parms[1]-M_PI)-SQR(theta-M_PI))
                    *exp(-(rab/parms[2]+rbc/parms[3]));
              break;
            case TRUNCATED_VESSAL:
              // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
              //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
              // ============================================================================
              // p_0/k_B [K/rad^(4+p_2)]
              // p_1     [degrees]
              // p_2     [-]
              // p_3     [A]
              U=parms[0]*(pow(theta,parms[2])*SQR(theta-parms[1])*SQR(theta+parms[1]-2.0*M_PI)
                    -0.5*parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*SQR(theta-parms[1])*pow(M_PI-parms[1],(REAL)3.0))
                    *exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8));
              break;
            default:
              fprintf(stderr, "Undefined Bond-Bend potential in routine 'CalculateFrameworkBondBendEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }

          // energy
          UHostBondBend+=U;
        }
      }
    }
  }
  return UHostBondBend;
}

// Bend/Bend cross term for a centered atom B
// the first angle is A-B-C
// the second angle is A-B-D
REAL CalculateFrameworkBendBendEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL U,rab,rbc,rbd;
  VECTOR Dab,Dbc,Dbd;
  REAL dot_abc,dot_abd,UHostBendBend;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL *parms;

  UHostBendBend=0.0;
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

        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
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

          parms=(REAL*)&Framework[CurrentSystem].BendBendArguments[f1][i];

          switch(Framework[CurrentSystem].BendBendType[f1][i])
          {
            case CVFF_BEND_BEND_CROSS:
            case CFF_BEND_BEND_CROSS:
              // p_0*(Theta1-p_1)*(Theta2-p_2)
              // ===================================
              // p_0/k_B [K/rad^2)]
              // p_1     [degrees]
              // p_2     [degrees]
              U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
              break;
            case MM3_BEND_BEND_CROSS:
              // -p_0*(Theta1-p_1)*(Theta2-p_2)
              // ===================================
              // p_0     [mdyne A/rad^2]
              // p_1     [degrees]
              // p_2     [degrees]
              U=-parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])*(Theta2-parms[2]);
              break;
            default:
              fprintf(stderr, "Undefined Bend-Bend potential in routine 'CalculateFrameworkBendBendEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
          // energy
          UHostBendBend+=U;
        }
      }
    }
  }
  return UHostBendBend;
}

REAL CalculateFrameworkBondTorsionEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL U,rab,rbc,rcd;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,temp;
  REAL CosPhi,CosPhi2;
  REAL *parms,UHostBondTorsion;

  UHostBondTorsion=0.0;
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

        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
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

          parms=(REAL*)&Framework[CurrentSystem].BondTorsionArguments[f1][i];

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
              U=parms[0]*temp*CosPhi+parms[1]*temp*(2.0*CosPhi2-1.0)+parms[2]*temp*(4.0*CosPhi2*CosPhi-3.0*CosPhi);
              break;
            default:
              fprintf(stderr, "Undefined Bond-Torsion potential in routine 'CalculateFrameworkBondTorsionEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }

          // energy
          UHostBondTorsion+=U;
        }
      }
    }
  }
  return UHostBondTorsion;
}

REAL CalculateFrameworkBendTorsionEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL U,rab,rbc,rcd;
  VECTOR Dab,Dbc,Dcd,dr,ds,Pb,Pc;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL sign,Phi,SinPhi;
  REAL *parms,UHostBendTorsion;

  UHostBendTorsion=0.0;
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

        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
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

          parms=(REAL*)&Framework[CurrentSystem].BendTorsionArguments[f1][i];

          switch(Framework[CurrentSystem].BendTorsionType[f1][i])
          {
            case CVFF_BEND_TORSION_CROSS:
            case CFF_BEND_TORSION_CROSS:
              // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
              // =====================================================================================
              // p_0/k_B [K/rad^3]
              // p_1     [degrees]
              // p_2     [degrees]
              U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi;
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
              U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2);
              break;
            case SMOOTHED_THREE_COSINE_DIHEDRAL:
              // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                     Smoothing(Theta1)*Smoothing(Theta2);
              break;
            case NICHOLAS_DIHEDRAL:
              // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                     Smoothing(Theta1);
              break;
            case SMOOTHED_CFF_DIHEDRAL:
              // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
              break;
            case SMOOTHED_CFF_DIHEDRAL2:
              // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*Smoothing(Theta1)*Smoothing(Theta2);
              break;
            case SMOOTHED_CFF_BEND_TORSION_CROSS:
              // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K/rad^3]
              // p_1     [degrees]
              // p_2     [degrees]
              U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi*Smoothing(Theta1)*Smoothing(Theta2);
              break;
            default:
              fprintf(stderr, "Undefined Bend-Torsion potential in routine 'CalculateFrameworkBendTorsionEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }

          // energy
          UHostBendTorsion+=U;
        }
      }
    }
  }
  return UHostBendTorsion;
}

int CalculateFrameworkIntraVDWEnergy(void)
{
  int i,j,typeA,typeB,A,B;
  REAL rr;
  VECTOR posA,posB,dr;
  int f1,f2;
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

        for(j=i+1;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],0))
          {
            typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
            posB=Framework[CurrentSystem].Atoms[f1][j].AnisotropicPosition;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
              UHostHostVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr,1.0);
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

        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
        {
          typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
          posB=Framework[CurrentSystem].Atoms[f2][j].AnisotropicPosition;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffVDWSquared)
            UHostHostVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr,1.0);
        }
      }
    }
  }

  // contributions from interactions 1-4 torsions
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
          posA=Framework[CurrentSystem].Atoms[f1][A].AnisotropicPosition;

          B=Framework[CurrentSystem].Torsions[f1][i].D;
          typeB=Framework[CurrentSystem].Atoms[f1][B].Type;
          posB=Framework[CurrentSystem].Atoms[f1][B].AnisotropicPosition;


          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          UHostHostVDW[CurrentSystem]+=parms[6]*PotentialValue(typeA,typeB,rr,1.0);
        }
      }
    }
  }

  return 0;
}

int CalculateFrameworkIntraChargeChargeEnergy(void)
{
  int i,j,typeA,typeB;
  REAL chargeA,chargeB;
  REAL r,rr;
  VECTOR posA,posB,dr;
  int f1,f2,A,B;
  REAL *parms;
  REAL alpha;

  // Framework-Framework energy
  UHostHostChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

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

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              r=sqrt(rr);
              chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

              UHostHostChargeChargeReal[CurrentSystem]+=PotentialValueCoulombic(chargeA,chargeB,r);
            }
          }
        }
      }
    }
  }


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

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            r=sqrt(rr);
            chargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;

            UHostHostChargeChargeReal[CurrentSystem]+=PotentialValueCoulombic(chargeA,chargeB,r);
          }
        }
      }
    }
  }

  // contributions from interactions 1-4 torsions
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
          posA=Framework[CurrentSystem].Atoms[f1][A].AnisotropicPosition;

          B=Framework[CurrentSystem].Torsions[f1][i].D;
          chargeB=Framework[CurrentSystem].Atoms[f1][B].Charge;
          posB=Framework[CurrentSystem].Atoms[f1][B].AnisotropicPosition;


          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

          UHostHostChargeChargeReal[CurrentSystem]+=parms[7]*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraChargeBondDipoleEnergy(void)
{
  int i,j,f1,f2;
  int A1,A2;
  int Type;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL rr,ri2,temp,length,chargeB;
  VECTOR dipoleA;
  REAL Bt1;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;

  Bt1=0.0;
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  UHostHostChargeBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

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
                UHostHostChargeBondDipoleReal[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA,dr,sqrt(rr));
            }
          }
        }
      }
    }
  }

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
                UHostHostChargeBondDipoleReal[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA,dr,sqrt(rr));
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraBondDipoleBondDipoleEnergy(void)
{
  int i,j,f1,f2;
  int A1,A2,B1,B2;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL rr,ri2,rk2,temp,length;
  VECTOR dipoleA,dipoleB;
  REAL Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL_MATRIX3x3 v;

  Bt1=Bt2=0.0;
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  UHostHostBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

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
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffBondDipoleBondDipoleSquared)
              UHostHostBondDipoleBondDipoleReal[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA,dipoleB,dr,sqrt(rr));
          }
        }
      }
    }
  }

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

          // if framework are different, or if they are the same but not excluded within the framework
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
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffBondDipoleBondDipoleSquared)
            UHostHostBondDipoleBondDipoleReal[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA,dipoleB,dr,sqrt(rr));
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateVDWEnergy(void)
{
  int i,j,k,l,f1;
  int typeA,typeB,type;
  REAL rr,scalingA,energy;
  VECTOR posA,posB,dr;
  VECTOR s;
  int icell0,icell;

  UHostAdsorbateVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
      scalingA=Adsorbates[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA])&&(!IsFractionalAdsorbateMolecule(i)))
      {
        UHostAdsorbateVDW[CurrentSystem]+=InterpolateVDWGrid(typeA,posA);
      }
      else
      {
        if(UseCellLists[CurrentSystem])
        {
          // convert from xyz to abc
          s.x=InverseBox[CurrentSystem].ax*posA.x+InverseBox[CurrentSystem].bx*posA.y+InverseBox[CurrentSystem].cx*posA.z;
          s.y=InverseBox[CurrentSystem].ay*posA.x+InverseBox[CurrentSystem].by*posA.y+InverseBox[CurrentSystem].cy*posA.z;
          s.z=InverseBox[CurrentSystem].az*posA.x+InverseBox[CurrentSystem].bz*posA.y+InverseBox[CurrentSystem].cz*posA.z;

          // apply boundary condition
          s.x-=(REAL)NINT(s.x);
          s.y-=(REAL)NINT(s.y);
          s.z-=(REAL)NINT(s.z);

          // s between 0 and 1
          s.x+=0.5;
          s.y+=0.5;
          s.z+=0.5;

          // compute the corresponding cell-id
          icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
                 ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
                 ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            // loop over cells
            for(l=0;l<27;l++)
            {
              icell=CellListMap[CurrentSystem][icell0][l];

              k=Framework[CurrentSystem].CellListHead[f1][icell];

              while(k>=0)
              {
                posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
                typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                  UHostAdsorbateVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr,1.0);

                k=Framework[CurrentSystem].CellList[f1][k];
              }
            }
          }
        }
        else
        {
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
                energy=PotentialValue(typeA,typeB,rr,scalingA);;
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UHostAdsorbateVDW[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationVDWEnergy(void)
{
  int i,j,k,l,f1;
  int typeA,typeB,type;
  REAL rr,scalingA,energy;
  VECTOR posA,posB,dr;
  VECTOR s;
  int icell0,icell;

  UHostCationVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
      scalingA=Cations[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA])&&(!IsFractionalCationMolecule(i)))
      {
        UHostCationVDW[CurrentSystem]+=scalingA*InterpolateVDWGrid(typeA,posA);
      }
      else
      {
        if(UseCellLists[CurrentSystem])
        {
          // convert from xyz to abc
          s.x=InverseBox[CurrentSystem].ax*posA.x+InverseBox[CurrentSystem].bx*posA.y+InverseBox[CurrentSystem].cx*posA.z;
          s.y=InverseBox[CurrentSystem].ay*posA.x+InverseBox[CurrentSystem].by*posA.y+InverseBox[CurrentSystem].cy*posA.z;
          s.z=InverseBox[CurrentSystem].az*posA.x+InverseBox[CurrentSystem].bz*posA.y+InverseBox[CurrentSystem].cz*posA.z;

          // apply boundary condition
          s.x-=(REAL)NINT(s.x);
          s.y-=(REAL)NINT(s.y);
          s.z-=(REAL)NINT(s.z);

          // s between 0 and 1
          s.x+=0.5;
          s.y+=0.5;
          s.z+=0.5;

          // compute the corresponding cell-id
          icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
                 ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
                 ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            // loop over cells
            for(l=0;l<27;l++)
            {
              icell=CellListMap[CurrentSystem][icell0][l];

              k=Framework[CurrentSystem].CellListHead[f1][icell];

              while(k>=0)
              {
                posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
                typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                  UHostCationVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr,1.0);

                k=Framework[CurrentSystem].CellList[f1][k];
              }
            }
          }
        }
        else
        {
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
                energy=PotentialValue(typeA,typeB,rr,scalingA);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UHostCationVDW[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateChargeChargeEnergy(void)
{
  int i,j,k;
  int typeA,typeB,type;
  REAL r,rr;
  REAL chargeA,chargeB;
  VECTOR posA,posB,dr;
  REAL scalingA;

  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      scalingA=Adsorbates[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
      chargeA=scalingA*Adsorbates[CurrentSystem][i].Atoms[j].Charge;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid)&&(!IsFractionalAdsorbateMolecule(i)))
      {
        UHostAdsorbateChargeChargeReal[CurrentSystem]+=InterpolateCoulombGrid(typeA,posA);
      }
      else
      {
        for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
          {
            posB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;
            typeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              r=sqrt(rr);
              chargeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Charge;

              UHostAdsorbateChargeChargeReal[CurrentSystem]+=PotentialValueCoulombic(chargeA,chargeB,r);
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationChargeChargeEnergy(void)
{
  int i,j,k;
  int typeA,typeB,type;
  REAL r,rr;
  REAL chargeA,chargeB;
  VECTOR posA,posB,dr;
  REAL scalingA;

  UHostCationChargeChargeReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      scalingA=Cations[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
      chargeA=scalingA*Cations[CurrentSystem][i].Atoms[j].Charge;
      posA=Cations[CurrentSystem][i].Atoms[j].Position;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid)&&(!IsFractionalCationMolecule(i)))
      {
        UHostCationChargeChargeReal[CurrentSystem]+=InterpolateCoulombGrid(typeA,posA);
      }
      else
      {
        for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
          {
            posB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;
            typeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared[CurrentSystem])
            {
              r=sqrt(rr);
              chargeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Charge;

              UHostCationChargeChargeReal[CurrentSystem]+=PotentialValueCoulombic(chargeA,chargeB,r);
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateChargeBondDipoleEnergy(void)
{
  int i,j,k,f1;
  int A1,A2;
  int Type,TypeA;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL rr,ri2,temp,length,chargeB;
  VECTOR dipoleA;
  REAL Bt0,Bt1;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;

  UHostAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

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
              UHostAdsorbateChargeBondDipoleReal[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA,dr,sqrt(rr));
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

      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
        {
          posB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;
          Type=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;
          chargeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
            UHostAdsorbateChargeBondDipoleReal[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA,dr,sqrt(rr));
        }
      }
    }
  }

  return 0;
}

int CalculateFrameworkCationChargeBondDipoleEnergy(void)
{
  int i,j,k,f1;
  int A1,A2;
  int Type,TypeA;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL rr,ri2,temp,length,chargeB;
  VECTOR dipoleA;
  REAL Bt0,Bt1;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;

  UHostCationChargeBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

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
              UHostCationChargeBondDipoleReal[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA,dr,sqrt(rr));
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

      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
        {
          posB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;
          Type=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;
          chargeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
            UHostCationChargeBondDipoleReal[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA,dr,sqrt(rr));
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergy(void)
{
  int i,k,l,f1;
  int A1,A2,B1,B2;
  int TypeB;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r2,ri2,rk2,temp,length;
  VECTOR dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;

  UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

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
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(r2<CutOffBondDipoleBondDipoleSquared)
            UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA,dipoleB,dr,sqrt(r2));
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationBondDipoleBondDipoleEnergy(void)
{
  int i,k,l,f1;
  int A1,A2,B1,B2;
  int TypeB;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r2,ri2,rk2,temp,length;
  VECTOR dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;

  UHostCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

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
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(r2<CutOffBondDipoleBondDipoleSquared)
            UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA,dipoleB,dr,sqrt(r2));
        }
      }
    }
  }
  return 0;
}

void CalculateFrameworkShiftEnergyDifferenceAdsorbateVDW(void)
{
  int i,j,m,typeA,typeB;
  POINT posB,posA_new,posA_old;
  REAL rr;
  VECTOR dr;

  // Framework-Adsorbate energy
  UAdsorbateVDWDelta[CurrentSystem]=0.0;
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][i].AnisotropicPosition;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferenceAnisotropicPosition;

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Adsorbates[CurrentSystem][m].Atoms[j].Type;
        posB=Adsorbates[CurrentSystem][m].Atoms[j].AnisotropicPosition;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<CutOffVDWSquared)
          UAdsorbateVDWDelta[CurrentSystem]+=PotentialValue(typeA,typeB,rr,1.0);

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<CutOffVDWSquared)
          UAdsorbateVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,1.0);

      }
    }
  }
}

void CalculateFrameworkShiftEnergyDifferenceAdsorbateCharge(void)
{
  int i,j,m,typeA,typeB;
  int A1,A2,B1,B2,TypeMolB;
  POINT posB,posA_new,posA_old;
  REAL rr,r,chargeA,chargeB;
  VECTOR posA1,posA2,dipoleA_new,dipoleA_old;
  VECTOR posB1,posB2,dipoleB;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,temp;
  VECTOR dr;

  // Framework-Adsorbate energy
  UAdsorbateChargeChargeRealDelta[CurrentSystem]=0.0;
  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;


  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferencePosition;
    chargeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge;

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Adsorbates[CurrentSystem][m].Atoms[j].Type;
        posB=Adsorbates[CurrentSystem][m].Atoms[j].Position;
        chargeB=Adsorbates[CurrentSystem][m].Atoms[j].Charge;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeChargeSquared[CurrentSystem])
        {
          r=sqrt(rr);

          UAdsorbateChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA,chargeB,r);
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeChargeSquared[CurrentSystem])
        {
          r=sqrt(rr);

          UAdsorbateChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA,chargeB,r);
        }
      }

      TypeMolB=Adsorbates[CurrentSystem][m].Type;
      for(j=0;j<Components[TypeMolB].NumberOfBondDipoles;j++)
      {
        B1=Components[TypeMolB].BondDipoles[j].A;
        B2=Components[TypeMolB].BondDipoles[j].B;
        posB1=Adsorbates[CurrentSystem][m].Atoms[B1].Position;
        posB2=Adsorbates[CurrentSystem][m].Atoms[B2].Position;
        DipoleMagnitudeB=Components[TypeMolB].BondDipoleMagnitude[j];
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA,dipoleB,dr,sqrt(rr));

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA,dipoleB,dr,sqrt(rr));
      }
    }
  }

  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
    dipoleA_new.x=posA2.x-posA1.x;
    dipoleA_new.y=posA2.y-posA1.y;
    dipoleA_new.z=posA2.z-posA1.z;
    dipoleA_new=ApplyBoundaryCondition(dipoleA_new);
    posA_new.x=posA1.x+0.5*dipoleA_new.x;
    posA_new.y=posA1.y+0.5*dipoleA_new.y;
    posA_new.z=posA1.z+0.5*dipoleA_new.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z));
    dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
    dipoleA_old.x=posA2.x-posA1.x;
    dipoleA_old.y=posA2.y-posA1.y;
    dipoleA_old.z=posA2.z-posA1.z;
    dipoleA_old=ApplyBoundaryCondition(dipoleA_old);
    posA_old.x=posA1.x+0.5*dipoleA_old.x;
    posA_old.y=posA1.y+0.5*dipoleA_old.y;
    posA_old.z=posA1.z+0.5*dipoleA_old.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z));
    dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      TypeMolB=Adsorbates[CurrentSystem][m].Type;
      for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Adsorbates[CurrentSystem][m].Atoms[j].Type;
        posB=Adsorbates[CurrentSystem][m].Atoms[j].Position;
        chargeB=Adsorbates[CurrentSystem][m].Atoms[j].Charge;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
      }

      for(j=0;j<Components[TypeMolB].NumberOfBondDipoles;j++)
      {
        B1=Components[TypeMolB].BondDipoles[j].A;
        B2=Components[TypeMolB].BondDipoles[j].B;
        posB1=Adsorbates[CurrentSystem][m].Atoms[B1].Position;
        posB2=Adsorbates[CurrentSystem][m].Atoms[B2].Position;
        DipoleMagnitudeB=Components[TypeMolB].BondDipoleMagnitude[j];
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffBondDipoleBondDipoleSquared)
          UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(rr));

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffBondDipoleBondDipoleSquared)
          UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(rr));
      }
    }
  }
}

void CalculateFrameworkShiftEnergyDifferenceCationVDW(void)
{
  int i,j,m,typeA,typeB;
  POINT posB,posA_new,posA_old;
  REAL rr;
  VECTOR dr;

  // Framework-Cation energy
  UCationVDWDelta[CurrentSystem]=0.0;
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][i].AnisotropicPosition;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferenceAnisotropicPosition;

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Cations[CurrentSystem][m].Atoms[j].Type;
        posB=Cations[CurrentSystem][m].Atoms[j].AnisotropicPosition;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<CutOffVDWSquared)
          UCationVDWDelta[CurrentSystem]+=PotentialValue(typeA,typeB,rr,1.0);

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<CutOffVDWSquared)
          UCationVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,1.0);
      }
    }
  }
}


void CalculateFrameworkShiftEnergyDifferenceCationCharge(void)
{
  int i,j,m,typeA,typeB;
  int A1,A2,B1,B2,TypeMolB;
  POINT posB,posA_new,posA_old;
  REAL rr,r,chargeA,chargeB;
  VECTOR posA1,posA2,dipoleA_new,dipoleA_old;
  VECTOR posB1,posB2,dipoleB;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,temp;
  VECTOR dr;

  // Framework-Cation energy
  UCationChargeChargeRealDelta[CurrentSystem]=0.0;
  UCationChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;


  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferencePosition;
    chargeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Cations[CurrentSystem][m].Atoms[j].Type;
        posB=Cations[CurrentSystem][m].Atoms[j].Position;
        chargeB=Cations[CurrentSystem][m].Atoms[j].Charge;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeChargeSquared[CurrentSystem])
        {
          r=sqrt(rr);
          UCationChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA,chargeB,r);
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeChargeSquared[CurrentSystem])
        {
          r=sqrt(rr);

          UCationChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA,chargeB,r);
        }
      }

      TypeMolB=Cations[CurrentSystem][m].Type;
      for(j=0;j<Components[TypeMolB].NumberOfBondDipoles;j++)
      {
        B1=Components[TypeMolB].BondDipoles[j].A;
        B2=Components[TypeMolB].BondDipoles[j].B;
        posB1=Cations[CurrentSystem][m].Atoms[B1].Position;
        posB2=Cations[CurrentSystem][m].Atoms[B2].Position;
        DipoleMagnitudeB=Components[TypeMolB].BondDipoleMagnitude[j];
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
          UCationChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA,dipoleB,dr,sqrt(rr));

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
          UCationChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA,dipoleB,dr,sqrt(rr));
      }
    }
  }

  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
    dipoleA_new.x=posA2.x-posA1.x;
    dipoleA_new.y=posA2.y-posA1.y;
    dipoleA_new.z=posA2.z-posA1.z;
    dipoleA_new=ApplyBoundaryCondition(dipoleA_new);
    posA_new.x=posA1.x+0.5*dipoleA_new.x;
    posA_new.y=posA1.y+0.5*dipoleA_new.y;
    posA_new.z=posA1.z+0.5*dipoleA_new.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z));
    dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
    dipoleA_old.x=posA2.x-posA1.x;
    dipoleA_old.y=posA2.y-posA1.y;
    dipoleA_old.z=posA2.z-posA1.z;
    dipoleA_old=ApplyBoundaryCondition(dipoleA_old);
    posA_old.x=posA1.x+0.5*dipoleA_old.x;
    posA_old.y=posA1.y+0.5*dipoleA_old.y;
    posA_old.z=posA1.z+0.5*dipoleA_old.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z));
    dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      TypeMolB=Cations[CurrentSystem][m].Type;
      for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Cations[CurrentSystem][m].Atoms[j].Type;
        posB=Cations[CurrentSystem][m].Atoms[j].Position;
        chargeB=Cations[CurrentSystem][m].Atoms[j].Charge;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
          UCationChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
          UCationChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
      }

      for(j=0;j<Components[TypeMolB].NumberOfBondDipoles;j++)
      {
        B1=Components[TypeMolB].BondDipoles[j].A;
        B2=Components[TypeMolB].BondDipoles[j].B;
        posB1=Cations[CurrentSystem][m].Atoms[B1].Position;
        posB2=Cations[CurrentSystem][m].Atoms[B2].Position;
        DipoleMagnitudeB=Components[TypeMolB].BondDipoleMagnitude[j];
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffBondDipoleBondDipoleSquared)
          UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(rr));

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffBondDipoleBondDipoleSquared)
          UCationBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(rr));
      }
    }
  }
}


REAL CalculateEnergyDifferenceFrameworkMoveVDW(int atom_id,VECTOR posA,int typeA)
{
  int j,f1,typeB;
  REAL rr;
  VECTOR dr,posB;
  REAL UHostVDW;
  int ncell;
  int index_j;

  UHostVDW=0.0;
  if(!InternalFrameworkLennardJonesInteractions) return 0.0;

  // loop over all frameworks
  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
    {
      for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
      {
        if(j!=atom_id) // no self-contribution for a single atom
        {
          index_j=j+ncell*Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[CurrentFramework][atom_id][index_j],0))
          {
            posB=Framework[CurrentSystem].Atoms[CurrentFramework][j].AnisotropicPosition;
            typeB=Framework[CurrentSystem].Atoms[CurrentFramework][j].Type;

            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UHostVDW+=PotentialValue(typeA,typeB,rr,1.0);
          }
        }
      }
    }
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      if(f1!=CurrentFramework)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][j].AnisotropicPosition;
          typeB=Framework[CurrentSystem].Atoms[f1][j].Type;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UHostVDW+=PotentialValue(typeA,typeB,rr,1.0);
          }
        }
      }
    }
  }
  else
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
      {
        if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][atom_id][j],0))
        {
          posB=Framework[CurrentSystem].Atoms[f1][j].AnisotropicPosition;
          typeB=Framework[CurrentSystem].Atoms[f1][j].Type;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
            UHostVDW+=PotentialValue(typeA,typeB,rr,1.0);
        }
      }
    }
  }
  return UHostVDW;
}

void CalculateEnergyDifferenceFrameworkMoveCharge(int atom_id)
{
  int i,j,f1,typeA,typeB;
  int A1,A2,B1,B2;
  REAL rr,r;
  VECTOR dr,posB;
  VECTOR posA_new,posA_old;
  REAL chargeA,chargeB;
  VECTOR posA1,posA2,posB1,posB2,posB_new,posB_old;
  VECTOR dipoleB,dipoleA_new,dipoleA_old,dipoleB_new,dipoleB_old;
  REAL temp;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  int ncell;
  int index_j;

  UHostChargeChargeRealDelta[CurrentSystem]=0.0;
  UHostChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UHostBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Position;
  posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].ReferencePosition;
  typeA=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Type;
  chargeA=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Charge;

  // loop over all frameworks
  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
    {
      for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
      {
        if(j!=atom_id) // no self-contribution for a single atom
        {
          index_j=j+ncell*Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[CurrentFramework][atom_id][index_j],1))
          {
            posB=Framework[CurrentSystem].Atoms[CurrentFramework][j].Position;
            typeB=Framework[CurrentSystem].Atoms[CurrentFramework][j].Type;
            chargeB=Framework[CurrentSystem].Atoms[CurrentFramework][j].Charge;

            dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffChargeChargeSquared[CurrentSystem])
              UHostChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA,chargeB,r);

            dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffChargeChargeSquared[CurrentSystem])
              UHostChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA,chargeB,r);
          }
        }
      }
    }
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      if(f1!=CurrentFramework)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][j].AnisotropicPosition;
          typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffChargeChargeSquared[CurrentSystem])
              UHostChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA,chargeB,r);

            dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffChargeChargeSquared[CurrentSystem])
              UHostChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA,chargeB,r);
          }
        }
      }
    }
  }
  else
  {
    // loop over all frameworks
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
      {
        if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][atom_id][j],1))
        {
          posB=Framework[CurrentSystem].Atoms[f1][j].Position;
          typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
            UHostChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA,chargeB,r);

          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
            UHostChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA,chargeB,r);
        }
      }

      // compute the charge-bond-dipole interactions
      for(j=0;j<Framework[CurrentSystem].NumberOfBondDipoles[f1];j++)
      {
        if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][atom_id][j],2))
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][j];
          B1=Framework[CurrentSystem].BondDipoles[f1][j].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][j].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
            UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA,dipoleB,dr,sqrt(rr));

          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
            UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA,dipoleB,dr,sqrt(rr));
        }
      }
    }

    // loop over all the atoms of the currently selected framework
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
      A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
      A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

      // check if the current-atom is part of the current bond-dipole i
      if((atom_id==A1)||(atom_id==A2))
      {
        posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
        posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
        dipoleA_new.x=posA2.x-posA1.x;
        dipoleA_new.y=posA2.y-posA1.y;
        dipoleA_new.z=posA2.z-posA1.z;
        dipoleA_new=ApplyBoundaryCondition(dipoleA_new);
        posA_new.x=posA1.x+0.5*dipoleA_new.x;
        posA_new.y=posA1.y+0.5*dipoleA_new.y;
        posA_new.z=posA1.z+0.5*dipoleA_new.z;
        temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z));
        dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;

        posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
        posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
        dipoleA_old.x=posA2.x-posA1.x;
        dipoleA_old.y=posA2.y-posA1.y;
        dipoleA_old.z=posA2.z-posA1.z;
        dipoleA_old=ApplyBoundaryCondition(dipoleA_old);
        posA_old.x=posA1.x+0.5*dipoleA_old.x;
        posA_old.y=posA1.y+0.5*dipoleA_old.y;
        posA_old.z=posA1.z+0.5*dipoleA_old.z;
        temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z));
        dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;

        // loop over all frameworks
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
          {
            if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][j][i],2))
            {
              posB_new=Framework[CurrentSystem].Atoms[f1][j].Position;
              posB_old=Framework[CurrentSystem].Atoms[f1][j].ReferencePosition;
              typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
              chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

              dr.x=posA_new.x-posB_new.x;
              dr.y=posA_new.y-posB_new.y;
              dr.z=posA_new.z-posB_new.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeBondDipoleSquared)
                UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));

              dr.x=posA_old.x-posB_old.x;
              dr.y=posA_old.y-posB_old.y;
              dr.z=posA_old.z-posB_old.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeBondDipoleSquared)
                UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
            }
          }

          // compute the bond-dipole-bond-dipole interactions
          for(j=0;j<Framework[CurrentSystem].NumberOfBondDipoles[f1];j++)
          {
            if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],3))
            {
              DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][j];
              B1=Framework[CurrentSystem].BondDipoles[f1][j].A;
              B2=Framework[CurrentSystem].BondDipoles[f1][j].B;

              posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
              posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
              dipoleB_new.x=posB2.x-posB1.x;
              dipoleB_new.y=posB2.y-posB1.y;
              dipoleB_new.z=posB2.z-posB1.z;
              dipoleB_new=ApplyBoundaryCondition(dipoleB_new);
              posB_new.x=posB1.x+0.5*dipoleB_new.x;
              posB_new.y=posB1.y+0.5*dipoleB_new.y;
              posB_new.z=posB1.z+0.5*dipoleB_new.z;
              temp=DipoleMagnitudeB/sqrt(SQR(dipoleB_new.x)+SQR(dipoleB_new.y)+SQR(dipoleB_new.z));
              dipoleB_new.x*=temp; dipoleB_new.y*=temp; dipoleB_new.z*=temp;

              posB1=Framework[CurrentSystem].Atoms[f1][B1].ReferencePosition;
              posB2=Framework[CurrentSystem].Atoms[f1][B2].ReferencePosition;
              dipoleB_old.x=posB2.x-posB1.x;
              dipoleB_old.y=posB2.y-posB1.y;
              dipoleB_old.z=posB2.z-posB1.z;
              dipoleB_old=ApplyBoundaryCondition(dipoleB_old);
              posB_old.x=posB1.x+0.5*dipoleB_old.x;
              posB_old.y=posB1.y+0.5*dipoleB_old.y;
              posB_old.z=posB1.z+0.5*dipoleB_old.z;
              temp=DipoleMagnitudeB/sqrt(SQR(dipoleB_old.x)+SQR(dipoleB_old.y)+SQR(dipoleB_old.z));
              dipoleB_old.x*=temp; dipoleB_old.y*=temp; dipoleB_old.z*=temp;

              dr.x=posA_new.x-posB_new.x;
              dr.y=posA_new.y-posB_new.y;
              dr.z=posA_new.z-posB_new.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffBondDipoleBondDipoleSquared)
                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB_new,dr,sqrt(rr));

              dr.x=posA_old.x-posB_old.x;
              dr.y=posA_old.y-posB_old.y;
              dr.z=posA_old.z-posB_old.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffBondDipoleBondDipoleSquared)
                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB_old,dr,sqrt(rr));
            }
          }
        }
      }
    }
  }
}



void CalculateFrameworkEnergyDifferenceShiftedFramework(void)
{
  int i,j,f1,typeA,typeB;
  int A1,A2,B1,B2;
  REAL rr,r,chargeA,chargeB;
  VECTOR dr,posA_new,posA_old,posB,posA1,posA2,posB1,posB2;
  VECTOR dipoleB,dipoleA_new,dipoleA_old;
  REAL temp;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;

  // Framework-Adsorbate energy
  UHostVDWDelta[CurrentSystem]=0.0;
  UHostChargeChargeRealDelta[CurrentSystem]=0.0;
  UHostChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UHostBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  // loop over all the atoms of the currently selected framework
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferencePosition;
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    chargeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge;

    // loop over all frameworks
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      // except the current one (no self-interaction, simple translation of the whole framework)
      if(f1!=CurrentFramework)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][j].Position;
          typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
            if(InternalFrameworkLennardJonesInteractions)
              UHostVDWDelta[CurrentSystem]+=PotentialValue(typeA,typeB,rr,1.0);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
            UHostChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA,chargeB,r);

          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
            if(InternalFrameworkLennardJonesInteractions)
              UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,1.0);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
            UHostChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA,chargeB,r);
        }

        // compute the charge-bond-dipole interactions
        for(j=0;j<Framework[CurrentSystem].NumberOfBondDipoles[f1];j++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][j];
          B1=Framework[CurrentSystem].BondDipoles[f1][j].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][j].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
            UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA,dipoleB,dr,sqrt(rr));

          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
            UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA,dipoleB,dr,sqrt(rr));
        }
      }
    }
  }

  // loop over all the atoms of the currently selected framework
  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
    dipoleA_new.x=posA2.x-posA1.x;
    dipoleA_new.y=posA2.y-posA1.y;
    dipoleA_new.z=posA2.z-posA1.z;
    dipoleA_new=ApplyBoundaryCondition(dipoleA_new);
    posA_new.x=posA1.x+0.5*dipoleA_new.x;
    posA_new.y=posA1.y+0.5*dipoleA_new.y;
    posA_new.z=posA1.z+0.5*dipoleA_new.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z));
    dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
    dipoleA_old.x=posA2.x-posA1.x;
    dipoleA_old.y=posA2.y-posA1.y;
    dipoleA_old.z=posA2.z-posA1.z;
    dipoleA_old=ApplyBoundaryCondition(dipoleA_old);
    posA_old.x=posA1.x+0.5*dipoleA_old.x;
    posA_old.y=posA1.y+0.5*dipoleA_old.y;
    posA_old.z=posA1.z+0.5*dipoleA_old.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z));
    dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;


    // loop over all frameworks
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      // except the current one (no self-interaction, simple translation of the whole framework)
      if(f1!=CurrentFramework)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][j].Position;
          typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
            UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));

          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
            UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
        }

        // compute the bond-dipole-bond-dipole interactions
        for(j=0;j<Framework[CurrentSystem].NumberOfBondDipoles[f1];j++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][j];
          B1=Framework[CurrentSystem].BondDipoles[f1][j].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][j].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffBondDipoleBondDipoleSquared)
            UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(rr));

          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffBondDipoleBondDipoleSquared)
            UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(rr));
        }
      }
    }
  }
}

int CalculateFrameworkIntraReplicaVDWEnergy(void)
{
  int i,j,typeA,typeB,start;
  REAL energy,rr;
  VECTOR posA,posB,dr;
  int f1,f2,ncell,index_j;
  int A,B;
  REAL r;
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

        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          if((f1==f2)&&(ncell==0)) start=i+1;
          else start=0;

          for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
          {
            index_j=j+ncell*Framework[CurrentSystem].NumberOfAtoms[f2];
            if(f1!=f2?TRUE:(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][index_j],0)))
            {
              typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
              posB=Framework[CurrentSystem].Atoms[f2][j].AnisotropicPosition;

              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr,1.0);

                if(ncell==0)
                  UHostHostVDW[CurrentSystem]+=energy;
                else
                  UHostHostVDW[CurrentSystem]+=0.5*energy;
              }
            }
          }
        }
      }
    }
  }

  // contributions from interactions 1-4 torsions
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
          posA=Framework[CurrentSystem].Atoms[f1][A].AnisotropicPosition;

          B=Framework[CurrentSystem].Torsions[f1][i].D;
          typeB=Framework[CurrentSystem].Atoms[f1][B].Type;
          posB=Framework[CurrentSystem].Atoms[f1][B].AnisotropicPosition;


          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          UHostHostVDW[CurrentSystem]+=parms[6]*PotentialValue(typeA,typeB,rr,1.0);
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraReplicaChargeChargeEnergy(void)
{
  int i,j,typeA,typeB,start;
  REAL chargeA,chargeB;
  REAL r,rr,energy;
  VECTOR posA,posB,dr;
  int f1,f2,ncell,index_j;
  int A,B;
  REAL alpha;
  REAL *parms;

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

        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
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

              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared[CurrentSystem])
              {
                energy=PotentialValueCoulombic(chargeA,chargeB,sqrt(rr));

                if(ncell==0)
                  UHostHostChargeChargeReal[CurrentSystem]+=energy;
                else
                  UHostHostChargeChargeReal[CurrentSystem]+=0.5*energy;
              }
            }
          }
        }
      }
    }
  }


  // contributions from interactions 1-4 torsions
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
          posA=Framework[CurrentSystem].Atoms[f1][A].AnisotropicPosition;

          B=Framework[CurrentSystem].Torsions[f1][i].D;
          chargeB=Framework[CurrentSystem].Atoms[f1][B].Charge;
          posB=Framework[CurrentSystem].Atoms[f1][B].AnisotropicPosition;


          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

          UHostHostChargeChargeReal[CurrentSystem]+=parms[7]*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
        }
      }
    }
  }
  return 0;
}


int CalculateFrameworkAdsorbateReplicaVDWEnergy(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,scaling;
  VECTOR posA,posB,dr;
  int ncell;

  UHostAdsorbateVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
      scaling=Adsorbates[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
      {
        UHostAdsorbateVDW[CurrentSystem]+=InterpolateVDWGrid(typeA,posA);
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
                UHostAdsorbateVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr,scaling);
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationReplicaVDWEnergy(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,scaling;
  VECTOR posA,posB,dr;
  int ncell;

  UHostCationVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
      scaling=Cations[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
      {
        UHostCationVDW[CurrentSystem]+=InterpolateVDWGrid(typeA,posA);
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
                UHostCationVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr,scaling);
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateReplicaChargeChargeEnergy(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,r,chargeA,chargeB;
  REAL energy,scaling;
  VECTOR posA,posB,dr;
  int ncell;

  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
      scaling=Adsorbates[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
      chargeA=scaling*Adsorbates[CurrentSystem][i].Atoms[j].Charge;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
      {
        UHostAdsorbateChargeChargeReal[CurrentSystem]+=InterpolateCoulombGrid(typeA,posA);
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared[CurrentSystem])
              {
                energy=PotentialValueCoulombic(chargeA,chargeB,r);

                // energy
                UHostAdsorbateChargeChargeReal[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationReplicaChargeChargeEnergy(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,r,chargeA,chargeB;
  REAL energy,scaling;
  VECTOR posA,posB,dr;
  int ncell;

  UHostCationChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].Position;
      scaling=Cations[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
      chargeA=scaling*Cations[CurrentSystem][i].Atoms[j].Charge;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
      {
        UHostCationChargeChargeReal[CurrentSystem]+=InterpolateCoulombGrid(typeA,posA);
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared[CurrentSystem])
              {
                energy=PotentialValueCoulombic(chargeA,chargeB,r);

                // energy
                UHostCationChargeChargeReal[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

// Monte Carlo routines
// ====================


/*********************************************************************************************************
 * Name       | CalculateFrameworkAdsorbateVDWEnergyDifference                                           *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes framework-molecule energy differences for a single adsorbate molecule:          *
 *            | 1) New=TRUE, Old=FAlse, computes the VDW energy of the new molecule (positions in        *
 *            |    'TrialAnisotropicPosition' and the molecule is of type 'comp'                         *
 *            | 2) New=FALSE, Old=TRUE, computes the VDW energy of molecule 'm' of type 'comp'           *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in VDW energy of the molecule 'm' with    *
 *            |    the new positions in 'TrialAnisotropicPosition' and the old positions                 *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the VDW energy of the new positions stored in             *
 *            |      'TrialAnisotropicPosition'                                                          *
 *            | False: whether or not to compute the VDW energy of the old positions of molecule 'm'     *
 * Note       |                                                                                          *
 * Used       | int TranslationMoveAdsorbate(void)                                                       *
 *            | int RandomTranslationMoveAdsorbate(void)                                                 *
 *            | int RotationMoveAdsorbate(void)                                                          *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CBCFSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateFrameworkAdsorbateVDWEnergyDifference(int m,int comp,int New,int Old,int CanUseGrid)
{
  int i,j,k,nr_atoms,typeA,typeB,f1;
  POINT posA_new,posA_old,posB;
  REAL rr,energy,scaling_new,scaling_old;
  VECTOR dr,s;
  int icell0,icell;
  int ncell;

  // Framework-Adsorbate energy
  OVERLAP=FALSE;
  UHostVDWDelta[CurrentSystem]=0.0;

  if(New)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      // grid interpolation; transfer to the correct coordinates
      nr_atoms=Components[comp].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        typeA=Components[comp].Type[j];
        posA_new=TrialAnisotropicPosition[CurrentSystem][j];
        scaling_new=CFVDWScaling[j];

        if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA])&&CanUseGrid)
          UHostVDWDelta[CurrentSystem]+=InterpolateVDWGrid(typeA,posA_new);
        else if(UseCellLists[CurrentSystem])
        {
          // convert from xyz to abc
          s.x=InverseBox[CurrentSystem].ax*posA_new.x+InverseBox[CurrentSystem].bx*posA_new.y+InverseBox[CurrentSystem].cx*posA_new.z;
          s.y=InverseBox[CurrentSystem].ay*posA_new.x+InverseBox[CurrentSystem].by*posA_new.y+InverseBox[CurrentSystem].cy*posA_new.z;
          s.z=InverseBox[CurrentSystem].az*posA_new.x+InverseBox[CurrentSystem].bz*posA_new.y+InverseBox[CurrentSystem].cz*posA_new.z;

          // apply boundary condition
          s.x-=(REAL)NINT(s.x);
          s.y-=(REAL)NINT(s.y);
          s.z-=(REAL)NINT(s.z);

          // s between 0 and 1
          s.x+=0.5;
          s.y+=0.5;
          s.z+=0.5;

          // compute the corresponding cell-id
          icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
                 ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
                 ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            // loop over cells
            for(i=0;i<27;i++)
            {
              icell=CellListMap[CurrentSystem][icell0][i];

              k=Framework[CurrentSystem].CellListHead[f1][icell];

              while(k>=0)
              {
                typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
                posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                {
                  energy=PotentialValue(typeA,typeB,rr,1.0);
                  if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                  UHostVDWDelta[CurrentSystem]+=energy;
                }

                k=Framework[CurrentSystem].CellList[f1][k];
              }
            }
          }
        }
        else if(UseReplicas[CurrentSystem])
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                {
                  energy=PotentialValue(typeA,typeB,rr,scaling_new);
                  if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                  UHostVDWDelta[CurrentSystem]+=energy;
                }
              }
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr,scaling_new);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UHostVDWDelta[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }


  if(Old)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      // grid interpolation; transfer to the correct coordinates
      nr_atoms=Components[comp].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        typeA=Components[comp].Type[j];
        posA_old=Adsorbates[CurrentSystem][m].Atoms[j].AnisotropicPosition;
        scaling_old=Adsorbates[CurrentSystem][m].Atoms[j].CFVDWScalingParameter;

        if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA])&&CanUseGrid)
          UHostVDWDelta[CurrentSystem]-=InterpolateVDWGrid(typeA,posA_old);
        else if(UseCellLists[CurrentSystem])
        {
          // convert from xyz to abc
          s.x=InverseBox[CurrentSystem].ax*posA_old.x+InverseBox[CurrentSystem].bx*posA_old.y+InverseBox[CurrentSystem].cx*posA_old.z;
          s.y=InverseBox[CurrentSystem].ay*posA_old.x+InverseBox[CurrentSystem].by*posA_old.y+InverseBox[CurrentSystem].cy*posA_old.z;
          s.z=InverseBox[CurrentSystem].az*posA_old.x+InverseBox[CurrentSystem].bz*posA_old.y+InverseBox[CurrentSystem].cz*posA_old.z;

          // apply boundary condition
          s.x-=(REAL)NINT(s.x);
          s.y-=(REAL)NINT(s.y);
          s.z-=(REAL)NINT(s.z);

          // s between 0 and 1
          s.x+=0.5;
          s.y+=0.5;
          s.z+=0.5;

          // compute the corresponding cell-id
          icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
                 ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
                 ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            // loop over cells
            for(i=0;i<27;i++)
            {
              icell=CellListMap[CurrentSystem][icell0][i];

              k=Framework[CurrentSystem].CellListHead[f1][icell];

              while(k>=0)
              {
                typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
                posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                  UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,1.0);

                k=Framework[CurrentSystem].CellList[f1][k];
              }
            }
          }
        }
        else if(UseReplicas[CurrentSystem])
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                  UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scaling_old);
              }
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scaling_old);
            }
          }
        }
      }
    }
  }

  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateFrameworkCationVDWEnergyDifference                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes framework-molecule energy differences for a single cation molecule:             *
 *            | 1) New=TRUE, Old=FAlse, computes the VDW energy of the new molecule (positions in        *
 *            |    'TrialAnisotropicPosition' and the molecule is of type 'comp'                         *
 *            | 2) New=FALSE, Old=TRUE, computes the VDW energy of molecule 'm' of type 'comp'           *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in VDW energy of the molecule 'm' with    *
 *            |    the new positions in 'TrialAnisotropicPosition' and the old positions                 *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the VDW energy of the new positions stored in             *
 *            |      'TrialAnisotropicPosition'                                                          *
 *            | False: whether or not to compute the VDW energy of the old positions of molecule 'm'     *
 * Note       |                                                                                          *
 * Used       | int TranslationMoveCation(void)                                                          *
 *            | int RandomTranslationMoveCation(void)                                                    *
 *            | int RotationMoveCation(void)                                                             *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CBCFSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateFrameworkCationVDWEnergyDifference(int m,int comp,int New,int Old,int CanUseGrid)
{
  int i,j,k,nr_atoms,typeA,typeB,f1;
  POINT posA_new,posA_old,posB;
  REAL rr,energy,scaling_new,scaling_old;
  VECTOR dr,s;
  int icell0,icell;
  int ncell;

  // Framework-Cation energy
  OVERLAP=FALSE;
  UHostVDWDelta[CurrentSystem]=0.0;
  if(New)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      // grid interpolation; transfer to the correct coordinates
      nr_atoms=Components[comp].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        typeA=Components[comp].Type[j];
        posA_new=TrialAnisotropicPosition[CurrentSystem][j];
        scaling_new=CFVDWScaling[j];

        if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA])&&CanUseGrid)
          UHostVDWDelta[CurrentSystem]+=InterpolateVDWGrid(typeA,posA_new);
        else if(UseCellLists[CurrentSystem])
        {
          // convert from xyz to abc
          s.x=InverseBox[CurrentSystem].ax*posA_new.x+InverseBox[CurrentSystem].bx*posA_new.y+InverseBox[CurrentSystem].cx*posA_new.z;
          s.y=InverseBox[CurrentSystem].ay*posA_new.x+InverseBox[CurrentSystem].by*posA_new.y+InverseBox[CurrentSystem].cy*posA_new.z;
          s.z=InverseBox[CurrentSystem].az*posA_new.x+InverseBox[CurrentSystem].bz*posA_new.y+InverseBox[CurrentSystem].cz*posA_new.z;

          // apply boundary condition
          s.x-=(REAL)NINT(s.x);
          s.y-=(REAL)NINT(s.y);
          s.z-=(REAL)NINT(s.z);

          // s between 0 and 1
          s.x+=0.5;
          s.y+=0.5;
          s.z+=0.5;

          // compute the corresponding cell-id
          icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
                 ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
                 ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            // loop over cells
            for(i=0;i<27;i++)
            {
              icell=CellListMap[CurrentSystem][icell0][i];

              k=Framework[CurrentSystem].CellListHead[f1][icell];

              while(k>=0)
              {
                typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
                posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                {
                  energy=PotentialValue(typeA,typeB,rr,1.0);
                  if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                  UHostVDWDelta[CurrentSystem]+=energy;
                }

                k=Framework[CurrentSystem].CellList[f1][k];
              }
            }
          }
        }
        else if(UseReplicas[CurrentSystem])
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                {
                  energy=PotentialValue(typeA,typeB,rr,scaling_new);
                  if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                  UHostVDWDelta[CurrentSystem]+=energy;
                }
              }
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr,scaling_new);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UHostVDWDelta[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }

  if(Old)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      // grid interpolation; transfer to the correct coordinates
      nr_atoms=Components[comp].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        typeA=Components[comp].Type[j];
        posA_old=Cations[CurrentSystem][m].Atoms[j].AnisotropicPosition;
        scaling_old=Cations[CurrentSystem][m].Atoms[j].CFVDWScalingParameter;

        if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA])&&CanUseGrid)
          UHostVDWDelta[CurrentSystem]-=InterpolateVDWGrid(typeA,posA_old);
        else if(UseCellLists[CurrentSystem])
        {
          // convert from xyz to abc
          s.x=InverseBox[CurrentSystem].ax*posA_old.x+InverseBox[CurrentSystem].bx*posA_old.y+InverseBox[CurrentSystem].cx*posA_old.z;
          s.y=InverseBox[CurrentSystem].ay*posA_old.x+InverseBox[CurrentSystem].by*posA_old.y+InverseBox[CurrentSystem].cy*posA_old.z;
          s.z=InverseBox[CurrentSystem].az*posA_old.x+InverseBox[CurrentSystem].bz*posA_old.y+InverseBox[CurrentSystem].cz*posA_old.z;

          // apply boundary condition
          s.x-=(REAL)NINT(s.x);
          s.y-=(REAL)NINT(s.y);
          s.z-=(REAL)NINT(s.z);

          // s between 0 and 1
          s.x+=0.5;
          s.y+=0.5;
          s.z+=0.5;

          // compute the corresponding cell-id
          icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
                 ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
                 ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            // loop over cells
            for(i=0;i<27;i++)
            {
              icell=CellListMap[CurrentSystem][icell0][i];

              k=Framework[CurrentSystem].CellListHead[f1][icell];

              while(k>=0)
              {
                typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
                posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                  UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,1.0);

                k=Framework[CurrentSystem].CellList[f1][k];
              }
            }
          }
        }
        else if(UseReplicas[CurrentSystem])
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                  UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scaling_old);
              }
              
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scaling_old);
              
            }
          }
        }
      }
    }
  }

  return 0;
}

/*********************************************************************************************************
 * Name       | CalculateFrameworkAdsorbateChargeChargeEnergyDifference                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes framework-molecule energy differences for a single adsorbate molecule:          *
 *            | 1) New=TRUE, Old=FAlse, computes the charge-charge energy of the new molecule (positions *
 *            |    in 'TrialPosition' and the molecule is of type 'comp'                                 *
 *            | 2) New=FALSE, Old=TRUE, computes the charge-charge energy of molecule 'm' of type 'comp' *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in charge-charge energy of the molecule   *
 *            |    'm' with the new positions in 'TrialPosition' and the old positions                   *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the charge-charge energy of the new positions stored in   *
 *            |      'TrialPosition'                                                                     *
 *            | False: whether or not to compute the charge-charge energy of the old positions of        *
 *            |        molecule 'm'                                                                      *
 * Note       |                                                                                          *
 * Used       | int TranslationMoveAdsorbate(void)                                                       *
 *            | int RandomTranslationMoveAdsorbate(void)                                                 *
 *            | int RotationMoveAdsorbate(void)                                                          *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CBCFSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateFrameworkAdsorbateChargeChargeEnergyDifference(int m,int comp,int New,int Old,int CanUseGrid)
{
  int j,k,f1,typeA;
  POINT posA_new,posA_old,posB;
  REAL rr,chargeB,energy;
  REAL chargeA_old,chargeA_new;
  VECTOR dr;
  int ncell;

  // Framework-Adsorbate energy
  UHostChargeChargeRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  if(New)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      // grid interpolation; transfer to the correct coordinates
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
      {
        posA_new=TrialPosition[CurrentSystem][j];
        typeA=Components[comp].Type[j];
        chargeA_new=CFChargeScaling[j]*PseudoAtoms[typeA].Charge1;

        if((ChargeMethod!=NONE)&&(Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid)&&CanUseGrid)
          UHostChargeChargeRealDelta[CurrentSystem]+=InterpolateCoulombGrid(typeA,posA_new);
        else if(UseReplicas[CurrentSystem])
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeChargeSquared[CurrentSystem])
                {
                  energy=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(rr));
                  if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                  UHostChargeChargeRealDelta[CurrentSystem]+=energy;
                }
              }
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared[CurrentSystem])
              {
                energy=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(rr));
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UHostChargeChargeRealDelta[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }

  if(Old)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      // grid interpolation; transfer to the correct coordinates
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
      {
        posA_old=Adsorbates[CurrentSystem][m].Atoms[j].Position;
        typeA=Adsorbates[CurrentSystem][m].Atoms[j].Type;
        chargeA_old=Adsorbates[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Adsorbates[CurrentSystem][m].Atoms[j].Charge;
        
        if((ChargeMethod!=NONE)&&(Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid)&&CanUseGrid)
          UHostChargeChargeRealDelta[CurrentSystem]-=InterpolateCoulombGrid(typeA,posA_old);
        else if(UseReplicas[CurrentSystem])
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  UHostChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(rr));
              }
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared[CurrentSystem])
                UHostChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(rr));
            }
          }
        }
      }
    }
  }

  return 0;
}

/*********************************************************************************************************
 * Name       | CalculateFrameworkCationChargeChargeEnergyDifference                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes framework-molecule energy differences for a single cation molecule:             *
 *            | 1) New=TRUE, Old=FAlse, computes the charge-charge energy of the new molecule (positions *
 *            |    in 'TrialPosition' and the molecule is of type 'comp'                                 *
 *            | 2) New=FALSE, Old=TRUE, computes the charge-charge energy of molecule 'm' of type 'comp' *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in charge-charge energy of the molecule   *
 *            |    'm' with the new positions in 'TrialPosition' and the old positions                   *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the charge-charge energy of the new positions stored in   *
 *            |      'TrialPosition'                                                                     *
 *            | False: whether or not to compute the charge-charge energy of the old positions of        *
 *            |        molecule 'm'                                                                      *
 * Note       |                                                                                          *
 * Used       | int TranslationMoveAdsorbate(void)                                                       *
 *            | int RandomTranslationMoveAdsorbate(void)                                                 *
 *            | int RotationMoveAdsorbate(void)                                                          *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CBCFSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateFrameworkCationChargeChargeEnergyDifference(int m,int comp,int New,int Old,int CanUseGrid)
{
  int j,k,f1,typeA;
  POINT posA_new,posA_old,posB;
  REAL rr,chargeB;
  REAL chargeA_old,chargeA_new;
  VECTOR dr;
  int ncell;

  // Framework-Cation energy
  UHostChargeChargeRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  if(New)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      // grid interpolation; transfer to the correct coordinates
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
      {
        posA_new=TrialPosition[CurrentSystem][j];
        typeA=Components[comp].Type[j];
        chargeA_new=CFChargeScaling[j]*PseudoAtoms[typeA].Charge1;;
        
        if((ChargeMethod!=NONE)&&(Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid)&&CanUseGrid)
          UHostChargeChargeRealDelta[CurrentSystem]+=InterpolateCoulombGrid(typeA,posA_new);
        else if(UseReplicas[CurrentSystem])
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  UHostChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(rr));
              }
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared[CurrentSystem])
                UHostChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(rr));
            }
          }
        }
      }
    }
  }

  if(Old)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      // grid interpolation; transfer to the correct coordinates
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
      {
        posA_old=Cations[CurrentSystem][m].Atoms[j].Position;
        typeA=Cations[CurrentSystem][m].Atoms[j].Type;
        chargeA_old=Cations[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Cations[CurrentSystem][m].Atoms[j].Charge;
        
        if((ChargeMethod!=NONE)&&(Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid)&&CanUseGrid)
          UHostChargeChargeRealDelta[CurrentSystem]-=InterpolateCoulombGrid(typeA,posA_old);
        else if(UseReplicas[CurrentSystem])
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  UHostChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(rr));
              }
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared[CurrentSystem])
                UHostChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(rr));
              
            }
          }
        }
      }
    }
  }

  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes framework-molecule energy differences for a single adsorbate molecule:          *
 *            | 1) New=TRUE, Old=FAlse, computes the charge-bonddipole energy of the new molecule        *
 *            |    (positions in 'TrialPosition' and the molecule is of type 'comp'                      *
 *            | 2) New=FALSE, Old=TRUE, computes the charge-bonddipole energy of molecule 'm' of type    *
 *            |    'comp'                                                                                *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in charge-bonddipole energy of the        *
 *            |    molecule 'm' with the new positions in 'TrialPosition' and the old positions          *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the charge-bonddipole energy of the new positions stored  *
 *            |      in 'TrialPosition'                                                                  *
 *            | False: whether or not to compute the charge-bonddipole energy of the old positions of    *
 *            |        molecule 'm'                                                                      *
 * Note       |                                                                                          *
 * Used       | int TranslationMoveAdsorbate(void)                                                       *
 *            | int RandomTranslationMoveAdsorbate(void)                                                 *
 *            | int RotationMoveAdsorbate(void)                                                          *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CBCFSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference(int m,int comp,int New,int Old)
{
  int j,k,f1;
  int A1,A2,B1,B2;
  int Type;
  VECTOR posA_new,posA_old,posB,posA1,posA2,posB1,posB2,dr;
  REAL rr,ri2,rk2,temp,length,chargeB,chargeA_new,chargeA_old;
  VECTOR dipoleA_new,dipoleA_old,dipoleB;
  REAL scalingA1,scalingA2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  int ncell;

  UHostChargeBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if(UseReplicas[CurrentSystem])
    {
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
      {
        if(New)
        {
          posA_new=TrialPosition[CurrentSystem][j];
          chargeA_new=CFChargeScaling[j]*PseudoAtoms[Components[comp].Type[j]].Charge1;
        }

        if(Old)
        {
          chargeA_old=Adsorbates[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Adsorbates[CurrentSystem][m].Atoms[j].Charge;
          posA_old=Adsorbates[CurrentSystem][m].Atoms[j].Position;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
          {
            DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
            B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
            B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
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

            if(New)
            {
              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);;
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeBondDipoleSquared)
                  UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
              }
            }

            if(Old)
            {
              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeBondDipoleSquared)
                  UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
              }
            }
          }
        }
      }

      for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
      {
        A1=Components[comp].BondDipoles[j].A;
        A2=Components[comp].BondDipoles[j].B;
        DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];

        if(New)
        {
          posA1=TrialPosition[CurrentSystem][A1];
          posA2=TrialPosition[CurrentSystem][A2];
          scalingA1=CFChargeScaling[A1];
          scalingA2=CFChargeScaling[A2];
          dipoleA_new.x=posA2.x-posA1.x;
          dipoleA_new.y=posA2.y-posA1.y;
          dipoleA_new.z=posA2.z-posA1.z;
          posA_new.x=posA1.x+0.5*dipoleA_new.x;
          posA_new.y=posA1.y+0.5*dipoleA_new.y;
          posA_new.z=posA1.z+0.5*dipoleA_new.z;
          ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1*scalingA2)*DipoleMagnitudeA/length;
          dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
        }

        if(Old)
        {
          posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
          posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
          scalingA1=Adsorbates[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
          scalingA2=Adsorbates[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
          dipoleA_old.x=posA2.x-posA1.x;
          dipoleA_old.y=posA2.y-posA1.y;
          dipoleA_old.z=posA2.z-posA1.z;
          posA_old.x=posA1.x+0.5*dipoleA_old.x;
          posA_old.y=posA1.y+0.5*dipoleA_old.y;
          posA_old.z=posA1.z+0.5*dipoleA_old.z;
          ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            Type=Framework[CurrentSystem].Atoms[f1][k].Type;

            if(PseudoAtoms[Type].HasCharges)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              if(New)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                {
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffChargeBondDipoleSquared)
                    UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
                }
              }

              if(Old)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                {
                  dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffChargeBondDipoleSquared)
                    UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
                }
              }
            }
          }
        }
      }
    }
    else
    {
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
      {
        if(New)
        {
          posA_new=TrialPosition[CurrentSystem][j];
          chargeA_new=CFChargeScaling[j]*PseudoAtoms[Components[comp].Type[j]].Charge1;
        }

        if(Old)
        {
          chargeA_old=Adsorbates[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Adsorbates[CurrentSystem][m].Atoms[j].Charge;
          posA_old=Adsorbates[CurrentSystem][m].Atoms[j].Position;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
          {
            DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
            B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
            B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
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

            if(New)
            {
              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeBondDipoleSquared)
                UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
            }

            if(Old)
            {
              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeBondDipoleSquared)
                UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
            }
          }
        }
      }

      for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
      {
        A1=Components[comp].BondDipoles[j].A;
        A2=Components[comp].BondDipoles[j].B;
        DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];

        if(New)
        {
          posA1=TrialPosition[CurrentSystem][A1];
          posA2=TrialPosition[CurrentSystem][A2];
          scalingA1=CFChargeScaling[A1];
          scalingA2=CFChargeScaling[A2];
          dipoleA_new.x=posA2.x-posA1.x;
          dipoleA_new.y=posA2.y-posA1.y;
          dipoleA_new.z=posA2.z-posA1.z;
          posA_new.x=posA1.x+0.5*dipoleA_new.x;
          posA_new.y=posA1.y+0.5*dipoleA_new.y;
          posA_new.z=posA1.z+0.5*dipoleA_new.z;
          ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
        }

        if(Old)
        {
          posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
          posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
          scalingA1=Adsorbates[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
          scalingA2=Adsorbates[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
          dipoleA_old.x=posA2.x-posA1.x;
          dipoleA_old.y=posA2.y-posA1.y;
          dipoleA_old.z=posA2.z-posA1.z;
          posA_old.x=posA1.x+0.5*dipoleA_old.x;
          posA_old.y=posA1.y+0.5*dipoleA_old.y;
          posA_old.z=posA1.z+0.5*dipoleA_old.z;
          ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            Type=Framework[CurrentSystem].Atoms[f1][k].Type;

            if(PseudoAtoms[Type].HasCharges)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              if(New)
              {
                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeBondDipoleSquared)
                  UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
              }

              if(Old)
              {
                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeBondDipoleSquared)
                  UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
              }
            }
          }
        }
      }
    }
  }
  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateFrameworkCationChargeBondDipoleEnergyDifference                                 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes framework-molecule energy differences for a single cation molecule:             *
 *            | 1) New=TRUE, Old=FAlse, computes the charge-bonddipole energy of the new molecule        *
 *            |    (positions in 'TrialPosition' and the molecule is of type 'comp'                      *
 *            | 2) New=FALSE, Old=TRUE, computes the charge-bonddipole energy of molecule 'm' of type    *
 *            |    'comp'                                                                                *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in charge-bonddipole energy of the        *
 *            |    molecule 'm' with the new positions in 'TrialPosition' and the old positions          *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the charge-bonddipole energy of the new positions stored  *
 *            |      in 'TrialPosition'                                                                  *
 *            | False: whether or not to compute the charge-bonddipole energy of the old positions of    *
 *            |        molecule 'm'                                                                      *
 * Note       |                                                                                          *
 * Used       | int TranslationMoveCation(void)                                                          *
 *            | int RandomTranslationMoveCation(void)                                                    *
 *            | int RotationMoveCation(void)                                                             *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CBCFSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateFrameworkCationChargeBondDipoleEnergyDifference(int m,int comp,int New,int Old)
{
  int j,k,f1;
  int A1,A2,B1,B2;
  int Type;
  VECTOR posA_new,posA_old,posB,posA1,posA2,posB1,posB2,dr;
  REAL rr,ri2,rk2,temp,length,chargeB,chargeA_new,chargeA_old;
  VECTOR dipoleA_new,dipoleA_old,dipoleB;
  REAL scalingA1,scalingA2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  int ncell;

  UHostChargeBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if(UseReplicas[CurrentSystem])
    {
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
      {
        if(New)
        {
          posA_new=TrialPosition[CurrentSystem][j];
          chargeA_new=CFChargeScaling[j]*PseudoAtoms[Components[comp].Type[j]].Charge1;
        }

        if(Old)
        {
          chargeA_old=Cations[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Cations[CurrentSystem][m].Atoms[j].Charge;
          posA_old=Cations[CurrentSystem][m].Atoms[j].Position;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
          {
            DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
            B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
            B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
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

            if(New)
            {
              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);;
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeBondDipoleSquared)
                  UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
              }
            }

            if(Old)
            {
              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeBondDipoleSquared)
                  UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
              }
            }
          }
        }
      }

      for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
      {
        A1=Components[comp].BondDipoles[j].A;
        A2=Components[comp].BondDipoles[j].B;
        DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];

        if(New)
        {
          posA1=TrialPosition[CurrentSystem][A1];
          posA2=TrialPosition[CurrentSystem][A2];
          scalingA1=CFChargeScaling[A1];
          scalingA2=CFChargeScaling[A2];
          dipoleA_new.x=posA2.x-posA1.x;
          dipoleA_new.y=posA2.y-posA1.y;
          dipoleA_new.z=posA2.z-posA1.z;
          posA_new.x=posA1.x+0.5*dipoleA_new.x;
          posA_new.y=posA1.y+0.5*dipoleA_new.y;
          posA_new.z=posA1.z+0.5*dipoleA_new.z;
          ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1*scalingA2)*DipoleMagnitudeA/length;
          dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
        }

        if(Old)
        {
          posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
          posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
          scalingA1=Cations[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
          scalingA2=Cations[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
          dipoleA_old.x=posA2.x-posA1.x;
          dipoleA_old.y=posA2.y-posA1.y;
          dipoleA_old.z=posA2.z-posA1.z;
          posA_old.x=posA1.x+0.5*dipoleA_old.x;
          posA_old.y=posA1.y+0.5*dipoleA_old.y;
          posA_old.z=posA1.z+0.5*dipoleA_old.z;
          ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            Type=Framework[CurrentSystem].Atoms[f1][k].Type;

            if(PseudoAtoms[Type].HasCharges)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              if(New)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                {
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffChargeBondDipoleSquared)
                    UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
                }
              }

              if(Old)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                {
                  dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffChargeBondDipoleSquared)
                    UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
                }
              }
            }
          }
        }
      }
    }
    else
    {
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
      {
        if(New)
        {
          posA_new=TrialPosition[CurrentSystem][j];
          chargeA_new=CFChargeScaling[j]*PseudoAtoms[Components[comp].Type[j]].Charge1;
        }

        if(Old)
        {
          chargeA_old=Cations[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Cations[CurrentSystem][m].Atoms[j].Charge;
          posA_old=Cations[CurrentSystem][m].Atoms[j].Position;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
          {
            DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
            B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
            B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
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

            if(New)
            {
              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeBondDipoleSquared)
                UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
            }

            if(Old)
            {
              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeBondDipoleSquared)
                UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
            }
          }
        }
      }

      for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
      {
        A1=Components[comp].BondDipoles[j].A;
        A2=Components[comp].BondDipoles[j].B;
        DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];

        if(New)
        {
          posA1=TrialPosition[CurrentSystem][A1];
          posA2=TrialPosition[CurrentSystem][A2];
          scalingA1=CFChargeScaling[A1];
          scalingA2=CFChargeScaling[A2];
          dipoleA_new.x=posA2.x-posA1.x;
          dipoleA_new.y=posA2.y-posA1.y;
          dipoleA_new.z=posA2.z-posA1.z;
          posA_new.x=posA1.x+0.5*dipoleA_new.x;
          posA_new.y=posA1.y+0.5*dipoleA_new.y;
          posA_new.z=posA1.z+0.5*dipoleA_new.z;
          ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
        }

        if(Old)
        {
          posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
          posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
          scalingA1=Cations[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
          scalingA2=Cations[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
          dipoleA_old.x=posA2.x-posA1.x;
          dipoleA_old.y=posA2.y-posA1.y;
          dipoleA_old.z=posA2.z-posA1.z;
          posA_old.x=posA1.x+0.5*dipoleA_old.x;
          posA_old.y=posA1.y+0.5*dipoleA_old.y;
          posA_old.z=posA1.z+0.5*dipoleA_old.z;
          ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            Type=Framework[CurrentSystem].Atoms[f1][k].Type;

            if(PseudoAtoms[Type].HasCharges)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

              if(New)
              {
                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeBondDipoleSquared)
                  UHostChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
              }

              if(Old)
              {
                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeBondDipoleSquared)
                  UHostChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
              }
            }
          }
        }
      }
    }
  }
  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes framework-molecule energy differences for a single adsorbate molecule:          *
 *            | 1) New=TRUE, Old=FAlse, computes the bonddipole-bonddipole energy of the new molecule    *
 *            |    (positions in 'TrialPosition' and the molecule is of type 'comp'                      *
 *            | 2) New=FALSE, Old=TRUE, computes the bonddipole-bonddipole energy of molecule 'm' of     *
 *            |    type 'comp'                                                                           *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in bonddipole-bonddipole energy of the    *
 *            |    molecule 'm' with the new positions in 'TrialPosition' and the old positions          *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the bonddipole-bonddipole energy of the new positions     *
 *            |      stored in 'TrialPosition'                                                           *
 *            | False: whether or not to compute the bonddipole-bonddipole energy of the old positions   *
 *            |        of molecule 'm'                                                                   *
 * Note       |                                                                                          *
 * Used       | int TranslationMoveAdsorbate(void)                                                       *
 *            | int RandomTranslationMoveAdsorbate(void)                                                 *
 *            | int RotationMoveAdsorbate(void)                                                          *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CBCFSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference(int m,int comp,int New,int Old)
{
  int j,k,f1;
  int A1,A2,B1,B2;
  REAL r;
  VECTOR posA_new,posA_old,posB,posA1,posA2,posB1,posB2,dr;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL ri2,rk2,rr,length,temp;
  REAL scalingA1,scalingA2;
  VECTOR dipoleA_new,dipoleA_old,dipoleB;
  int ncell;

  OVERLAP=FALSE;
  UHostBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if(UseReplicas[CurrentSystem])
    {
      for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
      {
        A1=Components[comp].BondDipoles[j].A;
        A2=Components[comp].BondDipoles[j].B;
        DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];

        if(New)
        {
          posA1=TrialPosition[CurrentSystem][A1];
          posA2=TrialPosition[CurrentSystem][A2];
          scalingA1=CFChargeScaling[A1];
          scalingA2=CFChargeScaling[A2];
          dipoleA_new.x=posA2.x-posA1.x;
          dipoleA_new.y=posA2.y-posA1.y;
          dipoleA_new.z=posA2.z-posA1.z;
          posA_new.x=posA1.x+0.5*dipoleA_new.x;
          posA_new.y=posA1.y+0.5*dipoleA_new.y;
          posA_new.z=posA1.z+0.5*dipoleA_new.z;
          ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
        }

        if(Old)
        {
          posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
          posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
          scalingA1=Adsorbates[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
          scalingA2=Adsorbates[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
          dipoleA_old.x=posA2.x-posA1.x;
          dipoleA_old.y=posA2.y-posA1.y;
          dipoleA_old.z=posA2.z-posA1.z;
          posA_old.x=posA1.x+0.5*dipoleA_old.x;
          posA_old.y=posA1.y+0.5*dipoleA_old.y;
          posA_old.z=posA1.z+0.5*dipoleA_old.z;
          ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
          {
            DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
            B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
            B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
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

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffBondDipoleBondDipoleSquared)
                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));

              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffBondDipoleBondDipoleSquared)
                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
            }
          }
        }
      }
    }
    else
    {
      for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
      {
        A1=Components[comp].BondDipoles[j].A;
        A2=Components[comp].BondDipoles[j].B;
        DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];

        if(New)
        {
          posA1=TrialPosition[CurrentSystem][A1];
          posA2=TrialPosition[CurrentSystem][A2];
          scalingA1=CFChargeScaling[A1];
          scalingA2=CFChargeScaling[A2];
          dipoleA_new.x=posA2.x-posA1.x;
          dipoleA_new.y=posA2.y-posA1.y;
          dipoleA_new.z=posA2.z-posA1.z;
          posA_new.x=posA1.x+0.5*dipoleA_new.x;
          posA_new.y=posA1.y+0.5*dipoleA_new.y;
          posA_new.z=posA1.z+0.5*dipoleA_new.z;
          ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
        }

        if(Old)
        {
          posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
          posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
          scalingA1=Adsorbates[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
          scalingA2=Adsorbates[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
          dipoleA_old.x=posA2.x-posA1.x;
          dipoleA_old.y=posA2.y-posA1.y;
          dipoleA_old.z=posA2.z-posA1.z;
          posA_old.x=posA1.x+0.5*dipoleA_old.x;
          posA_old.y=posA1.y+0.5*dipoleA_old.y;
          posA_old.z=posA1.z+0.5*dipoleA_old.z;
          ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
          {
            DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
            B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
            B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
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

            if(New)
            {
              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffBondDipoleBondDipoleSquared)
                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));
            }

            if(Old)
            {
              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffBondDipoleBondDipoleSquared)
                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
            }
          }
        }
      }
    }
  }
  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes framework-molecule energy differences for a single cation molecule:             *
 *            | 1) New=TRUE, Old=FAlse, computes the bonddipole-bonddipole energy of the new molecule    *
 *            |    (positions in 'TrialPosition' and the molecule is of type 'comp'                      *
 *            | 2) New=FALSE, Old=TRUE, computes the bonddipole-bonddipole energy of molecule 'm' of     *
 *            |    type 'comp'                                                                           *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in bonddipole-bonddipole energy of the    *
 *            |    molecule 'm' with the new positions in 'TrialPosition' and the old positions          *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the bonddipole-bonddipole energy of the new positions     *
 *            |      stored in 'TrialPosition'                                                           *
 *            | False: whether or not to compute the bonddipole-bonddipole energy of the old positions   *
 *            |        of molecule 'm'                                                                   *
 * Note       |                                                                                          *
 * Used       | int TranslationMoveAdsorbate(void)                                                       *
 *            | int RandomTranslationMoveAdsorbate(void)                                                 *
 *            | int RotationMoveAdsorbate(void)                                                          *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CBCFSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference(int m,int comp,int New,int Old)
{
  int j,k,f1;
  int A1,A2,B1,B2;
  REAL r;
  VECTOR posA_new,posA_old,posB,posA1,posA2,posB1,posB2,dr;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL ri2,rk2,rr,length,temp;
  REAL scalingA1,scalingA2;
  VECTOR dipoleA_new,dipoleA_old,dipoleB;
  int ncell;

  OVERLAP=FALSE;
  UHostBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if(UseReplicas[CurrentSystem])
    {
      for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
      {
        A1=Components[comp].BondDipoles[j].A;
        A2=Components[comp].BondDipoles[j].B;
        DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];

        if(New)
        {
          posA1=TrialPosition[CurrentSystem][A1];
          posA2=TrialPosition[CurrentSystem][A2];
          scalingA1=CFChargeScaling[A1];
          scalingA2=CFChargeScaling[A2];
          dipoleA_new.x=posA2.x-posA1.x;
          dipoleA_new.y=posA2.y-posA1.y;
          dipoleA_new.z=posA2.z-posA1.z;
          posA_new.x=posA1.x+0.5*dipoleA_new.x;
          posA_new.y=posA1.y+0.5*dipoleA_new.y;
          posA_new.z=posA1.z+0.5*dipoleA_new.z;
          ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
        }

        if(Old)
        {
          posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
          posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
          scalingA1=Cations[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
          scalingA2=Cations[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
          dipoleA_old.x=posA2.x-posA1.x;
          dipoleA_old.y=posA2.y-posA1.y;
          dipoleA_old.z=posA2.z-posA1.z;
          posA_old.x=posA1.x+0.5*dipoleA_old.x;
          posA_old.y=posA1.y+0.5*dipoleA_old.y;
          posA_old.z=posA1.z+0.5*dipoleA_old.z;
          ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
          {
            DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
            B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
            B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
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

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffBondDipoleBondDipoleSquared)
                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));

              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffBondDipoleBondDipoleSquared)
                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
            }
          }
        }
      }
    }
    else
    {
      for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
      {
        A1=Components[comp].BondDipoles[j].A;
        A2=Components[comp].BondDipoles[j].B;
        DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];

        if(New)
        {
          posA1=TrialPosition[CurrentSystem][A1];
          posA2=TrialPosition[CurrentSystem][A2];
          scalingA1=CFChargeScaling[A1];
          scalingA2=CFChargeScaling[A2];
          dipoleA_new.x=posA2.x-posA1.x;
          dipoleA_new.y=posA2.y-posA1.y;
          dipoleA_new.z=posA2.z-posA1.z;
          posA_new.x=posA1.x+0.5*dipoleA_new.x;
          posA_new.y=posA1.y+0.5*dipoleA_new.y;
          posA_new.z=posA1.z+0.5*dipoleA_new.z;
          ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
        }

        if(Old)
        {
          posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
          posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
          scalingA1=Cations[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
          scalingA2=Cations[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
          dipoleA_old.x=posA2.x-posA1.x;
          dipoleA_old.y=posA2.y-posA1.y;
          dipoleA_old.z=posA2.z-posA1.z;
          posA_old.x=posA1.x+0.5*dipoleA_old.x;
          posA_old.y=posA1.y+0.5*dipoleA_old.y;
          posA_old.z=posA1.z+0.5*dipoleA_old.z;
          ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
          length=sqrt(ri2);
          temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
          dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
        }

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
          {
            DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
            B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
            B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
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

            if(New)
            {
              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffBondDipoleBondDipoleSquared)
                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));
            }

            if(Old)
            {
              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffBondDipoleBondDipoleSquared)
                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
            }
          }
        }
      }
    }
  }
  return 0;
}




//HERE
int CalculateFrameworkAdsorbateVDWEnergyDifferenceRXCM(int reaction,REAL Lambda1,REAL Lambda2,REAL LambdaNew,int **LambdaRetraceMolecules,int direction)
{
  int i,j,k,l,m,f1;
  int typeA,typeB;
  REAL rr,energy,scalingA_new,scalingA_old;
  VECTOR posA,posB,dr;
  int TypeMolA,TypeMolB;
  int ncell;
  int NumberOfRXMCMolecules;
  int RXMCMoleculeNumbers[256];
  int RXMCMoleculeType[256];

  OVERLAP=FALSE;
  UHostVDWDelta[CurrentSystem]=0.0;

  NumberOfRXMCMolecules=0;
  for(j=0;j<NumberOfComponents;j++)
  {
    for(k=0;k<ReactantsStoichiometry[reaction][j];k++)
    {
      RXMCMoleculeNumbers[NumberOfRXMCMolecules]=Components[j].ReactantFractionalMolecules[CurrentSystem][reaction][k];
      RXMCMoleculeType[NumberOfRXMCMolecules]=j;

      for(l=0;l<Components[j].NumberOfAtoms;l++)
        CFVDWScalingRXMC[NumberOfRXMCMolecules][l]=Lambda1;

      NumberOfRXMCMolecules++;
    }
  }
  for(j=0;j<NumberOfComponents;j++)
  {
    for(k=0;k<ProductsStoichiometry[reaction][j];k++)
    {
      RXMCMoleculeNumbers[NumberOfRXMCMolecules]=Components[j].ProductFractionalMolecules[CurrentSystem][reaction][k];
      RXMCMoleculeType[NumberOfRXMCMolecules]=j;

      for(l=0;l<Components[j].NumberOfAtoms;l++)
        CFVDWScalingRXMC[NumberOfRXMCMolecules][l]=Lambda2;
      NumberOfRXMCMolecules++;
    }
  }

  // add retrace molecules to list
  switch(direction)
  {
    case FORWARD:
      for(i=0;i<NumberOfComponents;i++)
      {
        for(k=0;k<ReactantsStoichiometry[reaction][i];k++)
        {
          RXMCMoleculeNumbers[NumberOfRXMCMolecules]=LambdaRetraceMolecules[i][k];
          RXMCMoleculeType[NumberOfRXMCMolecules]=Adsorbates[CurrentSystem][LambdaRetraceMolecules[i][k]].Type;

          for(l=0;l<Components[i].NumberOfAtoms;l++)
            CFVDWScalingRXMC[NumberOfRXMCMolecules][l]=1.0-LambdaNew;
          NumberOfRXMCMolecules++;
        }
      }
      break;
    case BACKWARD:
      for(i=0;i<NumberOfComponents;i++)
      {
        for(k=0;k<ProductsStoichiometry[reaction][i];k++)
        {
          RXMCMoleculeNumbers[NumberOfRXMCMolecules]=LambdaRetraceMolecules[i][k];
          RXMCMoleculeType[NumberOfRXMCMolecules]=Adsorbates[CurrentSystem][LambdaRetraceMolecules[i][k]].Type;

          for(l=0;l<Components[i].NumberOfAtoms;l++)
            CFVDWScalingRXMC[NumberOfRXMCMolecules][l]=1.0-LambdaNew;
          NumberOfRXMCMolecules++;
        }
      }
      break;
    case NO_FORWARD_OR_BACKWARD:
    default:
      break;
  }

  for(m=0;m<NumberOfRXMCMolecules;m++)
  {
    TypeMolA=RXMCMoleculeType[m];
    for(k=0;k<Components[TypeMolA].NumberOfAtoms;k++)
    {
      typeA=Components[TypeMolA].Type[k];

      posA=Adsorbates[CurrentSystem][RXMCMoleculeNumbers[m]].Atoms[k].AnisotropicPosition;
      scalingA_new=CFVDWScalingRXMC[m][k];
      scalingA_old=Adsorbates[CurrentSystem][RXMCMoleculeNumbers[m]].Atoms[k].CFVDWScalingParameter;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(l=0;l<Framework[CurrentSystem].NumberOfAtoms[f1];l++)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][l].Type;
          posB=Framework[CurrentSystem].Atoms[f1][l].AnisotropicPosition;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
          {
            energy=PotentialValue(typeA,typeB,rr,scalingA_new);
            if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
            UHostVDWDelta[CurrentSystem]+=energy;

            energy=PotentialValue(typeA,typeB,rr,scalingA_old);
            if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
            UHostVDWDelta[CurrentSystem]-=energy;
          }
        }
      }
    }
  }
}

int CalculateFrameworkAdsorbateVDWEnergyDifferenceNewRXCM(int reaction,REAL LambdaNew,int direction)
{
  int i,j,k,l,f1;
  int indexA,indexB;
  int typeA,typeB;
  int TypeMolA,TypeMolB;
  REAL rr,energy;
  VECTOR posA,posB,dr;
  REAL scalingA_new,scalingB;
  REAL scalingB_new;
  int NumberOfRXMCMolecules;
  int RXMCMoleculeNumbers[256];
  int RXMCMoleculeType[256];


  OVERLAP=FALSE;
  UHostVDWDelta[CurrentSystem]=0.0;

  NumberOfRXMCMolecules=0;
  switch(direction)
  {
    case FORWARD:
      for(j=0;j<NumberOfComponents;j++)
      {
        for(k=0;k<ProductsStoichiometry[reaction][j];k++)
        {
          RXMCMoleculeNumbers[NumberOfRXMCMolecules]=Components[j].ProductFractionalMolecules[CurrentSystem][reaction][k];
          RXMCMoleculeType[NumberOfRXMCMolecules]=j;

          for(l=0;l<Components[j].NumberOfAtoms;l++)
            CFVDWScalingRXMC[NumberOfRXMCMolecules][l]=LambdaNew;
          NumberOfRXMCMolecules++;
        }
      }
      break;
    case BACKWARD:
      for(j=0;j<NumberOfComponents;j++)
      {
        for(k=0;k<ReactantsStoichiometry[reaction][j];k++)
        {
          RXMCMoleculeNumbers[NumberOfRXMCMolecules]=Components[j].ReactantFractionalMolecules[CurrentSystem][reaction][k];
          RXMCMoleculeType[NumberOfRXMCMolecules]=j;

          for(l=0;l<Components[j].NumberOfAtoms;l++)
            CFVDWScalingRXMC[NumberOfRXMCMolecules][l]=LambdaNew;
          NumberOfRXMCMolecules++;
        }
      }
      break;
    default:
      break;
  }

  for(i=0;i<NumberOfRXMCMolecules;i++)
  {
    TypeMolA=RXMCMoleculeType[i];
    for(k=0;k<Components[TypeMolA].NumberOfAtoms;k++)
    {
      typeA=Components[TypeMolA].Type[k];
      posA=RXMCTrialAnisotropicPositions[CurrentSystem][i][k];
      scalingA_new=CFVDWScalingRXMC[i][k];

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(l=0;l<Framework[CurrentSystem].NumberOfAtoms[f1];l++)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][l].Type;
          posB=Framework[CurrentSystem].Atoms[f1][l].AnisotropicPosition;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
          {
            energy=PotentialValue(typeA,typeB,rr,scalingA_new);
            if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
            UHostVDWDelta[CurrentSystem]+=energy;
          }
        }
      }
    }
  }
}

int CalculateFrameworkAdsorbateChargeChargeEnergyDifferenceRXCM(int reaction,REAL Lambda1,REAL Lambda2,REAL LambdaNew,int **LambdaRetraceMolecules,int direction)
{
  int i,j,k,l,m,f1;
  int typeA,typeB;
  REAL rr,energy,scalingA_new,scalingA_old;
  VECTOR posA,posB,dr;
  int TypeMolA,TypeMolB;
  int ncell;
  int NumberOfRXMCMolecules;
  int RXMCMoleculeNumbers[256];
  int RXMCMoleculeType[256];
  REAL chargeA,chargeB,chargeA_new,chargeA_old;

  OVERLAP=FALSE;
  UHostChargeChargeRealDelta[CurrentSystem]=0.0;

  NumberOfRXMCMolecules=0;
  for(j=0;j<NumberOfComponents;j++)
  {
    for(k=0;k<ReactantsStoichiometry[reaction][j];k++)
    {
      RXMCMoleculeNumbers[NumberOfRXMCMolecules]=Components[j].ReactantFractionalMolecules[CurrentSystem][reaction][k];
      RXMCMoleculeType[NumberOfRXMCMolecules]=j;

      for(l=0;l<Components[j].NumberOfAtoms;l++)
        CFChargeScalingRXMC[NumberOfRXMCMolecules][l]=pow(Lambda1,5);

      NumberOfRXMCMolecules++;
    }
  }
  for(j=0;j<NumberOfComponents;j++)
  {
    for(k=0;k<ProductsStoichiometry[reaction][j];k++)
    {
      RXMCMoleculeNumbers[NumberOfRXMCMolecules]=Components[j].ProductFractionalMolecules[CurrentSystem][reaction][k];
      RXMCMoleculeType[NumberOfRXMCMolecules]=j;

      for(l=0;l<Components[j].NumberOfAtoms;l++)
        CFChargeScalingRXMC[NumberOfRXMCMolecules][l]=pow(Lambda2,5);
      NumberOfRXMCMolecules++;
    }
  }

  // add retrace molecules to list
  switch(direction)
  {
    case FORWARD:
      for(i=0;i<NumberOfComponents;i++)
      {
        for(k=0;k<ReactantsStoichiometry[reaction][i];k++)
        {
          RXMCMoleculeNumbers[NumberOfRXMCMolecules]=LambdaRetraceMolecules[i][k];
          RXMCMoleculeType[NumberOfRXMCMolecules]=Adsorbates[CurrentSystem][LambdaRetraceMolecules[i][k]].Type;

          for(l=0;l<Components[i].NumberOfAtoms;l++)
            CFChargeScalingRXMC[NumberOfRXMCMolecules][l]=pow(1.0-LambdaNew,5);
          NumberOfRXMCMolecules++;
        }
      }
      break;
    case BACKWARD:
      for(i=0;i<NumberOfComponents;i++)
      {
        for(k=0;k<ProductsStoichiometry[reaction][i];k++)
        {
          RXMCMoleculeNumbers[NumberOfRXMCMolecules]=LambdaRetraceMolecules[i][k];
          RXMCMoleculeType[NumberOfRXMCMolecules]=Adsorbates[CurrentSystem][LambdaRetraceMolecules[i][k]].Type;

          for(l=0;l<Components[i].NumberOfAtoms;l++)
            CFChargeScalingRXMC[NumberOfRXMCMolecules][l]=pow(1.0-LambdaNew,5);
          NumberOfRXMCMolecules++;
        }
      }
      break;
    case NO_FORWARD_OR_BACKWARD:
    default:
      break;
  }

  for(m=0;m<NumberOfRXMCMolecules;m++)
  {
    TypeMolA=RXMCMoleculeType[m];
    for(k=0;k<Components[TypeMolA].NumberOfAtoms;k++)
    {
      chargeA=Components[TypeMolA].Charge[k];

      posA=Adsorbates[CurrentSystem][RXMCMoleculeNumbers[m]].Atoms[k].AnisotropicPosition;
      scalingA_new=CFChargeScalingRXMC[m][k];
      scalingA_old=Adsorbates[CurrentSystem][RXMCMoleculeNumbers[m]].Atoms[k].CFChargeScalingParameter;
      chargeA_new=scalingA_new*chargeA;
      chargeA_old=scalingA_old*chargeA;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(l=0;l<Framework[CurrentSystem].NumberOfAtoms[f1];l++)
        {
          chargeB=Framework[CurrentSystem].Atoms[f1][l].Charge;
          posB=Framework[CurrentSystem].Atoms[f1][l].AnisotropicPosition;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
          {
            energy=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(rr));
            if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
            UHostChargeChargeRealDelta[CurrentSystem]+=energy;

            energy=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(rr));
            if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
            UHostChargeChargeRealDelta[CurrentSystem]-=energy;
          }
        }
      }
    }
  }
}

int CalculateFrameworkAdsorbateChargeChargeEnergyDifferenceNewRXCM(int reaction,REAL LambdaNew,int direction)
{
  int i,j,k,l,f1;
  int indexA,indexB;
  int typeA,typeB;
  int TypeMolA,TypeMolB;
  REAL rr,energy;
  VECTOR posA,posB,dr;
  REAL scalingA_new,scalingB;
  REAL scalingB_new;
  int NumberOfRXMCMolecules;
  int RXMCMoleculeNumbers[256];
  int RXMCMoleculeType[256];
  REAL chargeA_new,chargeB;

  OVERLAP=FALSE;
  UHostChargeChargeRealDelta[CurrentSystem]=0.0;

  NumberOfRXMCMolecules=0;
  switch(direction)
  {
    case FORWARD:
      for(j=0;j<NumberOfComponents;j++)
      {
        for(k=0;k<ProductsStoichiometry[reaction][j];k++)
        {
          RXMCMoleculeNumbers[NumberOfRXMCMolecules]=Components[j].ProductFractionalMolecules[CurrentSystem][reaction][k];
          RXMCMoleculeType[NumberOfRXMCMolecules]=j;

          for(l=0;l<Components[j].NumberOfAtoms;l++)
            CFChargeScalingRXMC[NumberOfRXMCMolecules][l]=pow(LambdaNew,5);
          NumberOfRXMCMolecules++;
        }
      }
      break;
    case BACKWARD:
      for(j=0;j<NumberOfComponents;j++)
      {
        for(k=0;k<ReactantsStoichiometry[reaction][j];k++)
        {
          RXMCMoleculeNumbers[NumberOfRXMCMolecules]=Components[j].ReactantFractionalMolecules[CurrentSystem][reaction][k];
          RXMCMoleculeType[NumberOfRXMCMolecules]=j;

          for(l=0;l<Components[j].NumberOfAtoms;l++)
            CFChargeScalingRXMC[NumberOfRXMCMolecules][l]=pow(LambdaNew,5);
          NumberOfRXMCMolecules++;
        }
      }
      break;
    default:
      break;
  }

  for(i=0;i<NumberOfRXMCMolecules;i++)
  {
    TypeMolA=RXMCMoleculeType[i];
    for(k=0;k<Components[TypeMolA].NumberOfAtoms;k++)
    {
      posA=RXMCTrialAnisotropicPositions[CurrentSystem][i][k];
      scalingA_new=CFChargeScalingRXMC[i][k];
      chargeA_new=scalingA_new*Components[TypeMolA].Charge[k];

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(l=0;l<Framework[CurrentSystem].NumberOfAtoms[f1];l++)
        {
          chargeB=Framework[CurrentSystem].Atoms[f1][l].Charge;
          posB=Framework[CurrentSystem].Atoms[f1][l].AnisotropicPosition;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
          {
            energy=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(rr));
            if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
            UHostChargeChargeRealDelta[CurrentSystem]+=energy;
          }
        }
      }
    }
  }
}


REAL CalculateFrameworkElectrostaticPotential(POINT posA)
{
  int i,typeB,f;
  REAL r,rr,chargeB;
  VECTOR posB,dr;
  REAL ElectrostaticPotential;

  ElectrostaticPotential=0.0;

  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if(ChargeMethod!=NONE)
    {
      for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f];i++)
        {
          typeB=Framework[CurrentSystem].Atoms[f][i].Type;
          posB=Framework[CurrentSystem].Atoms[f][i].Position;
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared[CurrentSystem])
          {
            r=sqrt(rr);
            chargeB=Framework[CurrentSystem].Atoms[f][i].Charge;

            ElectrostaticPotential+=PotentialValueCoulombic(1.0,chargeB,r);
          }
        }
      }
    }
  }
  return ElectrostaticPotential;
}

REAL CalculateFrameworkVDWEnergyCorrection(VECTOR* Positions,VECTOR *AnisotropicPositions,REAL *scaling)
{
  int k;
  int typeA;
  REAL UVDWDelta;
  VECTOR posA;

  UVDWDelta=0.0;
  for(k=0;k<Components[CurrentComponent].NumberOfAtoms;k++)
  {
    typeA=Components[CurrentComponent].Type[k];

    if(PseudoAtoms[typeA].AnisotropicCorrection)
    {
      posA=AnisotropicPositions[k];
      UVDWDelta+=CalculateFrameworkVDWEnergyAtPosition(posA,typeA,scaling[k]);

      posA=Positions[k];
      UVDWDelta-=CalculateFrameworkVDWEnergyAtPosition(posA,typeA,scaling[k]);
    }
  }
  return UVDWDelta;
}

