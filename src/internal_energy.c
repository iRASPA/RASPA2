/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'internal_energy.c' is part of RASPA-2.0

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "potentials.h"
#include "molecule.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "utils.h"
#include "ewald.h"
#include "internal_energy.h"
#include "minimization.h"

// Generate bondlength using the acceptence/rejection scheme
REAL GenerateBondlength(int i)
{
  REAL BondLength,r1,ran1,ran2,*parms,temp,temp2;
  REAL energy,exp_term;

  parms=Components[CurrentComponent].BondArguments[i];
  switch(Components[CurrentComponent].BondType[i])
  {
    case HARMONIC_BOND:
      // 0.5*p0*SQR(r-p1);
      // ===============================================
      // p_0/k_B [K/A^2]   force constant
      // p_1     [A]       reference bond distance
      ran1=1.0/sqrt(Beta[CurrentSystem]*parms[0]);
      ran2=1.0/SQR(parms[1]+3.0*ran1);
      do  BondLength=parms[1]+ran1*RandomGaussianNumber();
      while((RandomNumber()>(SQR(BondLength)*ran2))||(BondLength<=0.0));
      break;
    case CORE_SHELL_SPRING:
      ran1=1.0/sqrt(Beta[CurrentSystem]*parms[0]);
      ran2=1.0/SQR(parms[1]+3.0*ran1);
      do  BondLength=ran1*RandomGaussianNumber();
      while((RandomNumber()>(SQR(BondLength)*ran2))||(BondLength<=0.0));
      break;
    case MORSE_BOND:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference bond distance
      do
      {
        BondLength=RandomNumber()*3.0;
        energy=parms[0]*(SQR(1.0-exp(-parms[1]*(BondLength-parms[2])))-1.0);
      }while(RandomNumber()>SQR(BondLength)*exp(-Beta[CurrentSystem]*energy));
      break;
    case LJ_12_6_BOND:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      do
      {
        BondLength=RandomNumber()*3.0;
        temp=CUBE(1.0/SQR(BondLength));
        energy=parms[0]*SQR(temp)-parms[1]*temp;
      }while(RandomNumber()>SQR(BondLength)*exp(-Beta[CurrentSystem]*energy));
      break;
    case LENNARD_JONES_BOND:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A]
      do
      {
        BondLength=RandomNumber()*3.0;
        temp=CUBE(parms[1]/SQR(BondLength));
        energy=4.0*parms[0]*(temp*(temp-1.0));
      }while(RandomNumber()>SQR(BondLength)*exp(-Beta[CurrentSystem]*energy));
      break;
    case BUCKINGHAM_BOND:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      do
      {
        BondLength=RandomNumber()*3.0+0.8;
        temp=parms[2]*CUBE(1.0/SQR(BondLength));
        exp_term=parms[0]*exp(-parms[1]*BondLength);
        energy=exp_term-temp;
      }while(RandomNumber()>SQR(BondLength)*exp(-Beta[CurrentSystem]*energy));
      break;
    case RESTRAINED_HARMONIC_BOND:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      do
      {
        BondLength=RandomNumber()*3.0;
        r1=BondLength-parms[1];
        energy=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
                +parms[0]*parms[2]*MAX2(fabs(r1)-parms[2],(REAL)0.0);
      }while(RandomNumber()>SQR(BondLength)*exp(-Beta[CurrentSystem]*energy));
      break;
    case QUARTIC_BOND:
      // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
      // ===========================================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2/k_B [K/A^3]
      // p_3/k_B [K/A^4]
      do
      {
        BondLength=RandomNumber()*3.0;
        temp=BondLength-parms[1];
        temp2=SQR(temp);
        energy=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
      }while(RandomNumber()>SQR(BondLength)*exp(-Beta[CurrentSystem]*energy));
      break;
    case CFF_QUARTIC_BOND:
      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2/k_B [K/A^3]
      // p_3/k_B [K/A^4]
      do
      {
        BondLength=RandomNumber()*3.0;
        temp=BondLength-parms[1];
        temp2=SQR(temp);
        energy=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
      }while(RandomNumber()>SQR(BondLength)*exp(-Beta[CurrentSystem]*energy));

    case MM3_BOND:
      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
      // =================================================================
      // p_0     [mdyne/A molecule]
      // p_1     [A]
      do
      {
        BondLength=RandomNumber()*3.0;
        temp=BondLength-parms[1];
        temp2=SQR(BondLength-parms[1]);
        energy=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
      }while(RandomNumber()>SQR(BondLength)*exp(-Beta[CurrentSystem]*energy));
      break;
    case RIGID_BOND:
      BondLength=parms[0];
      fprintf(stderr, "Error !!! in generate RIGID_BOND\n");
      exit(0);
      break;
    case FIXED_BOND:
      BondLength=parms[0];
      break;
    default:
      fprintf(stderr, "Undefined Bond potential in routine 'GenerateBondlength' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return BondLength;
}

REAL GenerateBendAngle(int i)
{
  REAL theta,*parms;
  REAL energy;
  REAL temp,temp2;

  parms=Components[CurrentComponent].BendArguments[i];
  switch(Components[CurrentComponent].BendType[i])
  {
    case HARMONIC_BEND:
      // (1/2)p_0*(theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      do
      {
        theta=M_PI*RandomNumber();
        energy=0.5*parms[0]*SQR(theta-parms[1]);
      }while(RandomNumber()>SQR(sin(theta))*exp(-Beta[CurrentSystem]*energy));
      break;
    case QUARTIC_BEND:
      // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
      // ======================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      do
      {
        theta=M_PI*RandomNumber();
        energy=0.5*parms[0]*pow(theta-parms[1],2)+(1.0/3.0)*parms[2]*pow(theta-parms[1],3)+
               0.25*parms[3]*pow(theta-parms[1],4);
      }while(RandomNumber()>SQR(sin(theta))*exp(-Beta[CurrentSystem]*energy));
      break;
    case CFF_QUARTIC_BEND:
      // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
      // =====================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      do
      {
        theta=M_PI*RandomNumber();
        energy=parms[0]*pow(theta-parms[1],2)+parms[2]*pow(theta-parms[1],3)+
               parms[3]*pow(theta-parms[1],4);
      }while(RandomNumber()>SQR(sin(theta))*exp(-Beta[CurrentSystem]*energy));
      break;
    case HARMONIC_COSINE_BEND:
      // (1/2)*p_0*(cos(theta)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      do
      {
        theta=M_PI*RandomNumber();
        energy=0.5*parms[0]*SQR(cos(theta)-cos(parms[1]));
      }while(RandomNumber()>SQR(sin(theta))*exp(-Beta[CurrentSystem]*energy));
      break;
    case COSINE_BEND:
      // p_0*(1+cos(p_1*theta-p_2))
      // ===============================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      do
      {
        theta=M_PI*RandomNumber();
        energy=parms[0]*(1.0+cos(parms[1]*theta-parms[2]));
      }while(RandomNumber()>SQR(sin(theta))*exp(-Beta[CurrentSystem]*energy));
      break;
    case TAFIPOLSKY_BEND:
      // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
      // ===============================================
      // p_0/k_B [K]
      do
      {
        theta=M_PI*RandomNumber();
        energy=0.5*parms[0]*(1.0+cos(theta))*(1.0+cos(2.0*theta));
      }while(RandomNumber()>SQR(sin(theta))*exp(-Beta[CurrentSystem]*energy));
      break;
    case MM3_BEND:
    case MM3_IN_PLANE_BEND:
      // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      do
      {
        theta=M_PI*RandomNumber();
        temp=RAD2DEG*(theta-parms[1]);
        temp2=SQR(temp);
        energy=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
      }while(RandomNumber()>SQR(sin(theta))*exp(-Beta[CurrentSystem]*energy));
      break;
    case MEASURE_BEND:
      theta=109.5*RAD2DEG;
      break;
    case FIXED_BEND:
      theta=parms[0];
      break;
    default:
      fprintf(stderr, "Undefined Bend potential in routine 'GenerateBendAngle' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return theta;
}


REAL CalculateAngle(int A,int B,int C,int Iu)
{
  REAL theta,rab,rbc;
  VECTOR Rab,Rbc;

  Rab.x=TrialPositions[Iu][A].x-TrialPositions[Iu][B].x;
  Rab.y=TrialPositions[Iu][A].y-TrialPositions[Iu][B].y;
  Rab.z=TrialPositions[Iu][A].z-TrialPositions[Iu][B].z;
  rab=1.0/sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
  Rab.x*=rab;
  Rab.y*=rab;
  Rab.z*=rab;

  Rbc.x=TrialPositions[Iu][C].x-TrialPositions[Iu][B].x;
  Rbc.y=TrialPositions[Iu][C].y-TrialPositions[Iu][B].y;
  Rbc.z=TrialPositions[Iu][C].z-TrialPositions[Iu][B].z;
  rbc=1.0/sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x*=rbc;
  Rbc.y*=rbc;
  Rbc.z*=rbc;

  theta=acos(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
  return theta;
}

REAL CalculateAngle2(int A,int B,int group_nr,int Iu)
{
  int i;
  REAL theta,rab,rbc;
  VECTOR Rab,Rbc;
  VECTOR com;
  int atom_nr;

  com.x=com.y=com.z=0.0;
  for(i=0;i<Components[CurrentComponent].Groups[group_nr].NumberOfGroupAtoms;i++)
  {
    atom_nr=Components[CurrentComponent].Groups[group_nr].Atoms[i];
    com.x+=TrialPositions[Iu][atom_nr].x;
    com.y+=TrialPositions[Iu][atom_nr].y;
    com.z+=TrialPositions[Iu][atom_nr].z;
  }
  com.x/=Components[CurrentComponent].Groups[group_nr].NumberOfGroupAtoms;
  com.y/=Components[CurrentComponent].Groups[group_nr].NumberOfGroupAtoms;
  com.z/=Components[CurrentComponent].Groups[group_nr].NumberOfGroupAtoms;

  Rab.x=TrialPositions[Iu][A].x-TrialPositions[Iu][B].x;
  Rab.y=TrialPositions[Iu][A].y-TrialPositions[Iu][B].y;
  Rab.z=TrialPositions[Iu][A].z-TrialPositions[Iu][B].z;
  rab=1.0/sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
  Rab.x*=rab;
  Rab.y*=rab;
  Rab.z*=rab;

  Rbc.x=com.x-TrialPositions[Iu][B].x;
  Rbc.y=com.y-TrialPositions[Iu][B].y;
  Rbc.z=com.z-TrialPositions[Iu][B].z;
  rbc=1.0/sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x*=rbc;
  Rbc.y*=rbc;
  Rbc.z*=rbc;

  theta=acos(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
  return theta;
}


REAL CalculateBondEnergy(int Itype,int Iu)
{
  int A,B;
  REAL r,rr,r1,*parms,temp,temp2,U;
  VECTOR dr;

  A=Components[CurrentComponent].Bonds[Itype].A;
  B=Components[CurrentComponent].Bonds[Itype].B;

  dr.x=TrialPositions[Iu][A].x-TrialPositions[Iu][B].x;
  dr.y=TrialPositions[Iu][A].y-TrialPositions[Iu][B].y;
  dr.z=TrialPositions[Iu][A].z-TrialPositions[Iu][B].z;
  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  r=sqrt(rr);

  parms=Components[CurrentComponent].BondArguments[Itype];
  switch(Components[CurrentComponent].BondType[Itype])
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
      U=parms[0]*exp(-parms[1]*r)-parms[2]*CUBE(1.0/rr);
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
    case RIGID_BOND:
      U=0.0;
      break;
    case FIXED_BOND:
      U=0.0;
      break;
    case MEASURE_BOND:
      U=0.0;
      break;
    default:
      fprintf(stderr, "Undefined Bond potential in routine 'CalculateBondEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return U;
}

REAL CalculateBondEnergyAdsorbate(int m)
{
  int i;
  int NumberOfBonds,Type,A,B;
  VECTOR posA,posB,dr;
  REAL UBond,temp,temp2,r,rr,r1;
  REAL *parms;

  UBond=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBonds=Components[Type].NumberOfBonds;
  for(i=0;i<NumberOfBonds;i++)
  {
    A=Components[Type].Bonds[i].A;
    B=Components[Type].Bonds[i].B;

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    parms=Components[Type].BondArguments[i];

    switch(Components[Type].BondType[i])
    {
      case HARMONIC_BOND:
        // 0.5*p0*SQR(r-p1);
        // ===============================================
        // p_0/k_B [K/A^2]   force constant
        // p_1     [A]       reference bond distance
        UBond+=0.5*parms[0]*SQR(r-parms[1]);
        break;
      case CORE_SHELL_SPRING:
        UBond+=0.5*parms[0]*SQR(r);
        break;
      case MORSE_BOND:
        // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
        // ===============================================
        // p_0/k_B [K]       force constant
        // p_1     [A^-1]    parameter
        // p_2     [A]       reference bond distance
        temp=exp(parms[1]*(parms[2]-r));
        UBond+=parms[0]*(SQR(1.0-temp)-1.0);
        break;
      case LJ_12_6_BOND:
        // A/r_ij^12-B/r_ij^6
        // ===============================================
        // p_0/k_B [K A^12]
        // p_1/k_B [K A^6]
        temp=CUBE(1.0/rr);
        UBond+=parms[0]*SQR(temp)-parms[1]*temp;
        break;
      case LENNARD_JONES_BOND:
        // 4*p_0*((p_1/r)^12-(p_1/r)^6)
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A]
        temp=CUBE(parms[1]/rr);
        UBond+=4.0*parms[0]*(temp*(temp-1.0));
        break;
      case BUCKINGHAM_BOND:
        // p_0*exp(-p_1 r)-p_2/r^6
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A^-1]
        // p_2/k_B [K A^6]
        UBond+=parms[0]*exp(-parms[1]*r)-parms[2]*CUBE(1.0/rr);
        break;
      case RESTRAINED_HARMONIC_BOND:
        // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
        // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        r1=r-parms[1];
        UBond+=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
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
        UBond+=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
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
        UBond+=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
        break;
      case MM3_BOND:
        // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
        // =================================================================
        // p_0     [mdyne/A molecule]
        // p_1     [A]
        temp=r-parms[1];
        temp2=SQR(temp);
        UBond+=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
        break;
      case RIGID_BOND:
        UBond+=0.0;
        break;
      case FIXED_BOND:
        UBond+=0.0;
        break;
      case MEASURE_BOND:
        UBond+=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Bond potential in routine 'CalculateBondEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBond;
}

void CalculateBondEnergyAdsorbates(void)
{
  int i;

  UAdsorbateBond[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateBond[CurrentSystem]+=CalculateBondEnergyAdsorbate(i);
}

REAL CalculateBondEnergyCation(int m)
{
  int i;
  int NumberOfBonds,Type,A,B;
  VECTOR posA,posB,dr;
  REAL UBond,temp,temp2,r,rr,r1;
  REAL *parms;

  UBond=0.0;
  Type=Cations[CurrentSystem][m].Type;
  NumberOfBonds=Components[Type].NumberOfBonds;
  for(i=0;i<NumberOfBonds;i++)
  {
    A=Components[Type].Bonds[i].A;
    B=Components[Type].Bonds[i].B;

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    parms=Components[Type].BondArguments[i];

    switch(Components[Type].BondType[i])
    {
      case HARMONIC_BOND:
        // 0.5*p0*SQR(r-p1);
        // ===============================================
        // p_0/k_B [K/A^2]   force constant
        // p_1     [A]       reference bond distance
        UBond+=0.5*parms[0]*SQR(r-parms[1]);
        break;
      case CORE_SHELL_SPRING:
        UBond+=0.5*parms[0]*SQR(r);
        break;
      case MORSE_BOND:
        // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
        // ===============================================
        // p_0/k_B [K]       force constant
        // p_1     [A^-1]    parameter
        // p_2     [A]       reference bond distance
        temp=exp(parms[1]*(parms[2]-r));
        UBond+=parms[0]*(SQR(1.0-temp)-1.0);
        break;
      case LJ_12_6_BOND:
        // A/r_ij^12-B/r_ij^6
        // ===============================================
        // p_0/k_B [K A^12]
        // p_1/k_B [K A^6]
        temp=CUBE(1.0/rr);
        UBond+=parms[0]*SQR(temp)-parms[1]*temp;
        break;
      case LENNARD_JONES_BOND:
        // 4*p_0*((p_1/r)^12-(p_1/r)^6)
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A]
        temp=CUBE(parms[1]/rr);
        UBond+=4.0*parms[0]*(temp*(temp-1.0));
        break;
      case BUCKINGHAM_BOND:
        // p_0*exp(-p_1 r)-p_2/r^6
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A^-1]
        // p_2/k_B [K A^6]
        UBond+=parms[0]*exp(-parms[1]*r)-parms[2]*CUBE(1.0/rr);
        break;
      case RESTRAINED_HARMONIC_BOND:
        // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
        // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        r1=r-parms[1];
        UBond+=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
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
        UBond+=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
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
        UBond+=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
        break;
      case MM3_BOND:
        // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
        // =================================================================
        // p_0     [mdyne/A molecule]
        // p_1     [A]
        temp=r-parms[1];
        temp2=SQR(temp);
        UBond+=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
        break;
      case RIGID_BOND:
        UBond+=0.0;
        break;
      case FIXED_BOND:
        UBond+=0.0;
        break;
      case MEASURE_BOND:
        UBond+=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Bond potential in routine 'CalculateBondEnergyCation' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBond;
}

void CalculateBondEnergyCations(void)
{
  int i;

  UCationBond[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationBond[CurrentSystem]+=CalculateBondEnergyCation(i);
}


REAL CalculateUreyBradleyEnergy(int Itype,int Iu)
{
  int A,B;
  REAL r,rr,r1,*parms,energy,temp,temp2;
  VECTOR dr;

  A=Components[CurrentComponent].UreyBradleys[Itype].A;
  B=Components[CurrentComponent].UreyBradleys[Itype].C;
  dr.x=TrialPositions[Iu][A].x-TrialPositions[Iu][B].x;
  dr.y=TrialPositions[Iu][A].y-TrialPositions[Iu][B].y;
  dr.z=TrialPositions[Iu][A].z-TrialPositions[Iu][B].z;
  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  r=sqrt(rr);

  energy=0.0;
  parms=Components[CurrentComponent].UreyBradleyArguments[Itype];
  switch(Components[CurrentComponent].UreyBradleyType[Itype])
  {
    case HARMONIC_UREYBRADLEY:
      // 0.5*p0*SQR(r-p1);
      // ===============================================
      // p_0/k_B [K/A^2]   force constant
      // p_1     [A]       reference bond distance
      energy=0.5*parms[0]*SQR(r-parms[1]);
      break;
    case MORSE_UREYBRADLEY:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference bond distance
      temp=exp(parms[1]*(parms[2]-r));
      energy=parms[0]*(SQR(1.0-temp)-1.0);
      break;
    case LJ_12_6_UREYBRADLEY:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      temp=CUBE(1.0/rr);
      energy=parms[0]*SQR(temp)-parms[1]*temp;
      break;
    case LENNARD_JONES_UREYBRADLEY:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A]
      temp=CUBE(parms[1]/rr);
      energy=4.0*parms[0]*(temp*(temp-1.0));
      break;
    case BUCKINGHAM_UREYBRADLEY:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      energy=parms[0]*exp(-parms[1]*r)-parms[2]*CUBE(1.0/rr);
      break;
    case RESTRAINED_HARMONIC_UREYBRADLEY:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      r1=r-parms[1];
      energy=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
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
      energy=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
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
      energy=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
      break;
    case MM3_UREYBRADLEY:
      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
      // =================================================================
      // p_0     [mdyne/A molecule]
      // p_1     [A]
      temp=r-parms[1];
      temp2=SQR(temp);
      energy=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
      break;
    case RIGID_UREYBRADLEY:
      energy=0.0;
      break;
    case FIXED_UREYBRADLEY:
      energy=0.0;
      break;
    case MEASURE_UREYBRADLEY:
      energy=0.0;
      break;
    default:
      fprintf(stderr, "Undefined Urey-Bradley potential in routine 'CalculateUreyBradleyEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return energy;
}

REAL CalculateUreyBradleyEnergyAdsorbate(int m)
{
  int i;
  int NumberOfUreyBradleys,Type,A,B;
  VECTOR posA,posB,dr;
  REAL UUreyBradley,r,rr,r1,temp,temp2;
  REAL *parms;

  UUreyBradley=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfUreyBradleys=Components[Type].NumberOfUreyBradleys;
  for(i=0;i<NumberOfUreyBradleys;i++)
  {
    A=Components[Type].UreyBradleys[i].A;
    B=Components[Type].UreyBradleys[i].C;

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    parms=Components[Type].UreyBradleyArguments[i];

    switch(Components[Type].UreyBradleyType[i])
    {
      case HARMONIC_UREYBRADLEY:
        // 0.5*p0*SQR(r-p1);
        // ===============================================
        // p_0/k_B [K/A^2]   force constant
        // p_1     [A]       reference bond distance
        UUreyBradley+=0.5*parms[0]*SQR(r-parms[1]);
        break;
      case MORSE_UREYBRADLEY:
        // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
        // ===============================================
        // p_0/k_B [K]       force constant
        // p_1     [A^-1]    parameter
        // p_2     [A]       reference bond distance
        temp=exp(parms[1]*(parms[2]-r));
        UUreyBradley+=parms[0]*(SQR(1.0-temp)-1.0);
        break;
      case LJ_12_6_UREYBRADLEY:
        // A/r_ij^12-B/r_ij^6
        // ===============================================
        // p_0/k_B [K A^12]
        // p_1/k_B [K A^6]
        temp=CUBE(1.0/rr);
        UUreyBradley+=parms[0]*SQR(temp)-parms[1]*temp;
        break;
      case LENNARD_JONES_UREYBRADLEY:
        // 4*p_0*((p_1/r)^12-(p_1/r)^6)
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A]
        temp=CUBE(parms[1]/rr);
        UUreyBradley+=4.0*parms[0]*(temp*(temp-1.0));
        break;
      case BUCKINGHAM_UREYBRADLEY:
        // p_0*exp(-p_1 r)-p_2/r^6
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A^-1]
        // p_2/k_B [K A^6]
        UUreyBradley+=parms[0]*exp(-parms[1]*r)-parms[2]*CUBE(1.0/rr);
        break;
      case RESTRAINED_HARMONIC_UREYBRADLEY:
        // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
        // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        r1=r-parms[1];
        UUreyBradley+=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
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
        UUreyBradley+=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
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
        UUreyBradley+=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
        break;
      case MM3_UREYBRADLEY:
        // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
        // =================================================================
        // p_0     [mdyne/A molecule]
        // p_1     [A]
        temp=r-parms[1];
        temp2=SQR(temp);
        UUreyBradley+=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
        break;
      case RIGID_UREYBRADLEY:
        UUreyBradley+=0.0;
        break;
      case FIXED_UREYBRADLEY:
        UUreyBradley+=0.0;
        break;
      case MEASURE_UREYBRADLEY:
        UUreyBradley+=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Urey-Bradley potential in routine 'CalculateUreyBradleyEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UUreyBradley;
}

void CalculateUreyBradleyEnergyAdsorbates(void)
{
  int i;

  UAdsorbateUreyBradley[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateUreyBradley[CurrentSystem]+=CalculateUreyBradleyEnergyAdsorbate(i);
}

REAL CalculateUreyBradleyEnergyCation(int m)
{
  int i;
  int NumberOfUreyBradleys,Type,A,B;
  VECTOR posA,posB,dr;
  REAL UUreyBradley,r1,r,rr,temp,temp2;
  REAL *parms;

  UUreyBradley=0.0;
  Type=Cations[CurrentSystem][m].Type;
  NumberOfUreyBradleys=Components[Type].NumberOfUreyBradleys;
  for(i=0;i<NumberOfUreyBradleys;i++)
  {
    A=Components[Type].UreyBradleys[i].A;
    B=Components[Type].UreyBradleys[i].C;

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    parms=Components[Type].UreyBradleyArguments[i];

    switch(Components[Type].UreyBradleyType[i])
    {
      case HARMONIC_UREYBRADLEY:
        // 0.5*p0*SQR(r-p1);
        // ===============================================
        // p_0/k_B [K/A^2]   force constant
        // p_1     [A]       reference bond distance
        UUreyBradley+=0.5*parms[0]*SQR(r-parms[1]);
        break;
      case MORSE_UREYBRADLEY:
        // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
        // ===============================================
        // p_0/k_B [K]       force constant
        // p_1     [A^-1]    parameter
        // p_2     [A]       reference bond distance
        temp=exp(parms[1]*(parms[2]-r));
        UUreyBradley+=parms[0]*(SQR(1.0-temp)-1.0);
        break;
      case LJ_12_6_UREYBRADLEY:
        // A/r_ij^12-B/r_ij^6
        // ===============================================
        // p_0/k_B [K A^12]
        // p_1/k_B [K A^6]
        temp=CUBE(1.0/rr);
        UUreyBradley+=parms[0]*SQR(temp)-parms[1]*temp;
        break;
      case LENNARD_JONES_UREYBRADLEY:
        // 4*p_0*((p_1/r)^12-(p_1/r)^6)
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A]
        temp=CUBE(parms[1]/rr);
        UUreyBradley+=4.0*parms[0]*(temp*(temp-1.0));
        break;
      case BUCKINGHAM_UREYBRADLEY:
        // p_0*exp(-p_1 r)-p_2/r^6
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A^-1]
        // p_2/k_B [K A^6]
        UUreyBradley+=parms[0]*exp(-parms[1]*r)-parms[2]*CUBE(1.0/rr);
        break;
      case RESTRAINED_HARMONIC_UREYBRADLEY:
        // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
        // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        r1=r-parms[1];
        UUreyBradley+=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
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
        UUreyBradley+=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
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
        UUreyBradley+=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
        break;
      case MM3_UREYBRADLEY:
        // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
        // =================================================================
        // p_0     [mdyne/A molecule]
        // p_1     [A]
        temp=r-parms[1];
        temp2=SQR(temp);
        UUreyBradley+=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
        break;
      case RIGID_UREYBRADLEY:
        UUreyBradley+=0.0;
        break;
      case FIXED_UREYBRADLEY:
        UUreyBradley+=0.0;
        break;
      case MEASURE_UREYBRADLEY:
        UUreyBradley+=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Urey-Bradley potential in routine 'CalculateUreyBradleyEnergyCations' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UUreyBradley;
}

void CalculateUreyBradleyEnergyCations(void)
{
  int i;

  UCationUreyBradley[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationUreyBradley[CurrentSystem]+=CalculateUreyBradleyEnergyCation(i);
}

REAL CalculateBendEnergy(int Itype,int Iu)
{
  int A,B,C,D;
  REAL *parms,U,temp,temp2;
  REAL CosTheta,Theta;
  REAL rab,rbc,rac,energy;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL delta,rt2,rap2,rcp2;

  energy=0.0;
  A=Components[CurrentComponent].Bends[Itype].A;
  B=Components[CurrentComponent].Bends[Itype].B;
  C=Components[CurrentComponent].Bends[Itype].C;

  posA=TrialPositions[Iu][A];
  posB=TrialPositions[Iu][B];
  posC=TrialPositions[Iu][C];

  parms=Components[CurrentComponent].BendArguments[Itype];
  switch(Components[CurrentComponent].BendType[Itype])
  {
    case MM3_IN_PLANE_BEND:
      D=Components[CurrentComponent].Bends[Itype].D;
      posD=TrialPositions[Iu][D];
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

  switch(Components[CurrentComponent].BendType[Itype])
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
      fprintf(stderr, "Undefined Bend potential in routine 'CalculateBendEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return U;
}


REAL CalculateBendEnergyAdsorbate(int m)
{
  int i,A,B,C,D,Type,NumberOfBends;
  REAL *parms,U,temp,temp2;
  REAL CosTheta,Theta;
  REAL rab,rbc,rac,UBend;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL delta,rt2,rap2,rcp2;

  UBend=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBends=Components[Type].NumberOfBends;
  for(i=0;i<NumberOfBends;i++)
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

    switch(Components[Type].BendType[i])
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
        fprintf(stderr, "Undefined Bend potential in routine 'CalculateBendEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }

    // add contribution to the energy
    UBend+=U;

  }
  return UBend;
}

void CalculateBendEnergyAdsorbates(void)
{
  int i;

  UAdsorbateBend[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateBend[CurrentSystem]+=CalculateBendEnergyAdsorbate(i);
}

REAL CalculateBendEnergyCation(int m)
{
  int i,A,B,C,D,Type,NumberOfBends;
  REAL *parms,U,temp,temp2;
  REAL CosTheta,Theta;
  REAL rab,rbc,rac,UBend;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL delta,rt2,rap2,rcp2;

  UBend=0.0;
  Type=Cations[CurrentSystem][m].Type;
  NumberOfBends=Components[Type].NumberOfBends;
  for(i=0;i<NumberOfBends;i++)
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

    switch(Components[Type].BendType[i])
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
        fprintf(stderr, "Undefined Bend potential in routine 'CalculateBendEnergyCation' ('internal_energy.c')\n");
        exit(0);
        break;
    }

    // add contribution to the energy
    UBend+=U;
  }
  return UBend;
}

void CalculateBendEnergyCations(void)
{
  int i;

  UCationBend[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationBend[CurrentSystem]+=CalculateBendEnergyCation(i);
}

REAL CalculateInversionBendEnergy(int Itype,int Iu)
{
  int A,B,C,D;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2;
  REAL CosChi,Chi,energy;
  REAL temp,temp2;
  VECTOR Rab,Rbc,Rbd,Rcd,Rad;
  POINT posA,posB,posC,posD;

  A=Components[CurrentComponent].InversionBends[Itype].A;
  B=Components[CurrentComponent].InversionBends[Itype].B;
  C=Components[CurrentComponent].InversionBends[Itype].C;
  D=Components[CurrentComponent].InversionBends[Itype].D;

  posA=TrialPositions[Iu][A];
  posB=TrialPositions[Iu][B];
  posC=TrialPositions[Iu][C];
  posD=TrialPositions[Iu][D];

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
  Rad=ApplyBoundaryCondition(Rad);

  parms=Components[CurrentComponent].InversionBendArguments[Itype];

  switch(Components[CurrentComponent].InversionBendType[Itype])
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
      fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateInversionBendEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }

  e=Rab.x*(Rbd.y*Rbc.z-Rbd.z*Rbc.y)+Rab.y*(Rbd.z*Rbc.x-Rbd.x*Rbc.z)+Rab.z*(Rbd.x*Rbc.y-Rbd.y*Rbc.x);
  CosChi=sqrt((Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z)-SQR(e)/c)/rrab;

  // Ensure CosChi is between -1 and 1.
  CosChi=SIGN(MIN2(fabs(CosChi),(REAL)1.0),CosChi);


  switch(Components[CurrentComponent].InversionBendType[Itype])
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
      fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateInversionBendEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return energy;
}

REAL CalculateInversionBendEnergyAdsorbate(int m)
{
  int i,A,B,C,D;
  int Type,NumberOfInversionBends;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2;
  REAL CosChi,Chi;
  REAL temp,temp2;
  VECTOR Rab,Rbc,Rbd,Rcd,Rad;
  POINT posA,posB,posC,posD;
  REAL UAdsorbateInversionBend;


  UAdsorbateInversionBend=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfInversionBends=Components[Type].NumberOfInversionBends;
  for(i=0;i<NumberOfInversionBends;i++)
  {
    A=Components[Type].InversionBends[i].A;
    B=Components[Type].InversionBends[i].B;
    C=Components[Type].InversionBends[i].C;
    D=Components[Type].InversionBends[i].D;
    parms=Components[Type].InversionBendArguments[i];

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

    switch(Components[Type].InversionBendType[i])
    {
      case HARMONIC_INVERSION:
      case HARMONIC_INVERSION2:
        // (1/2)*p_0*(chi-p_1)^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        Chi=acos(CosChi);
        UAdsorbateInversionBend+=0.5*parms[0]*SQR(Chi-parms[1]);
        break;
      case HARMONIC_COSINE_INVERSION:
      case HARMONIC_COSINE_INVERSION2:
        // (1/2)*p_0*(cos(phi)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        UAdsorbateInversionBend+=0.5*parms[0]*SQR(CosChi-parms[1]);
        break;
      case PLANAR_INVERSION:
      case PLANAR_INVERSION2:
        // (1/2)*p_0*(1-cos(phi))
        // ===============================================
        // p_0/k_B [K]
        UAdsorbateInversionBend+=parms[0]*(1.0-CosChi);
        break;
      case MM3_INVERSION:
        // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
        // =================================================================================================
        // p_0/k_B [mdyne A/rad^2]
        // p_1     [degrees]
        Chi=acos(CosChi);
        temp=RAD2DEG*(Chi-parms[1]);
        temp2=SQR(temp);
        UAdsorbateInversionBend+=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
        break;
      default:
        fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateInversionBendEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }

  return UAdsorbateInversionBend;
}

void CalculateInversionBendEnergyAdsorbates(void)
{
  int i;

  UAdsorbateInversionBend[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateInversionBend[CurrentSystem]+=CalculateInversionBendEnergyAdsorbate(i);
}

REAL CalculateInversionBendEnergyCation(int m)
{
  int i,A,B,C,D;
  int Type,NumberOfInversionBends;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2;
  REAL CosChi,Chi;
  REAL temp,temp2;
  VECTOR Rab,Rbc,Rbd,Rcd,Rad;
  POINT posA,posB,posC,posD;
  REAL UCationInversionBend;


  UCationInversionBend=0.0;
  Type=Cations[CurrentSystem][m].Type;
  NumberOfInversionBends=Components[Type].NumberOfInversionBends;
  for(i=0;i<NumberOfInversionBends;i++)
  {
    A=Components[Type].InversionBends[i].A;
    B=Components[Type].InversionBends[i].B;
    C=Components[Type].InversionBends[i].C;
    D=Components[Type].InversionBends[i].D;
    parms=Components[Type].InversionBendArguments[i];

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

    switch(Components[Type].InversionBendType[i])
    {
      case HARMONIC_INVERSION:
      case HARMONIC_INVERSION2:
        // (1/2)*p_0*(chi-p_1)^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        Chi=acos(CosChi);
        UCationInversionBend+=0.5*parms[0]*SQR(Chi-parms[1]);
        break;
      case HARMONIC_COSINE_INVERSION:
      case HARMONIC_COSINE_INVERSION2:
        // (1/2)*p_0*(cos(phi)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        UCationInversionBend+=0.5*parms[0]*SQR(CosChi-parms[1]);
        break;
      case PLANAR_INVERSION:
      case PLANAR_INVERSION2:
        // (1/2)*p_0*(1-cos(phi))
        // ===============================================
        // p_0/k_B [K]
        UCationInversionBend+=parms[0]*(1.0-CosChi);
        break;
      case MM3_INVERSION:
        // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
        // =================================================================================================
        // p_0/k_B [mdyne A/rad^2]
        // p_1     [degrees]
        Chi=acos(CosChi);
        temp=RAD2DEG*(Chi-parms[1]);
        temp2=SQR(temp);
        UCationInversionBend+=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
        break;
      default:
        fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateInversionBendEnergyCation' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }

  return UCationInversionBend;
}

void CalculateInversionBendEnergyCations(void)
{
  int i;

  UCationInversionBend[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationInversionBend[CurrentSystem]+=CalculateInversionBendEnergyCation(i);
}

REAL ComputeTorsionAngle(VECTOR posA,VECTOR posB,VECTOR posC, VECTOR posD)
{
  REAL rbc;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi;
  VECTOR Pb,Pc;

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

  if(Phi<0.0) Phi+=2.0*M_PI;
  return Phi;
}

REAL CalculateTorsionEnergy(int Itype,int Iu)
{
  int A,B,C,D;
  REAL rbc,energy;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2;
  REAL ShiftedCosPhi,ShiftedCosPhi2;
  VECTOR Pb,Pc;
  REAL *parms;

  A=Components[CurrentComponent].Torsions[Itype].A;
  B=Components[CurrentComponent].Torsions[Itype].B;
  C=Components[CurrentComponent].Torsions[Itype].C;
  D=Components[CurrentComponent].Torsions[Itype].D;
  parms=Components[CurrentComponent].TorsionArguments[Itype];

  Dab.x=TrialPositions[Iu][A].x-TrialPositions[Iu][B].x;
  Dab.y=TrialPositions[Iu][A].y-TrialPositions[Iu][B].y;
  Dab.z=TrialPositions[Iu][A].z-TrialPositions[Iu][B].z;

  Dbc.x=TrialPositions[Iu][C].x-TrialPositions[Iu][B].x;
  Dbc.y=TrialPositions[Iu][C].y-TrialPositions[Iu][B].y;
  Dbc.z=TrialPositions[Iu][C].z-TrialPositions[Iu][B].z;
  rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
  Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

  Dcd.x=TrialPositions[Iu][D].x-TrialPositions[Iu][C].x;
  Dcd.y=TrialPositions[Iu][D].y-TrialPositions[Iu][C].y;
  Dcd.z=TrialPositions[Iu][D].z-TrialPositions[Iu][C].z;

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


  // Ensure CosPhi is between -1 and 1.
  CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
  CosPhi2=SQR(CosPhi);

  // Note that we can rewrite various terms into CosPhi-contributions using
  // cos(2*Theta)=2*cos(Theta)**2-1
  // cos(3*Theta)=-3*cos(Theta)*(1-SQR(cos(Theta))+CUBE(cos(Theta))
  // ..........

  switch(Components[CurrentComponent].TorsionType[Itype])
  {
    case HARMONIC_DIHEDRAL:
      // (1/2)*p_0*(phi-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // potential defined in terms of 'phi' and therefore contains a singularity
      // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
      // same direction as Rbc, and negative otherwise
      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);
      Phi-=parms[1];
      Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
      energy=0.5*parms[0]*SQR(Phi);
      break;
    case HARMONIC_COSINE_DIHEDRAL:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      energy=0.5*parms[0]*SQR(CosPhi-parms[1]);
      break;
    case THREE_COSINE_DIHEDRAL:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      energy=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
      break;
    case MM3_DIHEDRAL:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0     [kcal/mol]
      // p_1     [kcal/mol]
      // p_2     [kcal/mol]
      energy=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
      break;
    case CVFF_BLOCKED_DIHEDRAL:
      //
      // ========================================================================
      // p_0     [rad]
      // p_1     [K]
      // p_2     [-]
      // p_3     [rad]
      // p_4     [rad]
      energy=0.0;
      break;
    case CFF_DIHEDRAL:
      // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      energy=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
      break;
    case CFF_DIHEDRAL2:
      // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      energy=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
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
      // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
      // polymer convention Phi'(trans)=0.
      energy=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
                parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
      break;
    case TRAPPE_DIHEDRAL:
      // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
      // ==========================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      energy=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
      break;
   case TRAPPE_DIHEDRAL_EXTENDED:
      // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
      // =============================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      energy=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
      break;
   case MOD_TRAPPE_DIHEDRAL:
     /* Salvador modification: 16/08/2016
      add phase in cos function:
      p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
     */
      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);
      Phi-=parms[4];           // shift Phi as Phi+parms[4]
      Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
      ShiftedCosPhi=cos(Phi);
      ShiftedCosPhi2=SQR(ShiftedCosPhi);
      energy=parms[0]+parms[1]+parms[3]+(parms[1]-3.0*parms[3])*ShiftedCosPhi -2.0*parms[2]*ShiftedCosPhi2 +4.0*parms[3]*ShiftedCosPhi*ShiftedCosPhi2;
      break;
    case CVFF_DIHEDRAL:
      // p_0*(1+cos(p_1*phi-p_2))
      // ========================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      // potential defined in terms of 'phi' and therefore contains a singularity
      // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
      // same direction as Rbc, and negative otherwise
      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);
      energy=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
      break;
    case OPLS_DIHEDRAL:
      // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
      // =================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      energy=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
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
      energy=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
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
      energy=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
        2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
        2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
      break;
    case FIXED_DIHEDRAL:
      energy=0.0;
      break;
    default:
      fprintf(stderr, "Undefined Torsion potential in routine 'CalculateTorsionEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return energy;
}

REAL CalculateTorsionEnergyAdsorbate(int m)
{
  int i,Type,NumberOfTorsions,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL rbc,UTorsion;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2;
  REAL ShiftedCosPhi,ShiftedCosPhi2;
  VECTOR Pb,Pc;
  REAL *parms;

  UTorsion=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfTorsions=Components[Type].NumberOfTorsions;
  for(i=0;i<NumberOfTorsions;i++)
  {
    A=Components[Type].Torsions[i].A;
    B=Components[Type].Torsions[i].B;
    C=Components[Type].Torsions[i].C;
    D=Components[Type].Torsions[i].D;
    parms=(REAL*)&Components[Type].TorsionArguments[i];

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

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].TorsionType[i])
    {
      case HARMONIC_DIHEDRAL:
        // (1/2)*p_0*(phi-p_1)^2
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        Phi-=parms[1];
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        UTorsion+=0.5*parms[0]*SQR(Phi);
        break;
      case HARMONIC_COSINE_DIHEDRAL:
        // (1/2)*p_0*(cos(phi)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        UTorsion+=0.5*parms[0]*SQR(CosPhi-parms[1]);
        break;
      case THREE_COSINE_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UTorsion+=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        break;
      case MM3_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0     [kcal/mol]
        // p_1     [kcal/mol]
        // p_2     [kcal/mol]
        UTorsion+=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        break;
      case CVFF_BLOCKED_DIHEDRAL:
        //
        // ========================================================================
        // p_0     [rad]
        // p_1     [K]
        // p_2     [-]
        // p_3     [rad]
        // p_4     [rad]
        UTorsion+=0.0;
        break;
      case CFF_DIHEDRAL:
        // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UTorsion+=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
        break;
      case CFF_DIHEDRAL2:
        // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UTorsion+=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
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
        // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
        // polymer convention Phi'(trans)=0.
        UTorsion+=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
                  parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
        break;
      case TRAPPE_DIHEDRAL:
        // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
        // ==========================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        UTorsion+=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
        break;
      case TRAPPE_DIHEDRAL_EXTENDED:
        // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
        // =============================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        UTorsion+=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
        break;
      case MOD_TRAPPE_DIHEDRAL:
        /* Salvador modification: 16/08/2016
         add phase in cos function:
         p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
        */
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        Phi-=parms[4];           // shift Phi as Phi+parms[4]
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        ShiftedCosPhi=cos(Phi);
        ShiftedCosPhi2=SQR(ShiftedCosPhi);
        UTorsion+=parms[0]+parms[1]+parms[3]+(parms[1]-3.0*parms[3])*ShiftedCosPhi -2.0*parms[2]*ShiftedCosPhi2 +4.0*parms[3]*ShiftedCosPhi*ShiftedCosPhi2;
        break;
      case CVFF_DIHEDRAL:
        // p_0*(1+cos(p_1*phi-p_2))
        // ========================
        // p_0/k_B [K]
        // p_1     [-]
        // p_2     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        UTorsion+=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
        break;
      case OPLS_DIHEDRAL:
        // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
        // =================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        UTorsion+=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
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
        UTorsion+=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
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
        UTorsion+=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
          2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
          2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
        break;
      case FIXED_DIHEDRAL:
        break;
      default:
        fprintf(stderr, "Undefined Torsion potential in routine 'CalculateTorsionEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UTorsion;
}

void CalculateTorsionEnergyAdsorbates(void)
{
  int i;

  UAdsorbateTorsion[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateTorsion[CurrentSystem]+=CalculateTorsionEnergyAdsorbate(i);
}

REAL CalculateTorsionEnergyCation(int m)
{
  int i,Type,NumberOfTorsions,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL rbc,UTorsion;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2;
  REAL ShiftedCosPhi,ShiftedCosPhi2;
  VECTOR Pb,Pc;
  REAL *parms;

  UTorsion=0.0;
  Type=Cations[CurrentSystem][m].Type;
  NumberOfTorsions=Components[Type].NumberOfTorsions;
  for(i=0;i<NumberOfTorsions;i++)
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

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    parms=Components[Type].TorsionArguments[i];
    switch(Components[Type].TorsionType[i])
    {
      case HARMONIC_DIHEDRAL:
        // (1/2)*p_0*(phi-p_1)^2
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        Phi-=parms[1];
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        UTorsion+=0.5*parms[0]*SQR(Phi);
        break;
      case HARMONIC_COSINE_DIHEDRAL:
        // (1/2)*p_0*(cos(phi)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        UTorsion+=0.5*parms[0]*SQR(CosPhi-parms[1]);
        break;
      case THREE_COSINE_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UTorsion+=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        break;
      case MM3_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0     [kcal/mol]
        // p_1     [kcal/mol]
        // p_2     [kcal/mol]
        UTorsion+=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        break;
      case CVFF_BLOCKED_DIHEDRAL:
        //
        // ========================================================================
        // p_0     [rad]
        // p_1     [K]
        // p_2     [-]
        // p_3     [rad]
        // p_4     [rad]
        UTorsion+=0.0;
        break;
      case CFF_DIHEDRAL:
        // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UTorsion+=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
        break;
      case CFF_DIHEDRAL2:
        // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UTorsion+=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
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
        // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
        // polymer convention Phi'(trans)=0.
        UTorsion+=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
                  parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
        break;
      case TRAPPE_DIHEDRAL:
        // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
        // ==========================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        UTorsion+=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
        break;
      case TRAPPE_DIHEDRAL_EXTENDED:
        // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
        // =============================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        UTorsion+=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
        break;
      case MOD_TRAPPE_DIHEDRAL:
        /* Salvador modification: 16/08/2016
         add phase in cos function:
         p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
        */
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        Phi-=parms[4];           // shift Phi as Phi+parms[4]
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        ShiftedCosPhi=cos(Phi);
        ShiftedCosPhi2=SQR(ShiftedCosPhi);
        UTorsion+=parms[0]+parms[1]+parms[3]+(parms[1]-3.0*parms[3])*ShiftedCosPhi -2.0*parms[2]*ShiftedCosPhi2 +4.0*parms[3]*ShiftedCosPhi*ShiftedCosPhi2;
        break;
      case CVFF_DIHEDRAL:
        // p_0*(1+cos(p_1*phi-p_2))
        // ========================
        // p_0/k_B [K]
        // p_1     [-]
        // p_2     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        UTorsion+=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
        break;
      case OPLS_DIHEDRAL:
        // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
        // =================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        UTorsion+=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
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
        UTorsion+=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
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
        UTorsion+=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
          2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
          2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
        break;
      case FIXED_DIHEDRAL:
        break;
      default:
        fprintf(stderr, "Undefined Torsion potential in routine 'CalculateTorsionEnergyCation' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UTorsion;
}

void CalculateTorsionEnergyCations(void)
{
  int i;

  UCationTorsion[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationTorsion[CurrentSystem]+=CalculateTorsionEnergyCation(i);
}


REAL CalculateImproperTorsionEnergy(int Itype,int Iu)
{
  int A,B,C,D;
  REAL rbc,U;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2;
  VECTOR Pb,Pc;
  REAL *parms;

  A=Components[CurrentComponent].ImproperTorsions[Itype].A;
  B=Components[CurrentComponent].ImproperTorsions[Itype].B;
  C=Components[CurrentComponent].ImproperTorsions[Itype].C;
  D=Components[CurrentComponent].ImproperTorsions[Itype].D;
  parms=Components[CurrentComponent].ImproperTorsionArguments[Itype];

  Dab.x=TrialPositions[Iu][A].x-TrialPositions[Iu][B].x;
  Dab.y=TrialPositions[Iu][A].y-TrialPositions[Iu][B].y;
  Dab.z=TrialPositions[Iu][A].z-TrialPositions[Iu][B].z;

  Dbc.x=TrialPositions[Iu][C].x-TrialPositions[Iu][B].x;
  Dbc.y=TrialPositions[Iu][C].y-TrialPositions[Iu][B].y;
  Dbc.z=TrialPositions[Iu][C].z-TrialPositions[Iu][B].z;
  rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
  Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

  Dcd.x=TrialPositions[Iu][D].x-TrialPositions[Iu][C].x;
  Dcd.y=TrialPositions[Iu][D].y-TrialPositions[Iu][C].y;
  Dcd.z=TrialPositions[Iu][D].z-TrialPositions[Iu][C].z;

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

  // Ensure CosPhi is between -1 and 1.
  CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
  CosPhi2=SQR(CosPhi);

  switch(Components[CurrentComponent].ImproperTorsionType[Itype])
  {
    case HARMONIC_IMPROPER_DIHEDRAL:
      // (1/2)*p_0*(phi-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // potential defined in terms of 'phi' and therefore contains a singularity
      // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
      // same direction as Rbc, and negative otherwise
      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);
      Phi-=parms[1];
      Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
      U=0.5*parms[0]*SQR(Phi);
      break;
    case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      U=0.5*parms[0]*SQR(CosPhi-parms[1]);
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
      // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
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
    case CVFF_IMPROPER_DIHEDRAL:
      // p_0*(1+cos(p_1*phi-p_2))
      // ========================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      // potential defined in terms of 'phi' and therefore contains a singularity
      // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
      // same direction as Rbc, and negative otherwise
      Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
      Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
      Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
      Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
      Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
      Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
      sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
            +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);
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
      U=0;
      break;
    default:
      fprintf(stderr, "Undefined Imporper-Torsion potential in routine 'CalculateImproperTorsionEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return U;
}

REAL CalculateImproperTorsionEnergyAdsorbate(int m)
{
  int i,Type,NumberOfImproperTorsions,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL rbc,UImproperTorsion;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2;
  VECTOR Pb,Pc;
  REAL *parms;

  UImproperTorsion=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfImproperTorsions=Components[Type].NumberOfImproperTorsions;
  for(i=0;i<NumberOfImproperTorsions;i++)
  {
    A=Components[Type].ImproperTorsions[i].A;
    B=Components[Type].ImproperTorsions[i].B;
    C=Components[Type].ImproperTorsions[i].C;
    D=Components[Type].ImproperTorsions[i].D;
    parms=(REAL*)&Components[Type].ImproperTorsionArguments[i];

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

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].ImproperTorsionType[i])
    {
      case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
        UImproperTorsion+=0.5*parms[0]*SQR(CosPhi-cos(parms[1]));
        break;
      case HARMONIC_IMPROPER_DIHEDRAL:
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        Phi-=parms[1];
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        UImproperTorsion+=0.5*parms[0]*SQR(Phi);
        break;
      case THREE_COSINE_IMPROPER_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UImproperTorsion+=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        break;
      case MM3_IMPROPER_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0     [kcal/mol]
        // p_1     [kcal/mol]
        // p_2     [kcal/mol]
        UImproperTorsion+=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        break;
      case CFF_IMPROPER_DIHEDRAL:
        // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UImproperTorsion+=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
        break;
      case CFF_IMPROPER_DIHEDRAL2:
        // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UImproperTorsion+=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
        break;
      case SIX_COSINE_IMPROPER_DIHEDRAL:
        // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
        // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
        // polymer convention Phi'(trans)=0.
        UImproperTorsion+=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
                          parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
        break;
      case TRAPPE_IMPROPER_DIHEDRAL:
        // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
        // ==========================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        UImproperTorsion+=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
        break;
      case TRAPPE_IMPROPER_DIHEDRAL_EXTENDED:
        // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
        // =============================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        UImproperTorsion+=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
        break;
      case CVFF_IMPROPER_DIHEDRAL:
        // p_0*(1+cos(p_1*phi-p_2))
        // ========================
        // p_0/k_B [K]
        // p_1     [-]
        // p_2     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        UImproperTorsion+=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
        break;
      case OPLS_IMPROPER_DIHEDRAL:
        // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
        // =================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        UImproperTorsion+=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
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
        UImproperTorsion+=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
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
        UImproperTorsion+=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
          2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
          2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
        break;
      case FIXED_IMPROPER_DIHEDRAL:
        break;
      default:
        fprintf(stderr, "Undefined Imporper-Torsion potential in routine 'CalculateImproperTorsionEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UImproperTorsion;
}

void CalculateImproperTorsionEnergyAdsorbates(void)
{
  int i;

  UAdsorbateImproperTorsion[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateImproperTorsion[CurrentSystem]+=CalculateImproperTorsionEnergyAdsorbate(i);
}

REAL CalculateImproperTorsionEnergyCation(int m)
{
  int i,Type,NumberOfImproperTorsions,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL rbc,UImproperTorsion;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2;
  VECTOR Pb,Pc;
  REAL *parms;

  UImproperTorsion=0.0;
  Type=Cations[CurrentSystem][m].Type;
  NumberOfImproperTorsions=Components[Type].NumberOfImproperTorsions;
  for(i=0;i<NumberOfImproperTorsions;i++)
  {
    A=Components[Type].ImproperTorsions[i].A;
    B=Components[Type].ImproperTorsions[i].B;
    C=Components[Type].ImproperTorsions[i].C;
    D=Components[Type].ImproperTorsions[i].D;
    parms=(REAL*)&Components[Type].ImproperTorsionArguments[i];

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

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].ImproperTorsionType[i])
    {
      case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
        UImproperTorsion+=0.5*parms[0]*SQR(CosPhi-cos(parms[1]));
        break;
      case HARMONIC_IMPROPER_DIHEDRAL:
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        Phi-=parms[1];
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        UImproperTorsion+=0.5*parms[0]*SQR(Phi);
        break;
      case THREE_COSINE_IMPROPER_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UImproperTorsion+=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        break;
      case MM3_IMPROPER_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0     [kcal/mol]
        // p_1     [kcal/mol]
        // p_2     [kcal/mol]
        UImproperTorsion+=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        break;
      case CFF_IMPROPER_DIHEDRAL:
        // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UImproperTorsion+=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
        break;
      case CFF_IMPROPER_DIHEDRAL2:
        // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UImproperTorsion+=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
        break;
      case SIX_COSINE_IMPROPER_DIHEDRAL:
        // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
        // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
        // polymer convention Phi'(trans)=0.
        UImproperTorsion+=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
                          parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
        break;
      case TRAPPE_IMPROPER_DIHEDRAL:
        // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
        // ==========================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        UImproperTorsion+=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
        break;
      case TRAPPE_IMPROPER_DIHEDRAL_EXTENDED:
        // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
        // =============================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        UImproperTorsion+=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
        break;
      case CVFF_IMPROPER_DIHEDRAL:
        // p_0*(1+cos(p_1*phi-p_2))
        // ========================
        // p_0/k_B [K]
        // p_1     [-]
        // p_2     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        UImproperTorsion+=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
        break;
      case OPLS_IMPROPER_DIHEDRAL:
        // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
        // =================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        UImproperTorsion+=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
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
        UImproperTorsion+=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
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
        UImproperTorsion+=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
          2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
          2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
        break;
      case FIXED_IMPROPER_DIHEDRAL:
        break;
      default:
        fprintf(stderr, "Undefined Imporper-Torsion potential in routine 'CalculateImproperTorsionEnergyCation' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UImproperTorsion;
}

void CalculateImproperTorsionEnergyCations(void)
{
  int i;

  UCationImproperTorsion[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationImproperTorsion[CurrentSystem]+=CalculateImproperTorsionEnergyCation(i);
}

REAL CalculateBondBondEnergy(int Itype,int Iu)
{
  int A,B,C;
  REAL *parms;
  REAL rab,rbc;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc;
  REAL UBondBond;

  UBondBond=0.0;

  A=Components[CurrentComponent].BondBonds[Itype].A;
  B=Components[CurrentComponent].BondBonds[Itype].B;
  C=Components[CurrentComponent].BondBonds[Itype].C;
  parms=(REAL*)&Components[CurrentComponent].BondBondArguments[Itype];

  posA=TrialPositions[Iu][A];
  posB=TrialPositions[Iu][B];
  posC=TrialPositions[Iu][C];

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;
  rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));

  switch(Components[CurrentComponent].BondBondType[Itype])
  {
    case CVFF_BOND_BOND_CROSS:
    case CFF_BOND_BOND_CROSS:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      UBondBond+=parms[0]*(rab-parms[1])*(rbc-parms[2]);
      break;
    default:
      fprintf(stderr, "Undefined Bond-Bond potential in routine 'CalculateBondBondEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return UBondBond;
}

REAL CalculateBondBondEnergyAdsorbate(int m)
{
  int i,A,B,C;
  int NumberOfBondBonds,Type;
  REAL *parms;
  REAL rab,rbc;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc;
  REAL UBondBond;

  UBondBond=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBondBonds=Components[Type].NumberOfBondBonds;
  for(i=0;i<NumberOfBondBonds;i++)
  {
    A=Components[Type].BondBonds[i].A;
    B=Components[Type].BondBonds[i].B;
    C=Components[Type].BondBonds[i].C;
    parms=(REAL*)&Components[Type].BondBondArguments[i];

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

    switch(Components[Type].BondBondType[i])
    {
      case CVFF_BOND_BOND_CROSS:
      case CFF_BOND_BOND_CROSS:
        // p_0*(rab-p_1)*(rbc-p_2)
        // =======================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        UBondBond+=parms[0]*(rab-parms[1])*(rbc-parms[2]);
        break;
      default:
        fprintf(stderr, "Undefined Bond-Bond potential in routine 'CalculateBondBondEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBondBond;
}

void CalculateBondBondEnergyAdsorbates(void)
{
  int i;

  UAdsorbateBondBond[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateBondBond[CurrentSystem]+=CalculateBondBondEnergyAdsorbate(i);
}

REAL CalculateBondBondEnergyCation(int m)
{
  int i,A,B,C;
  int NumberOfBondBonds,Type;
  REAL *parms;
  REAL rab,rbc;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc;
  REAL UBondBond;

  UBondBond=0.0;
  Type=Cations[CurrentSystem][m].Type;
  NumberOfBondBonds=Components[Type].NumberOfBondBonds;
  for(i=0;i<NumberOfBondBonds;i++)
  {
    A=Components[Type].BondBonds[i].A;
    B=Components[Type].BondBonds[i].B;
    C=Components[Type].BondBonds[i].C;
    parms=(REAL*)&Components[Type].BondBondArguments[i];

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

    switch(Components[Type].BondBondType[i])
    {
      case CVFF_BOND_BOND_CROSS:
      case CFF_BOND_BOND_CROSS:
        // p_0*(rab-p_1)*(rbc-p_2)
        // =======================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        UBondBond+=parms[0]*(rab-parms[1])*(rbc-parms[2]);
        break;
      default:
        fprintf(stderr, "Undefined Bond-Bond potential in routine 'CalculateBondBondEnergyCation' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBondBond;
}

void CalculateBondBondEnergyCations(void)
{
  int i;

  UCationBondBond[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationBondBond[CurrentSystem]+=CalculateBondBondEnergyCation(i);
}


REAL CalculateBondBendEnergy(int Itype,int Iu)
{
  int A,B,C;
  REAL *parms;
  REAL rab,rbc;
  VECTOR posA,posB,posC;
  REAL cost,theta;
  VECTOR Rab,Rbc;
  REAL UBondBend;

  UBondBend=0.0;

  A=Components[CurrentComponent].BondBends[Itype].A;
  B=Components[CurrentComponent].BondBends[Itype].B;
  C=Components[CurrentComponent].BondBends[Itype].C;
  parms=(REAL*)&Components[CurrentComponent].BondBendArguments[Itype];

  posA=TrialPositions[Iu][A];
  posB=TrialPositions[Iu][B];
  posC=TrialPositions[Iu][C];

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

  cost=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
  theta=acos(cost);

  parms=Components[CurrentComponent].BondBendArguments[Itype];

  switch(Components[CurrentComponent].BondBendType[Itype])
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
      UBondBend=(theta-parms[0])*(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]));
      break;
    case MM3_BOND_BEND_CROSS:
      // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
      // =====================================
      // p_0     [mdyne/rad]
      // p_1     [A]
      // p_2     [A]
      // p_3     [degrees]
      UBondBend=parms[0]*((rab-parms[1])+(rbc-parms[2]))*RAD2DEG*(theta-parms[3]);
      break;
    case TRUNCATED_HARMONIC:
      // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
      // ================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      UBondBend=0.5*parms[0]*SQR(theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8));
      break;
    case SCREENED_HARMONIC:
      // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      // p_3     [A]
      UBondBend=0.5*parms[0]*SQR(theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3]));
      break;
    case SCREENED_VESSAL:
      // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
      // ============================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      // p_3     [A]
      UBondBend=(parms[0]/(8.0*SQR(theta-M_PI)))*SQR(SQR(parms[1]-M_PI)-SQR(theta-M_PI))
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
      UBondBend=parms[0]*(pow(theta,parms[2])*SQR(theta-parms[1])*SQR(theta+parms[1]-2.0*M_PI)
            -0.5*parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*SQR(theta-parms[1])*pow(M_PI-parms[1],(REAL)3.0))
            *exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8));
      break;
    default:
      fprintf(stderr, "Undefined Bond-Bend potential in routine 'CalculateBondBendEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return UBondBend;

}


REAL CalculateBondBendEnergyAdsorbate(int m)
{
  int i;
  int NumberOfBondBends,Type,A,B,C;
  VECTOR posA,posB,posC,Rab,Rbc;
  REAL UBondBend,rab,rbc,cost,theta;
  REAL *parms;

  UBondBend=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBondBends=Components[Type].NumberOfBondBends;
  for(i=0;i<NumberOfBondBends;i++)
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

    cost=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
    theta=acos(cost);

    parms=Components[Type].BondBendArguments[i];

    switch(Components[Type].BondBendType[i])
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
        UBondBend+=(theta-parms[0])*(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]));
        break;
      case MM3_BOND_BEND_CROSS:
        // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
        // =====================================
        // p_0     [mdyne/rad]
        // p_1     [A]
        // p_2     [A]
        // p_3     [degrees]
        UBondBend+=parms[0]*((rab-parms[1])+(rbc-parms[2]))*RAD2DEG*(theta-parms[3]);
        break;
      case TRUNCATED_HARMONIC:
       // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
        // ================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        UBondBend+=0.5*parms[0]*SQR(theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8));
        break;
      case SCREENED_HARMONIC:
        // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        // p_3     [A]
        UBondBend+=0.5*parms[0]*SQR(theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3]));
        break;
      case SCREENED_VESSAL:
        // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
        // ============================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        // p_3     [A]
        UBondBend+=(parms[0]/(8.0*SQR(theta-M_PI)))*SQR(SQR(parms[1]-M_PI)-SQR(theta-M_PI))
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
        UBondBend+=parms[0]*(pow(theta,parms[2])*SQR(theta-parms[1])*SQR(theta+parms[1]-2.0*M_PI)
              -0.5*parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*SQR(theta-parms[1])*pow(M_PI-parms[1],(REAL)3.0))
              *exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8));
        break;
      default:
        fprintf(stderr, "Undefined Bond-Bend potential in routine 'CalculateBondBendEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBondBend;
}

void CalculateBondBendEnergyAdsorbates(void)
{
  int i;

  UAdsorbateBondBend[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateBondBend[CurrentSystem]+=CalculateBondBendEnergyAdsorbate(i);
}

REAL CalculateBondBendEnergyCation(int m)
{
  int i;
  int NumberOfBondBends,Type,A,B,C;
  VECTOR posA,posB,posC,Rab,Rbc;
  REAL UBondBend,rab,rbc,cost,theta;
  REAL *parms;

  UBondBend=0.0;
  Type=Cations[CurrentSystem][m].Type;
  NumberOfBondBends=Components[Type].NumberOfBondBends;
  for(i=0;i<NumberOfBondBends;i++)
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

    cost=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
    theta=acos(cost);

    parms=Components[Type].BondBendArguments[i];

    switch(Components[Type].BondBendType[i])
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
        UBondBend+=(theta-parms[0])*(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]));
        break;
      case MM3_BOND_BEND_CROSS:
        // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
        // =====================================
        // p_0     [mdyne/rad]
        // p_1     [A]
        // p_2     [A]
        // p_3     [degrees]
        UBondBend+=parms[0]*((rab-parms[1])+(rbc-parms[2]))*RAD2DEG*(theta-parms[3]);
        break;
      case TRUNCATED_HARMONIC:
       // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
        // ================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        UBondBend+=0.5*parms[0]*SQR(theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8));
        break;
      case SCREENED_HARMONIC:
        // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        // p_3     [A]
        UBondBend+=0.5*parms[0]*SQR(theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3]));
        break;
      case SCREENED_VESSAL:
        // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
        // ============================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        // p_3     [A]
        UBondBend+=(parms[0]/(8.0*SQR(theta-M_PI)))*SQR(SQR(parms[1]-M_PI)-SQR(theta-M_PI))
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
        UBondBend+=parms[0]*(pow(theta,parms[2])*SQR(theta-parms[1])*SQR(theta+parms[1]-2.0*M_PI)
              -0.5*parms[2]*pow((REAL)M_PI,(REAL)(parms[2]-1.0))*SQR(theta-parms[1])*pow(M_PI-parms[1],(REAL)3.0))
              *exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8));
        break;
      default:
        fprintf(stderr, "Undefined Bond-Bend potential in routine 'CalculateBondBendEnergyCation' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBondBend;
}

void CalculateBondBendEnergyCations(void)
{
  int i;

  UCationBondBend[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationBondBend[CurrentSystem]+=CalculateBondBendEnergyCation(i);
}

REAL CalculateBendBendEnergy(int Itype,int Iu)
{
  int A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL rab,rbc,rbd;
  VECTOR Dab,Dbc,Dbd;
  REAL dot_abc,dot_abd;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL energy;
  REAL *parms;

  energy=0.0;

  A=Components[CurrentComponent].BendBends[Itype].A;
  B=Components[CurrentComponent].BendBends[Itype].B;
  C=Components[CurrentComponent].BendBends[Itype].C;
  D=Components[CurrentComponent].BendBends[Itype].D;

  posA=TrialPositions[Iu][A];
  posB=TrialPositions[Iu][B];
  posC=TrialPositions[Iu][C];
  posD=TrialPositions[Iu][D];

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
  SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

  dot_abd=Dab.x*Dbd.x+Dab.y*Dbd.y+Dab.z*Dbd.z;
  CosTheta2=dot_abd;
  CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
  Theta2=acos(CosTheta2);
  SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

  parms=Components[CurrentComponent].BendBendArguments[Itype];
  switch(Components[CurrentComponent].BendBendType[Itype])
  {
    case CVFF_BEND_BEND_CROSS:
    case CFF_BEND_BEND_CROSS:
      // p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0/k_B [K/rad^2)]
      // p_1     [degrees]
      // p_2     [degrees]
      energy=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
      break;
    case MM3_BEND_BEND_CROSS:
      // -p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0     [mdyne A/rad^2]
      // p_1     [degrees]
      // p_2     [degrees]
      energy=-parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])*(Theta2-parms[2]);
      break;
    default:
      fprintf(stderr, "Undefined Bend-Bend potential in routine 'CalculateBendBendEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }

  return energy;
}

REAL CalculateBendBendEnergyAdsorbate(int m)
{
  int i,Type,NumberOfBendBends,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL rab,rbc,rbd;
  VECTOR Dab,Dbc,Dbd;
  REAL dot_abc,dot_abd;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL UBendBend;
  REAL *parms;

  UBendBend=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBendBends=Components[Type].NumberOfBendBends;
  for(i=0;i<NumberOfBendBends;i++)
  {
    A=Components[Type].BendBends[i].A;
    B=Components[Type].BendBends[i].B;
    C=Components[Type].BendBends[i].C;
    D=Components[Type].BendBends[i].D;
    parms=(REAL*)&Components[Type].BendBendArguments[i];

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
    SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

    dot_abd=Dab.x*Dbd.x+Dab.y*Dbd.y+Dab.z*Dbd.z;
    CosTheta2=dot_abd;
    CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
    Theta2=acos(CosTheta2);
    SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

    switch(Components[Type].BendBendType[i])
    {
      case CVFF_BEND_BEND_CROSS:
      case CFF_BEND_BEND_CROSS:
        // p_0*(Theta1-p_1)*(Theta2-p_2)
        // ===================================
        // p_0/k_B [K/rad^2)]
        // p_1     [degrees]
        // p_2     [degrees]
        UBendBend+=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
        break;
      case MM3_BEND_BEND_CROSS:
        // -p_0*(Theta1-p_1)*(Theta2-p_2)
        // ===================================
        // p_0     [mdyne A/rad^2]
        // p_1     [degrees]
        // p_2     [degrees]
        UBendBend-=parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])*(Theta2-parms[2]);
        break;
      default:
        fprintf(stderr, "Undefined Bend-Bend potential in routine 'CalculateBendBendEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBendBend;
}

void CalculateBendBendEnergyAdsorbates(void)
{
  int i;

  UAdsorbateBendBend[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateBendBend[CurrentSystem]+=CalculateBendBendEnergyAdsorbate(i);
}

REAL CalculateBendBendEnergyCation(int m)
{
  int i,Type,NumberOfBendBends,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL rab,rbc,rbd;
  VECTOR Dab,Dbc,Dbd;
  REAL dot_abc,dot_abd;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL UBendBend;
  REAL *parms;

  UBendBend=0.0;
  Type=Cations[CurrentSystem][m].Type;
  NumberOfBendBends=Components[Type].NumberOfBendBends;
  for(i=0;i<NumberOfBendBends;i++)
  {
    A=Components[Type].BendBends[i].A;
    B=Components[Type].BendBends[i].B;
    C=Components[Type].BendBends[i].C;
    D=Components[Type].BendBends[i].D;
    parms=(REAL*)&Components[Type].BendBendArguments[i];

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
    SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

    dot_abd=Dab.x*Dbd.x+Dab.y*Dbd.y+Dab.z*Dbd.z;
    CosTheta2=dot_abd;
    CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
    Theta2=acos(CosTheta2);
    SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

    switch(Components[Type].BendBendType[i])
    {
      case CVFF_BEND_BEND_CROSS:
      case CFF_BEND_BEND_CROSS:
        // p_0*(Theta1-p_1)*(Theta2-p_2)
        // ===================================
        // p_0/k_B [K/rad^2)]
        // p_1     [degrees]
        // p_2     [degrees]
        UBendBend+=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
        break;
      case MM3_BEND_BEND_CROSS:
        // -p_0*(Theta1-p_1)*(Theta2-p_2)
        // ===================================
        // p_0     [mdyne A/rad^2]
        // p_1     [degrees]
        // p_2     [degrees]
        UBendBend-=parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])*(Theta2-parms[2]);
        break;
      default:
        fprintf(stderr, "Undefined Bend-Bend potential in routine 'CalculateBendBendEnergyCation' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBendBend;
}

void CalculateBendBendEnergyCations(void)
{
  int i;

  UCationBendBend[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationBendBend[CurrentSystem]+=CalculateBendBendEnergyCation(i);
}

REAL CalculateBondTorsionEnergy(int Itype,int Iu)
{
  int A,B,C,D;
  REAL rab,rbc,rcd,temp;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  VECTOR posA,posB,posC,posD;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2,UBondTorsion;
  REAL *parms;

  UBondTorsion=0.0;

  A=Components[CurrentComponent].BondTorsions[Itype].A;
  B=Components[CurrentComponent].BondTorsions[Itype].B;
  C=Components[CurrentComponent].BondTorsions[Itype].C;
  D=Components[CurrentComponent].BondTorsions[Itype].D;
  parms=Components[CurrentComponent].BondTorsionArguments[Itype];

  posA=TrialPositions[Iu][A];
  posB=TrialPositions[Iu][B];
  posC=TrialPositions[Iu][C];
  posD=TrialPositions[Iu][D];

  Dab.x=posA.x-posB.x;
  Dab.y=posA.y-posB.y;
  Dab.z=posA.z-posB.z;
  rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

  Dcb.x=posC.x-posB.x;
  Dcb.y=posC.y-posB.y;
  Dcb.z=posC.z-posB.z;
  rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
  Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

  Ddc.x=posD.x-posC.x;
  Ddc.y=posD.y-posC.y;
  Ddc.z=posD.z-posC.z;
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

  // cos(2*Theta)=2*cos(Theta)**2-1
  // cos(3*Theta)=4*cos(Theta)**3-3*cos(Theta)
  switch(Components[CurrentComponent].BondTorsionType[Itype])
  {
    case MM3_BOND_TORSION_CROSS:
      // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
      // =====================================================================================
      // p_0     [kcal/A mole]
      // p_1     [kcal/A mole]
      // p_2     [kcal/A mole]
      // p_3     [A]
      temp=(rbc-parms[3]);
      UBondTorsion=parms[0]*temp*CosPhi+parms[1]*temp*(2.0*CosPhi2-1.0)+parms[2]*temp*(4.0*CosPhi2*CosPhi-3.0*CosPhi);
      break;
    default:
      fprintf(stderr, "Undefined Bond-Torsion potential in routine 'CalculateBondTorsionEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return UBondTorsion;

}

REAL CalculateBondTorsionEnergyAdsorbate(int m)
{
  int i,A,B,C,D,Type;
  POINT posA,posB,posC,posD;
  REAL UBondTorsion,rab,rbc,rcd,temp;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2;
  REAL *parms;

  UBondTorsion=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  for(i=0;i<Components[Type].NumberOfBondTorsions;i++)
  {
    A=Components[Type].BondTorsions[i].A;
    B=Components[Type].BondTorsions[i].B;
    C=Components[Type].BondTorsions[i].C;
    D=Components[Type].BondTorsions[i].D;
    parms=(REAL*)&Components[Type].BondTorsionArguments[i];

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
    posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;
    rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

    Dcb.x=posC.x-posB.x;
    Dcb.y=posC.y-posB.y;
    Dcb.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
    Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

    Ddc.x=posD.x-posC.x;
    Ddc.y=posD.y-posC.y;
    Ddc.z=posD.z-posC.z;
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

    switch(Components[Type].BondTorsionType[i])
    {
      case MM3_BOND_TORSION_CROSS:
        // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
        // =====================================================================================
        // p_0     [kcal/A mole]
        // p_1     [kcal/A mole]
        // p_2     [kcal/A mole]
        // p_3     [A]
        temp=(rbc-parms[3]);
        UBondTorsion+=parms[0]*temp*CosPhi+parms[1]*temp*(2.0*CosPhi2-1.0)+parms[2]*temp*(4.0*CosPhi2*CosPhi-3.0*CosPhi);
        break;
      default:
        fprintf(stderr, "Undefined Bond-Torsion potential in routine 'CalculateBondTorsionEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBondTorsion;
}


void CalculateBondTorsionEnergyAdsorbates(void)
{
  int i;

  UAdsorbateBondTorsion[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateBondTorsion[CurrentSystem]+=CalculateBondTorsionEnergyAdsorbate(i);
}

REAL CalculateBondTorsionEnergyCation(int m)
{
  int i,A,B,C,D,Type;
  POINT posA,posB,posC,posD;
  REAL UBondTorsion,rab,rbc,rcd,temp;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2;
  REAL *parms;

  UBondTorsion=0.0;
  Type=Cations[CurrentSystem][m].Type;
  for(i=0;i<Components[Type].NumberOfBondTorsions;i++)
  {
    A=Components[Type].BondTorsions[i].A;
    B=Components[Type].BondTorsions[i].B;
    C=Components[Type].BondTorsions[i].C;
    D=Components[Type].BondTorsions[i].D;
    parms=(REAL*)&Components[Type].BondTorsionArguments[i];

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;
    posD=Cations[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;
    rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

    Dcb.x=posC.x-posB.x;
    Dcb.y=posC.y-posB.y;
    Dcb.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
    Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

    Ddc.x=posD.x-posC.x;
    Ddc.y=posD.y-posC.y;
    Ddc.z=posD.z-posC.z;
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

    switch(Components[Type].BondTorsionType[i])
    {
      case MM3_BOND_TORSION_CROSS:
        // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
        // =====================================================================================
        // p_0     [kcal/A mole]
        // p_1     [kcal/A mole]
        // p_2     [kcal/A mole]
        // p_3     [A]
        temp=(rbc-parms[3]);
        UBondTorsion+=parms[0]*temp*CosPhi+parms[1]*temp*(2.0*CosPhi2-1.0)+parms[2]*temp*(4.0*CosPhi2*CosPhi-3.0*CosPhi);
        break;
      default:
        fprintf(stderr, "Undefined Bond-Torsion potential in routine 'CalculateBondTorsionEnergyCation' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBondTorsion;
}

void CalculateBondTorsionEnergyCations(void)
{
  int i;

  UCationBondTorsion[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationBondTorsion[CurrentSystem]+=CalculateBondTorsionEnergyCation(i);
}

REAL CalculateBendTorsionEnergy(int Itype,int Iu)
{
  int A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL UBendTorsion,rab,rbc,rcd;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2;
  VECTOR Pb,Pc;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL sign,Phi,SinPhi;
  REAL *parms;

  UBendTorsion=0.0;

  A=Components[CurrentComponent].BendTorsions[Itype].A;
  B=Components[CurrentComponent].BendTorsions[Itype].B;
  C=Components[CurrentComponent].BendTorsions[Itype].C;
  D=Components[CurrentComponent].BendTorsions[Itype].D;

  posA=TrialPositions[Iu][A];
  posB=TrialPositions[Iu][B];
  posC=TrialPositions[Iu][C];
  posD=TrialPositions[Iu][D];

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
  SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

  dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;
  CosTheta2=-dot_cd/rcd;
  CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
  Theta2=acos(CosTheta2);
  SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

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

  // Ensure CosPhi is between -1 and 1.
  CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
  CosPhi2=SQR(CosPhi);

  parms=Components[CurrentComponent].BendTorsionArguments[Itype];
  switch(Components[CurrentComponent].BendTorsionType[Itype])
  {
    case CVFF_BEND_TORSION_CROSS:
    case CFF_BEND_TORSION_CROSS:
      // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
      // =====================================================================================
      // p_0/k_B [K/rad^3]
      // p_1     [degrees]
      // p_2     [degrees]
      UBendTorsion=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi;
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
      UBendTorsion=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2);
      break;
    case SMOOTHED_THREE_COSINE_DIHEDRAL:
      // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      UBendTorsion=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
             Smoothing(Theta1)*Smoothing(Theta2);
      break;
    case NICHOLAS_DIHEDRAL:
      // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      UBendTorsion=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
             Smoothing(Theta1);
      break;
    case SMOOTHED_CFF_DIHEDRAL:
      // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      UBendTorsion=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
      break;
    case SMOOTHED_CFF_DIHEDRAL2:
      // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      UBendTorsion=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*Smoothing(Theta1)*Smoothing(Theta2);
      break;
    case SMOOTHED_CFF_BEND_TORSION_CROSS:
      // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K/rad^3]
      // p_1     [degrees]
      // p_2     [degrees]
      UBendTorsion=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi*Smoothing(Theta1)*Smoothing(Theta2);
      break;
    default:
      fprintf(stderr, "Undefined Bend-Torsion potential in routine 'CalculateBendTorsionEnergy' ('internal_energy.c')\n");
      exit(0);
      break;
  }
  return UBendTorsion;


}

REAL CalculateBendTorsionEnergyAdsorbate(int m)
{
  int i,A,B,C,D,Type,NumberOfBendTorsions;
  POINT posA,posB,posC,posD;
  REAL UBendTorsion,rab,rbc,rcd;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2;
  VECTOR Pb,Pc;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL sign,Phi,SinPhi;
  REAL *parms;

  UBendTorsion=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBendTorsions=Components[Type].NumberOfBendTorsions;
  for(i=0;i<NumberOfBendTorsions;i++)
  {
    A=Components[Type].BendTorsions[i].A;
    B=Components[Type].BendTorsions[i].B;
    C=Components[Type].BendTorsions[i].C;
    D=Components[Type].BendTorsions[i].D;
    parms=Components[Type].BendTorsionArguments[i];

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
    SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

    dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;
    CosTheta2=-dot_cd/rcd;
    CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
    Theta2=acos(CosTheta2);
    SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

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

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].BendTorsionType[i])
    {
      case CVFF_BEND_TORSION_CROSS:
      case CFF_BEND_TORSION_CROSS:
        // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
        // =====================================================================================
        // p_0/k_B [K/rad^3]
        // p_1     [degrees]
        // p_2     [degrees]
        UBendTorsion+=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi;
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
        UBendTorsion+=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2);
        break;
      case SMOOTHED_THREE_COSINE_DIHEDRAL:
        // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UBendTorsion+=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
               Smoothing(Theta1)*Smoothing(Theta2);
        break;
      case NICHOLAS_DIHEDRAL:
        // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UBendTorsion+=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
               Smoothing(Theta1);
        break;
      case SMOOTHED_CFF_DIHEDRAL:
        // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UBendTorsion+=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
        break;
      case SMOOTHED_CFF_DIHEDRAL2:
        // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UBendTorsion+=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*Smoothing(Theta1)*Smoothing(Theta2);
        break;
      case SMOOTHED_CFF_BEND_TORSION_CROSS:
        // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K/rad^3]
        // p_1     [degrees]
        // p_2     [degrees]
        UBendTorsion=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi*Smoothing(Theta1)*Smoothing(Theta2);
        break;
      default:
        fprintf(stderr, "Undefined Bend-Torsion potential in routine 'CalculateBendTorsionEnergyAdsorbate' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBendTorsion;
}

void CalculateBendTorsionEnergyAdsorbates(void)
{
  int i;

  UAdsorbateBendTorsion[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateBendTorsion[CurrentSystem]+=CalculateBendTorsionEnergyAdsorbate(i);
}

REAL CalculateBendTorsionEnergyCation(int m)
{
  int i,A,B,C,D,Type,NumberOfBendTorsions;
  POINT posA,posB,posC,posD;
  REAL UBendTorsion,rab,rbc,rcd;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2;
  VECTOR Pb,Pc;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL sign,Phi,SinPhi;
  REAL *parms;

  UBendTorsion=0.0;
  Type=Cations[CurrentSystem][m].Type;
  NumberOfBendTorsions=Components[Type].NumberOfBendTorsions;
  for(i=0;i<NumberOfBendTorsions;i++)
  {
    A=Components[Type].BendTorsions[i].A;
    B=Components[Type].BendTorsions[i].B;
    C=Components[Type].BendTorsions[i].C;
    D=Components[Type].BendTorsions[i].D;
    parms=Components[Type].BendTorsionArguments[i];

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
    SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

    dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;
    CosTheta2=-dot_cd/rcd;
    CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
    Theta2=acos(CosTheta2);
    SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

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

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].BendTorsionType[i])
    {
      case CVFF_BEND_TORSION_CROSS:
      case CFF_BEND_TORSION_CROSS:
        // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
        // =====================================================================================
        // p_0/k_B [K/rad^3]
        // p_1     [degrees]
        // p_2     [degrees]
        UBendTorsion+=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi;
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
        UBendTorsion+=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2);
        break;
      case SMOOTHED_THREE_COSINE_DIHEDRAL:
        // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UBendTorsion+=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
               Smoothing(Theta1)*Smoothing(Theta2);
        break;
      case NICHOLAS_DIHEDRAL:
        // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UBendTorsion+=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
               Smoothing(Theta1);
        break;
      case SMOOTHED_CFF_DIHEDRAL:
        // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UBendTorsion+=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
        break;
      case SMOOTHED_CFF_DIHEDRAL2:
        // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        UBendTorsion+=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*Smoothing(Theta1)*Smoothing(Theta2);
        break;
      case SMOOTHED_CFF_BEND_TORSION_CROSS:
        // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K/rad^3]
        // p_1     [degrees]
        // p_2     [degrees]
        UBendTorsion=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi*Smoothing(Theta1)*Smoothing(Theta2);
        break;
      default:
        fprintf(stderr, "Undefined Bend-Torsion potential in routine 'CalculateBendTorsionEnergyCation' ('internal_energy.c')\n");
        exit(0);
        break;
    }
  }
  return UBendTorsion;
}

void CalculateBendTorsionEnergyCations(void)
{
  int i;

  UCationBendTorsion[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationBendTorsion[CurrentSystem]+=CalculateBendTorsionEnergyCation(i);
}


REAL CalculateIntraVDWEnergyAdsorbate(int m)
{
  int i,NumberOfIntraLJ,Type,TypeA,TypeB,A,B;
  REAL UIntraVDW,rr,Scaling;
  VECTOR dr;
  POINT posA,posB;

  UIntraVDW=0.0;
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
      UIntraVDW+=Scaling*PotentialValue(TypeA,TypeB,rr,1.0);
    }
  }
  return UIntraVDW;
}

void CalculateIntraVDWEnergyAdsorbates(void)
{
  int i;

  UAdsorbateIntraVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateIntraVDW[CurrentSystem]+=CalculateIntraVDWEnergyAdsorbate(i);
}

REAL CalculateIntraVDWEnergyCation(int m)
{
  int i,NumberOfIntraLJ,Type,TypeA,TypeB,A,B;
  REAL UIntraVDW,rr,Scaling;
  VECTOR dr;
  POINT posA,posB;

  UIntraVDW=0.0;
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
      UIntraVDW+=Scaling*PotentialValue(TypeA,TypeB,rr,1.0);
    }
  }
  return UIntraVDW;
}

void CalculateIntraVDWEnergyCations(void)
{
  int i;

  UCationIntraVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationIntraVDW[CurrentSystem]+=CalculateIntraVDWEnergyCation(i);
}

REAL CalculateIntraChargeChargeEnergyAdsorbate(int m)
{
  int i,NumberOfIntraChargeChargeCoulomb,Type,TypeA,TypeB,A,B;
  REAL UIntraCoulomb,r,rr,chargeA,chargeB,Scaling;
  VECTOR dr;
  POINT posA,posB;

  UIntraCoulomb=0.0;
  if(ChargeMethod==NONE) return UIntraCoulomb;

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

    UIntraCoulomb+=Scaling*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
  }
  return UIntraCoulomb;
}

void CalculateIntraChargeChargeEnergyAdsorbates(void)
{
  int i;

  UAdsorbateIntraChargeCharge[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateIntraChargeCharge[CurrentSystem]+=CalculateIntraChargeChargeEnergyAdsorbate(i);
}

REAL CalculateIntraChargeChargeEnergyCation(int m)
{
  int i,NumberOfIntraChargeChargeCoulomb,Type,TypeA,TypeB,A,B;
  REAL UIntraCoulomb,r,rr,chargeA,chargeB,Scaling;
  VECTOR dr;
  POINT posA,posB;

  UIntraCoulomb=0.0;
  if(ChargeMethod==NONE) return UIntraCoulomb;

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

    UIntraCoulomb+=Scaling*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
  }
  return UIntraCoulomb;
}

void CalculateIntraChargeChargeEnergyCations(void)
{
  int i;

  UCationIntraChargeCharge[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationIntraChargeCharge[CurrentSystem]+=CalculateIntraChargeChargeEnergyCation(i);
}

REAL CalculateIntraChargeBondDipoleEnergyAdsorbate(int m)
{
  int i,NumberOfIntraChargeBondDipoleCoulomb,Type,A,B,B1,B2;
  REAL r,rr,ri2,length;
  REAL Bt1,cosB;
  VECTOR dr,dipoleB;
  POINT posA,posB,posB1,posB2;
  REAL ChargeA,temp,UIntraCoulomb;
  REAL DipoleMagnitudeB,energy;

  UIntraCoulomb=0.0;
  if(ChargeMethod==NONE) return UIntraCoulomb;

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
    energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
    UIntraCoulomb-=energy;
  }
  return UIntraCoulomb;
}

void CalculateIntraChargeBondDipoleEnergyAdsorbates(void)
{
  int i;

  UAdsorbateIntraChargeBondDipole[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateIntraChargeBondDipole[CurrentSystem]+=CalculateIntraChargeBondDipoleEnergyAdsorbate(i);
}

REAL CalculateIntraChargeBondDipoleEnergyCation(int m)
{
  int i,NumberOfIntraChargeBondDipoleCoulomb,Type,A,B,B1,B2;
  REAL r,rr,ri2,length;
  REAL Bt1,cosB;
  VECTOR dr,dipoleB;
  POINT posA,posB,posB1,posB2;
  REAL ChargeA,temp,UIntraCoulomb;
  REAL DipoleMagnitudeB,energy;

  UIntraCoulomb=0.0;
  if(ChargeMethod==NONE) return UIntraCoulomb;

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
    energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
    UIntraCoulomb-=energy;
  }
  return UIntraCoulomb;
}

void CalculateIntraChargeBondDipoleEnergyCations(void)
{
  int i;

  UCationIntraChargeBondDipole[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationIntraChargeBondDipole[CurrentSystem]+=CalculateIntraChargeBondDipoleEnergyCation(i);
}

REAL CalculateIntraBondDipoleBondDipoleEnergyAdsorbate(int m)
{
  int i,NumberOfIntraBondDipoleBondDipoleCoulomb,Type,A,B,A1,A2,B1,B2;
  REAL r,rr,ri2,rk2,length;
  REAL Bt0,Bt1,Bt2,Bt3,cosA,cosB,cosAB;
  VECTOR dr,dipoleA,dipoleB;
  POINT posA,posB,posA1,posA2,posB1,posB2;
  REAL temp,UIntraCoulomb;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,energy;

  UIntraCoulomb=0.0;
  if(ChargeMethod==NONE) return UIntraCoulomb;

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
    UIntraCoulomb+=energy;
  }
  return UIntraCoulomb;
}

void CalculateIntraBondDipoleBondDipoleEnergyAdsorbates(void)
{
  int i;

  UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=CalculateIntraBondDipoleBondDipoleEnergyAdsorbate(i);
}

REAL CalculateIntraBondDipoleBondDipoleEnergyCation(int m)
{
  int i,NumberOfIntraBondDipoleBondDipoleCoulomb,Type,A,B,A1,A2,B1,B2;
  REAL r,rr,ri2,rk2,length;
  REAL Bt0,Bt1,Bt2,Bt3,cosA,cosB,cosAB;
  VECTOR dr,dipoleA,dipoleB;
  POINT posA,posB,posA1,posA2,posB1,posB2;
  REAL temp,UIntraCoulomb;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,energy;

  UIntraCoulomb=0.0;
  if(ChargeMethod==NONE) return UIntraCoulomb;

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
    dr=ApplyBoundaryCondition(dr);
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
    UIntraCoulomb+=energy;
  }
  return UIntraCoulomb;
}

void CalculateIntraBondDipoleBondDipoleEnergyCations(void)
{
  int i;

  UCationIntraBondDipoleBondDipole[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    UCationIntraBondDipoleBondDipole[CurrentSystem]+=CalculateIntraBondDipoleBondDipoleEnergyCation(i);
}


void CalculateHarmonicBondConstraintEnergy(void)
{
  int m;
  REAL U,r,rr;
  REAL parms0,parms1;
  POINT posA,posB;
  VECTOR dr;

  UDistanceConstraints[CurrentSystem]=0.0;
  for(m=0;m<NumberOfHarmonicDistanceConstraints[CurrentSystem];m++)
  {
    posA=HarmonicDistanceConstraints[CurrentSystem][m][0]->Position;
    posB=HarmonicDistanceConstraints[CurrentSystem][m][1]->Position;

    parms0=HarmonicDistanceConstraintParameters[CurrentSystem][m][0];
    parms1=HarmonicDistanceConstraintParameters[CurrentSystem][m][1];

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    U=0.5*parms0*SQR(r-parms1);

    // add contribution to the energy
    UDistanceConstraints[CurrentSystem]+=U;
  }
}

void CalculateHarmonicAngleConstraintEnergy(void)
{
  int m;
  REAL U;
  REAL CosTheta,Theta;
  REAL rab,rac,rbc;
  POINT posA,posB,posC;
  VECTOR Rab,Rac,Rbc;
  REAL parms0,parms1;

  UAngleConstraints[CurrentSystem]=0.0;
  for(m=0;m<NumberOfHarmonicAngleConstraints[CurrentSystem];m++)
  {
    posA=HarmonicAngleConstraints[CurrentSystem][m][0]->Position;
    posB=HarmonicAngleConstraints[CurrentSystem][m][1]->Position;
    posC=HarmonicAngleConstraints[CurrentSystem][m][2]->Position;

    parms0=HarmonicAngleConstraintParameters[CurrentSystem][m][0];
    parms1=HarmonicAngleConstraintParameters[CurrentSystem][m][1];

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
    Theta=acos(CosTheta);

    U=0.5*parms0*SQR(Theta-parms1);

    // add contribution to the energy
    UAngleConstraints[CurrentSystem]+=U;
  }
}

void CalculateHarmonicDihedralConstraintEnergy(void)
{
  int m;
  POINT posA,posB,posC,posD;
  REAL rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,U;
  VECTOR Pb,Pc;
  REAL parms0,parms1;

  UDihedralConstraints[CurrentSystem]=0.0;
  for(m=0;m<NumberOfHarmonicDihedralConstraints[CurrentSystem];m++)
  {
    posA=HarmonicDihedralConstraints[CurrentSystem][m][0]->Position;
    posB=HarmonicDihedralConstraints[CurrentSystem][m][1]->Position;
    posC=HarmonicDihedralConstraints[CurrentSystem][m][2]->Position;
    posD=HarmonicDihedralConstraints[CurrentSystem][m][3]->Position;

    parms0=HarmonicDihedralConstraintParameters[CurrentSystem][m][0];
    parms1=HarmonicDihedralConstraintParameters[CurrentSystem][m][1];

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;

    Dcb.x=posC.x-posB.x;
    Dcb.y=posC.y-posB.y;
    Dcb.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
    Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

    Ddc.x=posD.x-posC.x;
    Ddc.y=posD.y-posC.y;
    Ddc.z=posD.z-posC.z;
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
    Phi-=parms1;
    Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
    U=0.5*parms0*SQR(Phi);

    UDihedralConstraints[CurrentSystem]+=U;
  }
}
