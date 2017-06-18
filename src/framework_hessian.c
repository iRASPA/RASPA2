/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_hessian.c' is part of RASPA-2.0

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
#include <float.h>
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
#include "rigid.h"
#include "minimization.h"
#include "internal_hessian.h"
#include "numerical.h"


void CalculateDerivativesAtPositionVDW(VECTOR pos,int typeA,REAL *value,VECTOR *first_derivative,
                                       REAL_MATRIX3x3 *second_derivative,REAL *third_derivative)
{
  int i,f;
  int typeB;
  VECTOR dr;
  REAL rr,F,DF,DDF,DDDF;

  *value=0;
  first_derivative->x=0.0;
  first_derivative->y=0.0;
  first_derivative->z=0.0;
  second_derivative->ax=second_derivative->bx=second_derivative->cx=0.0;
  second_derivative->ay=second_derivative->by=second_derivative->cy=0.0;
  second_derivative->az=second_derivative->bz=second_derivative->cz=0.0;
  *third_derivative=0.0;

  for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f];i++)
    {
      typeB=Framework[CurrentSystem].Atoms[f][i].Type;
      dr.x=pos.x-Framework[CurrentSystem].Atoms[f][i].AnisotropicPosition.x;
      dr.y=pos.y-Framework[CurrentSystem].Atoms[f][i].AnisotropicPosition.y;
      dr.z=pos.z-Framework[CurrentSystem].Atoms[f][i].AnisotropicPosition.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      if(rr<CutOffVDWSquared)
      {
        PotentialThirdDerivative(typeA,typeB,rr,&F,&DF,&DDF,&DDDF);

        *value+=F;

        first_derivative->x+=dr.x*DF;
        first_derivative->y+=dr.y*DF;
        first_derivative->z+=dr.z*DF;

        // add contribution to the second derivatives (Hessian matrix)
        second_derivative->ax+=DDF*dr.x*dr.x+DF; second_derivative->bx+=DDF*dr.y*dr.x;    second_derivative->cx+=DDF*dr.z*dr.x;
        second_derivative->ay+=DDF*dr.x*dr.y;    second_derivative->by+=DDF*dr.y*dr.y+DF; second_derivative->cy+=DDF*dr.z*dr.y;
        second_derivative->az+=DDF*dr.x*dr.z;    second_derivative->bz+=DDF*dr.y*dr.z;    second_derivative->cz+=DDF*dr.z*dr.z+DF;

        *third_derivative+=DDDF*dr.x*dr.y*dr.z;
      }
    }
  }
}

REAL CalculateDerivativesAtPositionReal(VECTOR pos,int typeA,REAL *value,VECTOR *first_derivative,
                                       REAL_MATRIX3x3 *second_derivative,REAL *third_derivative)
{
  int i,f;
  int typeB;
  VECTOR dr;
  REAL ChargeA,ChargeB;
  REAL r,rr,F,DF,DDF,DDDF;
  REAL smallest_r;

  *value=0;
  first_derivative->x=0.0;
  first_derivative->y=0.0;
  first_derivative->z=0.0;
  second_derivative->ax=second_derivative->bx=second_derivative->cx=0.0;
  second_derivative->ay=second_derivative->by=second_derivative->cy=0.0;
  second_derivative->az=second_derivative->bz=second_derivative->cz=0.0;
  *third_derivative=0.0;

  ChargeA=1.0;
  smallest_r=DBL_MAX;
  for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f];i++)
    {
      typeB=Framework[CurrentSystem].Atoms[f][i].Type;
      ChargeB=Framework[CurrentSystem].Atoms[f][i].Charge;
      dr.x=pos.x-Framework[CurrentSystem].Atoms[f][i].Position.x;
      dr.y=pos.y-Framework[CurrentSystem].Atoms[f][i].Position.y;
      dr.z=pos.z-Framework[CurrentSystem].Atoms[f][i].Position.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      if(rr<CutOffChargeChargeSquared[CurrentSystem])
      {
        r=sqrt(rr);

        if(r<smallest_r) smallest_r=r;

        F=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(erfc(Alpha[CurrentSystem]*r)/r);

        DF=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
            (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)+erfc(Alpha[CurrentSystem]*r))/
            (r*rr);

        DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
             (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*r*(3.0+2.0*SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)+
              3.0*erfc(Alpha[CurrentSystem]*r))/(rr*rr*r);

        DDDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
             (-2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*r*(15.0+10.0*SQR(Alpha[CurrentSystem])*rr+4.0*SQR(SQR(Alpha[CurrentSystem])*rr))/sqrt(M_PI)
              -15.0*erfc(Alpha[CurrentSystem]*r))/(rr*rr*rr*r);

        *value+=F;

        first_derivative->x+=dr.x*DF;
        first_derivative->y+=dr.y*DF;
        first_derivative->z+=dr.z*DF;

        // add contribution to the second derivatives (Hessian matrix)
        second_derivative->ax+=DDF*dr.x*dr.x+DF; second_derivative->bx+=DDF*dr.y*dr.x;    second_derivative->cx+=DDF*dr.z*dr.x;
        second_derivative->ay+=DDF*dr.x*dr.y;    second_derivative->by+=DDF*dr.y*dr.y+DF; second_derivative->cy+=DDF*dr.z*dr.y;
        second_derivative->az+=DDF*dr.x*dr.z;    second_derivative->bz+=DDF*dr.y*dr.z;    second_derivative->cz+=DDF*dr.z*dr.z+DF;

        *third_derivative+=DDDF*dr.x*dr.y*dr.z;
      }
    }
  }
  return smallest_r;
}


// Hessian: Center of mass - Center of mass
// ========================================
static inline void HessianAtomicPositionPosition(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor)
{
  REAL_MATRIX3x3 Hessian;

  Hessian.ax=f2*dr.x*dr.x+f1; Hessian.bx=f2*dr.y*dr.x;    Hessian.cx=f2*dr.z*dr.x;
  Hessian.ay=f2*dr.x*dr.y;    Hessian.by=f2*dr.y*dr.y+f1; Hessian.cy=f2*dr.z*dr.y;
  Hessian.az=f2*dr.x*dr.z;    Hessian.bz=f2*dr.y*dr.z;    Hessian.cz=f2*dr.z*dr.z+f1;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x]+=ReplicaFactor*Hessian.ax;
    if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y]+=ReplicaFactor*Hessian.ay;
    if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z]+=ReplicaFactor*Hessian.az;
    if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y]+=ReplicaFactor*Hessian.by;
    if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z]+=ReplicaFactor*Hessian.bz;
    if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z]+=ReplicaFactor*Hessian.cz;
  }

  if((index_j.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_j.x][index_j.x]+=ReplicaFactor*Hessian.ax;
  if((index_j.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.x][index_j.y]+=ReplicaFactor*Hessian.ay;
  if((index_j.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.x][index_j.z]+=ReplicaFactor*Hessian.az;
  if((index_j.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.y][index_j.y]+=ReplicaFactor*Hessian.by;
  if((index_j.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.y][index_j.z]+=ReplicaFactor*Hessian.bz;
  if((index_j.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.z][index_j.z]+=ReplicaFactor*Hessian.cz;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]-=Hessian.ax;
    if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]-=Hessian.ay;
    if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]-=Hessian.az;
    if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]-=Hessian.ay;
    if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]-=Hessian.by;
    if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]-=Hessian.bz;
    if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]-=Hessian.az;
    if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]-=Hessian.bz;
    if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]-=Hessian.cz;
  }
}


// Hessian: Center of mass - Orientation
// =====================================
static inline void HessianCenterOfMassOrientationJ(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,INT_VECTOR3 index_j2,int index2,REAL f1,REAL f2,VECTOR dr)
{
  REAL_MATRIX3x3 Hessian;
  VECTOR vecj1,vecj2,vecj3;

  vecj1=DVecX[index2];
  vecj2=DVecY[index2];
  vecj3=DVecZ[index2];

  Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
  Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
  Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

  Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
  Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
  Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

  Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
  Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
  Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

  if((index_j.x>=0)&&(index_j2.x>=0))  HessianMatrix.element[index_j.x][index_j2.x]+=Hessian.ax;
  if((index_j.x>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_j.x][index_j2.y]+=Hessian.ay;
  if((index_j.x>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_j.x][index_j2.z]+=Hessian.az;

  if((index_j.y>=0)&&(index_j2.x>=0))  HessianMatrix.element[index_j.y][index_j2.x]+=Hessian.bx;
  if((index_j.y>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_j.y][index_j2.y]+=Hessian.by;
  if((index_j.y>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_j.y][index_j2.z]+=Hessian.bz;

  if((index_j.z>=0)&&(index_j2.x>=0))  HessianMatrix.element[index_j.z][index_j2.x]+=Hessian.cx;
  if((index_j.z>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_j.z][index_j2.y]+=Hessian.cy;
  if((index_j.z>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_j.z][index_j2.z]+=Hessian.cz;

  Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
  Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
  Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

  Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
  Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
  Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

  Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
  Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
  Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    if((index_i.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.x][index_j2.x]-=Hessian.ax;
    if((index_i.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.x][index_j2.y]-=Hessian.ay;
    if((index_i.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.x][index_j2.z]-=Hessian.az;

    if((index_i.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.y][index_j2.x]-=Hessian.bx;
    if((index_i.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.y][index_j2.y]-=Hessian.by;
    if((index_i.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.y][index_j2.z]-=Hessian.bz;

    if((index_i.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.z][index_j2.x]-=Hessian.cx;
    if((index_i.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.z][index_j2.y]-=Hessian.cy;
    if((index_i.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.z][index_j2.z]-=Hessian.cz;
  }
}

// Hessian: Orientation - Orientation
// ==================================
static inline void HessianOrientationOrientationJ(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_j2,int index2,REAL f1,REAL f2,VECTOR dr)
{
  REAL_MATRIX3x3 Hessian;
  VECTOR vecj1,vecj2,vecj3;
  VECTOR DDvecJAX,DDvecJBY,DDvecJCZ,DDvecJAY,DDvecJAZ,DDvecJBZ;

  vecj1=DVecX[index2];
  vecj2=DVecY[index2];
  vecj3=DVecZ[index2];

  DDvecJAX=DDVecAX[index2];
  DDvecJBY=DDVecBY[index2];
  DDvecJCZ=DDVecCZ[index2];
  DDvecJAY=DDVecAY[index2];
  DDvecJAZ=DDVecAZ[index2];
  DDvecJBZ=DDVecBZ[index2];

  Hessian.ax=-f1*(dr.x*DDvecJAX.x+dr.y*DDvecJAX.y+dr.z*DDvecJAX.z)
             +f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)
             +f1*(vecj1.x*vecj1.x+vecj1.y*vecj1.y+vecj1.z*vecj1.z);

  Hessian.by=-f1*(dr.x*DDvecJBY.x+dr.y*DDvecJBY.y+dr.z*DDvecJBY.z)+
             f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)+
             f1*(vecj2.x*vecj2.x+vecj2.y*vecj2.y+vecj2.z*vecj2.z);

  Hessian.cz=-f1*(dr.x*DDvecJCZ.x+dr.y*DDvecJCZ.y+dr.z*DDvecJCZ.z)+
             f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
             f1*(vecj3.x*vecj3.x+vecj3.y*vecj3.y+vecj3.z*vecj3.z);

  Hessian.ay=-f1*(dr.x*DDvecJAY.x+dr.y*DDvecJAY.y+dr.z*DDvecJAY.z)+
             f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)+
             f1*(vecj1.x*vecj2.x+vecj1.y*vecj2.y+vecj1.z*vecj2.z);

  Hessian.az=-f1*(dr.x*DDvecJAZ.x+dr.y*DDvecJAZ.y+dr.z*DDvecJAZ.z)+
             f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
             f1*(vecj1.x*vecj3.x+vecj1.y*vecj3.y+vecj1.z*vecj3.z);

  Hessian.bz=-f1*(dr.x*DDvecJBZ.x+dr.y*DDvecJBZ.y+dr.z*DDvecJBZ.z)+
             f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
             f1*(vecj2.x*vecj3.x+vecj2.y*vecj3.y+vecj2.z*vecj3.z);

  if((index_j2.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_j2.x][index_j2.x]+=Hessian.ax;
  if((index_j2.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_j2.y][index_j2.y]+=Hessian.by;
  if((index_j2.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_j2.z][index_j2.z]+=Hessian.cz;
  if((index_j2.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_j2.x][index_j2.y]+=Hessian.ay;
  if((index_j2.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_j2.x][index_j2.z]+=Hessian.az;
  if((index_j2.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_j2.y][index_j2.z]+=Hessian.bz;
}

// Hessian: Center of mass - Strain (part I)
// =========================================
static inline void HessianAtomicPositionStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,REAL f1,REAL f2,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;   // xx x + yy x + zz x
        if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;   // xx y + yy y + zz y
        if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;   // xx z + yy z + zz z
      }

      if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;
      if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;
      if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.y*dr.y*dr.x;
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.z*dr.z*dr.x;

            if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.z*dr.z*dr.y;

            if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.y*dr.y*dr.z;
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
          }

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.y*dr.y*dr.x;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.z*dr.z*dr.x;

          if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.z*dr.z*dr.y;

          if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.y*dr.y*dr.z;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            if(index_i.x>=0)  HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
            if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;
            if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+2]+=f2*dr.x*dr.z*dr.x+f1*dr.z;
            if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+3]+=f2*dr.y*dr.y*dr.x;
            if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+4]+=f2*dr.y*dr.z*dr.x;
            if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+5]+=f2*dr.z*dr.z*dr.x;

            if(index_i.y>=0)  HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;
            if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;
            if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+2]+=f2*dr.x*dr.z*dr.y;
            if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+3]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
            if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+4]+=f2*dr.y*dr.z*dr.y+f1*dr.z;
            if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+5]+=f2*dr.z*dr.z*dr.y;

            if(index_i.z>=0)  HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;
            if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.y*dr.z;
            if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+2]+=f2*dr.x*dr.z*dr.z+f1*dr.x;
            if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+3]+=f2*dr.y*dr.y*dr.z;
            if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+4]+=f2*dr.y*dr.z*dr.z+f1*dr.y;
            if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+5]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
          }

          if(index_j.x>=0)  HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
          if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;
          if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+2]-=f2*dr.x*dr.z*dr.x+f1*dr.z;
          if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+3]-=f2*dr.y*dr.y*dr.x;
          if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+4]-=f2*dr.y*dr.z*dr.x;
          if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+5]-=f2*dr.z*dr.z*dr.x;

          if(index_j.y>=0)  HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
          if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;
          if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+2]-=f2*dr.x*dr.z*dr.y;
          if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+3]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
          if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+4]-=f2*dr.y*dr.z*dr.y+f1*dr.z;
          if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+5]-=f2*dr.z*dr.z*dr.y;

          if(index_j.z>=0)  HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
          if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.y*dr.z;
          if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+2]-=f2*dr.x*dr.z*dr.z+f1*dr.x;
          if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+3]-=f2*dr.y*dr.y*dr.z;
          if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+4]-=f2*dr.y*dr.z*dr.z+f1*dr.y;
          if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+5]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.y*dr.y*dr.x;
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.z*dr.x;
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.z*dr.y+f1*dr.z;
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.y*dr.y*dr.z;
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.z*dr.z+f1*dr.y;
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
              }

              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+1]-=f2*dr.y*dr.y*dr.x;
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+2]-=f2*dr.y*dr.z*dr.x;
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+3]-=f2*dr.z*dr.z*dr.x;

              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+2]-=f2*dr.y*dr.z*dr.y+f1*dr.z;
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+3]-=f2*dr.z*dr.z*dr.y;

              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+1]-=f2*dr.y*dr.y*dr.z;
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+2]-=f2*dr.y*dr.z*dr.z+f1*dr.y;
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.z*dr.x+f1*dr.z;
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.y*dr.x;
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;

                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;
                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.z*dr.y;
                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;

                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;
                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.z*dr.z+f1*dr.x;
                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.y*dr.z;
                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
              }

              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.z*dr.x+f1*dr.z;
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+2]-=f2*dr.y*dr.y*dr.x;
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+3]-=f2*dr.z*dr.z*dr.x;

              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.z*dr.y;
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+2]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+3]-=f2*dr.z*dr.z*dr.y;

              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.z*dr.z+f1*dr.x;
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+2]-=f2*dr.y*dr.y*dr.z;
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.y*dr.x;
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.y*dr.z;
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.y*dr.z;
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
              }

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.y*dr.y*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=f2*dr.z*dr.z*dr.x;

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=f2*dr.z*dr.z*dr.y;

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.y*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}

// Hessian: Center of mass - Strain (part II)
// ==========================================
static inline void HessianCenterOfMassStrainJ(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,REAL f1,REAL f2,VECTOR dr,VECTOR posB,VECTOR comB)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1)+(posB.y-comB.y)*(f2*dr.x*dr.y)+(posB.z-comB.z)*(f2*dr.x*dr.z);
        if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x)+(posB.y-comB.y)*(f2*dr.y*dr.y+f1)+(posB.z-comB.z)*(f2*dr.y*dr.z);
        if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x)+(posB.y-comB.y)*(f2*dr.z*dr.y)+(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
      }

      if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1)+(posB.y-comB.y)*(f2*dr.x*dr.y)+(posB.z-comB.z)*(f2*dr.x*dr.z);
      if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x)+(posB.y-comB.y)*(f2*dr.y*dr.y+f1)+(posB.z-comB.z)*(f2*dr.y*dr.z);
      if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x)+(posB.y-comB.y)*(f2*dr.z*dr.y)+(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            if(index_i.x>=0)  HessianMatrix.element[index_i.x][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
            if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+1]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
            if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+2]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

            if(index_i.y>=0)  HessianMatrix.element[index_i.y][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
            if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+1]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
            if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+2]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

            if(index_i.z>=0)  HessianMatrix.element[index_i.z][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
            if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+1]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
            if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+2]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
          }

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

          if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

          if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
          break;
        case REGULAR:
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.y));
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.z));
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=0.5*((posB.z-comB.z)*(f2*dr.x*dr.y)+(posB.y-comB.y)*(f2*dr.x*dr.z));
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

            if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.y+f1));
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.z));
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=0.5*((posB.z-comB.z)*(f2*dr.y*dr.y+f1)+(posB.y-comB.y)*(f2*dr.y*dr.z));
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

            if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.y));
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.z+f1));
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=0.5*((posB.z-comB.z)*(f2*dr.z*dr.y)+(posB.y-comB.y)*(f2*dr.z*dr.z+f1));
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
          }

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.y));
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.z));
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]-=0.5*((posB.z-comB.z)*(f2*dr.x*dr.y)+(posB.y-comB.y)*(f2*dr.x*dr.z));
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

          if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.y+f1));
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.z));
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]-=0.5*((posB.z-comB.z)*(f2*dr.y*dr.y+f1)+(posB.y-comB.y)*(f2*dr.y*dr.z));
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

          if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.y));
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.z+f1));
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]-=0.5*((posB.z-comB.z)*(f2*dr.z*dr.y)+(posB.y-comB.y)*(f2*dr.z*dr.z+f1));
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
          break;
        case REGULAR_UPPER_TRIANGLE:
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=(posB.y-comB.y)*(f2*dr.x*dr.x+f1);
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=(posB.z-comB.z)*(f2*dr.x*dr.x+f1);
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=(posB.z-comB.z)*(f2*dr.x*dr.y);
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

            if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=(posB.y-comB.y)*(f2*dr.y*dr.x);
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=(posB.z-comB.z)*(f2*dr.y*dr.x);
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=(posB.z-comB.z)*(f2*dr.y*dr.y+f1);
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

            if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=(posB.y-comB.y)*(f2*dr.z*dr.x);
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=(posB.z-comB.z)*(f2*dr.z*dr.x);
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=(posB.z-comB.z)*(f2*dr.z*dr.y);
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
          }

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n  ]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=(posB.y-comB.y)*(f2*dr.x*dr.x+f1);
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=(posB.z-comB.z)*(f2*dr.x*dr.x+f1);
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]-=(posB.z-comB.z)*(f2*dr.x*dr.y);
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

          if(index_j.y>=0) HessianMatrix.element[index_j.y][n  ]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=(posB.y-comB.y)*(f2*dr.y*dr.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=(posB.z-comB.z)*(f2*dr.y*dr.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]-=(posB.z-comB.z)*(f2*dr.y*dr.y+f1);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

          if(index_j.z>=0) HessianMatrix.element[index_j.z][n  ]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=(posB.y-comB.y)*(f2*dr.z*dr.x);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=(posB.z-comB.z)*(f2*dr.z*dr.x);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]-=(posB.z-comB.z)*(f2*dr.z*dr.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.x*dr.y)+(posB.y-comB.y)*(f2*dr.x*dr.z));
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.y*dr.y+f1)+(posB.y-comB.y)*(f2*dr.y*dr.z));
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.z*dr.y)+(posB.y-comB.y)*(f2*dr.z*dr.z+f1));
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.x*dr.y)+(posB.y-comB.y)*(f2*dr.x*dr.z));
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.y*dr.y+f1)+(posB.y-comB.y)*(f2*dr.y*dr.z));
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.z*dr.y)+(posB.y-comB.y)*(f2*dr.z*dr.z+f1));
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=0.5*((posB.z-comB.z)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.z));
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=0.5*((posB.z-comB.z)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.z));
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=0.5*((posB.z-comB.z)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.z+f1));
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=0.5*((posB.z-comB.z)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.z));
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=0.5*((posB.z-comB.z)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.z));
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=0.5*((posB.z-comB.z)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.z+f1));
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.y));
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.y+f1));
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.y));
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.y));
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.y+f1));
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.y));
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n  ]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+1]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+2]+=(posB.z-comB.z)*(f2*dr.x*dr.y);
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n  ]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+1]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+2]+=(posB.z-comB.z)*(f2*dr.y*dr.y+f1);
                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n  ]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+1]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+2]+=(posB.z-comB.z)*(f2*dr.z*dr.y);
                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }

              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n  ]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+1]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+2]-=(posB.z-comB.z)*(f2*dr.x*dr.y);
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n  ]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+1]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+2]-=(posB.z-comB.z)*(f2*dr.y*dr.y+f1);
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n  ]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+1]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+2]-=(posB.z-comB.z)*(f2*dr.z*dr.y);
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n  ]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+1]+=(posB.z-comB.z)*(f2*dr.x*dr.x+f1);
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+2]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n  ]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+1]+=(posB.z-comB.z)*(f2*dr.y*dr.x);
                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+2]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n  ]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+1]+=(posB.z-comB.z)*(f2*dr.z*dr.x);
                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+2]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }

              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n  ]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+1]-=(posB.z-comB.z)*(f2*dr.x*dr.x+f1);
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+2]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
              if(index_j.x>=0)  HessianMatrix.element[index_j.x][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n  ]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+1]-=(posB.z-comB.z)*(f2*dr.y*dr.x);
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+2]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
              if(index_j.y>=0)  HessianMatrix.element[index_j.y][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n  ]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+1]-=(posB.z-comB.z)*(f2*dr.z*dr.x);
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+2]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
              if(index_j.z>=0)  HessianMatrix.element[index_j.z][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=(posB.y-comB.y)*(f2*dr.x*dr.x+f1);
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=(posB.y-comB.y)*(f2*dr.y*dr.x);
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=(posB.y-comB.y)*(f2*dr.z*dr.x);
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n  ]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=(posB.y-comB.y)*(f2*dr.x*dr.x+f1);
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n  ]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=(posB.y-comB.y)*(f2*dr.y*dr.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n  ]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=(posB.y-comB.y)*(f2*dr.z*dr.x);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }


}

// Hessian: Orientation - Strain (part I)
// ======================================
static inline void HessianOrientationStrainJ(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_j2,int index2,REAL f1,REAL f2,VECTOR dr)
{
  int n;
  REAL_MATRIX3x3 Hessian;
  VECTOR vecj1,vecj2,vecj3;

  n=NumberOfCoordinatesMinimizationVariables;

  vecj1=DVecX[index2];
  vecj2=DVecY[index2];
  vecj3=DVecZ[index2];

  Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
  Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
  Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

  Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
  Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
  Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

  Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
  Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
  Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n]-=Hessian.ax*dr.x+Hessian.bx*dr.y+Hessian.cx*dr.z;
      if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n]-=Hessian.ay*dr.x+Hessian.by*dr.y+Hessian.cy*dr.z;
      if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n]-=Hessian.az*dr.x+Hessian.bz*dr.y+Hessian.cz*dr.z;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n]-=Hessian.ax*dr.x;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=Hessian.bx*dr.y;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=Hessian.cx*dr.z;

          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n]-=Hessian.ay*dr.x;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=Hessian.by*dr.y;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=Hessian.cy*dr.z;

          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n]-=Hessian.az*dr.x;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=Hessian.bz*dr.y;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=Hessian.cz*dr.z;
          break;
        case REGULAR:
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n]-=Hessian.ax*dr.x;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=0.5*(Hessian.ax*dr.y+Hessian.bx*dr.x);
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=0.5*(Hessian.ax*dr.z+Hessian.cx*dr.x);
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+3]-=Hessian.bx*dr.y;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+4]-=0.5*(Hessian.bx*dr.z+Hessian.cx*dr.y);
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+5]-=Hessian.cx*dr.z;

          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n]-=Hessian.ay*dr.x;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=0.5*(Hessian.ay*dr.y+Hessian.by*dr.x);
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=0.5*(Hessian.ay*dr.z+Hessian.cy*dr.x);
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+3]-=Hessian.by*dr.y;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+4]-=0.5*(Hessian.by*dr.z+Hessian.cy*dr.y);
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+5]-=Hessian.cy*dr.z;

          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n]-=Hessian.az*dr.x;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=0.5*(Hessian.az*dr.y+Hessian.bz*dr.x);
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=0.5*(Hessian.az*dr.z+Hessian.cz*dr.x);
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+3]-=Hessian.bz*dr.y;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+4]-=0.5*(Hessian.bz*dr.z+Hessian.cz*dr.y);
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+5]-=Hessian.cz*dr.z;
          break;
        case REGULAR_UPPER_TRIANGLE:
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n  ]-=Hessian.ax*dr.x;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=Hessian.ax*dr.y;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=Hessian.ax*dr.z;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+3]-=Hessian.bx*dr.y;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+4]-=Hessian.bx*dr.z;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+5]-=Hessian.cx*dr.z;

          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n  ]-=Hessian.ay*dr.x;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=Hessian.ay*dr.y;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=Hessian.ay*dr.z;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+3]-=Hessian.by*dr.y;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+4]-=Hessian.by*dr.z;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+5]-=Hessian.cy*dr.z;

          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n  ]-=Hessian.az*dr.x;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=Hessian.az*dr.y;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=Hessian.az*dr.z;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+3]-=Hessian.bz*dr.y;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+4]-=Hessian.bz*dr.z;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+5]-=Hessian.cz*dr.z;
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n]-=Hessian.ax*dr.x;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=Hessian.bx*dr.y;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=0.5*(Hessian.bx*dr.z+Hessian.cx*dr.y);
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+3]-=Hessian.cx*dr.z;

              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n]-=Hessian.ay*dr.x;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=Hessian.by*dr.y;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=0.5*(Hessian.by*dr.z+Hessian.cy*dr.y);
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+3]-=Hessian.cy*dr.z;

              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n]-=Hessian.az*dr.x;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=Hessian.bz*dr.y;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=0.5*(Hessian.bz*dr.z+Hessian.cz*dr.y);
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+3]-=Hessian.cz*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n]-=Hessian.ax*dr.x;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=0.5*(Hessian.ax*dr.z+Hessian.cx*dr.x);
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=Hessian.bx*dr.y;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+3]-=Hessian.cx*dr.z;

              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n]-=Hessian.ay*dr.x;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=0.5*(Hessian.ay*dr.z+Hessian.cy*dr.x);
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=Hessian.by*dr.y;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+3]-=Hessian.cy*dr.z;

              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n]-=Hessian.az*dr.x;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=0.5*(Hessian.az*dr.z+Hessian.cz*dr.x);
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=Hessian.bz*dr.y;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+3]-=Hessian.cz*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n]-=Hessian.ax*dr.x;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=0.5*(Hessian.ax*dr.y+Hessian.bx*dr.x);
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=Hessian.bx*dr.y;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+3]-=Hessian.cx*dr.z;

              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n]-=Hessian.ay*dr.x;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=0.5*(Hessian.ay*dr.y+Hessian.by*dr.x);
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=Hessian.by*dr.y;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+3]-=Hessian.cy*dr.z;

              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n]-=Hessian.az*dr.x;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=0.5*(Hessian.az*dr.y+Hessian.bz*dr.x);
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=Hessian.bz*dr.y;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+3]-=Hessian.cz*dr.z;
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n  ]-=Hessian.ax*dr.x;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=Hessian.bx*dr.y;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=Hessian.bx*dr.z;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+3]-=Hessian.cx*dr.z;

              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n  ]-=Hessian.ay*dr.x;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=Hessian.by*dr.y;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=Hessian.by*dr.z;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+3]-=Hessian.cy*dr.z;

              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n  ]-=Hessian.az*dr.x;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=Hessian.bz*dr.y;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=Hessian.bz*dr.z;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+3]-=Hessian.cz*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n  ]-=Hessian.ax*dr.x;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=Hessian.ax*dr.z;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=Hessian.bx*dr.y;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+3]-=Hessian.cx*dr.z;

              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n  ]-=Hessian.ay*dr.x;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=Hessian.ay*dr.z;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=Hessian.by*dr.y;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+3]-=Hessian.cy*dr.z;

              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n  ]-=Hessian.az*dr.x;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=Hessian.az*dr.z;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=Hessian.bz*dr.y;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+3]-=Hessian.cz*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n  ]-=Hessian.ax*dr.x;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=Hessian.ax*dr.y;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=Hessian.bx*dr.y;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+3]-=Hessian.cx*dr.z;

              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n  ]-=Hessian.ay*dr.x;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=Hessian.ay*dr.y;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=Hessian.by*dr.y;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+3]-=Hessian.cy*dr.z;

              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n  ]-=Hessian.az*dr.x;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=Hessian.az*dr.y;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=Hessian.bz*dr.y;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+3]-=Hessian.cz*dr.z;
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}

// Hessian: Orientation - Strain (part II)
// =======================================
static inline void HessianOrientationStrainJ_FR(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_j2,int index2,REAL f1,REAL f2,
                                                VECTOR posB,VECTOR comB,VECTOR dr)
{
  int n;
  REAL_MATRIX3x3 Hessian;
  VECTOR vecj1,vecj2,vecj3;

  n=NumberOfCoordinatesMinimizationVariables;

  vecj1=DVecX[index2];
  vecj2=DVecY[index2];
  vecj3=DVecZ[index2];

  Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
  Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
  Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

  Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
  Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
  Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

  Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
  Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
  Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n]-=(posB.x-comB.x)*Hessian.ax+(posB.y-comB.y)*Hessian.bx+(posB.z-comB.z)*Hessian.cx;
      if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n]-=(posB.x-comB.x)*Hessian.ay+(posB.y-comB.y)*Hessian.by+(posB.z-comB.z)*Hessian.cy;
      if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n]-=(posB.x-comB.x)*Hessian.az+(posB.y-comB.y)*Hessian.bz+(posB.z-comB.z)*Hessian.cz;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=(posB.x-comB.x)*Hessian.ax;
          if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=(posB.y-comB.y)*Hessian.bx;
          if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=(posB.z-comB.z)*Hessian.cx;

          if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=(posB.x-comB.x)*Hessian.ay;
          if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=(posB.y-comB.y)*Hessian.by;
          if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=(posB.z-comB.z)*Hessian.cy;

          if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=(posB.x-comB.x)*Hessian.az;
          if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=(posB.y-comB.y)*Hessian.bz;
          if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=(posB.z-comB.z)*Hessian.cz;
          break;
        case REGULAR:
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n]-=(posB.x-comB.x)*Hessian.ax;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=0.5*((posB.y-comB.y)*Hessian.ax+(posB.x-comB.x)*Hessian.bx);
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=0.5*((posB.z-comB.z)*Hessian.ax+(posB.x-comB.x)*Hessian.cx);
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+3]-=(posB.y-comB.y)*Hessian.bx;
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+4]-=0.5*((posB.z-comB.z)*Hessian.bx+(posB.y-comB.y)*Hessian.cx);
          if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+5]-=(posB.z-comB.z)*Hessian.cx;

          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n]-=(posB.x-comB.x)*Hessian.ay;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=0.5*((posB.y-comB.y)*Hessian.ay+(posB.x-comB.x)*Hessian.by);
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=0.5*((posB.z-comB.z)*Hessian.ay+(posB.x-comB.x)*Hessian.cy);
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+3]-=(posB.y-comB.y)*Hessian.by;
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+4]-=0.5*((posB.z-comB.z)*Hessian.by+(posB.y-comB.y)*Hessian.cy);
          if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+5]-=(posB.z-comB.z)*Hessian.cy;

          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n]-=(posB.x-comB.x)*Hessian.az;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=0.5*((posB.y-comB.y)*Hessian.az+(posB.x-comB.x)*Hessian.bz);
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=0.5*((posB.z-comB.z)*Hessian.az+(posB.x-comB.x)*Hessian.cz);
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+3]-=(posB.y-comB.y)*Hessian.bz;
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+4]-=0.5*((posB.z-comB.z)*Hessian.bz+(posB.y-comB.y)*Hessian.cz);
          if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+5]-=(posB.z-comB.z)*Hessian.cz;
          break;
        case REGULAR_UPPER_TRIANGLE:
          if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n  ]-=(posB.x-comB.x)*Hessian.ax;
          if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=(posB.y-comB.y)*Hessian.ax;
          if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=(posB.z-comB.z)*Hessian.ax;
          if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=(posB.y-comB.y)*Hessian.bx;
          if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+4]-=(posB.z-comB.z)*Hessian.bx;
          if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+5]-=(posB.z-comB.z)*Hessian.cx;

          if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n  ]-=(posB.x-comB.x)*Hessian.ay;
          if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=(posB.y-comB.y)*Hessian.ay;
          if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=(posB.z-comB.z)*Hessian.ay;
          if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=(posB.y-comB.y)*Hessian.by;
          if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+4]-=(posB.z-comB.z)*Hessian.by;
          if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+5]-=(posB.z-comB.z)*Hessian.cy;

          if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n  ]-=(posB.x-comB.x)*Hessian.az;
          if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=(posB.y-comB.y)*Hessian.az;
          if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=(posB.z-comB.z)*Hessian.az;
          if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=(posB.y-comB.y)*Hessian.bz;
          if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+4]-=(posB.z-comB.z)*Hessian.bz;
          if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+5]-=(posB.z-comB.z)*Hessian.cz;
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=(posB.x-comB.x)*Hessian.ax;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=(posB.y-comB.y)*Hessian.bx;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=0.5*((posB.z-comB.z)*Hessian.bx+(posB.y-comB.y)*Hessian.cx);
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=(posB.z-comB.z)*Hessian.cx;

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=(posB.x-comB.x)*Hessian.ay;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=(posB.y-comB.y)*Hessian.by;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=0.5*((posB.z-comB.z)*Hessian.by+(posB.y-comB.y)*Hessian.cy);
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=(posB.z-comB.z)*Hessian.cy;

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=(posB.x-comB.x)*Hessian.az;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=(posB.y-comB.y)*Hessian.bz;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=0.5*((posB.z-comB.z)*Hessian.bz+(posB.y-comB.y)*Hessian.cz);
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=(posB.z-comB.z)*Hessian.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=(posB.x-comB.x)*Hessian.ax;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=0.5*((posB.z-comB.z)*Hessian.ax+(posB.x-comB.x)*Hessian.cx);
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=(posB.y-comB.y)*Hessian.bx;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=(posB.z-comB.z)*Hessian.cx;

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=(posB.x-comB.x)*Hessian.ay;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=0.5*((posB.z-comB.z)*Hessian.ay+(posB.x-comB.x)*Hessian.cy);
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=(posB.y-comB.y)*Hessian.by;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=(posB.z-comB.z)*Hessian.cy;

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=(posB.x-comB.x)*Hessian.az;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=0.5*((posB.z-comB.z)*Hessian.az+(posB.x-comB.x)*Hessian.cz);
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=(posB.y-comB.y)*Hessian.bz;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=(posB.z-comB.z)*Hessian.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n]-=(posB.x-comB.x)*Hessian.ax;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+1]-=0.5*((posB.y-comB.y)*Hessian.ax+(posB.x-comB.x)*Hessian.bx);
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+2]-=(posB.y-comB.y)*Hessian.bx;
              if(index_j2.x>=0)  HessianMatrix.element[index_j2.x][n+3]-=(posB.z-comB.z)*Hessian.cx;

              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n]-=(posB.x-comB.x)*Hessian.ay;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+1]-=0.5*((posB.y-comB.y)*Hessian.ay+(posB.x-comB.x)*Hessian.by);
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+2]-=(posB.y-comB.y)*Hessian.by;
              if(index_j2.y>=0)  HessianMatrix.element[index_j2.y][n+3]-=(posB.z-comB.z)*Hessian.cy;

              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n]-=(posB.x-comB.x)*Hessian.az;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+1]-=0.5*((posB.y-comB.y)*Hessian.az+(posB.x-comB.x)*Hessian.bz);
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+2]-=(posB.y-comB.y)*Hessian.bz;
              if(index_j2.z>=0)  HessianMatrix.element[index_j2.z][n+3]-=(posB.z-comB.z)*Hessian.cz;
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n  ]-=(posB.x-comB.x)*Hessian.ax;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=(posB.y-comB.y)*Hessian.bx;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=(posB.z-comB.z)*Hessian.bx;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=(posB.z-comB.z)*Hessian.cx;

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n  ]-=(posB.x-comB.x)*Hessian.ay;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=(posB.y-comB.y)*Hessian.by;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=(posB.z-comB.z)*Hessian.by;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=(posB.z-comB.z)*Hessian.cy;

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n  ]-=(posB.x-comB.x)*Hessian.az;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=(posB.y-comB.y)*Hessian.bz;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=(posB.z-comB.z)*Hessian.bz;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=(posB.z-comB.z)*Hessian.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n  ]-=(posB.x-comB.x)*Hessian.ax;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=(posB.z-comB.z)*Hessian.ax;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=(posB.y-comB.y)*Hessian.bx;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=(posB.z-comB.z)*Hessian.cx;

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n  ]-=(posB.x-comB.x)*Hessian.ay;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=(posB.z-comB.z)*Hessian.ay;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=(posB.y-comB.y)*Hessian.by;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=(posB.z-comB.z)*Hessian.cy;

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n  ]-=(posB.x-comB.x)*Hessian.az;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=(posB.z-comB.z)*Hessian.az;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=(posB.y-comB.y)*Hessian.bz;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=(posB.z-comB.z)*Hessian.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n  ]-=(posB.x-comB.x)*Hessian.ax;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=(posB.y-comB.y)*Hessian.ax;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=(posB.y-comB.y)*Hessian.bx;
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=(posB.z-comB.z)*Hessian.cx;

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n  ]-=(posB.x-comB.x)*Hessian.ay;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=(posB.y-comB.y)*Hessian.ay;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=(posB.y-comB.y)*Hessian.by;
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=(posB.z-comB.z)*Hessian.cy;

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n  ]-=(posB.x-comB.x)*Hessian.az;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=(posB.y-comB.y)*Hessian.az;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=(posB.y-comB.y)*Hessian.bz;
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=(posB.z-comB.z)*Hessian.cz;
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}


// Hessian: Strain - Strain
// ========================
static inline void HessianAtomicStrainStrain(REAL_MATRIX HessianMatrix,REAL f1,REAL f2,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x+2.0*f2*dr.x*dr.x*dr.y*dr.y+2.0*f2*dr.x*dr.x*dr.z*dr.z+
                            f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y+2.0*f2*dr.y*dr.y*dr.z*dr.z+f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x; // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                  // xxyy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.z*dr.z;                  // xxzz
          HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
          HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz
          HessianMatrix.element[n+2][n+2]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
          break;
        case REGULAR:
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;              // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;                  // xxxy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;                  // xxxz
          HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.y*dr.y;                               // xxyy
          HessianMatrix.element[n  ][n+4]+=f2*dr.x*dr.x*dr.y*dr.z;                               // xxyz
          HessianMatrix.element[n  ][n+5]+=f2*dr.x*dr.x*dr.z*dr.z;                               // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+0.5*f1*(dr.x*dr.x+dr.y*dr.y);  // xyxy
          HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.x*dr.z+0.5*f1*dr.y*dr.z;              // xyxz
          HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.y*dr.y+f1*dr.x*dr.y;                  // xyyy
          HessianMatrix.element[n+1][n+4]+=f2*dr.x*dr.y*dr.y*dr.z+0.5*f1*(dr.x*dr.z);            // xyyz
          HessianMatrix.element[n+1][n+5]+=f2*dr.x*dr.y*dr.z*dr.z;                               // xyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.x*dr.z*dr.x*dr.z+0.5*f1*(dr.x*dr.x+dr.z*dr.z);  // xzxz
          HessianMatrix.element[n+2][n+3]+=f2*dr.x*dr.z*dr.y*dr.y;                               // xzyy
          HessianMatrix.element[n+2][n+4]+=f2*dr.x*dr.z*dr.y*dr.z+0.5*f1*dr.x*dr.y;              // xzyz
          HessianMatrix.element[n+2][n+5]+=f2*dr.x*dr.z*dr.z*dr.z+f1*dr.x*dr.z;                  // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;                  // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

          HessianMatrix.element[n+4][n+4]+=f2*dr.y*dr.z*dr.y*dr.z+0.5*f1*(dr.y*dr.y+dr.z*dr.z);  // yzyz
          HessianMatrix.element[n+4][n+5]+=f2*dr.y*dr.z*dr.z*dr.z+f1*dr.y*dr.z;                  // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
          break;
        case REGULAR_UPPER_TRIANGLE:
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x; // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;     // xxxy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;     // xxxz
          HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.y*dr.y;                  // xxyy
          HessianMatrix.element[n  ][n+4]+=f2*dr.x*dr.x*dr.y*dr.z;                  // xxyz
          HessianMatrix.element[n  ][n+5]+=f2*dr.x*dr.x*dr.z*dr.z;                  // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+f1*dr.y*dr.y;     // xyxy
          HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.x*dr.z+f1*dr.y*dr.z;     // xyxz
          HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.y*dr.y;                  // xyyy
          HessianMatrix.element[n+1][n+4]+=f2*dr.x*dr.y*dr.y*dr.z;                  // xyyz
          HessianMatrix.element[n+1][n+5]+=f2*dr.x*dr.y*dr.z*dr.z;                  // xyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.x*dr.z*dr.x*dr.z+f1*dr.z*dr.z;     // xzxz
          HessianMatrix.element[n+2][n+3]+=f2*dr.x*dr.z*dr.y*dr.y;                  // xzyy
          HessianMatrix.element[n+2][n+4]+=f2*dr.x*dr.z*dr.y*dr.z;                  // xzyz
          HessianMatrix.element[n+2][n+5]+=f2*dr.x*dr.z*dr.z*dr.z;                  // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;     // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

          HessianMatrix.element[n+4][n+4]+=f2*dr.y*dr.z*dr.y*dr.z+f1*dr.z*dr.z;     // yzyz
          HessianMatrix.element[n+4][n+5]+=f2*dr.y*dr.z*dr.z*dr.z;                  // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;              // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                               // xxyy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.z;                               // xxyz
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;                  // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.z*dr.y*dr.z+0.5*f1*(dr.y*dr.y+dr.z*dr.z);  // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.z*dr.z*dr.z+f1*dr.y*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;              // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;                  // xxxz
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                               // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.z*dr.x*dr.z+0.5*f1*(dr.x*dr.x+dr.z*dr.z);  // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.z*dr.y*dr.y;                               // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.z*dr.z*dr.z+f1*dr.x*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;              // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;                  // xxxy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                               // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+0.5*f1*(dr.x*dr.x+dr.y*dr.y);  // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.y*dr.y+f1*dr.x*dr.y;                  // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.z*dr.z;                               // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.z;                  // xxyz
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;     // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.z*dr.y*dr.z+f1*dr.z*dr.z;     // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.z*dr.z*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;     // xxxz
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.z*dr.x*dr.z+f1*dr.z*dr.z;     // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.z*dr.y*dr.y;                  // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.z*dr.z*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;     // xxxy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+f1*dr.y*dr.y;     // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.y*dr.y;                  // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.z*dr.z;                  // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}

// Hessian: Strain - Strain (part II)
// ==================================
static inline void HessianCorrectionStrainStrain(REAL_MATRIX HessianMatrix,REAL f1,REAL f2,VECTOR dr,VECTOR posB,VECTOR comB)
{
  int n;
  VECTOR dB,drJ;

  n=NumberOfCoordinatesMinimizationVariables;
  dB.x=posB.x-comB.x;
  dB.y=posB.y-comB.y;
  dB.z=posB.z-comB.z;

  drJ.x=dr.x+dB.x;
  drJ.y=dr.y+dB.y;
  drJ.z=dr.z+dB.z;


  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // atomic correction
      HessianMatrix.element[n][n]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxy
      HessianMatrix.element[n][n]+=2.0*f2*dr.x*dB.x*dr.y*dr.y;              // xxyy
      HessianMatrix.element[n][n]+=2.0*f2*dr.x*dB.x*dr.z*dr.z;              // xxzz

      HessianMatrix.element[n][n]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
      HessianMatrix.element[n][n]+=2.0*f2*dr.y*dB.y*dr.z*dr.z;              // yyzz

      HessianMatrix.element[n][n]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

      // rigid correction
      HessianMatrix.element[n][n]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x;   // xxxy
      HessianMatrix.element[n][n]+=2.0*f2*dr.x*drJ.x*dr.y*dB.y;             // xxyy
      HessianMatrix.element[n][n]+=2.0*f2*dr.x*drJ.x*dr.z*dB.z;             // xxzz

      HessianMatrix.element[n][n]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y;   // yyyy
      HessianMatrix.element[n][n]+=2.0*f2*dr.y*drJ.y*dr.z*dB.z;             // yyzz

      HessianMatrix.element[n][n]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z;   // zzzz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          // atomic correction
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxy
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dB.x*dr.y*dr.y;                  // xxyy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dB.x*dr.z*dr.z;                  // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
          HessianMatrix.element[n+1][n+2]+=f2*dr.y*dB.y*dr.z*dr.z;                  // yyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

          // rigid correction
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x;   // xxxy
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*drJ.x*dr.y*dB.y;                 // xxyy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*drJ.x*dr.z*dB.z;                 // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y;   // yyyy
          HessianMatrix.element[n+1][n+2]+=f2*dr.y*drJ.y*dr.z*dB.z;                 // yyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z;   // zzzz
          break;
        case REGULAR:
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dB.x*dr.x+2.0*f1*dB.x*dr.x; // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dB.x*dr.y+f1*dB.x*dr.y;     // xxxy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dB.x*dr.z+f1*dB.x*dr.z;     // xxxz
          HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.y*dB.x*dr.y;                  // xxyy
          HessianMatrix.element[n  ][n+4]+=f2*dr.x*dr.y*dB.x*dr.z;                  // xxyz
          HessianMatrix.element[n  ][n+5]+=f2*dr.x*dr.z*dB.x*dr.z;                  // xxzz

          HessianMatrix.element[n+1][n+1]+=0.5*(f2*(dr.x*dB.y*dr.x*dr.y+dr.y*dB.x*dr.x*dr.y)+f1*(dB.x*dr.x+dB.y*dr.y));   // xyxy
          HessianMatrix.element[n+1][n+2]+=0.5*(f2*(dr.x*dB.y*dr.x*dr.z+dr.y*dB.x*dr.x*dr.z)+f1*(dB.y*dr.z));             // xyxz
          HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.x*dB.y*dr.y*dr.y+dr.y*dB.x*dr.y*dr.y)+f1*(dB.x*dr.y);               // xyyy
          HessianMatrix.element[n+1][n+4]+=0.5*f2*(dr.x*dB.y*dr.y*dr.z+dr.y*dB.x*dr.y*dr.z)+0.5*f1*(dB.x*dr.z);           // xyyz
          HessianMatrix.element[n+1][n+5]+=0.5*f2*(dr.x*dB.y*dr.z*dr.z+dr.y*dB.x*dr.z*dr.z);                              // xyzz

          HessianMatrix.element[n+2][n+2]+=0.5*(f2*(dr.x*dB.z*dr.x*dr.z+dr.z*dB.x*dr.x*dr.z)+f1*(dB.x*dr.x+dB.z*dr.z));   // xzxz
          HessianMatrix.element[n+2][n+3]+=0.5*f2*(dr.x*dB.z*dr.y*dr.y+dB.x*dr.z*dr.y*dr.y);                              // xzyy
          HessianMatrix.element[n+2][n+4]+=0.5*(f2*(dr.x*dB.z*dr.y*dr.z+dB.x*dr.z*dr.y*dr.z)+f1*dB.x*dr.y);               // xzyz
          HessianMatrix.element[n+2][n+5]+=0.5*f2*(dr.x*dB.z*dr.z*dr.z+dB.x*dr.z*dr.z*dr.z)+f1*dB.x*dr.z;                 // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dr.y*dB.y*dr.y+2.0*f1*dB.y*dr.y;   // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dB.y*dr.y*dr.z+f1*dB.y*dr.z;       // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dB.y*dr.z*dr.z;                    // yyzz

          HessianMatrix.element[n+4][n+4]+=0.5*(f2*(dr.y*dB.z*dr.y*dr.z+dB.y*dr.z*dr.y*dr.z)+f1*(dB.y*dr.y+dB.z*dr.z));   // yzyz
          HessianMatrix.element[n+4][n+5]+=0.5*f2*(dr.y*dB.z*dr.z*dr.z+dB.y*dr.z*dr.z*dr.z)+f1*dB.y*dr.z;                 // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z;

          HessianMatrix.element[n][n]+=f2*dr.x*dB.x*drJ.x*dr.x+f1*dB.x*drJ.x;
          HessianMatrix.element[n][n+1]+=0.5*f2*(dr.x*dB.y+dr.y*dB.x)*drJ.x*dr.x+0.5*f1*dB.y*(drJ.x);
          HessianMatrix.element[n][n+2]+=0.5*f2*(dr.x*dB.z+dr.z*dB.x)*drJ.x*dr.x+0.5*f1*dB.z*(drJ.x);
          HessianMatrix.element[n][n+3]+=0.5*(dr.y*dB.y+dr.y*dB.y)*f2*drJ.x*dr.x;
          HessianMatrix.element[n][n+4]+=0.5*(dr.y*dB.z+dr.z*dB.y)*f2*drJ.x*dr.x;
          HessianMatrix.element[n][n+5]+=0.5*(dr.z*dB.z+dr.z*dB.z)*f2*drJ.x*dr.x;

          HessianMatrix.element[n+1][n+1]+=0.25*f2*(dr.x*dB.y+dr.y*dB.x)*(drJ.x*dr.y+drJ.y*dr.x)+
                                           0.25*f1*(dB.y*drJ.y+dB.x*drJ.x);
          HessianMatrix.element[n+1][n+2]+=0.25*f2*(dr.x*dB.z+dr.z*dB.x)*(drJ.x*dr.y+drJ.y*dr.x)+0.25*f1*dB.z*drJ.y;
          HessianMatrix.element[n+1][n+3]+=0.25*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.x*dr.y+drJ.y*dr.x)
                                           +0.25*f1*(dB.y*drJ.x+dB.y*drJ.x);
          HessianMatrix.element[n+1][n+4]+=0.25*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.x*dr.y+drJ.y*dr.x)+0.25*f1*dB.z*drJ.x;
          HessianMatrix.element[n+1][n+5]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.x*dr.y+drJ.y*dr.x);

          HessianMatrix.element[n+2][n+2]+=0.25*f2*(dr.x*dB.z+dr.z*dB.x)*(drJ.x*dr.z+drJ.z*dr.x)+
                                           0.25*f1*(dB.z*drJ.z+dB.x*drJ.x);
          HessianMatrix.element[n+2][n+3]+=0.25*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.x*dr.z+drJ.z*dr.x);
          HessianMatrix.element[n+2][n+4]+=0.25*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.x*dr.z+drJ.z*dr.x)+0.25*f1*dB.y*drJ.x;
          HessianMatrix.element[n+2][n+5]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.x*dr.z+drJ.z*dr.x)+
                                           0.25*f1*(dB.z*drJ.x+dB.z*drJ.x);

          HessianMatrix.element[n+3][n+3]+=0.5*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.y*dr.y)+0.5*f1*(dB.y*drJ.y+dB.y*drJ.y);
          HessianMatrix.element[n+3][n+4]+=0.5*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.y*dr.y)+0.5*f1*dB.z*drJ.y;
          HessianMatrix.element[n+3][n+5]+=0.5*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.y);

          HessianMatrix.element[n+4][n+4]+=0.25*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.y*dr.z+drJ.z*dr.y)+
                                           0.25*f1*(dB.z*drJ.z+dB.y*drJ.y);
          HessianMatrix.element[n+4][n+5]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.z+drJ.z*dr.y)+
                                           0.25*f1*(dB.z*drJ.y+dB.z*drJ.y);

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dB.z*drJ.z*dr.z+f1*dB.z*drJ.z;
          break;
        case REGULAR_UPPER_TRIANGLE:
          // atomic correction
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dB.x*dr.x*dr.y+f1*dB.x*dr.y;     // xxxy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dB.x*dr.x*dr.z+f1*dB.x*dr.z;     // xxxz
          HessianMatrix.element[n  ][n+3]+=f2*dr.x*dB.x*dr.y*dr.y;                  // xxyy
          HessianMatrix.element[n  ][n+4]+=f2*dr.x*dB.x*dr.y*dr.z;                  // xxyz
          HessianMatrix.element[n  ][n+5]+=f2*dr.x*dB.x*dr.z*dr.z;                  // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.x*dB.y*dr.x*dr.y+f1*dB.y*dr.y;     // xyxy
          HessianMatrix.element[n+1][n+2]+=f2*dr.x*dB.y*dr.x*dr.z+f1*dB.y*dr.z;     // xyxz
          HessianMatrix.element[n+1][n+3]+=f2*dr.x*dB.y*dr.y*dr.y;                  // xyyy
          HessianMatrix.element[n+1][n+4]+=f2*dr.x*dB.y*dr.y*dr.z;                  // xyyz
          HessianMatrix.element[n+1][n+5]+=f2*dr.x*dB.y*dr.z*dr.z;                  // xyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.x*dB.z*dr.x*dr.z+f1*dB.z*dr.z;     // xzxz
          HessianMatrix.element[n+2][n+3]+=f2*dr.x*dB.z*dr.y*dr.y;                  // xzyy
          HessianMatrix.element[n+2][n+4]+=f2*dr.x*dB.z*dr.y*dr.z;                  // xzyz
          HessianMatrix.element[n+2][n+5]+=f2*dr.x*dB.z*dr.z*dr.z;                  // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dB.y*dr.y*dr.z+f1*dB.y*dr.z;     // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dB.y*dr.z*dr.z;                  // yyzz

          HessianMatrix.element[n+4][n+4]+=f2*dr.y*dB.z*dr.y*dr.z+f1*dB.z*dr.z;     // yzyz
          HessianMatrix.element[n+4][n+5]+=f2*dr.y*dB.z*dr.z*dr.z;                  // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

          // rigid correction
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x; // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*drJ.x*dr.x*dB.y+f1*drJ.x*dB.y; // xxxy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*drJ.x*dr.x*dB.z+f1*drJ.x*dB.z; // xxxz
          HessianMatrix.element[n  ][n+3]+=f2*dr.x*drJ.x*dr.y*dB.y;               // xxyy
          HessianMatrix.element[n  ][n+4]+=f2*dr.x*drJ.x*dr.y*dB.z;               // xxyz
          HessianMatrix.element[n  ][n+5]+=f2*dr.x*drJ.x*dr.z*dB.z;               // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.x*drJ.y*dr.x*dB.y+f1*drJ.y*dB.y; // xyxy
          HessianMatrix.element[n+1][n+2]+=f2*dr.x*drJ.y*dr.x*dB.z+f1*drJ.y*dB.z; // xyxz
          HessianMatrix.element[n+1][n+3]+=f2*dr.x*drJ.y*dr.y*dB.y;               // xyyy
          HessianMatrix.element[n+1][n+4]+=f2*dr.x*drJ.y*dr.y*dB.z;               // xyyz
          HessianMatrix.element[n+1][n+5]+=f2*dr.x*drJ.y*dr.z*dB.z;               // xyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.x*drJ.z*dr.x*dB.z+f1*drJ.z*dB.z; // xzxz
          HessianMatrix.element[n+2][n+3]+=f2*dr.x*drJ.z*dr.y*dB.y;               // xzyy
          HessianMatrix.element[n+2][n+4]+=f2*dr.x*drJ.z*dr.y*dB.z;               // xzyz
          HessianMatrix.element[n+2][n+5]+=f2*dr.x*drJ.z*dr.z*dB.z;               // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y; // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*drJ.y*dr.y*dB.z+f1*drJ.y*dB.z; // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*drJ.y*dr.z*dB.z;               // yyzz

          HessianMatrix.element[n+4][n+4]+=f2*dr.y*drJ.z*dr.y*dB.z+f1*drJ.z*dB.z; // yzyz
          HessianMatrix.element[n+4][n+5]+=f2*dr.y*drJ.z*dr.z*dB.z;               // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z; // zzzz
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dB.x*dr.x+2.0*f1*dB.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.y*dB.x*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.y*dB.x*dr.z;                  // xxyz
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.z*dB.x*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dB.y*dr.y+2.0*f1*dB.y*dr.y;   // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dB.y*dr.y*dr.z+f1*dB.y*dr.z;       // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+2][n+2]+=0.5*(f2*(dr.y*dB.z*dr.y*dr.z+dB.y*dr.z*dr.y*dr.z)+f1*(dB.y*dr.y+dB.z*dr.z));   // yzyz
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(dr.y*dB.z*dr.z*dr.z+dB.y*dr.z*dr.z*dr.z)+f1*dB.y*dr.z;                 // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z;

              HessianMatrix.element[n][n]+=f2*dr.x*dB.x*drJ.x*dr.x+f1*dB.x*drJ.x;
              HessianMatrix.element[n][n+1]+=0.5*(dr.y*dB.y+dr.y*dB.y)*f2*drJ.x*dr.x;
              HessianMatrix.element[n][n+2]+=0.5*(dr.y*dB.z+dr.z*dB.y)*f2*drJ.x*dr.x;
              HessianMatrix.element[n][n+3]+=0.5*(dr.z*dB.z+dr.z*dB.z)*f2*drJ.x*dr.x;

              HessianMatrix.element[n+1][n+1]+=0.5*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.y*dr.y)+0.5*f1*(dB.y*drJ.y+dB.y*drJ.y);
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.y*dr.y)+0.5*f1*dB.z*drJ.y;
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.y);

              HessianMatrix.element[n+2][n+2]+=0.25*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.y*dr.z+drJ.z*dr.y)+
                                               0.25*f1*(dB.z*drJ.z+dB.y*drJ.y);
              HessianMatrix.element[n+2][n+3]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.z+drJ.z*dr.y)+
                                               0.25*f1*(dB.z*drJ.y+dB.z*drJ.y);

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*drJ.z*dr.z+f1*dB.z*drJ.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dB.x*dr.x+2.0*f1*dB.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dB.x*dr.z+f1*dB.x*dr.z;     // xxxz
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.y*dB.x*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.z*dB.x*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=0.5*(f2*(dr.x*dB.z*dr.x*dr.z+dr.z*dB.x*dr.x*dr.z)+f1*(dB.x*dr.x+dB.z*dr.z));   // xzxz
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(dr.x*dB.z*dr.y*dr.y+dB.x*dr.z*dr.y*dr.y);                              // xzyy
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.x*dB.z*dr.z*dr.z+dB.x*dr.z*dr.z*dr.z)+f1*dB.x*dr.z;                 // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dB.y*dr.y+2.0*f1*dB.y*dr.y;   // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z;

              HessianMatrix.element[n][n]+=f2*dr.x*dB.x*drJ.x*dr.x+f1*dB.x*drJ.x;
              HessianMatrix.element[n][n+1]+=0.5*f2*(dr.x*dB.z+dr.z*dB.x)*drJ.x*dr.x+0.5*f1*dB.z*(drJ.x);
              HessianMatrix.element[n][n+2]+=0.5*(dr.y*dB.y+dr.y*dB.y)*f2*drJ.x*dr.x;
              HessianMatrix.element[n][n+3]+=0.5*(dr.z*dB.z+dr.z*dB.z)*f2*drJ.x*dr.x;

              HessianMatrix.element[n+1][n+1]+=0.25*f2*(dr.x*dB.z+dr.z*dB.x)*(drJ.x*dr.z+drJ.z*dr.x)+
                                               0.25*f1*(dB.z*drJ.z+dB.x*drJ.x);
              HessianMatrix.element[n+1][n+2]+=0.25*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.x*dr.z+drJ.z*dr.x);
              HessianMatrix.element[n+1][n+3]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.x*dr.z+drJ.z*dr.x)+
                                               0.25*f1*(dB.z*drJ.x+dB.z*drJ.x);

              HessianMatrix.element[n+2][n+2]+=0.5*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.y*dr.y)+0.5*f1*(dB.y*drJ.y+dB.y*drJ.y);
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.y);

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*drJ.z*dr.z+f1*dB.z*drJ.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dB.x*dr.x+2.0*f1*dB.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dB.x*dr.y+f1*dB.x*dr.y;     // xxxy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.y*dB.x*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.z*dB.x*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=0.5*(f2*(dr.x*dB.y*dr.x*dr.y+dr.y*dB.x*dr.x*dr.y)+f1*(dB.x*dr.x+dB.y*dr.y));   // xyxy
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(dr.x*dB.y*dr.y*dr.y+dr.y*dB.x*dr.y*dr.y)+f1*(dB.x*dr.y);               // xyyy
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.x*dB.y*dr.z*dr.z+dr.y*dB.x*dr.z*dr.z);                              // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dB.y*dr.y+2.0*f1*dB.y*dr.y;   // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z;

              HessianMatrix.element[n][n]+=f2*dr.x*dB.x*drJ.x*dr.x+f1*dB.x*drJ.x;
              HessianMatrix.element[n][n+1]+=0.5*f2*(dr.x*dB.y+dr.y*dB.x)*drJ.x*dr.x+0.5*f1*dB.y*(drJ.x);
              HessianMatrix.element[n][n+2]+=0.5*(dr.y*dB.y+dr.y*dB.y)*f2*drJ.x*dr.x;
              HessianMatrix.element[n][n+3]+=0.5*(dr.z*dB.z+dr.z*dB.z)*f2*drJ.x*dr.x;

              HessianMatrix.element[n+1][n+1]+=0.25*f2*(dr.x*dB.y+dr.y*dB.x)*(drJ.x*dr.y+drJ.y*dr.x)+
                                               0.25*f1*(dB.y*drJ.y+dB.x*drJ.x);
              HessianMatrix.element[n+1][n+2]+=0.25*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.x*dr.y+drJ.y*dr.x)
                                               +0.25*f1*(dB.y*drJ.x+dB.y*drJ.x);
              HessianMatrix.element[n+1][n+3]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.x*dr.y+drJ.y*dr.x);

              HessianMatrix.element[n+2][n+2]+=0.5*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.y*dr.y)+0.5*f1*(dB.y*drJ.y+dB.y*drJ.y);
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.y);

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*drJ.z*dr.z+f1*dB.z*drJ.z;
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // atomic correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dB.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dB.x*dr.y*dr.z;                  // xxyz
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dB.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dB.y*dr.y*dr.z+f1*dB.y*dr.z;     // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dB.z*dr.y*dr.z+f1*dB.z*dr.z;     // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dB.z*dr.z*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

              // rigid correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*drJ.x*dr.y*dB.y;               // xxyy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*drJ.x*dr.y*dB.z;               // xxyz
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*drJ.x*dr.z*dB.z;               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y; // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*drJ.y*dr.y*dB.z+f1*drJ.y*dB.z; // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*drJ.y*dr.z*dB.z;               // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*drJ.z*dr.y*dB.z+f1*drJ.z*dB.z; // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*drJ.z*dr.z*dB.z;               // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z; // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // atomic correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dB.x*dr.x*dr.z+f1*dB.x*dr.z;     // xxxz
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dB.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dB.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dB.z*dr.x*dr.z+f1*dB.z*dr.z;     // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dB.z*dr.y*dr.y;                  // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dB.z*dr.z*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

              // rigid correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*drJ.x*dr.x*dB.z+f1*drJ.x*dB.z; // xxxz
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*drJ.x*dr.y*dB.y;               // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*drJ.x*dr.z*dB.z;               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*drJ.z*dr.x*dB.z+f1*drJ.z*dB.z; // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*drJ.z*dr.y*dB.y;               // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*drJ.z*dr.z*dB.z;               // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*drJ.y*dr.z*dB.z;               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z; // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // atomic correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dB.x*dr.x*dr.y+f1*dB.x*dr.y;     // xxxy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dB.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dB.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dB.y*dr.x*dr.y+f1*dB.y*dr.y;     // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dB.y*dr.y*dr.y;                  // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dB.y*dr.z*dr.z;                  // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

              // rigid correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*drJ.x*dr.x*dB.y+f1*drJ.x*dB.y; // xxxy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*drJ.x*dr.y*dB.y;               // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*drJ.x*dr.z*dB.z;               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*drJ.y*dr.x*dB.y+f1*drJ.y*dB.y; // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*drJ.y*dr.y*dB.y;               // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*drJ.y*dr.z*dB.z;               // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*drJ.y*dr.z*dB.z;               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z; // zzzz
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

static inline void GradientStrain(REAL *Gradient,REAL f1,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]+=f1*(dr.x*dr.x+dr.y*dr.y+dr.z*dr.z); // xx + yy + zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n  ]+=f1*dr.x*dr.x; // xx
          Gradient[n+1]+=f1*dr.y*dr.y; // yy
          Gradient[n+2]+=f1*dr.z*dr.z; // zz
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n  ]+=f1*dr.x*dr.x;  // xx
          Gradient[n+1]+=f1*dr.x*dr.y;  // xy
          Gradient[n+2]+=f1*dr.x*dr.z;  // xz
          Gradient[n+3]+=f1*dr.y*dr.y;  // yy
          Gradient[n+4]+=f1*dr.y*dr.z;  // yz
          Gradient[n+5]+=f1*dr.z*dr.z;  // zz
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]+=f1*dr.x*dr.x;  // xx
              Gradient[n+1]+=f1*dr.y*dr.y;  // yy
              Gradient[n+2]+=f1*dr.y*dr.z;  // yz
              Gradient[n+3]+=f1*dr.z*dr.z;  // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]+=f1*dr.x*dr.x;  // xx
              Gradient[n+1]+=f1*dr.x*dr.z;  // xz
              Gradient[n+2]+=f1*dr.y*dr.y;  // yy
              Gradient[n+3]+=f1*dr.z*dr.z;  // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]+=f1*dr.x*dr.x;  // xx
              Gradient[n+1]+=f1*dr.x*dr.y;  // xy
              Gradient[n+2]+=f1*dr.y*dr.y;  // yy
              Gradient[n+3]+=f1*dr.z*dr.z;  // zz
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}

static inline void GradientStrainJ(REAL *Gradient,REAL f1,VECTOR dr,VECTOR posB,VECTOR comB)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n  ]+=f1*(dr.x*(posB.x-comB.x)+dr.y*(posB.y-comB.y)+dr.z*(posB.z-comB.z)); // xx + yy + zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n  ]+=f1*dr.x*(posB.x-comB.x); // xx
          Gradient[n+1]+=f1*dr.y*(posB.y-comB.y); // yy
          Gradient[n+2]+=f1*dr.z*(posB.z-comB.z); // zz
          break;
        case REGULAR:
          Gradient[n  ]+=f1*dr.x*(posB.x-comB.x);                            // xx
          Gradient[n+1]+=f1*0.5*(dr.x*(posB.y-comB.y)+dr.y*(posB.x-comB.x)); // xy
          Gradient[n+2]+=f1*0.5*(dr.x*(posB.z-comB.z)+dr.z*(posB.x-comB.x)); // xz
          Gradient[n+3]+=f1*(posB.y-comB.y)*dr.y;                            // yy
          Gradient[n+4]+=f1*0.5*(dr.y*(posB.z-comB.z)+dr.z*(posB.y-comB.y)); // yz
          Gradient[n+5]+=f1*dr.z*(posB.z-comB.z);                            // zz
          break;
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n  ]+=f1*dr.x*(posB.x-comB.x); // xx
          Gradient[n+1]+=f1*dr.x*(posB.y-comB.y); // xy
          Gradient[n+2]+=f1*dr.x*(posB.z-comB.z); // xz
          Gradient[n+3]+=f1*dr.y*(posB.y-comB.y); // yy
          Gradient[n+4]+=f1*dr.y*(posB.z-comB.z); // yz
          Gradient[n+5]+=f1*dr.z*(posB.z-comB.z); // zz
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x);                            // xx
              Gradient[n+1]+=f1*(posB.y-comB.y)*dr.y;                            // yy
              Gradient[n+2]+=f1*0.5*(dr.y*(posB.z-comB.z)+dr.z*(posB.y-comB.y)); // yz
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z);                            // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x);                            // xx
              Gradient[n+1]+=f1*0.5*(dr.x*(posB.z-comB.z)+dr.z*(posB.x-comB.x)); // xz
              Gradient[n+2]+=f1*(posB.y-comB.y)*dr.y;                            // yy
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z);                            // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x);                            // xx
              Gradient[n+1]+=f1*0.5*(dr.x*(posB.y-comB.y)+dr.y*(posB.x-comB.x)); // xy
              Gradient[n+2]+=f1*(posB.y-comB.y)*dr.y;                            // yy
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z);                            // zz
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x); // xx
              Gradient[n+1]+=f1*dr.y*(posB.y-comB.y); // yy
              Gradient[n+2]+=f1*dr.y*(posB.z-comB.z); // yz
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z); // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x); // xx
              Gradient[n+1]+=f1*dr.x*(posB.z-comB.z); // xz
              Gradient[n+2]+=f1*dr.y*(posB.y-comB.y); // yy
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z); // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x); // xx
              Gradient[n+1]+=f1*dr.x*(posB.y-comB.y); // xy
              Gradient[n+2]+=f1*dr.y*(posB.y-comB.y); // yy
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z); // zz
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

/*********************************************************************************************************
 * Name       | ComputeFrameworkAdsorbateVDWHessian                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the contributions of the framework-adsorbate Van der Waals interactions to the  *
 *            | total energy, gradient, and Hessian-matrix.                                              *
 *            | interactions.                                                                            *
 * Parameters | Energy: pointer to the total energy                                                      *
 *            | Gradient: array of the gradients                                                         *
 *            | HessianMatrix: the Hessian matrix                                                        *
 *            | StrainDerivative: pointer to a 3x3x matrix of the strain derivative                      *
 *            | ComputeGradient: Boolean whether or not to compute the gradient                          *
 *            | ComputeHessian: Boolean whether or not to compute the Hessian matrix                     *
 * Note       |                                                                                          *
 * Used in    | void EvaluateDerivatives(...)                                                            *
 *********************************************************************************************************/

void ComputeFrameworkAdsorbateVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                         REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr;
  REAL energy,f1,f2;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR pos,comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;
  REAL scalingB;

  f1=f2=0.0;
  index2=0;
  // first loop over adsorbate molecules

  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
    {
      index_i=Framework[CurrentSystem].Atoms[fr][i].HessianIndex;
      posA=Framework[CurrentSystem].Atoms[fr][i].Position;
      typeA=Framework[CurrentSystem].Atoms[fr][i].Type;
      ChargeA=Framework[CurrentSystem].Atoms[fr][i].Charge;

        // second loop over adsorbates
      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
            {
              TypeMolB=Adsorbates[CurrentSystem][J].Type;
              for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
              {
                for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                {
                  j=Components[TypeMolB].Groups[jg].Atoms[ja];

                  if(Components[TypeMolB].Groups[jg].Rigid)
                  {
                    index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                    index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                    index2=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                  }
                  else
                  {
                    index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                    index_j2=UNDEFINED_INT_VECTOR3;
                    index2=-1;
                  }

                  typeB=Adsorbates[CurrentSystem][J].Atoms[j].Type;
                  posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;
                  ChargeB=Adsorbates[CurrentSystem][J].Atoms[j].Charge;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffVDWSquared)
                  {
                    scalingB=Adsorbates[CurrentSystem][J].Atoms[j].CFVDWScalingParameter;

                    PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2,scalingB);

                    // add contribution to the energy
                    *Energy+=energy;

                    //if((index_i<0)&&(index_j<0)&&(index_j2<0)) continue;

                    StrainDerivative->ax+=f1*dr.x*dr.x;
                    StrainDerivative->bx+=f1*dr.x*dr.y;
                    StrainDerivative->cx+=f1*dr.x*dr.z;

                    StrainDerivative->ay+=f1*dr.y*dr.x;
                    StrainDerivative->by+=f1*dr.y*dr.y;
                    StrainDerivative->cy+=f1*dr.y*dr.z;

                    StrainDerivative->az+=f1*dr.z*dr.x;
                    StrainDerivative->bz+=f1*dr.z*dr.y;
                    StrainDerivative->cz+=f1*dr.z*dr.z;

                    if(Components[TypeMolB].Groups[jg].Rigid)
                    {
                      comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                      comB.x+=ReplicaShift[ncell].x;
                      comB.y+=ReplicaShift[ncell].y;
                      comB.z+=ReplicaShift[ncell].z;

                      pos=Components[TypeMolB].Positions[j];

                      temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                      temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                      temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                      StrainDerivative->ax+=(posB.x-comB.x)*f1*dr.x;
                      StrainDerivative->ay+=temp1;
                      StrainDerivative->az+=temp2;
                      StrainDerivative->bx+=temp1;
                      StrainDerivative->by+=(posB.y-comB.y)*f1*dr.y;
                      StrainDerivative->bz+=temp3;
                      StrainDerivative->cx+=temp2;
                      StrainDerivative->cy+=temp3;
                      StrainDerivative->cz+=(posB.z-comB.z)*f1*dr.z;
                    }

                    // add contribution to the first derivatives
                    if(ComputeGradient)
                    {
                      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      {
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;
                      }

                      if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                      if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                      if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                        // add contribution to the first derivatives
                        if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                        if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                        if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,1.0);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,f1,f2,dr);
                      HessianAtomicStrainStrain(HessianMatrix,f1,f2,dr);


                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianCenterOfMassStrainJ(HessianMatrix,index_i,index_j,f1,f2,dr,posB,comB);
                        HessianOrientationStrainJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianOrientationStrainJ_FR(HessianMatrix,index_j2,index2,f1,f2,posB,comB,dr);

                        HessianCorrectionStrainStrain(HessianMatrix,f1,f2,dr,posB,comB);
                      }
                    }
                  }
                }
              }
            }
            ncell++;
          }
    }
  }


}

void ComputeFrameworkAdsorbateChargeChargeHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                                  REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr,U;
  REAL f1,f2;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;
  REAL scalingB;

  f1=f2=0.0;
  index2=0;
  // first loop over adsorbate molecules
  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
    {
      index_i=Framework[CurrentSystem].Atoms[fr][i].HessianIndex;
      posA=Framework[CurrentSystem].Atoms[fr][i].Position;
      typeA=Framework[CurrentSystem].Atoms[fr][i].Type;
      ChargeA=Framework[CurrentSystem].Atoms[fr][i].Charge;

      // second loop over adsorbates
      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
            {
              TypeMolB=Adsorbates[CurrentSystem][J].Type;
              for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
              {
                for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                {
                  j=Components[TypeMolB].Groups[jg].Atoms[ja];

                  if(Components[TypeMolB].Groups[jg].Rigid)
                  {
                    index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                    index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                    index2=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                  }
                  else
                  {
                    index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                    index_j2=UNDEFINED_INT_VECTOR3;
                    index2=-1;
                  }

                  typeB=Adsorbates[CurrentSystem][J].Atoms[j].Type;
                  posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  {
                    scalingB=Adsorbates[CurrentSystem][J].Atoms[j].CFChargeScalingParameter;
                    ChargeB=scalingB*Adsorbates[CurrentSystem][J].Atoms[j].Charge;

                    PotentialSecondDerivativeCoulombic(ChargeA,ChargeB,rr,&U,&f1,&f2);

                    *Energy+=U;

                    StrainDerivative->ax+=f1*dr.x*dr.x;
                    StrainDerivative->bx+=f1*dr.x*dr.y;
                    StrainDerivative->cx+=f1*dr.x*dr.z;

                    StrainDerivative->ay+=f1*dr.y*dr.x;
                    StrainDerivative->by+=f1*dr.y*dr.y;
                    StrainDerivative->cy+=f1*dr.y*dr.z;

                    StrainDerivative->az+=f1*dr.z*dr.x;
                    StrainDerivative->bz+=f1*dr.z*dr.y;
                    StrainDerivative->cz+=f1*dr.z*dr.z;

                    if(Components[TypeMolB].Groups[jg].Rigid)
                    {
                      comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                      comB.x+=ReplicaShift[ncell].x;
                      comB.y+=ReplicaShift[ncell].y;
                      comB.z+=ReplicaShift[ncell].z;

                      temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                      temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                      temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                      StrainDerivative->ax+=(posB.x-comB.x)*f1*dr.x;
                      StrainDerivative->ay+=temp1;
                      StrainDerivative->az+=temp2;
                      StrainDerivative->bx+=temp1;
                      StrainDerivative->by+=(posB.y-comB.y)*f1*dr.y;
                      StrainDerivative->bz+=temp3;
                      StrainDerivative->cx+=temp2;
                      StrainDerivative->cy+=temp3;
                      StrainDerivative->cz+=(posB.z-comB.z)*f1*dr.z;
                    }

                    // add contribution to the first derivatives
                    if(ComputeGradient)
                    {
                      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      {
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;
                      }

                      if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                      if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                      if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                         // add contribution to the first derivatives
                        if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                        if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                        if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,1.0);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,f1,f2,dr);
                      HessianAtomicStrainStrain(HessianMatrix,f1,f2,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianCenterOfMassStrainJ(HessianMatrix,index_i,index_j,f1,f2,dr,posB,comB);
                        HessianOrientationStrainJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianOrientationStrainJ_FR(HessianMatrix,index_j2,index2,f1,f2,posB,comB,dr);

                        HessianCorrectionStrainStrain(HessianMatrix,f1,f2,dr,posB,comB);
                      }
                    }
                  }
                }
              }
            }
            ncell++;
          }
    }
  }
}

/*********************************************************************************************************
 * Name       | ComputeFrameworkCationVDWHessian                                                         *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the contributions of the framework-cation Van der Waals interactions to the     *
 *            | total energy, gradient, and Hessian-matrix.                                              *
 *            | interactions.                                                                            *
 * Parameters | Energy: pointer to the total energy                                                      *
 *            | Gradient: array of the gradients                                                         *
 *            | HessianMatrix: the Hessian matrix                                                        *
 *            | StrainDerivative: pointer to a 3x3x matrix of the strain derivative                      *
 *            | ComputeGradient: Boolean whether or not to compute the gradient                          *
 *            | ComputeHessian: Boolean whether or not to compute the Hessian matrix                     *
 * Note       |                                                                                          *
 * Used in    | void EvaluateDerivatives(...)                                                            *
 *********************************************************************************************************/

void ComputeFrameworkCationVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                      REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr;
  REAL energy,f1,f2;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR pos,comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;
  REAL scalingB;

  f1=f2=0.0;
  index2=0;

  // first loop over framework atoms
  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
    {
      index_i=Framework[CurrentSystem].Atoms[fr][i].HessianIndex;
      posA=Framework[CurrentSystem].Atoms[fr][i].Position;
      typeA=Framework[CurrentSystem].Atoms[fr][i].Type;
      ChargeA=Framework[CurrentSystem].Atoms[fr][i].Charge;

        // second loop over adsorbates
      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
            {
              TypeMolB=Cations[CurrentSystem][J].Type;
              for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
              {
                for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                {
                  j=Components[TypeMolB].Groups[jg].Atoms[ja];

                  if(Components[TypeMolB].Groups[jg].Rigid)
                  {
                    index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                    index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                    index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                  }
                  else
                  {
                    index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                    index_j2=UNDEFINED_INT_VECTOR3;
                    index2=-1;
                  }

                  typeB=Cations[CurrentSystem][J].Atoms[j].Type;
                  posB=Cations[CurrentSystem][J].Atoms[j].Position;
                  ChargeB=Cations[CurrentSystem][J].Atoms[j].Charge;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffVDWSquared)
                  {
                    scalingB=Cations[CurrentSystem][J].Atoms[j].CFVDWScalingParameter;

                    PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2,scalingB);

                    // add contribution to the energy
                    *Energy+=energy;

                    //if((index_i<0)&&(index_j<0)&&(index_j2<0)) continue;

                    StrainDerivative->ax+=f1*dr.x*dr.x;
                    StrainDerivative->bx+=f1*dr.x*dr.y;
                    StrainDerivative->cx+=f1*dr.x*dr.z;

                    StrainDerivative->ay+=f1*dr.y*dr.x;
                    StrainDerivative->by+=f1*dr.y*dr.y;
                    StrainDerivative->cy+=f1*dr.y*dr.z;

                    StrainDerivative->az+=f1*dr.z*dr.x;
                    StrainDerivative->bz+=f1*dr.z*dr.y;
                    StrainDerivative->cz+=f1*dr.z*dr.z;

                    if(Components[TypeMolB].Groups[jg].Rigid)
                    {
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                      comB.x+=ReplicaShift[ncell].x;
                      comB.y+=ReplicaShift[ncell].y;
                      comB.z+=ReplicaShift[ncell].z;

                      pos=Components[TypeMolB].Positions[j];

                      temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                      temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                      temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                      StrainDerivative->ax+=(posB.x-comB.x)*f1*dr.x;
                      StrainDerivative->ay+=temp1;
                      StrainDerivative->az+=temp2;
                      StrainDerivative->bx+=temp1;
                      StrainDerivative->by+=(posB.y-comB.y)*f1*dr.y;
                      StrainDerivative->bz+=temp3;
                      StrainDerivative->cx+=temp2;
                      StrainDerivative->cy+=temp3;
                      StrainDerivative->cz+=(posB.z-comB.z)*f1*dr.z;
                    }

                    // add contribution to the first derivatives
                    if(ComputeGradient)
                    {
                      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      {
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;
                      }

                      if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                      if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                      if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                        // add contribution to the first derivatives
                        if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                        if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                        if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,1.0);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,f1,f2,dr);
                      HessianAtomicStrainStrain(HessianMatrix,f1,f2,dr);


                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianCenterOfMassStrainJ(HessianMatrix,index_i,index_j,f1,f2,dr,posB,comB);
                        HessianOrientationStrainJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianOrientationStrainJ_FR(HessianMatrix,index_j2,index2,f1,f2,posB,comB,dr);

                        HessianCorrectionStrainStrain(HessianMatrix,f1,f2,dr,posB,comB);
                      }
                    }
                  }
                }
              }
            }
            ncell++;
          }
    }
  }
}

void ComputeFrameworkCationChargeChargeHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                               REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr,U;
  REAL f1,f2;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;
  REAL scalingB;

  f1=f2=0.0;
  index2=0;
  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
    {
      index_i=Framework[CurrentSystem].Atoms[fr][i].HessianIndex;
      posA=Framework[CurrentSystem].Atoms[fr][i].Position;
      typeA=Framework[CurrentSystem].Atoms[fr][i].Type;
      ChargeA=Framework[CurrentSystem].Atoms[fr][i].Charge;

      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
            {
              TypeMolB=Cations[CurrentSystem][J].Type;
              for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
              {
                for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                {
                  j=Components[TypeMolB].Groups[jg].Atoms[ja];

                  if(Components[TypeMolB].Groups[jg].Rigid)
                  {
                    index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                    index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                    index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                  }
                  else
                  {
                    index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                    index_j2=UNDEFINED_INT_VECTOR3;
                    index2=-1;
                  }

                  typeB=Cations[CurrentSystem][J].Atoms[j].Type;
                  posB=Cations[CurrentSystem][J].Atoms[j].Position;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  {
                    scalingB=Cations[CurrentSystem][J].Atoms[j].CFChargeScalingParameter;
                    ChargeB=scalingB*Cations[CurrentSystem][J].Atoms[j].Charge;

                    PotentialSecondDerivativeCoulombic(ChargeA,ChargeB,rr,&U,&f1,&f2);

                    *Energy+=U;

                    StrainDerivative->ax+=f1*dr.x*dr.x;
                    StrainDerivative->bx+=f1*dr.x*dr.y;
                    StrainDerivative->cx+=f1*dr.x*dr.z;

                    StrainDerivative->ay+=f1*dr.y*dr.x;
                    StrainDerivative->by+=f1*dr.y*dr.y;
                    StrainDerivative->cy+=f1*dr.y*dr.z;

                    StrainDerivative->az+=f1*dr.z*dr.x;
                    StrainDerivative->bz+=f1*dr.z*dr.y;
                    StrainDerivative->cz+=f1*dr.z*dr.z;

                    if(Components[TypeMolB].Groups[jg].Rigid)
                    {
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                      comB.x+=ReplicaShift[ncell].x;
                      comB.y+=ReplicaShift[ncell].y;
                      comB.z+=ReplicaShift[ncell].z;

                      temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                      temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                      temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                      StrainDerivative->ax+=(posB.x-comB.x)*f1*dr.x;
                      StrainDerivative->ay+=temp1;
                      StrainDerivative->az+=temp2;
                      StrainDerivative->bx+=temp1;
                      StrainDerivative->by+=(posB.y-comB.y)*f1*dr.y;
                      StrainDerivative->bz+=temp3;
                      StrainDerivative->cx+=temp2;
                      StrainDerivative->cy+=temp3;
                      StrainDerivative->cz+=(posB.z-comB.z)*f1*dr.z;
                    }

                    // add contribution to the first derivatives
                    if(ComputeGradient)
                    {
                      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      {
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;
                      }

                      if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                      if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                      if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                        // add contribution to the first derivatives
                        if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                        if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                        if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,1.0);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,f1,f2,dr);
                      HessianAtomicStrainStrain(HessianMatrix,f1,f2,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianCenterOfMassStrainJ(HessianMatrix,index_i,index_j,f1,f2,dr,posB,comB);
                        HessianOrientationStrainJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianOrientationStrainJ_FR(HessianMatrix,index_j2,index2,f1,f2,posB,comB,dr);

                        HessianCorrectionStrainStrain(HessianMatrix,f1,f2,dr,posB,comB);
                      }
                    }
                  }
                }
              }
            }
            ncell++;
          }
    }
  }
}


void ComputeFrameworkIntraVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                     REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,j,typeA,typeB,f1,f2,start;
  REAL energy,DF,DDF;
  REAL rr;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j;
  int A,B;
  int ncell,k1,k2,k3,indj;
  REAL ReplicaFactor;
  REAL *parms;

  if(!InternalFrameworkLennardJonesInteractions) return;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if((f1==f2)&&(ncell==0)) start=i+1;
              else start=0;
              for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
              {
                indj=j+ncell*Framework[CurrentSystem].NumberOfAtoms[f2];
                if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][indj],0))
                {
                  typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
                  posB=Framework[CurrentSystem].Atoms[f2][j].Position;
                  index_j=Framework[CurrentSystem].Atoms[f2][j].HessianIndex;

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
                    if(ncell==0) ReplicaFactor=1.0;
                    else ReplicaFactor=0.5;

                    PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,1.0);

                    // add contribution to the energy
                    *Energy+=ReplicaFactor*energy;

                    if(ComputeGradient)
                    {
                      // add contribution to the first derivatives
                      if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
                      if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
                      if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

                      if(ncell==0)
                      {
                        if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
                        if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
                        if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;
                      }

                      GradientStrain(Gradient,ReplicaFactor*DF,dr);
                    }

                    // add contribution to the strain derivative tensor
                    StrainDerivative->ax+=ReplicaFactor*DF*dr.x*dr.x;
                    StrainDerivative->bx+=ReplicaFactor*DF*dr.y*dr.x;
                    StrainDerivative->cx+=ReplicaFactor*DF*dr.z*dr.x;

                    StrainDerivative->ay+=ReplicaFactor*DF*dr.x*dr.y;
                    StrainDerivative->by+=ReplicaFactor*DF*dr.y*dr.y;
                    StrainDerivative->cy+=ReplicaFactor*DF*dr.z*dr.y;

                    StrainDerivative->az+=ReplicaFactor*DF*dr.x*dr.z;
                    StrainDerivative->bz+=ReplicaFactor*DF*dr.y*dr.z;
                    StrainDerivative->cz+=ReplicaFactor*DF*dr.z*dr.z;

                    if(ComputeHessian)
                    {
                      // add contribution to the second derivatives (Hessian matrix)
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,ReplicaFactor);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,ReplicaFactor*DF,ReplicaFactor*DDF,dr);
                      HessianAtomicStrainStrain(HessianMatrix,ReplicaFactor*DF,ReplicaFactor*DDF,dr);
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
          A=MIN2(Framework[CurrentSystem].Torsions[f1][i].A,Framework[CurrentSystem].Torsions[f1][i].D);
          B=MAX2(Framework[CurrentSystem].Torsions[f1][i].A,Framework[CurrentSystem].Torsions[f1][i].D);

          typeA=Framework[CurrentSystem].Atoms[f1][A].Type;
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;

          typeB=Framework[CurrentSystem].Atoms[f1][B].Type;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,1.0);

          energy*=parms[6];
          DF*=parms[6];
          DDF*=parms[6];

          // add contribution to the energy
          *Energy+=energy;

          if(ComputeGradient)
          {
            // add contribution to the first derivatives
            if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
            if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
            if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

            if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
            if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
            if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

            GradientStrain(Gradient,DF,dr);
          }

          // add contribution to the strain derivative tensor
          StrainDerivative->ax+=DF*dr.x*dr.x;
          StrainDerivative->bx+=DF*dr.y*dr.x;
          StrainDerivative->cx+=DF*dr.z*dr.x;

          StrainDerivative->ay+=DF*dr.x*dr.y;
          StrainDerivative->by+=DF*dr.y*dr.y;
          StrainDerivative->cy+=DF*dr.z*dr.y;

          StrainDerivative->az+=DF*dr.x*dr.z;
          StrainDerivative->bz+=DF*dr.y*dr.z;
          StrainDerivative->cz+=DF*dr.z*dr.z;

          if(ComputeHessian)
          {
            // add contribution to the second derivatives (Hessian matrix)
            HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
            HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,DF,DDF,dr);
            HessianAtomicStrainStrain(HessianMatrix,DF,DDF,dr);
          }
        }
      }
    }
  }
}

void ComputeFrameworkIntraChargeChargeHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                              REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,j,typeA,typeB,f1,f2,start;
  REAL ChargeA,ChargeB,DF,DDF;
  REAL rr,r,U;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j;
  int A,B;
  int ncell,k1,k2,k3,indj;
  REAL ReplicaFactor;
  REAL *parms;

  if(ChargeMethod==NONE) return;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        ChargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;
        index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if((f1==f2)&&(ncell==0)) start=i+1;
              else start=0;
              for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
              {
                indj=j+ncell*Framework[CurrentSystem].NumberOfAtoms[f2];
                if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][indj],1))
                {
                  typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
                  posB=Framework[CurrentSystem].Atoms[f2][j].Position;
                  ChargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;
                  index_j=Framework[CurrentSystem].Atoms[f2][j].HessianIndex;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  {
                    if(ncell==0) ReplicaFactor=1.0;
                    else ReplicaFactor=0.5;

                    PotentialSecondDerivativeCoulombic(ChargeA,ChargeB,rr,&U,&DF,&DDF);

                    *Energy+=ReplicaFactor*U;

                    if(ComputeGradient)
                    {
                      // add contribution to the first derivatives
                      if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
                      if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
                      if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

                      if(ncell==0)
                      {
                        if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
                        if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
                        if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;
                      }

                      GradientStrain(Gradient,ReplicaFactor*DF,dr);
                    }

                    // add contribution to the strain derivative tensor
                    StrainDerivative->ax+=ReplicaFactor*DF*dr.x*dr.x;
                    StrainDerivative->bx+=ReplicaFactor*DF*dr.y*dr.x;
                    StrainDerivative->cx+=ReplicaFactor*DF*dr.z*dr.x;

                    StrainDerivative->ay+=ReplicaFactor*DF*dr.x*dr.y;
                    StrainDerivative->by+=ReplicaFactor*DF*dr.y*dr.y;
                    StrainDerivative->cy+=ReplicaFactor*DF*dr.z*dr.y;

                    StrainDerivative->az+=ReplicaFactor*DF*dr.x*dr.z;
                    StrainDerivative->bz+=ReplicaFactor*DF*dr.y*dr.z;
                    StrainDerivative->cz+=ReplicaFactor*DF*dr.z*dr.z;

                    if(ComputeHessian)
                    {
                      // add contribution to the second derivatives (Hessian matrix)
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,ReplicaFactor);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,ReplicaFactor*DF,ReplicaFactor*DDF,dr);
                      HessianAtomicStrainStrain(HessianMatrix,ReplicaFactor*DF,ReplicaFactor*DDF,dr);
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

        if(fabs(parms[7])>1e-8)
        {
          A=MIN2(Framework[CurrentSystem].Torsions[f1][i].A,Framework[CurrentSystem].Torsions[f1][i].D);
          B=MAX2(Framework[CurrentSystem].Torsions[f1][i].A,Framework[CurrentSystem].Torsions[f1][i].D);

          typeA=Framework[CurrentSystem].Atoms[f1][A].Type;
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          ChargeA=Framework[CurrentSystem].Atoms[f1][A].Charge;
          index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;

          typeB=Framework[CurrentSystem].Atoms[f1][B].Type;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          ChargeB=Framework[CurrentSystem].Atoms[f1][B].Charge;
          index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          r=sqrt(rr);


          // add contribution to the energy
          *Energy+=parms[7]*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r);
          DF=-parms[7]*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
          DDF=parms[7]*3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);

          //if((index_i<0)&&(index_j)<0) continue;

          if(ComputeGradient)
          {
            // add contribution to the first derivatives
            if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
            if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
            if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

            if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
            if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
            if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

            GradientStrain(Gradient,DF,dr);
          }

          // add contribution to the strain derivative tensor
          StrainDerivative->ax+=DF*dr.x*dr.x;
          StrainDerivative->bx+=DF*dr.y*dr.x;
          StrainDerivative->cx+=DF*dr.z*dr.x;

          StrainDerivative->ay+=DF*dr.x*dr.y;
          StrainDerivative->by+=DF*dr.y*dr.y;
          StrainDerivative->cy+=DF*dr.z*dr.y;

          StrainDerivative->az+=DF*dr.x*dr.z;
          StrainDerivative->bz+=DF*dr.y*dr.z;
          StrainDerivative->cz+=DF*dr.z*dr.z;

          if(ComputeHessian)
          {
            // add contribution to the second derivatives (Hessian matrix)
            HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
            HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,DF,DDF,dr);
            HessianAtomicStrainStrain(HessianMatrix,DF,DDF,dr);
          }
        }
      }
    }
  }
}


void ComputeFrameworkBondHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i;   // loop variable
  int f1;  // loop over all frameworks
  REAL r;  // distance
  REAL rr; // distance squared
  REAL temp,temp2; // temporary
  REAL exp_term;   // temporary
  REAL U;  // energy of a specific interaction
  REAL DF; // first derivative
  REAL DDF;  // second derivative
  VECTOR dr; // atoms separation vector
  int A,B;   // atom-indices
  INT_VECTOR3 index_i,index_j; // indices of the Hessian
  REAL *parms;  // pointer to potential parameter

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBonds[f1];i++)
    {
      A=Framework[CurrentSystem].Bonds[f1][i].A;
      B=Framework[CurrentSystem].Bonds[f1][i].B;

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;

      dr.x=Framework[CurrentSystem].Atoms[f1][A].Position.x-
           Framework[CurrentSystem].Atoms[f1][B].Position.x;
      dr.y=Framework[CurrentSystem].Atoms[f1][A].Position.y-
           Framework[CurrentSystem].Atoms[f1][B].Position.y;
      dr.z=Framework[CurrentSystem].Atoms[f1][A].Position.z-
           Framework[CurrentSystem].Atoms[f1][B].Position.z;

      // apply boundary condition
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
          DDF=(parms[0]*parms[1])/(r*rr);
          //DDF=(parms[0]/r-DF)/rr;
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
          fprintf(stderr, "Undefined Bond potential in routine 'CalculateFrameworkBondHessian' ('framework_hessian.c')\n");
          exit(0);
          break;
      }

      // add contribution to the energy
      *Energy+=U;

      //if((index_i<0)&&(index_j<0)) continue;

      if(ComputeGradient)
      {
        // add contribution to the first derivatives
        if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
        if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
        if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

        if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
        if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
        if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

        GradientStrain(Gradient,DF,dr);
      }

      // add contribution to the strain derivative tensor
      StrainDerivative->ax+=dr.x*DF*dr.x;
      StrainDerivative->bx+=dr.y*DF*dr.x;
      StrainDerivative->cx+=dr.z*DF*dr.x;

      StrainDerivative->ay+=dr.x*DF*dr.y;
      StrainDerivative->by+=dr.y*DF*dr.y;
      StrainDerivative->cy+=dr.z*DF*dr.y;

      StrainDerivative->az+=dr.x*DF*dr.z;
      StrainDerivative->bz+=dr.y*DF*dr.z;
      StrainDerivative->cz+=dr.z*DF*dr.z;

      if(ComputeHessian)
      {
        // add contribution to the second derivatives (Hessian matrix)
        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
        HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,DF,DDF,dr);
        HessianAtomicStrainStrain(HessianMatrix,DF,DDF,dr);
      }
    }
  }
}


void ComputeFrameworkUreyBradleyHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i;   // loop variable
  int f1;  // loop over all frameworks
  REAL r;  // distance
  REAL rr; // distance squared
  REAL temp,temp2; // temporary
  REAL exp_term;   // temporary
  REAL U;  // energy of a specific interaction
  REAL DF; // first derivative
  REAL DDF;  // second derivative
  VECTOR dr; // atoms separation vector
  int A,C;   // atom-indices
  INT_VECTOR3 index_i,index_j; // indices of the Hessian
  REAL *parms;  // pointer to potential parameter

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfUreyBradleys[f1];i++)
    {
      A=Framework[CurrentSystem].UreyBradleys[f1][i].A;
      C=Framework[CurrentSystem].UreyBradleys[f1][i].C;

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][C].HessianIndex;

      dr.x=Framework[CurrentSystem].Atoms[f1][A].Position.x-
           Framework[CurrentSystem].Atoms[f1][C].Position.x;
      dr.y=Framework[CurrentSystem].Atoms[f1][A].Position.y-
           Framework[CurrentSystem].Atoms[f1][C].Position.y;
      dr.z=Framework[CurrentSystem].Atoms[f1][A].Position.z-
           Framework[CurrentSystem].Atoms[f1][C].Position.z;

      // apply boundary condition
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
          DDF=(parms[0]*parms[1])/(r*rr);
          //DDF=(parms[0]/r-DF)/rr;
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
          fprintf(stderr, "Undefined UreyBradley potential in routine 'CalculateFrameworkUreyBradleyHessian' ('framework_hessian.c')\n");
          exit(0);
          break;
      }

      // add contribution to the energy
      *Energy+=U;

      //if((index_i<0)&&(index_j<0)) continue;

      if(ComputeGradient)
      {
        // add contribution to the first derivatives
        if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
        if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
        if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

        if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
        if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
        if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

        GradientStrain(Gradient,DF,dr);
      }

      // add contribution to the strain derivative tensor
      StrainDerivative->ax+=dr.x*DF*dr.x;
      StrainDerivative->bx+=dr.y*DF*dr.x;
      StrainDerivative->cx+=dr.z*DF*dr.x;

      StrainDerivative->ay+=dr.x*DF*dr.y;
      StrainDerivative->by+=dr.y*DF*dr.y;
      StrainDerivative->cy+=dr.z*DF*dr.y;

      StrainDerivative->az+=dr.x*DF*dr.z;
      StrainDerivative->bz+=dr.y*DF*dr.z;
      StrainDerivative->cz+=dr.z*DF*dr.z;

      if(ComputeHessian)
      {
        // add contribution to the second derivatives (Hessian matrix)
        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
        HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,DF,DDF,dr);
        HessianAtomicStrainStrain(HessianMatrix,DF,DDF,dr);
      }
    }
  }
}


static inline void GradientStrainBend(REAL *Gradient,REAL_MATRIX3x3 S)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]+=S.ax+S.by+S.cz;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n]+=S.ax;
          Gradient[n+1]+=S.by;
          Gradient[n+2]+=S.cz;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n]+=S.ax;
          Gradient[n+1]+=S.bx;
          Gradient[n+2]+=S.cx;
          Gradient[n+3]+=S.by;
          Gradient[n+4]+=S.cy;
          Gradient[n+5]+=S.cz;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.by;
              Gradient[n+2]+=S.cy;
              Gradient[n+3]+=S.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.cx;
              Gradient[n+2]+=S.by;
              Gradient[n+3]+=S.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.bx;
              Gradient[n+2]+=S.by;
              Gradient[n+3]+=S.cz;
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}



void ComputeFrameworkBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,f1,A,B,C;
  REAL *parms,U;
  REAL CosTheta,Theta,SinTheta,temp,temp2;
  REAL rab,rbc,rac;
  POINT posA,posB,posC;
  VECTOR Rab,Rbc,Rac;
  INT_VECTOR3 index_i,index_j,index_k;
  int index1,index2,index3;
  REAL DTDX,DF,DDF;
  REAL_MATRIX3x3 D2I,D2K,D2IK;
  VECTOR dtA,dtB,dtC;
  REAL_MATRIX3x3 S;
  VECTOR vec_u,vec_v;
  REAL u,v;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBends[f1];i++)
    {
      A=Framework[CurrentSystem].Bends[f1][i].A;
      B=Framework[CurrentSystem].Bends[f1][i].B;
      C=Framework[CurrentSystem].Bends[f1][i].C;
      parms=Framework[CurrentSystem].BendArguments[f1][i];

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;
      index_k=Framework[CurrentSystem].Atoms[f1][C].HessianIndex;

      index1=Framework[CurrentSystem].Atoms[f1][A].HessianAtomIndex;
      index2=Framework[CurrentSystem].Atoms[f1][B].HessianAtomIndex;
      index3=Framework[CurrentSystem].Atoms[f1][C].HessianAtomIndex;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;

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


      switch(Framework[CurrentSystem].BendType[f1][i])
      {
        case HARMONIC_BEND:
          // (1/2)p_0*(theta-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          temp=Theta-parms[1];
          U=0.5*parms[0]*SQR(temp);
          DF=parms[0]*temp*DTDX;
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
          fprintf(stderr, "Undefined Bend potential\n");
          exit(0);
          break;
      }

      // add contribution to the energy
      *Energy+=U;

      //if((index_i<0)&&(index_j<0)&&(index_k<0)) continue;

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

      S.ax=rab*Rab.x*DF*dtA.x+rbc*Rbc.x*DF*dtC.x;
      S.bx=rab*Rab.y*DF*dtA.x+rbc*Rbc.y*DF*dtC.x;
      S.cx=rab*Rab.z*DF*dtA.x+rbc*Rbc.z*DF*dtC.x;

      S.ay=rab*Rab.x*DF*dtA.y+rbc*Rbc.x*DF*dtC.y;
      S.by=rab*Rab.y*DF*dtA.y+rbc*Rbc.y*DF*dtC.y;
      S.cy=rab*Rab.z*DF*dtA.y+rbc*Rbc.z*DF*dtC.y;

      S.az=rab*Rab.x*DF*dtA.z+rbc*Rbc.x*DF*dtC.z;
      S.bz=rab*Rab.y*DF*dtA.z+rbc*Rbc.y*DF*dtC.z;
      S.cz=rab*Rab.z*DF*dtA.z+rbc*Rbc.z*DF*dtC.z;

      StrainDerivative->ax+=S.ax;
      StrainDerivative->bx+=S.bx;
      StrainDerivative->cx+=S.cx;

      StrainDerivative->ay+=S.ay;
      StrainDerivative->by+=S.by;
      StrainDerivative->cy+=S.cy;

      StrainDerivative->az+=S.az;
      StrainDerivative->bz+=S.bz;
      StrainDerivative->cz+=S.cz;

      // add contribution to the first derivatives
      if(ComputeGradient)
      {
        if(index_i.x>=0) Gradient[index_i.x]+=DF*dtA.x;
        if(index_i.y>=0) Gradient[index_i.y]+=DF*dtA.y;
        if(index_i.z>=0) Gradient[index_i.z]+=DF*dtA.z;

        if(index_j.x>=0) Gradient[index_j.x]+=DF*dtB.x;
        if(index_j.y>=0) Gradient[index_j.y]+=DF*dtB.y;
        if(index_j.z>=0) Gradient[index_j.z]+=DF*dtB.z;

        if(index_k.x>=0) Gradient[index_k.x]+=DF*dtC.x;
        if(index_k.y>=0) Gradient[index_k.y]+=DF*dtC.y;
        if(index_k.z>=0) Gradient[index_k.z]+=DF*dtC.z;

        GradientStrainBend(Gradient,S);
      }

      if(ComputeHessian)
      {
        // Calculate some diagonal Hessian terms for A
        D2I.ax=-DF*(2.0*dtA.x*Rab.x+CosTheta*(1.0-Rab.x*Rab.x)/rab)/rab;
        D2I.by=-DF*(2.0*dtA.y*Rab.y+CosTheta*(1.0-Rab.y*Rab.y)/rab)/rab;
        D2I.cz=-DF*(2.0*dtA.z*Rab.z+CosTheta*(1.0-Rab.z*Rab.z)/rab)/rab;
        D2I.ay=DF*(CosTheta*Rab.x*Rab.y/rab-dtA.x*Rab.y-dtA.y*Rab.x)/rab;
        D2I.az=DF*(CosTheta*Rab.x*Rab.z/rab-dtA.x*Rab.z-dtA.z*Rab.x)/rab;
        D2I.bz=DF*(CosTheta*Rab.y*Rab.z/rab-dtA.y*Rab.z-dtA.z*Rab.y)/rab;

        // Calculate some diagonal Hessian terms for C
        D2K.ax=-DF*(2.0*dtC.x*Rbc.x+CosTheta*(1.0-Rbc.x*Rbc.x)/rbc)/rbc;
        D2K.by=-DF*(2.0*dtC.y*Rbc.y+CosTheta*(1.0-Rbc.y*Rbc.y)/rbc)/rbc;
        D2K.cz=-DF*(2.0*dtC.z*Rbc.z+CosTheta*(1.0-Rbc.z*Rbc.z)/rbc)/rbc;
        D2K.ay=DF*(CosTheta*Rbc.x*Rbc.y/rbc-dtC.x*Rbc.y-dtC.y*Rbc.x)/rbc;
        D2K.az=DF*(CosTheta*Rbc.x*Rbc.z/rbc-dtC.x*Rbc.z-dtC.z*Rbc.x)/rbc;
        D2K.bz=DF*(CosTheta*Rbc.y*Rbc.z/rbc-dtC.y*Rbc.z-dtC.z*Rbc.y)/rbc;

        // Calculate the AC off-diagonal terms.
        D2IK.ax=DF*(CosTheta*Rab.x*Rbc.x-Rab.x*Rab.x-Rbc.x*Rbc.x+1.0)/(rab*rbc);
        D2IK.ay=DF*(CosTheta*Rab.x*Rbc.y-Rab.x*Rab.y-Rbc.x*Rbc.y)/(rab*rbc);
        D2IK.az=DF*(CosTheta*Rab.x*Rbc.z-Rab.x*Rab.z-Rbc.x*Rbc.z)/(rab*rbc);
        D2IK.bx=DF*(CosTheta*Rab.y*Rbc.x-Rab.y*Rab.x-Rbc.y*Rbc.x)/(rab*rbc);
        D2IK.by=DF*(CosTheta*Rab.y*Rbc.y-Rab.y*Rab.y-Rbc.y*Rbc.y+1.0)/(rab*rbc);
        D2IK.bz=DF*(CosTheta*Rab.y*Rbc.z-Rab.y*Rab.z-Rbc.y*Rbc.z)/(rab*rbc);
        D2IK.cx=DF*(CosTheta*Rab.z*Rbc.x-Rab.z*Rab.x-Rbc.z*Rbc.x)/(rab*rbc);
        D2IK.cy=DF*(CosTheta*Rab.z*Rbc.y-Rab.z*Rab.y-Rbc.z*Rbc.y)/(rab*rbc);
        D2IK.cz=DF*(CosTheta*Rab.z*Rbc.z-Rab.z*Rab.z-Rbc.z*Rbc.z+1.0)/(rab*rbc);

        // Calculate AA-block of the Hessian
        if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x]+=DDF*dtA.x*dtA.x+D2I.ax;
        if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y]+=DDF*dtA.x*dtA.y+D2I.ay;
        if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y]+=DDF*dtA.y*dtA.y+D2I.by;
        if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z]+=DDF*dtA.x*dtA.z+D2I.az;
        if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z]+=DDF*dtA.y*dtA.z+D2I.bz;
        if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z]+=DDF*dtA.z*dtA.z+D2I.cz;

        // Calculate BB-block of the Hessian
        if((index_j.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_j.x][index_j.x]+=DDF*dtB.x*dtB.x+(D2I.ax+D2K.ax+2.0*D2IK.ax);
        if((index_j.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.x][index_j.y]+=DDF*dtB.x*dtB.y+(D2I.ay+D2K.ay+D2IK.ay+D2IK.bx);
        if((index_j.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.y][index_j.y]+=DDF*dtB.y*dtB.y+(D2I.by+D2K.by+2.0*D2IK.by);
        if((index_j.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.x][index_j.z]+=DDF*dtB.x*dtB.z+(D2I.az+D2K.az+D2IK.az+D2IK.cx);
        if((index_j.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.y][index_j.z]+=DDF*dtB.y*dtB.z+(D2I.bz+D2K.bz+D2IK.bz+D2IK.cy);
        if((index_j.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.z][index_j.z]+=DDF*dtB.z*dtB.z+(D2I.cz+D2K.cz+2.0*D2IK.cz);

        // Calculate the CC-block of the Hessian
        if((index_k.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_k.x][index_k.x]+=DDF*dtC.x*dtC.x+D2K.ax;
        if((index_k.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_k.x][index_k.y]+=DDF*dtC.x*dtC.y+D2K.ay;
        if((index_k.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_k.y][index_k.y]+=DDF*dtC.y*dtC.y+D2K.by;
        if((index_k.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.x][index_k.z]+=DDF*dtC.x*dtC.z+D2K.az;
        if((index_k.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.y][index_k.z]+=DDF*dtC.y*dtC.z+D2K.bz;
        if((index_k.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.z][index_k.z]+=DDF*dtC.z*dtC.z+D2K.cz;

        // Calculate the AB-block of the Hessian
        if(index1<index2)
        {
          if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]+=DDF*dtA.x*dtB.x-D2I.ax-D2IK.ax;
          if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]+=DDF*dtA.x*dtB.y-D2I.ay-D2IK.ay;
          if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]+=DDF*dtA.x*dtB.z-D2I.az-D2IK.az;
          if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]+=DDF*dtA.y*dtB.x-D2I.ay-D2IK.bx;
          if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]+=DDF*dtA.y*dtB.y-D2I.by-D2IK.by;
          if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]+=DDF*dtA.y*dtB.z-D2I.bz-D2IK.bz;
          if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]+=DDF*dtA.z*dtB.x-D2I.az-D2IK.cx;
          if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]+=DDF*dtA.z*dtB.y-D2I.bz-D2IK.cy;
          if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]+=DDF*dtA.z*dtB.z-D2I.cz-D2IK.cz;
        }
        else
        {
          if((index_j.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.x][index_i.x]+=DDF*dtA.x*dtB.x-D2I.ax-D2IK.ax;
          if((index_j.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.x][index_i.y]+=DDF*dtA.y*dtB.x-D2I.ay-D2IK.bx;
          if((index_j.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.x][index_i.z]+=DDF*dtA.z*dtB.x-D2I.az-D2IK.cx;
          if((index_j.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.y][index_i.x]+=DDF*dtA.x*dtB.y-D2I.ay-D2IK.ay;
          if((index_j.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.y][index_i.y]+=DDF*dtA.y*dtB.y-D2I.by-D2IK.by;
          if((index_j.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.y][index_i.z]+=DDF*dtA.z*dtB.y-D2I.bz-D2IK.cy;
          if((index_j.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.z][index_i.x]+=DDF*dtA.x*dtB.z-D2I.az-D2IK.az;
          if((index_j.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.z][index_i.y]+=DDF*dtA.y*dtB.z-D2I.bz-D2IK.bz;
          if((index_j.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.z][index_i.z]+=DDF*dtA.z*dtB.z-D2I.cz-D2IK.cz;
        }

        // Calculate the AC-block of the Hessian
        if(index1<index3)
        {
          if((index_i.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.x][index_k.x]+=DDF*dtA.x*dtC.x+D2IK.ax;
          if((index_i.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.x][index_k.y]+=DDF*dtA.x*dtC.y+D2IK.ay;
          if((index_i.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.x][index_k.z]+=DDF*dtA.x*dtC.z+D2IK.az;
          if((index_i.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.y][index_k.x]+=DDF*dtA.y*dtC.x+D2IK.bx;
          if((index_i.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.y][index_k.y]+=DDF*dtA.y*dtC.y+D2IK.by;
          if((index_i.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.y][index_k.z]+=DDF*dtA.y*dtC.z+D2IK.bz;
          if((index_i.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.z][index_k.x]+=DDF*dtA.z*dtC.x+D2IK.cx;
          if((index_i.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.z][index_k.y]+=DDF*dtA.z*dtC.y+D2IK.cy;
          if((index_i.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.z][index_k.z]+=DDF*dtA.z*dtC.z+D2IK.cz;
        }
        else
        {
          if((index_k.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.x][index_i.x]+=DDF*dtA.x*dtC.x+D2IK.ax;
          if((index_k.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.x][index_i.y]+=DDF*dtA.y*dtC.x+D2IK.bx;
          if((index_k.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.x][index_i.z]+=DDF*dtA.z*dtC.x+D2IK.cx;
          if((index_k.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.y][index_i.x]+=DDF*dtA.x*dtC.y+D2IK.ay;
          if((index_k.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.y][index_i.y]+=DDF*dtA.y*dtC.y+D2IK.by;
          if((index_k.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.y][index_i.z]+=DDF*dtA.z*dtC.y+D2IK.cy;
          if((index_k.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.z][index_i.x]+=DDF*dtA.x*dtC.z+D2IK.az;
          if((index_k.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.z][index_i.y]+=DDF*dtA.y*dtC.z+D2IK.bz;
          if((index_k.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.z][index_i.z]+=DDF*dtA.z*dtC.z+D2IK.cz;
        }

        // Calculate the BC-block of the Hessian
        if(index3<index2)
        {
          if((index_k.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.x][index_j.x]+=DDF*dtC.x*dtB.x-D2K.ax-D2IK.ax;
          if((index_k.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.x][index_j.y]+=DDF*dtC.x*dtB.y-D2K.ay-D2IK.bx;
          if((index_k.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.x][index_j.z]+=DDF*dtC.x*dtB.z-D2K.az-D2IK.cx;
          if((index_k.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.y][index_j.x]+=DDF*dtC.y*dtB.x-D2K.ay-D2IK.ay;
          if((index_k.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.y][index_j.y]+=DDF*dtC.y*dtB.y-D2K.by-D2IK.by;
          if((index_k.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.y][index_j.z]+=DDF*dtC.y*dtB.z-D2K.bz-D2IK.cy;
          if((index_k.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.z][index_j.x]+=DDF*dtC.z*dtB.x-D2K.az-D2IK.az;
          if((index_k.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.z][index_j.y]+=DDF*dtC.z*dtB.y-D2K.bz-D2IK.bz;
          if((index_k.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.z][index_j.z]+=DDF*dtC.z*dtB.z-D2K.cz-D2IK.cz;
        }
        else
        {
          if((index_j.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.x][index_k.x]+=DDF*dtC.x*dtB.x-D2K.ax-D2IK.ax;
          if((index_j.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.x][index_k.y]+=DDF*dtC.y*dtB.x-D2K.ay-D2IK.ay;
          if((index_j.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.x][index_k.z]+=DDF*dtC.z*dtB.x-D2K.az-D2IK.az;
          if((index_j.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.y][index_k.x]+=DDF*dtC.x*dtB.y-D2K.ay-D2IK.bx;
          if((index_j.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.y][index_k.y]+=DDF*dtC.y*dtB.y-D2K.by-D2IK.by;
          if((index_j.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.y][index_k.z]+=DDF*dtC.z*dtB.y-D2K.bz-D2IK.bz;
          if((index_j.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.z][index_k.x]+=DDF*dtC.x*dtB.z-D2K.az-D2IK.cx;
          if((index_j.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.z][index_k.y]+=DDF*dtC.y*dtB.z-D2K.bz-D2IK.cy;
          if((index_j.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.z][index_k.z]+=DDF*dtC.z*dtB.z-D2K.cz-D2IK.cz;
        }

        HessianBendStrainPosition(index_i,index_j,index_k,HessianMatrix,vec_u,vec_v,u,v,
                                  rab,rbc,Rab,Rbc,dtA,dtC,DF,DDF,S,CosTheta);

        HessianBendStrainStrain(HessianMatrix,vec_u,vec_v,u,v,DF,DDF,S,CosTheta);
      }
    }
  }
}

static inline void GradientStrainTorsion(REAL *Gradient,REAL_MATRIX3x3 S)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]+=S.ax+S.by+S.cz;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n]+=S.ax;
          Gradient[n+1]+=S.by;
          Gradient[n+2]+=S.cz;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n]+=S.ax;
          Gradient[n+1]+=S.bx;
          Gradient[n+2]+=S.cx;
          Gradient[n+3]+=S.by;
          Gradient[n+4]+=S.cy;
          Gradient[n+5]+=S.cz;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.by;
              Gradient[n+2]+=S.cy;
              Gradient[n+3]+=S.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.cx;
              Gradient[n+2]+=S.by;
              Gradient[n+3]+=S.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.bx;
              Gradient[n+2]+=S.by;
              Gradient[n+3]+=S.cz;
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}



void ComputeFrameworkTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  REAL ShiftedCosPhi,ShiftedCosPhi2,ShiftedSinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL *parms;
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  int index1,index2,index3,index4;
  REAL_MATRIX3x3 D2I,D2L,D2IL,D2IJ,D2IK,D2JL;
  REAL_MATRIX3x3 D2KL,D2J,D2K,D2JK,S;
  VECTOR DIL,DDJ,DDK,DEJ,DEK;
  REAL DF,DDF;
  VECTOR vec_u,vec_v,vec_w;
  REAL u,w,v;
  VECTOR fa,fb,fc,fd;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
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

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;
      index_k=Framework[CurrentSystem].Atoms[f1][C].HessianIndex;
      index_l=Framework[CurrentSystem].Atoms[f1][D].HessianIndex;

      index1=Framework[CurrentSystem].Atoms[f1][A].HessianAtomIndex;
      index2=Framework[CurrentSystem].Atoms[f1][B].HessianAtomIndex;
      index3=Framework[CurrentSystem].Atoms[f1][C].HessianAtomIndex;
      index4=Framework[CurrentSystem].Atoms[f1][D].HessianAtomIndex;

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
          fprintf(stderr, "Undefined Torsion potential in 'framework_hessian.c'\n");
          U=DF=DDF=0.0;
          exit(0);
          break;
      }

      *Energy+=U;

      //if((index_i<0)&&(index_j<0)&&(index_k<0)&&(index_l<0)) continue;

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
      fa.x=DF*dtA.x;
      fa.y=DF*dtA.y;
      fa.z=DF*dtA.z;

      fb.x=DF*dtB.x;
      fb.y=DF*dtB.y;
      fb.z=DF*dtB.z;

      fc.x=DF*dtC.x;
      fc.y=DF*dtC.y;
      fc.z=DF*dtC.z;

      fd.x=DF*dtD.x;
      fd.y=DF*dtD.y;
      fd.z=DF*dtD.z;

      // add contribution to the strain derivative tensor
      // Note: rbc is here because the vector was normalized before
      S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
      S.bx=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
      S.cx=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

      S.ay=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
      S.by=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
      S.cy=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

      S.az=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
      S.bz=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
      S.cz=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;

      StrainDerivative->ax+=S.ax;
      StrainDerivative->bx+=S.bx;
      StrainDerivative->cx+=S.cx;

      StrainDerivative->ay+=S.ay;
      StrainDerivative->by+=S.by;
      StrainDerivative->cy+=S.cy;

      StrainDerivative->az+=S.az;
      StrainDerivative->bz+=S.bz;
      StrainDerivative->cz+=S.cz;

      if(ComputeGradient)
      {
        if(index_i.x>=0) Gradient[index_i.x]+=fa.x;
        if(index_i.y>=0) Gradient[index_i.y]+=fa.y;
        if(index_i.z>=0) Gradient[index_i.z]+=fa.z;

        if(index_j.x>=0) Gradient[index_j.x]+=fb.x;
        if(index_j.y>=0) Gradient[index_j.y]+=fb.y;
        if(index_j.z>=0) Gradient[index_j.z]+=fb.z;

        if(index_k.x>=0) Gradient[index_k.x]+=fc.x;
        if(index_k.y>=0) Gradient[index_k.y]+=fc.y;
        if(index_k.z>=0) Gradient[index_k.z]+=fc.z;

        if(index_l.x>=0) Gradient[index_l.x]+=fd.x;
        if(index_l.y>=0) Gradient[index_l.y]+=fd.y;
        if(index_l.z>=0) Gradient[index_l.z]+=fd.z;

        GradientStrainTorsion(Gradient,S);
      }

      if(ComputeHessian)
      {
        // Calculate the derivatives of DOTIJ and DOTLK.
        DIL.x=DF*Dcb.x/rbc;
        DIL.y=DF*Dcb.y/rbc;
        DIL.z=DF*Dcb.z/rbc;
        DDJ.x=DF*((2.0*d-1.0)*Dcb.x-Dab.x/rbc)/rbc;
        DDJ.y=DF*((2.0*d-1.0)*Dcb.y-Dab.y/rbc)/rbc;
        DDJ.z=DF*((2.0*d-1.0)*Dcb.z-Dab.z/rbc)/rbc;
        DDK.x=-DF*(2.0*d*Dcb.x-Dab.x/rbc)/rbc;
        DDK.y=-DF*(2.0*d*Dcb.y-Dab.y/rbc)/rbc;
        DDK.z=-DF*(2.0*d*Dcb.z-Dab.z/rbc)/rbc;
        DEJ.x=DF*(2.0*e*Dcb.x-Ddc.x/rbc)/rbc;
        DEJ.y=DF*(2.0*e*Dcb.y-Ddc.y/rbc)/rbc;
        DEJ.z=DF*(2.0*e*Dcb.z-Ddc.z/rbc)/rbc;
        DEK.x=-DF*((2.0*e+1.0)*Dcb.x-Ddc.x/rbc)/rbc;
        DEK.y=-DF*((2.0*e+1.0)*Dcb.y-Ddc.y/rbc)/rbc;
        DEK.z=-DF*((2.0*e+1.0)*Dcb.z-Ddc.z/rbc)/rbc;

        // Calculate some diagonal Hessian terms for I.
        D2I.ax=DF*(CosPhi*(dr.x*dr.x+Dcb.x*Dcb.x-1.0)/r-2.0*dr.x*dtA.x)/r;
        D2I.by=DF*(CosPhi*(dr.y*dr.y+Dcb.y*Dcb.y-1.0)/r-2.0*dr.y*dtA.y)/r;
        D2I.cz=DF*(CosPhi*(dr.z*dr.z+Dcb.z*Dcb.z-1.0)/r-2.0*dr.z*dtA.z)/r;
        D2I.ay=DF*(CosPhi*(dr.x*dr.y+Dcb.x*Dcb.y)/r-dr.x*dtA.y-dr.y*dtA.x)/r;
        D2I.az=DF*(CosPhi*(dr.x*dr.z+Dcb.x*Dcb.z)/r-dr.x*dtA.z-dr.z*dtA.x)/r;
        D2I.bz=DF*(CosPhi*(dr.y*dr.z+Dcb.y*Dcb.z)/r-dr.y*dtA.z-dr.z*dtA.y)/r;

        // Calculate some diagonal Hessian terms for L.
        D2L.ax=DF*(CosPhi*(ds.x*ds.x+Dcb.x*Dcb.x-1.0)/s-2.0*ds.x*dtD.x)/s;
        D2L.by=DF*(CosPhi*(ds.y*ds.y+Dcb.y*Dcb.y-1.0)/s-2.0*ds.y*dtD.y)/s;
        D2L.cz=DF*(CosPhi*(ds.z*ds.z+Dcb.z*Dcb.z-1.0)/s-2.0*ds.z*dtD.z)/s;
        D2L.ay=DF*(CosPhi*(ds.x*ds.y+Dcb.x*Dcb.y)/s-ds.x*dtD.y-ds.y*dtD.x)/s;
        D2L.az=DF*(CosPhi*(ds.x*ds.z+Dcb.x*Dcb.z)/s-ds.x*dtD.z-ds.z*dtD.x)/s;
        D2L.bz=DF*(CosPhi*(ds.y*ds.z+Dcb.y*Dcb.z)/s-ds.y*dtD.z-ds.z*dtD.y)/s;

        // Calculate the IL off-diagonal terms.
        D2IL.ax=DF*(CosPhi*dr.x*ds.x-dr.x*dr.x-ds.x*ds.x-Dcb.x*Dcb.x+1.0)/(r*s);
        D2IL.ay=DF*(CosPhi*dr.x*ds.y-dr.x*dr.y-ds.x*ds.y-Dcb.x*Dcb.y)/(r*s);
        D2IL.az=DF*(CosPhi*dr.x*ds.z-dr.x*dr.z-ds.x*ds.z-Dcb.x*Dcb.z)/(r*s);
        D2IL.bx=DF*(CosPhi*dr.y*ds.x-dr.y*dr.x-ds.y*ds.x-Dcb.y*Dcb.x)/(r*s);
        D2IL.by=DF*(CosPhi*dr.y*ds.y-dr.y*dr.y-ds.y*ds.y-Dcb.y*Dcb.y+1.0)/(r*s);
        D2IL.bz=DF*(CosPhi*dr.y*ds.z-dr.y*dr.z-ds.y*ds.z-Dcb.y*Dcb.z)/(r*s);
        D2IL.cx=DF*(CosPhi*dr.z*ds.x-dr.z*dr.x-ds.z*ds.x-Dcb.z*Dcb.x)/(r*s);
        D2IL.cy=DF*(CosPhi*dr.z*ds.y-dr.z*dr.y-ds.z*ds.y-Dcb.z*Dcb.y)/(r*s);
        D2IL.cz=DF*(CosPhi*dr.z*ds.z-dr.z*dr.z-ds.z*ds.z-Dcb.z*Dcb.z+1.0)/(r*s);

        // Calculate the IJ off-diagonal terms.
        D2IJ.ax=D2I.ax*(d-1.0)+D2IL.ax*e+DIL.x*dtA.x;
        D2IJ.ay=D2I.ay*(d-1.0)+D2IL.ay*e+DIL.x*dtA.y;
        D2IJ.az=D2I.az*(d-1.0)+D2IL.az*e+DIL.x*dtA.z;
        D2IJ.bx=D2I.ay*(d-1.0)+D2IL.bx*e+DIL.y*dtA.x;
        D2IJ.by=D2I.by*(d-1.0)+D2IL.by*e+DIL.y*dtA.y;
        D2IJ.bz=D2I.bz*(d-1.0)+D2IL.bz*e+DIL.y*dtA.z;
        D2IJ.cx=D2I.az*(d-1.0)+D2IL.cx*e+DIL.z*dtA.x;
        D2IJ.cy=D2I.bz*(d-1.0)+D2IL.cy*e+DIL.z*dtA.y;
        D2IJ.cz=D2I.cz*(d-1.0)+D2IL.cz*e+DIL.z*dtA.z;

        // Calculate the IK off-diagonal terms.
        D2IK.ax=-D2I.ax*d-D2IL.ax*(e+1.0)-DIL.x*dtA.x;
        D2IK.ay=-D2I.ay*d-D2IL.ay*(e+1.0)-DIL.x*dtA.y;
        D2IK.az=-D2I.az*d-D2IL.az*(e+1.0)-DIL.x*dtA.z;
        D2IK.bx=-D2I.ay*d-D2IL.bx*(e+1.0)-DIL.y*dtA.x;
        D2IK.by=-D2I.by*d-D2IL.by*(e+1.0)-DIL.y*dtA.y;
        D2IK.bz=-D2I.bz*d-D2IL.bz*(e+1.0)-DIL.y*dtA.z;
        D2IK.cx=-D2I.az*d-D2IL.cx*(e+1.0)-DIL.z*dtA.x;
        D2IK.cy=-D2I.bz*d-D2IL.cy*(e+1.0)-DIL.z*dtA.y;
        D2IK.cz=-D2I.cz*d-D2IL.cz*(e+1.0)-DIL.z*dtA.z;

        // Calculate the JL off-diagonal terms.
        D2JL.ax=D2IL.ax*(d-1.0)+D2L.ax*e+DIL.x*dtD.x;
        D2JL.ay=D2IL.ay*(d-1.0)+D2L.ay*e+DIL.y*dtD.x;
        D2JL.az=D2IL.az*(d-1.0)+D2L.az*e+DIL.z*dtD.x;
        D2JL.bx=D2IL.bx*(d-1.0)+D2L.ay*e+DIL.x*dtD.y;
        D2JL.by=D2IL.by*(d-1.0)+D2L.by*e+DIL.y*dtD.y;
        D2JL.bz=D2IL.bz*(d-1.0)+D2L.bz*e+DIL.z*dtD.y;
        D2JL.cx=D2IL.cx*(d-1.0)+D2L.az*e+DIL.x*dtD.z;
        D2JL.cy=D2IL.cy*(d-1.0)+D2L.bz*e+DIL.y*dtD.z;
        D2JL.cz=D2IL.cz*(d-1.0)+D2L.cz*e+DIL.z*dtD.z;

        // Calculate the KL off-diagonal terms.
        D2KL.ax=-D2IL.ax*d-D2L.ax*(e+1.0)-DIL.x*dtD.x;
        D2KL.ay=-D2IL.ay*d-D2L.ay*(e+1.0)-DIL.y*dtD.x;
        D2KL.az=-D2IL.az*d-D2L.az*(e+1.0)-DIL.z*dtD.x;
        D2KL.bx=-D2IL.bx*d-D2L.ay*(e+1.0)-DIL.x*dtD.y;
        D2KL.by=-D2IL.by*d-D2L.by*(e+1.0)-DIL.y*dtD.y;
        D2KL.bz=-D2IL.bz*d-D2L.bz*(e+1.0)-DIL.z*dtD.y;
        D2KL.cx=-D2IL.cx*d-D2L.az*(e+1.0)-DIL.x*dtD.z;
        D2KL.cy=-D2IL.cy*d-D2L.bz*(e+1.0)-DIL.y*dtD.z;
        D2KL.cz=-D2IL.cz*d-D2L.cz*(e+1.0)-DIL.z*dtD.z;

        // Calculate the JJ diagonal terms.
        D2J.ax=D2IJ.ax*(d-1.0)+D2JL.ax*e+DDJ.x*dtA.x+DEJ.x*dtD.x;
        D2J.by=D2IJ.by*(d-1.0)+D2JL.by*e+DDJ.y*dtA.y+DEJ.y*dtD.y;
        D2J.cz=D2IJ.cz*(d-1.0)+D2JL.cz*e+DDJ.z*dtA.z+DEJ.z*dtD.z;
        D2J.ay=D2IJ.ay*(d-1.0)+D2JL.bx*e+DDJ.y*dtA.x+DEJ.y*dtD.x;
        D2J.az=D2IJ.az*(d-1.0)+D2JL.cx*e+DDJ.z*dtA.x+DEJ.z*dtD.x;
        D2J.bz=D2IJ.bz*(d-1.0)+D2JL.cy*e+DDJ.z*dtA.y+DEJ.z*dtD.y;

        // Calculate the KK diagonal terms.
        D2K.ax=-D2KL.ax*(e+1.0)-D2IK.ax*d-DDK.x*dtA.x-DEK.x*dtD.x;
        D2K.by=-D2KL.by*(e+1.0)-D2IK.by*d-DDK.y*dtA.y-DEK.y*dtD.y;
        D2K.cz=-D2KL.cz*(e+1.0)-D2IK.cz*d-DDK.z*dtA.z-DEK.z*dtD.z;
        D2K.ay=-D2KL.ay*(e+1.0)-D2IK.bx*d-DDK.x*dtA.y-DEK.x*dtD.y;
        D2K.az=-D2KL.az*(e+1.0)-D2IK.cx*d-DDK.x*dtA.z-DEK.x*dtD.z;
        D2K.bz=-D2KL.bz*(e+1.0)-D2IK.cy*d-DDK.y*dtA.z-DEK.y*dtD.z;

        // Calculate the JK off-diagonal terms.
        D2JK.ax=-D2IJ.ax*d-D2JL.ax*(e+1.0)-DDJ.x*dtA.x-DEJ.x*dtD.x;
        D2JK.ay=-D2IJ.bx*d-D2JL.ay*(e+1.0)-DDJ.x*dtA.y-DEJ.x*dtD.y;
        D2JK.az=-D2IJ.cx*d-D2JL.az*(e+1.0)-DDJ.x*dtA.z-DEJ.x*dtD.z;
        D2JK.bx=-D2IJ.ay*d-D2JL.bx*(e+1.0)-DDJ.y*dtA.x-DEJ.y*dtD.x;
        D2JK.by=-D2IJ.by*d-D2JL.by*(e+1.0)-DDJ.y*dtA.y-DEJ.y*dtD.y;
        D2JK.bz=-D2IJ.cy*d-D2JL.bz*(e+1.0)-DDJ.y*dtA.z-DEJ.y*dtD.z;
        D2JK.cx=-D2IJ.az*d-D2JL.cx*(e+1.0)-DDJ.z*dtA.x-DEJ.z*dtD.x;
        D2JK.cy=-D2IJ.bz*d-D2JL.cy*(e+1.0)-DDJ.z*dtA.y-DEJ.z*dtD.y;
        D2JK.cz=-D2IJ.cz*d-D2JL.cz*(e+1.0)-DDJ.z*dtA.z-DEJ.z*dtD.z;

        // Calculate the diagonal blocks of the hessian.
        if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x]+=DDF*dtA.x*dtA.x+D2I.ax;
        if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y]+=DDF*dtA.x*dtA.y+D2I.ay;
        if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y]+=DDF*dtA.y*dtA.y+D2I.by;
        if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z]+=DDF*dtA.x*dtA.z+D2I.az;
        if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z]+=DDF*dtA.y*dtA.z+D2I.bz;
        if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z]+=DDF*dtA.z*dtA.z+D2I.cz;

        if((index_j.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_j.x][index_j.x]+=DDF*dtB.x*dtB.x+D2J.ax;
        if((index_j.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.x][index_j.y]+=DDF*dtB.x*dtB.y+D2J.ay;
        if((index_j.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.y][index_j.y]+=DDF*dtB.y*dtB.y+D2J.by;
        if((index_j.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.x][index_j.z]+=DDF*dtB.x*dtB.z+D2J.az;
        if((index_j.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.y][index_j.z]+=DDF*dtB.y*dtB.z+D2J.bz;
        if((index_j.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.z][index_j.z]+=DDF*dtB.z*dtB.z+D2J.cz;

        if((index_k.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_k.x][index_k.x]+=DDF*dtC.x*dtC.x+D2K.ax;
        if((index_k.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_k.x][index_k.y]+=DDF*dtC.x*dtC.y+D2K.ay;
        if((index_k.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_k.y][index_k.y]+=DDF*dtC.y*dtC.y+D2K.by;
        if((index_k.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.x][index_k.z]+=DDF*dtC.x*dtC.z+D2K.az;
        if((index_k.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.y][index_k.z]+=DDF*dtC.y*dtC.z+D2K.bz;
        if((index_k.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.z][index_k.z]+=DDF*dtC.z*dtC.z+D2K.cz;

        if((index_l.x>=0)&&(index_l.x>=0)) HessianMatrix.element[index_l.x][index_l.x]+=DDF*dtD.x*dtD.x+D2L.ax;
        if((index_l.x>=0)&&(index_l.y>=0)) HessianMatrix.element[index_l.x][index_l.y]+=DDF*dtD.x*dtD.y+D2L.ay;
        if((index_l.y>=0)&&(index_l.y>=0)) HessianMatrix.element[index_l.y][index_l.y]+=DDF*dtD.y*dtD.y+D2L.by;
        if((index_l.x>=0)&&(index_l.z>=0)) HessianMatrix.element[index_l.x][index_l.z]+=DDF*dtD.x*dtD.z+D2L.az;
        if((index_l.y>=0)&&(index_l.z>=0)) HessianMatrix.element[index_l.y][index_l.z]+=DDF*dtD.y*dtD.z+D2L.bz;
        if((index_l.z>=0)&&(index_l.z>=0)) HessianMatrix.element[index_l.z][index_l.z]+=DDF*dtD.z*dtD.z+D2L.cz;

        if(index1<index2)
        {
          if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]+=DDF*dtA.x*dtB.x+D2IJ.ax;
          if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]+=DDF*dtA.x*dtB.y+D2IJ.ay;
          if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]+=DDF*dtA.x*dtB.z+D2IJ.az;
          if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]+=DDF*dtA.y*dtB.x+D2IJ.bx;
          if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]+=DDF*dtA.y*dtB.y+D2IJ.by;
          if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]+=DDF*dtA.y*dtB.z+D2IJ.bz;
          if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]+=DDF*dtA.z*dtB.x+D2IJ.cx;
          if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]+=DDF*dtA.z*dtB.y+D2IJ.cy;
          if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]+=DDF*dtA.z*dtB.z+D2IJ.cz;
        }
        else
        {
          if((index_j.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.x][index_i.x]+=DDF*dtA.x*dtB.x+D2IJ.ax;
          if((index_j.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.x][index_i.y]+=DDF*dtA.y*dtB.x+D2IJ.bx;
          if((index_j.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.x][index_i.z]+=DDF*dtA.z*dtB.x+D2IJ.cx;
          if((index_j.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.y][index_i.x]+=DDF*dtA.x*dtB.y+D2IJ.ay;
          if((index_j.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.y][index_i.y]+=DDF*dtA.y*dtB.y+D2IJ.by;
          if((index_j.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.y][index_i.z]+=DDF*dtA.z*dtB.y+D2IJ.cy;
          if((index_j.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.z][index_i.x]+=DDF*dtA.x*dtB.z+D2IJ.az;
          if((index_j.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.z][index_i.y]+=DDF*dtA.y*dtB.z+D2IJ.bz;
          if((index_j.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.z][index_i.z]+=DDF*dtA.z*dtB.z+D2IJ.cz;
        }

        if(index1<index3)
        {
          if((index_i.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.x][index_k.x]+=DDF*dtA.x*dtC.x+D2IK.ax;
          if((index_i.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.x][index_k.y]+=DDF*dtA.x*dtC.y+D2IK.ay;
          if((index_i.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.x][index_k.z]+=DDF*dtA.x*dtC.z+D2IK.az;
          if((index_i.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.y][index_k.x]+=DDF*dtA.y*dtC.x+D2IK.bx;
          if((index_i.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.y][index_k.y]+=DDF*dtA.y*dtC.y+D2IK.by;
          if((index_i.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.y][index_k.z]+=DDF*dtA.y*dtC.z+D2IK.bz;
          if((index_i.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.z][index_k.x]+=DDF*dtA.z*dtC.x+D2IK.cx;
          if((index_i.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.z][index_k.y]+=DDF*dtA.z*dtC.y+D2IK.cy;
          if((index_i.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.z][index_k.z]+=DDF*dtA.z*dtC.z+D2IK.cz;
        }
        else
        {
          if((index_k.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.x][index_i.x]+=DDF*dtA.x*dtC.x+D2IK.ax;
          if((index_k.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.x][index_i.y]+=DDF*dtA.y*dtC.x+D2IK.bx;
          if((index_k.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.x][index_i.z]+=DDF*dtA.z*dtC.x+D2IK.cx;
          if((index_k.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.y][index_i.x]+=DDF*dtA.x*dtC.y+D2IK.ay;
          if((index_k.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.y][index_i.y]+=DDF*dtA.y*dtC.y+D2IK.by;
          if((index_k.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.y][index_i.z]+=DDF*dtA.z*dtC.y+D2IK.cy;
          if((index_k.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.z][index_i.x]+=DDF*dtA.x*dtC.z+D2IK.az;
          if((index_k.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.z][index_i.y]+=DDF*dtA.y*dtC.z+D2IK.bz;
          if((index_k.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.z][index_i.z]+=DDF*dtA.z*dtC.z+D2IK.cz;
        }

        if(index1<index4)
        {
          if((index_i.x>=0)&&(index_l.x>=0)) HessianMatrix.element[index_i.x][index_l.x]+=DDF*dtA.x*dtD.x+D2IL.ax;
          if((index_i.x>=0)&&(index_l.y>=0)) HessianMatrix.element[index_i.x][index_l.y]+=DDF*dtA.x*dtD.y+D2IL.ay;
          if((index_i.x>=0)&&(index_l.z>=0)) HessianMatrix.element[index_i.x][index_l.z]+=DDF*dtA.x*dtD.z+D2IL.az;
          if((index_i.y>=0)&&(index_l.x>=0)) HessianMatrix.element[index_i.y][index_l.x]+=DDF*dtA.y*dtD.x+D2IL.bx;
          if((index_i.y>=0)&&(index_l.y>=0)) HessianMatrix.element[index_i.y][index_l.y]+=DDF*dtA.y*dtD.y+D2IL.by;
          if((index_i.y>=0)&&(index_l.z>=0)) HessianMatrix.element[index_i.y][index_l.z]+=DDF*dtA.y*dtD.z+D2IL.bz;
          if((index_i.z>=0)&&(index_l.x>=0)) HessianMatrix.element[index_i.z][index_l.x]+=DDF*dtA.z*dtD.x+D2IL.cx;
          if((index_i.z>=0)&&(index_l.y>=0)) HessianMatrix.element[index_i.z][index_l.y]+=DDF*dtA.z*dtD.y+D2IL.cy;
          if((index_i.z>=0)&&(index_l.z>=0)) HessianMatrix.element[index_i.z][index_l.z]+=DDF*dtA.z*dtD.z+D2IL.cz;
        }
        else
        {
          if((index_l.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_l.x][index_i.x]+=DDF*dtA.x*dtD.x+D2IL.ax;
          if((index_l.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_l.x][index_i.y]+=DDF*dtA.y*dtD.x+D2IL.bx;
          if((index_l.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_l.x][index_i.z]+=DDF*dtA.z*dtD.x+D2IL.cx;
          if((index_l.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_l.y][index_i.x]+=DDF*dtA.x*dtD.y+D2IL.ay;
          if((index_l.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_l.y][index_i.y]+=DDF*dtA.y*dtD.y+D2IL.by;
          if((index_l.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_l.y][index_i.z]+=DDF*dtA.z*dtD.y+D2IL.cy;
          if((index_l.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_l.z][index_i.x]+=DDF*dtA.x*dtD.z+D2IL.az;
          if((index_l.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_l.z][index_i.y]+=DDF*dtA.y*dtD.z+D2IL.bz;
          if((index_l.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_l.z][index_i.z]+=DDF*dtA.z*dtD.z+D2IL.cz;
        }

        if(index2<index3)
        {
          if((index_j.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.x][index_k.x]+=DDF*dtB.x*dtC.x+D2JK.ax;
          if((index_j.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.x][index_k.y]+=DDF*dtB.x*dtC.y+D2JK.ay;
          if((index_j.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.x][index_k.z]+=DDF*dtB.x*dtC.z+D2JK.az;
          if((index_j.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.y][index_k.x]+=DDF*dtB.y*dtC.x+D2JK.bx;
          if((index_j.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.y][index_k.y]+=DDF*dtB.y*dtC.y+D2JK.by;
          if((index_j.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.y][index_k.z]+=DDF*dtB.y*dtC.z+D2JK.bz;
          if((index_j.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.z][index_k.x]+=DDF*dtB.z*dtC.x+D2JK.cx;
          if((index_j.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.z][index_k.y]+=DDF*dtB.z*dtC.y+D2JK.cy;
          if((index_j.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.z][index_k.z]+=DDF*dtB.z*dtC.z+D2JK.cz;
        }
        else
        {
          if((index_k.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.x][index_j.x]+=DDF*dtB.x*dtC.x+D2JK.ax;
          if((index_k.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.x][index_j.y]+=DDF*dtB.y*dtC.x+D2JK.bx;
          if((index_k.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.x][index_j.z]+=DDF*dtB.z*dtC.x+D2JK.cx;
          if((index_k.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.y][index_j.x]+=DDF*dtB.x*dtC.y+D2JK.ay;
          if((index_k.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.y][index_j.y]+=DDF*dtB.y*dtC.y+D2JK.by;
          if((index_k.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.y][index_j.z]+=DDF*dtB.z*dtC.y+D2JK.cy;
          if((index_k.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.z][index_j.x]+=DDF*dtB.x*dtC.z+D2JK.az;
          if((index_k.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.z][index_j.y]+=DDF*dtB.y*dtC.z+D2JK.bz;
          if((index_k.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.z][index_j.z]+=DDF*dtB.z*dtC.z+D2JK.cz;
        }

        if(index2<index4)
        {
          if((index_j.x>=0)&&(index_l.x>=0)) HessianMatrix.element[index_j.x][index_l.x]+=DDF*dtB.x*dtD.x+D2JL.ax;
          if((index_j.x>=0)&&(index_l.y>=0)) HessianMatrix.element[index_j.x][index_l.y]+=DDF*dtB.x*dtD.y+D2JL.ay;
          if((index_j.x>=0)&&(index_l.z>=0)) HessianMatrix.element[index_j.x][index_l.z]+=DDF*dtB.x*dtD.z+D2JL.az;
          if((index_j.y>=0)&&(index_l.x>=0)) HessianMatrix.element[index_j.y][index_l.x]+=DDF*dtB.y*dtD.x+D2JL.bx;
          if((index_j.y>=0)&&(index_l.y>=0)) HessianMatrix.element[index_j.y][index_l.y]+=DDF*dtB.y*dtD.y+D2JL.by;
          if((index_j.y>=0)&&(index_l.z>=0)) HessianMatrix.element[index_j.y][index_l.z]+=DDF*dtB.y*dtD.z+D2JL.bz;
          if((index_j.z>=0)&&(index_l.x>=0)) HessianMatrix.element[index_j.z][index_l.x]+=DDF*dtB.z*dtD.x+D2JL.cx;
          if((index_j.z>=0)&&(index_l.y>=0)) HessianMatrix.element[index_j.z][index_l.y]+=DDF*dtB.z*dtD.y+D2JL.cy;
          if((index_j.z>=0)&&(index_l.z>=0)) HessianMatrix.element[index_j.z][index_l.z]+=DDF*dtB.z*dtD.z+D2JL.cz;
        }
        else
        {
          if((index_l.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_l.x][index_j.x]+=DDF*dtB.x*dtD.x+D2JL.ax;
          if((index_l.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_l.x][index_j.y]+=DDF*dtB.y*dtD.x+D2JL.bx;
          if((index_l.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_l.x][index_j.z]+=DDF*dtB.z*dtD.x+D2JL.cx;
          if((index_l.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_l.y][index_j.x]+=DDF*dtB.x*dtD.y+D2JL.ay;
          if((index_l.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_l.y][index_j.y]+=DDF*dtB.y*dtD.y+D2JL.by;
          if((index_l.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_l.y][index_j.z]+=DDF*dtB.z*dtD.y+D2JL.cy;
          if((index_l.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_l.z][index_j.x]+=DDF*dtB.x*dtD.z+D2JL.az;
          if((index_l.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_l.z][index_j.y]+=DDF*dtB.y*dtD.z+D2JL.bz;
          if((index_l.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_l.z][index_j.z]+=DDF*dtB.z*dtD.z+D2JL.cz;
        }

        if(index3<index4)
        {
          if((index_k.x>=0)&&(index_l.x>=0)) HessianMatrix.element[index_k.x][index_l.x]+=DDF*dtC.x*dtD.x+D2KL.ax;
          if((index_k.x>=0)&&(index_l.y>=0)) HessianMatrix.element[index_k.x][index_l.y]+=DDF*dtC.x*dtD.y+D2KL.ay;
          if((index_k.x>=0)&&(index_l.z>=0)) HessianMatrix.element[index_k.x][index_l.z]+=DDF*dtC.x*dtD.z+D2KL.az;
          if((index_k.y>=0)&&(index_l.x>=0)) HessianMatrix.element[index_k.y][index_l.x]+=DDF*dtC.y*dtD.x+D2KL.bx;
          if((index_k.y>=0)&&(index_l.y>=0)) HessianMatrix.element[index_k.y][index_l.y]+=DDF*dtC.y*dtD.y+D2KL.by;
          if((index_k.y>=0)&&(index_l.z>=0)) HessianMatrix.element[index_k.y][index_l.z]+=DDF*dtC.y*dtD.z+D2KL.bz;
          if((index_k.z>=0)&&(index_l.x>=0)) HessianMatrix.element[index_k.z][index_l.x]+=DDF*dtC.z*dtD.x+D2KL.cx;
          if((index_k.z>=0)&&(index_l.y>=0)) HessianMatrix.element[index_k.z][index_l.y]+=DDF*dtC.z*dtD.y+D2KL.cy;
          if((index_k.z>=0)&&(index_l.z>=0)) HessianMatrix.element[index_k.z][index_l.z]+=DDF*dtC.z*dtD.z+D2KL.cz;
        }
        else
        {
          if((index_l.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_l.x][index_k.x]+=DDF*dtC.x*dtD.x+D2KL.ax;
          if((index_l.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_l.x][index_k.y]+=DDF*dtC.y*dtD.x+D2KL.bx;
          if((index_l.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_l.x][index_k.z]+=DDF*dtC.z*dtD.x+D2KL.cx;
          if((index_l.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_l.y][index_k.x]+=DDF*dtC.x*dtD.y+D2KL.ay;
          if((index_l.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_l.y][index_k.y]+=DDF*dtC.y*dtD.y+D2KL.by;
          if((index_l.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_l.y][index_k.z]+=DDF*dtC.z*dtD.y+D2KL.cy;
          if((index_l.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_l.z][index_k.x]+=DDF*dtC.x*dtD.z+D2KL.az;
          if((index_l.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_l.z][index_k.y]+=DDF*dtC.y*dtD.z+D2KL.bz;
          if((index_l.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_l.z][index_k.z]+=DDF*dtC.z*dtD.z+D2KL.cz;
        }

        HessianTorsionStrainPosition(HessianMatrix,index_i,index_j,index_k,index_l,
             vec_u,vec_v,vec_w,u,v,w,fa,fc,fd,Dab,Dcb,Ddc,rbc,D2I,D2K,D2L,
             D2IJ,D2IK,D2IL,D2JK,D2JL,D2KL,dtA,dtB,dtC,dtD,DDF,S,CosPhi);

        HessianTorsionStrainStrain(HessianMatrix,vec_u,vec_v,vec_w,u,v,w,DF,DDF,S,CosPhi);
      }
    }
  }
}

void ComputeFrameworkImproperTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL *parms;
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  int index1,index2,index3,index4;
  REAL_MATRIX3x3 D2I,D2L,D2IL,D2IJ,D2IK,D2JL;
  REAL_MATRIX3x3 D2KL,D2J,D2K,D2JK,S;
  VECTOR DIL,DDJ,DDK,DEJ,DEK;
  REAL DF,DDF;
  REAL u,w,v;
  VECTOR vec_u,vec_v,vec_w;
  VECTOR fa,fb,fc,fd;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
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

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;
      index_k=Framework[CurrentSystem].Atoms[f1][C].HessianIndex;
      index_l=Framework[CurrentSystem].Atoms[f1][D].HessianIndex;

      index1=Framework[CurrentSystem].Atoms[f1][A].HessianAtomIndex;
      index2=Framework[CurrentSystem].Atoms[f1][B].HessianAtomIndex;
      index3=Framework[CurrentSystem].Atoms[f1][C].HessianAtomIndex;
      index4=Framework[CurrentSystem].Atoms[f1][D].HessianAtomIndex;

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
        case FIXED_IMPROPER_DIHEDRAL:
          U=DF=DDF=0.0;
          break;
        default:
          fprintf(stderr, "Undefined Improper Torsion potential in 'framework_hessian.c'\n");
          U=DF=DDF=0.0;
          exit(0);
          break;
      }

      *Energy+=U;

      //if((index_i<0)&&(index_j<0)&&(index_k<0)&&(index_l<0)) continue;

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
      fa.x=DF*dtA.x;
      fa.y=DF*dtA.y;
      fa.z=DF*dtA.z;

      fb.x=DF*dtB.x;
      fb.y=DF*dtB.y;
      fb.z=DF*dtB.z;

      fc.x=DF*dtC.x;
      fc.y=DF*dtC.y;
      fc.z=DF*dtC.z;

      fd.x=DF*dtD.x;
      fd.y=DF*dtD.y;
      fd.z=DF*dtD.z;

      // add contribution to the strain derivative tensor
      // Note: rbc is here because the vector was normalized before
      S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
      S.bx=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
      S.cx=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

      S.ay=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
      S.by=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
      S.cy=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

      S.az=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
      S.bz=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
      S.cz=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;

      StrainDerivative->ax+=S.ax;
      StrainDerivative->bx+=S.bx;
      StrainDerivative->cx+=S.cx;

      StrainDerivative->ay+=S.ay;
      StrainDerivative->by+=S.by;
      StrainDerivative->cy+=S.cy;

      StrainDerivative->az+=S.az;
      StrainDerivative->bz+=S.bz;
      StrainDerivative->cz+=S.cz;

      if(ComputeGradient)
      {
        if(index_i.x>=0) Gradient[index_i.x]+=fa.x;
        if(index_i.y>=0) Gradient[index_i.y]+=fa.y;
        if(index_i.z>=0) Gradient[index_i.z]+=fa.z;

        if(index_j.x>=0) Gradient[index_j.x]+=fb.x;
        if(index_j.y>=0) Gradient[index_j.y]+=fb.y;
        if(index_j.z>=0) Gradient[index_j.z]+=fb.z;

        if(index_k.x>=0) Gradient[index_k.x]+=fc.x;
        if(index_k.y>=0) Gradient[index_k.y]+=fc.y;
        if(index_k.z>=0) Gradient[index_k.z]+=fc.z;

        if(index_l.x>=0) Gradient[index_l.x]+=fd.x;
        if(index_l.y>=0) Gradient[index_l.y]+=fd.y;
        if(index_l.z>=0) Gradient[index_l.z]+=fd.z;

        GradientStrainTorsion(Gradient,S);
      }

      if(ComputeHessian)
      {
        // Calculate the derivatives of DOTIJ and DOTLK.
        DIL.x=DF*Dcb.x/rbc;
        DIL.y=DF*Dcb.y/rbc;
        DIL.z=DF*Dcb.z/rbc;
        DDJ.x=DF*((2.0*d-1.0)*Dcb.x-Dab.x/rbc)/rbc;
        DDJ.y=DF*((2.0*d-1.0)*Dcb.y-Dab.y/rbc)/rbc;
        DDJ.z=DF*((2.0*d-1.0)*Dcb.z-Dab.z/rbc)/rbc;
        DDK.x=-DF*(2.0*d*Dcb.x-Dab.x/rbc)/rbc;
        DDK.y=-DF*(2.0*d*Dcb.y-Dab.y/rbc)/rbc;
        DDK.z=-DF*(2.0*d*Dcb.z-Dab.z/rbc)/rbc;
        DEJ.x=DF*(2.0*e*Dcb.x-Ddc.x/rbc)/rbc;
        DEJ.y=DF*(2.0*e*Dcb.y-Ddc.y/rbc)/rbc;
        DEJ.z=DF*(2.0*e*Dcb.z-Ddc.z/rbc)/rbc;
        DEK.x=-DF*((2.0*e+1.0)*Dcb.x-Ddc.x/rbc)/rbc;
        DEK.y=-DF*((2.0*e+1.0)*Dcb.y-Ddc.y/rbc)/rbc;
        DEK.z=-DF*((2.0*e+1.0)*Dcb.z-Ddc.z/rbc)/rbc;

        // Calculate some diagonal Hessian terms for I.
        D2I.ax=DF*(CosPhi*(dr.x*dr.x+Dcb.x*Dcb.x-1.0)/r-2.0*dr.x*dtA.x)/r;
        D2I.by=DF*(CosPhi*(dr.y*dr.y+Dcb.y*Dcb.y-1.0)/r-2.0*dr.y*dtA.y)/r;
        D2I.cz=DF*(CosPhi*(dr.z*dr.z+Dcb.z*Dcb.z-1.0)/r-2.0*dr.z*dtA.z)/r;
        D2I.ay=DF*(CosPhi*(dr.x*dr.y+Dcb.x*Dcb.y)/r-dr.x*dtA.y-dr.y*dtA.x)/r;
        D2I.az=DF*(CosPhi*(dr.x*dr.z+Dcb.x*Dcb.z)/r-dr.x*dtA.z-dr.z*dtA.x)/r;
        D2I.bz=DF*(CosPhi*(dr.y*dr.z+Dcb.y*Dcb.z)/r-dr.y*dtA.z-dr.z*dtA.y)/r;

        // Calculate some diagonal Hessian terms for L.
        D2L.ax=DF*(CosPhi*(ds.x*ds.x+Dcb.x*Dcb.x-1.0)/s-2.0*ds.x*dtD.x)/s;
        D2L.by=DF*(CosPhi*(ds.y*ds.y+Dcb.y*Dcb.y-1.0)/s-2.0*ds.y*dtD.y)/s;
        D2L.cz=DF*(CosPhi*(ds.z*ds.z+Dcb.z*Dcb.z-1.0)/s-2.0*ds.z*dtD.z)/s;
        D2L.ay=DF*(CosPhi*(ds.x*ds.y+Dcb.x*Dcb.y)/s-ds.x*dtD.y-ds.y*dtD.x)/s;
        D2L.az=DF*(CosPhi*(ds.x*ds.z+Dcb.x*Dcb.z)/s-ds.x*dtD.z-ds.z*dtD.x)/s;
        D2L.bz=DF*(CosPhi*(ds.y*ds.z+Dcb.y*Dcb.z)/s-ds.y*dtD.z-ds.z*dtD.y)/s;

        // Calculate the IL off-diagonal terms.
        D2IL.ax=DF*(CosPhi*dr.x*ds.x-dr.x*dr.x-ds.x*ds.x-Dcb.x*Dcb.x+1.0)/(r*s);
        D2IL.ay=DF*(CosPhi*dr.x*ds.y-dr.x*dr.y-ds.x*ds.y-Dcb.x*Dcb.y)/(r*s);
        D2IL.az=DF*(CosPhi*dr.x*ds.z-dr.x*dr.z-ds.x*ds.z-Dcb.x*Dcb.z)/(r*s);
        D2IL.bx=DF*(CosPhi*dr.y*ds.x-dr.y*dr.x-ds.y*ds.x-Dcb.y*Dcb.x)/(r*s);
        D2IL.by=DF*(CosPhi*dr.y*ds.y-dr.y*dr.y-ds.y*ds.y-Dcb.y*Dcb.y+1.0)/(r*s);
        D2IL.bz=DF*(CosPhi*dr.y*ds.z-dr.y*dr.z-ds.y*ds.z-Dcb.y*Dcb.z)/(r*s);
        D2IL.cx=DF*(CosPhi*dr.z*ds.x-dr.z*dr.x-ds.z*ds.x-Dcb.z*Dcb.x)/(r*s);
        D2IL.cy=DF*(CosPhi*dr.z*ds.y-dr.z*dr.y-ds.z*ds.y-Dcb.z*Dcb.y)/(r*s);
        D2IL.cz=DF*(CosPhi*dr.z*ds.z-dr.z*dr.z-ds.z*ds.z-Dcb.z*Dcb.z+1.0)/(r*s);

        // Calculate the IJ off-diagonal terms.
        D2IJ.ax=D2I.ax*(d-1.0)+D2IL.ax*e+DIL.x*dtA.x;
        D2IJ.ay=D2I.ay*(d-1.0)+D2IL.ay*e+DIL.x*dtA.y;
        D2IJ.az=D2I.az*(d-1.0)+D2IL.az*e+DIL.x*dtA.z;
        D2IJ.bx=D2I.ay*(d-1.0)+D2IL.bx*e+DIL.y*dtA.x;
        D2IJ.by=D2I.by*(d-1.0)+D2IL.by*e+DIL.y*dtA.y;
        D2IJ.bz=D2I.bz*(d-1.0)+D2IL.bz*e+DIL.y*dtA.z;
        D2IJ.cx=D2I.az*(d-1.0)+D2IL.cx*e+DIL.z*dtA.x;
        D2IJ.cy=D2I.bz*(d-1.0)+D2IL.cy*e+DIL.z*dtA.y;
        D2IJ.cz=D2I.cz*(d-1.0)+D2IL.cz*e+DIL.z*dtA.z;

        // Calculate the IK off-diagonal terms.
        D2IK.ax=-D2I.ax*d-D2IL.ax*(e+1.0)-DIL.x*dtA.x;
        D2IK.ay=-D2I.ay*d-D2IL.ay*(e+1.0)-DIL.x*dtA.y;
        D2IK.az=-D2I.az*d-D2IL.az*(e+1.0)-DIL.x*dtA.z;
        D2IK.bx=-D2I.ay*d-D2IL.bx*(e+1.0)-DIL.y*dtA.x;
        D2IK.by=-D2I.by*d-D2IL.by*(e+1.0)-DIL.y*dtA.y;
        D2IK.bz=-D2I.bz*d-D2IL.bz*(e+1.0)-DIL.y*dtA.z;
        D2IK.cx=-D2I.az*d-D2IL.cx*(e+1.0)-DIL.z*dtA.x;
        D2IK.cy=-D2I.bz*d-D2IL.cy*(e+1.0)-DIL.z*dtA.y;
        D2IK.cz=-D2I.cz*d-D2IL.cz*(e+1.0)-DIL.z*dtA.z;

        // Calculate the JL off-diagonal terms.
        D2JL.ax=D2IL.ax*(d-1.0)+D2L.ax*e+DIL.x*dtD.x;
        D2JL.ay=D2IL.ay*(d-1.0)+D2L.ay*e+DIL.y*dtD.x;
        D2JL.az=D2IL.az*(d-1.0)+D2L.az*e+DIL.z*dtD.x;
        D2JL.bx=D2IL.bx*(d-1.0)+D2L.ay*e+DIL.x*dtD.y;
        D2JL.by=D2IL.by*(d-1.0)+D2L.by*e+DIL.y*dtD.y;
        D2JL.bz=D2IL.bz*(d-1.0)+D2L.bz*e+DIL.z*dtD.y;
        D2JL.cx=D2IL.cx*(d-1.0)+D2L.az*e+DIL.x*dtD.z;
        D2JL.cy=D2IL.cy*(d-1.0)+D2L.bz*e+DIL.y*dtD.z;
        D2JL.cz=D2IL.cz*(d-1.0)+D2L.cz*e+DIL.z*dtD.z;

        // Calculate the KL off-diagonal terms.
        D2KL.ax=-D2IL.ax*d-D2L.ax*(e+1.0)-DIL.x*dtD.x;
        D2KL.ay=-D2IL.ay*d-D2L.ay*(e+1.0)-DIL.y*dtD.x;
        D2KL.az=-D2IL.az*d-D2L.az*(e+1.0)-DIL.z*dtD.x;
        D2KL.bx=-D2IL.bx*d-D2L.ay*(e+1.0)-DIL.x*dtD.y;
        D2KL.by=-D2IL.by*d-D2L.by*(e+1.0)-DIL.y*dtD.y;
        D2KL.bz=-D2IL.bz*d-D2L.bz*(e+1.0)-DIL.z*dtD.y;
        D2KL.cx=-D2IL.cx*d-D2L.az*(e+1.0)-DIL.x*dtD.z;
        D2KL.cy=-D2IL.cy*d-D2L.bz*(e+1.0)-DIL.y*dtD.z;
        D2KL.cz=-D2IL.cz*d-D2L.cz*(e+1.0)-DIL.z*dtD.z;

        // Calculate the JJ diagonal terms.
        D2J.ax=D2IJ.ax*(d-1.0)+D2JL.ax*e+DDJ.x*dtA.x+DEJ.x*dtD.x;
        D2J.by=D2IJ.by*(d-1.0)+D2JL.by*e+DDJ.y*dtA.y+DEJ.y*dtD.y;
        D2J.cz=D2IJ.cz*(d-1.0)+D2JL.cz*e+DDJ.z*dtA.z+DEJ.z*dtD.z;
        D2J.ay=D2IJ.ay*(d-1.0)+D2JL.bx*e+DDJ.y*dtA.x+DEJ.y*dtD.x;
        D2J.az=D2IJ.az*(d-1.0)+D2JL.cx*e+DDJ.z*dtA.x+DEJ.z*dtD.x;
        D2J.bz=D2IJ.bz*(d-1.0)+D2JL.cy*e+DDJ.z*dtA.y+DEJ.z*dtD.y;

        // Calculate the KK diagonal terms.
        D2K.ax=-D2KL.ax*(e+1.0)-D2IK.ax*d-DDK.x*dtA.x-DEK.x*dtD.x;
        D2K.by=-D2KL.by*(e+1.0)-D2IK.by*d-DDK.y*dtA.y-DEK.y*dtD.y;
        D2K.cz=-D2KL.cz*(e+1.0)-D2IK.cz*d-DDK.z*dtA.z-DEK.z*dtD.z;
        D2K.ay=-D2KL.ay*(e+1.0)-D2IK.bx*d-DDK.x*dtA.y-DEK.x*dtD.y;
        D2K.az=-D2KL.az*(e+1.0)-D2IK.cx*d-DDK.x*dtA.z-DEK.x*dtD.z;
        D2K.bz=-D2KL.bz*(e+1.0)-D2IK.cy*d-DDK.y*dtA.z-DEK.y*dtD.z;

        // Calculate the JK off-diagonal terms.
        D2JK.ax=-D2IJ.ax*d-D2JL.ax*(e+1.0)-DDJ.x*dtA.x-DEJ.x*dtD.x;
        D2JK.ay=-D2IJ.bx*d-D2JL.ay*(e+1.0)-DDJ.x*dtA.y-DEJ.x*dtD.y;
        D2JK.az=-D2IJ.cx*d-D2JL.az*(e+1.0)-DDJ.x*dtA.z-DEJ.x*dtD.z;
        D2JK.bx=-D2IJ.ay*d-D2JL.bx*(e+1.0)-DDJ.y*dtA.x-DEJ.y*dtD.x;
        D2JK.by=-D2IJ.by*d-D2JL.by*(e+1.0)-DDJ.y*dtA.y-DEJ.y*dtD.y;
        D2JK.bz=-D2IJ.cy*d-D2JL.bz*(e+1.0)-DDJ.y*dtA.z-DEJ.y*dtD.z;
        D2JK.cx=-D2IJ.az*d-D2JL.cx*(e+1.0)-DDJ.z*dtA.x-DEJ.z*dtD.x;
        D2JK.cy=-D2IJ.bz*d-D2JL.cy*(e+1.0)-DDJ.z*dtA.y-DEJ.z*dtD.y;
        D2JK.cz=-D2IJ.cz*d-D2JL.cz*(e+1.0)-DDJ.z*dtA.z-DEJ.z*dtD.z;

        // Calculate the diagonal blocks of the hessian.
        if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x]+=DDF*dtA.x*dtA.x+D2I.ax;
        if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y]+=DDF*dtA.x*dtA.y+D2I.ay;
        if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y]+=DDF*dtA.y*dtA.y+D2I.by;
        if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z]+=DDF*dtA.x*dtA.z+D2I.az;
        if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z]+=DDF*dtA.y*dtA.z+D2I.bz;
        if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z]+=DDF*dtA.z*dtA.z+D2I.cz;

        if((index_j.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_j.x][index_j.x]+=DDF*dtB.x*dtB.x+D2J.ax;
        if((index_j.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.x][index_j.y]+=DDF*dtB.x*dtB.y+D2J.ay;
        if((index_j.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.y][index_j.y]+=DDF*dtB.y*dtB.y+D2J.by;
        if((index_j.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.x][index_j.z]+=DDF*dtB.x*dtB.z+D2J.az;
        if((index_j.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.y][index_j.z]+=DDF*dtB.y*dtB.z+D2J.bz;
        if((index_j.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.z][index_j.z]+=DDF*dtB.z*dtB.z+D2J.cz;

        if((index_k.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_k.x][index_k.x]+=DDF*dtC.x*dtC.x+D2K.ax;
        if((index_k.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_k.x][index_k.y]+=DDF*dtC.x*dtC.y+D2K.ay;
        if((index_k.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_k.y][index_k.y]+=DDF*dtC.y*dtC.y+D2K.by;
        if((index_k.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.x][index_k.z]+=DDF*dtC.x*dtC.z+D2K.az;
        if((index_k.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.y][index_k.z]+=DDF*dtC.y*dtC.z+D2K.bz;
        if((index_k.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.z][index_k.z]+=DDF*dtC.z*dtC.z+D2K.cz;

        if((index_l.x>=0)&&(index_l.x>=0)) HessianMatrix.element[index_l.x][index_l.x]+=DDF*dtD.x*dtD.x+D2L.ax;
        if((index_l.x>=0)&&(index_l.y>=0)) HessianMatrix.element[index_l.x][index_l.y]+=DDF*dtD.x*dtD.y+D2L.ay;
        if((index_l.y>=0)&&(index_l.y>=0)) HessianMatrix.element[index_l.y][index_l.y]+=DDF*dtD.y*dtD.y+D2L.by;
        if((index_l.x>=0)&&(index_l.z>=0)) HessianMatrix.element[index_l.x][index_l.z]+=DDF*dtD.x*dtD.z+D2L.az;
        if((index_l.y>=0)&&(index_l.z>=0)) HessianMatrix.element[index_l.y][index_l.z]+=DDF*dtD.y*dtD.z+D2L.bz;
        if((index_l.z>=0)&&(index_l.z>=0)) HessianMatrix.element[index_l.z][index_l.z]+=DDF*dtD.z*dtD.z+D2L.cz;

        if(index1<index2)
        {
          if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]+=DDF*dtA.x*dtB.x+D2IJ.ax;
          if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]+=DDF*dtA.x*dtB.y+D2IJ.ay;
          if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]+=DDF*dtA.x*dtB.z+D2IJ.az;
          if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]+=DDF*dtA.y*dtB.x+D2IJ.bx;
          if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]+=DDF*dtA.y*dtB.y+D2IJ.by;
          if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]+=DDF*dtA.y*dtB.z+D2IJ.bz;
          if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]+=DDF*dtA.z*dtB.x+D2IJ.cx;
          if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]+=DDF*dtA.z*dtB.y+D2IJ.cy;
          if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]+=DDF*dtA.z*dtB.z+D2IJ.cz;
        }
        else
        {
          if((index_j.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.x][index_i.x]+=DDF*dtA.x*dtB.x+D2IJ.ax;
          if((index_j.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.x][index_i.y]+=DDF*dtA.y*dtB.x+D2IJ.bx;
          if((index_j.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.x][index_i.z]+=DDF*dtA.z*dtB.x+D2IJ.cx;
          if((index_j.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.y][index_i.x]+=DDF*dtA.x*dtB.y+D2IJ.ay;
          if((index_j.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.y][index_i.y]+=DDF*dtA.y*dtB.y+D2IJ.by;
          if((index_j.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.y][index_i.z]+=DDF*dtA.z*dtB.y+D2IJ.cy;
          if((index_j.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.z][index_i.x]+=DDF*dtA.x*dtB.z+D2IJ.az;
          if((index_j.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.z][index_i.y]+=DDF*dtA.y*dtB.z+D2IJ.bz;
          if((index_j.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.z][index_i.z]+=DDF*dtA.z*dtB.z+D2IJ.cz;
        }

        if(index1<index3)
        {
          if((index_i.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.x][index_k.x]+=DDF*dtA.x*dtC.x+D2IK.ax;
          if((index_i.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.x][index_k.y]+=DDF*dtA.x*dtC.y+D2IK.ay;
          if((index_i.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.x][index_k.z]+=DDF*dtA.x*dtC.z+D2IK.az;
          if((index_i.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.y][index_k.x]+=DDF*dtA.y*dtC.x+D2IK.bx;
          if((index_i.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.y][index_k.y]+=DDF*dtA.y*dtC.y+D2IK.by;
          if((index_i.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.y][index_k.z]+=DDF*dtA.y*dtC.z+D2IK.bz;
          if((index_i.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.z][index_k.x]+=DDF*dtA.z*dtC.x+D2IK.cx;
          if((index_i.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.z][index_k.y]+=DDF*dtA.z*dtC.y+D2IK.cy;
          if((index_i.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.z][index_k.z]+=DDF*dtA.z*dtC.z+D2IK.cz;
        }
        else
        {
          if((index_k.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.x][index_i.x]+=DDF*dtA.x*dtC.x+D2IK.ax;
          if((index_k.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.x][index_i.y]+=DDF*dtA.y*dtC.x+D2IK.bx;
          if((index_k.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.x][index_i.z]+=DDF*dtA.z*dtC.x+D2IK.cx;
          if((index_k.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.y][index_i.x]+=DDF*dtA.x*dtC.y+D2IK.ay;
          if((index_k.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.y][index_i.y]+=DDF*dtA.y*dtC.y+D2IK.by;
          if((index_k.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.y][index_i.z]+=DDF*dtA.z*dtC.y+D2IK.cy;
          if((index_k.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.z][index_i.x]+=DDF*dtA.x*dtC.z+D2IK.az;
          if((index_k.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.z][index_i.y]+=DDF*dtA.y*dtC.z+D2IK.bz;
          if((index_k.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.z][index_i.z]+=DDF*dtA.z*dtC.z+D2IK.cz;
        }

        if(index1<index4)
        {
          if((index_i.x>=0)&&(index_l.x>=0)) HessianMatrix.element[index_i.x][index_l.x]+=DDF*dtA.x*dtD.x+D2IL.ax;
          if((index_i.x>=0)&&(index_l.y>=0)) HessianMatrix.element[index_i.x][index_l.y]+=DDF*dtA.x*dtD.y+D2IL.ay;
          if((index_i.x>=0)&&(index_l.z>=0)) HessianMatrix.element[index_i.x][index_l.z]+=DDF*dtA.x*dtD.z+D2IL.az;
          if((index_i.y>=0)&&(index_l.x>=0)) HessianMatrix.element[index_i.y][index_l.x]+=DDF*dtA.y*dtD.x+D2IL.bx;
          if((index_i.y>=0)&&(index_l.y>=0)) HessianMatrix.element[index_i.y][index_l.y]+=DDF*dtA.y*dtD.y+D2IL.by;
          if((index_i.y>=0)&&(index_l.z>=0)) HessianMatrix.element[index_i.y][index_l.z]+=DDF*dtA.y*dtD.z+D2IL.bz;
          if((index_i.z>=0)&&(index_l.x>=0)) HessianMatrix.element[index_i.z][index_l.x]+=DDF*dtA.z*dtD.x+D2IL.cx;
          if((index_i.z>=0)&&(index_l.y>=0)) HessianMatrix.element[index_i.z][index_l.y]+=DDF*dtA.z*dtD.y+D2IL.cy;
          if((index_i.z>=0)&&(index_l.z>=0)) HessianMatrix.element[index_i.z][index_l.z]+=DDF*dtA.z*dtD.z+D2IL.cz;
        }
        else
        {
          if((index_l.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_l.x][index_i.x]+=DDF*dtA.x*dtD.x+D2IL.ax;
          if((index_l.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_l.x][index_i.y]+=DDF*dtA.y*dtD.x+D2IL.bx;
          if((index_l.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_l.x][index_i.z]+=DDF*dtA.z*dtD.x+D2IL.cx;
          if((index_l.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_l.y][index_i.x]+=DDF*dtA.x*dtD.y+D2IL.ay;
          if((index_l.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_l.y][index_i.y]+=DDF*dtA.y*dtD.y+D2IL.by;
          if((index_l.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_l.y][index_i.z]+=DDF*dtA.z*dtD.y+D2IL.cy;
          if((index_l.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_l.z][index_i.x]+=DDF*dtA.x*dtD.z+D2IL.az;
          if((index_l.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_l.z][index_i.y]+=DDF*dtA.y*dtD.z+D2IL.bz;
          if((index_l.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_l.z][index_i.z]+=DDF*dtA.z*dtD.z+D2IL.cz;
        }

        if(index2<index3)
        {
          if((index_j.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.x][index_k.x]+=DDF*dtB.x*dtC.x+D2JK.ax;
          if((index_j.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.x][index_k.y]+=DDF*dtB.x*dtC.y+D2JK.ay;
          if((index_j.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.x][index_k.z]+=DDF*dtB.x*dtC.z+D2JK.az;
          if((index_j.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.y][index_k.x]+=DDF*dtB.y*dtC.x+D2JK.bx;
          if((index_j.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.y][index_k.y]+=DDF*dtB.y*dtC.y+D2JK.by;
          if((index_j.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.y][index_k.z]+=DDF*dtB.y*dtC.z+D2JK.bz;
          if((index_j.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.z][index_k.x]+=DDF*dtB.z*dtC.x+D2JK.cx;
          if((index_j.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.z][index_k.y]+=DDF*dtB.z*dtC.y+D2JK.cy;
          if((index_j.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.z][index_k.z]+=DDF*dtB.z*dtC.z+D2JK.cz;
        }
        else
        {
          if((index_k.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.x][index_j.x]+=DDF*dtB.x*dtC.x+D2JK.ax;
          if((index_k.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.x][index_j.y]+=DDF*dtB.y*dtC.x+D2JK.bx;
          if((index_k.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.x][index_j.z]+=DDF*dtB.z*dtC.x+D2JK.cx;
          if((index_k.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.y][index_j.x]+=DDF*dtB.x*dtC.y+D2JK.ay;
          if((index_k.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.y][index_j.y]+=DDF*dtB.y*dtC.y+D2JK.by;
          if((index_k.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.y][index_j.z]+=DDF*dtB.z*dtC.y+D2JK.cy;
          if((index_k.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.z][index_j.x]+=DDF*dtB.x*dtC.z+D2JK.az;
          if((index_k.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.z][index_j.y]+=DDF*dtB.y*dtC.z+D2JK.bz;
          if((index_k.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.z][index_j.z]+=DDF*dtB.z*dtC.z+D2JK.cz;
        }

        if(index2<index4)
        {
          if((index_j.x>=0)&&(index_l.x>=0)) HessianMatrix.element[index_j.x][index_l.x]+=DDF*dtB.x*dtD.x+D2JL.ax;
          if((index_j.x>=0)&&(index_l.y>=0)) HessianMatrix.element[index_j.x][index_l.y]+=DDF*dtB.x*dtD.y+D2JL.ay;
          if((index_j.x>=0)&&(index_l.z>=0)) HessianMatrix.element[index_j.x][index_l.z]+=DDF*dtB.x*dtD.z+D2JL.az;
          if((index_j.y>=0)&&(index_l.x>=0)) HessianMatrix.element[index_j.y][index_l.x]+=DDF*dtB.y*dtD.x+D2JL.bx;
          if((index_j.y>=0)&&(index_l.y>=0)) HessianMatrix.element[index_j.y][index_l.y]+=DDF*dtB.y*dtD.y+D2JL.by;
          if((index_j.y>=0)&&(index_l.z>=0)) HessianMatrix.element[index_j.y][index_l.z]+=DDF*dtB.y*dtD.z+D2JL.bz;
          if((index_j.z>=0)&&(index_l.x>=0)) HessianMatrix.element[index_j.z][index_l.x]+=DDF*dtB.z*dtD.x+D2JL.cx;
          if((index_j.z>=0)&&(index_l.y>=0)) HessianMatrix.element[index_j.z][index_l.y]+=DDF*dtB.z*dtD.y+D2JL.cy;
          if((index_j.z>=0)&&(index_l.z>=0)) HessianMatrix.element[index_j.z][index_l.z]+=DDF*dtB.z*dtD.z+D2JL.cz;
        }
        else
        {
          if((index_l.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_l.x][index_j.x]+=DDF*dtB.x*dtD.x+D2JL.ax;
          if((index_l.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_l.x][index_j.y]+=DDF*dtB.y*dtD.x+D2JL.bx;
          if((index_l.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_l.x][index_j.z]+=DDF*dtB.z*dtD.x+D2JL.cx;
          if((index_l.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_l.y][index_j.x]+=DDF*dtB.x*dtD.y+D2JL.ay;
          if((index_l.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_l.y][index_j.y]+=DDF*dtB.y*dtD.y+D2JL.by;
          if((index_l.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_l.y][index_j.z]+=DDF*dtB.z*dtD.y+D2JL.cy;
          if((index_l.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_l.z][index_j.x]+=DDF*dtB.x*dtD.z+D2JL.az;
          if((index_l.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_l.z][index_j.y]+=DDF*dtB.y*dtD.z+D2JL.bz;
          if((index_l.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_l.z][index_j.z]+=DDF*dtB.z*dtD.z+D2JL.cz;
        }

        if(index3<index4)
        {
          if((index_k.x>=0)&&(index_l.x>=0)) HessianMatrix.element[index_k.x][index_l.x]+=DDF*dtC.x*dtD.x+D2KL.ax;
          if((index_k.x>=0)&&(index_l.y>=0)) HessianMatrix.element[index_k.x][index_l.y]+=DDF*dtC.x*dtD.y+D2KL.ay;
          if((index_k.x>=0)&&(index_l.z>=0)) HessianMatrix.element[index_k.x][index_l.z]+=DDF*dtC.x*dtD.z+D2KL.az;
          if((index_k.y>=0)&&(index_l.x>=0)) HessianMatrix.element[index_k.y][index_l.x]+=DDF*dtC.y*dtD.x+D2KL.bx;
          if((index_k.y>=0)&&(index_l.y>=0)) HessianMatrix.element[index_k.y][index_l.y]+=DDF*dtC.y*dtD.y+D2KL.by;
          if((index_k.y>=0)&&(index_l.z>=0)) HessianMatrix.element[index_k.y][index_l.z]+=DDF*dtC.y*dtD.z+D2KL.bz;
          if((index_k.z>=0)&&(index_l.x>=0)) HessianMatrix.element[index_k.z][index_l.x]+=DDF*dtC.z*dtD.x+D2KL.cx;
          if((index_k.z>=0)&&(index_l.y>=0)) HessianMatrix.element[index_k.z][index_l.y]+=DDF*dtC.z*dtD.y+D2KL.cy;
          if((index_k.z>=0)&&(index_l.z>=0)) HessianMatrix.element[index_k.z][index_l.z]+=DDF*dtC.z*dtD.z+D2KL.cz;
        }
        else
        {
          if((index_l.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_l.x][index_k.x]+=DDF*dtC.x*dtD.x+D2KL.ax;
          if((index_l.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_l.x][index_k.y]+=DDF*dtC.y*dtD.x+D2KL.bx;
          if((index_l.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_l.x][index_k.z]+=DDF*dtC.z*dtD.x+D2KL.cx;
          if((index_l.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_l.y][index_k.x]+=DDF*dtC.x*dtD.y+D2KL.ay;
          if((index_l.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_l.y][index_k.y]+=DDF*dtC.y*dtD.y+D2KL.by;
          if((index_l.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_l.y][index_k.z]+=DDF*dtC.z*dtD.y+D2KL.cy;
          if((index_l.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_l.z][index_k.x]+=DDF*dtC.x*dtD.z+D2KL.az;
          if((index_l.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_l.z][index_k.y]+=DDF*dtC.y*dtD.z+D2KL.bz;
          if((index_l.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_l.z][index_k.z]+=DDF*dtC.z*dtD.z+D2KL.cz;
        }

        HessianTorsionStrainPosition(HessianMatrix,index_i,index_j,index_k,index_l,
               vec_u,vec_v,vec_w,u,v,w,fa,fc,fd,Dab,Dcb,Ddc,rbc,D2I,D2K,D2L,
               D2IJ,D2IK,D2IL,D2JK,D2JL,D2KL,dtA,dtB,dtC,dtD,DDF,S,CosPhi);

        HessianTorsionStrainStrain(HessianMatrix,vec_u,vec_v,vec_w,u,v,w,DF,DDF,S,CosPhi);
      }
    }
  }
}


// Numerical routines
// TODO: convert them to analytical expressions

void CalculateFrameworkInversionBendForces(int f1,int i,VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,REAL *Energy,VECTOR *fa,VECTOR *fb,VECTOR *fc,VECTOR *fd,REAL_MATRIX3x3 *strain_derivative)
{
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2,rrbc,rbc2,rbd2,rad2,rac2,dot;
  REAL CosChi,Chi,energy;
  REAL temp,temp2;
  VECTOR Rab,Rbc,Rbd,Rad,Rac;
  REAL term,dedcos;
  VECTOR dccdia,dccdic,dccdid;
  VECTOR deedia,deedic,deedid;

  parms=Framework[CurrentSystem].InversionBendArguments[f1][i];

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
        *Energy=energy;

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

        fa->x=-dedcos*(dccdia.x+deedia.x);
        fa->y=-dedcos*(dccdia.y+deedia.y);
        fa->z=-dedcos*(dccdia.z+deedia.z);
        fc->x=-dedcos*(dccdic.x+deedic.x);
        fc->y=-dedcos*(dccdic.y+deedic.y);
        fc->z=-dedcos*(dccdic.z+deedic.z);
        fd->x=-dedcos*(dccdid.x+deedid.x);
        fd->y=-dedcos*(dccdid.y+deedid.y);
        fd->z=-dedcos*(dccdid.z+deedid.z);

        fb->x=-(fa->x+fc->x+fd->x);
        fb->y=-(fa->y+fc->y+fd->y);
        fb->z=-(fa->z+fc->z+fd->z);

        strain_derivative->ax=Rab.x*fa->x+Rbc.x*fc->x+Rbd.x*fd->x;
        strain_derivative->ay=Rab.x*fa->y+Rbc.x*fc->y+Rbd.x*fd->y;
        strain_derivative->az=Rab.x*fa->z+Rbc.x*fc->z+Rbd.x*fd->z;

        strain_derivative->bx=Rab.y*fa->x+Rbc.y*fc->x+Rbd.y*fd->x;
        strain_derivative->by=Rab.y*fa->y+Rbc.y*fc->y+Rbd.y*fd->y;
        strain_derivative->bz=Rab.y*fa->z+Rbc.y*fc->z+Rbd.y*fd->z;

        strain_derivative->cx=Rab.z*fa->x+Rbc.z*fc->x+Rbd.z*fd->x;
        strain_derivative->cy=Rab.z*fa->y+Rbc.z*fc->y+Rbd.z*fd->y;
        strain_derivative->cz=Rab.z*fa->z+Rbc.z*fc->z+Rbd.z*fd->z;

}


void ComputeFrameworkInversionBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivativeTensor,int ComputeGradient,int ComputeHessian)
{
  int i,f1;
  int A,B,C,D;
  const REAL Delta=1e-7;
  VECTOR posA0,posB0,posC0,posD0;
  VECTOR posA,posB,posC,posD;
  VECTOR fa0,fb0,fc0,fd0,fa[4],fb[4],fc[4],fd[4];
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  REAL U;
  REAL_MATRIX3x3 strain_derivative;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfInversionBends[f1];i++)
    {
      A=Framework[CurrentSystem].InversionBends[f1][i].A;
      B=Framework[CurrentSystem].InversionBends[f1][i].B;
      C=Framework[CurrentSystem].InversionBends[f1][i].C;
      D=Framework[CurrentSystem].InversionBends[f1][i].D;

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;
      index_k=Framework[CurrentSystem].Atoms[f1][C].HessianIndex;
      index_l=Framework[CurrentSystem].Atoms[f1][D].HessianIndex;

      posA0=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB0=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC0=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD0=Framework[CurrentSystem].Atoms[f1][D].Position;

      CalculateFrameworkInversionBendForces(f1,i,posA0,posB0,posC0,posD0,&U,&fa0,&fb0,&fc0,&fd0,&strain_derivative);

      *Energy+=U;

      StrainDerivativeTensor->ax+=strain_derivative.ax;
      StrainDerivativeTensor->ay+=strain_derivative.ay;
      StrainDerivativeTensor->az+=strain_derivative.az;

      StrainDerivativeTensor->bx+=strain_derivative.bx;
      StrainDerivativeTensor->by+=strain_derivative.by;
      StrainDerivativeTensor->bz+=strain_derivative.bz;

      StrainDerivativeTensor->cx+=strain_derivative.cx;
      StrainDerivativeTensor->cy+=strain_derivative.cy;
      StrainDerivativeTensor->cz+=strain_derivative.cz;

      // add contribution to the first derivatives
      if(ComputeGradient)
      {
        if(index_i.x>=0) Gradient[index_i.x]-=fa0.x;
        if(index_i.y>=0) Gradient[index_i.y]-=fa0.y;
        if(index_i.z>=0) Gradient[index_i.z]-=fa0.z;

        if(index_j.x>=0) Gradient[index_j.x]-=fb0.x;
        if(index_j.y>=0) Gradient[index_j.y]-=fb0.y;
        if(index_j.z>=0) Gradient[index_j.z]-=fb0.z;

        if(index_k.x>=0) Gradient[index_k.x]-=fc0.x;
        if(index_k.y>=0) Gradient[index_k.y]-=fc0.y;
        if(index_k.z>=0) Gradient[index_k.z]-=fc0.z;

        if(index_l.x>=0) Gradient[index_l.x]-=fd0.x;
        if(index_l.y>=0) Gradient[index_l.y]-=fd0.y;
        if(index_l.z>=0) Gradient[index_l.z]-=fd0.z;
      }

      if(ComputeHessian)
      {
        // Atom A

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_i.x>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_i.x][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_i.x][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_i.x][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_i.x][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_i.x][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_i.x][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_i.x][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_i.x][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_i.x][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_i.x][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_i.x][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_i.x][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_i.y>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_i.y][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_i.y][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_i.y][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_i.y][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_i.y][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_i.y][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_i.y][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_i.y][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_i.y][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_i.y][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_i.y][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_i.y][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_i.z>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_i.z][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_i.z][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_i.z][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_i.z][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_i.z][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_i.z][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_i.z][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_i.z][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_i.z][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_i.z][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_i.z][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_i.z][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        // Atom B

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_j.x>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_j.x][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_j.x][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_j.x][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_j.x][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_j.x][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_j.x][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_j.x][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_j.x][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_j.x][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_j.x][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_j.x][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_j.x][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_j.y>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_j.y][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_j.y][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_j.y][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_j.y][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_j.y][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_j.y][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_j.y][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_j.y][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_j.y][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_j.y][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_j.y][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_j.z>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_j.z][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_j.z][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_j.z][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_j.z][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_j.z][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_j.z][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_j.z][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_j.z][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_j.z][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_j.z][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_j.z][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        // Atom C

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_k.x>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_k.x][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_k.x][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_k.x][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_k.x][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_k.x][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_k.x][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_k.x][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_k.x][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_k.x][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_k.x][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_k.x][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_k.y>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_k.y][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_k.y][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_k.y][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_k.y][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_k.y][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_k.y][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_k.y][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_k.y][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_k.y][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_k.y][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_k.y][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_k.y][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_k.z>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_k.z][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_k.z][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_k.z][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_k.z][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_k.z][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_k.z][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_k.z][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_k.z][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_k.z][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_k.z][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_k.z][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_k.z][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        // Atom D

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_l.x>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_l.x][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_l.x][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_l.x][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_l.x][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_l.x][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_l.x][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_l.x][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_l.x][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_l.x][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_l.x][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_l.x][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_l.x][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_l.y>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_l.y][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_l.y][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_l.y][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_l.y][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_l.y][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_l.y][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_l.y][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_l.y][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_l.y][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_l.y][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_l.y][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_l.y][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z+=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z+=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z-=0.5*Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z-=Delta;
        CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_l.z>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_l.z][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_l.z][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_l.z][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_l.z][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_l.z][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_l.z][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_l.z][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_l.z][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_l.z][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_l.z][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_l.z][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_l.z][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }
      }
    }
  }
}

void CalculateFrameworkBendTorsionForces(int f1,int i,VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,REAL *Energy,VECTOR *fan,VECTOR *fbn,VECTOR *fcn,VECTOR *fdn,REAL_MATRIX3x3 *strain_derivativen)
{
  REAL d,e, energy,rab,rbc,rcd;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2,DCos;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL DTheta1,DTheta2,sign,Phi,SinPhi;
  REAL *parms;
  VECTOR fa,fb,fc,fd;

  fan->x=fbn->x=fcn->x=fdn->x=0.0;
  fan->y=fbn->y=fcn->y=fdn->y=0.0;
  fan->z=fbn->z=fcn->z=fdn->z=0.0;

  strain_derivativen->ax=strain_derivativen->bx=strain_derivativen->cx=0.0;
  strain_derivativen->ay=strain_derivativen->by=strain_derivativen->cy=0.0;
  strain_derivativen->az=strain_derivativen->bz=strain_derivativen->cz=0.0;

  parms=(REAL*)&Framework[CurrentSystem].BendTorsionArguments[f1][i];

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
  *Energy=energy;

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

  // forces torsion are oppositely directed to the gradient
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

  fan->x+=fa.x; fan->y+=fa.y; fan->z+=fa.z;
  fbn->x+=fb.x; fbn->y+=fb.y; fbn->z+=fb.z;
  fcn->x+=fc.x; fcn->y+=fc.y; fcn->z+=fc.z;
  fdn->x+=fd.x; fdn->y+=fd.y; fdn->z+=fd.z;

  // add contribution to the stress tensor
  // Note: rbc is here because the vector was normalized before
  strain_derivativen->ax-=Dab.x*fa.x+Dbc.x*rbc*(fc.x+fd.x)+Dcd.x*fd.x;
  strain_derivativen->bx-=Dab.y*fa.x+Dbc.y*rbc*(fc.x+fd.x)+Dcd.y*fd.x;
  strain_derivativen->cx-=Dab.z*fa.x+Dbc.z*rbc*(fc.x+fd.x)+Dcd.z*fd.x;

  strain_derivativen->ay-=Dab.x*fa.y+Dbc.x*rbc*(fc.y+fd.y)+Dcd.x*fd.y;
  strain_derivativen->by-=Dab.y*fa.y+Dbc.y*rbc*(fc.y+fd.y)+Dcd.y*fd.y;
  strain_derivativen->cy-=Dab.z*fa.y+Dbc.z*rbc*(fc.y+fd.y)+Dcd.z*fd.y;

  strain_derivativen->az-=Dab.x*fa.z+Dbc.x*rbc*(fc.z+fd.z)+Dcd.x*fd.z;
  strain_derivativen->bz-=Dab.y*fa.z+Dbc.y*rbc*(fc.z+fd.z)+Dcd.y*fd.z;
  strain_derivativen->cz-=Dab.z*fa.z+Dbc.z*rbc*(fc.z+fd.z)+Dcd.z*fd.z;

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

  fan->x+=fa.x;
  fan->y+=fa.y;
  fan->z+=fa.z;

  fbn->x+=fb.x-fa.x-fc.x;
  fbn->y+=fb.y-fa.y-fc.y;
  fbn->z+=fb.z-fa.z-fc.z;

  fcn->x+=fc.x-fb.x-fd.x;
  fcn->y+=fc.y-fb.y-fd.y;
  fcn->z+=fc.z-fb.z-fd.z;

  fdn->x+=fd.x;
  fdn->y+=fd.y;
  fdn->z+=fd.z;

  strain_derivativen->ax-=rab*Dab.x*fa.x+rbc*Dbc.x*fc.x+rbc*Dbc.x*fb.x+rcd*Dcd.x*fd.x;
  strain_derivativen->bx-=rab*Dab.y*fa.x+rbc*Dbc.y*fc.x+rbc*Dbc.y*fb.x+rcd*Dcd.y*fd.x;
  strain_derivativen->cx-=rab*Dab.z*fa.x+rbc*Dbc.z*fc.x+rbc*Dbc.z*fb.x+rcd*Dcd.z*fd.x;

  strain_derivativen->ay-=rab*Dab.x*fa.y+rbc*Dbc.x*fc.y+rbc*Dbc.x*fb.y+rcd*Dcd.x*fd.y;
  strain_derivativen->by-=rab*Dab.y*fa.y+rbc*Dbc.y*fc.y+rbc*Dbc.y*fb.y+rcd*Dcd.y*fd.y;
  strain_derivativen->cy-=rab*Dab.z*fa.y+rbc*Dbc.z*fc.y+rbc*Dbc.z*fb.y+rcd*Dcd.z*fd.y;

  strain_derivativen->az-=rab*Dab.x*fa.z+rbc*Dbc.x*fc.z+rbc*Dbc.x*fb.z+rcd*Dcd.x*fd.z;
  strain_derivativen->bz-=rab*Dab.y*fa.z+rbc*Dbc.y*fc.z+rbc*Dbc.y*fb.z+rcd*Dcd.y*fd.z;
  strain_derivativen->cz-=rab*Dab.z*fa.z+rbc*Dbc.z*fc.z+rbc*Dbc.z*fb.z+rcd*Dcd.z*fd.z;

}

void ComputeFrameworkBendTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivativeTensor,int ComputeGradient,int ComputeHessian)
{
  int i,f1;
  int A,B,C,D;
  const REAL Delta=1e-7;
  VECTOR posA0,posB0,posC0,posD0;
  VECTOR posA,posB,posC,posD;
  VECTOR fa0,fb0,fc0,fd0,fa[4],fb[4],fc[4],fd[4];
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  REAL U;
  REAL_MATRIX3x3 strain_derivative;
  int n;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBendTorsions[f1];i++)
    {
      A=Framework[CurrentSystem].BendTorsions[f1][i].A;
      B=Framework[CurrentSystem].BendTorsions[f1][i].B;
      C=Framework[CurrentSystem].BendTorsions[f1][i].C;
      D=Framework[CurrentSystem].BendTorsions[f1][i].D;

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;
      index_k=Framework[CurrentSystem].Atoms[f1][C].HessianIndex;
      index_l=Framework[CurrentSystem].Atoms[f1][D].HessianIndex;

      posA0=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB0=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC0=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD0=Framework[CurrentSystem].Atoms[f1][D].Position;

      CalculateFrameworkBendTorsionForces(f1,i,posA0,posB0,posC0,posD0,&U,&fa0,&fb0,&fc0,&fd0,&strain_derivative);

      *Energy+=U;

      StrainDerivativeTensor->ax+=strain_derivative.ax;
      StrainDerivativeTensor->ay+=strain_derivative.ay;
      StrainDerivativeTensor->az+=strain_derivative.az;

      StrainDerivativeTensor->bx+=strain_derivative.bx;
      StrainDerivativeTensor->by+=strain_derivative.by;
      StrainDerivativeTensor->bz+=strain_derivative.bz;

      StrainDerivativeTensor->cx+=strain_derivative.cx;
      StrainDerivativeTensor->cy+=strain_derivative.cy;
      StrainDerivativeTensor->cz+=strain_derivative.cz;

      // add contribution to the first derivatives
      if(ComputeGradient)
      {
        if(index_i.x>=0) Gradient[index_i.x]-=fa0.x;
        if(index_i.y>=0) Gradient[index_i.y]-=fa0.y;
        if(index_i.z>=0) Gradient[index_i.z]-=fa0.z;

        if(index_j.x>=0) Gradient[index_j.x]-=fb0.x;
        if(index_j.y>=0) Gradient[index_j.y]-=fb0.y;
        if(index_j.z>=0) Gradient[index_j.z]-=fb0.z;

        if(index_k.x>=0) Gradient[index_k.x]-=fc0.x;
        if(index_k.y>=0) Gradient[index_k.y]-=fc0.y;
        if(index_k.z>=0) Gradient[index_k.z]-=fc0.z;

        if(index_l.x>=0) Gradient[index_l.x]-=fd0.x;
        if(index_l.y>=0) Gradient[index_l.y]-=fd0.y;
        if(index_l.z>=0) Gradient[index_l.z]-=fd0.z;

        n=NumberOfCoordinatesMinimizationVariables;
        switch(Ensemble[CurrentSystem])
        {
          case NPT:
            Gradient[n]+=strain_derivative.ax;
            break;
          case NPTPR:
          case NPHPR:
            switch(NPTPRCellType[CurrentSystem])
            {
              case ISOTROPIC:
                Gradient[n  ]+=(strain_derivative.ax+strain_derivative.by+strain_derivative.cz)/3.0; // xx
                Gradient[n+1]+=(strain_derivative.ax+strain_derivative.by+strain_derivative.cz)/3.0; // xx
                Gradient[n+2]+=(strain_derivative.ax+strain_derivative.by+strain_derivative.cz)/3.0; // xx
                break;
              case ANISOTROPIC:
                Gradient[n  ]+=strain_derivative.ax; // xx
                Gradient[n+1]+=strain_derivative.by; // yy
                Gradient[n+2]+=strain_derivative.cz; // zz
                break;
              case REGULAR:
              case REGULAR_UPPER_TRIANGLE:
                Gradient[n  ]+=strain_derivative.ax;  // xx
                Gradient[n+1]+=strain_derivative.ay;  // xy
                Gradient[n+2]+=strain_derivative.az;  // xz
                Gradient[n+3]+=strain_derivative.by;  // yy
                Gradient[n+4]+=strain_derivative.bz;  // yz
                Gradient[n+5]+=strain_derivative.cz;  // zz
                break;
              case MONOCLINIC:
              case MONOCLINIC_UPPER_TRIANGLE:
                switch(MonoclinicAngleType[CurrentSystem])
                {
                  case MONOCLINIC_ALPHA_ANGLE:
                    Gradient[n]+=strain_derivative.ax;
                    Gradient[n+1]+=strain_derivative.by;
                    Gradient[n+2]+=strain_derivative.bz;
                    Gradient[n+3]+=strain_derivative.cz;
                    break;
                  case MONOCLINIC_BETA_ANGLE:
                    Gradient[n]+=strain_derivative.ax;
                    Gradient[n+1]+=strain_derivative.az;
                    Gradient[n+2]+=strain_derivative.by;
                    Gradient[n+3]+=strain_derivative.cz;
                    break;
                  case MONOCLINIC_GAMMA_ANGLE:
                  default:
                    Gradient[n]+=strain_derivative.ax;
                    Gradient[n+1]+=strain_derivative.ay;
                    Gradient[n+2]+=strain_derivative.by;
                    Gradient[n+3]+=strain_derivative.cz;
                    break;
                }
                break;
              default:
                fprintf(stderr, "Unknown NPTPRCellType\n");
                exit(0);
                break;
            }
            break;
          case NVT:
          case NVE:
            break;
        }
      }

      if(ComputeHessian)
      {
        // Atom A

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_i.x>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_i.x][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_i.x][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_i.x][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_i.x][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_i.x][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_i.x][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_i.x][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_i.x][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_i.x][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_i.x][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_i.x][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_i.x][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_i.y>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_i.y][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_i.y][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_i.y][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_i.y][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_i.y][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_i.y][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_i.y][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_i.y][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_i.y][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_i.y][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_i.y][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_i.y][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_i.z>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_i.z][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_i.z][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_i.z][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_i.z][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_i.z][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_i.z][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_i.z][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_i.z][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_i.z][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_i.z][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_i.z][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_i.z][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        // Atom B

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_j.x>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_j.x][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_j.x][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_j.x][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_j.x][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_j.x][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_j.x][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_j.x][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_j.x][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_j.x][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_j.x][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_j.x][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_j.x][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_j.y>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_j.y][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_j.y][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_j.y][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_j.y][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_j.y][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_j.y][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_j.y][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_j.y][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_j.y][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_j.y][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_j.y][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_j.z>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_j.z][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_j.z][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_j.z][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_j.z][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_j.z][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_j.z][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_j.z][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_j.z][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_j.z][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_j.z][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_j.z][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        // Atom C

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_k.x>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_k.x][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_k.x][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_k.x][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_k.x][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_k.x][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_k.x][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_k.x][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_k.x][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_k.x][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_k.x][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_k.x][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_k.y>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_k.y][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_k.y][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_k.y][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_k.y][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_k.y][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_k.y][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_k.y][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_k.y][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_k.y][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_k.y][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_k.y][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_k.y][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_k.z>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_k.z][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_k.z][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_k.z][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_k.z][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_k.z][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_k.z][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_k.z][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_k.z][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_k.z][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_k.z][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_k.z][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_k.z][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        // Atom D

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_l.x>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_l.x][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_l.x][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_l.x][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_l.x][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_l.x][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_l.x][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_l.x][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_l.x][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_l.x][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_l.x][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_l.x][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_l.x][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_l.y>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_l.y][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_l.y][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_l.y][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_l.y][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_l.y][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_l.y][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_l.y][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_l.y][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_l.y][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_l.y][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_l.y][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_l.y][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z+=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z+=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z-=0.5*Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z-=Delta;
        CalculateFrameworkBendTorsionForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

        if(index_l.z>=0)
        {
          if(index_i.x>=0) HessianMatrix.element[index_l.z][index_i.x]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
          if(index_i.y>=0) HessianMatrix.element[index_l.z][index_i.y]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
          if(index_i.z>=0) HessianMatrix.element[index_l.z][index_i.z]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);

          if(index_j.x>=0) HessianMatrix.element[index_l.z][index_j.x]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
          if(index_j.y>=0) HessianMatrix.element[index_l.z][index_j.y]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
          if(index_j.z>=0) HessianMatrix.element[index_l.z][index_j.z]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);

          if(index_k.x>=0) HessianMatrix.element[index_l.z][index_k.x]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
          if(index_k.y>=0) HessianMatrix.element[index_l.z][index_k.y]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
          if(index_k.z>=0) HessianMatrix.element[index_l.z][index_k.z]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);

          if(index_l.x>=0) HessianMatrix.element[index_l.z][index_l.x]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
          if(index_l.y>=0) HessianMatrix.element[index_l.z][index_l.y]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
          if(index_l.z>=0) HessianMatrix.element[index_l.z][index_l.z]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
        }
      }
    }
  }
}
