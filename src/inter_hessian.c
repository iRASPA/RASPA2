/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'inter_hessian.c' is part of RASPA-2.0

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
#include <sys/stat.h>
#include "constants.h"
#include "utils.h"
#include "simulation.h"
#include "potentials.h"
#include "output.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "grids.h"
#include "ewald.h"
#include "inter_energy.h"
#include "inter_force.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "spectra.h"
#include "rigid.h"
#include "minimization.h"

// Hessian: Center of mass - Center of mass
// ========================================
void HessianAtomicPositionPosition(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,
                                                 REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor)
{
  REAL_MATRIX3x3 Hessian;

  // the first and second term of Eq. S54 of Ref. Dubbeldam, Krishna, Snurr, 2009
  // ============================================================================
  // I==J: f_2 r_ij,alpha r_ij,beta + f_1 delta(alpha,beta)
  // I!=J: -f_2 r_ij,alpha r_ij,beta - f_1 delta(alpha,beta)

  Hessian.ax=f2*dr.x*dr.x+f1; Hessian.bx=f2*dr.y*dr.x;    Hessian.cx=f2*dr.z*dr.x;
  Hessian.ay=f2*dr.x*dr.y;    Hessian.by=f2*dr.y*dr.y+f1; Hessian.cy=f2*dr.z*dr.y;
  Hessian.az=f2*dr.x*dr.z;    Hessian.bz=f2*dr.y*dr.z;    Hessian.cz=f2*dr.z*dr.z+f1;

  // case [I,I]: Center of mass - Center of mass
  switch(Dimension)
  {
    case 3:
      if(index_i.z>=0)
      {
        if(index_i.x>=0) HessianMatrix.element[index_i.x][index_i.z]+=ReplicaFactor*Hessian.az;
        if(index_i.y>=0) HessianMatrix.element[index_i.y][index_i.z]+=ReplicaFactor*Hessian.bz;
        if(index_i.z>=0) HessianMatrix.element[index_i.z][index_i.z]+=ReplicaFactor*Hessian.cz;
      }
    case 2:
      if(index_i.y>=0)
      {
        if(index_i.x>=0) HessianMatrix.element[index_i.x][index_i.y]+=ReplicaFactor*Hessian.ay;
        if(index_i.y>=0) HessianMatrix.element[index_i.y][index_i.y]+=ReplicaFactor*Hessian.by;
      }
    case 1:
      if(index_i.x>=0) HessianMatrix.element[index_i.x][index_i.x]+=ReplicaFactor*Hessian.ax;
      break;
  }


  // case [J,J]: Center of mass - Center of mass
  switch(Dimension)
  {
    case 3:
      if(index_j.z>=0)
      {
        if(index_j.x>=0) HessianMatrix.element[index_j.x][index_j.z]+=ReplicaFactor*Hessian.az;
        if(index_j.y>=0) HessianMatrix.element[index_j.y][index_j.z]+=ReplicaFactor*Hessian.bz;
        if(index_j.z>=0) HessianMatrix.element[index_j.z][index_j.z]+=ReplicaFactor*Hessian.cz;
      }
    case 2:
      if(index_j.y>=0)
      {
        if(index_j.x>=0) HessianMatrix.element[index_j.x][index_j.y]+=ReplicaFactor*Hessian.ay;
        if(index_j.y>=0) HessianMatrix.element[index_j.y][index_j.y]+=ReplicaFactor*Hessian.by;
      }
    case 1:
      if(index_j.x>=0) HessianMatrix.element[index_j.x][index_j.x]+=ReplicaFactor*Hessian.ax;
      break;
  }

  // case [I,J]: Center of mass - Center of mass
  switch(Dimension)
  {
    case 3:
      if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]-=Hessian.az;
      if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]-=Hessian.bz;
      if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]-=Hessian.az;
      if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]-=Hessian.bz;
      if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]-=Hessian.cz;
    case 2:
      if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]-=Hessian.ay;
      if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]-=Hessian.ay;
      if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]-=Hessian.by;
    case 1:
      if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]-=Hessian.ax;
      break;
  }
}


// Hessian: Center of mass - Orientation
// =====================================
void HessianCenterOfMassOrientation(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_i2,INT_VECTOR3 index_j,INT_VECTOR3 index_j2,
              int index1,int index2,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor)
{
  REAL_MATRIX3x3 Hessian;
  VECTOR veci1,veci2,veci3;
  VECTOR vecj1,vecj2,vecj3;

  // the first and second term of Eq. S56 of Ref. Dubbeldam, Krishna, Snurr, 2009
  // ============================================================================
  // I==J: f_2 (r_ij . R^I_beta r_i^0) r_ij,alpha + f1 [R^I_beta r_i^0]_alpha
  // I!=J: -f_2 (r_ij . R^J_beta r_j^0) r_ij,alpha - f1 [R^J_beta r_j^0]_alpha

  // case [I,I]: Center of mass I - Orientation J
  if(index1>=0)
  {
    veci1=DVecX[index1];
    veci2=DVecY[index1];
    veci3=DVecZ[index1];
  }

  if(index2>=0)
  {
    vecj1=DVecX[index2];
    vecj2=DVecY[index2];
    vecj3=DVecZ[index2];
  }

  Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
  Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
  Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

  Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
  Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
  Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

  Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
  Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
  Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

  switch(Dimension)
  {
    case 3:
      if((index_i.x>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i.x][index_i2.z]+=ReplicaFactor*Hessian.az;  // xz
      if((index_i.y>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i.y][index_i2.z]+=ReplicaFactor*Hessian.bz;  // yz
      if((index_i.z>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i.z][index_i2.x]+=ReplicaFactor*Hessian.cx;  // zx
      if((index_i.z>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i.z][index_i2.y]+=ReplicaFactor*Hessian.cy;  // zy
      if((index_i.z>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i.z][index_i2.z]+=ReplicaFactor*Hessian.cz;  // zz
    case 2:
      if((index_i.x>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i.x][index_i2.y]+=ReplicaFactor*Hessian.ay;  // xy
      if((index_i.y>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i.y][index_i2.x]+=ReplicaFactor*Hessian.bx;  // yx
      if((index_i.y>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i.y][index_i2.y]+=ReplicaFactor*Hessian.by;  // yy
    case 1:
      if((index_i.x>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i.x][index_i2.x]+=ReplicaFactor*Hessian.ax;  // xx
      break;
  }

  // case [I,J]: Orientation I - Center of mass J
  if(index1>=0)
  {
    veci1=DVecX[index1];
    veci2=DVecY[index1];
    veci3=DVecZ[index1];
  }

  if(index2>=0)
  {
    vecj1=DVecX[index2];
    vecj2=DVecY[index2];
    vecj3=DVecZ[index2];
  }

  Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
  Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
  Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

  Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
  Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
  Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

  Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
  Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
  Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

  switch(Dimension)
  {
    case 3:
      if((index_i2.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.z][index_j.x]-=Hessian.az;  // zx
      if((index_i2.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.z][index_j.y]-=Hessian.bz;  // zy
      if((index_i2.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.x][index_j.z]-=Hessian.cx;  // xz
      if((index_i2.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.y][index_j.z]-=Hessian.cy;  // yz
      if((index_i2.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.z][index_j.z]-=Hessian.cz;  // zz
    case 2:
      if((index_i2.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.y][index_j.x]-=Hessian.ay;  // yx
      if((index_i2.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.x][index_j.y]-=Hessian.bx;  // xy
      if((index_i2.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.y][index_j.y]-=Hessian.by;  // yy
    case 1:
      if((index_i2.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.x][index_j.x]-=Hessian.ax;  // xx
      break;
  }

  // case [J,J]: Center of mass J - Orientation J
  if(index1>=0)
  {
    veci1=DVecX[index1];
    veci2=DVecY[index1];
    veci3=DVecZ[index1];
  }

  if(index2>=0)
  {
    vecj1=DVecX[index2];
    vecj2=DVecY[index2];
    vecj3=DVecZ[index2];
  }

  Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
  Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
  Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

  Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
  Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
  Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

  Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
  Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
  Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

  switch(Dimension)
  {
    case 3:
      if((index_j.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_j.x][index_j2.z]+=ReplicaFactor*Hessian.az;
      if((index_j.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_j.y][index_j2.z]+=ReplicaFactor*Hessian.bz;
      if((index_j.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_j.z][index_j2.x]+=ReplicaFactor*Hessian.cx;
      if((index_j.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_j.z][index_j2.y]+=ReplicaFactor*Hessian.cy;
      if((index_j.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_j.z][index_j2.z]+=ReplicaFactor*Hessian.cz;
    case 2:
      if((index_j.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_j.x][index_j2.y]+=ReplicaFactor*Hessian.ay;
      if((index_j.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_j.y][index_j2.x]+=ReplicaFactor*Hessian.bx;
      if((index_j.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_j.y][index_j2.y]+=ReplicaFactor*Hessian.by;
    case 1:
      if((index_j.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_j.x][index_j2.x]+=ReplicaFactor*Hessian.ax;
      break;
  }

  // case [I,J]: Center of mass I - Orientation J
  if(index1>=0)
  {
    veci1=DVecX[index1];
    veci2=DVecY[index1];
    veci3=DVecZ[index1];
  }

  if(index2>=0)
  {
    vecj1=DVecX[index2];
    vecj2=DVecY[index2];
    vecj3=DVecZ[index2];
  }

  Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
  Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
  Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

  Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
  Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
  Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

  Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
  Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
  Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

  switch(Dimension)
  {
    case 3:
      if((index_i.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.x][index_j2.z]-=Hessian.az;  // xz
      if((index_i.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.y][index_j2.z]-=Hessian.bz;  // yz
      if((index_i.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.z][index_j2.x]-=Hessian.cx;  // zx
      if((index_i.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.z][index_j2.y]-=Hessian.cy;  // zy
      if((index_i.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.z][index_j2.z]-=Hessian.cz;  // zz
    case 2:
      if((index_i.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.x][index_j2.y]-=Hessian.ay;  // xy
      if((index_i.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.y][index_j2.x]-=Hessian.bx;  // yx
      if((index_i.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.y][index_j2.y]-=Hessian.by;  // yy
    case 1:
      if((index_i.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.x][index_j2.x]-=Hessian.ax;  // xx
      break;
  }
}


// Hessian: Orientation - Orientation
// ==================================
void HessianOrientationOrientation(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_i2,INT_VECTOR3 index_j,INT_VECTOR3 index_j2,
             int index1,int index2,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor)
{
  REAL_MATRIX3x3 Hessian;
  VECTOR veci1,veci2,veci3;
  VECTOR vecj1,vecj2,vecj3;
  VECTOR DDvecIAX,DDvecIBY,DDvecICZ,DDvecIAY,DDvecIAZ,DDvecIBZ;
  VECTOR DDvecJAX,DDvecJBY,DDvecJCZ,DDvecJAY,DDvecJAZ,DDvecJBZ;

  // the first, second and third terms of Eq. S55 of Ref. Dubbeldam, Krishna, Snurr, 2009
  // ====================================================================================
  // I==J: f_2 (r_ij . R^I_alpha r_i^0) (r_ij . R^I_beta r_i^0)+
  //       f_1 (R^I_alpha r_i^0).(R^I_beta r_i^0)+
  //       f_1 (r_ij . R^I_alpha,beta r_i^0)
  // I!=J: -f_2 (r_ij . R^I_alpha r_i^0) (r_ij . R^J_beta r_j^0)-
  //       f_1 (R^I_alpha r_i^0).(R^J_beta r_j^0)

  // case [I,I]: Orientation I - Orientation I
  veci1=DVecX[index1]; vecj1=DVecX[index2];
  veci2=DVecY[index1]; vecj2=DVecY[index2];
  veci3=DVecZ[index1]; vecj3=DVecZ[index2];

  DDvecIAX=DDVecAX[index1]; DDvecJAX=DDVecAX[index2];
  DDvecIBY=DDVecBY[index1]; DDvecJBY=DDVecBY[index2];
  DDvecICZ=DDVecCZ[index1]; DDvecJCZ=DDVecCZ[index2];
  DDvecIAY=DDVecAY[index1]; DDvecJAY=DDVecAY[index2];
  DDvecIAZ=DDVecAZ[index1]; DDvecJAZ=DDVecAZ[index2];
  DDvecIBZ=DDVecBZ[index1]; DDvecJBZ=DDVecBZ[index2];

  Hessian.ax=f1*(dr.x*DDvecIAX.x+dr.y*DDvecIAX.y+dr.z*DDvecIAX.z)+
             f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)+
             f1*(veci1.x*veci1.x+veci1.y*veci1.y+veci1.z*veci1.z);

  Hessian.by=f1*(dr.x*DDvecIBY.x+dr.y*DDvecIBY.y+dr.z*DDvecIBY.z)+
             f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)+
             f1*(veci2.x*veci2.x+veci2.y*veci2.y+veci2.z*veci2.z);

  Hessian.cz=f1*(dr.x*DDvecICZ.x+dr.y*DDvecICZ.y+dr.z*DDvecICZ.z)+
             f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)+
             f1*(veci3.x*veci3.x+veci3.y*veci3.y+veci3.z*veci3.z);


  Hessian.ay=f1*(dr.x*DDvecIAY.x+dr.y*DDvecIAY.y+dr.z*DDvecIAY.z)+
             f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)+
             f1*(veci1.x*veci2.x+veci1.y*veci2.y+veci1.z*veci2.z);

  Hessian.az=f1*(dr.x*DDvecIAZ.x+dr.y*DDvecIAZ.y+dr.z*DDvecIAZ.z)+
             f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)+
             f1*(veci1.x*veci3.x+veci1.y*veci3.y+veci1.z*veci3.z);

  Hessian.bz=f1*(dr.x*DDvecIBZ.x+dr.y*DDvecIBZ.y+dr.z*DDvecIBZ.z)+
             f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)+
             f1*(veci2.x*veci3.x+veci2.y*veci3.y+veci2.z*veci3.z);

  switch(Dimension)
  {
    case 3:
      if((index_i2.z>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.z][index_i2.z]+=ReplicaFactor*Hessian.cz;
      if((index_i2.x>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.x][index_i2.z]+=ReplicaFactor*Hessian.az;
      if((index_i2.y>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.y][index_i2.z]+=ReplicaFactor*Hessian.bz;
    case 2:
      if((index_i2.y>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i2.y][index_i2.y]+=ReplicaFactor*Hessian.by;
      if((index_i2.x>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i2.x][index_i2.y]+=ReplicaFactor*Hessian.ay;
    case 1:
      if((index_i2.x>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i2.x][index_i2.x]+=ReplicaFactor*Hessian.ax;
      break;
  }

  // case [J,J]: Orientation J - Orientation J
  veci1=DVecX[index1]; vecj1=DVecX[index2];
  veci2=DVecY[index1]; vecj2=DVecY[index2];
  veci3=DVecZ[index1]; vecj3=DVecZ[index2];

  DDvecIAX=DDVecAX[index1]; DDvecJAX=DDVecAX[index2];
  DDvecIBY=DDVecBY[index1]; DDvecJBY=DDVecBY[index2];
  DDvecICZ=DDVecCZ[index1]; DDvecJCZ=DDVecCZ[index2];
  DDvecIAY=DDVecAY[index1]; DDvecJAY=DDVecAY[index2];
  DDvecIAZ=DDVecAZ[index1]; DDvecJAZ=DDVecAZ[index2];
  DDvecIBZ=DDVecBZ[index1]; DDvecJBZ=DDVecBZ[index2];

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

  switch(Dimension)
  {
    case 3:
      if((index_j2.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_j2.z][index_j2.z]+=ReplicaFactor*Hessian.cz;
      if((index_j2.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_j2.x][index_j2.z]+=ReplicaFactor*Hessian.az;
      if((index_j2.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_j2.y][index_j2.z]+=ReplicaFactor*Hessian.bz;
    case 2:
      if((index_j2.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_j2.y][index_j2.y]+=ReplicaFactor*Hessian.by;
      if((index_j2.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_j2.x][index_j2.y]+=ReplicaFactor*Hessian.ay;
    case 1:
      if((index_j2.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_j2.x][index_j2.x]+=ReplicaFactor*Hessian.ax;
      break;
  }

  // case [I,J]: Orientation I - Orientation J
  veci1=DVecX[index1]; vecj1=DVecX[index2];
  veci2=DVecY[index1]; vecj2=DVecY[index2];
  veci3=DVecZ[index1]; vecj3=DVecZ[index2];

  DDvecIAX=DDVecAX[index1]; DDvecJAX=DDVecAX[index2];
  DDvecIBY=DDVecBY[index1]; DDvecJBY=DDVecBY[index2];
  DDvecICZ=DDVecCZ[index1]; DDvecJCZ=DDVecCZ[index2];
  DDvecIAY=DDVecAY[index1]; DDvecJAY=DDVecAY[index2];
  DDvecIAZ=DDVecAZ[index1]; DDvecJAZ=DDVecAZ[index2];
  DDvecIBZ=DDVecBZ[index1]; DDvecJBZ=DDVecBZ[index2];

  Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)+
             f1*(veci1.x*vecj1.x+veci1.y*vecj1.y+veci1.z*vecj1.z);
  Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)+
             f1*(veci2.x*vecj2.x+veci2.y*vecj2.y+veci2.z*vecj2.z);
  Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
             f1*(veci3.x*vecj3.x+veci3.y*vecj3.y+veci3.z*vecj3.z);
  Hessian.ay=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)+
             f1*(veci1.x*vecj2.x+veci1.y*vecj2.y+veci1.z*vecj2.z);
  Hessian.bx=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)+
             f1*(veci2.x*vecj1.x+veci2.y*vecj1.y+veci2.z*vecj1.z);
  Hessian.az=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
             f1*(veci1.x*vecj3.x+veci1.y*vecj3.y+veci1.z*vecj3.z);
  Hessian.cx=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)+
             f1*(veci3.x*vecj1.x+veci3.y*vecj1.y+veci3.z*vecj1.z);
  Hessian.bz=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
             f1*(veci2.x*vecj3.x+veci2.y*vecj3.y+veci2.z*vecj3.z);
  Hessian.cy=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)+
             f1*(veci3.x*vecj2.x+veci3.y*vecj2.y+veci3.z*vecj2.z);

  switch(Dimension)
  {
    case 3:
      if((index_i2.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.z][index_j2.z]-=Hessian.cz;
      if((index_i2.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.x][index_j2.z]-=Hessian.az;
      if((index_i2.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.z][index_j2.x]-=Hessian.cx;
      if((index_i2.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.y][index_j2.z]-=Hessian.bz;
      if((index_i2.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.z][index_j2.y]-=Hessian.cy;
    case 2:
      if((index_i2.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.y][index_j2.y]-=Hessian.by;
      if((index_i2.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.x][index_j2.y]-=Hessian.ay;
      if((index_i2.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.y][index_j2.x]-=Hessian.bx;
    case 1:
      if((index_i2.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.x][index_j2.x]-=Hessian.ax;
      break;
  }
}



// Hessian: Center of mass - Strain
// ================================
void HessianCenterOfMassStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,
      REAL f1,REAL f2,VECTOR dr,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB,int RigidA,int RigidB)
{
  int n;
  VECTOR dI,dJ;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // the first and second term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
      // ===========================================================================
      switch(Dimension)
      {
        case 3:
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;   // xx z + yy z + zz z
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
        case 2:
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;   // xx y + yy y + zz y
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;
        case 1:
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;   // xx x + yy x + zz x
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;
          break;
      }

      // the third and fourth term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
      // ===========================================================================
      if(RigidA)
      {
        dI.x=posA.x-comA.x;
        dI.y=posA.y-comA.y;
        dI.z=posA.z-comA.z;

        switch(Dimension)
        {
          case 3:
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=dI.x*(f2*dr.z*dr.x)+dI.y*(f2*dr.z*dr.y)+dI.z*(f2*dr.z*dr.z+f1); // xx z + yy z + zz z
            if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=dI.x*(f2*dr.z*dr.x)+dI.y*(f2*dr.z*dr.y)+dI.z*(f2*dr.z*dr.z+f1);
          case 2:
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=dI.x*(f2*dr.y*dr.x)+dI.y*(f2*dr.y*dr.y+f1)+dI.z*(f2*dr.y*dr.z); // xx y + yy y + zz y
            if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=dI.x*(f2*dr.y*dr.x)+dI.y*(f2*dr.y*dr.y+f1)+dI.z*(f2*dr.y*dr.z);
          case 1:
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=dI.x*(f2*dr.x*dr.x+f1)+dI.y*(f2*dr.x*dr.y)+dI.z*(f2*dr.x*dr.z); // xx x + yy x + zz x
            if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=dI.x*(f2*dr.x*dr.x+f1)+dI.y*(f2*dr.x*dr.y)+dI.z*(f2*dr.x*dr.z);
            break;
        }
      }
      if(RigidB)
      {
        dJ.x=posB.x-comB.x;
        dJ.y=posB.y-comB.y;
        dJ.z=posB.z-comB.z;

        switch(Dimension)
        {
          case 3:
            if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=dJ.x*(f2*dr.z*dr.x)+dJ.y*(f2*dr.z*dr.y)+dJ.z*(f2*dr.z*dr.z+f1); // xx z + yy z + zz z
            if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=dJ.x*(f2*dr.z*dr.x)+dJ.y*(f2*dr.z*dr.y)+dJ.z*(f2*dr.z*dr.z+f1);
          case 2:
            if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=dJ.x*(f2*dr.y*dr.x)+dJ.y*(f2*dr.y*dr.y+f1)+dJ.z*(f2*dr.y*dr.z); // xx y + yy y + zz y
            if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=dJ.x*(f2*dr.y*dr.x)+dJ.y*(f2*dr.y*dr.y+f1)+dJ.z*(f2*dr.y*dr.z);
          case 1:
            if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=dJ.x*(f2*dr.x*dr.x+f1)+dJ.y*(f2*dr.x*dr.y)+dJ.z*(f2*dr.x*dr.z); // xx x + yy x + zz x
            if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=dJ.x*(f2*dr.x*dr.x+f1)+dJ.y*(f2*dr.x*dr.y)+dJ.z*(f2*dr.x*dr.z);
            break;
        }
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          // the first and second term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ===========================================================================
          switch(Dimension)
          {
            case 3:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.z*dr.z*dr.x;             // zz x
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.z*dr.z*dr.y;             // zz y
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;               // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.z*dr.z*dr.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.z*dr.z*dr.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.y*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
            case 2:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;               // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.y*dr.y*dr.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
            case 1:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;   // xx x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              break;
          }

          // the third and fourth term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ===========================================================================
          if(RigidA)
          {
            dI.x=posA.x-comA.x;
            dI.y=posA.y-comA.y;
            dI.z=posA.z-comA.z;

            switch(Dimension)
            {
              case 3:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=dI.z*(f2*dr.x*dr.z);    // zz x
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=dI.z*(f2*dr.y*dr.z);    // zz y
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=dI.x*(f2*dr.z*dr.x);      // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=dI.y*(f2*dr.z*dr.y);    // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=dI.z*(f2*dr.z*dr.z+f1); // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=dI.z*(f2*dr.x*dr.z);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=dI.z*(f2*dr.y*dr.z);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=dI.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=dI.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=dI.z*(f2*dr.z*dr.z+f1);
              case 2:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=dI.y*(f2*dr.x*dr.y);    // yy x
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=dI.x*(f2*dr.y*dr.x);      // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=dI.y*(f2*dr.y*dr.y+f1); // yy y

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=dI.y*(f2*dr.x*dr.y);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=dI.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=dI.y*(f2*dr.y*dr.y+f1);
              case 1:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=dI.x*(f2*dr.x*dr.x+f1);   // xx x
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=dI.x*(f2*dr.x*dr.x+f1);
                break;
            }

          }
          if(RigidB)
          {
            dJ.x=posB.x-comB.x;
            dJ.y=posB.y-comB.y;
            dJ.z=posB.z-comB.z;

            switch(Dimension)
            {
              case 3:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=dJ.z*(f2*dr.x*dr.z);    // zz x
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=dJ.z*(f2*dr.y*dr.z);    // zz y
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=dJ.x*(f2*dr.z*dr.x);      // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=dJ.y*(f2*dr.z*dr.y);    // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=dJ.z*(f2*dr.z*dr.z+f1); // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=dJ.z*(f2*dr.x*dr.z);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=dJ.z*(f2*dr.y*dr.z);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=dJ.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=dJ.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=dJ.z*(f2*dr.z*dr.z+f1);
              case 2:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=dJ.y*(f2*dr.x*dr.y);    // yy x
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=dJ.x*(f2*dr.y*dr.x);      // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=dJ.y*(f2*dr.y*dr.y+f1); // yy y

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=dJ.y*(f2*dr.x*dr.y);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=dJ.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=dJ.y*(f2*dr.y*dr.y+f1);
              case 1:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=dJ.x*(f2*dr.x*dr.x+f1);   // xx x
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=dJ.x*(f2*dr.x*dr.x+f1);
                break;
            }
          }
          break;
        case REGULAR:
          // the first and second term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ===========================================================================
          switch(Dimension)
          {
            case 3:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=f2*dr.y*dr.z*dr.x;             // yz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=f2*dr.z*dr.z*dr.x;             // zz x
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.x*dr.z*dr.y;             // xz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=f2*dr.z*dr.z*dr.y;             // zz y
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;               // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.y*dr.z;             // xy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.x*dr.z*dr.x+f1*dr.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]-=f2*dr.y*dr.z*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]-=f2*dr.z*dr.z*dr.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.x*dr.z*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]-=f2*dr.y*dr.z*dr.y+f1*dr.z;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]-=f2*dr.z*dr.z*dr.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.x*dr.z*dr.z+f1*dr.x;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.y*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]-=f2*dr.y*dr.z*dr.z+f1*dr.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
            case 2:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;             // xy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+Dimension]+=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;                       // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;             // xy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+Dimension]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+Dimension]-=f2*dr.y*dr.y*dr.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+Dimension]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
            case 1:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;     // xx x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              break;
          }

          // the third and fourth term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ===========================================================================
          if(RigidA)
          {
            dI.x=posA.x-comA.x;
            dI.y=posA.y-comA.y;
            dI.z=posA.z-comA.z;

            switch(Dimension)
            {
              case 3:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=0.5*(dI.z*(f2*dr.x*dr.x+f1)+dI.x*(f2*dr.x*dr.z));  // xz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]-=0.5*(dI.z*(f2*dr.x*dr.y)+dI.y*(f2*dr.x*dr.z));     // yz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]-=dI.z*(f2*dr.x*dr.z);                               // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=0.5*(dI.z*(f2*dr.y*dr.x)+dI.x*(f2*dr.y*dr.z));     // xz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]-=0.5*(dI.z*(f2*dr.y*dr.y+f1)+dI.y*(f2*dr.y*dr.z));  // yz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]-=dI.z*(f2*dr.y*dr.z);                               // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=dI.x*(f2*dr.z*dr.x);                                 // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=0.5*(dI.y*(f2*dr.z*dr.x)+dI.x*(f2*dr.z*dr.y));     // xy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=0.5*(dI.z*(f2*dr.z*dr.x)+dI.x*(f2*dr.z*dr.z+f1));  // xz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=dI.y*(f2*dr.z*dr.y);                               // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]-=0.5*(dI.z*(f2*dr.z*dr.y)+dI.y*(f2*dr.z*dr.z+f1));  // yz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]-=dI.z*(f2*dr.z*dr.z+f1);                            // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=0.5*(dI.z*(f2*dr.x*dr.x+f1)+dI.x*(f2*dr.x*dr.z));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]+=0.5*(dI.z*(f2*dr.x*dr.y)+dI.y*(f2*dr.x*dr.z));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]+=dI.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=0.5*(dI.z*(f2*dr.y*dr.x)+dI.x*(f2*dr.y*dr.z));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]+=0.5*(dI.z*(f2*dr.y*dr.y+f1)+dI.y*(f2*dr.y*dr.z));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]+=dI.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=dI.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=0.5*(dI.y*(f2*dr.z*dr.x)+dI.x*(f2*dr.z*dr.y));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=0.5*(dI.z*(f2*dr.z*dr.x)+dI.x*(f2*dr.z*dr.z+f1));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=dI.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]+=0.5*(dI.z*(f2*dr.z*dr.y)+dI.y*(f2*dr.z*dr.z+f1));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]+=dI.z*(f2*dr.z*dr.z+f1);
              case 2:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=0.5*(dI.y*(f2*dr.x*dr.x+f1)+dI.x*(f2*dr.x*dr.y));  // xy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+Dimension]-=dI.y*(f2*dr.x*dr.y);                       // yy x
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=dI.x*(f2*dr.y*dr.x);                                 // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=0.5*(dI.y*(f2*dr.y*dr.x)+dI.x*(f2*dr.y*dr.y+f1));  // xy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+Dimension]-=dI.y*(f2*dr.y*dr.y+f1);                    // yy y

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=0.5*(dI.y*(f2*dr.x*dr.x+f1)+dI.x*(f2*dr.x*dr.y));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+Dimension]+=dI.y*(f2*dr.x*dr.y);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=dI.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=0.5*(dI.y*(f2*dr.y*dr.x)+dI.x*(f2*dr.y*dr.y+f1));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+Dimension]+=dI.y*(f2*dr.y*dr.y+f1);
              case 1:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=dI.x*(f2*dr.x*dr.x+f1);                            // xx x
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=dI.x*(f2*dr.x*dr.x+f1);
                break;
            }
          }
          if(RigidB)
          {
            dJ.x=posB.x-comB.x;
            dJ.y=posB.y-comB.y;
            dJ.z=posB.z-comB.z;

            switch(Dimension)
            {
              case 3:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=0.5*(dJ.z*(f2*dr.x*dr.x+f1)+dJ.x*(f2*dr.x*dr.z)); // xz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=0.5*(dJ.z*(f2*dr.x*dr.y)+dJ.y*(f2*dr.x*dr.z));    // yz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=dJ.z*(f2*dr.x*dr.z);                              // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=0.5*(dJ.z*(f2*dr.y*dr.x)+dJ.x*(f2*dr.y*dr.z));    // xz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=0.5*(dJ.z*(f2*dr.y*dr.y+f1)+dJ.y*(f2*dr.y*dr.z)); // yz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=dJ.z*(f2*dr.y*dr.z);                              // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=dJ.x*(f2*dr.z*dr.x);                                // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=0.5*(dJ.y*(f2*dr.z*dr.x)+dJ.x*(f2*dr.z*dr.y));    // xy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=0.5*(dJ.z*(f2*dr.z*dr.x)+dJ.x*(f2*dr.z*dr.z+f1)); // xz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=dJ.y*(f2*dr.z*dr.y);                              // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=0.5*(dJ.z*(f2*dr.z*dr.y)+dJ.y*(f2*dr.z*dr.z+f1)); // yz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=dJ.z*(f2*dr.z*dr.z+f1);                           // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=0.5*(dJ.z*(f2*dr.x*dr.x+f1)+dJ.x*(f2*dr.x*dr.z));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]-=0.5*(dJ.z*(f2*dr.x*dr.y)+dJ.y*(f2*dr.x*dr.z));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]-=dJ.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=0.5*(dJ.z*(f2*dr.y*dr.x)+dJ.x*(f2*dr.y*dr.z));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]-=0.5*(dJ.z*(f2*dr.y*dr.y+f1)+dJ.y*(f2*dr.y*dr.z));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]-=dJ.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=dJ.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=0.5*(dJ.y*(f2*dr.z*dr.x)+dJ.x*(f2*dr.z*dr.y));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=0.5*(dJ.z*(f2*dr.z*dr.x)+dJ.x*(f2*dr.z*dr.z+f1));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=dJ.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]-=0.5*(dJ.z*(f2*dr.z*dr.y)+dJ.y*(f2*dr.z*dr.z+f1));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]-=dJ.z*(f2*dr.z*dr.z+f1);
              case 2:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=0.5*(dJ.y*(f2*dr.x*dr.x+f1)+dJ.x*(f2*dr.x*dr.y)); // xy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+Dimension]+=dJ.y*(f2*dr.x*dr.y);                      // yy x
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=dJ.x*(f2*dr.y*dr.x);                                // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=0.5*(dJ.y*(f2*dr.y*dr.x)+dJ.x*(f2*dr.y*dr.y+f1)); // xy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+Dimension]+=dJ.y*(f2*dr.y*dr.y+f1);                   // yy y

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=0.5*(dJ.y*(f2*dr.x*dr.x+f1)+dJ.x*(f2*dr.x*dr.y));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+Dimension]-=dJ.y*(f2*dr.x*dr.y);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=dJ.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=0.5*(dJ.y*(f2*dr.y*dr.x)+dJ.x*(f2*dr.y*dr.y+f1));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+Dimension]-=dJ.y*(f2*dr.y*dr.y+f1);
              case 1:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=dJ.x*(f2*dr.x*dr.x+f1);                             // xx x
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=dJ.x*(f2*dr.x*dr.x+f1);
                break;
            }
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          // the first and second term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ===========================================================================
          switch(Dimension)
          {
            case 3:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=f2*dr.y*dr.z*dr.x;             // yz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.x*dr.z*dr.y;             // xz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;               // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.y*dr.z;             // xy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.x*dr.z*dr.x+f1*dr.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]-=f2*dr.y*dr.z*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]-=f2*dr.z*dr.z*dr.x;

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.x*dr.z*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]-=f2*dr.y*dr.z*dr.y+f1*dr.z;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]-=f2*dr.z*dr.z*dr.y;

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.x*dr.z*dr.z+f1*dr.x;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.y*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]-=f2*dr.y*dr.z*dr.z+f1*dr.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
            case 2:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;             // xy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+Dimension]+=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;                       // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;             // xy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+Dimension]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+Dimension]-=f2*dr.y*dr.y*dr.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+Dimension]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
            case 1:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;     // xx x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              break;
          }

          // the third and fourth term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ===========================================================================
          if(RigidA)
          {
            dI.x=posA.x-comA.x;
            dI.y=posA.y-comA.y;
            dI.z=posA.z-comA.z;

            switch(Dimension)
            {
              case 3:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=dI.z*(f2*dr.x*dr.x+f1);  // xz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]-=dI.z*(f2*dr.x*dr.y);     // yz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]-=dI.z*(f2*dr.x*dr.z);     // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=dI.z*(f2*dr.y*dr.x);     // xz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]-=dI.z*(f2*dr.y*dr.y+f1);  // yz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]-=dI.z*(f2*dr.y*dr.z);     // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=dI.x*(f2*dr.z*dr.x);       // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=dI.y*(f2*dr.z*dr.x);     // xy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=dI.z*(f2*dr.z*dr.x);     // xz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=dI.y*(f2*dr.z*dr.y);     // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]-=dI.z*(f2*dr.z*dr.y);     // yz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]-=dI.z*(f2*dr.z*dr.z+f1);  // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=dI.z*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]+=dI.z*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]+=dI.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=dI.z*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]+=dI.z*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]+=dI.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=dI.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=dI.y*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=dI.z*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=dI.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]+=dI.z*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]+=dI.z*(f2*dr.z*dr.z+f1);
              case 2:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=dI.y*(f2*dr.x*dr.x+f1);  // xy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=dI.y*(f2*dr.x*dr.y);     // yy x
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=dI.x*(f2*dr.y*dr.x);       // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=dI.y*(f2*dr.y*dr.x);     // xy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=dI.y*(f2*dr.y*dr.y+f1);  // yy y

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=dI.y*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=dI.y*(f2*dr.x*dr.y);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=dI.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=dI.y*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=dI.y*(f2*dr.y*dr.y+f1);
              case 1:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=dI.x*(f2*dr.x*dr.x+f1);    // xx x
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=dI.x*(f2*dr.x*dr.x+f1);
                break;
            }
          }
          if(RigidB)
          {
            dJ.x=posB.x-comB.x;
            dJ.y=posB.y-comB.y;
            dJ.z=posB.z-comB.z;

            switch(Dimension)
            {
              case 3:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=dJ.z*(f2*dr.x*dr.x+f1); // xz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=dJ.z*(f2*dr.x*dr.y);    // yz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=dJ.z*(f2*dr.x*dr.z);    // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=dJ.z*(f2*dr.y*dr.x);    // xz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=dJ.z*(f2*dr.y*dr.y+f1); // yz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=dJ.z*(f2*dr.y*dr.z);    // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=dJ.x*(f2*dr.z*dr.x);      // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=dJ.y*(f2*dr.z*dr.x);    // xy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=dJ.z*(f2*dr.z*dr.x);    // xz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=dJ.y*(f2*dr.z*dr.y);    // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=dJ.z*(f2*dr.z*dr.y);    // yz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=dJ.z*(f2*dr.z*dr.z+f1); // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=dJ.z*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]-=dJ.z*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]-=dJ.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=dJ.z*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]-=dJ.z*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]-=dJ.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=dJ.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=dJ.y*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=dJ.z*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=dJ.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]-=dJ.z*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]-=dJ.z*(f2*dr.z*dr.z+f1);
              case 2:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=dJ.y*(f2*dr.x*dr.x+f1);   // xy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=dJ.y*(f2*dr.x*dr.y);      // yy x
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=dJ.x*(f2*dr.y*dr.x);      // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=dJ.y*(f2*dr.y*dr.x);    // xy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=dJ.y*(f2*dr.y*dr.y+f1); // yy y

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=dJ.y*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=dJ.y*(f2*dr.x*dr.y);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=dJ.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=dJ.y*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=dJ.y*(f2*dr.y*dr.y+f1);
              case 1:
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=dJ.x*(f2*dr.x*dr.x+f1);     // xx x
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=dJ.x*(f2*dr.x*dr.x+f1);
                break;
            }
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first and second term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;     // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.y*dr.y*dr.x;               // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.z*dr.x;               // yz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;               // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;               // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;               // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.y*dr.y*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.y*dr.z*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=f2*dr.z*dr.z*dr.x;

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.y*dr.z*dr.y+f1*dr.z;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=f2*dr.z*dr.z*dr.y;

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.y*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.y*dr.z*dr.z+f1*dr.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;

              // the third and fourth term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(RigidA)
              {
                dI.x=posA.x-comA.x;
                dI.y=posA.y-comA.y;
                dI.z=posA.z-comA.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=dI.x*(f2*dr.x*dr.x+f1);                              // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=dI.y*(f2*dr.x*dr.y);                               // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=0.5*(dI.z*(f2*dr.x*dr.y)+dI.y*(f2*dr.x*dr.z));     // yz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=dI.z*(f2*dr.x*dr.z);                               // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=dI.x*(f2*dr.y*dr.x);                                 // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=dI.y*(f2*dr.y*dr.y+f1);                            // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=0.5*(dI.z*(f2*dr.y*dr.y+f1)+dI.y*(f2*dr.y*dr.z));  // yz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=dI.z*(f2*dr.y*dr.z);                               // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=dI.x*(f2*dr.z*dr.x);                                 // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=dI.y*(f2*dr.z*dr.y);                               // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=0.5*(dI.z*(f2*dr.z*dr.y)+dI.y*(f2*dr.z*dr.z+f1));  // yz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=dI.z*(f2*dr.z*dr.z+f1);                            // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=dI.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=dI.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=0.5*(dI.z*(f2*dr.x*dr.y)+dI.y*(f2*dr.x*dr.z));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=dI.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=dI.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=dI.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=0.5*(dI.z*(f2*dr.y*dr.y+f1)+dI.y*(f2*dr.y*dr.z));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=dI.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=dI.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=dI.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=0.5*(dI.z*(f2*dr.z*dr.y)+dI.y*(f2*dr.z*dr.z+f1));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=dI.z*(f2*dr.z*dr.z+f1);
              }
              if(RigidB)
              {
                dJ.x=posB.x-comB.x;
                dJ.y=posB.y-comB.y;
                dJ.z=posB.z-comB.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=dJ.x*(f2*dr.x*dr.x+f1);                             // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=dJ.y*(f2*dr.x*dr.y);                              // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=0.5*(dJ.z*(f2*dr.x*dr.y)+dJ.y*(f2*dr.x*dr.z));    // yz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=dJ.z*(f2*dr.x*dr.z);                              // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=dJ.x*(f2*dr.y*dr.x);                                // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=dJ.y*(f2*dr.y*dr.y+f1);                           // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=0.5*(dJ.z*(f2*dr.y*dr.y+f1)+dJ.y*(f2*dr.y*dr.z)); // yz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=dJ.z*(f2*dr.y*dr.z);                              // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=dJ.x*(f2*dr.z*dr.x);                                // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=dJ.y*(f2*dr.z*dr.y);                              // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=0.5*(dJ.z*(f2*dr.z*dr.y)+dJ.y*(f2*dr.z*dr.z+f1)); // yz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=dJ.z*(f2*dr.z*dr.z+f1);                           // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=dJ.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=dJ.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=0.5*(dJ.z*(f2*dr.x*dr.y)+dJ.y*(f2*dr.x*dr.z));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=dJ.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=dJ.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=dJ.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=0.5*(dJ.z*(f2*dr.y*dr.y+f1)+dJ.y*(f2*dr.y*dr.z));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=dJ.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=dJ.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=dJ.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=0.5*(dJ.z*(f2*dr.z*dr.y)+dJ.y*(f2*dr.z*dr.z+f1));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=dJ.z*(f2*dr.z*dr.z+f1);
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the first and second term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;   // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;               // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.z*dr.y;             // xz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;               // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.z*dr.x+f1*dr.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.y*dr.y*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=f2*dr.z*dr.z*dr.x;

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.z*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=f2*dr.z*dr.z*dr.y;

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.z*dr.z+f1*dr.x;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.y*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;

              // the third and fourth term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(RigidA)
              {
                dI.x=posA.x-comA.x;
                dI.y=posA.y-comA.y;
                dI.z=posA.z-comA.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=dI.x*(f2*dr.x*dr.x+f1);                              // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=0.5*(dI.z*(f2*dr.x*dr.x+f1)+dI.x*(f2*dr.x*dr.z));  // xz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=dI.y*(f2*dr.x*dr.y);                               // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=dI.z*(f2*dr.x*dr.z);                               // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=dI.x*(f2*dr.y*dr.x);                                 // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=0.5*(dI.z*(f2*dr.y*dr.x)+dI.x*(f2*dr.y*dr.z));     // xz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=dI.y*(f2*dr.y*dr.y+f1);                            // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=dI.z*(f2*dr.y*dr.z);                               // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=dI.x*(f2*dr.z*dr.x);                                 // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=0.5*(dI.z*(f2*dr.z*dr.x)+dI.x*(f2*dr.z*dr.z+f1));  // xz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=dI.y*(f2*dr.z*dr.y);                               // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=dI.z*(f2*dr.z*dr.z+f1);                            // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=dI.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=0.5*(dI.z*(f2*dr.x*dr.x+f1)+dI.x*(f2*dr.x*dr.z));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=dI.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=dI.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=dI.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=0.5*(dI.z*(f2*dr.y*dr.x)+dI.x*(f2*dr.y*dr.z));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=dI.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=dI.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=dI.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=0.5*(dI.z*(f2*dr.z*dr.x)+dI.x*(f2*dr.z*dr.z+f1));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=dI.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=dI.z*(f2*dr.z*dr.z+f1);
              }
              if(RigidB)
              {
                dJ.x=posB.x-comB.x;
                dJ.y=posB.y-comB.y;
                dJ.z=posB.z-comB.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=dJ.x*(f2*dr.x*dr.x+f1);                             // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=0.5*(dJ.z*(f2*dr.x*dr.x+f1)+dJ.x*(f2*dr.x*dr.z)); // xz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=dJ.y*(f2*dr.x*dr.y);                              // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=dJ.z*(f2*dr.x*dr.z);                              // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=dJ.x*(f2*dr.y*dr.x);                                // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=0.5*(dJ.z*(f2*dr.y*dr.x)+dJ.x*(f2*dr.y*dr.z));    // xz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=dJ.y*(f2*dr.y*dr.y+f1);                           // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=dJ.z*(f2*dr.y*dr.z);                              // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=dJ.x*(f2*dr.z*dr.x);                                // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=0.5*(dJ.z*(f2*dr.z*dr.x)+dJ.x*(f2*dr.z*dr.z+f1)); // xz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=dJ.y*(f2*dr.z*dr.y);                              // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=dJ.z*(f2*dr.z*dr.z+f1);                           // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=dJ.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=0.5*(dJ.z*(f2*dr.x*dr.x+f1)+dJ.x*(f2*dr.x*dr.z));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=dJ.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=dJ.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=dJ.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=0.5*(dJ.z*(f2*dr.y*dr.x)+dJ.x*(f2*dr.y*dr.z));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=dJ.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=dJ.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=dJ.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=0.5*(dJ.z*(f2*dr.z*dr.x)+dJ.x*(f2*dr.z*dr.z+f1));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=dJ.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=dJ.z*(f2*dr.z*dr.z+f1);
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first and second term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;   // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;     // xy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;               // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;     // xy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;               // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.y*dr.z;             // xy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

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

              // the third and fourth term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(RigidA)
              {
                dI.x=posA.x-comA.x;
                dI.y=posA.y-comA.y;
                dI.z=posA.z-comA.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=dI.x*(f2*dr.x*dr.x+f1);                              // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=0.5*(dI.y*(f2*dr.x*dr.x+f1)+dI.x*(f2*dr.x*dr.y));  // xy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=dI.y*(f2*dr.x*dr.y);                               // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=dI.z*(f2*dr.x*dr.z);                               // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=dI.x*(f2*dr.y*dr.x);                                 // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=0.5*(dI.y*(f2*dr.y*dr.x)+dI.x*(f2*dr.y*dr.y+f1));  // xy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=dI.y*(f2*dr.y*dr.y+f1);                            // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=dI.z*(f2*dr.y*dr.z);                               // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=dI.x*(f2*dr.z*dr.x);                                 // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=0.5*(dI.y*(f2*dr.z*dr.x)+dI.x*(f2*dr.z*dr.y));     // xy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=dI.y*(f2*dr.z*dr.y);                               // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=dI.z*(f2*dr.z*dr.z+f1);                            // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=dI.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=0.5*(dI.y*(f2*dr.x*dr.x+f1)+dI.x*(f2*dr.x*dr.y));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=dI.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=dI.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=dI.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=0.5*(dI.y*(f2*dr.y*dr.x)+dI.x*(f2*dr.y*dr.y+f1));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=dI.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=dI.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=dI.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=0.5*(dI.y*(f2*dr.z*dr.x)+dI.x*(f2*dr.z*dr.y));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=dI.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=dI.z*(f2*dr.z*dr.z+f1);
              }
              if(RigidB)
              {
                dJ.x=posB.x-comB.x;
                dJ.y=posB.y-comB.y;
                dJ.z=posB.z-comB.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=dJ.x*(f2*dr.x*dr.x+f1);                             // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=0.5*(dJ.y*(f2*dr.x*dr.x+f1)+dJ.x*(f2*dr.x*dr.y)); // xy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=dJ.y*(f2*dr.x*dr.y);                              // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=dJ.z*(f2*dr.x*dr.z);                              // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=dJ.x*(f2*dr.y*dr.x);                                // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=0.5*(dJ.y*(f2*dr.y*dr.x)+dJ.x*(f2*dr.y*dr.y+f1)); // xy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=dJ.y*(f2*dr.y*dr.y+f1);                           // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=dJ.z*(f2*dr.y*dr.z);                              // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=dJ.x*(f2*dr.z*dr.x);                                // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=0.5*(dJ.y*(f2*dr.z*dr.x)+dJ.x*(f2*dr.z*dr.y));    // xy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=dJ.y*(f2*dr.z*dr.y);                              // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=dJ.z*(f2*dr.z*dr.z+f1);                           // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=dJ.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=0.5*(dJ.y*(f2*dr.x*dr.x+f1)+dJ.x*(f2*dr.x*dr.y));
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=dJ.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=dJ.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=dJ.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=0.5*(dJ.y*(f2*dr.y*dr.x)+dJ.x*(f2*dr.y*dr.y+f1));
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=dJ.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=dJ.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=dJ.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=0.5*(dJ.y*(f2*dr.z*dr.x)+dJ.x*(f2*dr.z*dr.y));
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=dJ.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=dJ.z*(f2*dr.z*dr.z+f1);
              }
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first and second term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;   // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.z*dr.x;             // yz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;               // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;               // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.y*dr.y*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.y*dr.z*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=f2*dr.z*dr.z*dr.x;

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.y*dr.z*dr.y+f1*dr.z;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=f2*dr.z*dr.z*dr.y;

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.y*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.y*dr.z*dr.z+f1*dr.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;

              // the third and fourth term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(RigidA)
              {
                dI.x=posA.x-comA.x;
                dI.y=posA.y-comA.y;
                dI.z=posA.z-comA.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=dI.x*(f2*dr.x*dr.x+f1);    // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=dI.y*(f2*dr.x*dr.y);     // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=dI.z*(f2*dr.x*dr.y);     // yz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=dI.z*(f2*dr.x*dr.z);     // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=dI.x*(f2*dr.y*dr.x);       // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=dI.y*(f2*dr.y*dr.y+f1);  // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=dI.z*(f2*dr.y*dr.y+f1);  // yz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=dI.z*(f2*dr.y*dr.z);     // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=dI.x*(f2*dr.z*dr.x);       // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=dI.y*(f2*dr.z*dr.y);     // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=dI.z*(f2*dr.z*dr.y);     // yz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=dI.z*(f2*dr.z*dr.z+f1);  // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=dI.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=dI.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=dI.z*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=dI.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=dI.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=dI.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=dI.z*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=dI.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=dI.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=dI.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=dI.z*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=dI.z*(f2*dr.z*dr.z+f1);
              }
              if(RigidB)
              {
                dJ.x=posB.x-comB.x;
                dJ.y=posB.y-comB.y;
                dJ.z=posB.z-comB.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=dJ.x*(f2*dr.x*dr.x+f1);   // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=dJ.y*(f2*dr.x*dr.y);    // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=dJ.z*(f2*dr.x*dr.y);    // yz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=dJ.z*(f2*dr.x*dr.z);    // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=dJ.x*(f2*dr.y*dr.x);      // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=dJ.y*(f2*dr.y*dr.y+f1); // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=dJ.z*(f2*dr.y*dr.y+f1); // yz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=dJ.z*(f2*dr.y*dr.z);    // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=dJ.x*(f2*dr.z*dr.x);      // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=dJ.y*(f2*dr.z*dr.y);    // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=dJ.z*(f2*dr.z*dr.y);    // yz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=dJ.z*(f2*dr.z*dr.z+f1); // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=dJ.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=dJ.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=dJ.z*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=dJ.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=dJ.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=dJ.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=dJ.z*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=dJ.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=dJ.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=dJ.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=dJ.z*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=dJ.z*(f2*dr.z*dr.z+f1);
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the first and second term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;     // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.z*dr.x+f1*dr.z;       // xz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.y*dr.x;               // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;               // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;               // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.z*dr.y;             // xz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;               // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.z*dr.x+f1*dr.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.y*dr.y*dr.x;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=f2*dr.z*dr.z*dr.x;

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.z*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=f2*dr.z*dr.z*dr.y;

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.z*dr.z+f1*dr.x;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.y*dr.y*dr.z;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;

              // the third and fourth term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(RigidA)
              {
                dI.x=posA.x-comA.x;
                dI.y=posA.y-comA.y;
                dI.z=posA.z-comA.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=dI.x*(f2*dr.x*dr.x+f1);    // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=dI.z*(f2*dr.x*dr.x+f1);  // xz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=dI.y*(f2*dr.x*dr.y);     // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=dI.z*(f2*dr.x*dr.z);     // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=dI.x*(f2*dr.y*dr.x);       // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=dI.z*(f2*dr.y*dr.x);     // xz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=dI.y*(f2*dr.y*dr.y+f1);  // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=dI.z*(f2*dr.y*dr.z);     // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=dI.x*(f2*dr.z*dr.x);       // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=dI.z*(f2*dr.z*dr.x);     // xz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=dI.y*(f2*dr.z*dr.y);     // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=dI.z*(f2*dr.z*dr.z+f1);  // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=dI.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=dI.z*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=dI.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=dI.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=dI.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=dI.z*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=dI.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=dI.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=dI.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=dI.z*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=dI.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=dI.z*(f2*dr.z*dr.z+f1);
              }
              if(RigidB)
              {
                dJ.x=posB.x-comB.x;
                dJ.y=posB.y-comB.y;
                dJ.z=posB.z-comB.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=dJ.x*(f2*dr.x*dr.x+f1);   // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=dJ.z*(f2*dr.x*dr.x+f1); // xz x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=dJ.y*(f2*dr.x*dr.y);    // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=dJ.z*(f2*dr.x*dr.z);    // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=dJ.x*(f2*dr.y*dr.x);      // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=dJ.z*(f2*dr.y*dr.x);    // xz y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=dJ.y*(f2*dr.y*dr.y+f1); // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=dJ.z*(f2*dr.y*dr.z);    // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=dJ.x*(f2*dr.z*dr.x);      // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=dJ.z*(f2*dr.z*dr.x);    // xz z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=dJ.y*(f2*dr.z*dr.y);    // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=dJ.z*(f2*dr.z*dr.z+f1); // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=dJ.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=dJ.z*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=dJ.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=dJ.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=dJ.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=dJ.z*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=dJ.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=dJ.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=dJ.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=dJ.z*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=dJ.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=dJ.z*(f2*dr.z*dr.z+f1);
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first and second term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;   // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;     // xy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;               // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;     // xy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;               // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.y*dr.z;             // xy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

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

              // the third and fourth term of Eq. 39 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===========================================================================
              if(RigidA)
              {
                dI.x=posA.x-comA.x;
                dI.y=posA.y-comA.y;
                dI.z=posA.z-comA.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=dI.x*(f2*dr.x*dr.x+f1);    // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=dI.y*(f2*dr.x*dr.x+f1);  // xy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=dI.y*(f2*dr.x*dr.y);     // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=dI.z*(f2*dr.x*dr.z);     // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=dI.x*(f2*dr.y*dr.x);       // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=dI.y*(f2*dr.y*dr.x);     // xy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=dI.y*(f2*dr.y*dr.y+f1);  // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=dI.z*(f2*dr.y*dr.z);     // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=dI.x*(f2*dr.z*dr.x);       // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=dI.y*(f2*dr.z*dr.x);     // xy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=dI.y*(f2*dr.z*dr.y);     // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=dI.z*(f2*dr.z*dr.z+f1);  // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=dI.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=dI.y*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=dI.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=dI.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=dI.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=dI.y*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=dI.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=dI.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=dI.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=dI.y*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=dI.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=dI.z*(f2*dr.z*dr.z+f1);
              }
              if(RigidB)
              {
                dJ.x=posB.x-comB.x;
                dJ.y=posB.y-comB.y;
                dJ.z=posB.z-comB.z;

                if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=dJ.x*(f2*dr.x*dr.x+f1);   // xx x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=dJ.y*(f2*dr.x*dr.x+f1); // xy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=dJ.y*(f2*dr.x*dr.y);    // yy x
                if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=dJ.z*(f2*dr.x*dr.z);    // zz x

                if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=dJ.x*(f2*dr.y*dr.x);      // xx y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=dJ.y*(f2*dr.y*dr.x);    // xy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=dJ.y*(f2*dr.y*dr.y+f1); // yy y
                if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=dJ.z*(f2*dr.y*dr.z);    // zz y

                if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=dJ.x*(f2*dr.z*dr.x);      // xx z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=dJ.y*(f2*dr.z*dr.x);    // xy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=dJ.y*(f2*dr.z*dr.y);    // yy z
                if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=dJ.z*(f2*dr.z*dr.z+f1); // zz z

                if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=dJ.x*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=dJ.y*(f2*dr.x*dr.x+f1);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=dJ.y*(f2*dr.x*dr.y);
                if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=dJ.z*(f2*dr.x*dr.z);

                if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=dJ.x*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=dJ.y*(f2*dr.y*dr.x);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=dJ.y*(f2*dr.y*dr.y+f1);
                if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=dJ.z*(f2*dr.y*dr.z);

                if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=dJ.x*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=dJ.y*(f2*dr.z*dr.x);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=dJ.y*(f2*dr.z*dr.y);
                if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=dJ.z*(f2*dr.z*dr.z+f1);
              }
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


// Hessian: Orientation - Strain
// =============================
void HessianOrientationStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i2,INT_VECTOR3 index_j2,int index1,int index2,
                  REAL f1,REAL f2,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB,VECTOR dr)
{
  int n;
  REAL_MATRIX3x3 Hessian;
  VECTOR drJI;
  VECTOR veci1,veci2,veci3;
  VECTOR vecj1,vecj2,vecj3;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // Eq. 40 of Ref. Dubbeldam, Krishna, Snurr, 2009
      // ==============================================
      veci1=DVecX[index1]; vecj1=DVecX[index2];
      veci2=DVecY[index1]; vecj2=DVecY[index2];
      veci3=DVecZ[index1]; vecj3=DVecZ[index2];

      drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
      drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
      drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

      Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
      Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
      Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

      Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
      Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
      Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

      Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
      Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
      Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

      switch(Dimension)
      {
        case 3:
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.az+drJI.y*Hessian.bz+drJI.z*Hessian.cz; // (xx z + yy z + zz z)/3
        case 2:
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=drJI.x*Hessian.ay+drJI.y*Hessian.by+drJI.z*Hessian.cy; // (xx y + yy y + zz y)/3
        case 1:
          if(index_i2.x>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.ax+drJI.y*Hessian.bx+drJI.z*Hessian.cx;  // (xx x + yy x + zz x)/3
          break;
      }

      veci1=DVecX[index1]; vecj1=DVecX[index2];
      veci2=DVecY[index1]; vecj2=DVecY[index2];
      veci3=DVecZ[index1]; vecj3=DVecZ[index2];

      drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
      drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
      drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

      Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
      Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
      Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

      Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
      Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
      Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

      Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
      Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
      Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

      switch(Dimension)
      {
        case 3:
          if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=drJI.x*Hessian.az+drJI.y*Hessian.bz+drJI.z*Hessian.cz; // (xx z + yy z + zz z)/3
        case 2:
          if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=drJI.x*Hessian.ay+drJI.y*Hessian.by+drJI.z*Hessian.cy; // (xx y + yy y + zz y)/3
        case 1:
          if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=drJI.x*Hessian.ax+drJI.y*Hessian.bx+drJI.z*Hessian.cx;   // (xx x + yy x + zz x)/3
          break;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          // Eq. 40 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ==============================================
          veci1=DVecX[index1]; vecj1=DVecX[index2];
          veci2=DVecY[index1]; vecj2=DVecY[index2];
          veci3=DVecZ[index1]; vecj3=DVecZ[index2];

          drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
          drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
          drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

          Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
          Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
          Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

          Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
          Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
          Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

          Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
          Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
          Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

          switch(Dimension)
          {
            case 3:
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=drJI.z*Hessian.cx;       // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=drJI.z*Hessian.cy;       // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.az;         // xx z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=drJI.y*Hessian.bz;       // yy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=drJI.z*Hessian.cz;       // zz z
            case 2:
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=drJI.y*Hessian.bx;       // yy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=drJI.x*Hessian.ay;         // xx y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=drJI.y*Hessian.by;       // yy y
            case 1:
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=drJI.x*Hessian.ax;         // xx x
              break;
          }

          veci1=DVecX[index1]; vecj1=DVecX[index2];
          veci2=DVecY[index1]; vecj2=DVecY[index2];
          veci3=DVecZ[index1]; vecj3=DVecZ[index2];

          drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
          drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
          drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

          Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
          Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
          Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

          Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
          Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
          Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

          Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
          Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
          Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

          switch(Dimension)
          {
            case 3:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=drJI.z*Hessian.cx;       // zz x
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=drJI.z*Hessian.cy;       // zz y
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=drJI.x*Hessian.az;         // xx z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=drJI.y*Hessian.bz;       // yy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=drJI.z*Hessian.cz;       // zz z
            case 2:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=drJI.y*Hessian.bx;       // yy x
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=drJI.x*Hessian.ay;         // xx y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=drJI.y*Hessian.by;       // yy y
            case 1:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=drJI.x*Hessian.ax;         // xx x
              break;
          }
          break;
        case REGULAR:
          // Eq. 40 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ==============================================
          veci1=DVecX[index1]; vecj1=DVecX[index2];
          veci2=DVecY[index1]; vecj2=DVecY[index2];
          veci3=DVecZ[index1]; vecj3=DVecZ[index2];

          drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
          drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
          drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

          Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
          Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
          Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

          Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
          Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
          Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

          Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
          Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
          Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

          switch(Dimension)
          {
            case 3:
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=0.5*(drJI.z*Hessian.ax+drJI.x*Hessian.cx); // xz x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+4]+=0.5*(drJI.z*Hessian.bx+drJI.y*Hessian.cx); // yz x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+5]+=drJI.z*Hessian.cx;                         // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=0.5*(drJI.z*Hessian.ay+drJI.x*Hessian.cy); // xz y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+4]+=0.5*(drJI.z*Hessian.by+drJI.y*Hessian.cy); // yz y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+5]+=drJI.z*Hessian.cy;                         // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.az;                           // xx z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=0.5*(drJI.y*Hessian.az+drJI.x*Hessian.bz); // xy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=0.5*(drJI.z*Hessian.az+drJI.x*Hessian.cz); // xz z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=drJI.y*Hessian.bz;                         // yy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+4]+=0.5*(drJI.z*Hessian.bz+drJI.y*Hessian.cz); // yz z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+5]+=drJI.z*Hessian.cz;                         // zz z
            case 2:
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=0.5*(drJI.y*Hessian.ax+drJI.x*Hessian.bx); // xy x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+Dimension]+=drJI.y*Hessian.bx;                 // yy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=drJI.x*Hessian.ay;                           // xx y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=0.5*(drJI.y*Hessian.ay+drJI.x*Hessian.by); // xy y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+Dimension]+=drJI.y*Hessian.by;                 // yy y
            case 1:
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=drJI.x*Hessian.ax;                           // xx x
              break;
          }

          veci1=DVecX[index1]; vecj1=DVecX[index2];
          veci2=DVecY[index1]; vecj2=DVecY[index2];
          veci3=DVecZ[index1]; vecj3=DVecZ[index2];

          drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
          drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
          drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

          Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
          Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
          Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

          Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
          Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
          Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

          Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
          Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
          Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

          switch(Dimension)
          {
            case 3:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=0.5*(drJI.z*Hessian.ax+drJI.x*Hessian.cx); // xz x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+4]-=0.5*(drJI.z*Hessian.bx+drJI.y*Hessian.cx); // yz x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+5]-=drJI.z*Hessian.cx;                         // zz x
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=0.5*(drJI.z*Hessian.ay+drJI.x*Hessian.cy); // xz y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+4]-=0.5*(drJI.z*Hessian.by+drJI.y*Hessian.cy); // yz y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+5]-=drJI.z*Hessian.cy;                         // zz y
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=drJI.x*Hessian.az;                           // xx z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=0.5*(drJI.y*Hessian.az+drJI.x*Hessian.bz); // xy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=0.5*(drJI.z*Hessian.az+drJI.x*Hessian.cz); // xz z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=drJI.y*Hessian.bz;                         // yy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+4]-=0.5*(drJI.z*Hessian.bz+drJI.y*Hessian.cz); // yz z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+5]-=drJI.z*Hessian.cz;                         // zz z
            case 2:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=0.5*(drJI.y*Hessian.ax+drJI.x*Hessian.bx); // xy x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+Dimension]-=drJI.y*Hessian.bx;                 // yy x
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=drJI.x*Hessian.ay;                           // xx y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=0.5*(drJI.y*Hessian.ay+drJI.x*Hessian.by); // xy y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+Dimension]-=drJI.y*Hessian.by;                 // yy y
            case 1:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=drJI.x*Hessian.ax;                           // xx x
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          // Eq. 40 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ==============================================
          veci1=DVecX[index1]; vecj1=DVecX[index2];
          veci2=DVecY[index1]; vecj2=DVecY[index2];
          veci3=DVecZ[index1]; vecj3=DVecZ[index2];

          drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
          drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
          drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

          Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
          Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
          Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

          Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
          Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
          Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

          Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
          Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
          Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

          switch(Dimension)
          {
            case 3:
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=drJI.z*Hessian.ax; // xz x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+4]+=drJI.z*Hessian.bx; // yz x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+5]+=drJI.z*Hessian.cx; // zz x

              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=drJI.z*Hessian.ay; // xz y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+4]+=drJI.z*Hessian.by; // yz y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+5]+=drJI.z*Hessian.cy; // zz y

              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.az;   // xx z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=drJI.y*Hessian.az; // xy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=drJI.z*Hessian.az; // xz z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=drJI.y*Hessian.bz; // yy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+4]+=drJI.z*Hessian.bz; // yz z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+5]+=drJI.z*Hessian.cz; // zz z
            case 2:
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=drJI.y*Hessian.ax;         // xy x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+Dimension]+=drJI.y*Hessian.bx; // yy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=drJI.x*Hessian.ay;           // xx y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=drJI.y*Hessian.ay;         // xy y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+Dimension]+=drJI.y*Hessian.by; // yy y
            case 1:
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=drJI.x*Hessian.ax;     // xx x
              break;
          }

          veci1=DVecX[index1]; vecj1=DVecX[index2];
          veci2=DVecY[index1]; vecj2=DVecY[index2];
          veci3=DVecZ[index1]; vecj3=DVecZ[index2];

          drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
          drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
          drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

          Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
          Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
          Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

          Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
          Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
          Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

          Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
          Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
          Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

          switch(Dimension)
          {
            case 3:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=drJI.z*Hessian.ax; // xz x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+4]-=drJI.z*Hessian.bx; // yz x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+5]-=drJI.z*Hessian.cx; // zz x

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=drJI.z*Hessian.ay; // xz y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+4]-=drJI.z*Hessian.by; // yz y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+5]-=drJI.z*Hessian.cy; // zz y

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=drJI.x*Hessian.az;   // xx z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=drJI.y*Hessian.az; // xy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=drJI.z*Hessian.az; // xz z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=drJI.y*Hessian.bz; // yy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+4]-=drJI.z*Hessian.bz; // yz z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+5]-=drJI.z*Hessian.cz; // zz z
            case 2:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=drJI.y*Hessian.ax;         // xy x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+Dimension]-=drJI.y*Hessian.bx; // yy x
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=drJI.x*Hessian.ay;           // xx y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=drJI.y*Hessian.ay;         // xy y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+Dimension]-=drJI.y*Hessian.by; // yy y
            case 1:
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=drJI.x*Hessian.ax;     // xx x
              break;
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // Eq. 40 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ==============================================
              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
              Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
              Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

              Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
              Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
              Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

              Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
              Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
              Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=drJI.x*Hessian.ax;                           // xx x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=drJI.y*Hessian.bx;                         // yy x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=0.5*(drJI.z*Hessian.bx+drJI.y*Hessian.cx); // yz x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=drJI.z*Hessian.cx;                         // zz x

              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=drJI.x*Hessian.ay;                           // xx y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=drJI.y*Hessian.by;                         // yy y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=0.5*(drJI.z*Hessian.by+drJI.y*Hessian.cy); // yz y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=drJI.z*Hessian.cy;                         // zz y

              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.az;                           // xx z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=drJI.y*Hessian.bz;                         // yy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=0.5*(drJI.z*Hessian.bz+drJI.y*Hessian.cz); // yz z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=drJI.z*Hessian.cz;                         // zz z

              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
              Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
              Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

              Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
              Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
              Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

              Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
              Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
              Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=drJI.x*Hessian.ax;                           // xx x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=drJI.y*Hessian.bx;                         // yy x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=0.5*(drJI.z*Hessian.bx+drJI.y*Hessian.cx); // yz x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=drJI.z*Hessian.cx;                         // zz x

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=drJI.x*Hessian.ay;                           // xx y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=drJI.y*Hessian.by;                         // yy y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=0.5*(drJI.z*Hessian.by+drJI.y*Hessian.cy); // yz y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=drJI.z*Hessian.cy;                         // zz y

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=drJI.x*Hessian.az;                           // xx z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=drJI.y*Hessian.bz;                         // yy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=0.5*(drJI.z*Hessian.bz+drJI.y*Hessian.cz); // yz z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=drJI.z*Hessian.cz;                         // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // Eq. 40 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ==============================================
              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
              Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
              Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

              Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
              Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
              Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

              Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
              Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
              Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=drJI.x*Hessian.ax;                           // xx x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=0.5*(drJI.y*Hessian.ax+drJI.x*Hessian.bx); // xy x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=drJI.y*Hessian.bx;                         // yy x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=drJI.z*Hessian.cx;                         // zz x

              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=drJI.x*Hessian.ay;                           // xx y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=0.5*(drJI.y*Hessian.ay+drJI.x*Hessian.by); // xy y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=drJI.y*Hessian.by;                         // yy y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=drJI.z*Hessian.cy;                         // zz y

              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.az;                           // xx z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=0.5*(drJI.y*Hessian.az+drJI.x*Hessian.bz); // xy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=drJI.y*Hessian.bz;                         // yy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=drJI.z*Hessian.cz;                         // zz z

              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
              Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
              Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

              Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
              Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
              Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

              Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
              Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
              Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=drJI.x*Hessian.ax;                             // xx x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=0.5*(drJI.y*Hessian.ax+drJI.x*Hessian.bx);   // xy x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=drJI.y*Hessian.bx;                           // yy x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=drJI.z*Hessian.cx;                           // zz x

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=drJI.x*Hessian.ay;                           // xx y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=0.5*(drJI.y*Hessian.ay+drJI.x*Hessian.by); // xy y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=drJI.y*Hessian.by;                         // yy y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=drJI.z*Hessian.cy;                         // zz y

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=drJI.x*Hessian.az;                           // xx z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=0.5*(drJI.y*Hessian.az+drJI.x*Hessian.bz); // xy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=drJI.y*Hessian.bz;                         // yy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=drJI.z*Hessian.cz;                         // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
            default:
              // Eq. 40 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ==============================================
              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
              Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
              Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

              Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
              Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
              Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

              Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
              Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
              Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=drJI.x*Hessian.ax;                             // xx x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=0.5*(drJI.z*Hessian.ax+drJI.x*Hessian.cx);   // xz x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=drJI.y*Hessian.bx;                           // yy x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=drJI.z*Hessian.cx;                           // zz x

              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=drJI.x*Hessian.ay;                           // xx y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=0.5*(drJI.z*Hessian.ay+drJI.x*Hessian.cy); // xz y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=drJI.y*Hessian.by;                         // yy y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=drJI.z*Hessian.cy;                         // zz y

              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.az;                           // xx z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=0.5*(drJI.z*Hessian.az+drJI.x*Hessian.cz); // xz z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=drJI.y*Hessian.bz;                         // yy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=drJI.z*Hessian.cz;                         // zz z

              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
              Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
              Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

              Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
              Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
              Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

              Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
              Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
              Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=drJI.x*Hessian.ax;                             // xx x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=0.5*(drJI.z*Hessian.ax+drJI.x*Hessian.cx);   // xz x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=drJI.y*Hessian.bx;                           // yy x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=drJI.z*Hessian.cx;                           // zz x

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=drJI.x*Hessian.ay;                           // xx y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=0.5*(drJI.z*Hessian.ay+drJI.x*Hessian.cy); // xz y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=drJI.y*Hessian.by;                         // yy y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=drJI.z*Hessian.cy;                         // zz y

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=drJI.x*Hessian.az;                           // xx z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=0.5*(drJI.z*Hessian.az+drJI.x*Hessian.cz); // xz z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=drJI.y*Hessian.bz;                         // yy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=drJI.z*Hessian.cz;                         // zz z
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // Eq. 40 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ==============================================
              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
              Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
              Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

              Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
              Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
              Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

              Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
              Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
              Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=drJI.x*Hessian.ax;   // xx x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=drJI.y*Hessian.bx; // yy x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=drJI.z*Hessian.bx; // yz x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=drJI.z*Hessian.cx; // zz x

              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=drJI.x*Hessian.ay;   // xx y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=drJI.y*Hessian.by; // yy y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=drJI.z*Hessian.by; // yz y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=drJI.z*Hessian.cy; // zz y

              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.az;   // xx z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=drJI.y*Hessian.bz; // yy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=drJI.z*Hessian.bz; // yz z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=drJI.z*Hessian.cz; // zz z

              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
              Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
              Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

              Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
              Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
              Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

              Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
              Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
              Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=drJI.x*Hessian.ax;   // xx x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=drJI.y*Hessian.bx; // yy x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=drJI.z*Hessian.bx; // yz x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=drJI.z*Hessian.cx; // zz x

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=drJI.x*Hessian.ay;   // xx y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=drJI.y*Hessian.by; // yy y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=drJI.z*Hessian.by; // yz y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=drJI.z*Hessian.cy; // zz y

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=drJI.x*Hessian.az;   // xx z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=drJI.y*Hessian.bz; // yy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=drJI.z*Hessian.bz; // yz z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=drJI.z*Hessian.cz; // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              // Eq. 40 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ==============================================
              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
              Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
              Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

              Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
              Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
              Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

              Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
              Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
              Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=drJI.x*Hessian.ax;   // xx x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=drJI.z*Hessian.ax; // xz x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=drJI.y*Hessian.bx; // yy x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=drJI.z*Hessian.cx; // zz x

              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=drJI.x*Hessian.ay;   // xx y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=drJI.z*Hessian.ay; // xz y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=drJI.y*Hessian.by; // yy y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=drJI.z*Hessian.cy; // zz y

              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.az;   // xx z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=drJI.z*Hessian.az; // xz z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=drJI.y*Hessian.bz; // yy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=drJI.z*Hessian.cz; // zz z

              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
              Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
              Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

              Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
              Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
              Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

              Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
              Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
              Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=drJI.x*Hessian.ax;   // xx x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=drJI.z*Hessian.ax; // xz x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=drJI.y*Hessian.bx; // yy x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=drJI.z*Hessian.cx; // zz x

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=drJI.x*Hessian.ay;   // xx y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=drJI.z*Hessian.ay; // xz y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=drJI.y*Hessian.by; // yy y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=drJI.z*Hessian.cy; // zz y

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=drJI.x*Hessian.az;   // xx z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=drJI.z*Hessian.az; // xz z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=drJI.y*Hessian.bz; // yy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=drJI.z*Hessian.cz; // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // Eq. 40 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ==============================================
              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
              Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
              Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

              Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
              Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
              Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

              Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
              Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
              Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=drJI.x*Hessian.ax;     // xx x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=drJI.y*Hessian.ax;   // xy x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=drJI.y*Hessian.bx;   // yy x
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=drJI.z*Hessian.cx;   // zz x

              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=drJI.x*Hessian.ay;   // xx y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=drJI.y*Hessian.ay; // xy y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=drJI.y*Hessian.by; // yy y
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=drJI.z*Hessian.cy; // zz y

              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=drJI.x*Hessian.az;   // xx z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=drJI.y*Hessian.az; // xy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=drJI.y*Hessian.bz; // yy z
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=drJI.z*Hessian.cz; // zz z

              veci1=DVecX[index1]; vecj1=DVecX[index2];
              veci2=DVecY[index1]; vecj2=DVecY[index2];
              veci3=DVecZ[index1]; vecj3=DVecZ[index2];

              drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
              drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
              drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

              Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
              Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
              Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

              Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
              Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
              Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

              Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
              Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
              Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n]-=drJI.x*Hessian.ax;   // xx x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+1]-=drJI.y*Hessian.ax; // xy x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+2]-=drJI.y*Hessian.bx; // yy x
              if(index_j2.x>=0) HessianMatrix.element[index_j2.x][n+3]-=drJI.z*Hessian.cx; // zz x

              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n]-=drJI.x*Hessian.ay;   // xx y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+1]-=drJI.y*Hessian.ay; // xy y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+2]-=drJI.y*Hessian.by; // yy y
              if(index_j2.y>=0) HessianMatrix.element[index_j2.y][n+3]-=drJI.z*Hessian.cy; // zz y

              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n]-=drJI.x*Hessian.az;   // xx z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+1]-=drJI.y*Hessian.az; // xy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+2]-=drJI.y*Hessian.bz; // yy z
              if(index_j2.z>=0) HessianMatrix.element[index_j2.z][n+3]-=drJI.z*Hessian.cz; // zz z
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
void HessianAtomicStrainStrain(REAL_MATRIX HessianMatrix,
                       REAL f1,REAL f2,VECTOR dr,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB)
{
  int n;
  VECTOR dJI,drJI;
  REAL temp;

  n=NumberOfCoordinatesMinimizationVariables;

  dJI.x=(posB.x-comB.x)-(posA.x-comA.x);
  dJI.y=(posB.y-comB.y)-(posA.y-comA.y);
  dJI.z=(posB.z-comB.z)-(posA.z-comA.z);

  drJI.x=dr.x+dJI.x;
  drJI.y=dr.y+dJI.y;
  drJI.z=dr.z+dJI.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      temp=0.0;
      // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
      // ============================================================================
      switch(Dimension)
      {
        case 3:
          temp+=f2*(drJI.x*dr.x*dr.z*dr.z+dr.x*drJI.x*dr.z*dr.z);
          temp+=f2*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z);
          temp+=f2*drJI.z*dr.z*dr.z*dr.z+2.0*f1*(drJI.z*dr.z);
        case 2:
          temp+=f2*(drJI.x*dr.x*dr.y*dr.y+dr.x*drJI.x*dr.y*dr.y);
          temp+=f2*drJI.y*dr.y*dr.y*dr.y+2.0*f1*(drJI.y*dr.y);
        case 1:
          temp+=f2*drJI.x*dr.x*dr.x*dr.x+2.0*f1*drJI.x*dr.x;
          break;
      }

      // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
      // zero for non-rigid particles
      // ============================================================================
      switch(Dimension)
      {
        case 3:
          temp+=0.5*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);
          temp+=0.5*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.z*dr.z+dJI.z*dr.z);
          temp+=0.25*f2*(drJI.z*dr.z+drJI.z*dr.z)*(dJI.z*dr.z+dJI.z*dr.z)+f1*(drJI.z*dJI.z);
        case 2:
          temp+=0.5*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
          temp+=0.25*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.y*dr.y+dJI.y*dr.y)+f1*(drJI.y*dJI.y);
        case 1:
          temp+=f2*drJI.x*dr.x*dJI.x*dr.x+f1*drJI.x*dJI.x;
          break;
      }
      HessianMatrix.element[n][n]+=temp;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          temp=0.0;
          // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ============================================================================
          switch(Dimension)
          {
            case 3:
              temp+=f2*(drJI.x*dr.x*dr.z*dr.z+dr.x*drJI.x*dr.z*dr.z);
              temp+=f2*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z);
              temp+=f2*drJI.z*dr.z*dr.z*dr.z+2.0*f1*(drJI.z*dr.z);
            case 2:
              temp+=f2*(drJI.x*dr.x*dr.y*dr.y+dr.x*drJI.x*dr.y*dr.y);
              temp+=f2*drJI.y*dr.y*dr.y*dr.y+2.0*f1*(drJI.y*dr.y);
            case 1:
              temp+=f2*drJI.x*dr.x*dr.x*dr.x+2.0*f1*drJI.x*dr.x;
              break;
          }

          // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // zero for non-rigid particles
          // ============================================================================
          switch(Dimension)
          {
            case 3:
              temp+=0.5*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);
              temp+=0.5*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.z*dr.z+dJI.z*dr.z);
              temp+=0.25*f2*(drJI.z*dr.z+drJI.z*dr.z)*(dJI.z*dr.z+dJI.z*dr.z)+f1*(drJI.z*dJI.z);
              HessianMatrix.element[n][n+2]+=temp/9.0;
              HessianMatrix.element[n+1][n+2]+=temp/9.0;
              HessianMatrix.element[n+2][n+2]+=temp/9.0;
            case 2:
              temp+=0.5*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
              temp+=0.25*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.y*dr.y+dJI.y*dr.y)+f1*(drJI.y*dJI.y);
              HessianMatrix.element[n][n+1]+=temp/9.0;
              HessianMatrix.element[n+1][n+1]+=temp/9.0;
            case 1:
              temp+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.x*dr.x+dJI.x*dr.x)+0.25*f1*(drJI.x*dJI.x+drJI.x*dJI.x+dJI.x*drJI.x+dJI.x*drJI.x);
              HessianMatrix.element[n][n]+=temp/9.0;
              break;
          }
          break;
        case ANISOTROPIC:
          // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ============================================================================
          switch(Dimension)
          {
            case 3:
              HessianMatrix.element[n][n+2]+=f2*drJI.x*dr.x*dr.z*dr.z;                      // xxzz
              HessianMatrix.element[n+1][n+2]+=f2*drJI.y*dr.y*dr.z*dr.z;                    // yyzz
              HessianMatrix.element[n+2][n+2]+=f2*drJI.z*dr.z*dr.z*dr.z+2.0*f1*drJI.z*dr.z; // zzzz
            case 2:
              HessianMatrix.element[n][n+1]+=f2*drJI.x*dr.x*dr.y*dr.y;                      // xxyy
              HessianMatrix.element[n+1][n+1]+=f2*drJI.y*dr.y*dr.y*dr.y+2.0*f1*drJI.y*dr.y; // yyyy
            case 1:
              HessianMatrix.element[n][n]+=f2*drJI.x*dr.x*dr.x*dr.x+2.0*f1*drJI.x*dr.x;     // xxxx
              break;
          }

          // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // zero for non-rigid particles
          // ============================================================================
          switch(Dimension)
          {
            case 3:
              HessianMatrix.element[n][n+2]+=f2*(drJI.x*dr.x*dJI.z*dr.z);                      // xxzz
              HessianMatrix.element[n+1][n+2]+=f2*(drJI.y*dr.y*dJI.z*dr.z);                    // yyzz
              HessianMatrix.element[n+2][n+2]+=f2*(drJI.z*dr.z*dJI.z*dr.z)+f1*(drJI.z*dJI.z);  // zzzz
            case 2:
              HessianMatrix.element[n][n+1]+=f2*(drJI.x*dr.x*dJI.y*dr.y);                      // xxyy
              HessianMatrix.element[n+1][n+1]+=f2*(drJI.y*dr.y*dJI.y*dr.y)+f1*(drJI.y*dJI.y);  // yyyy
            case 1:
              HessianMatrix.element[n][n]+=f2*(drJI.x*dr.x*dJI.x*dr.x)+f1*(drJI.x*dJI.x);      // xxxx
              break;
          }
          break;
        case REGULAR:
          // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ============================================================================
          switch(Dimension)
          {
            case 3:
              HessianMatrix.element[n][n+2]+=f2*drJI.x*dr.x*dr.x*dr.z+f1*(drJI.x*dr.z);                                               // xxxz
              HessianMatrix.element[n][n+4]+=f2*drJI.x*dr.x*dr.y*dr.z;                                                                // xxyz
              HessianMatrix.element[n][n+5]+=f2*drJI.x*dr.x*dr.z*dr.z;                                                                // xxzz
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(drJI.x*dr.y*dr.x*dr.z+dr.x*drJI.y*dr.x*dr.z)+0.5*f1*(drJI.y*dr.z);             // xyxz
              HessianMatrix.element[n+1][n+4]+=0.5*f2*(drJI.x*dr.y*dr.y*dr.z+dr.x*drJI.y*dr.y*dr.z)+0.5*f1*(drJI.x*dr.z);             // xyyz
              HessianMatrix.element[n+1][n+5]+=0.5*f2*(drJI.x*dr.y*dr.z*dr.z+dr.x*drJI.y*dr.z*dr.z);                                  // xyzz
              HessianMatrix.element[n+2][n+2]+=0.5*f2*(drJI.x*dr.z*dr.x*dr.z+dr.x*drJI.z*dr.x*dr.z)+0.5*f1*(drJI.z*dr.z+drJI.x*dr.x); // xzxz
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(drJI.x*dr.z*dr.y*dr.y+dr.x*drJI.z*dr.y*dr.y);                                  // xzyy
              HessianMatrix.element[n+2][n+4]+=0.5*f2*(drJI.x*dr.z*dr.y*dr.z+dr.x*drJI.z*dr.y*dr.z)+0.5*f1*(drJI.x*dr.y);             // xzyz
              HessianMatrix.element[n+2][n+5]+=0.5*f2*(drJI.x*dr.z*dr.z*dr.z+dr.x*drJI.z*dr.z*dr.z)+f1*(drJI.x*dr.z);                 // xzzz
              HessianMatrix.element[n+3][n+4]+=0.5*f2*(drJI.y*dr.y*dr.y*dr.z+dr.y*drJI.y*dr.y*dr.z)+f1*(drJI.y*dr.z);                 // yyyz
              HessianMatrix.element[n+3][n+5]+=0.5*f2*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z);                                  // yyzz
              HessianMatrix.element[n+4][n+4]+=0.5*f2*(drJI.y*dr.z*dr.y*dr.z+dr.y*drJI.z*dr.y*dr.z)+0.5*f1*(drJI.z*dr.z+drJI.y*dr.y); // yzyz
              HessianMatrix.element[n+4][n+5]+=0.5*f2*(drJI.y*dr.z*dr.z*dr.z+dr.y*drJI.z*dr.z*dr.z)+f1*(drJI.y*dr.z);                 // yzzz
              HessianMatrix.element[n+5][n+5]+=f2*drJI.z*dr.z*dr.z*dr.z+2.0*f1*drJI.z*dr.z;                                           // zzzz
            case 2:
              HessianMatrix.element[n][n+1]+=f2*drJI.x*dr.x*dr.x*dr.y+f1*drJI.x*dr.y;                                                 // xxxy
              HessianMatrix.element[n][n+Dimension]+=f2*drJI.x*dr.x*dr.y*dr.y;                                                        // xxyy
              HessianMatrix.element[n+1][n+1]+=0.5*f2*(drJI.x*dr.y*dr.x*dr.y+dr.x*drJI.y*dr.x*dr.y)+0.5*f1*(drJI.y*dr.y+drJI.x*dr.x); // xyxy
              HessianMatrix.element[n+1][n+Dimension]+=0.5*f2*(drJI.x*dr.y*dr.y*dr.y+dr.x*drJI.y*dr.y*dr.y)+f1*(drJI.x*dr.y);         // xyyy
              HessianMatrix.element[n+Dimension][n+Dimension]+=f2*drJI.y*dr.y*dr.y*dr.y+2.0*f1*drJI.y*dr.y;                           // yyyy
            case 1:
              HessianMatrix.element[n][n]+=f2*drJI.x*dr.x*dr.x*dr.x+2.0*f1*drJI.x*dr.x; // xxxx
              break;
          }

          // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ============================================================================
          switch(Dimension)
          {
            case 3:
              HessianMatrix.element[n][n+2]+=0.5*f2*drJI.x*dr.x*(dJI.x*dr.z+dJI.z*dr.x)+0.5*f1*(drJI.x*dJI.z);                                // xxxz
              HessianMatrix.element[n][n+4]+=0.5*f2*drJI.x*dr.x*(dJI.y*dr.z+dJI.z*dr.y);                                                      // xxyz
              HessianMatrix.element[n][n+5]+=0.5*f2*drJI.x*dr.x*(dJI.z*dr.z+dJI.z*dr.z);                                                      // xxzz
              HessianMatrix.element[n+1][n+2]+=0.25*f2*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*f1*(drJI.y*dJI.z);              // xyxz
              HessianMatrix.element[n+1][n+4]+=0.25*f2*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*f1*(drJI.x*dJI.z);              // xyyz
              HessianMatrix.element[n+1][n+5]+=0.25*f2*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);                                     // xyzz
              HessianMatrix.element[n+2][n+2]+=0.25*f2*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*f1*(drJI.z*dJI.z+drJI.x*dJI.x); // xzxz
              HessianMatrix.element[n+2][n+3]+=0.25*f2*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);                                     // xzyy
              HessianMatrix.element[n+2][n+4]+=0.25*f2*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*f1*(drJI.x*dJI.y);              // xzyz
              HessianMatrix.element[n+2][n+5]+=0.25*f2*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*f1*(drJI.x*dJI.z);               // xzzz
              HessianMatrix.element[n+3][n+4]+=0.5*f2*drJI.y*dr.y*(dJI.y*dr.z+dJI.z*dr.y)+0.5*f1*(drJI.y*dJI.z);                              // yyyz
              HessianMatrix.element[n+3][n+5]+=0.5*f2*drJI.y*dr.y*(dJI.z*dr.z+dJI.z*dr.z);                                                    // yyzz
              HessianMatrix.element[n+4][n+4]+=0.25*f2*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*f1*(drJI.z*dJI.z+drJI.y*dJI.y); // yzyz
              HessianMatrix.element[n+4][n+5]+=0.25*f2*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*f1*(drJI.y*dJI.z);               // yzzz
              HessianMatrix.element[n+5][n+5]+=f2*drJI.z*dr.z*dJI.z*dr.z+f1*drJI.z*dJI.z;                                                     // zzzz
            case 2:
              HessianMatrix.element[n][n+1]+=0.5*f2*drJI.x*dr.x*(dJI.x*dr.y+dJI.y*dr.x)+0.5*f1*(drJI.x*dJI.y);                                // xxxy
              HessianMatrix.element[n][n+Dimension]+=0.5*f2*drJI.x*dr.x*(dJI.y*dr.y+dJI.y*dr.y);                                              // xxyy
              HessianMatrix.element[n+1][n+1]+=0.25*f2*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.y+dJI.y*dr.x)+0.25*f1*(drJI.y*dJI.y+drJI.x*dJI.x); // xyxy
              HessianMatrix.element[n+1][n+Dimension]+=0.25*f2*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.y+dJI.y*dr.y)+0.5*f1*(drJI.x*dJI.y);       // xyyy
              HessianMatrix.element[n+Dimension][n+Dimension]+=f2*drJI.y*dr.y*dJI.y*dr.y+f1*drJI.y*dJI.y;                                     // yyyy
            case 1:
              HessianMatrix.element[n][n]+=f2*drJI.x*dr.x*dJI.x*dr.x+f1*drJI.x*dJI.x;                                                         // xxxx
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ============================================================================
          switch(Dimension)
          {
            case 3:
              HessianMatrix.element[n][n+2]+=f2*drJI.x*dr.x*dr.x*dr.z+f1*drJI.x*dr.z;       // xxxz
              HessianMatrix.element[n][n+4]+=f2*drJI.x*dr.x*dr.y*dr.z;                      // xxyz
              HessianMatrix.element[n][n+5]+=f2*drJI.x*dr.x*dr.z*dr.z;                      // xxzz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*drJI.y*dr.x*dr.z+f1*drJI.y*dr.z;     // xyxz
              HessianMatrix.element[n+1][n+4]+=f2*dr.x*drJI.y*dr.y*dr.z;                    // xyyz
              HessianMatrix.element[n+1][n+5]+=f2*dr.x*drJI.y*dr.z*dr.z;                    // xyzz
              HessianMatrix.element[n+2][n+2]+=f2*dr.x*drJI.z*dr.x*dr.z+f1*drJI.z*dr.z;     // xzxz
              HessianMatrix.element[n+2][n+3]+=f2*dr.x*drJI.z*dr.y*dr.y;                    // xzyy
              HessianMatrix.element[n+2][n+4]+=f2*dr.x*drJI.z*dr.y*dr.z;                    // xzyz
              HessianMatrix.element[n+2][n+5]+=f2*dr.x*drJI.z*dr.z*dr.z;                    // xzzz
              HessianMatrix.element[n+3][n+4]+=f2*drJI.y*dr.y*dr.y*dr.z+f1*drJI.y*dr.z;     // yyyz
              HessianMatrix.element[n+3][n+5]+=f2*drJI.y*dr.y*dr.z*dr.z;                    // yyzz
              HessianMatrix.element[n+4][n+4]+=f2*dr.y*drJI.z*dr.y*dr.z+f1*drJI.z*dr.z;     // yzyz
              HessianMatrix.element[n+4][n+5]+=f2*dr.y*drJI.z*dr.z*dr.z;                    // yzzz
              HessianMatrix.element[n+5][n+5]+=f2*drJI.z*dr.z*dr.z*dr.z+2.0*f1*drJI.z*dr.z; // zzzz
            case 2:
              HessianMatrix.element[n][n+1]+=f2*drJI.x*dr.x*dr.x*dr.y+f1*drJI.x*dr.y;                       // xxxy
              HessianMatrix.element[n][n+Dimension]+=f2*drJI.x*dr.x*dr.y*dr.y;                              // xxyy
              HessianMatrix.element[n+1][n+1]+=f2*dr.x*drJI.y*dr.x*dr.y+f1*drJI.y*dr.y;                     // xyxy
              HessianMatrix.element[n+1][n+Dimension]+=f2*dr.x*drJI.y*dr.y*dr.y;                            // xyyy
              HessianMatrix.element[n+Dimension][n+Dimension]+=f2*drJI.y*dr.y*dr.y*dr.y+2.0*f1*drJI.y*dr.y; // yyyy
            case 1:
              HessianMatrix.element[n][n  ]+=f2*drJI.x*dr.x*dr.x*dr.x+2.0*f1*drJI.x*dr.x;   // xxxx
              break;
          }

          // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ============================================================================
          switch(Dimension)
          {
            case 3:
              HessianMatrix.element[n][n+2]+=f2*drJI.x*dr.x*dJI.z*dr.x+f1*drJI.x*dJI.z;    // xxxz
              HessianMatrix.element[n][n+4]+=f2*drJI.x*dr.x*dJI.z*dr.y;                    // xxyz
              HessianMatrix.element[n][n+5]+=f2*drJI.x*dr.x*dJI.z*dr.z;                    // xxzz
              HessianMatrix.element[n+1][n+2]+=f2*drJI.y*dr.x*dJI.z*dr.x+f1*drJI.y*dJI.z;  // xyxz
              HessianMatrix.element[n+1][n+4]+=f2*drJI.y*dr.x*dJI.z*dr.y;                  // xyyz
              HessianMatrix.element[n+1][n+5]+=f2*drJI.y*dr.x*dJI.z*dr.z;                  // xyzz
              HessianMatrix.element[n+2][n+2]+=f2*drJI.z*dr.x*dJI.z*dr.x+f1*drJI.z*dJI.z;  // xzxz
              HessianMatrix.element[n+2][n+3]+=f2*drJI.z*dr.x*dJI.y*dr.y;                  // xzyy
              HessianMatrix.element[n+2][n+4]+=f2*drJI.z*dr.x*dJI.z*dr.y;                  // xzyz
              HessianMatrix.element[n+2][n+5]+=f2*drJI.z*dr.x*dJI.z*dr.z;                  // xzzz
              HessianMatrix.element[n+3][n+4]+=f2*drJI.y*dr.y*dJI.z*dr.y+f1*drJI.y*dJI.z;  // yyyz
              HessianMatrix.element[n+3][n+5]+=f2*drJI.y*dr.y*dJI.z*dr.z;                  // yyzz
              HessianMatrix.element[n+4][n+4]+=f2*drJI.z*dr.y*dJI.z*dr.y+f1*drJI.z*dJI.z;  // yzyz
              HessianMatrix.element[n+4][n+5]+=f2*drJI.z*dr.y*dJI.z*dr.z;                  // yzzz
              HessianMatrix.element[n+5][n+5]+=f2*drJI.z*dr.z*dJI.z*dr.z+f1*drJI.z*dJI.z;  // zzzz
            case 2:
              HessianMatrix.element[n][n+1]+=f2*drJI.x*dr.x*dJI.y*dr.x+f1*drJI.x*dJI.y;                    // xxxy
              HessianMatrix.element[n][n+Dimension]+=f2*drJI.x*dr.x*dJI.y*dr.y;                            // xxyy
              HessianMatrix.element[n+1][n+1]+=f2*drJI.y*dr.x*dJI.y*dr.x+f1*drJI.y*dJI.y;                  // xyxy
              HessianMatrix.element[n+1][n+Dimension]+=f2*drJI.y*dr.x*dJI.y*dr.y;                          // xyyy
              HessianMatrix.element[n+Dimension][n+Dimension]+=f2*drJI.y*dr.y*dJI.y*dr.y+f1*drJI.y*dJI.y;  // yyyy
            case 1:
              HessianMatrix.element[n][n  ]+=f2*drJI.x*dr.x*dJI.x*dr.x+f1*drJI.x*dJI.x;                    // xxxx
              break;
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              // xxxx
              HessianMatrix.element[n][n]+=f2*drJI.x*dr.x*dr.x*dr.x+
                                           +2.0*f1*drJI.x*dr.x;
              // xxyy
              HessianMatrix.element[n][n+1]+=0.5*f2*(drJI.x*dr.x*dr.y*dr.y+dr.x*drJI.x*dr.y*dr.y);
              // xxyz
              HessianMatrix.element[n][n+2]+=0.5*f2*(drJI.x*dr.x*dr.y*dr.z+dr.x*drJI.x*dr.y*dr.z);
              // xxzz
              HessianMatrix.element[n][n+3]+=0.5*f2*(drJI.x*dr.x*dr.z*dr.z+dr.x*drJI.x*dr.z*dr.z);

              // yyyy
              HessianMatrix.element[n+1][n+1]+=0.5*f2*(drJI.y*dr.y*dr.y*dr.y+dr.y*drJI.y*dr.y*dr.y)
                                               +2.0*f1*(drJI.y*dr.y);
              // yyyz
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(drJI.y*dr.y*dr.y*dr.z+dr.y*drJI.y*dr.y*dr.z)
                                               +f1*(drJI.y*dr.z);
              // yyzz
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z);

              // yzyz
              HessianMatrix.element[n+2][n+2]+=0.5*f2*(drJI.y*dr.z*dr.y*dr.z+dr.y*drJI.z*dr.y*dr.z)
                                               +0.5*f1*(drJI.z*dr.z+drJI.y*dr.y);
              // yzzz
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(drJI.y*dr.z*dr.z*dr.z+dr.y*drJI.z*dr.z*dr.z)
                                               +f1*(drJI.y*dr.z);

              // zzzz
              HessianMatrix.element[n+3][n+3]+=0.5*f2*(drJI.z*dr.z*dr.z*dr.z+dr.z*drJI.z*dr.z*dr.z)
                                               +2.0*f1*(drJI.z*dr.z);


              // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              // xxxx
              HessianMatrix.element[n][n]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.x*dr.x+dJI.x*dr.x)
                                           +0.25*f1*(drJI.x*dJI.x+drJI.x*dJI.x+dJI.x*drJI.x+dJI.x*drJI.x);
              // xxyy
              HessianMatrix.element[n][n+1]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
              // xxyz
              HessianMatrix.element[n][n+2]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.y*dr.z+dJI.z*dr.y);
              // xxzz
              HessianMatrix.element[n][n+3]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);

              // yyyy
              HessianMatrix.element[n+1][n+1]+=0.25*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.y*dr.y+dJI.y*dr.y)
                                               +f1*(drJI.y*dJI.y);
              // yyyz
              HessianMatrix.element[n+1][n+2]+=0.25*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.y*dr.z+dJI.z*dr.y)
                                               +0.5*f1*(drJI.y*dJI.z);
              // yyzz
              HessianMatrix.element[n+1][n+3]+=0.25*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.z*dr.z+dJI.z*dr.z);

              // yzyz
              HessianMatrix.element[n+2][n+2]+=0.25*f2*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.y*dr.z+dJI.z*dr.y)
                                               +0.25*f1*(drJI.z*dJI.z+drJI.y*dJI.y);
              // yzzz
              HessianMatrix.element[n+2][n+3]+=0.25*f2*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.z*dr.z+dJI.z*dr.z)
                                               +0.5*f1*(drJI.y*dJI.z);

              // zzzz
              HessianMatrix.element[n+3][n+3]+=0.25*f2*(drJI.z*dr.z+drJI.z*dr.z)*(dJI.z*dr.z+dJI.z*dr.z)
                                               +f1*(drJI.z*dJI.z);
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              // xxxx
              HessianMatrix.element[n][n]+=f2*drJI.x*dr.x*dr.x*dr.x+
                                           +2.0*f1*drJI.x*dr.x;
              // xxxy
              HessianMatrix.element[n][n+1]+=0.5*f2*(drJI.x*dr.x*dr.x*dr.y+dr.x*drJI.x*dr.x*dr.y)
                                            +0.5*f1*(drJI.x*dr.y+drJI.x*dr.y);
              // xxyy
              HessianMatrix.element[n][n+2]+=0.5*f2*(drJI.x*dr.x*dr.y*dr.y+dr.x*drJI.x*dr.y*dr.y);
              // xxzz
              HessianMatrix.element[n][n+3]+=0.5*f2*(drJI.x*dr.x*dr.z*dr.z+dr.x*drJI.x*dr.z*dr.z);

              // xyxy
              HessianMatrix.element[n+1][n+1]+=0.5*f2*(drJI.x*dr.y*dr.x*dr.y+dr.x*drJI.y*dr.x*dr.y)
                                               +0.5*f1*(drJI.y*dr.y+drJI.x*dr.x);
              // xyyy
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(drJI.x*dr.y*dr.y*dr.y+dr.x*drJI.y*dr.y*dr.y)
                                               +f1*(drJI.x*dr.y);
              // xyzz
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(drJI.x*dr.y*dr.z*dr.z+dr.x*drJI.y*dr.z*dr.z);

              // yyyy
              HessianMatrix.element[n+2][n+2]+=0.5*f2*(drJI.y*dr.y*dr.y*dr.y+dr.y*drJI.y*dr.y*dr.y)
                                               +2.0*f1*(drJI.y*dr.y);
              // yyzz
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z);

              // zzzz
              HessianMatrix.element[n+3][n+3]+=0.5*f2*(drJI.z*dr.z*dr.z*dr.z+dr.z*drJI.z*dr.z*dr.z)
                                               +2.0*f1*(drJI.z*dr.z);


              // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              // xxxx
              HessianMatrix.element[n][n]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.x*dr.x+dJI.x*dr.x)
                                           +0.25*f1*(drJI.x*dJI.x+drJI.x*dJI.x+dJI.x*drJI.x+dJI.x*drJI.x);
              // xxxy
              HessianMatrix.element[n][n+1]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.x*dr.y+dJI.y*dr.x)
                                             +0.25*f1*(drJI.x*dJI.y+drJI.x*dJI.y);
              // xxyy
              HessianMatrix.element[n][n+2]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
              // xxzz
              HessianMatrix.element[n][n+3]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);

              // xyxy
              HessianMatrix.element[n+1][n+1]+=0.25*f2*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.y+dJI.y*dr.x)
                                               +0.25*f1*(drJI.y*dJI.y+drJI.x*dJI.x);
              // xyyy
              HessianMatrix.element[n+1][n+2]+=0.25*f2*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.y+dJI.y*dr.y)
                                               +0.5*f1*(drJI.x*dJI.y);
              // xyzz
              HessianMatrix.element[n+1][n+3]+=0.25*f2*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);

              // yyyy
              HessianMatrix.element[n+2][n+2]+=0.25*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.y*dr.y+dJI.y*dr.y)
                                               +f1*(drJI.y*dJI.y);
              // yyzz
              HessianMatrix.element[n+2][n+3]+=0.25*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.z*dr.z+dJI.z*dr.z);

              // zzzz
              HessianMatrix.element[n+3][n+3]+=0.25*f2*(drJI.z*dr.z+drJI.z*dr.z)*(dJI.z*dr.z+dJI.z*dr.z)
                                               +f1*(drJI.z*dJI.z);
              break;
            case MONOCLINIC_BETA_ANGLE:
            default:
              // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              // xxxx
              HessianMatrix.element[n][n]+=f2*drJI.x*dr.x*dr.x*dr.x+
                                           +2.0*f1*drJI.x*dr.x;
              // xxxz
              HessianMatrix.element[n][n+1]+=0.5*f2*(drJI.x*dr.x*dr.x*dr.z+dr.x*drJI.x*dr.x*dr.z)
                                            +0.5*f1*(drJI.x*dr.z+drJI.x*dr.z);
              // xxyy
              HessianMatrix.element[n][n+2]+=0.5*f2*(drJI.x*dr.x*dr.y*dr.y+dr.x*drJI.x*dr.y*dr.y);
              // xxzz
              HessianMatrix.element[n][n+3]+=0.5*f2*(drJI.x*dr.x*dr.z*dr.z+dr.x*drJI.x*dr.z*dr.z);

              // xzxz
              HessianMatrix.element[n+1][n+1]+=0.5*f2*(drJI.x*dr.z*dr.x*dr.z+dr.x*drJI.z*dr.x*dr.z)
                                               +0.5*f1*(drJI.z*dr.z+drJI.x*dr.x);
              // xzyy
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(drJI.x*dr.z*dr.y*dr.y+dr.x*drJI.z*dr.y*dr.y);
              // xzzz
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(drJI.x*dr.z*dr.z*dr.z+dr.x*drJI.z*dr.z*dr.z)
                                               +f1*(drJI.x*dr.z);

              // yyyy
              HessianMatrix.element[n+2][n+2]+=0.5*f2*(drJI.y*dr.y*dr.y*dr.y+dr.y*drJI.y*dr.y*dr.y)
                                               +2.0*f1*(drJI.y*dr.y);
              // yyzz
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z);

              // zzzz
              HessianMatrix.element[n+3][n+3]+=0.5*f2*(drJI.z*dr.z*dr.z*dr.z+dr.z*drJI.z*dr.z*dr.z)
                                               +2.0*f1*(drJI.z*dr.z);


              // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              // xxxx
              HessianMatrix.element[n][n]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.x*dr.x+dJI.x*dr.x)
                                           +0.25*f1*(drJI.x*dJI.x+drJI.x*dJI.x+dJI.x*drJI.x+dJI.x*drJI.x);
              // xxxz
              HessianMatrix.element[n][n+1]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)
                                             +0.25*f1*(drJI.x*dJI.z+drJI.x*dJI.z);
              // xxyy
              HessianMatrix.element[n][n+2]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
              // xxzz
              HessianMatrix.element[n][n+3]+=0.25*f2*(drJI.x*dr.x+drJI.x*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);

              // xzxz
              HessianMatrix.element[n+1][n+1]+=0.25*f2*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)
                                               +0.25*f1*(drJI.z*dJI.z+drJI.x*dJI.x);
              // xzyy
              HessianMatrix.element[n+1][n+2]+=0.25*f2*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
              // xzzz
              HessianMatrix.element[n+1][n+3]+=0.25*f2*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.z*dr.z+dJI.z*dr.z)
                                               +0.5*f1*(drJI.x*dJI.z);

              // yyyy
              HessianMatrix.element[n+2][n+2]+=0.25*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.y*dr.y+dJI.y*dr.y)
                                               +f1*(drJI.y*dJI.y);
              // yyzz
              HessianMatrix.element[n+2][n+3]+=0.25*f2*(drJI.y*dr.y+drJI.y*dr.y)*(dJI.z*dr.z+dJI.z*dr.z);

              // zzzz
              HessianMatrix.element[n+3][n+3]+=0.25*f2*(drJI.z*dr.z+drJI.z*dr.z)*(dJI.z*dr.z+dJI.z*dr.z)
                                               +f1*(drJI.z*dJI.z);

              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              HessianMatrix.element[n][n  ]+=f2*drJI.x*dr.x*dr.x*dr.x+2.0*f1*drJI.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*drJI.x*dr.x*dr.y*dr.y;                      // xxyy
              HessianMatrix.element[n][n+2]+=f2*drJI.x*dr.x*dr.y*dr.z;                      // xxyz
              HessianMatrix.element[n][n+3]+=f2*drJI.x*dr.x*dr.z*dr.z;                      // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*drJI.y*dr.y*dr.y*dr.y+2.0*f1*drJI.y*dr.y; // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*drJI.y*dr.y*dr.y*dr.z+f1*drJI.y*dr.z;     // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*drJI.y*dr.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*drJI.z*dr.y*dr.z+f1*drJI.z*dr.z;     // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*drJI.z*dr.z*dr.z;                    // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*drJI.z*dr.z*dr.z*dr.z+2.0*f1*drJI.z*dr.z; // zzzz


              // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              HessianMatrix.element[n][n  ]+=f2*drJI.x*dr.x*dJI.x*dr.x+f1*drJI.x*dJI.x;    // xxxx
              HessianMatrix.element[n][n+1]+=f2*drJI.x*dr.x*dJI.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+2]+=f2*drJI.x*dr.x*dJI.z*dr.y;                    // xxyz
              HessianMatrix.element[n][n+3]+=f2*drJI.x*dr.x*dJI.z*dr.z;                    // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*drJI.y*dr.y*dJI.y*dr.y+f1*drJI.y*dJI.y;  // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*drJI.y*dr.y*dJI.z*dr.y+f1*drJI.y*dJI.z;  // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*drJI.y*dr.y*dJI.z*dr.z;                  // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*drJI.z*dr.y*dJI.z*dr.y+f1*drJI.z*dJI.z;  // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*drJI.z*dr.y*dJI.z*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*drJI.z*dr.z*dJI.z*dr.z+f1*drJI.z*dJI.z;  // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              HessianMatrix.element[n][n  ]+=f2*drJI.x*dr.x*dr.x*dr.x+2.0*f1*drJI.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*drJI.x*dr.x*dr.x*dr.z+f1*drJI.x*dr.z;       // xxxz
              HessianMatrix.element[n][n+2]+=f2*drJI.x*dr.x*dr.y*dr.y;                      // xxyy
              HessianMatrix.element[n][n+3]+=f2*drJI.x*dr.x*dr.z*dr.z;                      // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*drJI.z*dr.x*dr.z+f1*drJI.z*dr.z;     // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*drJI.z*dr.y*dr.y;                    // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*drJI.z*dr.z*dr.z;                    // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*drJI.y*dr.y*dr.y*dr.y+2.0*f1*drJI.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*drJI.y*dr.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*drJI.z*dr.z*dr.z*dr.z+2.0*f1*drJI.z*dr.z; // zzzz


              // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              HessianMatrix.element[n][n  ]+=f2*drJI.x*dr.x*dJI.x*dr.x+f1*drJI.x*dJI.x;    // xxxx
              HessianMatrix.element[n][n+1]+=f2*drJI.x*dr.x*dJI.z*dr.x+f1*drJI.x*dJI.z;    // xxxz
              HessianMatrix.element[n][n+2]+=f2*drJI.x*dr.x*dJI.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+3]+=f2*drJI.x*dr.x*dJI.z*dr.z;                    // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*drJI.z*dr.x*dJI.z*dr.x+f1*drJI.z*dJI.z;  // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*drJI.z*dr.x*dJI.y*dr.y;                  // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*drJI.z*dr.x*dJI.z*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*drJI.y*dr.y*dJI.y*dr.y+f1*drJI.y*dJI.y;  // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*drJI.y*dr.y*dJI.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*drJI.z*dr.z*dJI.z*dr.z+f1*drJI.z*dJI.z;  // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first and second terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              HessianMatrix.element[n][n  ]+=f2*drJI.x*dr.x*dr.x*dr.x+2.0*f1*drJI.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*drJI.x*dr.x*dr.x*dr.y+f1*drJI.x*dr.y;       // xxxy
              HessianMatrix.element[n][n+2]+=f2*drJI.x*dr.x*dr.y*dr.y;                      // xxyy
              HessianMatrix.element[n][n+3]+=f2*drJI.x*dr.x*dr.z*dr.z;                      // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*drJI.y*dr.x*dr.y+f1*drJI.y*dr.y;     // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*drJI.y*dr.y*dr.y;                    // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*drJI.y*dr.z*dr.z;                    // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*drJI.y*dr.y*dr.y*dr.y+2.0*f1*drJI.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*drJI.y*dr.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*drJI.z*dr.z*dr.z*dr.z+2.0*f1*drJI.z*dr.z; // zzzz


              // the third and last terms of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ============================================================================
              HessianMatrix.element[n][n  ]+=f2*drJI.x*dr.x*dJI.x*dr.x+f1*drJI.x*dJI.x;    // xxxx
              HessianMatrix.element[n][n+1]+=f2*drJI.x*dr.x*dJI.y*dr.x+f1*drJI.x*dJI.y;    // xxxy
              HessianMatrix.element[n][n+2]+=f2*drJI.x*dr.x*dJI.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+3]+=f2*drJI.x*dr.x*dJI.z*dr.z;                    // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*drJI.y*dr.x*dJI.y*dr.x+f1*drJI.y*dJI.y;  // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*drJI.y*dr.x*dJI.y*dr.y;                  // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*drJI.y*dr.x*dJI.z*dr.z;                  // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*drJI.y*dr.y*dJI.y*dr.y+f1*drJI.y*dJI.y;  // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*drJI.y*dr.y*dJI.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*drJI.z*dr.z*dJI.z*dr.z+f1*drJI.z*dJI.z;  // zzzz
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
      switch(Dimension)
      {
        case 3:
          Gradient[n]+=f1*dr.z*dr.z;
        case 2:
          Gradient[n]+=f1*dr.y*dr.y;
        case 1:
          Gradient[n]+=f1*dr.x*dr.x;
          break;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=f1*dr.z*dr.z/Dimension;
            case 2:
              Gradient[n+1]+=f1*dr.y*dr.y/Dimension;
            case 1:
              Gradient[n]+=f1*dr.x*dr.x/Dimension;
              break;
          }
          break;
        case ANISOTROPIC:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=f1*dr.z*dr.z;
            case 2:
              Gradient[n+1]+=f1*dr.y*dr.y;
            case 1:
              Gradient[n]+=f1*dr.x*dr.x;
              break;
          }
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=f1*dr.x*dr.z;
              Gradient[n+4]+=f1*dr.y*dr.z;
              Gradient[n+5]+=f1*dr.z*dr.z;
            case 2:
              Gradient[n+1]+=f1*dr.x*dr.y;
              Gradient[n+Dimension]+=f1*dr.y*dr.y;
            case 1:
             Gradient[n]+=f1*dr.x*dr.x;
             break;
          }
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]+=f1*dr.x*dr.x;
              Gradient[n+1]+=f1*dr.y*dr.y;
              Gradient[n+2]+=f1*dr.y*dr.z;
              Gradient[n+3]+=f1*dr.z*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]+=f1*dr.x*dr.x;
              Gradient[n+1]+=f1*dr.x*dr.y;
              Gradient[n+2]+=f1*dr.y*dr.y;
              Gradient[n+3]+=f1*dr.z*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
            default:
              Gradient[n]+=f1*dr.x*dr.x;
              Gradient[n+1]+=f1*dr.x*dr.z;
              Gradient[n+2]+=f1*dr.y*dr.y;
              Gradient[n+3]+=f1*dr.z*dr.z;
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

void GradientStrainI(REAL *Gradient,REAL f1,VECTOR dr,VECTOR posA,VECTOR comA)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      switch(Dimension)
      {
        case 3:
          Gradient[n]-=f1*(posA.z-comA.z)*dr.z;
        case 2:
          Gradient[n]-=f1*(posA.y-comA.y)*dr.y;
        case 1:
          Gradient[n]-=f1*(posA.x-comA.x)*dr.x;
          break;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]-=(posA.z-comA.z)*f1*dr.z;
            case 2:
              Gradient[n+1]-=(posA.y-comA.y)*f1*dr.y;
            case 1:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              break;
          }
          break;
        case REGULAR:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]-=0.5*f1*((posA.z-comA.z)*dr.x+(posA.x-comA.x)*dr.z);
              Gradient[n+4]-=0.5*f1*((posA.z-comA.z)*dr.y+(posA.y-comA.y)*dr.z);
              Gradient[n+5]-=f1*(posA.z-comA.z)*dr.z;
            case 2:
              Gradient[n+1]-=0.5*f1*((posA.y-comA.y)*dr.x+(posA.x-comA.x)*dr.y);
              Gradient[n+3]-=f1*(posA.y-comA.y)*dr.y;
            case 1:
              Gradient[n]-=f1*(posA.x-comA.x)*dr.x;
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]-=f1*(posA.z-comA.z)*dr.x;
              Gradient[n+4]-=f1*(posA.z-comA.z)*dr.y;
              Gradient[n+5]-=f1*(posA.z-comA.z)*dr.z;
            case 2:
              Gradient[n+1]-=f1*(posA.y-comA.y)*dr.x;
              Gradient[n+Dimension]-=f1*(posA.y-comA.y)*dr.y;
            case 1:
              Gradient[n]-=f1*(posA.x-comA.x)*dr.x;
              break;
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+2]-=0.5*f1*((posA.z-comA.z)*dr.y+(posA.y-comA.y)*dr.z);
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=0.5*f1*((posA.y-comA.y)*dr.x+(posA.x-comA.x)*dr.y);
              Gradient[n+2]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=0.5*f1*((posA.z-comA.z)*dr.x+(posA.x-comA.x)*dr.z);
              Gradient[n+2]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            default:
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+2]-=f1*(posA.z-comA.z)*dr.y;
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=f1*(posA.z-comA.z)*dr.x;
              Gradient[n+2]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=f1*(posA.y-comA.y)*dr.x;
              Gradient[n+2]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            default:
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

void GradientStrainJ(REAL *Gradient,REAL f1,VECTOR dr,VECTOR posB,VECTOR comB)
{
  int n;
  REAL temp1,temp2,temp3;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      switch(Dimension)
      {
        case 3:
          Gradient[n]+=f1*(posB.z-comB.z)*dr.z;
        case 2:
          Gradient[n]+=f1*(posB.y-comB.y)*dr.y;
        case 1:
          Gradient[n]+=f1*(posB.x-comB.x)*dr.x;
          break;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=(posB.z-comB.z)*f1*dr.z;
            case 2:
              Gradient[n+1]+=(posB.y-comB.y)*f1*dr.y;
            case 1:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              break;
          }
          break;
        case REGULAR:
          temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
          temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
          temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=temp2;
              Gradient[n+4]+=temp3;
              Gradient[n+5]+=(posB.z-comB.z)*f1*dr.z;
            case 2:
              Gradient[n+1]+=temp1;
              Gradient[n+3]+=(posB.y-comB.y)*f1*dr.y;
            case 1:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          temp1=(posB.y-comB.y)*f1*dr.x;
          temp2=(posB.z-comB.z)*f1*dr.x;
          temp3=(posB.z-comB.z)*f1*dr.y;
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=temp2;
              Gradient[n+4]+=temp3;
              Gradient[n+5]+=(posB.z-comB.z)*f1*dr.z;
            case 2:
              Gradient[n+1]+=temp1;
              Gradient[n+Dimension]+=(posB.y-comB.y)*f1*dr.y;
            case 1:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              break;
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+2]+=temp3;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=temp1;
              Gradient[n+2]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=temp2;
              Gradient[n+2]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            default:
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+2]+=(posB.z-comB.z)*f1*dr.y;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=(posB.z-comB.z)*f1*dr.x;
              Gradient[n+2]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=(posB.y-comB.y)*f1*dr.x;
              Gradient[n+2]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            default:
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

void PreComputeRotationDerivatives(void)
{
  int i,m,l,A,Type,index;
  VECTOR com,pos,p;
  REAL_MATRIX3x3 RotationMatrix;
  REAL_MATRIX3x3 DRotationMatrixIX,DRotationMatrixIY,DRotationMatrixIZ;
  REAL_MATRIX3x3 DDRotationMatrixIAX,DDRotationMatrixIBY,DDRotationMatrixICZ;
  REAL_MATRIX3x3 DDRotationMatrixIAY,DDRotationMatrixIAZ,DDRotationMatrixIBZ;

  index=0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        p=Adsorbates[CurrentSystem][m].Groups[l].EulerAxis;
        RotationMatrix=ComputeRotationMatrix(p);
        com=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition;

        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          Adsorbates[CurrentSystem][m].Atoms[A].HessianAtomIndex=index;
          //index=Adsorbates[CurrentSystem][m].Atoms[A].HessianAtomIndex;

          pos=Components[Type].Positions[A];

          DRotationMatrixIX=ComputeRotationMatrixDerivativeX(Adsorbates[CurrentSystem][m].Groups[l].EulerAxis);
          DRotationMatrixIY=ComputeRotationMatrixDerivativeY(Adsorbates[CurrentSystem][m].Groups[l].EulerAxis);
          DRotationMatrixIZ=ComputeRotationMatrixDerivativeZ(Adsorbates[CurrentSystem][m].Groups[l].EulerAxis);

          DVecX[index].x=(DRotationMatrixIX.ax*pos.x+DRotationMatrixIX.bx*pos.y+DRotationMatrixIX.cx*pos.z);
          DVecX[index].y=(DRotationMatrixIX.ay*pos.x+DRotationMatrixIX.by*pos.y+DRotationMatrixIX.cy*pos.z);
          DVecX[index].z=(DRotationMatrixIX.az*pos.x+DRotationMatrixIX.bz*pos.y+DRotationMatrixIX.cz*pos.z);
          DVecY[index].x=(DRotationMatrixIY.ax*pos.x+DRotationMatrixIY.bx*pos.y+DRotationMatrixIY.cx*pos.z);
          DVecY[index].y=(DRotationMatrixIY.ay*pos.x+DRotationMatrixIY.by*pos.y+DRotationMatrixIY.cy*pos.z);
          DVecY[index].z=(DRotationMatrixIY.az*pos.x+DRotationMatrixIY.bz*pos.y+DRotationMatrixIY.cz*pos.z);
          DVecZ[index].x=(DRotationMatrixIZ.ax*pos.x+DRotationMatrixIZ.bx*pos.y+DRotationMatrixIZ.cx*pos.z);
          DVecZ[index].y=(DRotationMatrixIZ.ay*pos.x+DRotationMatrixIZ.by*pos.y+DRotationMatrixIZ.cy*pos.z);
          DVecZ[index].z=(DRotationMatrixIZ.az*pos.x+DRotationMatrixIZ.bz*pos.y+DRotationMatrixIZ.cz*pos.z);

          DDRotationMatrixIAX=ComputeRotationMatrixSecondDerivativeAX(Adsorbates[CurrentSystem][m].Groups[l].EulerAxis);
          DDRotationMatrixIBY=ComputeRotationMatrixSecondDerivativeBY(Adsorbates[CurrentSystem][m].Groups[l].EulerAxis);
          DDRotationMatrixICZ=ComputeRotationMatrixSecondDerivativeCZ(Adsorbates[CurrentSystem][m].Groups[l].EulerAxis);

          DDRotationMatrixIAY=ComputeRotationMatrixSecondDerivativeAY(Adsorbates[CurrentSystem][m].Groups[l].EulerAxis);
          DDRotationMatrixIAZ=ComputeRotationMatrixSecondDerivativeAZ(Adsorbates[CurrentSystem][m].Groups[l].EulerAxis);
          DDRotationMatrixIBZ=ComputeRotationMatrixSecondDerivativeBZ(Adsorbates[CurrentSystem][m].Groups[l].EulerAxis);

          DDVecAX[index].x=(DDRotationMatrixIAX.ax*pos.x+DDRotationMatrixIAX.bx*pos.y+DDRotationMatrixIAX.cx*pos.z);
          DDVecAX[index].y=(DDRotationMatrixIAX.ay*pos.x+DDRotationMatrixIAX.by*pos.y+DDRotationMatrixIAX.cy*pos.z);
          DDVecAX[index].z=(DDRotationMatrixIAX.az*pos.x+DDRotationMatrixIAX.bz*pos.y+DDRotationMatrixIAX.cz*pos.z);

          DDVecBY[index].x=(DDRotationMatrixIBY.ax*pos.x+DDRotationMatrixIBY.bx*pos.y+DDRotationMatrixIBY.cx*pos.z);
          DDVecBY[index].y=(DDRotationMatrixIBY.ay*pos.x+DDRotationMatrixIBY.by*pos.y+DDRotationMatrixIBY.cy*pos.z);
          DDVecBY[index].z=(DDRotationMatrixIBY.az*pos.x+DDRotationMatrixIBY.bz*pos.y+DDRotationMatrixIBY.cz*pos.z);

          DDVecCZ[index].x=(DDRotationMatrixICZ.ax*pos.x+DDRotationMatrixICZ.bx*pos.y+DDRotationMatrixICZ.cx*pos.z);
          DDVecCZ[index].y=(DDRotationMatrixICZ.ay*pos.x+DDRotationMatrixICZ.by*pos.y+DDRotationMatrixICZ.cy*pos.z);
          DDVecCZ[index].z=(DDRotationMatrixICZ.az*pos.x+DDRotationMatrixICZ.bz*pos.y+DDRotationMatrixICZ.cz*pos.z);

          DDVecAY[index].x=(DDRotationMatrixIAY.ax*pos.x+DDRotationMatrixIAY.bx*pos.y+DDRotationMatrixIAY.cx*pos.z);
          DDVecAY[index].y=(DDRotationMatrixIAY.ay*pos.x+DDRotationMatrixIAY.by*pos.y+DDRotationMatrixIAY.cy*pos.z);
          DDVecAY[index].z=(DDRotationMatrixIAY.az*pos.x+DDRotationMatrixIAY.bz*pos.y+DDRotationMatrixIAY.cz*pos.z);

          DDVecAZ[index].x=(DDRotationMatrixIAZ.ax*pos.x+DDRotationMatrixIAZ.bx*pos.y+DDRotationMatrixIAZ.cx*pos.z);
          DDVecAZ[index].y=(DDRotationMatrixIAZ.ay*pos.x+DDRotationMatrixIAZ.by*pos.y+DDRotationMatrixIAZ.cy*pos.z);
          DDVecAZ[index].z=(DDRotationMatrixIAZ.az*pos.x+DDRotationMatrixIAZ.bz*pos.y+DDRotationMatrixIAZ.cz*pos.z);

          DDVecBZ[index].x=(DDRotationMatrixIBZ.ax*pos.x+DDRotationMatrixIBZ.bx*pos.y+DDRotationMatrixIBZ.cx*pos.z);
          DDVecBZ[index].y=(DDRotationMatrixIBZ.ay*pos.x+DDRotationMatrixIBZ.by*pos.y+DDRotationMatrixIBZ.cy*pos.z);
          DDVecBZ[index].z=(DDRotationMatrixIBZ.az*pos.x+DDRotationMatrixIBZ.bz*pos.y+DDRotationMatrixIBZ.cz*pos.z);

          index++;
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        p=Cations[CurrentSystem][m].Groups[l].EulerAxis;
        RotationMatrix=ComputeRotationMatrix(p);
        com=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition;

        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          Cations[CurrentSystem][m].Atoms[A].HessianAtomIndex=index;
          //index=Cations[CurrentSystem][m].Atoms[A].HessianAtomIndex;

          pos=Components[Type].Positions[A];

          DRotationMatrixIX=ComputeRotationMatrixDerivativeX(Cations[CurrentSystem][m].Groups[l].EulerAxis);
          DRotationMatrixIY=ComputeRotationMatrixDerivativeY(Cations[CurrentSystem][m].Groups[l].EulerAxis);
          DRotationMatrixIZ=ComputeRotationMatrixDerivativeZ(Cations[CurrentSystem][m].Groups[l].EulerAxis);

          DVecX[index].x=(DRotationMatrixIX.ax*pos.x+DRotationMatrixIX.bx*pos.y+DRotationMatrixIX.cx*pos.z);
          DVecX[index].y=(DRotationMatrixIX.ay*pos.x+DRotationMatrixIX.by*pos.y+DRotationMatrixIX.cy*pos.z);
          DVecX[index].z=(DRotationMatrixIX.az*pos.x+DRotationMatrixIX.bz*pos.y+DRotationMatrixIX.cz*pos.z);
          DVecY[index].x=(DRotationMatrixIY.ax*pos.x+DRotationMatrixIY.bx*pos.y+DRotationMatrixIY.cx*pos.z);
          DVecY[index].y=(DRotationMatrixIY.ay*pos.x+DRotationMatrixIY.by*pos.y+DRotationMatrixIY.cy*pos.z);
          DVecY[index].z=(DRotationMatrixIY.az*pos.x+DRotationMatrixIY.bz*pos.y+DRotationMatrixIY.cz*pos.z);
          DVecZ[index].x=(DRotationMatrixIZ.ax*pos.x+DRotationMatrixIZ.bx*pos.y+DRotationMatrixIZ.cx*pos.z);
          DVecZ[index].y=(DRotationMatrixIZ.ay*pos.x+DRotationMatrixIZ.by*pos.y+DRotationMatrixIZ.cy*pos.z);
          DVecZ[index].z=(DRotationMatrixIZ.az*pos.x+DRotationMatrixIZ.bz*pos.y+DRotationMatrixIZ.cz*pos.z);

          DDRotationMatrixIAX=ComputeRotationMatrixSecondDerivativeAX(Cations[CurrentSystem][m].Groups[l].EulerAxis);
          DDRotationMatrixIBY=ComputeRotationMatrixSecondDerivativeBY(Cations[CurrentSystem][m].Groups[l].EulerAxis);
          DDRotationMatrixICZ=ComputeRotationMatrixSecondDerivativeCZ(Cations[CurrentSystem][m].Groups[l].EulerAxis);

          DDRotationMatrixIAY=ComputeRotationMatrixSecondDerivativeAY(Cations[CurrentSystem][m].Groups[l].EulerAxis);
          DDRotationMatrixIAZ=ComputeRotationMatrixSecondDerivativeAZ(Cations[CurrentSystem][m].Groups[l].EulerAxis);
          DDRotationMatrixIBZ=ComputeRotationMatrixSecondDerivativeBZ(Cations[CurrentSystem][m].Groups[l].EulerAxis);

          DDVecAX[index].x=(DDRotationMatrixIAX.ax*pos.x+DDRotationMatrixIAX.bx*pos.y+DDRotationMatrixIAX.cx*pos.z);
          DDVecAX[index].y=(DDRotationMatrixIAX.ay*pos.x+DDRotationMatrixIAX.by*pos.y+DDRotationMatrixIAX.cy*pos.z);
          DDVecAX[index].z=(DDRotationMatrixIAX.az*pos.x+DDRotationMatrixIAX.bz*pos.y+DDRotationMatrixIAX.cz*pos.z);

          DDVecBY[index].x=(DDRotationMatrixIBY.ax*pos.x+DDRotationMatrixIBY.bx*pos.y+DDRotationMatrixIBY.cx*pos.z);
          DDVecBY[index].y=(DDRotationMatrixIBY.ay*pos.x+DDRotationMatrixIBY.by*pos.y+DDRotationMatrixIBY.cy*pos.z);
          DDVecBY[index].z=(DDRotationMatrixIBY.az*pos.x+DDRotationMatrixIBY.bz*pos.y+DDRotationMatrixIBY.cz*pos.z);

          DDVecCZ[index].x=(DDRotationMatrixICZ.ax*pos.x+DDRotationMatrixICZ.bx*pos.y+DDRotationMatrixICZ.cx*pos.z);
          DDVecCZ[index].y=(DDRotationMatrixICZ.ay*pos.x+DDRotationMatrixICZ.by*pos.y+DDRotationMatrixICZ.cy*pos.z);
          DDVecCZ[index].z=(DDRotationMatrixICZ.az*pos.x+DDRotationMatrixICZ.bz*pos.y+DDRotationMatrixICZ.cz*pos.z);

          DDVecAY[index].x=(DDRotationMatrixIAY.ax*pos.x+DDRotationMatrixIAY.bx*pos.y+DDRotationMatrixIAY.cx*pos.z);
          DDVecAY[index].y=(DDRotationMatrixIAY.ay*pos.x+DDRotationMatrixIAY.by*pos.y+DDRotationMatrixIAY.cy*pos.z);
          DDVecAY[index].z=(DDRotationMatrixIAY.az*pos.x+DDRotationMatrixIAY.bz*pos.y+DDRotationMatrixIAY.cz*pos.z);

          DDVecAZ[index].x=(DDRotationMatrixIAZ.ax*pos.x+DDRotationMatrixIAZ.bx*pos.y+DDRotationMatrixIAZ.cx*pos.z);
          DDVecAZ[index].y=(DDRotationMatrixIAZ.ay*pos.x+DDRotationMatrixIAZ.by*pos.y+DDRotationMatrixIAZ.cy*pos.z);
          DDVecAZ[index].z=(DDRotationMatrixIAZ.az*pos.x+DDRotationMatrixIAZ.bz*pos.y+DDRotationMatrixIAZ.cz*pos.z);

          DDVecBZ[index].x=(DDRotationMatrixIBZ.ax*pos.x+DDRotationMatrixIBZ.bx*pos.y+DDRotationMatrixIBZ.cx*pos.z);
          DDVecBZ[index].y=(DDRotationMatrixIBZ.ay*pos.x+DDRotationMatrixIBZ.by*pos.y+DDRotationMatrixIBZ.cy*pos.z);
          DDVecBZ[index].z=(DDRotationMatrixIBZ.az*pos.x+DDRotationMatrixIBZ.bz*pos.y+DDRotationMatrixIBZ.cz*pos.z);

          index++;
        }
      }
    }
  }
}


void ComputeInterVDWMolecularHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                     REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int I,J,i,j,ig,jg,ia,ja;
  int typeA,typeB;
  int TypeMolA,TypeMolB;
  REAL rr;
  REAL energy,f1,f2;
  VECTOR posA,posB,dr;
  int index1,index2;
  INT_VECTOR3 index_i,index_j;
  INT_VECTOR3 index_i2,index_j2;
  VECTOR comA,comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3,start;
  REAL ReplicaFactor;
  int RigidI,RigidJ;
  REAL scalingA,scalingB;

  f1=f2=0.0;
  index1=index2=0;

  // first loop over adsorbate molecules
  for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
  {
    TypeMolA=Adsorbates[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      RigidI=Components[TypeMolA].Groups[ig].Rigid;
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(RigidI)
        {
          index_i=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndex;
          index_i2=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndexOrientation;
          comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;
          index1=Adsorbates[CurrentSystem][I].Atoms[i].HessianAtomIndex;
        }
        else
        {
          index_i=Adsorbates[CurrentSystem][I].Atoms[i].HessianIndex;
          index_i2=UNDEFINED_INT_VECTOR3;
          comA=Adsorbates[CurrentSystem][I].Atoms[i].Position;
          index1=-1;
        }

        typeA=Adsorbates[CurrentSystem][I].Atoms[i].Type;
        posA=Adsorbates[CurrentSystem][I].Atoms[i].Position;
        scalingA=Adsorbates[CurrentSystem][I].Atoms[i].CFVDWScalingParameter;

        // second loop over adsorbates
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if(ncell==0) start=I+1;
              else start=0;

              for(J=start;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
              {
                TypeMolB=Adsorbates[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];

                    if(RigidJ)
                    {
                      index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                      index2=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=UNDEFINED_INT_VECTOR3;
                      comB=Adsorbates[CurrentSystem][J].Atoms[j].Position;
                      index2=-1;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

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


                    if(rr<CutOffVDWSquared)
                    {
                      if(ncell==0) ReplicaFactor=1.0;
                      else ReplicaFactor=0.5;

                      scalingB=Adsorbates[CurrentSystem][J].Atoms[j].CFVDWScalingParameter;

                      PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2,scalingA*scalingB);

                      *Energy+=ReplicaFactor*energy;

                      //if((index_i<0)&&(index_i2<0)&&(index_j<0)&&(index_j2<0)) continue;

                      StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
                      StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
                      StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

                      StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
                      StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
                      StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

                      StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
                      StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
                      StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=ReplicaFactor*(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=ReplicaFactor*temp1;
                        StrainDerivative->az-=ReplicaFactor*temp2;
                        StrainDerivative->bx-=ReplicaFactor*temp1;
                        StrainDerivative->by-=ReplicaFactor*(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=ReplicaFactor*temp3;
                        StrainDerivative->cx-=ReplicaFactor*temp2;
                        StrainDerivative->cy-=ReplicaFactor*temp3;
                        StrainDerivative->cz-=ReplicaFactor*(posA.z-comA.z)*f1*dr.z;
                      }

                      if(RigidJ)
                      {
                        temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                        temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                        temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                        StrainDerivative->ax+=ReplicaFactor*(posB.x-comB.x)*f1*dr.x;
                        StrainDerivative->ay+=ReplicaFactor*temp1;
                        StrainDerivative->az+=ReplicaFactor*temp2;
                        StrainDerivative->bx+=ReplicaFactor*temp1;
                        StrainDerivative->by+=ReplicaFactor*(posB.y-comB.y)*f1*dr.y;
                        StrainDerivative->bz+=ReplicaFactor*temp3;
                        StrainDerivative->cx+=ReplicaFactor*temp2;
                        StrainDerivative->cy+=ReplicaFactor*temp3;
                        StrainDerivative->cz+=ReplicaFactor*(posB.z-comB.z)*f1*dr.z;
                      }

                      // add contribution to the first derivatives
                      if(ComputeGradient)
                      {
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2.x>=0) Gradient[index_i2.x]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                          if(index_i2.y>=0) Gradient[index_i2.y]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                          if(index_i2.z>=0) Gradient[index_i2.z]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                        }

                        GradientStrain(Gradient,ReplicaFactor*f1,dr);

                        if(ncell==0)
                        {
                          if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                          if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                          if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                          if(RigidJ)
                          {
                            GradientStrainJ(Gradient,f1,dr,posB,comB);

                            // add contribution to the first derivatives
                            if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                            if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                            if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                          }
                        }
                      }

                      if(ComputeHessian)
                      {
                        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                           f1,f2,dr,ReplicaFactor);
                        HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                           f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,ReplicaFactor*f1,ReplicaFactor*f2,dr,posA,comA,posB,comB,
                                                  RigidI,RigidJ);

                        HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,
                                            ReplicaFactor*f1,ReplicaFactor*f2,posA,comA,posB,comB,dr);

                        HessianAtomicStrainStrain(HessianMatrix,ReplicaFactor*f1,ReplicaFactor*f2,dr,posA,comA,posB,comB);
                      }
                    }
                  }
                }
              }

              for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
              {
                TypeMolB=Cations[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];

                    if(RigidJ)
                    {
                      index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                      index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=UNDEFINED_INT_VECTOR3;
                      comB=Cations[CurrentSystem][J].Atoms[j].Position;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

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


                    if(rr<CutOffVDWSquared)
                    {
                      scalingB=Adsorbates[CurrentSystem][J].Atoms[j].CFVDWScalingParameter;

                      PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2,scalingA*scalingB);

                      *Energy+=energy;

                      //if((index_i<0)&&(index_i2<0)&&(index_j<0)&&(index_j2<0)) continue;

                      StrainDerivative->ax+=f1*dr.x*dr.x;
                      StrainDerivative->ay+=f1*dr.x*dr.y;
                      StrainDerivative->az+=f1*dr.x*dr.z;

                      StrainDerivative->bx+=f1*dr.y*dr.x;
                      StrainDerivative->by+=f1*dr.y*dr.y;
                      StrainDerivative->bz+=f1*dr.y*dr.z;

                      StrainDerivative->cx+=f1*dr.z*dr.x;
                      StrainDerivative->cy+=f1*dr.z*dr.y;
                      StrainDerivative->cz+=f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=temp1;
                        StrainDerivative->az-=temp2;
                        StrainDerivative->bx-=temp1;
                        StrainDerivative->by-=(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=temp3;
                        StrainDerivative->cx-=temp2;
                        StrainDerivative->cy-=temp3;
                        StrainDerivative->cz-=(posA.z-comA.z)*f1*dr.z;
                      }

                      if(RigidJ)
                      {
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
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2.x>=0) Gradient[index_i2.x]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                          if(index_i2.y>=0) Gradient[index_i2.y]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                          if(index_i2.z>=0) Gradient[index_i2.z]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                        }

                        GradientStrain(Gradient,f1,dr);

                        if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                        if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                        if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                        if(RigidJ)
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

                        HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                           f1,f2,dr,1.0);
                        HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                           f1,f2,dr,1.0);

                        HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,f1,f2,dr,posA,comA,posB,comB,
                                                  RigidI,RigidJ);

                        HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,
                                            f1,f2,posA,comA,posB,comB,dr);

                        HessianAtomicStrainStrain(HessianMatrix,f1,f2,dr,posA,comA,posB,comB);
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

  // first loop over cation molecules
  for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
  {
    TypeMolA=Cations[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      RigidI=Components[TypeMolA].Groups[ig].Rigid;
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(RigidI)
        {
          index_i=Cations[CurrentSystem][I].Groups[ig].HessianIndex;
          index_i2=Cations[CurrentSystem][I].Groups[ig].HessianIndexOrientation;
          comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

          index1=Cations[CurrentSystem][I].Atoms[i].HessianAtomIndex;
        }
        else
        {
          index_i=Cations[CurrentSystem][I].Atoms[i].HessianIndex;
          index_i2=UNDEFINED_INT_VECTOR3;
          comA=Cations[CurrentSystem][I].Atoms[i].Position;
        }

        typeA=Cations[CurrentSystem][I].Atoms[i].Type;
        posA=Cations[CurrentSystem][I].Atoms[i].Position;

        scalingA=Cations[CurrentSystem][I].Atoms[i].CFVDWScalingParameter;

        // second loop over cations
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if(ncell==0) start=I+1;
              else start=0;

              for(J=start;J<NumberOfCationMolecules[CurrentSystem];J++)
              {
                TypeMolB=Cations[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];

                    if(RigidJ)
                    {
                      index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                      index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=UNDEFINED_INT_VECTOR3;
                      comB=Cations[CurrentSystem][J].Atoms[j].Position;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

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


                    if(rr<CutOffVDWSquared)
                    {
                      if(ncell==0) ReplicaFactor=1.0;
                      else ReplicaFactor=0.5;

                      scalingB=Cations[CurrentSystem][J].Atoms[j].CFVDWScalingParameter;

                      PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2,scalingA*scalingB);

                      *Energy+=ReplicaFactor*energy;

                      //if((index_i<0)&&(index_i2<0)&&(index_j<0)&&(index_j2<0)) continue;

                      StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
                      StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
                      StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

                      StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
                      StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
                      StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

                      StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
                      StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
                      StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=ReplicaFactor*(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=ReplicaFactor*temp1;
                        StrainDerivative->az-=ReplicaFactor*temp2;
                        StrainDerivative->bx-=ReplicaFactor*temp1;
                        StrainDerivative->by-=ReplicaFactor*(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=ReplicaFactor*temp3;
                        StrainDerivative->cx-=ReplicaFactor*temp2;
                        StrainDerivative->cy-=ReplicaFactor*temp3;
                        StrainDerivative->cz-=ReplicaFactor*(posA.z-comA.z)*f1*dr.z;
                      }

                      if(RigidJ)
                      {
                        temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                        temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                        temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                        StrainDerivative->ax+=ReplicaFactor*(posB.x-comB.x)*f1*dr.x;
                        StrainDerivative->ay+=ReplicaFactor*temp1;
                        StrainDerivative->az+=ReplicaFactor*temp2;
                        StrainDerivative->bx+=ReplicaFactor*temp1;
                        StrainDerivative->by+=ReplicaFactor*(posB.y-comB.y)*f1*dr.y;
                        StrainDerivative->bz+=ReplicaFactor*temp3;
                        StrainDerivative->cx+=ReplicaFactor*temp2;
                        StrainDerivative->cy+=ReplicaFactor*temp3;
                        StrainDerivative->cz+=ReplicaFactor*(posB.z-comB.z)*f1*dr.z;
                      }

                      // add contribution to the first derivatives
                      if(ComputeGradient)
                      {
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2.x>=0) Gradient[index_i2.x]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                          if(index_i2.y>=0) Gradient[index_i2.y]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                          if(index_i2.z>=0) Gradient[index_i2.z]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                        }

                        GradientStrain(Gradient,ReplicaFactor*f1,dr);

                        if(ncell==0)
                        {
                          if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                          if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                          if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                          if(RigidJ)
                          {
                            GradientStrainJ(Gradient,f1,dr,posB,comB);

                            // add contribution to the first derivatives
                            if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                            if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                            if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                          }
                        }
                      }

                      if(ComputeHessian)
                      {
                        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                           f1,f2,dr,ReplicaFactor);
                        HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                           f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,ReplicaFactor*f1,ReplicaFactor*f2,dr,posA,comA,posB,comB,
                                                  RigidI,RigidJ);

                        HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,
                                            ReplicaFactor*f1,ReplicaFactor*f2,posA,comA,posB,comB,dr);

                        HessianAtomicStrainStrain(HessianMatrix,ReplicaFactor*f1,ReplicaFactor*f2,dr,posA,comA,posB,comB);
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
}

void ComputeInterChargeChargeMolecularHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                              REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int I,J,i,j,ig,jg,ia,ja;
  int typeA,typeB;
  int TypeMolA,TypeMolB;
  REAL ChargeA,ChargeB,rr;
  REAL U,f1,f2;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j;
  INT_VECTOR3 index_i2,index_j2;
  VECTOR comA,comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3,start;
  REAL ReplicaFactor;
  int index1,index2;
  int RigidI,RigidJ;
  REAL scalingA,scalingB;

  if(ChargeMethod==NONE) return;

  f1=f2=0.0;
  index1=index2=0;
  // first loop over adsorbate molecules
  for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
  {
    TypeMolA=Adsorbates[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      RigidI=Components[TypeMolA].Groups[ig].Rigid;
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(RigidI)
        {
          index_i=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndex;
          index_i2=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndexOrientation;
          comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

          index1=Adsorbates[CurrentSystem][I].Atoms[i].HessianAtomIndex;
        }
        else
        {
          index_i=Adsorbates[CurrentSystem][I].Atoms[i].HessianIndex;
          index_i2=UNDEFINED_INT_VECTOR3;
          comA=Adsorbates[CurrentSystem][I].Atoms[i].Position;
          index1=-1;
        }

        typeA=Adsorbates[CurrentSystem][I].Atoms[i].Type;
        posA=Adsorbates[CurrentSystem][I].Atoms[i].Position;
        scalingA=Adsorbates[CurrentSystem][I].Atoms[i].CFChargeScalingParameter;
        ChargeA=scalingA*Adsorbates[CurrentSystem][I].Atoms[i].Charge;

        // second loop over adsorbates
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if(ncell==0) start=I+1;
              else start=0;

              for(J=start;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
              {
                TypeMolB=Adsorbates[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];

                    if(RigidJ)
                    {
                      index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                      index2=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=UNDEFINED_INT_VECTOR3;
                      comB=Adsorbates[CurrentSystem][J].Atoms[j].Position;
                      index2=-1;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

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

                      if(ncell==0) ReplicaFactor=1.0;
                      else ReplicaFactor=0.5;

                      PotentialSecondDerivativeCoulombic(ChargeA,ChargeB,rr,&U,&f1,&f2);

                      *Energy+=ReplicaFactor*U;

                      StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
                      StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
                      StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

                      StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
                      StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
                      StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

                      StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
                      StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
                      StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=ReplicaFactor*(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=ReplicaFactor*temp1;
                        StrainDerivative->az-=ReplicaFactor*temp2;
                        StrainDerivative->bx-=ReplicaFactor*temp1;
                        StrainDerivative->by-=ReplicaFactor*(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=ReplicaFactor*temp3;
                        StrainDerivative->cx-=ReplicaFactor*temp2;
                        StrainDerivative->cy-=ReplicaFactor*temp3;
                        StrainDerivative->cz-=ReplicaFactor*(posA.z-comA.z)*f1*dr.z;
                      }


                      if(RigidJ)
                      {
                        temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                        temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                        temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                        StrainDerivative->ax+=ReplicaFactor*(posB.x-comB.x)*f1*dr.x;
                        StrainDerivative->ay+=ReplicaFactor*temp1;
                        StrainDerivative->az+=ReplicaFactor*temp2;
                        StrainDerivative->bx+=ReplicaFactor*temp1;
                        StrainDerivative->by+=ReplicaFactor*(posB.y-comB.y)*f1*dr.y;
                        StrainDerivative->bz+=ReplicaFactor*temp3;
                        StrainDerivative->cx+=ReplicaFactor*temp2;
                        StrainDerivative->cy+=ReplicaFactor*temp3;
                        StrainDerivative->cz+=ReplicaFactor*(posB.z-comB.z)*f1*dr.z;
                      }

                      // add contribution to the first derivatives
                      if(ComputeGradient)
                      {
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2.x>=0) Gradient[index_i2.x]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                          if(index_i2.y>=0) Gradient[index_i2.y]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                          if(index_i2.z>=0) Gradient[index_i2.z]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                        }

                        GradientStrain(Gradient,ReplicaFactor*f1,dr);

                        if(ncell==0)
                        {
                          if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                          if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                          if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                          if(RigidJ)
                          {
                            GradientStrainJ(Gradient,f1,dr,posB,comB);

                            // add contribution to the first derivatives
                            if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                            if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                            if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                          }
                        }
                      }

                      if(ComputeHessian)
                      {
                        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);
                        HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,ReplicaFactor*f1,ReplicaFactor*f2,dr,posA,comA,posB,comB,
                                                   RigidI,RigidJ);

                        HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,
                                         ReplicaFactor*f1,ReplicaFactor*f2,posA,comA,posB,comB,dr);

                        HessianAtomicStrainStrain(HessianMatrix,ReplicaFactor*f1,ReplicaFactor*f2,dr,posA,comA,posB,comB);
                      }
                    }
                  }
                }
              }
              for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
              {
                TypeMolB=Cations[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];

                    if(RigidJ)
                    {
                      index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                      index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=UNDEFINED_INT_VECTOR3;
                      comB=Cations[CurrentSystem][J].Atoms[j].Position;
                      index2=-1;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

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
                      StrainDerivative->ay+=f1*dr.x*dr.y;
                      StrainDerivative->az+=f1*dr.x*dr.z;

                      StrainDerivative->bx+=f1*dr.y*dr.x;
                      StrainDerivative->by+=f1*dr.y*dr.y;
                      StrainDerivative->bz+=f1*dr.y*dr.z;

                      StrainDerivative->cx+=f1*dr.z*dr.x;
                      StrainDerivative->cy+=f1*dr.z*dr.y;
                      StrainDerivative->cz+=f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=temp1;
                        StrainDerivative->az-=temp2;
                        StrainDerivative->bx-=temp1;
                        StrainDerivative->by-=(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=temp3;
                        StrainDerivative->cx-=temp2;
                        StrainDerivative->cy-=temp3;
                        StrainDerivative->cz-=(posA.z-comA.z)*f1*dr.z;
                      }


                      if(RigidJ)
                      {
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
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2.x>=0) Gradient[index_i2.x]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                          if(index_i2.y>=0) Gradient[index_i2.y]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                          if(index_i2.z>=0) Gradient[index_i2.z]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                        }

                        GradientStrain(Gradient,f1,dr);

                        if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                        if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                        if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                        if(RigidJ)
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

                        HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,1.0);
                        HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,1.0);

                        HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,f1,f2,dr,posA,comA,posB,comB,
                                                   RigidI,RigidJ);

                        HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,
                             f1,f2,posA,comA,posB,comB,dr);

                        HessianAtomicStrainStrain(HessianMatrix,f1,f2,dr,posA,comA,posB,comB);
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


  // Cation-Cation
  for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
  {
    TypeMolA=Cations[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      RigidI=Components[TypeMolA].Groups[ig].Rigid;
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(RigidI)
        {
          index_i=Cations[CurrentSystem][I].Groups[ig].HessianIndex;
          index_i2=Cations[CurrentSystem][I].Groups[ig].HessianIndexOrientation;
          comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

          index1=Cations[CurrentSystem][I].Atoms[i].HessianAtomIndex;
        }
        else
        {
          index_i=Cations[CurrentSystem][I].Atoms[i].HessianIndex;
          index_i2=UNDEFINED_INT_VECTOR3;
          comA=Cations[CurrentSystem][I].Atoms[i].Position;
          index1=-1;
        }

        typeA=Cations[CurrentSystem][I].Atoms[i].Type;
        posA=Cations[CurrentSystem][I].Atoms[i].Position;
        scalingA=Cations[CurrentSystem][I].Atoms[i].CFChargeScalingParameter;
        ChargeA=scalingA*Cations[CurrentSystem][I].Atoms[i].Charge;

        // second loop over adsorbates
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if(ncell==0) start=I+1;
              else start=0;

              for(J=start;J<NumberOfCationMolecules[CurrentSystem];J++)
              {
                TypeMolB=Cations[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];

                    if(RigidJ)
                    {
                      index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                      index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=UNDEFINED_INT_VECTOR3;
                      comB=Cations[CurrentSystem][J].Atoms[j].Position;
                      index2=-1;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

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
                      if(ncell==0) ReplicaFactor=1.0;
                      else ReplicaFactor=0.5;

                      scalingB=Cations[CurrentSystem][J].Atoms[j].CFChargeScalingParameter;
                      ChargeB=scalingB*Cations[CurrentSystem][J].Atoms[j].Charge;

                      PotentialSecondDerivativeCoulombic(ChargeA,ChargeB,rr,&U,&f1,&f2);

                      *Energy+=ReplicaFactor*U;

                      StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
                      StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
                      StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

                      StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
                      StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
                      StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

                      StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
                      StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
                      StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=ReplicaFactor*(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=ReplicaFactor*temp1;
                        StrainDerivative->az-=ReplicaFactor*temp2;
                        StrainDerivative->bx-=ReplicaFactor*temp1;
                        StrainDerivative->by-=ReplicaFactor*(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=ReplicaFactor*temp3;
                        StrainDerivative->cx-=ReplicaFactor*temp2;
                        StrainDerivative->cy-=ReplicaFactor*temp3;
                        StrainDerivative->cz-=ReplicaFactor*(posA.z-comA.z)*f1*dr.z;
                      }


                      if(RigidJ)
                      {
                        temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                        temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                        temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                        StrainDerivative->ax+=ReplicaFactor*(posB.x-comB.x)*f1*dr.x;
                        StrainDerivative->ay+=ReplicaFactor*temp1;
                        StrainDerivative->az+=ReplicaFactor*temp2;
                        StrainDerivative->bx+=ReplicaFactor*temp1;
                        StrainDerivative->by+=ReplicaFactor*(posB.y-comB.y)*f1*dr.y;
                        StrainDerivative->bz+=ReplicaFactor*temp3;
                        StrainDerivative->cx+=ReplicaFactor*temp2;
                        StrainDerivative->cy+=ReplicaFactor*temp3;
                        StrainDerivative->cz+=ReplicaFactor*(posB.z-comB.z)*f1*dr.z;
                      }

                      // add contribution to the first derivatives
                      if(ComputeGradient)
                      {
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2.x>=0) Gradient[index_i2.x]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                          if(index_i2.y>=0) Gradient[index_i2.y]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                          if(index_i2.z>=0) Gradient[index_i2.z]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                        }

                        GradientStrain(Gradient,ReplicaFactor*f1,dr);

                        if(ncell==0)
                        {
                          if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                          if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                          if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                          if(RigidJ)
                          {
                            GradientStrainJ(Gradient,f1,dr,posB,comB);

                            // add contribution to the first derivatives
                            if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                            if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                            if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                          }
                        }
                      }

                      if(ComputeHessian)
                      {
                        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);
                        HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,ReplicaFactor*f1,ReplicaFactor*f2,dr,posA,comA,posB,comB,
                                                   RigidI,RigidJ);

                        HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,
                                  ReplicaFactor*f1,ReplicaFactor*f2,posA,comA,posB,comB,dr);

                        HessianAtomicStrainStrain(HessianMatrix,ReplicaFactor*f1,ReplicaFactor*f2,dr,posA,comA,posB,comB);
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
}


void CalculateBondConstraintExclusionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int m;
  REAL r,rr;
  POINT posA,posB;
  VECTOR dr;
  INT_VECTOR3 index_i,index_j;
  int ncell,k1,k2,k3;
  int typeA,typeB;
  REAL f1,f2,ReplicaFactor,energy;

  for(m=0;m<NumberOfDistanceConstraints[CurrentSystem];m++)
  {
    index_i=DistanceConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=DistanceConstraints[CurrentSystem][m][1]->HessianIndex;

    posA=DistanceConstraints[CurrentSystem][m][0]->Position;
    posB=DistanceConstraints[CurrentSystem][m][1]->Position;

    typeA=DistanceConstraints[CurrentSystem][m][0]->Type;
    typeB=DistanceConstraints[CurrentSystem][m][1]->Type;

    // TODO: FIX CHARGES
    //ChargeA=PseudoAtoms[typeA].Charge;
    //ChargeB=PseudoAtoms[typeB].Charge;

    ncell=0;
    for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
      for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
        for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
        {
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
            energy=f1=f2=0.0;
            if(ncell==0) ReplicaFactor=1.0;
            else ReplicaFactor=0.5;

            PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2,1.0);

            (*Energy)-=ReplicaFactor*energy;
            f1=-f1;
            f2=-f2;

            //if((index_i<0)&&(index_j<0)) continue;

            StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
            StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
            StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

            StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
            StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
            StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

            StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
            StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
            StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

            // add contribution to the first derivatives
            if(ComputeGradient)
            {
              if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
              if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
              if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;

              if(ncell==0)
              {
                if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;
              }
            }

            if(ComputeHessian)
            {
              HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);
            }
          }

/*
          if(rr<CutOffChargeChargeSquared)
          {
            energy=f1=f2=0.0;
            if(ncell==0) ReplicaFactor=1.0;
            else ReplicaFactor=0.5;

            switch(ChargeMethod)
            {
              case SHIFTED_COULOMB:
                *Energy-=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r-InverseCutOffChargeCharge);
                f1=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                f2=-3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                break;
              case TRUNCATED_COULOMB:
                *Energy-=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r);
                f1=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                f2=-3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                break;
              case EWALD:
              default:
                (*Energy)-=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*erfc(Alpha[CurrentSystem]*r)*
                                                    ChargeA*ChargeB/r;

                f1+=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                    (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                    (r*rr);

                f2-=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                    (((3.0*erfc(Alpha[CurrentSystem]*r)/(r*rr))+
                    (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                    (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));

                if(!OmitEwaldFourier)
                {
                  (*Energy)-=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*erf(Alpha[CurrentSystem]*r)*
                                   ChargeA*ChargeB/r;

                  f1+=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                      (erf(Alpha[CurrentSystem]*r)-2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                      (r*rr);

                  f2+=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                       (((-3.0*erf(Alpha[CurrentSystem]*r)/(r*rr))+
                       (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                       (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));
                }
                break;
            }

            if((index_i<0)&&(index_j<0)) continue;

            StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
            StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
            StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

            StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
            StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
            StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

            StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
            StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
            StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

            // add contribution to the first derivatives
            if(ComputeGradient)
            {
              if(index_i>=0)
              {
                Gradient[index_i]+=f1*dr.x;
                Gradient[index_i+1]+=f1*dr.y;
                Gradient[index_i+2]+=f1*dr.z;
              }

              if(ncell==0)
              {
                if(index_j>=0)
                {
                  Gradient[index_j]-=f1*dr.x;
                  Gradient[index_j+1]-=f1*dr.y;
                  Gradient[index_j+2]-=f1*dr.z;
                }
              }
            }

            if(ComputeHessian)
            {
              HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);
            }
          }
          ncell++;
*/
        }
  }
}

