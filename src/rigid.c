/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'rigid.c' is part of RASPA-2.0

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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "molecule.h"
#include "potentials.h"
#include "utils.h"
#include "integrate.h"
#include "simulation.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "internal_force.h"
#include "recrossing.h"
#include "rigid.h"
#include "ewald.h"
#include "cbmc.h"
#include "thermo_baro_stats.h"
#include "minimization.h"
#include "mc_moves.h"
#include "spacegroup.h"

int DistanceConstraintType;
int BendConstraintType;
int DihedralConstraintType;
int ImproperDihedralConstraintType;
int InversionBendConstraintType;
int OutOfPlaneConstraintType;

int ComputeRattleSteps;
REAL *NumberOfRattleCyclesStage1;
int *MaximumNumberOfRattleCyclesStage1;
REAL *NumberOfRattleCyclesStage2;
int *MaximumNumberOfRattleCyclesStage2;

// The space-fixed coordinate system is chosen in an inertial reference frame in which Newton's equations of motion are valid.
// The body-fixed coordinate system is an origin and three orthogonal axes (unit vectors) that are fixed to the body and rotate,
// tumble, spin and twist along with it.

// The unit quaternion q is introduced in order to generate a minimal, nonsingular, representation of the rotation matrix from a
// space-fixed denoted ‘s’ to a body-fixed coordinate system denoted ‘b’.
// The Quaternion Rotation Operator is L_q(v)=q* v q; this operator represents a rotation through an angle alpha about a vector q as its axis.
// We apply the quaternion rotation operator to a 3D vector v (a pure quaternion defined in the space-fixed frame), and express it as w in the
// body-fixed frame w_b=R(q) v_s
void BuildRotationMatrix(REAL_MATRIX3x3 *R,QUATERNION q)
{
  R->ax=2.0*(SQR(q.r)+SQR(q.i))-1.0; R->bx=2.0*(q.i*q.j+q.r*q.k);       R->cx=2.0*(q.i*q.k-q.r*q.j);
  R->ay=2.0*(q.i*q.j-q.r*q.k);       R->by=2.0*(SQR(q.j)+SQR(q.r))-1.0; R->cy=2.0*(q.j*q.k+q.r*q.i);
  R->az=2.0*(q.i*q.k+q.r*q.j);       R->bz=2.0*(q.j*q.k-q.r*q.i);       R->cz=2.0*(SQR(q.r)+SQR(q.k))-1.0;
}

// To do the opposite, namely apply the quaternion rotation operator to a 3D vector v (a pure quaternion defined in the body-fixed frame),
// and express it as w in the space-fixed frame, we have w_s=R(q)^-1 v_b=R(q)^T v_b
void BuildRotationMatrixInverse(REAL_MATRIX3x3 *R,QUATERNION q)
{
  R->ax=2.0*(SQR(q.r)+SQR(q.i))-1.0; R->bx=2.0*(q.i*q.j-q.r*q.k);       R->cx=2.0*(q.i*q.k+q.r*q.j);
  R->ay=2.0*(q.i*q.j+q.r*q.k);       R->by=2.0*(SQR(q.j)+SQR(q.r))-1.0; R->cy=2.0*(q.j*q.k-q.r*q.i);
  R->az=2.0*(q.i*q.k-q.r*q.j);       R->bz=2.0*(q.j*q.k+q.r*q.i);       R->cz=2.0*(SQR(q.r)+SQR(q.k))-1.0;
}

/***********************************************************************************************************
 * Name       | ComputeQuaternions                                                                         *
 * ------------------------------------------------------------------------------------------------------- *
 * Function   | Computes the quaternions for all rigid units in the system                                 *
 * Outline    | First it needs to be determined whether the rigid unit is 1D, 2D or 3D. Next, several      *
 *            | atoms need to be selected that defined the orientation. By comparing the current positions *
 *            | to the defined body-frame coordinates a rotation matrix can be computed. This rotation     *
 *            | matrix is converted to a quaternion.                                                       *
 * Parameters | -                                                                                          *
 ***********************************************************************************************************/

void ComputeQuaternions(void)
{
  int i,j,k,dimensionality;
  int A,B,C;
  REAL rotall,rotmin;
  REAL length;
  REAL_MATRIX3x3 aa,bb;
  int atom1,atom2,atom3;
  REAL dettest,det,rsq;

  atom1=atom2=atom3=0;
  for(j=0;j<NumberOfComponents;j++)
  {
    for(k=0;k<Components[j].NumberOfGroups;k++)
    {
      rotall=Components[j].Groups[k].InertiaVector.x+Components[j].Groups[k].InertiaVector.y+Components[j].Groups[k].InertiaVector.z;
      if(rotall<1e-5) rotall=1.0;

      rotmin=MIN3(Components[j].Groups[k].InertiaVector.x,Components[j].Groups[k].InertiaVector.y,Components[j].Groups[k].InertiaVector.z)/rotall;

      Components[j].Groups[k].rot_min=rotmin;

      dimensionality=0;
      if(Components[j].Groups[k].InertiaVector.x/rotall<1.0e-5) dimensionality++;
      if(Components[j].Groups[k].InertiaVector.y/rotall<1.0e-5) dimensionality++;
      if(Components[j].Groups[k].InertiaVector.z/rotall<1.0e-5) dimensionality++;

      switch(dimensionality)
      {
        // non-linear molecule
        case 0:
          atom1=0;
          atom2=0;
          atom3=0;
          do
          {
            atom2++;
            atom3=atom2;
            do
            {
              atom3++;
              A=Components[j].Groups[k].Atoms[atom1];
              B=Components[j].Groups[k].Atoms[atom2];
              C=Components[j].Groups[k].Atoms[atom3];

              // first row vector is A-B in body-fixed coordinates
              aa.ax=Components[j].Positions[A].x-Components[j].Positions[B].x;
              aa.bx=Components[j].Positions[A].y-Components[j].Positions[B].y;
              aa.cx=Components[j].Positions[A].z-Components[j].Positions[B].z;

              // second row vector is A-C in body-fixed coordinates
              aa.ay=Components[j].Positions[A].x-Components[j].Positions[C].x;
              aa.by=Components[j].Positions[A].y-Components[j].Positions[C].y;
              aa.cy=Components[j].Positions[A].z-Components[j].Positions[C].z;

              // third row vector is perpendicular to (A-B) and (A-C), i.e. the crossproduct
              aa.az=aa.bx*aa.cy-aa.cx*aa.by;
              aa.bz=aa.cx*aa.ay-aa.ax*aa.cy;
              aa.cz=aa.ax*aa.by-aa.bx*aa.ay;

              // invert matrix
              Invert3x3Matrix(&aa,&bb,&det);

              // check on size of determinant - to see if the 3 sites are
              // too close to being linear for safety.
              dettest=0.01;
            }
            while((fabs(det)<dettest)&&(atom3<Components[j].NumberOfAtoms));
          }
          while((fabs(det)<dettest)&&(atom2<Components[j].NumberOfAtoms-1));

          Components[j].Groups[k].InverseOriginalCoordinateSystem=bb;
          Components[j].Groups[k].Type=NONLINEAR_MOLECULE;
          break;
        // linear molecule
        case 1:
          atom1=0;
          atom2=1;
          atom3=0;
          A=Components[j].Groups[k].Atoms[atom1];
          B=Components[j].Groups[k].Atoms[atom2];
          C=Components[j].Groups[k].Atoms[atom3];

          // first row vector is A-B in body-fixed coordinates
          aa.ax=Components[j].Positions[A].x-Components[j].Positions[B].x;
          aa.bx=Components[j].Positions[A].y-Components[j].Positions[B].y;
          aa.cx=Components[j].Positions[A].z-Components[j].Positions[B].z;

          // length of the first vector
          length=sqrt(SQR(aa.ax)+SQR(aa.bx)+SQR(aa.cx));

          if(fabs(aa.cx/length)>0.5)
          {
            // set vector perpendicular in the (y,z)-plane
            rsq=sqrt(SQR(aa.bx)+SQR(aa.cx));
            aa.ay=0.0;
            aa.by=aa.cx/rsq;
            aa.cy=-aa.bx/rsq;
          }
          else if(fabs(aa.bx/length)>0.5)
          {
            // set vector perpendicular in the (x,y)-plane
            rsq=sqrt(SQR(aa.bx)+SQR(aa.ax));
            aa.ay=-aa.bx/rsq;
            aa.by=aa.ax/rsq;
            aa.cy=0.0;
          }
          else if(fabs(aa.ax/length)>0.5)
          {
            // set vector perpendicular in the (x,z)-plane
            rsq=sqrt(SQR(aa.ax)+SQR(aa.cx));
            aa.ay=-aa.cx/rsq;
            aa.by=0.0;
            aa.cy=aa.ax/rsq;
          }

          // third row vector is set perpendicular to the previous two vectors using the crossproduct
          aa.az=aa.bx*aa.cy-aa.cx*aa.by;
          aa.bz=aa.cx*aa.ay-aa.ax*aa.cy;
          aa.cz=aa.ax*aa.by-aa.bx*aa.ay;

          // invert matrix
          Invert3x3Matrix(&aa,&bb,&det);

          // check for singularity
          if(fabs(det)<1e-5)
          {
            fprintf(stderr, "failed to find principal axis system in routine 'void ComputeQuaternions(void)' in file 'rigid.c')\n");
            exit(0);
          }

          Components[j].Groups[k].InverseOriginalCoordinateSystem=bb;
          Components[j].Groups[k].Type=LINEAR_MOLECULE;
          break;
        // point particle
        default:
          atom1=atom2=atom3=0;
          bb.ax=bb.bx=bb.cx=0.0;
          bb.ay=bb.by=bb.cy=0.0;
          bb.az=bb.bz=bb.cz=0.0;
          Components[j].Groups[k].InverseOriginalCoordinateSystem=bb;
          Components[j].Groups[k].Type=POINT_PARTICLE;
          break;
      }

      Components[j].Groups[k].RotationalDegreesOfFreedom=3;
      rotall=1.0/MAX2((REAL)1.0e-5,Components[j].Groups[k].InertiaVector.x+Components[j].Groups[k].InertiaVector.y+Components[j].Groups[k].InertiaVector.z);

      if(rotall*Components[j].Groups[k].InertiaVector.x<1.0e-5)
        Components[j].Groups[k].RotationalDegreesOfFreedom--;

      if(rotall*Components[j].Groups[k].InertiaVector.y<1.0e-5)
        Components[j].Groups[k].RotationalDegreesOfFreedom--;

      if(rotall*Components[j].Groups[k].InertiaVector.z<1.0e-5)
        Components[j].Groups[k].RotationalDegreesOfFreedom--;

      A=Components[j].Groups[k].Atoms[atom1];
      B=Components[j].Groups[k].Atoms[atom2];
      C=Components[j].Groups[k].Atoms[atom3];

      // store the 3 atoms defining the orientation
      Components[j].Groups[k].orientation.A=A;
      Components[j].Groups[k].orientation.B=B;
      Components[j].Groups[k].orientation.C=C;
    }
  }

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      ComputeQuaternionAdsorbate(i);
      ComputeEulerAxisFromQuaternionAdsorbate(i);
    }

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      ComputeQuaternionCation(i);
      ComputeEulerAxisFromQuaternionCation(i);
    }
  }
}


/*********************************************************************************************************
 * Name       | RotationalMatrixToQuaternion                                                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Transform a rotation matrix to a quaternion                                              *
 * Parameters | The 3x3 rotation matrix R                                                                *
 * Note       | The algorithm returns a single quaternion q, but -q represents the same rotation         *
 * Ref.       | Page 305 "3-D Computer Graphics A Mathematical Introduction with OpenGL", Samuel R. Buss *
 *********************************************************************************************************/

QUATERNION RotationalMatrixToQuaternion(REAL_MATRIX3x3 R)
{
  int Switch;
  REAL trace,max;
  QUATERNION q;

  trace=R.ax+R.by+R.cz;

  Switch=0;
  max=trace;
  if(R.ax>max) {max=R.ax; Switch=1;}
  if(R.by>max) {max=R.by; Switch=2;}
  if(R.cz>max) {max=R.cz; Switch=3;}

  switch(Switch)
  {
    case 0:
      q.r=0.5*sqrt(trace+1.0);
      q.i=(R.cy-R.bz)/(4.0*q.r);
      q.j=(R.az-R.cx)/(4.0*q.r);
      q.k=(R.bx-R.ay)/(4.0*q.r);
      break;
    case 1:
      q.i=0.5*sqrt(2.0*R.ax-trace+1.0);
      q.r=(R.cy-R.bz)/(4.0*q.i);
      q.j=(R.bx+R.ay)/(4.0*q.i);
      q.k=(R.az+R.cx)/(4.0*q.i);
      break;
    case 2:
      q.j=0.5*sqrt(2.0*R.by-trace+1.0);
      q.r=(R.az-R.cx)/(4.0*q.j);
      q.i=(R.bx+R.ay)/(4.0*q.j);
      q.k=(R.cy+R.bz)/(4.0*q.j);
      break;
    case 3:
      q.k=0.5*sqrt(2.0*R.cz-trace+1.0);
      q.r=(R.bx-R.ay)/(4.0*q.k);
      q.i=(R.az+R.cx)/(4.0*q.k);
      q.j=(R.cy+R.bz)/(4.0*q.k);
      break;
  }

  return q;
}

/*********************************************************************************************************
 * Name       | ComputeQuaternionAdsorbate                                                               *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the quaternions for the adsorbates from their atomic space-fixed positions      *
 * Parameters | The adsorbate identifier 'm'                                                             *
 * Used in    | 'void ComputeQuaternions(void)'                                                          *
 *********************************************************************************************************/

void ComputeQuaternionAdsorbate(int m)
{
/*
  int k;
  int Type,A,B,C;
  REAL_MATRIX3x3 rot,aa,bb;
  REAL rotmin,rsq;
  VECTOR dr;

  Type=Adsorbates[CurrentSystem][m].Type;
  for(k=0;k<Components[Type].NumberOfGroups;k++)
  {
    if(Components[Type].Groups[k].Rigid)
    {
      rotmin=Components[Type].Groups[k].rot_min;
      A=Components[Type].Groups[k].orientation.A;
      B=Components[Type].Groups[k].orientation.B;
      C=Components[Type].Groups[k].orientation.C;

      // construct the final coordinate system matrix
      aa.ax=aa.ay=aa.cy=aa.by=0.0;

      // first row vector is A-B in space-fixed coordinates
      dr.x=Adsorbates[CurrentSystem][m].Atoms[A].Position.x-Adsorbates[CurrentSystem][m].Atoms[B].Position.x;
      dr.y=Adsorbates[CurrentSystem][m].Atoms[A].Position.y-Adsorbates[CurrentSystem][m].Atoms[B].Position.y;
      dr.z=Adsorbates[CurrentSystem][m].Atoms[A].Position.z-Adsorbates[CurrentSystem][m].Atoms[B].Position.z;
      dr=ApplyBoundaryCondition(dr);
      aa.ax=dr.x;
      aa.bx=dr.y;
      aa.cx=dr.z;

      // second row vector is A-C in space-fixed coordinates
      if(rotmin>1.0e-5)
      {
        dr.x=Adsorbates[CurrentSystem][m].Atoms[A].Position.x-Adsorbates[CurrentSystem][m].Atoms[C].Position.x;
        dr.y=Adsorbates[CurrentSystem][m].Atoms[A].Position.y-Adsorbates[CurrentSystem][m].Atoms[C].Position.y;
        dr.z=Adsorbates[CurrentSystem][m].Atoms[A].Position.z-Adsorbates[CurrentSystem][m].Atoms[C].Position.z;
        dr=ApplyBoundaryCondition(dr);
        aa.ay=dr.x;
        aa.by=dr.y;
        aa.cy=dr.z;
      }
      else
      {
        rsq=sqrt(SQR(aa.ax)+SQR(aa.bx)+SQR(aa.cx));
        if(fabs(aa.cx/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.bx)+SQR(aa.cx));
          aa.ay=0.0;
          aa.by=aa.cx/rsq;
          aa.cy=-aa.bx/rsq;
        }
        else if(fabs(aa.bx/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.bx)+SQR(aa.ax));
          aa.ay=-aa.bx/rsq;
          aa.by=aa.ax/rsq;
          aa.cy=0.0;
        }
        else if(fabs(aa.ax/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.ax)+SQR(aa.cx));
          aa.ay=-aa.cx/rsq;
          aa.by=0.0;
          aa.cy=aa.ax/rsq;
        }
      }
      // third row vector is (A-B)x(A-C) in space-fixed coordinates
      aa.az=aa.bx*aa.cy-aa.cx*aa.by;
      aa.bz=aa.cx*aa.ay-aa.ax*aa.cy;
      aa.cz=aa.ax*aa.by-aa.bx*aa.ay;

      // get the inverse of the original coordinate system
      bb=Components[Type].Groups[k].InverseOriginalCoordinateSystem;

      // the desired rotation matrix can be found by multiplying the inverse of
      // the original coordinate system with the final coordinate system matrix.
      // M_rot=M_final.M_orig^-1
      rot.ax=bb.ax*aa.ax+bb.bx*aa.ay+bb.cx*aa.az;
      rot.ay=bb.ay*aa.ax+bb.by*aa.ay+bb.cy*aa.az;
      rot.az=bb.az*aa.ax+bb.bz*aa.ay+bb.cz*aa.az;

      rot.bx=bb.ax*aa.bx+bb.bx*aa.by+bb.cx*aa.bz;
      rot.by=bb.ay*aa.bx+bb.by*aa.by+bb.cy*aa.bz;
      rot.bz=bb.az*aa.bx+bb.bz*aa.by+bb.cz*aa.bz;

      rot.cx=bb.ax*aa.cx+bb.bx*aa.cy+bb.cx*aa.cz;
      rot.cy=bb.ay*aa.cx+bb.by*aa.cy+bb.cy*aa.cz;
      rot.cz=bb.az*aa.cx+bb.bz*aa.cy+bb.cz*aa.cz;

      // transform the rotation matrix to a quaternion
      Adsorbates[CurrentSystem][m].Groups[k].Quaternion=RotationalMatrixToQuaternion(rot);
    }
  }
*/
  int k;
  int Type,A,B,C;
  REAL_MATRIX3x3 rot,aa,bb;
  REAL rotmin,rsq;
  REAL aq,bq,cq,dq,eq,fq,gq,hq;
  QUATERNION q;

  Type=Adsorbates[CurrentSystem][m].Type;
  for(k=0;k<Components[Type].NumberOfGroups;k++)
  {
    if(Components[Type].Groups[k].Rigid)
    {
      rotmin=Components[Type].Groups[k].rot_min;
      A=Components[Type].Groups[k].orientation.A;
      B=Components[Type].Groups[k].orientation.B;
      C=Components[Type].Groups[k].orientation.C;

      aa.ax=aa.ay=aa.cy=aa.by=0.0;
      aa.ax=Adsorbates[CurrentSystem][m].Atoms[A].Position.x-Adsorbates[CurrentSystem][m].Atoms[B].Position.x;
      aa.bx=Adsorbates[CurrentSystem][m].Atoms[A].Position.y-Adsorbates[CurrentSystem][m].Atoms[B].Position.y;
      aa.cx=Adsorbates[CurrentSystem][m].Atoms[A].Position.z-Adsorbates[CurrentSystem][m].Atoms[B].Position.z;

      if(rotmin>1.0e-5)
      {
        aa.ay=Adsorbates[CurrentSystem][m].Atoms[A].Position.x-Adsorbates[CurrentSystem][m].Atoms[C].Position.x;
        aa.by=Adsorbates[CurrentSystem][m].Atoms[A].Position.y-Adsorbates[CurrentSystem][m].Atoms[C].Position.y;
        aa.cy=Adsorbates[CurrentSystem][m].Atoms[A].Position.z-Adsorbates[CurrentSystem][m].Atoms[C].Position.z;
      }
      else
      {
        rsq=sqrt(SQR(aa.ax)+SQR(aa.bx)+SQR(aa.cx));
        if(fabs(aa.cx/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.bx)+SQR(aa.cx));
          aa.ay=0.0;
          aa.by=aa.cx/rsq;
          aa.cy=-aa.bx/rsq;
        }
        else if(fabs(aa.bx/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.bx)+SQR(aa.ax));
          aa.ay=-aa.bx/rsq;
          aa.by=aa.ax/rsq;
          aa.cy=0.0;
        }
        else if(fabs(aa.ax/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.ax)+SQR(aa.cx));
          aa.ay=-aa.cx/rsq;
          aa.by=0.0;
          aa.cy=aa.ax/rsq;
        }
      }
      aa.az=aa.bx*aa.cy-aa.cx*aa.by;
      aa.bz=aa.cx*aa.ay-aa.ax*aa.cy;
      aa.cz=aa.ax*aa.by-aa.bx*aa.ay;

          // group rotational matrix
      bb=Components[Type].Groups[k].InverseOriginalCoordinateSystem;
      rot.ax=bb.ax*aa.ax+bb.bx*aa.ay+bb.cx*aa.az;
      rot.ay=bb.ay*aa.ax+bb.by*aa.ay+bb.cy*aa.az;
      rot.az=bb.az*aa.ax+bb.bz*aa.ay+bb.cz*aa.az;

      rot.bx=bb.ax*aa.bx+bb.bx*aa.by+bb.cx*aa.bz;
      rot.by=bb.ay*aa.bx+bb.by*aa.by+bb.cy*aa.bz;
      rot.bz=bb.az*aa.bx+bb.bz*aa.by+bb.cz*aa.bz;

      rot.cx=bb.ax*aa.cx+bb.bx*aa.cy+bb.cx*aa.cz;
      rot.cy=bb.ay*aa.cx+bb.by*aa.cy+bb.cy*aa.cz;
      rot.cz=bb.az*aa.cx+bb.bz*aa.cy+bb.cz*aa.cz;

      // determine quaternions from rotational matrix
      aq=rot.ax+rot.by;
      bq=rot.ay-rot.bx;
      cq=rot.bz-rot.cy;
      dq=rot.ay+rot.bx;
      eq=rot.az+rot.cx;
      fq=rot.bz+rot.cy;
      gq=rot.az-rot.cx;
      hq=rot.ax-rot.by;

      q.r=0.5*sqrt(aq+sqrt(aq*aq+bq*bq));

      if(q.r>1.0e-4)
      {
        q.i=-0.25*cq/q.r;
        q.j=0.25*gq/q.r;
        q.k=-0.25*bq/q.r;
      }
      else
      {
        q.i=0.5*sqrt(hq+sqrt(hq*hq+dq*dq));
        if(q.i>1.0e-4)
        {
          q.j=0.25*dq/q.i;
          q.k=0.25*eq/q.i;
        }
        else
        {
          q.j=0.5*sqrt(-hq+sqrt(hq*hq+dq*dq));
          if(q.j>1.0e-4)
            q.k=0.25*fq/q.j;
          else
            q.k=1.0;
        }
      }

      // normalise quaternions
      rsq=sqrt(SQR(q.r)+SQR(q.i)+SQR(q.j)+SQR(q.k));
      q.r/=rsq; q.i/=rsq;
      q.j/=rsq; q.k/=rsq;
      Adsorbates[CurrentSystem][m].Groups[k].Quaternion=q;
    }
  }

}


/*********************************************************************************************************
 * Name       | ComputeQuaternionCation                                                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the quaternions for the cations from their atomic space-fixed positions         *
 * Parameters | The cation identifier 'm'                                                                *
 * Used in    | 'void ComputeQuaternions(void)'                                                          *
 *********************************************************************************************************/

void ComputeQuaternionCation(int m)
{
/*
  int k;
  int Type,A,B,C;
  REAL_MATRIX3x3 rot,aa,bb;
  REAL rotmin,rsq;
  VECTOR dr;

  Type=Cations[CurrentSystem][m].Type;
  for(k=0;k<Components[Type].NumberOfGroups;k++)
  {
    if(Components[Type].Groups[k].Rigid)
    {
      rotmin=Components[Type].Groups[k].rot_min;
      A=Components[Type].Groups[k].orientation.A;
      B=Components[Type].Groups[k].orientation.B;
      C=Components[Type].Groups[k].orientation.C;

      // construct the final coordinate system matrix
      aa.ax=aa.ay=aa.cy=aa.by=0.0;

      // first row vector is A-B in space-fixed coordinates
      dr.x=Cations[CurrentSystem][m].Atoms[A].Position.x-Cations[CurrentSystem][m].Atoms[B].Position.x;
      dr.y=Cations[CurrentSystem][m].Atoms[A].Position.y-Cations[CurrentSystem][m].Atoms[B].Position.y;
      dr.z=Cations[CurrentSystem][m].Atoms[A].Position.z-Cations[CurrentSystem][m].Atoms[B].Position.z;
      dr=ApplyBoundaryCondition(dr);
      aa.ax=dr.x;
      aa.bx=dr.y;
      aa.cx=dr.z;

      // second row vector is A-C in space-fixed coordinates
      if(rotmin>1.0e-5)
      {
        dr.x=Cations[CurrentSystem][m].Atoms[A].Position.x-Cations[CurrentSystem][m].Atoms[C].Position.x;
        dr.y=Cations[CurrentSystem][m].Atoms[A].Position.y-Cations[CurrentSystem][m].Atoms[C].Position.y;
        dr.z=Cations[CurrentSystem][m].Atoms[A].Position.z-Cations[CurrentSystem][m].Atoms[C].Position.z;
        dr=ApplyBoundaryCondition(dr);
        aa.ay=dr.x;
        aa.by=dr.y;
        aa.cy=dr.z;
      }
      else
      {
        rsq=sqrt(SQR(aa.ax)+SQR(aa.bx)+SQR(aa.cx));
        if(fabs(aa.cx/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.bx)+SQR(aa.cx));
          aa.ay=0.0;
          aa.by=aa.cx/rsq;
          aa.cy=-aa.bx/rsq;
        }
        else if(fabs(aa.bx/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.bx)+SQR(aa.ax));
          aa.ay=-aa.bx/rsq;
          aa.by=aa.ax/rsq;
          aa.cy=0.0;
        }
        else if(fabs(aa.ax/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.ax)+SQR(aa.cx));
          aa.ay=-aa.cx/rsq;
          aa.by=0.0;
          aa.cy=aa.ax/rsq;
        }
      }
      // third row vector is (A-B)x(A-C) in space-fixed coordinates
      aa.az=aa.bx*aa.cy-aa.cx*aa.by;
      aa.bz=aa.cx*aa.ay-aa.ax*aa.cy;
      aa.cz=aa.ax*aa.by-aa.bx*aa.ay;

      // get the inverse of the original coordinate system
      bb=Components[Type].Groups[k].InverseOriginalCoordinateSystem;

      // the desired rotation matrix can be found by multiplying the inverse of
      // the original coordinate system with the final coordinate system matrix.
      // M_rot=M_final.M_orig^-1
      rot.ax=bb.ax*aa.ax+bb.bx*aa.ay+bb.cx*aa.az;
      rot.ay=bb.ay*aa.ax+bb.by*aa.ay+bb.cy*aa.az;
      rot.az=bb.az*aa.ax+bb.bz*aa.ay+bb.cz*aa.az;

      rot.bx=bb.ax*aa.bx+bb.bx*aa.by+bb.cx*aa.bz;
      rot.by=bb.ay*aa.bx+bb.by*aa.by+bb.cy*aa.bz;
      rot.bz=bb.az*aa.bx+bb.bz*aa.by+bb.cz*aa.bz;

      rot.cx=bb.ax*aa.cx+bb.bx*aa.cy+bb.cx*aa.cz;
      rot.cy=bb.ay*aa.cx+bb.by*aa.cy+bb.cy*aa.cz;
      rot.cz=bb.az*aa.cx+bb.bz*aa.cy+bb.cz*aa.cz;

      // transform the rotation matrix to a quaternion
      Cations[CurrentSystem][m].Groups[k].Quaternion=RotationalMatrixToQuaternion(rot);
    }
  }
*/
  int k;
  int Type,A,B,C;
  REAL_MATRIX3x3 rot,aa,bb;
  REAL rotmin,rsq;
  REAL aq,bq,cq,dq,eq,fq,gq,hq;
  QUATERNION q;

  Type=Cations[CurrentSystem][m].Type;
  for(k=0;k<Components[Type].NumberOfGroups;k++)
  {
    if(Components[Type].Groups[k].Rigid)
    {
      rotmin=Components[Type].Groups[k].rot_min;
      A=Components[Type].Groups[k].orientation.A;
      B=Components[Type].Groups[k].orientation.B;
      C=Components[Type].Groups[k].orientation.C;

      aa.ax=aa.ay=aa.cy=aa.by=0.0;
      aa.ax=Cations[CurrentSystem][m].Atoms[A].Position.x-Cations[CurrentSystem][m].Atoms[B].Position.x;
      aa.bx=Cations[CurrentSystem][m].Atoms[A].Position.y-Cations[CurrentSystem][m].Atoms[B].Position.y;
      aa.cx=Cations[CurrentSystem][m].Atoms[A].Position.z-Cations[CurrentSystem][m].Atoms[B].Position.z;
      if(rotmin>1.0e-5)
      {
        aa.ay=Cations[CurrentSystem][m].Atoms[A].Position.x-Cations[CurrentSystem][m].Atoms[C].Position.x;
        aa.by=Cations[CurrentSystem][m].Atoms[A].Position.y-Cations[CurrentSystem][m].Atoms[C].Position.y;
        aa.cy=Cations[CurrentSystem][m].Atoms[A].Position.z-Cations[CurrentSystem][m].Atoms[C].Position.z;
      }
      else
      {
        rsq=sqrt(SQR(aa.ax)+SQR(aa.bx)+SQR(aa.cx));
        if(fabs(aa.cx/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.bx)+SQR(aa.cx));
          aa.ay=0.0;
          aa.by=aa.cx/rsq;
          aa.cy=-aa.bx/rsq;
        }
        else if(fabs(aa.bx/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.bx)+SQR(aa.ax));
          aa.ay=-aa.bx/rsq;
          aa.by=aa.ax/rsq;
          aa.cy=0.0;
        }
        else if(fabs(aa.ax/rsq)>0.5)
        {
          rsq=sqrt(SQR(aa.ax)+SQR(aa.cx));
          aa.ay=-aa.cx/rsq;
          aa.by=0.0;
          aa.cy=aa.ax/rsq;
        }
      }
      aa.az=aa.bx*aa.cy-aa.cx*aa.by;
      aa.bz=aa.cx*aa.ay-aa.ax*aa.cy;
      aa.cz=aa.ax*aa.by-aa.bx*aa.ay;

          // group rotational matrix
      bb=Components[Type].Groups[k].InverseOriginalCoordinateSystem;
      rot.ax=bb.ax*aa.ax+bb.bx*aa.ay+bb.cx*aa.az;
      rot.ay=bb.ay*aa.ax+bb.by*aa.ay+bb.cy*aa.az;
      rot.az=bb.az*aa.ax+bb.bz*aa.ay+bb.cz*aa.az;

      rot.bx=bb.ax*aa.bx+bb.bx*aa.by+bb.cx*aa.bz;
      rot.by=bb.ay*aa.bx+bb.by*aa.by+bb.cy*aa.bz;
      rot.bz=bb.az*aa.bx+bb.bz*aa.by+bb.cz*aa.bz;

      rot.cx=bb.ax*aa.cx+bb.bx*aa.cy+bb.cx*aa.cz;
      rot.cy=bb.ay*aa.cx+bb.by*aa.cy+bb.cy*aa.cz;
      rot.cz=bb.az*aa.cx+bb.bz*aa.cy+bb.cz*aa.cz;

      // determine quaternions from rotational matrix
      aq=rot.ax+rot.by;
      bq=rot.ay-rot.bx;
      cq=rot.bz-rot.cy;
      dq=rot.ay+rot.bx;
      eq=rot.az+rot.cx;
      fq=rot.bz+rot.cy;
      gq=rot.az-rot.cx;
      hq=rot.ax-rot.by;

      q.r=0.5*sqrt(aq+sqrt(aq*aq+bq*bq));

      if(q.r>1.0e-4)
      {
        q.i=-0.25*cq/q.r;
        q.j=0.25*gq/q.r;
        q.k=-0.25*bq/q.r;
      }
      else
      {
        q.i=0.5*sqrt(hq+sqrt(hq*hq+dq*dq));
        if(q.i>1.0e-4)
        {
          q.j=0.25*dq/q.i;
          q.k=0.25*eq/q.i;
        }
        else
        {
          q.j=0.5*sqrt(-hq+sqrt(hq*hq+dq*dq));
          if(q.j>1.0e-4)
            q.k=0.25*fq/q.j;
          else
            q.k=1.0;
        }
      }

      // normalise quaternions
      rsq=sqrt(SQR(q.r)+SQR(q.i)+SQR(q.j)+SQR(q.k));
      q.r/=rsq; q.i/=rsq;
      q.j/=rsq; q.k/=rsq;
      Cations[CurrentSystem][m].Groups[k].Quaternion=q;
    }
  }
}

void ComputeEulerAxisFromQuaternionAdsorbate(int m)
{
  int l,MolType;
  VECTOR EulerAxis;
  REAL_MATRIX3x3 RotationMatrix;
  REAL theta;

  MolType=Adsorbates[CurrentSystem][m].Type;
  for(l=0;l<Components[MolType].NumberOfGroups;l++)
  {
    if(Components[MolType].Groups[l].Rigid) // rigid unit
    {
      BuildRotationMatrix(&RotationMatrix,Adsorbates[CurrentSystem][m].Groups[l].Quaternion);

      // use Euler-axis as generalized coordinates for rotation
      theta=acos(0.5*(RotationMatrix.ax+RotationMatrix.by+RotationMatrix.cz-1.0));
      if(fabs(theta)<1e-8)
      {
        EulerAxis.x=EulerAxis.y=EulerAxis.z=0.0;
      }
      else
      {
        EulerAxis.x=theta*(RotationMatrix.cy-RotationMatrix.bz)/(2.0*sin(theta));
        EulerAxis.y=theta*(RotationMatrix.az-RotationMatrix.cx)/(2.0*sin(theta));
        EulerAxis.z=theta*(RotationMatrix.bx-RotationMatrix.ay)/(2.0*sin(theta));
      }
      Adsorbates[CurrentSystem][m].Groups[l].EulerAxis=EulerAxis;
    }
  }
}

void ComputeEulerAxisFromQuaternionCation(int m)
{
  int l,MolType;
  VECTOR EulerAxis;
  REAL_MATRIX3x3 RotationMatrix;
  REAL theta;

  MolType=Cations[CurrentSystem][m].Type;
  for(l=0;l<Components[MolType].NumberOfGroups;l++)
  {
    if(Components[MolType].Groups[l].Rigid) // rigid unit
    {
      BuildRotationMatrix(&RotationMatrix,Cations[CurrentSystem][m].Groups[l].Quaternion);

      // use Euler-axis as generalized coordinates for rotation
      theta=acos(0.5*(RotationMatrix.ax+RotationMatrix.by+RotationMatrix.cz-1.0));
      if(fabs(theta)<1e-8)
      {
        EulerAxis.x=EulerAxis.y=EulerAxis.z=0.0;
      }
      else
      {
        EulerAxis.x=theta*(RotationMatrix.cy-RotationMatrix.bz)/(2.0*sin(theta));
        EulerAxis.y=theta*(RotationMatrix.az-RotationMatrix.cx)/(2.0*sin(theta));
        EulerAxis.z=theta*(RotationMatrix.bx-RotationMatrix.ay)/(2.0*sin(theta));
      }
      Cations[CurrentSystem][m].Groups[l].EulerAxis=EulerAxis;
    }
  }
}


// computes the center-of-mass positions and center-of-mass velocities for rigid units in the molecules
void ComputeMolecularVelocitiesAndPositions(void)
{
  int i,k,l;
  int A,Type;
  REAL Mass,TotalMass;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      TotalMass=0.0;

      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.x=0.0;
      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.y=0.0;
      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.z=0.0;

      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x=0.0;
      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y=0.0;
      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z=0.0;

      for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
      {
        A=Components[Type].Groups[l].Atoms[k];
        Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
        TotalMass+=Mass;

        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.x+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Position.x;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.y+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Position.y;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.z+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Position.z;

        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z;
      }

      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.x/=TotalMass;
      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.y/=TotalMass;
      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.z/=TotalMass;

      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x/=TotalMass;
      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y/=TotalMass;
      Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z/=TotalMass;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      TotalMass=0.0;

      Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition.x=0.0;
      Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition.y=0.0;
      Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition.z=0.0;

      Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x=0.0;
      Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y=0.0;
      Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z=0.0;

      for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
      {
        A=Components[Type].Groups[l].Atoms[k];
        Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[A].Type].Mass;
        TotalMass+=Mass;

        Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition.x+=Mass*Cations[CurrentSystem][i].Atoms[A].Position.x;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition.y+=Mass*Cations[CurrentSystem][i].Atoms[A].Position.y;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition.z+=Mass*Cations[CurrentSystem][i].Atoms[A].Position.z;

        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.x;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.y;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.z;
      }

      Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition.x/=TotalMass;
      Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition.y/=TotalMass;
      Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition.z/=TotalMass;

      Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x/=TotalMass;
      Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y/=TotalMass;
      Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z/=TotalMass;
    }
  }
}

// computes the angular velocities of the rigid units from the atomic Cartesian velocities
void ComputeAngularVelocities(void)
{
  int i,l;
  int Type;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
      Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity=AtomicVelocityToAngularVelocityAdsorbates(i,l);
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
      Cations[CurrentSystem][i].Groups[l].AngularVelocity=AtomicVelocityToAngularVelocityCations(i,l);
  }
}

VECTOR AtomicVelocityToAngularVelocityAdsorbates(int i,int g)
{
  int j;
  int Type,A;
  VECTOR omega,inv,vel,dr,AngularVelocity;
  REAL Mass;
  REAL_MATRIX3x3 M;

  Adsorbates[CurrentSystem][i].Groups[g].AngularVelocity.x=0.0;
  Adsorbates[CurrentSystem][i].Groups[g].AngularVelocity.y=0.0;
  Adsorbates[CurrentSystem][i].Groups[g].AngularVelocity.z=0.0;

  // angular momentum
  omega.x=0.0;
  omega.y=0.0;
  omega.z=0.0;

  Type=Adsorbates[CurrentSystem][i].Type;
  if(!Components[Type].Groups[g].Rigid) return omega;

  for(j=0;j<Components[Type].Groups[g].NumberOfGroupAtoms;j++)
  {
    A=Components[Type].Groups[g].Atoms[j];

    dr.x=Adsorbates[CurrentSystem][i].Atoms[A].Position.x-Adsorbates[CurrentSystem][i].Groups[g].CenterOfMassPosition.x;
    dr.y=Adsorbates[CurrentSystem][i].Atoms[A].Position.y-Adsorbates[CurrentSystem][i].Groups[g].CenterOfMassPosition.y;
    dr.z=Adsorbates[CurrentSystem][i].Atoms[A].Position.z-Adsorbates[CurrentSystem][i].Groups[g].CenterOfMassPosition.z;

    Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
    vel=Adsorbates[CurrentSystem][i].Atoms[A].Velocity;
    omega.x+=Mass*(dr.y*vel.z-dr.z*vel.y);
    omega.y+=Mass*(dr.z*vel.x-dr.x*vel.z);
    omega.z+=Mass*(dr.x*vel.y-dr.y*vel.x);
  }

  // build rotation matrix
  BuildRotationMatrix(&M,Adsorbates[CurrentSystem][i].Groups[g].Quaternion);
  inv=Components[Type].Groups[g].InverseInertiaVector;
  AngularVelocity.x=(M.ax*omega.x+M.bx*omega.y+M.cx*omega.z)*inv.x;
  AngularVelocity.y=(M.ay*omega.x+M.by*omega.y+M.cy*omega.z)*inv.y;
  AngularVelocity.z=(M.az*omega.x+M.bz*omega.y+M.cz*omega.z)*inv.z;
  return AngularVelocity;
}

VECTOR AtomicVelocityToAngularVelocityCations(int i,int g)
{
  int j;
  int Type,A;
  VECTOR omega,inv,vel,dr,AngularVelocity;
  REAL Mass;
  REAL_MATRIX3x3 M;

  Type=Cations[CurrentSystem][i].Type;
  Cations[CurrentSystem][i].Groups[g].AngularVelocity.x=0.0;
  Cations[CurrentSystem][i].Groups[g].AngularVelocity.y=0.0;
  Cations[CurrentSystem][i].Groups[g].AngularVelocity.z=0.0;

  // angular momentum
  omega.x=0.0;
  omega.y=0.0;
  omega.z=0.0;
  for(j=0;j<Components[Type].Groups[g].NumberOfGroupAtoms;j++)
  {
    A=Components[Type].Groups[g].Atoms[j];

    dr.x=Cations[CurrentSystem][i].Atoms[A].Position.x-Cations[CurrentSystem][i].Groups[g].CenterOfMassPosition.x;
    dr.y=Cations[CurrentSystem][i].Atoms[A].Position.y-Cations[CurrentSystem][i].Groups[g].CenterOfMassPosition.y;
    dr.z=Cations[CurrentSystem][i].Atoms[A].Position.z-Cations[CurrentSystem][i].Groups[g].CenterOfMassPosition.z;

    Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[A].Type].Mass;
    vel=Cations[CurrentSystem][i].Atoms[A].Velocity;
    omega.x+=Mass*(dr.y*vel.z-dr.z*vel.y);
    omega.y+=Mass*(dr.z*vel.x-dr.x*vel.z);
    omega.z+=Mass*(dr.x*vel.y-dr.y*vel.x);
  }

  // build rotation matrix
  BuildRotationMatrix(&M,Cations[CurrentSystem][i].Groups[g].Quaternion);
  inv=Components[Cations[CurrentSystem][i].Type].Groups[g].InverseInertiaVector;
  AngularVelocity.x=(M.ax*omega.x+M.bx*omega.y+M.cx*omega.z)*inv.x;
  AngularVelocity.y=(M.ay*omega.x+M.by*omega.y+M.cy*omega.z)*inv.y;
  AngularVelocity.z=(M.az*omega.x+M.bz*omega.y+M.cz*omega.z)*inv.z;
  return AngularVelocity;
}

void ConvertAngularVelocityToAtomicVelocity(void)
{
  int i,k,l;
  int A,Type;
  VECTOR pos,w,vel,omega;
  QUATERNION q;
  REAL_MATRIX3x3 M;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        q=Adsorbates[CurrentSystem][i].Groups[l].Quaternion;
        BuildRotationMatrixInverse(&M,q);
        omega=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity;

        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];

          // site position in body frame
          pos=Components[Type].Positions[A];
          w.x=omega.y*pos.z-omega.z*pos.y;
          w.y=omega.z*pos.x-omega.x*pos.z;
          w.z=omega.x*pos.y-omega.y*pos.x;

          // new atomic velocites in lab frame
          vel=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x=(M.ax*w.x+M.bx*w.y+M.cx*w.z)+vel.x;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y=(M.ay*w.x+M.by*w.y+M.cy*w.z)+vel.y;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z=(M.az*w.x+M.bz*w.y+M.cz*w.z)+vel.z;
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        q=Cations[CurrentSystem][i].Groups[l].Quaternion;
        BuildRotationMatrixInverse(&M,q);
        omega=Cations[CurrentSystem][i].Groups[l].AngularVelocity;

        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];

          // site position in body frame
          pos=Components[Type].Positions[A];
          w.x=omega.y*pos.z-omega.z*pos.y;
          w.y=omega.z*pos.x-omega.x*pos.z;
          w.z=omega.x*pos.y-omega.y*pos.x;

          // new atomic velocites in lab frame
          vel=Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
          Cations[CurrentSystem][i].Atoms[A].Velocity.x=(M.ax*w.x+M.bx*w.y+M.cx*w.z)+vel.x;
          Cations[CurrentSystem][i].Atoms[A].Velocity.y=(M.ay*w.x+M.by*w.y+M.cy*w.z)+vel.y;
          Cations[CurrentSystem][i].Atoms[A].Velocity.z=(M.az*w.x+M.bz*w.y+M.cz*w.z)+vel.z;
        }
      }
    }
  }
}



void ComputeQuaternionMomenta(void)
{
  int i,l;
  int Type;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
      Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum=AngularVelocityToQuaternionMomentumAdsorbates(i,l);
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
      Cations[CurrentSystem][i].Groups[l].QuaternionMomentum=AngularVelocityToQuaternionMomentumCations(i,l);
  }
}

QUATERNION AngularVelocityToQuaternionMomentumAdsorbates(int i,int g)
{
  QUATERNION p,q;
  VECTOR omega,I;

  omega=Adsorbates[CurrentSystem][i].Groups[g].AngularVelocity;
  I=Components[Adsorbates[CurrentSystem][i].Type].Groups[g].InertiaVector;
  q=Adsorbates[CurrentSystem][i].Groups[g].Quaternion;
  p.r=2.0*(-q.i*(I.x*omega.x)-q.j*(I.y*omega.y)-q.k*(I.z*omega.z));
  p.i=2.0*(q.r*(I.x*omega.x)-q.k*(I.y*omega.y)+q.j*(I.z*omega.z));
  p.j=2.0*(q.k*(I.x*omega.x)+q.r*(I.y*omega.y)-q.i*(I.z*omega.z));
  p.k=2.0*(-q.j*(I.x*omega.x)+q.i*(I.y*omega.y)+q.r*(I.z*omega.z));
  return p;
}

QUATERNION AngularVelocityToQuaternionMomentumCations(int i,int g)
{
  QUATERNION p,q;
  VECTOR omega,I;

  omega=Cations[CurrentSystem][i].Groups[g].AngularVelocity;
  I=Components[Cations[CurrentSystem][i].Type].Groups[g].InertiaVector;
  q=Cations[CurrentSystem][i].Groups[g].Quaternion;
  p.r=2.0*(-q.i*(I.x*omega.x)-q.j*(I.y*omega.y)-q.k*(I.z*omega.z));
  p.i=2.0*(q.r*(I.x*omega.x)-q.k*(I.y*omega.y)+q.j*(I.z*omega.z));
  p.j=2.0*(q.k*(I.x*omega.x)+q.r*(I.y*omega.y)-q.i*(I.z*omega.z));
  p.k=2.0*(-q.j*(I.x*omega.x)+q.i*(I.y*omega.y)+q.r*(I.z*omega.z));
  return p;
}


VECTOR QuaternionMomentumToAngularVelocityAdsorbates(int i,int g)
{
  QUATERNION p,q;
  VECTOR inv,omega;

  q=Adsorbates[CurrentSystem][i].Groups[g].Quaternion;
  p=Adsorbates[CurrentSystem][i].Groups[g].QuaternionMomentum;
  inv=Components[Adsorbates[CurrentSystem][i].Type].Groups[g].InverseInertiaVector;
  omega.x=0.5*(-q.i*p.r+q.r*p.i+q.k*p.j-q.j*p.k)*inv.x;
  omega.y=0.5*(-q.j*p.r-q.k*p.i+q.r*p.j+q.i*p.k)*inv.y;
  omega.z=0.5*(-q.k*p.r+q.j*p.i-q.i*p.j+q.r*p.k)*inv.z;
  return omega;
}

VECTOR QuaternionMomentumToAngularVelocityCations(int i,int g)
{
  QUATERNION p,q;
  VECTOR inv,omega;

  q=Cations[CurrentSystem][i].Groups[g].Quaternion;
  p=Cations[CurrentSystem][i].Groups[g].QuaternionMomentum;
  inv=Components[Cations[CurrentSystem][i].Type].Groups[g].InverseInertiaVector;
  omega.x=0.5*(-q.i*p.r+q.r*p.i+q.k*p.j-q.j*p.k)*inv.x;
  omega.y=0.5*(-q.j*p.r-q.k*p.i+q.r*p.j+q.i*p.k)*inv.y;
  omega.z=0.5*(-q.k*p.r+q.j*p.i-q.i*p.j+q.r*p.k)*inv.z;
  return omega;
}


VECTOR TotalAngularMomentum(void)
{
  int i,j;
  int Type;
  VECTOR ang,pos,vel,com;
  REAL TotalMass,Mass;

  // compute system center of mass
  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      TotalMass+=(Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[j].Type].Mass);
      com.x+=Mass*Adsorbates[CurrentSystem][i].Atoms[j].Position.x;
      com.y+=Mass*Adsorbates[CurrentSystem][i].Atoms[j].Position.y;
      com.z+=Mass*Adsorbates[CurrentSystem][i].Atoms[j].Position.z;
    }

  }
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      TotalMass+=(Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[j].Type].Mass);
      com.x+=Mass*Cations[CurrentSystem][i].Atoms[j].Position.x;
      com.y+=Mass*Cations[CurrentSystem][i].Atoms[j].Position.y;
      com.z+=Mass*Cations[CurrentSystem][i].Atoms[j].Position.z;
    }
  }
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
    {
      Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[0][i].Type].Mass;
      TotalMass+=Mass;
      com.x+=Mass*Framework[CurrentSystem].Atoms[0][i].Position.x;
      com.y+=Mass*Framework[CurrentSystem].Atoms[0][i].Position.y;
      com.z+=Mass*Framework[CurrentSystem].Atoms[0][i].Position.z;
    }
  }
  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;

  ang.x=ang.y=ang.z=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[j].Type].Mass;
      pos.x=Adsorbates[CurrentSystem][i].Atoms[j].Position.x-com.x;
      pos.y=Adsorbates[CurrentSystem][i].Atoms[j].Position.y-com.y;
      pos.z=Adsorbates[CurrentSystem][i].Atoms[j].Position.z-com.z;
      vel=Adsorbates[CurrentSystem][i].Atoms[j].Velocity;
      ang.x+=Mass*(pos.y*vel.z-pos.z*vel.y);
      ang.y+=Mass*(pos.z*vel.x-pos.x*vel.z);
      ang.z+=Mass*(pos.x*vel.y-pos.y*vel.x);
    }
  }
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[j].Type].Mass;
      pos.x=Cations[CurrentSystem][i].Atoms[j].Position.x-com.x;
      pos.y=Cations[CurrentSystem][i].Atoms[j].Position.y-com.y;
      pos.z=Cations[CurrentSystem][i].Atoms[j].Position.z-com.z;
      vel=Cations[CurrentSystem][i].Atoms[j].Velocity;
      ang.x+=Mass*(pos.y*vel.z-pos.z*vel.y);
      ang.y+=Mass*(pos.z*vel.x-pos.x*vel.z);
      ang.z+=Mass*(pos.x*vel.y-pos.y*vel.x);
    }
  }
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
    {
      Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[0][i].Type].Mass;
      pos.x=Framework[CurrentSystem].Atoms[0][i].Position.x-com.x;
      pos.y=Framework[CurrentSystem].Atoms[0][i].Position.y-com.y;
      pos.z=Framework[CurrentSystem].Atoms[0][i].Position.z-com.z;
      vel=Framework[CurrentSystem].Atoms[0][i].Velocity;
      ang.x+=Mass*(pos.y*vel.z-pos.z*vel.y);
      ang.y+=Mass*(pos.z*vel.x-pos.x*vel.z);
      ang.z+=Mass*(pos.x*vel.y-pos.y*vel.x);
    }
  }
  return ang;
}


/***********************************************************************************************************
 * Name       | CreateCartesianPositions                                                                   *
 * ------------------------------------------------------------------------------------------------------- *
 * Function   | Calculates the Cartesian position from the center-of-mass and the quaternions.             *
 * Parameters | -                                                                                          *
 ***********************************************************************************************************/

void CreateCartesianPositions(void)
{
  REAL_MATRIX3x3 M;
  VECTOR pos,t,com;
  int i,k,l;
  int Type,A;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        BuildRotationMatrixInverse(&M,Adsorbates[CurrentSystem][i].Groups[l].Quaternion);

        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          pos=Components[Type].Positions[A];

          t.x=M.ax*pos.x+M.bx*pos.y+M.cx*pos.z;
          t.y=M.ay*pos.x+M.by*pos.y+M.cy*pos.z;
          t.z=M.az*pos.x+M.bz*pos.y+M.cz*pos.z;
          com=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition;
          Adsorbates[CurrentSystem][i].Atoms[A].Position.x=com.x+t.x;
          Adsorbates[CurrentSystem][i].Atoms[A].Position.y=com.y+t.y;
          Adsorbates[CurrentSystem][i].Atoms[A].Position.z=com.z+t.z;
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Type=Cations[CurrentSystem][i].Type;
        BuildRotationMatrixInverse(&M,Cations[CurrentSystem][i].Groups[l].Quaternion);

        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          pos=Components[Type].Positions[A];

          t.x=M.ax*pos.x+M.bx*pos.y+M.cx*pos.z;
          t.y=M.ay*pos.x+M.by*pos.y+M.cy*pos.z;
          t.z=M.az*pos.x+M.bz*pos.y+M.cz*pos.z;
          com=Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition;
          Cations[CurrentSystem][i].Atoms[A].Position.x=com.x+t.x;
          Cations[CurrentSystem][i].Atoms[A].Position.y=com.y+t.y;
          Cations[CurrentSystem][i].Atoms[A].Position.z=com.z+t.z;
        }
      }
    }
  }
}


/***********************************************************************************************************
 * Name       | CreateCartesianVelocities                                                                  *
 * ------------------------------------------------------------------------------------------------------- *
 * Function   | Calculates the Cartesian velocities from the center-of-mass velocity and the quaternion    *
 *            | momenta.                                                                                   *
 * Outline    | in body-frame: v = w x r, convert to space-fixed and add the com-velocity                  *
 * Parameters | -                                                                                          *
 ***********************************************************************************************************/

void CreateCartesianVelocities(void)
{
  int i,k,l;
  int Type,A;
  REAL_MATRIX3x3 M;
  VECTOR pos,omega,w,vel;
  QUATERNION q;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        q=Adsorbates[CurrentSystem][i].Groups[l].Quaternion;
        BuildRotationMatrixInverse(&M,q);
        omega=QuaternionMomentumToAngularVelocityAdsorbates(i,l);

        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];

          // site position in body frame
          pos=Components[Type].Positions[A];
          w.x=omega.y*pos.z-omega.z*pos.y;
          w.y=omega.z*pos.x-omega.x*pos.z;
          w.z=omega.x*pos.y-omega.y*pos.x;

          // new atomic velocites in lab frame
          vel=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x=(M.ax*w.x+M.bx*w.y+M.cx*w.z)+vel.x;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y=(M.ay*w.x+M.by*w.y+M.cy*w.z)+vel.y;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z=(M.az*w.x+M.bz*w.y+M.cz*w.z)+vel.z;
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        q=Cations[CurrentSystem][i].Groups[l].Quaternion;
        BuildRotationMatrixInverse(&M,q);
        omega=QuaternionMomentumToAngularVelocityCations(i,l);

        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          // site position in body frame
          pos=Components[Type].Positions[A];
          w.x=omega.y*pos.z-omega.z*pos.y;
          w.y=omega.z*pos.x-omega.x*pos.z;
          w.z=omega.x*pos.y-omega.y*pos.x;

          // new atomic velocites in lab frame
          vel=Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
          Cations[CurrentSystem][i].Atoms[A].Velocity.x=(M.ax*w.x+M.bx*w.y+M.cx*w.z)+vel.x;
          Cations[CurrentSystem][i].Atoms[A].Velocity.y=(M.ay*w.x+M.by*w.y+M.cy*w.z)+vel.y;
          Cations[CurrentSystem][i].Atoms[A].Velocity.z=(M.az*w.x+M.bz*w.y+M.cz*w.z)+vel.z;
        }
      }
    }
  }
}

void NoSquishRotate(int k,REAL dt)
{
  int i,l;
  int Type;
  REAL zeta;
  QUATERNION p,q,pn,qn;
  VECTOR I;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        p=Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum;
        q=Adsorbates[CurrentSystem][i].Groups[l].Quaternion;

        I=Components[Type].Groups[l].InverseInertiaVector;

        switch(k)
        {
          case 1:
            zeta=dt*(-p.r*q.i+p.i*q.r+p.j*q.k-p.k*q.j)*I.x/4.0;
            pn.r=cos(zeta)*p.r-sin(zeta)*p.i;
            pn.i=cos(zeta)*p.i+sin(zeta)*p.r;
            pn.j=cos(zeta)*p.j+sin(zeta)*p.k;
            pn.k=cos(zeta)*p.k-sin(zeta)*p.j;

            qn.r=cos(zeta)*q.r-sin(zeta)*q.i;
            qn.i=cos(zeta)*q.i+sin(zeta)*q.r;
            qn.j=cos(zeta)*q.j+sin(zeta)*q.k;
            qn.k=cos(zeta)*q.k-sin(zeta)*q.j;
            break;
          case 2:
            zeta=dt*(-p.r*q.j-p.i*q.k+p.j*q.r+p.k*q.i)*I.y/4.0;
            pn.r=cos(zeta)*p.r-sin(zeta)*p.j;
            pn.i=cos(zeta)*p.i-sin(zeta)*p.k;
            pn.j=cos(zeta)*p.j+sin(zeta)*p.r;
            pn.k=cos(zeta)*p.k+sin(zeta)*p.i;

            qn.r=cos(zeta)*q.r-sin(zeta)*q.j;
            qn.i=cos(zeta)*q.i-sin(zeta)*q.k;
            qn.j=cos(zeta)*q.j+sin(zeta)*q.r;
            qn.k=cos(zeta)*q.k+sin(zeta)*q.i;
            break;
          case 3:
            zeta=dt*(-p.r*q.k+p.i*q.j-p.j*q.i+p.k*q.r)*I.z/4.0;
            pn.r=cos(zeta)*p.r-sin(zeta)*p.k;
            pn.i=cos(zeta)*p.i+sin(zeta)*p.j;
            pn.j=cos(zeta)*p.j-sin(zeta)*p.i;
            pn.k=cos(zeta)*p.k+sin(zeta)*p.r;

            qn.r=cos(zeta)*q.r-sin(zeta)*q.k;
            qn.i=cos(zeta)*q.i+sin(zeta)*q.j;
            qn.j=cos(zeta)*q.j-sin(zeta)*q.i;
            qn.k=cos(zeta)*q.k+sin(zeta)*q.r;
            break;
          default:
            fprintf(stderr, "error\n");
            exit(0);
            break;
        }
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum=pn;
        Adsorbates[CurrentSystem][i].Groups[l].Quaternion=qn;
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        p=Cations[CurrentSystem][i].Groups[l].QuaternionMomentum;
        q=Cations[CurrentSystem][i].Groups[l].Quaternion;

        I=Components[Type].Groups[l].InverseInertiaVector;

        switch(k)
        {
          case 1:
            zeta=dt*(-p.r*q.i+p.i*q.r+p.j*q.k-p.k*q.j)*I.x/4.0;
            pn.r=cos(zeta)*p.r-sin(zeta)*p.i;
            pn.i=cos(zeta)*p.i+sin(zeta)*p.r;
            pn.j=cos(zeta)*p.j+sin(zeta)*p.k;
            pn.k=cos(zeta)*p.k-sin(zeta)*p.j;

            qn.r=cos(zeta)*q.r-sin(zeta)*q.i;
            qn.i=cos(zeta)*q.i+sin(zeta)*q.r;
            qn.j=cos(zeta)*q.j+sin(zeta)*q.k;
            qn.k=cos(zeta)*q.k-sin(zeta)*q.j;
            break;
          case 2:
            zeta=dt*(-p.r*q.j-p.i*q.k+p.j*q.r+p.k*q.i)*I.y/4.0;
            pn.r=cos(zeta)*p.r-sin(zeta)*p.j;
            pn.i=cos(zeta)*p.i-sin(zeta)*p.k;
            pn.j=cos(zeta)*p.j+sin(zeta)*p.r;
            pn.k=cos(zeta)*p.k+sin(zeta)*p.i;

            qn.r=cos(zeta)*q.r-sin(zeta)*q.j;
            qn.i=cos(zeta)*q.i-sin(zeta)*q.k;
            qn.j=cos(zeta)*q.j+sin(zeta)*q.r;
            qn.k=cos(zeta)*q.k+sin(zeta)*q.i;
            break;
          case 3:
            zeta=dt*(-p.r*q.k+p.i*q.j-p.j*q.i+p.k*q.r)*I.z/4.0;
            pn.r=cos(zeta)*p.r-sin(zeta)*p.k;
            pn.i=cos(zeta)*p.i+sin(zeta)*p.j;
            pn.j=cos(zeta)*p.j-sin(zeta)*p.i;
            pn.k=cos(zeta)*p.k+sin(zeta)*p.r;

            qn.r=cos(zeta)*q.r-sin(zeta)*q.k;
            qn.i=cos(zeta)*q.i+sin(zeta)*q.j;
            qn.j=cos(zeta)*q.j-sin(zeta)*q.i;
            qn.k=cos(zeta)*q.k+sin(zeta)*q.r;
            break;
          default:
            fprintf(stderr, "error\n");
            exit(0);
            break;
        }
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum=pn;
        Cations[CurrentSystem][i].Groups[l].Quaternion=qn;
      }
    }
  }
}

void NoSquishFreeRotorOrderTwo(void)
{
  int i,m_rot;

  // second order
  m_rot=5;
  for(i=0;i<m_rot;i++)
  {
    NoSquishRotate(3,0.5*DeltaT/(REAL)m_rot);
    NoSquishRotate(2,0.5*DeltaT/(REAL)m_rot);
    NoSquishRotate(1,DeltaT/(REAL)m_rot);
    NoSquishRotate(2,0.5*DeltaT/(REAL)m_rot);
    NoSquishRotate(3,0.5*DeltaT/(REAL)m_rot);
  }
}

void NoSquishFreeRotorOrderFour(void)
{
  int i,m_rot;
  REAL w[2];

  w[0]=1.0/(2.0-pow(2.0,1.0/3.0));
  w[1]=1.0-2.0*w[0];

  m_rot=5;
  for(i=0;i<m_rot;i++)
  {
    NoSquishRotate(3,w[0]*0.5*DeltaT/(REAL)m_rot);
    NoSquishRotate(2,w[0]*0.5*DeltaT/(REAL)m_rot);
    NoSquishRotate(1,w[0]*DeltaT/(REAL)m_rot);

    NoSquishRotate(2,w[0]*0.5*DeltaT/(REAL)m_rot);
    NoSquishRotate(3,(w[0]+w[1])*0.5*DeltaT/(REAL)m_rot);
    NoSquishRotate(2,w[1]*0.5*DeltaT/(REAL)m_rot);

    NoSquishRotate(1,w[1]*DeltaT/(REAL)m_rot);

    NoSquishRotate(2,w[1]*0.5*DeltaT/(REAL)m_rot);
    NoSquishRotate(3,(w[0]+w[1])*0.5*DeltaT/(REAL)m_rot);
    NoSquishRotate(2,w[0]*0.5*DeltaT/(REAL)m_rot);

    NoSquishRotate(1,w[0]*DeltaT/(REAL)m_rot);
    NoSquishRotate(2,w[0]*0.5*DeltaT/(REAL)m_rot);
    NoSquishRotate(3,w[0]*0.5*DeltaT/(REAL)m_rot);
  }
}

void CreateComVelocities(void)
{
  int i,j,l;
  int Type,A;
  REAL Mass,TotalMass;
  VECTOR COMVel;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        TotalMass=COMVel.x=COMVel.y=COMVel.z=0.0;
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
          TotalMass+=Mass;
          COMVel.x+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x;
          COMVel.y+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y;
          COMVel.z+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z;
        }
        COMVel.x/=TotalMass;
        COMVel.y/=TotalMass;
        COMVel.z/=TotalMass;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity=COMVel;
      }
    }
  }
}

void CreateQuaternionAndComForces(void)
{
  int i,j,l;
  int Type,A;
  VECTOR COMForce,Force,F;
  VECTOR dr,Torque;
  REAL Mass,MMass;
  REAL_MATRIX3x3 M;
  QUATERNION q,F_qua;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        BuildRotationMatrix(&M,Adsorbates[CurrentSystem][i].Groups[l].Quaternion);

        COMForce.x=COMForce.y=COMForce.z=0.0;
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          COMForce.x+=Adsorbates[CurrentSystem][i].Atoms[A].Force.x;
          COMForce.y+=Adsorbates[CurrentSystem][i].Atoms[A].Force.y;
          COMForce.z+=Adsorbates[CurrentSystem][i].Atoms[A].Force.z;
        }
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce=COMForce;

        Torque.x=0.0;
        Torque.y=0.0;
        Torque.z=0.0;
        MMass=Components[Type].Groups[l].Mass;
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];

          Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
          Force.x=Adsorbates[CurrentSystem][i].Atoms[A].Force.x-COMForce.x*Mass/MMass;
          Force.y=Adsorbates[CurrentSystem][i].Atoms[A].Force.y-COMForce.y*Mass/MMass;
          Force.z=Adsorbates[CurrentSystem][i].Atoms[A].Force.z-COMForce.z*Mass/MMass;
          F.x=M.ax*Force.x+M.bx*Force.y+M.cx*Force.z;
          F.y=M.ay*Force.x+M.by*Force.y+M.cy*Force.z;
          F.z=M.az*Force.x+M.bz*Force.y+M.cz*Force.z;
          dr=Components[Adsorbates[CurrentSystem][i].Type].Positions[A];
          Torque.x+=dr.y*F.z-dr.z*F.y;
          Torque.y+=dr.z*F.x-dr.x*F.z;
          Torque.z+=dr.x*F.y-dr.y*F.x;
        }
        Adsorbates[CurrentSystem][i].Groups[l].Torque=Torque;

        q=Adsorbates[CurrentSystem][i].Groups[l].Quaternion;
        F_qua.r=2.0*(-q.i*Torque.x-q.j*Torque.y-q.k*Torque.z);
        F_qua.i=2.0*(q.r*Torque.x-q.k*Torque.y+q.j*Torque.z);
        F_qua.j=2.0*(q.k*Torque.x+q.r*Torque.y-q.i*Torque.z);
        F_qua.k=2.0*(-q.j*Torque.x+q.i*Torque.y+q.r*Torque.z);
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce=F_qua;
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        BuildRotationMatrix(&M,Cations[CurrentSystem][i].Groups[l].Quaternion);

        COMForce.x=COMForce.y=COMForce.z=0.0;
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          COMForce.x+=Cations[CurrentSystem][i].Atoms[A].Force.x;
          COMForce.y+=Cations[CurrentSystem][i].Atoms[A].Force.y;
          COMForce.z+=Cations[CurrentSystem][i].Atoms[A].Force.z;
        }
        Cations[CurrentSystem][i].Groups[l].CenterOfMassForce=COMForce;

        Torque.x=0.0;
        Torque.y=0.0;
        Torque.z=0.0;
        MMass=Components[Type].Groups[l].Mass;
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];

          Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[A].Type].Mass;
          Force.x=Cations[CurrentSystem][i].Atoms[A].Force.x-COMForce.x*Mass/MMass;
          Force.y=Cations[CurrentSystem][i].Atoms[A].Force.y-COMForce.y*Mass/MMass;
          Force.z=Cations[CurrentSystem][i].Atoms[A].Force.z-COMForce.z*Mass/MMass;
          F.x=M.ax*Force.x+M.bx*Force.y+M.cx*Force.z;
          F.y=M.ay*Force.x+M.by*Force.y+M.cy*Force.z;
          F.z=M.az*Force.x+M.bz*Force.y+M.cz*Force.z;
          dr=Components[Cations[CurrentSystem][i].Type].Positions[A];
          Torque.x+=dr.y*F.z-dr.z*F.y;
          Torque.y+=dr.z*F.x-dr.x*F.z;
          Torque.z+=dr.x*F.y-dr.y*F.x;
        }
        Cations[CurrentSystem][i].Groups[l].Torque=Torque;

        q=Cations[CurrentSystem][i].Groups[l].Quaternion;
        F_qua.r=2.0*(-q.i*Torque.x-q.j*Torque.y-q.k*Torque.z);
        F_qua.i=2.0*(q.r*Torque.x-q.k*Torque.y+q.j*Torque.z);
        F_qua.j=2.0*(q.k*Torque.x+q.r*Torque.y-q.i*Torque.z);
        F_qua.k=2.0*(-q.j*Torque.x+q.i*Torque.y+q.r*Torque.z);
        Cations[CurrentSystem][i].Groups[l].QuaternionForce=F_qua;
      }
    }
  }
}

void NoSquishEvolve0Todt2(void)
{
  int i,l;
  int Type;
  REAL Mass;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Groups[l].Mass;

        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.r+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.r;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.i+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.i;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.j+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.j;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.k+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.k;

        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.x/Mass;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.y/Mass;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.z/Mass;
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Groups[l].Mass;

        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.r+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.r;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.i+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.i;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.j+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.j;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.k+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.k;

        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.x/Mass;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.y/Mass;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.z/Mass;
      }
    }
  }
}

void NoSquishEvolvedt2Todt(void)
{
  int i,l;
  int Type;
  REAL Mass;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Groups[l].Mass;

        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.r+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.r;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.i+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.i;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.j+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.j;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.k+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.k;

        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.x/Mass;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.y/Mass;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce.z/Mass;
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Groups[l].Mass;

        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.r+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.r;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.i+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.i;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.j+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.j;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.k+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.k;

        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.x/Mass;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.y/Mass;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].CenterOfMassForce.z/Mass;
      }
    }
  }
}

void ComputeDegreesOfFreedom(void)
{
  int i,j,l,f1;
  int Type;

  for(j=0;j<NumberOfSystems;j++)
  {
    DegreesOfFreedomTranslation[j]=0;
    DegreesOfFreedomRotation[j]=0;

    DegreesOfFreedomAdsorbates[j]=0;
    DegreesOfFreedomTranslationalAdsorbates[j]=0;
    DegreesOfFreedomRotationalAdsorbates[j]=0;

    DegreesOfFreedomCations[j]=0;
    DegreesOfFreedomTranslationalCations[j]=0;
    DegreesOfFreedomRotationalCations[j]=0;

    for(i=0;i<NumberOfAdsorbateMolecules[j];i++)
    {
      Type=Adsorbates[j][i].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid)
        {
          DegreesOfFreedomAdsorbates[j]+=3;
          DegreesOfFreedomTranslation[j]+=3;
          DegreesOfFreedomTranslationalAdsorbates[j]+=3;

          DegreesOfFreedomRotation[j]+=Components[Type].Groups[l].RotationalDegreesOfFreedom;
          DegreesOfFreedomAdsorbates[j]+=Components[Type].Groups[l].RotationalDegreesOfFreedom;
          DegreesOfFreedomRotationalAdsorbates[j]+=Components[Type].Groups[l].RotationalDegreesOfFreedom;

        }
        else
        {
          DegreesOfFreedomTranslation[j]+=3*Components[Type].Groups[l].NumberOfGroupAtoms;
          DegreesOfFreedomAdsorbates[j]+=3*Components[Type].Groups[l].NumberOfGroupAtoms;
          DegreesOfFreedomTranslationalAdsorbates[j]+=3*Components[Type].Groups[l].NumberOfGroupAtoms;
        }
      }

      DegreesOfFreedomAdsorbates[j]-=Components[Type].NumberOfConstraintBonds;
      DegreesOfFreedomTranslation[j]-=Components[Type].NumberOfConstraintBonds;
      DegreesOfFreedomTranslationalAdsorbates[j]-=Components[Type].NumberOfConstraintBonds;

      DegreesOfFreedomAdsorbates[j]-=Components[Type].NumberOfConstraintBends;
      DegreesOfFreedomTranslation[j]-=Components[Type].NumberOfConstraintBends;
      DegreesOfFreedomTranslationalAdsorbates[j]-=Components[Type].NumberOfConstraintBends;

      DegreesOfFreedomAdsorbates[j]-=Components[Type].NumberOfConstraintInversionBends;
      DegreesOfFreedomTranslation[j]-=Components[Type].NumberOfConstraintInversionBends;
      DegreesOfFreedomTranslationalAdsorbates[j]-=Components[Type].NumberOfConstraintInversionBends;

      DegreesOfFreedomAdsorbates[j]-=Components[Type].NumberOfConstraintTorsions;
      DegreesOfFreedomTranslation[j]-=Components[Type].NumberOfConstraintTorsions;
      DegreesOfFreedomTranslationalAdsorbates[j]-=Components[Type].NumberOfConstraintTorsions;

      DegreesOfFreedomAdsorbates[j]-=Components[Type].NumberOfConstraintImproperTorsions;
      DegreesOfFreedomTranslation[j]-=Components[Type].NumberOfConstraintImproperTorsions;
      DegreesOfFreedomTranslationalAdsorbates[j]-=Components[Type].NumberOfConstraintImproperTorsions;
    }

    DegreesOfFreedomCations[j]=0;
    for(i=0;i<NumberOfCationMolecules[j];i++)
    {
      Type=Cations[j][i].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid)
        {
          DegreesOfFreedomCations[j]+=3;
          DegreesOfFreedomTranslation[j]+=3;
          DegreesOfFreedomTranslationalCations[j]+=3;

          DegreesOfFreedomRotation[j]+=Components[Type].Groups[l].RotationalDegreesOfFreedom;
          DegreesOfFreedomCations[j]+=Components[Type].Groups[l].RotationalDegreesOfFreedom;
          DegreesOfFreedomRotationalCations[j]+=Components[Type].Groups[l].RotationalDegreesOfFreedom;
        }
        else
        {
          DegreesOfFreedomTranslation[j]+=3*Components[Type].Groups[l].NumberOfGroupAtoms;
          DegreesOfFreedomCations[j]+=3*Components[Type].Groups[l].NumberOfGroupAtoms;
          DegreesOfFreedomTranslationalCations[j]+=3*Components[Type].Groups[l].NumberOfGroupAtoms;
        }
      }
      DegreesOfFreedomCations[j]-=Components[Type].NumberOfConstraintBonds;
      DegreesOfFreedomTranslation[j]-=Components[Type].NumberOfConstraintBonds;
      DegreesOfFreedomTranslationalCations[j]-=Components[Type].NumberOfConstraintBonds;

      DegreesOfFreedomCations[j]-=Components[Type].NumberOfConstraintBends;
      DegreesOfFreedomTranslation[j]-=Components[Type].NumberOfConstraintBends;
      DegreesOfFreedomTranslationalCations[j]-=Components[Type].NumberOfConstraintBends;

      DegreesOfFreedomCations[j]-=Components[Type].NumberOfConstraintInversionBends;
      DegreesOfFreedomTranslation[j]-=Components[Type].NumberOfConstraintInversionBends;
      DegreesOfFreedomTranslationalCations[j]-=Components[Type].NumberOfConstraintInversionBends;

      DegreesOfFreedomCations[j]-=Components[Type].NumberOfConstraintTorsions;
      DegreesOfFreedomTranslation[j]-=Components[Type].NumberOfConstraintTorsions;
      DegreesOfFreedomTranslationalCations[j]-=Components[Type].NumberOfConstraintTorsions;

      DegreesOfFreedomCations[j]-=Components[Type].NumberOfConstraintImproperTorsions;
      DegreesOfFreedomTranslation[j]-=Components[Type].NumberOfConstraintImproperTorsions;
      DegreesOfFreedomTranslationalCations[j]-=Components[Type].NumberOfConstraintImproperTorsions;

    }

    if(Framework[j].FrameworkModel==FLEXIBLE)
    {
      for(f1=0;f1<Framework[j].NumberOfFrameworks;f1++)
      {
        DegreesOfFreedomFramework[j]+=3*(Framework[j].NumberOfFreeAtoms[f1]-Framework[j].NumberOfCoreShells[f1]);
        DegreesOfFreedomTranslation[j]+=3*(Framework[j].NumberOfFreeAtoms[f1]-Framework[j].NumberOfCoreShells[f1]);
      }
    }
    else
      DegreesOfFreedomFramework[j]=0;

    DegreesOfFreedom[j]=DegreesOfFreedomTranslation[j]+DegreesOfFreedomRotation[j];
  }

}

void AdjustSystemAngularRotationToZero(void)
{
  int i,j,l;
  int Type;
  VECTOR com,am,pos,vel,w;
  REAL Mass,TotalMass,rsq,det;
  REAL_MATRIX3x3 roti,rotinv;

    com=MeasureVelocityDrift();
    RemoveVelocityDrift();
    com=MeasureVelocityDrift();
    //if(FrameworkModel==FLEXIBLE)
    //  DegreesOfFreedom[CurrentSystem]-=3;

    // compute system center of mass
    TotalMass=0.0;
    com.x=com.y=com.z=0.0;
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        TotalMass+=(Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[j].Type].Mass);
        com.x+=Mass*Adsorbates[CurrentSystem][i].Atoms[j].Position.x;
        com.y+=Mass*Adsorbates[CurrentSystem][i].Atoms[j].Position.y;
        com.z+=Mass*Adsorbates[CurrentSystem][i].Atoms[j].Position.z;
      }
    }
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        TotalMass+=(Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[j].Type].Mass);
        com.x+=Mass*Cations[CurrentSystem][i].Atoms[j].Position.x;
        com.y+=Mass*Cations[CurrentSystem][i].Atoms[j].Position.y;
        com.z+=Mass*Cations[CurrentSystem][i].Atoms[j].Position.z;
      }
    }
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
      {
        Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[0][i].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Framework[CurrentSystem].Atoms[0][i].Position.x;
        com.y+=Mass*Framework[CurrentSystem].Atoms[0][i].Position.y;
        com.z+=Mass*Framework[CurrentSystem].Atoms[0][i].Position.z;
      }
    }

    com.x/=TotalMass;
    com.y/=TotalMass;
    com.z/=TotalMass;

    // shift to center of mass origin
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Adsorbates[CurrentSystem][i].Atoms[j].Position.x-=com.x;
        Adsorbates[CurrentSystem][i].Atoms[j].Position.y-=com.y;
        Adsorbates[CurrentSystem][i].Atoms[j].Position.z-=com.z;
      }
    }
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Cations[CurrentSystem][i].Atoms[j].Position.x-=com.x;
        Cations[CurrentSystem][i].Atoms[j].Position.y-=com.y;
        Cations[CurrentSystem][i].Atoms[j].Position.z-=com.z;
      }
    }
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
      {
        Framework[CurrentSystem].Atoms[0][i].Position.x-=com.x;
        Framework[CurrentSystem].Atoms[0][i].Position.y-=com.y;
        Framework[CurrentSystem].Atoms[0][i].Position.z-=com.z;
      }
    }

    roti.ax=roti.bx=roti.cx=0.0;
    roti.ay=roti.by=roti.cy=0.0;
    roti.az=roti.bz=roti.cz=0.0;
    am.x=am.y=am.z=0.0;

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[j].Type].Mass;
        pos=Adsorbates[CurrentSystem][i].Atoms[j].Position;
        vel=Adsorbates[CurrentSystem][i].Atoms[j].Velocity;
        am.x+=Mass*(pos.y*vel.z-pos.z*vel.y);
        am.y+=Mass*(pos.z*vel.x-pos.x*vel.z);
        am.z+=Mass*(pos.x*vel.y-pos.y*vel.x);

        rsq=SQR(pos.x)+SQR(pos.y)+SQR(pos.z);
        roti.ax+=Mass*(SQR(pos.x)-rsq);
        roti.ay+=Mass*pos.x*pos.y;
        roti.az+=Mass*pos.x*pos.z;
        roti.bx+=Mass*pos.y*pos.x;
        roti.by+=Mass*(SQR(pos.y)-rsq);
        roti.bz+=Mass*pos.y*pos.z;
        roti.cx+=Mass*pos.z*pos.x;
        roti.cy+=Mass*pos.z*pos.y;
        roti.cz+=Mass*(SQR(pos.z)-rsq);
      }
    }

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[j].Type].Mass;
        pos=Cations[CurrentSystem][i].Atoms[j].Position;
        vel=Cations[CurrentSystem][i].Atoms[j].Velocity;
        am.x+=Mass*(pos.y*vel.z-pos.z*vel.y);
        am.y+=Mass*(pos.z*vel.x-pos.x*vel.z);
        am.z+=Mass*(pos.x*vel.y-pos.y*vel.x);
        rsq=SQR(pos.x)+SQR(pos.y)+SQR(pos.z);
        roti.ax+=Mass*(SQR(pos.x)-rsq);
        roti.ay+=Mass*pos.x*pos.y;
        roti.az+=Mass*pos.x*pos.z;
        roti.bx+=Mass*pos.y*pos.x;
        roti.by+=Mass*(SQR(pos.y)-rsq);
        roti.bz+=Mass*pos.y*pos.z;
        roti.cx+=Mass*pos.z*pos.x;
        roti.cy+=Mass*pos.z*pos.y;
        roti.cz+=Mass*(SQR(pos.z)-rsq);
      }
    }
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
      {
        Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[0][i].Type].Mass;
        pos=Framework[CurrentSystem].Atoms[0][i].Position;
        vel=Framework[CurrentSystem].Atoms[0][i].Velocity;
        am.x+=Mass*(pos.y*vel.z-pos.z*vel.y);
        am.y+=Mass*(pos.z*vel.x-pos.x*vel.z);
        am.z+=Mass*(pos.x*vel.y-pos.y*vel.x);
        rsq=SQR(pos.x)+SQR(pos.y)+SQR(pos.z);
        roti.ax+=Mass*(SQR(pos.x)-rsq);
        roti.ay+=Mass*pos.x*pos.y;
        roti.az+=Mass*pos.x*pos.z;
        roti.bx+=Mass*pos.y*pos.x;
        roti.by+=Mass*(SQR(pos.y)-rsq);
        roti.bz+=Mass*pos.y*pos.z;
        roti.cx+=Mass*pos.z*pos.x;
        roti.cy+=Mass*pos.z*pos.y;
        roti.cz+=Mass*(SQR(pos.z)-rsq);
      }
    }

    Invert3x3Matrix(&roti,&rotinv,&det);

    w.x=rotinv.ax*am.x+rotinv.ay*am.y+rotinv.az*am.z;
    w.y=rotinv.bx*am.x+rotinv.by*am.y+rotinv.bz*am.z;
    w.z=rotinv.cx*am.x+rotinv.cy*am.y+rotinv.cz*am.z;

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        pos=Adsorbates[CurrentSystem][i].Atoms[j].Position;
        vel=Adsorbates[CurrentSystem][i].Atoms[j].Velocity;
        vel.x+=w.y*pos.z-w.z*pos.y;
        vel.y+=w.z*pos.x-w.x*pos.z;
        vel.z+=w.x*pos.y-w.y*pos.x;
        Adsorbates[CurrentSystem][i].Atoms[j].Velocity=vel;
      }
    }

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        pos=Cations[CurrentSystem][i].Atoms[j].Position;
        vel=Cations[CurrentSystem][i].Atoms[j].Velocity;
        vel.x+=w.y*pos.z-w.z*pos.y;
        vel.y+=w.z*pos.x-w.x*pos.z;
        vel.z+=w.x*pos.y-w.y*pos.x;
        Cations[CurrentSystem][i].Atoms[j].Velocity=vel;
      }
    }
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
      {
        pos=Framework[CurrentSystem].Atoms[0][i].Position;
        vel=Framework[CurrentSystem].Atoms[0][i].Velocity;
        vel.x+=w.y*pos.z-w.z*pos.y;
        vel.y+=w.z*pos.x-w.x*pos.z;
        vel.z+=w.x*pos.y-w.y*pos.x;
        Framework[CurrentSystem].Atoms[0][i].Velocity=vel;
      }
    }

    // reset positions to original reference frame
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Adsorbates[CurrentSystem][i].Atoms[j].Position.x+=com.x;
        Adsorbates[CurrentSystem][i].Atoms[j].Position.y+=com.y;
        Adsorbates[CurrentSystem][i].Atoms[j].Position.z+=com.z;
      }
    }
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Cations[CurrentSystem][i].Atoms[j].Position.x+=com.x;
        Cations[CurrentSystem][i].Atoms[j].Position.y+=com.y;
        Cations[CurrentSystem][i].Atoms[j].Position.z+=com.z;
      }
    }
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
      {
        Framework[CurrentSystem].Atoms[0][i].Position.x+=com.x;
        Framework[CurrentSystem].Atoms[0][i].Position.y+=com.y;
        Framework[CurrentSystem].Atoms[0][i].Position.z+=com.z;
      }
    }

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      Type=Adsorbates[CurrentSystem][i].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity=AtomicVelocityToAngularVelocityAdsorbates(i,l);
    }

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      Type=Cations[CurrentSystem][i].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
        Cations[CurrentSystem][i].Groups[l].AngularVelocity=AtomicVelocityToAngularVelocityCations(i,l);
    }

    vel=TotalAngularMomentum();
}

// ==============================================================================================================================
// Routines for the use of SHAKE in minimization and RATTLE in molecular dynamics
// ==============================================================================================================================

// auxiliary routines for the distance-constraints
//------------------------------------------------

REAL ReturnConstraintDistance(VECTOR posA,VECTOR posB)
{
  VECTOR dr;
  REAL rr,r;

  dr.x=posA.x-posB.x;
  dr.y=posA.y-posB.y;
  dr.z=posA.z-posB.z;
  dr=ApplyBoundaryCondition(dr);
  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

  switch(DistanceConstraintType)
  {
    case DISTANCE_R:
      r=sqrt(rr);
      return r;
    default:
    case DISTANCE_R_SQUARED:
      return rr;
  }
}

REAL ReturnWilsonVectorsDistanceRATTLE(VECTOR posA,VECTOR posB,VECTOR *wa,VECTOR *wb)
{
  VECTOR Rba;
  REAL r2,rab;

  Rba.x=posA.x-posB.x;
  Rba.y=posA.y-posB.y;
  Rba.z=posA.z-posB.z;
  Rba=ApplyBoundaryCondition(Rba);
  r2=SQR(Rba.x)+SQR(Rba.y)+SQR(Rba.z);

  switch(DistanceConstraintType)
  {
    case DISTANCE_R:
      rab=sqrt(r2);
      wa->x=Rba.x/rab;
      wa->y=Rba.y/rab;
      wa->z=Rba.z/rab;
      wb->x=-Rba.x/rab;
      wb->y=-Rba.y/rab;
      wb->z=-Rba.z/rab;
      return rab;
    default:
    case DISTANCE_R_SQUARED:
      wa->x=2.0*Rba.x;
      wa->y=2.0*Rba.y;
      wa->z=2.0*Rba.z;
      wb->x=-2.0*Rba.x;
      wb->y=-2.0*Rba.y;
      wb->z=-2.0*Rba.z;
      return r2;
  }
}


REAL ReturnDistanceConstrainDerivative(VECTOR Rba,VECTOR Vba)
{
  REAL deriv;

  switch(DistanceConstraintType)
  {
    case DISTANCE_R:
      deriv=(Rba.x*Vba.x+Rba.y*Vba.y+Rba.z*Vba.z)/sqrt(SQR(Rba.x)+SQR(Rba.y)+SQR(Rba.z));
      return deriv;
    default:
    case DISTANCE_R_SQUARED:
      deriv=2.0*(Rba.x*Vba.x+Rba.y*Vba.y+Rba.z*Vba.z);
      return deriv;
  }
}


// auxiliary routines for the bend-angle-constraints
//--------------------------------------------------

REAL ReturnConstraintBendAngle(VECTOR posA,VECTOR posB,VECTOR posC)
{
  REAL CosTheta;
  VECTOR Rba,Rbc;
  REAL rab,rbc;

  Rba.x=posA.x-posB.x;
  Rba.y=posA.y-posB.y;
  Rba.z=posA.z-posB.z;
  Rba=ApplyBoundaryCondition(Rba);
  rab=sqrt(SQR(Rba.x)+SQR(Rba.y)+SQR(Rba.z));
  Rba.x/=rab;
  Rba.y/=rab;
  Rba.z/=rab;

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  Rbc=ApplyBoundaryCondition(Rbc);
  rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x/=rbc;
  Rbc.y/=rbc;
  Rbc.z/=rbc;

  CosTheta=(Rba.x*Rbc.x+Rba.y*Rbc.y+Rba.z*Rbc.z);

  switch(BendConstraintType)
  {
    case COS_THETA_SQUARED:
      return SQR(CosTheta);
    case COS_THETA:
      return CosTheta;
    case THETA:
    default:
      return acos(CosTheta);
  }

}

// the Wilson vectors are the conversion from internal theta-angle derivative to Cartesian derivatives: d\phi\dr_{1,2,3}
REAL ReturnWilsonVectorsBendRATTLE(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR *wa,VECTOR *wb,VECTOR *wc)
{
  VECTOR dtA,dtB,dtC;
  VECTOR Rba,Rbc;
  REAL rab,rbc;
  REAL CosTheta,SinTheta;

  Rba.x=posA.x-posB.x;
  Rba.y=posA.y-posB.y;
  Rba.z=posA.z-posB.z;
  Rba=ApplyBoundaryCondition(Rba);
  rab=sqrt(SQR(Rba.x)+SQR(Rba.y)+SQR(Rba.z));
  Rba.x/=rab;
  Rba.y/=rab;
  Rba.z/=rab;

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  Rbc=ApplyBoundaryCondition(Rbc);
  rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x/=rbc;
  Rbc.y/=rbc;
  Rbc.z/=rbc;

  CosTheta=(Rba.x*Rbc.x+Rba.y*Rbc.y+Rba.z*Rbc.z);
  SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));

  // Calculate the components of the derivatives.
  dtA.x=(CosTheta*Rba.x-Rbc.x)/rab;
  dtA.y=(CosTheta*Rba.y-Rbc.y)/rab;
  dtA.z=(CosTheta*Rba.z-Rbc.z)/rab;

  dtC.x=(CosTheta*Rbc.x-Rba.x)/rbc;
  dtC.y=(CosTheta*Rbc.y-Rba.y)/rbc;
  dtC.z=(CosTheta*Rbc.z-Rba.z)/rbc;

  dtB.x=-(dtA.x+dtC.x);
  dtB.y=-(dtA.y+dtC.y);
  dtB.z=-(dtA.z+dtC.z);


  switch(BendConstraintType)
  {
    case COS_THETA_SQUARED:
      wa->x=-2.0*dtA.x*CosTheta;
      wa->y=-2.0*dtA.y*CosTheta;
      wa->z=-2.0*dtA.z*CosTheta;

      wb->x=-2.0*dtB.x*CosTheta;
      wb->y=-2.0*dtB.y*CosTheta;
      wb->z=-2.0*dtB.z*CosTheta;

      wc->x=-2.0*dtC.x*CosTheta;
      wc->y=-2.0*dtC.y*CosTheta;
      wc->z=-2.0*dtC.z*CosTheta;
      return SQR(CosTheta);
    case COS_THETA:
      wa->x=-dtA.x;
      wa->y=-dtA.y;
      wa->z=-dtA.z;

      wb->x=-dtB.x;
      wb->y=-dtB.y;
      wb->z=-dtB.z;

      wc->x=-dtC.x;
      wc->y=-dtC.y;
      wc->z=-dtC.z;
      return CosTheta;
    default:
    case THETA:
      wa->x=dtA.x/SinTheta;
      wa->y=dtA.y/SinTheta;
      wa->z=dtA.z/SinTheta;

      wb->x=dtB.x/SinTheta;
      wb->y=dtB.y/SinTheta;
      wb->z=dtB.z/SinTheta;

      wc->x=dtC.x/SinTheta;
      wc->y=dtC.y/SinTheta;
      wc->z=dtC.z/SinTheta;
      return acos(CosTheta);
  }
}

REAL ReturnAngleConstrainDerivative(VECTOR Rba,VECTOR Rbc,VECTOR Vba,VECTOR Vbc)
{
  REAL deriv;
  REAL CosTheta,SinTheta;
  REAL dot_product1,dot_product2;
  REAL dot_product3,dot_product4;
  REAL length1,length2;
  REAL length_squared1,length_squared2;
  REAL term1,term2;

  length_squared1=SQR(Rba.x)+SQR(Rba.y)+SQR(Rba.z);
  length_squared2=SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z);

  length1=sqrt(length_squared1);
  length2=sqrt(length_squared2);

  CosTheta=(Rba.x*Rbc.x+Rba.y*Rbc.y+Rba.z*Rbc.z)/(length1*length2);
  SinTheta=sqrt(1.0-SQR(CosTheta));

  dot_product1=Rba.x*Vba.x+Rba.y*Vba.y+Rba.z*Vba.z;
  dot_product2=Rbc.x*Vbc.x+Rbc.y*Vbc.y+Rbc.z*Vbc.z;
  dot_product3=Rba.x*Vbc.x+Rba.y*Vbc.y+Rba.z*Vbc.z;
  dot_product4=Rbc.x*Vba.x+Rbc.y*Vba.y+Rbc.z*Vba.z;

  term1=(dot_product3+dot_product4)/(length1*length2);
  term2=(dot_product1/length_squared1+dot_product2/length_squared2)*CosTheta;

  switch(BendConstraintType)
  {
    case COS_THETA_SQUARED:
      deriv=2.0*CosTheta*(term1-term2);
      return deriv;
    case THETA:
      deriv=(term2-term1)/SinTheta;
      return deriv;
    case COS_THETA:
    default:
      deriv=(term1-term2);
      return deriv;
  }
}

// auxiliary routines for the dihedral-angle-constraints
//------------------------------------------------------

REAL ReturnConstraintDihedralAngle(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD)
{
  REAL rbc;
  VECTOR Rba,Rbc,Rcd,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi;
  VECTOR Pb,Pc;

  Rba.x=posA.x-posB.x;
  Rba.y=posA.y-posB.y;
  Rba.z=posA.z-posB.z;
  Rba=ApplyBoundaryCondition(Rba);

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  Rbc=ApplyBoundaryCondition(Rbc);
  rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x/=rbc; Rbc.y/=rbc; Rbc.z/=rbc;

  Rcd.x=posD.x-posC.x;
  Rcd.y=posD.y-posC.y;
  Rcd.z=posD.z-posC.z;
  Rcd=ApplyBoundaryCondition(Rcd);

  dot_ab=Rba.x*Rbc.x+Rba.y*Rbc.y+Rba.z*Rbc.z;
  dot_cd=Rcd.x*Rbc.x+Rcd.y*Rbc.y+Rcd.z*Rbc.z;

  dr.x=Rba.x-dot_ab*Rbc.x;
  dr.y=Rba.y-dot_ab*Rbc.y;
  dr.z=Rba.z-dot_ab*Rbc.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  dr.x/=r; dr.y/=r; dr.z/=r;

  ds.x=Rcd.x-dot_cd*Rbc.x;
  ds.y=Rcd.y-dot_cd*Rbc.y;
  ds.z=Rcd.z-dot_cd*Rbc.z;
  s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
  ds.x/=s; ds.y/=s; ds.z/=s;

  // compute Cos(Phi)
  // Phi is defined in protein convention Phi(trans)=Pi
  CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

  // Ensure CosPhi is between -1 and 1.
  CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);

  switch(DihedralConstraintType)
  {
    case COS_PHI_SQUARED:
      return SQR(CosPhi);
    case COS_PHI:
      return CosPhi;
    case PHI:
    default:
      // potential defined in terms of 'phi' and therefore contains a singularity
      // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
      // same direction as Rbc, and negative otherwise
      Pb.x=Rba.z*Rbc.y-Rba.y*Rbc.z;
      Pb.y=Rba.x*Rbc.z-Rba.z*Rbc.x;
      Pb.z=Rba.y*Rbc.x-Rba.x*Rbc.y;
      Pc.x=Rbc.y*Rcd.z-Rbc.z*Rcd.y;
      Pc.y=Rbc.z*Rcd.x-Rbc.x*Rcd.z;
      Pc.z=Rbc.x*Rcd.y-Rbc.y*Rcd.x;
      sign=(Rbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Rbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)+Rbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);
      if(Phi<=-M_PI) Phi+=2.0*M_PI;
      return Phi;
  }
}


// the Wilson vectors are the conversion from internal theta-angle derivative to Cartesian derivatives: d\phi\dr_{1,2,3}
REAL ReturnWilsonVectorsTorsionRATTLE(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR *wa,VECTOR *wb,VECTOR *wc,VECTOR *wd)
{
  REAL dot_ab,dot_cd;
  VECTOR Rba,Rbc,Rcd;
  VECTOR dr,ds;
  REAL r,s,rbc,CosPhi,SinPhi,d,e;
  VECTOR dtA,dtB,dtC,dtD;
  VECTOR Pb,Pc;
  REAL sign,Phi;

  Rba.x=posA.x-posB.x;
  Rba.y=posA.y-posB.y;
  Rba.z=posA.z-posB.z;
  Rba=ApplyBoundaryCondition(Rba);

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  Rbc=ApplyBoundaryCondition(Rbc);
  rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x/=rbc; Rbc.y/=rbc; Rbc.z/=rbc;

  Rcd.x=posD.x-posC.x;
  Rcd.y=posD.y-posC.y;
  Rcd.z=posD.z-posC.z;
  Rcd=ApplyBoundaryCondition(Rcd);

  dot_ab=Rba.x*Rbc.x+Rba.y*Rbc.y+Rba.z*Rbc.z;
  dot_cd=Rcd.x*Rbc.x+Rcd.y*Rbc.y+Rcd.z*Rbc.z;

  dr.x=Rba.x-dot_ab*Rbc.x;
  dr.y=Rba.y-dot_ab*Rbc.y;
  dr.z=Rba.z-dot_ab*Rbc.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  dr.x/=r; dr.y/=r; dr.z/=r;

  ds.x=Rcd.x-dot_cd*Rbc.x;
  ds.y=Rcd.y-dot_cd*Rbc.y;
  ds.z=Rcd.z-dot_cd*Rbc.z;
  s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
  ds.x/=s; ds.y/=s; ds.z/=s;

  // compute Cos(Phi)
  // Phi is defined in protein convention Phi(trans)=Pi
  CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

  // Ensure CosPhi is between -1 and 1.
  CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);

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

  switch(DihedralConstraintType)
  {
    case COS_PHI_SQUARED:
      wa->x=2.0*CosPhi*dtA.x;
      wa->y=2.0*CosPhi*dtA.y;
      wa->z=2.0*CosPhi*dtA.z;

      wb->x=2.0*CosPhi*dtB.x;
      wb->y=2.0*CosPhi*dtB.y;
      wb->z=2.0*CosPhi*dtB.z;

      wc->x=2.0*CosPhi*dtC.x;
      wc->y=2.0*CosPhi*dtC.y;
      wc->z=2.0*CosPhi*dtC.z;

      wd->x=2.0*CosPhi*dtD.x;
      wd->y=2.0*CosPhi*dtD.y;
      wd->z=2.0*CosPhi*dtD.z;
      return SQR(CosPhi);
    case COS_PHI:
      wa->x=dtA.x;
      wa->y=dtA.y;
      wa->z=dtA.z;

      wb->x=dtB.x;
      wb->y=dtB.y;
      wb->z=dtB.z;

      wc->x=dtC.x;
      wc->y=dtC.y;
      wc->z=dtC.z;

      wd->x=dtD.x;
      wd->y=dtD.y;
      wd->z=dtD.z;
      return CosPhi;
    default:
    case PHI:
      Pb.x=Rba.z*Rbc.y-Rba.y*Rbc.z;
      Pb.y=Rba.x*Rbc.z-Rba.z*Rbc.x;
      Pb.z=Rba.y*Rbc.x-Rba.x*Rbc.y;
      Pc.x=Rbc.y*Rcd.z-Rbc.z*Rcd.y;
      Pc.y=Rbc.z*Rcd.x-Rbc.x*Rcd.z;
      Pc.z=Rbc.x*Rcd.y-Rbc.y*Rcd.x;
      sign=(Rbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Rbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)+Rbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
      Phi=SIGN(acos(CosPhi),sign);
      if(Phi<=-M_PI) Phi+=2.0*M_PI;
      SinPhi=sin(Phi);

      wa->x=-dtA.x/SinPhi;
      wa->y=-dtA.y/SinPhi;
      wa->z=-dtA.z/SinPhi;

      wb->x=-dtB.x/SinPhi;
      wb->y=-dtB.y/SinPhi;
      wb->z=-dtB.z/SinPhi;

      wc->x=-dtC.x/SinPhi;
      wc->y=-dtC.y/SinPhi;
      wc->z=-dtC.z/SinPhi;

      wd->x=-dtD.x/SinPhi;
      wd->y=-dtD.y/SinPhi;
      wd->z=-dtD.z/SinPhi;
      return Phi;
  }
}

// d\sigma/d\lambda with sigma=acos[(r_ab x r_bc).(r_bc x r_cd)/(|r_ab x r_bc| |r_bc x r_cd|)]
REAL ReturnTorsionConstrainDerivative(VECTOR Rab,VECTOR Rbc,VECTOR Rcd,VECTOR Vab,VECTOR Vbc,VECTOR Vcd)
{
  REAL term1,term2,Phi,CosPhi,sign;
  VECTOR CrossProduct1,CrossProduct2;
  VECTOR CrossProduct3,CrossProduct4;
  REAL CrossProduct1LengthSquared,CrossProduct2LengthSquared;
  REAL CrossProduct1Length,CrossProduct2Length;

  CrossProduct1.x=Rab.y*Rbc.z-Rab.z*Rbc.y; CrossProduct2.x=Rbc.y*Rcd.z-Rbc.z*Rcd.y;
  CrossProduct1.y=Rab.z*Rbc.x-Rab.x*Rbc.z; CrossProduct2.y=Rbc.z*Rcd.x-Rbc.x*Rcd.z;
  CrossProduct1.z=Rab.x*Rbc.y-Rab.y*Rbc.x; CrossProduct2.z=Rbc.x*Rcd.y-Rbc.y*Rcd.x;
  CrossProduct1LengthSquared=SQR(CrossProduct1.x)+SQR(CrossProduct1.y)+SQR(CrossProduct1.z);
  CrossProduct2LengthSquared=SQR(CrossProduct2.x)+SQR(CrossProduct2.y)+SQR(CrossProduct2.z);
  CrossProduct1Length=sqrt(CrossProduct1LengthSquared);
  CrossProduct2Length=sqrt(CrossProduct2LengthSquared);

  CrossProduct3.x=Rcd.y*Vbc.z-Rcd.z*Vbc.y+Vcd.y*Rbc.z-Vcd.z*Rbc.y;
  CrossProduct3.y=Rcd.z*Vbc.x-Rcd.x*Vbc.z+Vcd.z*Rbc.x-Vcd.x*Rbc.z;
  CrossProduct3.z=Rcd.x*Vbc.y-Rcd.y*Vbc.x+Vcd.x*Rbc.y-Vcd.y*Rbc.x;

  CrossProduct4.x=Rbc.y*Vab.z-Rbc.z*Vab.y+Vbc.y*Rab.z-Vbc.z*Rab.y;
  CrossProduct4.y=Rbc.z*Vab.x-Rbc.x*Vab.z+Vbc.z*Rab.x-Vbc.x*Rab.z;
  CrossProduct4.z=Rbc.x*Vab.y-Rbc.y*Vab.x+Vbc.x*Rab.y-Vbc.y*Rab.x;

  term1=(CrossProduct1.x*CrossProduct3.x+CrossProduct1.y*CrossProduct3.y+CrossProduct1.z*CrossProduct3.z+
         CrossProduct2.x*CrossProduct4.x+CrossProduct2.y*CrossProduct4.y+CrossProduct2.z*CrossProduct4.z)/
        (CrossProduct1Length*CrossProduct2Length);
  term2=(CrossProduct1.x*CrossProduct4.x+CrossProduct1.y*CrossProduct4.y+CrossProduct1.z*CrossProduct4.z)/CrossProduct1LengthSquared+
        (CrossProduct2.x*CrossProduct3.x+CrossProduct2.y*CrossProduct3.y+CrossProduct2.z*CrossProduct3.z)/CrossProduct2LengthSquared;

  CosPhi=(CrossProduct1.x*CrossProduct2.x+CrossProduct1.y*CrossProduct2.y+CrossProduct1.z*CrossProduct2.z)/(CrossProduct1Length*CrossProduct2Length);
  CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);

  switch(DihedralConstraintType)
  {
    case COS_PHI_SQUARED:
      return -2.0*CosPhi*(term1-term2*CosPhi);
    case COS_PHI:
      return -(term1-term2*CosPhi);
    default:
    case PHI:
      // the proper sign for Phi is needed
      sign=Rbc.x*(CrossProduct1.y*CrossProduct2.z-CrossProduct1.z*CrossProduct2.y)
          +Rbc.y*(CrossProduct1.z*CrossProduct2.x-CrossProduct1.x*CrossProduct2.z)
          +Rbc.z*(CrossProduct1.x*CrossProduct2.y-CrossProduct1.y*CrossProduct2.x);
      Phi=SIGN(acos(CosPhi),sign);
      return (term1-term2*cos(Phi))/sin(Phi);
  }
}



// auxiliary routines for the Inversion-Bend-angle-constraints
//------------------------------------------------------------

// Note: the angle can be positive or negative for chiral molecules
// right-handed is positive, left-handed is negative
REAL ReturnConstraintInversionBendAngle(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD)
{
  REAL r;
  VECTOR cross_product_dc,cross_product_ad,cross_product_ca;
  REAL dot_product;
  VECTOR Rba,Rbc,Rbd;
  REAL SinChiA,SinChiC,SinChiD;

  Rba.x=posA.x-posB.x;
  Rba.y=posA.y-posB.y;
  Rba.z=posA.z-posB.z;
  Rba=ApplyBoundaryCondition(Rba);
  r=sqrt(SQR(Rba.x)+SQR(Rba.y)+SQR(Rba.z));
  Rba.x/=r;
  Rba.y/=r;
  Rba.z/=r;

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  Rbc=ApplyBoundaryCondition(Rbc);
  r=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x/=r;
  Rbc.y/=r;
  Rbc.z/=r;

  Rbd.x=posD.x-posB.x;
  Rbd.y=posD.y-posB.y;
  Rbd.z=posD.z-posB.z;
  Rbd=ApplyBoundaryCondition(Rbd);
  r=sqrt(SQR(Rbd.x)+SQR(Rbd.y)+SQR(Rbd.z));
  Rbd.x/=r;
  Rbd.y/=r;
  Rbd.z/=r;

  cross_product_dc.x=Rbd.y*Rbc.z-Rbd.z*Rbc.y;
  cross_product_dc.y=Rbd.z*Rbc.x-Rbd.x*Rbc.z;
  cross_product_dc.z=Rbd.x*Rbc.y-Rbd.y*Rbc.x;
  dot_product=cross_product_dc.x*Rba.x+cross_product_dc.y*Rba.y+cross_product_dc.z*Rba.z;
  SinChiA=dot_product/sin(ReturnBendAngle(posC,posB,posD));

  cross_product_ad.x=Rba.y*Rbd.z-Rba.z*Rbd.y;
  cross_product_ad.y=Rba.z*Rbd.x-Rba.x*Rbd.z;
  cross_product_ad.z=Rba.x*Rbd.y-Rba.y*Rbd.x;
  dot_product=cross_product_ad.x*Rbc.x+cross_product_ad.y*Rbc.y+cross_product_ad.z*Rbc.z;
  SinChiC=dot_product/sin(ReturnBendAngle(posA,posB,posD));

  cross_product_ca.x=Rbc.y*Rba.z-Rbc.z*Rba.y;
  cross_product_ca.y=Rbc.z*Rba.x-Rbc.x*Rba.z;
  cross_product_ca.z=Rbc.x*Rba.y-Rbc.y*Rba.x;
  dot_product=cross_product_ca.x*Rbd.x+cross_product_ca.y*Rbd.y+cross_product_ca.z*Rbd.z;
  SinChiD=dot_product/sin(ReturnBendAngle(posA,posB,posC));

  switch(InversionBendConstraintType)
  {
    case SIN_CHI_SQUARED:
      return (SQR(SinChiA)+SQR(SinChiC)+SQR(SinChiD))/3.0;
    case SIN_CHI:
      return (SinChiA+SinChiC+SinChiD)/3.0;
    default:
    case PHI:
      return (asin(SinChiA)+asin(SinChiC)+asin(SinChiD))/3.0;
  }
}

REAL ReturnWilsonVectorsInversionBendRATTLE(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR *wa,VECTOR *wb,VECTOR *wc,VECTOR *wd)
{
  REAL rab,rbc,rbd;
  VECTOR cross_product_dc,cross_product_ad,cross_product_ca;
  REAL dot_product;
  REAL ThetaA,SinThetaA,CosThetaA,ChiA,CosChiA,SinChiA,TanChiA;
  REAL ThetaC,SinThetaC,CosThetaC,ChiC,CosChiC,SinChiC,TanChiC;
  REAL ThetaD,SinThetaD,CosThetaD,ChiD,CosChiD,SinChiD,TanChiD;
  VECTOR Rba,Rbc,Rbd;

  Rba.x=posA.x-posB.x;
  Rba.y=posA.y-posB.y;
  Rba.z=posA.z-posB.z;
  Rba=ApplyBoundaryCondition(Rba);
  rab=sqrt(SQR(Rba.x)+SQR(Rba.y)+SQR(Rba.z));
  Rba.x/=rab;
  Rba.y/=rab;
  Rba.z/=rab;

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  Rbc=ApplyBoundaryCondition(Rbc);
  rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x/=rbc;
  Rbc.y/=rbc;
  Rbc.z/=rbc;

  Rbd.x=posD.x-posB.x;
  Rbd.y=posD.y-posB.y;
  Rbd.z=posD.z-posB.z;
  Rbd=ApplyBoundaryCondition(Rbd);
  rbd=sqrt(SQR(Rbd.x)+SQR(Rbd.y)+SQR(Rbd.z));
  Rbd.x/=rbd;
  Rbd.y/=rbd;
  Rbd.z/=rbd;

  ThetaA=ReturnBendAngle(posC,posB,posD);
  CosThetaA=cos(ReturnBendAngle(posC,posB,posD));
  SinThetaA=sin(ReturnBendAngle(posC,posB,posD));

  ThetaC=ReturnBendAngle(posA,posB,posD);
  CosThetaC=cos(ReturnBendAngle(posA,posB,posD));
  SinThetaC=sin(ReturnBendAngle(posA,posB,posD));

  ThetaD=ReturnBendAngle(posA,posB,posC);
  CosThetaD=cos(ReturnBendAngle(posA,posB,posC));
  SinThetaD=sin(ReturnBendAngle(posA,posB,posC));

  cross_product_dc.x=Rbd.y*Rbc.z-Rbd.z*Rbc.y;
  cross_product_dc.y=Rbd.z*Rbc.x-Rbd.x*Rbc.z;
  cross_product_dc.z=Rbd.x*Rbc.y-Rbd.y*Rbc.x;
  dot_product=cross_product_dc.x*Rba.x+cross_product_dc.y*Rba.y+cross_product_dc.z*Rba.z;
  ChiA=asin(dot_product/SinThetaA);
  CosChiA=cos(ChiA);
  SinChiA=sin(ChiA);
  TanChiA=tan(ChiA);

  cross_product_ad.x=Rba.y*Rbd.z-Rba.z*Rbd.y;
  cross_product_ad.y=Rba.z*Rbd.x-Rba.x*Rbd.z;
  cross_product_ad.z=Rba.x*Rbd.y-Rba.y*Rbd.x;
  dot_product=cross_product_ad.x*Rbc.x+cross_product_ad.y*Rbc.y+cross_product_ad.z*Rbc.z;
  ChiC=asin(dot_product/sin(ReturnBendAngle(posA,posB,posD)));
  CosChiC=cos(ChiC);
  SinChiC=sin(ChiC);
  TanChiC=tan(ChiC);

  cross_product_ca.x=Rbc.y*Rba.z-Rbc.z*Rba.y;
  cross_product_ca.y=Rbc.z*Rba.x-Rbc.x*Rba.z;
  cross_product_ca.z=Rbc.x*Rba.y-Rbc.y*Rba.x;
  dot_product=cross_product_ca.x*Rbd.x+cross_product_ca.y*Rbd.y+cross_product_ca.z*Rbd.z;
  ChiD=asin(dot_product/sin(ReturnBendAngle(posA,posB,posC)));
  CosChiD=cos(ChiD);
  SinChiD=sin(ChiD);
  TanChiD=tan(ChiD);

  switch(InversionBendConstraintType)
  {
    case SIN_CHI_SQUARED:
      wa->x=2.0*SinChiA*(cross_product_dc.x/SinThetaA-Rba.x*SinChiA)/(3.0*rab);
      wa->y=2.0*SinChiA*(cross_product_dc.y/SinThetaA-Rba.y*SinChiA)/(3.0*rab);
      wa->z=2.0*SinChiA*(cross_product_dc.z/SinThetaA-Rba.z*SinChiA)/(3.0*rab);

      wa->x+=2.0*SinChiD*(cross_product_dc.x/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rba.x-Rbc.x*CosThetaD))/(3.0*rab);
      wa->y+=2.0*SinChiD*(cross_product_dc.y/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rba.y-Rbc.y*CosThetaD))/(3.0*rab);
      wa->z+=2.0*SinChiD*(cross_product_dc.z/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rba.z-Rbc.z*CosThetaD))/(3.0*rab);

      wa->x+=2.0*SinChiC*(cross_product_dc.x/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rba.x-Rbd.x*CosThetaC))/(3.0*rab);
      wa->y+=2.0*SinChiC*(cross_product_dc.y/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rba.y-Rbd.y*CosThetaC))/(3.0*rab);
      wa->z+=2.0*SinChiC*(cross_product_dc.z/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rba.z-Rbd.z*CosThetaC))/(3.0*rab);


      wc->x=2.0*SinChiA*(cross_product_ad.x/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbc.x-Rbd.x*CosThetaA))/(3.0*rbc);
      wc->y=2.0*SinChiA*(cross_product_ad.y/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbc.y-Rbd.y*CosThetaA))/(3.0*rbc);
      wc->z=2.0*SinChiA*(cross_product_ad.z/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbc.z-Rbd.z*CosThetaA))/(3.0*rbc);

      wc->x+=2.0*SinChiD*(cross_product_ad.x/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rbc.x-Rba.x*CosThetaD))/(3.0*rbc);
      wc->y+=2.0*SinChiD*(cross_product_ad.y/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rbc.y-Rba.y*CosThetaD))/(3.0*rbc);
      wc->z+=2.0*SinChiD*(cross_product_ad.z/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rbc.z-Rba.z*CosThetaD))/(3.0*rbc);

      wc->x+=2.0*SinChiC*(cross_product_ad.x/SinThetaC-Rbc.x*SinChiC)/(3.0*rbc);
      wc->y+=2.0*SinChiC*(cross_product_ad.y/SinThetaC-Rbc.y*SinChiC)/(3.0*rbc);
      wc->z+=2.0*SinChiC*(cross_product_ad.z/SinThetaC-Rbc.z*SinChiC)/(3.0*rbc);


      wd->x=2.0*SinChiA*(cross_product_ca.x/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbd.x-Rbc.x*CosThetaA))/(3.0*rbd);
      wd->y=2.0*SinChiA*(cross_product_ca.y/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbd.y-Rbc.y*CosThetaA))/(3.0*rbd);
      wd->z=2.0*SinChiA*(cross_product_ca.z/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbd.z-Rbc.z*CosThetaA))/(3.0*rbd);

      wd->x+=2.0*SinChiD*(cross_product_ca.x/SinThetaD-Rbd.x*SinChiD)/(3.0*rbd);
      wd->y+=2.0*SinChiD*(cross_product_ca.y/SinThetaD-Rbd.y*SinChiD)/(3.0*rbd);
      wd->z+=2.0*SinChiD*(cross_product_ca.z/SinThetaD-Rbd.z*SinChiD)/(3.0*rbd);

      wd->x+=2.0*SinChiC*(cross_product_ca.x/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rbd.x-Rba.x*CosThetaC))/(3.0*rbd);
      wd->y+=2.0*SinChiC*(cross_product_ca.y/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rbd.y-Rba.y*CosThetaC))/(3.0*rbd);
      wd->z+=2.0*SinChiC*(cross_product_ca.z/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rbd.z-Rba.z*CosThetaC))/(3.0*rbd);
      break;
    case SIN_CHI:
      wa->x=(cross_product_dc.x/SinThetaA-Rba.x*SinChiA)/(3.0*rab);
      wa->y=(cross_product_dc.y/SinThetaA-Rba.y*SinChiA)/(3.0*rab);
      wa->z=(cross_product_dc.z/SinThetaA-Rba.z*SinChiA)/(3.0*rab);

      wa->x+=(cross_product_dc.x/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rba.x-Rbc.x*CosThetaD))/(3.0*rab);
      wa->y+=(cross_product_dc.y/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rba.y-Rbc.y*CosThetaD))/(3.0*rab);
      wa->z+=(cross_product_dc.z/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rba.z-Rbc.z*CosThetaD))/(3.0*rab);

      wa->x+=(cross_product_dc.x/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rba.x-Rbd.x*CosThetaC))/(3.0*rab);
      wa->y+=(cross_product_dc.y/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rba.y-Rbd.y*CosThetaC))/(3.0*rab);
      wa->z+=(cross_product_dc.z/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rba.z-Rbd.z*CosThetaC))/(3.0*rab);


      wc->x=(cross_product_ad.x/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbc.x-Rbd.x*CosThetaA))/(3.0*rbc);
      wc->y=(cross_product_ad.y/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbc.y-Rbd.y*CosThetaA))/(3.0*rbc);
      wc->z=(cross_product_ad.z/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbc.z-Rbd.z*CosThetaA))/(3.0*rbc);

      wc->x+=(cross_product_ad.x/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rbc.x-Rba.x*CosThetaD))/(3.0*rbc);
      wc->y+=(cross_product_ad.y/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rbc.y-Rba.y*CosThetaD))/(3.0*rbc);
      wc->z+=(cross_product_ad.z/SinThetaD-(SinChiD/(SQR(SinThetaD)))*(Rbc.z-Rba.z*CosThetaD))/(3.0*rbc);

      wc->x+=(cross_product_ad.x/SinThetaC-Rbc.x*SinChiC)/(3.0*rbc);
      wc->y+=(cross_product_ad.y/SinThetaC-Rbc.y*SinChiC)/(3.0*rbc);
      wc->z+=(cross_product_ad.z/SinThetaC-Rbc.z*SinChiC)/(3.0*rbc);


      wd->x=(cross_product_ca.x/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbd.x-Rbc.x*CosThetaA))/(3.0*rbd);
      wd->y=(cross_product_ca.y/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbd.y-Rbc.y*CosThetaA))/(3.0*rbd);
      wd->z=(cross_product_ca.z/SinThetaA-(SinChiA/(SQR(SinThetaA)))*(Rbd.z-Rbc.z*CosThetaA))/(3.0*rbd);

      wd->x+=(cross_product_ca.x/SinThetaD-Rbd.x*SinChiD)/(3.0*rbd);
      wd->y+=(cross_product_ca.y/SinThetaD-Rbd.y*SinChiD)/(3.0*rbd);
      wd->z+=(cross_product_ca.z/SinThetaD-Rbd.z*SinChiD)/(3.0*rbd);

      wd->x+=(cross_product_ca.x/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rbd.x-Rba.x*CosThetaC))/(3.0*rbd);
      wd->y+=(cross_product_ca.y/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rbd.y-Rba.y*CosThetaC))/(3.0*rbd);
      wd->z+=(cross_product_ca.z/SinThetaC-(SinChiC/(SQR(SinThetaC)))*(Rbd.z-Rba.z*CosThetaC))/(3.0*rbd);
      break;
    default:
    case PHI:
      wa->x=(cross_product_dc.x/(CosChiA*SinThetaA)-Rba.x*TanChiA)/(3.0*rab);
      wa->y=(cross_product_dc.y/(CosChiA*SinThetaA)-Rba.y*TanChiA)/(3.0*rab);
      wa->z=(cross_product_dc.z/(CosChiA*SinThetaA)-Rba.z*TanChiA)/(3.0*rab);

      wa->x+=(cross_product_dc.x/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Rba.x-Rbc.x*CosThetaD))/(3.0*rab);
      wa->y+=(cross_product_dc.y/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Rba.y-Rbc.y*CosThetaD))/(3.0*rab);
      wa->z+=(cross_product_dc.z/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Rba.z-Rbc.z*CosThetaD))/(3.0*rab);

      wa->x+=(cross_product_dc.x/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Rba.x-Rbd.x*CosThetaC))/(3.0*rab);
      wa->y+=(cross_product_dc.y/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Rba.y-Rbd.y*CosThetaC))/(3.0*rab);
      wa->z+=(cross_product_dc.z/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Rba.z-Rbd.z*CosThetaC))/(3.0*rab);


      wc->x=(cross_product_ad.x/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Rbc.x-Rbd.x*CosThetaA))/(3.0*rbc);
      wc->y=(cross_product_ad.y/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Rbc.y-Rbd.y*CosThetaA))/(3.0*rbc);
      wc->z=(cross_product_ad.z/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Rbc.z-Rbd.z*CosThetaA))/(3.0*rbc);

      wc->x+=(cross_product_ad.x/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Rbc.x-Rba.x*CosThetaD))/(3.0*rbc);
      wc->y+=(cross_product_ad.y/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Rbc.y-Rba.y*CosThetaD))/(3.0*rbc);
      wc->z+=(cross_product_ad.z/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Rbc.z-Rba.z*CosThetaD))/(3.0*rbc);

      wc->x+=(cross_product_ad.x/(CosChiC*SinThetaC)-Rbc.x*TanChiC)/(3.0*rbc);
      wc->y+=(cross_product_ad.y/(CosChiC*SinThetaC)-Rbc.y*TanChiC)/(3.0*rbc);
      wc->z+=(cross_product_ad.z/(CosChiC*SinThetaC)-Rbc.z*TanChiC)/(3.0*rbc);


      wd->x=(cross_product_ca.x/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Rbd.x-Rbc.x*CosThetaA))/(3.0*rbd);
      wd->y=(cross_product_ca.y/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Rbd.y-Rbc.y*CosThetaA))/(3.0*rbd);
      wd->z=(cross_product_ca.z/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Rbd.z-Rbc.z*CosThetaA))/(3.0*rbd);

      wd->x+=(cross_product_ca.x/(CosChiD*SinThetaD)-Rbd.x*TanChiD)/(3.0*rbd);
      wd->y+=(cross_product_ca.y/(CosChiD*SinThetaD)-Rbd.y*TanChiD)/(3.0*rbd);
      wd->z+=(cross_product_ca.z/(CosChiD*SinThetaD)-Rbd.z*TanChiD)/(3.0*rbd);

      wd->x+=(cross_product_ca.x/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Rbd.x-Rba.x*CosThetaC))/(3.0*rbd);
      wd->y+=(cross_product_ca.y/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Rbd.y-Rba.y*CosThetaC))/(3.0*rbd);
      wd->z+=(cross_product_ca.z/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Rbd.z-Rba.z*CosThetaC))/(3.0*rbd);
      break;
  }

  wb->x=-(wa->x+wc->x+wd->x);
  wb->y=-(wa->y+wc->y+wd->y);
  wb->z=-(wa->z+wc->z+wd->z);

  return (ChiA+ChiC+ChiD)/3.0;
}

REAL ReturnInversionBendConstrainDerivative(VECTOR Rba,VECTOR Rbc,VECTOR Rbd,VECTOR Vba,VECTOR Vbc,VECTOR Vbd)
{
  REAL deriv1,deriv2,deriv3;
  REAL SinThetaABC,SinThetaCBD,SinThetaDBA;
  REAL CosThetaABC,CosThetaCBD,CosThetaDBA;
  REAL SinChiA,SinChiC,SinChiD,CosChiA,CosChiC,CosChiD;
  REAL LengthRba,LengthRbc,LengthRbd;
  REAL LengthRbaSquared,LengthRbcSquared,LengthRbdSquared;
  REAL temp;

  LengthRbaSquared=SQR(Rba.x)+SQR(Rba.y)+SQR(Rba.z);
  LengthRbcSquared=SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z);
  LengthRbdSquared=SQR(Rbd.x)+SQR(Rbd.y)+SQR(Rbd.z);
  LengthRba=sqrt(LengthRbaSquared);
  LengthRbc=sqrt(LengthRbcSquared);
  LengthRbd=sqrt(LengthRbdSquared);

  CosThetaABC=(Rba.x*Rbc.x+Rba.y*Rbc.y+Rba.z*Rbc.z)/(LengthRba*LengthRbc);
  CosThetaCBD=(Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z)/(LengthRbc*LengthRbd);
  CosThetaDBA=(Rbd.x*Rba.x+Rbd.y*Rba.y+Rbd.z*Rba.z)/(LengthRbd*LengthRba);

  SinThetaABC=sqrt(1.0-SQR(Rba.x*Rbc.x+Rba.y*Rbc.y+Rba.z*Rbc.z)/(LengthRbaSquared*LengthRbcSquared));
  SinThetaCBD=sqrt(1.0-SQR(Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z)/(LengthRbcSquared*LengthRbdSquared));
  SinThetaDBA=sqrt(1.0-SQR(Rbd.x*Rba.x+Rbd.y*Rba.y+Rbd.z*Rba.z)/(LengthRbdSquared*LengthRbaSquared));

  SinChiA=(Rba.x*(Rbd.y*Rbc.z-Rbd.z*Rbc.y)+Rba.y*(Rbd.z*Rbc.x-Rbd.x*Rbc.z)+Rba.z*(Rbd.x*Rbc.y-Rbd.y*Rbc.x))/(LengthRba*SinThetaCBD*LengthRbc*LengthRbd);
  SinChiC=(Rbc.x*(Rba.y*Rbd.z-Rba.z*Rbd.y)+Rbc.y*(Rba.z*Rbd.x-Rba.x*Rbd.z)+Rbc.z*(Rba.x*Rbd.y-Rba.y*Rbd.x))/(LengthRba*SinThetaDBA*LengthRbc*LengthRbd);
  SinChiD=(Rbd.x*(Rbc.y*Rba.z-Rbc.z*Rba.y)+Rbd.y*(Rbc.z*Rba.x-Rbc.x*Rba.z)+Rbd.z*(Rbc.x*Rba.y-Rbc.y*Rba.x))/(LengthRba*SinThetaABC*LengthRbc*LengthRbd);

  CosChiA=sqrt(1.0-SQR(SinChiA));
  CosChiC=sqrt(1.0-SQR(SinChiC));
  CosChiD=sqrt(1.0-SQR(SinChiD));

  // derivative for Sin(Chi^a)
  temp=SQR(CosThetaCBD)*(Rbc.x*Vbc.x + Rbc.y*Vbc.y + Rbc.z*Vbc.z)/LengthRbcSquared
      +SQR(CosThetaCBD)*(Rbd.x*Vbd.x+Rbd.y*Vbd.y+Rbd.z*Vbd.z)/LengthRbdSquared
      -CosThetaCBD*(Rbd.x*Vbc.x+Rbd.y*Vbc.y+Rbd.z*Vbc.z + Rbc.x*Vbd.x+Rbc.y*Vbd.y+Rbc.z*Vbd.z)/(LengthRbc*LengthRbd);

  deriv1=(Vba.x*(Rbc.z*Rbd.y-Rbc.y*Rbd.z)+Rba.x*((Rbd.y*Vbc.z-Rbd.z*Vbc.y)+(Vbd.y*Rbc.z-Vbd.z*Rbc.y)))/(LengthRba*SinThetaCBD*LengthRbc*LengthRbd)
        +(Vba.y*(Rbc.x*Rbd.z-Rbc.z*Rbd.x)+Rba.y*((Rbd.z*Vbc.x-Rbd.x*Vbc.z)+(Vbd.z*Rbc.x-Vbd.x*Rbc.z)))/(LengthRba*SinThetaCBD*LengthRbc*LengthRbd)
        +(Vba.z*(Rbc.y*Rbd.x-Rbc.x*Rbd.y)+Rba.z*((Rbd.x*Vbc.y-Rbd.y*Vbc.x)+(Vbd.x*Rbc.y-Vbd.y*Rbc.x)))/(LengthRba*SinThetaCBD*LengthRbc*LengthRbd)
        -SinChiA*((Rba.x*Vba.x+Rba.y*Vba.y+Rba.z*Vba.z)/SQR(LengthRba)
                 +(Rbc.x*Vbc.x+Rbc.y*Vbc.y+Rbc.z*Vbc.z)/SQR(LengthRbc)
                 +(Rbd.x*Vbd.x+Rbd.y*Vbd.y+Rbd.z*Vbd.z)/SQR(LengthRbd)
                 +temp/SQR(SinThetaCBD));

  // derivative for Sin(Chi^c)
  temp=SQR(CosThetaDBA)*(Rbd.x*Vbd.x + Rbd.y*Vbd.y + Rbd.z*Vbd.z)/LengthRbdSquared
      +SQR(CosThetaDBA)*(Rba.x*Vba.x+Rba.y*Vba.y+Rba.z*Vba.z)/LengthRbaSquared
      -CosThetaDBA*(Rba.x*Vbd.x+Rba.y*Vbd.y+Rba.z*Vbd.z + Rbd.x*Vba.x+Rbd.y*Vba.y+Rbd.z*Vba.z)/(LengthRbd*LengthRba);

  deriv2=(Vbc.x*(Rbd.z*Rba.y-Rbd.y*Rba.z)+Rbc.x*((Rba.y*Vbd.z-Rba.z*Vbd.y)+(Vba.y*Rbd.z-Vba.z*Rbd.y)))/(LengthRba*SinThetaDBA*LengthRbc*LengthRbd)
        +(Vbc.y*(Rbd.x*Rba.z-Rbd.z*Rba.x)+Rbc.y*((Rba.z*Vbd.x-Rba.x*Vbd.z)+(Vba.z*Rbd.x-Vba.x*Rbd.z)))/(LengthRba*SinThetaDBA*LengthRbc*LengthRbd)
        +(Vbc.z*(Rbd.y*Rba.x-Rbd.x*Rba.y)+Rbc.z*((Rba.x*Vbd.y-Rba.y*Vbd.x)+(Vba.x*Rbd.y-Vba.y*Rbd.x)))/(LengthRba*SinThetaDBA*LengthRbc*LengthRbd)
        -SinChiC*((Rbc.x*Vbc.x+Rbc.y*Vbc.y+Rbc.z*Vbc.z)/SQR(LengthRbc)
                 +(Rbd.x*Vbd.x+Rbd.y*Vbd.y+Rbd.z*Vbd.z)/SQR(LengthRbd)
                 +(Rba.x*Vba.x+Rba.y*Vba.y+Rba.z*Vba.z)/SQR(LengthRba)
                 +temp/SQR(SinThetaDBA));


  // derivative for Sin(Chi^d)
  temp=SQR(CosThetaABC)*(Rba.x*Vba.x + Rba.y*Vba.y + Rba.z*Vba.z)/LengthRbaSquared
      +SQR(CosThetaABC)*(Rbc.x*Vbc.x+Rbc.y*Vbc.y+Rbc.z*Vbc.z)/LengthRbcSquared
      -CosThetaABC*(Rbc.x*Vba.x+Rbc.y*Vba.y+Rbc.z*Vba.z + Rba.x*Vbc.x+Rba.y*Vbc.y+Rba.z*Vbc.z)/(LengthRba*LengthRbc);

  deriv3=(Vbd.x*(Rba.z*Rbc.y-Rba.y*Rbc.z)+Rbd.x*((Rbc.y*Vba.z-Rbc.z*Vba.y)+(Vbc.y*Rba.z-Vbc.z*Rba.y)))/(LengthRba*SinThetaABC*LengthRbc*LengthRbd)
        +(Vbd.y*(Rba.x*Rbc.z-Rba.z*Rbc.x)+Rbd.y*((Rbc.z*Vba.x-Rbc.x*Vba.z)+(Vbc.z*Rba.x-Vbc.x*Rba.z)))/(LengthRba*SinThetaABC*LengthRbc*LengthRbd)
        +(Vbd.z*(Rba.y*Rbc.x-Rba.x*Rbc.y)+Rbd.z*((Rbc.x*Vba.y-Rbc.y*Vba.x)+(Vbc.x*Rba.y-Vbc.y*Rba.x)))/(LengthRba*SinThetaABC*LengthRbc*LengthRbd)
        -SinChiD*((Rbd.x*Vbd.x+Rbd.y*Vbd.y+Rbd.z*Vbd.z)/SQR(LengthRbd)
                 +(Rba.x*Vba.x+Rba.y*Vba.y+Rba.z*Vba.z)/SQR(LengthRba)
                 +(Rbc.x*Vbc.x+Rbc.y*Vbc.y+Rbc.z*Vbc.z)/SQR(LengthRbc)
                 +temp/SQR(SinThetaABC));


  switch(InversionBendConstraintType)
  {
    case SIN_CHI_SQUARED:
      return 2.0*(SinChiA*deriv1+SinChiC*deriv2+SinChiD*deriv3)/3.0;
    case SIN_CHI:
      return (deriv1+deriv2+deriv3)/3.0;
    default:
    case CHI:
      return (deriv1/CosChiA+deriv2/CosChiC+deriv3/CosChiD)/3.0;
  }

}

// auxiliary routines for the OutOfPlane-distance-constraints
//-----------------------------------------------------------

REAL ReturnConstraintOutOfPlaneDistance(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD)
{
  return 0.0;
}

REAL ReturnWilsonVectorsOutOfPlaneDistanceRATTLE(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR *wa,VECTOR *wb,VECTOR *wc,VECTOR *wd)
{
  return 0.0;
}

REAL ReturnOutOfPlaneDistanceConstrainDerivative(VECTOR Rij,VECTOR Rjk,VECTOR Rkl,VECTOR Vij,VECTOR Vjk,VECTOR Vkl)
{
  return 0.0;
}


// Note that shake makes no contribution to the stress (and the constraint forces are zero) when used in minimization.
// The routine is just used as a feasibility correction to satisfy the contraints
void ShakeInMinimization(void)
{
  int j,m,Type;
  int A,B,C,D;
  VECTOR posA,posB,posC,posD;
  REAL r,r0,gamma;
  REAL max_error,error;
  REAL MassA,MassB,MassC,MassD;
  VECTOR ha,hb,hc,hd;
  VECTOR Rab,Rbc,Rcd,Rba,Rbd, Vab,Vbc,Vcd,Vba,Vbd;
  REAL Phi,Phi0,Theta,Theta0,Chi,Chi0;
  REAL InverseMassA,InverseMassB,InverseMassC,InverseMassD;
  VECTOR Va,Vb,Vc,Vd;

  CurrentSystem=0;

  // store the initial positions at the start of the routine
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(j=0;j<Components[Type].NumberOfAtoms;j++)
      Adsorbates[CurrentSystem][m].Atoms[j].RattleReferencePosition=Adsorbates[CurrentSystem][m].Atoms[j].Position;
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(j=0;j<Components[Type].NumberOfAtoms;j++)
      Cations[CurrentSystem][m].Atoms[j].RattleReferencePosition=Cations[CurrentSystem][m].Atoms[j].Position;
  }


  // store the initial distance constraint gradients
  for(j=0;j<NumberOfDistanceConstraints[CurrentSystem];j++)
  {
    MassA=PseudoAtoms[DistanceConstraints[CurrentSystem][j][0]->Type].Mass;
    MassB=PseudoAtoms[DistanceConstraints[CurrentSystem][j][1]->Type].Mass;

    posA=DistanceConstraints[CurrentSystem][j][0]->Position;
    posB=DistanceConstraints[CurrentSystem][j][1]->Position;

    ReturnWilsonVectorsDistanceRATTLE(posA,posB,&ha,&hb);

    DistanceConstraintsDerivatives[CurrentSystem][j][0].x=ha.x*SQR(DeltaT)/MassA;
    DistanceConstraintsDerivatives[CurrentSystem][j][0].y=ha.y*SQR(DeltaT)/MassA;
    DistanceConstraintsDerivatives[CurrentSystem][j][0].z=ha.z*SQR(DeltaT)/MassA;

    DistanceConstraintsDerivatives[CurrentSystem][j][1].x=hb.x*SQR(DeltaT)/MassB;
    DistanceConstraintsDerivatives[CurrentSystem][j][1].y=hb.y*SQR(DeltaT)/MassB;
    DistanceConstraintsDerivatives[CurrentSystem][j][1].z=hb.z*SQR(DeltaT)/MassB;
  }

  // store the initial bend-angle constraint gradients
  for(j=0;j<NumberOfAngleConstraints[CurrentSystem];j++)
  {
    MassA=PseudoAtoms[AngleConstraints[CurrentSystem][j][0]->Type].Mass;
    MassB=PseudoAtoms[AngleConstraints[CurrentSystem][j][1]->Type].Mass;
    MassC=PseudoAtoms[AngleConstraints[CurrentSystem][j][2]->Type].Mass;

    posA=AngleConstraints[CurrentSystem][j][0]->Position;
    posB=AngleConstraints[CurrentSystem][j][1]->Position;
    posC=AngleConstraints[CurrentSystem][j][2]->Position;

    ReturnWilsonVectorsBendRATTLE(posA,posB,posC,&ha,&hb,&hc);

    AngleConstraintsDerivatives[CurrentSystem][j][0].x=ha.x*SQR(DeltaT)/MassA;
    AngleConstraintsDerivatives[CurrentSystem][j][0].y=ha.y*SQR(DeltaT)/MassA;
    AngleConstraintsDerivatives[CurrentSystem][j][0].z=ha.z*SQR(DeltaT)/MassA;

    AngleConstraintsDerivatives[CurrentSystem][j][1].x=hb.x*SQR(DeltaT)/MassB;
    AngleConstraintsDerivatives[CurrentSystem][j][1].y=hb.y*SQR(DeltaT)/MassB;
    AngleConstraintsDerivatives[CurrentSystem][j][1].z=hb.z*SQR(DeltaT)/MassB;

    AngleConstraintsDerivatives[CurrentSystem][j][2].x=hc.x*SQR(DeltaT)/MassC;
    AngleConstraintsDerivatives[CurrentSystem][j][2].y=hc.y*SQR(DeltaT)/MassC;
    AngleConstraintsDerivatives[CurrentSystem][j][2].z=hc.z*SQR(DeltaT)/MassC;
  }

  // store the initial inversion-bend-angle constraint gradients
  for(j=0;j<NumberOfInversionBendConstraints[CurrentSystem];j++)
  {
    MassA=PseudoAtoms[InversionBendConstraints[CurrentSystem][j][0]->Type].Mass;
    MassB=PseudoAtoms[InversionBendConstraints[CurrentSystem][j][1]->Type].Mass;
    MassC=PseudoAtoms[InversionBendConstraints[CurrentSystem][j][2]->Type].Mass;
    MassD=PseudoAtoms[InversionBendConstraints[CurrentSystem][j][3]->Type].Mass;

    posA=InversionBendConstraints[CurrentSystem][j][0]->Position;
    posB=InversionBendConstraints[CurrentSystem][j][1]->Position;
    posC=InversionBendConstraints[CurrentSystem][j][2]->Position;
    posD=InversionBendConstraints[CurrentSystem][j][3]->Position;

    ReturnWilsonVectorsInversionBendRATTLE(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

    InversionBendConstraintsDerivatives[CurrentSystem][j][0].x=ha.x*SQR(DeltaT)/MassA;
    InversionBendConstraintsDerivatives[CurrentSystem][j][0].y=ha.y*SQR(DeltaT)/MassA;
    InversionBendConstraintsDerivatives[CurrentSystem][j][0].z=ha.z*SQR(DeltaT)/MassA;

    InversionBendConstraintsDerivatives[CurrentSystem][j][1].x=hb.x*SQR(DeltaT)/MassB;
    InversionBendConstraintsDerivatives[CurrentSystem][j][1].y=hb.y*SQR(DeltaT)/MassB;
    InversionBendConstraintsDerivatives[CurrentSystem][j][1].z=hb.z*SQR(DeltaT)/MassB;

    InversionBendConstraintsDerivatives[CurrentSystem][j][2].x=hc.x*SQR(DeltaT)/MassC;
    InversionBendConstraintsDerivatives[CurrentSystem][j][2].y=hc.y*SQR(DeltaT)/MassC;
    InversionBendConstraintsDerivatives[CurrentSystem][j][2].z=hc.z*SQR(DeltaT)/MassC;

    InversionBendConstraintsDerivatives[CurrentSystem][j][3].x=hd.x*SQR(DeltaT)/MassD;
    InversionBendConstraintsDerivatives[CurrentSystem][j][3].y=hd.y*SQR(DeltaT)/MassD;
    InversionBendConstraintsDerivatives[CurrentSystem][j][3].z=hd.z*SQR(DeltaT)/MassD;
  }


  // store the initial dihedral-angle constraint gradients
  for(j=0;j<NumberOfDihedralConstraints[CurrentSystem];j++)
  {
    MassA=PseudoAtoms[DihedralConstraints[CurrentSystem][j][0]->Type].Mass;
    MassB=PseudoAtoms[DihedralConstraints[CurrentSystem][j][1]->Type].Mass;
    MassC=PseudoAtoms[DihedralConstraints[CurrentSystem][j][2]->Type].Mass;
    MassD=PseudoAtoms[DihedralConstraints[CurrentSystem][j][3]->Type].Mass;

    posA=DihedralConstraints[CurrentSystem][j][0]->Position;
    posB=DihedralConstraints[CurrentSystem][j][1]->Position;
    posC=DihedralConstraints[CurrentSystem][j][2]->Position;
    posD=DihedralConstraints[CurrentSystem][j][3]->Position;

    ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

    DihedralConstraintsDerivatives[CurrentSystem][j][0].x=ha.x*SQR(DeltaT)/MassA;
    DihedralConstraintsDerivatives[CurrentSystem][j][0].y=ha.y*SQR(DeltaT)/MassA;
    DihedralConstraintsDerivatives[CurrentSystem][j][0].z=ha.z*SQR(DeltaT)/MassA;

    DihedralConstraintsDerivatives[CurrentSystem][j][1].x=hb.x*SQR(DeltaT)/MassB;
    DihedralConstraintsDerivatives[CurrentSystem][j][1].y=hb.y*SQR(DeltaT)/MassB;
    DihedralConstraintsDerivatives[CurrentSystem][j][1].z=hb.z*SQR(DeltaT)/MassB;

    DihedralConstraintsDerivatives[CurrentSystem][j][2].x=hc.x*SQR(DeltaT)/MassC;
    DihedralConstraintsDerivatives[CurrentSystem][j][2].y=hc.y*SQR(DeltaT)/MassC;
    DihedralConstraintsDerivatives[CurrentSystem][j][2].z=hc.z*SQR(DeltaT)/MassC;

    DihedralConstraintsDerivatives[CurrentSystem][j][3].x=hd.x*SQR(DeltaT)/MassD;
    DihedralConstraintsDerivatives[CurrentSystem][j][3].y=hd.y*SQR(DeltaT)/MassD;
    DihedralConstraintsDerivatives[CurrentSystem][j][3].z=hd.z*SQR(DeltaT)/MassD;

  }

  // store the initial improper dihedral-angle constraint gradients
  for(j=0;j<NumberOfImproperDihedralConstraints[CurrentSystem];j++)
  {
    MassA=PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][j][0]->Type].Mass;
    MassB=PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][j][1]->Type].Mass;
    MassC=PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][j][2]->Type].Mass;
    MassD=PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][j][3]->Type].Mass;

    posA=ImproperDihedralConstraints[CurrentSystem][j][0]->Position;
    posB=ImproperDihedralConstraints[CurrentSystem][j][1]->Position;
    posC=ImproperDihedralConstraints[CurrentSystem][j][2]->Position;
    posD=ImproperDihedralConstraints[CurrentSystem][j][3]->Position;

    ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][0].x=ha.x*SQR(DeltaT)/MassA;
    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][0].y=ha.y*SQR(DeltaT)/MassA;
    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][0].z=ha.z*SQR(DeltaT)/MassA;

    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].x=hb.x*SQR(DeltaT)/MassB;
    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].y=hb.y*SQR(DeltaT)/MassB;
    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].z=hb.z*SQR(DeltaT)/MassB;

    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].x=hc.x*SQR(DeltaT)/MassC;
    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].y=hc.y*SQR(DeltaT)/MassC;
    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].z=hc.z*SQR(DeltaT)/MassC;

    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][3].x=hd.x*SQR(DeltaT)/MassD;
    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][3].y=hd.y*SQR(DeltaT)/MassD;
    ImproperDihedralConstraintsDerivatives[CurrentSystem][j][3].z=hd.z*SQR(DeltaT)/MassD;
  }


  int iter=0;
  do
  {
    iter++;
    max_error=0.0;

    if(iter==10000)
    {
      fprintf(stderr, "Shake in minimization failed, maximum number of iterations reached\n");
      exit(0);
    }

    // distance constraints in adsorbate molecules
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfBonds;j++)
      {
        if(Components[Type].BondType[j]==FIXED_BOND)
        {
          A=Components[Type].Bonds[j].A;
          B=Components[Type].Bonds[j].B;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;

          switch(DistanceConstraintType)
          {
            case DISTANCE_R:
              r0=Components[Type].BondArguments[j][0];
              break;
            default:
            case DISTANCE_R_SQUARED:
              r0=SQR(Components[Type].BondArguments[j][0]);
              break;
          }
          posA=Adsorbates[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          ReturnWilsonVectorsDistanceRATTLE(posA,posB,&Va,&Vb);
          Va.x*=SQR(DeltaT)*InverseMassA;  Vb.x*=SQR(DeltaT)*InverseMassB;
          Va.y*=SQR(DeltaT)*InverseMassA;  Vb.y*=SQR(DeltaT)*InverseMassB;
          Va.z*=SQR(DeltaT)*InverseMassA;  Vb.z*=SQR(DeltaT)*InverseMassB;
          Vba.x=Va.x-Vb.x;
          Vba.y=Va.y-Vb.y;
          Vba.z=Va.z-Vb.z;

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          Rba.x=posA.x-posB.x;
          Rba.y=posA.y-posB.y;
          Rba.z=posA.z-posB.z;

          r=ReturnConstraintDistance(posA,posB);

          gamma=(r-r0)/ReturnDistanceConstrainDerivative(Rba,Vba);

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x-=gamma*Va.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y-=gamma*Va.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z-=gamma*Va.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Position.x-=gamma*Vb.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.y-=gamma*Vb.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.z-=gamma*Vb.z;

          error=fabs(r-r0)/r0;
          max_error=MAX2(error,max_error);

/*
          switch(DistanceConstraintType)
          {
            case DISTANCE_R:
              fprintf(stderr, "Gamma %g distance: %g %g\n",gamma,r,r0);
              break;
            default:
            case DISTANCE_R_SQUARED:
              fprintf(stderr, "Gamma %g distance: %g %g\n",gamma,sqrt(r),sqrt(r0));
              break;
          }
*/
        }
      }
    }

    // distance constraints in cation molecules
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      Type=Cations[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfBonds;j++)
      {
        if(Components[Type].BondType[j]==FIXED_BOND)
        {
          A=Components[Type].Bonds[j].A;
          B=Components[Type].Bonds[j].B;

          InverseMassA=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[B].Type].Mass;

          switch(DistanceConstraintType)
          {
            case DISTANCE_R:
              r0=Components[Type].BondArguments[j][0];
              break;
            default:
            case DISTANCE_R_SQUARED:
              r0=SQR(Components[Type].BondArguments[j][0]);
              break;
          }
          posA=Cations[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Cations[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          ReturnWilsonVectorsDistanceRATTLE(posA,posB,&Va,&Vb);
          Va.x*=SQR(DeltaT)*InverseMassA;  Vb.x*=SQR(DeltaT)*InverseMassB;
          Va.y*=SQR(DeltaT)*InverseMassA;  Vb.y*=SQR(DeltaT)*InverseMassB;
          Va.z*=SQR(DeltaT)*InverseMassA;  Vb.z*=SQR(DeltaT)*InverseMassB;
          Vba.x=Va.x-Vb.x;
          Vba.y=Va.y-Vb.y;
          Vba.z=Va.z-Vb.z;

          posA=Cations[CurrentSystem][m].Atoms[A].Position;
          posB=Cations[CurrentSystem][m].Atoms[B].Position;
          Rba.x=posA.x-posB.x;
          Rba.y=posA.y-posB.y;
          Rba.z=posA.z-posB.z;

          r=ReturnConstraintDistance(posA,posB);

          gamma=(r-r0)/ReturnDistanceConstrainDerivative(Rba,Vba);

          Cations[CurrentSystem][m].Atoms[A].Position.x-=gamma*Va.x;
          Cations[CurrentSystem][m].Atoms[A].Position.y-=gamma*Va.y;
          Cations[CurrentSystem][m].Atoms[A].Position.z-=gamma*Va.z;

          Cations[CurrentSystem][m].Atoms[B].Position.x-=gamma*Vb.x;
          Cations[CurrentSystem][m].Atoms[B].Position.y-=gamma*Vb.y;
          Cations[CurrentSystem][m].Atoms[B].Position.z-=gamma*Vb.z;

          error=fabs(r-r0)/r0;
          max_error=MAX2(error,max_error);

/*
          switch(DistanceConstraintType)
          {
            case DISTANCE_R:
              fprintf(stderr, "Gamma %g distance: %g %g\n",gamma,r,r0);
              break;
            default:
            case DISTANCE_R_SQUARED:
              fprintf(stderr, "Gamma %g distance: %g %g\n",gamma,sqrt(r),sqrt(r0));
              break;
          }
*/
        }
      }
    }


    // bend-angle constraints in adsorbate molecules
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfBends;j++)
      {
        if(Components[Type].BendType[j]==FIXED_BEND)
        {
          A=Components[Type].Bends[j].A;
          B=Components[Type].Bends[j].B;
          C=Components[Type].Bends[j].C;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;

          switch(BendConstraintType)
          {
            case COS_THETA_SQUARED:
              Theta0=SQR(cos(Components[Type].BendArguments[j][0]));
              break;
            case COS_THETA:
              Theta0=cos(Components[Type].BendArguments[j][0]);
              break;
            case THETA:
            default:
              Theta0=Components[Type].BendArguments[j][0];
               break;
          }

          posA=Adsorbates[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          ReturnWilsonVectorsBendRATTLE(posA,posB,posC,&Va,&Vb,&Vc);
          Va.x*=SQR(DeltaT)*InverseMassA; Vb.x*=SQR(DeltaT)*InverseMassB; Vc.x*=SQR(DeltaT)*InverseMassC;
          Va.y*=SQR(DeltaT)*InverseMassA; Vb.y*=SQR(DeltaT)*InverseMassB; Vc.y*=SQR(DeltaT)*InverseMassC;
          Va.z*=SQR(DeltaT)*InverseMassA; Vb.z*=SQR(DeltaT)*InverseMassB; Vc.z*=SQR(DeltaT)*InverseMassC;
          Vba.x=Va.x-Vb.x; Vbc.x=Vc.x-Vb.x;
          Vba.y=Va.y-Vb.y; Vbc.y=Vc.y-Vb.y;
          Vba.z=Va.z-Vb.z; Vbc.z=Vc.z-Vb.z;

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

          Theta=ReturnConstraintBendAngle(posA,posB,posC);

          Rba.x=posA.x-posB.x;
          Rba.y=posA.y-posB.y;
          Rba.z=posA.z-posB.z;

          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;

          gamma=(Theta-Theta0)/ReturnAngleConstrainDerivative(Rba,Rbc,Vba,Vbc);

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x-=gamma*Va.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y-=gamma*Va.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z-=gamma*Va.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Position.x-=gamma*Vb.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.y-=gamma*Vb.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.z-=gamma*Vb.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Position.x-=gamma*Vc.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.y-=gamma*Vc.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.z-=gamma*Vc.z;

          error=fabs(Theta-Theta0);
          max_error=MAX2(error,max_error);

/*
          switch(BendConstraintType)
          {
            case COS_THETA_SQUARED:
              fprintf(stderr, "Gamma %g bend angle: %g %g\n",gamma,acos(sqrt(Theta))*RAD2DEG,acos(sqrt(Theta0))*RAD2DEG);
              break;
            case COS_THETA:
              fprintf(stderr, "Gamma %g bend angle: %g %g\n",gamma,acos(Theta)*RAD2DEG,acos(Theta0)*RAD2DEG);
              break;
            case THETA:
            default:
              fprintf(stderr, "Gamma %g bend angle: %g %g\n",gamma,(Theta)*RAD2DEG,(Theta0)*RAD2DEG);
              break;
          }
*/
        }
      }
    }

    // bend-angle constraints in cation molecules
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      Type=Cations[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfBends;j++)
      {
        if(Components[Type].BendType[j]==FIXED_BEND)
        {
          A=Components[Type].Bends[j].A;
          B=Components[Type].Bends[j].B;
          C=Components[Type].Bends[j].C;

          InverseMassA=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[C].Type].Mass;

          switch(BendConstraintType)
          {
            case COS_THETA_SQUARED:
              Theta0=SQR(cos(Components[Type].BendArguments[j][0]));
              break;
            case COS_THETA:
              Theta0=cos(Components[Type].BendArguments[j][0]);
              break;
            case THETA:
            default:
              Theta0=Components[Type].BendArguments[j][0];
               break;
          }

          posA=Cations[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Cations[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Cations[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          ReturnWilsonVectorsBendRATTLE(posA,posB,posC,&Va,&Vb,&Vc);
          Va.x*=SQR(DeltaT)*InverseMassA; Vb.x*=SQR(DeltaT)*InverseMassB; Vc.x*=SQR(DeltaT)*InverseMassC;
          Va.y*=SQR(DeltaT)*InverseMassA; Vb.y*=SQR(DeltaT)*InverseMassB; Vc.y*=SQR(DeltaT)*InverseMassC;
          Va.z*=SQR(DeltaT)*InverseMassA; Vb.z*=SQR(DeltaT)*InverseMassB; Vc.z*=SQR(DeltaT)*InverseMassC;
          Vba.x=Va.x-Vb.x; Vbc.x=Vc.x-Vb.x;
          Vba.y=Va.y-Vb.y; Vbc.y=Vc.y-Vb.y;
          Vba.z=Va.z-Vb.z; Vbc.z=Vc.z-Vb.z;

          posA=Cations[CurrentSystem][m].Atoms[A].Position;
          posB=Cations[CurrentSystem][m].Atoms[B].Position;
          posC=Cations[CurrentSystem][m].Atoms[C].Position;

          Theta=ReturnConstraintBendAngle(posA,posB,posC);

          Rba.x=posA.x-posB.x;
          Rba.y=posA.y-posB.y;
          Rba.z=posA.z-posB.z;

          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;

          gamma=(Theta-Theta0)/ReturnAngleConstrainDerivative(Rba,Rbc,Vba,Vbc);

          Cations[CurrentSystem][m].Atoms[A].Position.x-=gamma*Va.x;
          Cations[CurrentSystem][m].Atoms[A].Position.y-=gamma*Va.y;
          Cations[CurrentSystem][m].Atoms[A].Position.z-=gamma*Va.z;

          Cations[CurrentSystem][m].Atoms[B].Position.x-=gamma*Vb.x;
          Cations[CurrentSystem][m].Atoms[B].Position.y-=gamma*Vb.y;
          Cations[CurrentSystem][m].Atoms[B].Position.z-=gamma*Vb.z;

          Cations[CurrentSystem][m].Atoms[C].Position.x-=gamma*Vc.x;
          Cations[CurrentSystem][m].Atoms[C].Position.y-=gamma*Vc.y;
          Cations[CurrentSystem][m].Atoms[C].Position.z-=gamma*Vc.z;

          error=fabs(Theta-Theta0);
          max_error=MAX2(error,max_error);

/*
          switch(BendConstraintType)
          {
            case COS_THETA_SQUARED:
              fprintf(stderr, "Gamma %g bend angle: %g %g\n",gamma,acos(sqrt(Theta))*RAD2DEG,acos(sqrt(Theta0))*RAD2DEG);
              break;
            case COS_THETA:
              fprintf(stderr, "Gamma %g bend angle: %g %g\n",gamma,acos(Theta)*RAD2DEG,acos(Theta0)*RAD2DEG);
              break;
            case THETA:
            default:
              fprintf(stderr, "Gamma %g bend angle: %g %g\n",gamma,(Theta)*RAD2DEG,(Theta0)*RAD2DEG);
              break;
          }
*/
        }
      }
    }


    // dihedral-angle constraints in adsorbate molecules
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfTorsions;j++)
      {
        if(Components[Type].TorsionType[j]==FIXED_DIHEDRAL)
        {
          A=Components[Type].Torsions[j].A;
          B=Components[Type].Torsions[j].B;
          C=Components[Type].Torsions[j].C;
          D=Components[Type].Torsions[j].D;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[D].Type].Mass;

          switch(DihedralConstraintType)
          {
            case COS_PHI_SQUARED:
              Phi0=SQR(cos(Components[Type].TorsionArguments[j][0]));
              break;
            case COS_PHI:
              Phi0=cos(Components[Type].TorsionArguments[j][0]);
              break;
            case PHI:
            default:
              Phi0=Components[Type].TorsionArguments[j][0];
              break;
          }

          posA=Adsorbates[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].RattleReferencePosition;
          ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&Va,&Vb,&Vc,&Vd);
          Va.x*=SQR(DeltaT)*InverseMassA; Vb.x*=SQR(DeltaT)*InverseMassB; Vc.x*=SQR(DeltaT)*InverseMassC; Vd.x*=SQR(DeltaT)*InverseMassD;
          Va.y*=SQR(DeltaT)*InverseMassA; Vb.y*=SQR(DeltaT)*InverseMassB; Vc.y*=SQR(DeltaT)*InverseMassC; Vd.y*=SQR(DeltaT)*InverseMassD;
          Va.z*=SQR(DeltaT)*InverseMassA; Vb.z*=SQR(DeltaT)*InverseMassB; Vc.z*=SQR(DeltaT)*InverseMassC; Vd.z*=SQR(DeltaT)*InverseMassD;
          Vab.x=Vb.x-Va.x; Vbc.x=Vc.x-Vb.x; Vcd.x=Vd.x-Vc.x;
          Vab.y=Vb.y-Va.y; Vbc.y=Vc.y-Vb.y; Vcd.y=Vd.y-Vc.y;
          Vab.z=Vb.z-Va.z; Vbc.z=Vc.z-Vb.z; Vcd.z=Vd.z-Vc.z;

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;
          Phi=ReturnConstraintDihedralAngle(posA,posB,posC,posD);

          Rab.x=posB.x-posA.x;
          Rab.y=posB.y-posA.y;
          Rab.z=posB.z-posA.z;

          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;

          Rcd.x=posD.x-posC.x;
          Rcd.y=posD.y-posC.y;
          Rcd.z=posD.z-posC.z;

          gamma=(Phi-Phi0)/ReturnTorsionConstrainDerivative(Rab,Rbc,Rcd,Vab,Vbc,Vcd);

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x-=gamma*Va.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y-=gamma*Va.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z-=gamma*Va.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Position.x-=gamma*Vb.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.y-=gamma*Vb.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.z-=gamma*Vb.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Position.x-=gamma*Vc.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.y-=gamma*Vc.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.z-=gamma*Vc.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Position.x-=gamma*Vd.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.y-=gamma*Vd.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.z-=gamma*Vd.z;

          error=fabs(Phi-Phi0);
          max_error=MAX2(error,max_error);

          switch(DihedralConstraintType)
          {
            case COS_PHI_SQUARED:
              fprintf(stderr, "Cos^2(Phi): Gamma %g dihedral angle: %g %g\n",gamma,acos(sqrt(Phi))*RAD2DEG,acos(sqrt(Phi0))*RAD2DEG);
              break;
            case COS_PHI:
              fprintf(stderr, "Cos(Phi): Gamma %g dihedral angle: %g %g\n",gamma,acos(Phi)*RAD2DEG,acos(Phi0)*RAD2DEG);
              break;
            case PHI:
            default:
              fprintf(stderr, "Phi: Gamma %g dihedral angle: %g %g\n",gamma,(Phi)*RAD2DEG,(Phi0)*RAD2DEG);
              break;
          }
        }
      }
    }

    // dihedral-angle constraints in cation molecules
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      Type=Cations[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfTorsions;j++)
      {
        if(Components[Type].TorsionType[j]==FIXED_DIHEDRAL)
        {
          A=Components[Type].Torsions[j].A;
          B=Components[Type].Torsions[j].B;
          C=Components[Type].Torsions[j].C;
          D=Components[Type].Torsions[j].D;

          InverseMassA=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[D].Type].Mass;

          switch(DihedralConstraintType)
          {
            case COS_PHI_SQUARED:
              Phi0=SQR(cos(Components[Type].TorsionArguments[j][0]));
              break;
            case COS_PHI:
              Phi0=cos(Components[Type].TorsionArguments[j][0]);
              break;
            case PHI:
            default:
              Phi0=Components[Type].TorsionArguments[j][0];
              break;
          }

          posA=Cations[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Cations[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Cations[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          posD=Cations[CurrentSystem][m].Atoms[D].RattleReferencePosition;
          ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&Va,&Vb,&Vc,&Vd);
          Va.x*=SQR(DeltaT)*InverseMassA; Vb.x*=SQR(DeltaT)*InverseMassB; Vc.x*=SQR(DeltaT)*InverseMassC; Vd.x*=SQR(DeltaT)*InverseMassD;
          Va.y*=SQR(DeltaT)*InverseMassA; Vb.y*=SQR(DeltaT)*InverseMassB; Vc.y*=SQR(DeltaT)*InverseMassC; Vd.y*=SQR(DeltaT)*InverseMassD;
          Va.z*=SQR(DeltaT)*InverseMassA; Vb.z*=SQR(DeltaT)*InverseMassB; Vc.z*=SQR(DeltaT)*InverseMassC; Vd.z*=SQR(DeltaT)*InverseMassD;
          Vab.x=Vb.x-Va.x; Vbc.x=Vc.x-Vb.x; Vcd.x=Vd.x-Vc.x;
          Vab.y=Vb.y-Va.y; Vbc.y=Vc.y-Vb.y; Vcd.y=Vd.y-Vc.y;
          Vab.z=Vb.z-Va.z; Vbc.z=Vc.z-Vb.z; Vcd.z=Vd.z-Vc.z;

          posA=Cations[CurrentSystem][m].Atoms[A].Position;
          posB=Cations[CurrentSystem][m].Atoms[B].Position;
          posC=Cations[CurrentSystem][m].Atoms[C].Position;
          posD=Cations[CurrentSystem][m].Atoms[D].Position;
          Phi=ReturnConstraintDihedralAngle(posA,posB,posC,posD);

          Rab.x=posB.x-posA.x;
          Rab.y=posB.y-posA.y;
          Rab.z=posB.z-posA.z;

          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;

          Rcd.x=posD.x-posC.x;
          Rcd.y=posD.y-posC.y;
          Rcd.z=posD.z-posC.z;

          gamma=(Phi-Phi0)/ReturnTorsionConstrainDerivative(Rab,Rbc,Rcd,Vab,Vbc,Vcd);

          Cations[CurrentSystem][m].Atoms[A].Position.x-=gamma*Va.x;
          Cations[CurrentSystem][m].Atoms[A].Position.y-=gamma*Va.y;
          Cations[CurrentSystem][m].Atoms[A].Position.z-=gamma*Va.z;

          Cations[CurrentSystem][m].Atoms[B].Position.x-=gamma*Vb.x;
          Cations[CurrentSystem][m].Atoms[B].Position.y-=gamma*Vb.y;
          Cations[CurrentSystem][m].Atoms[B].Position.z-=gamma*Vb.z;

          Cations[CurrentSystem][m].Atoms[C].Position.x-=gamma*Vc.x;
          Cations[CurrentSystem][m].Atoms[C].Position.y-=gamma*Vc.y;
          Cations[CurrentSystem][m].Atoms[C].Position.z-=gamma*Vc.z;

          Cations[CurrentSystem][m].Atoms[D].Position.x-=gamma*Vd.x;
          Cations[CurrentSystem][m].Atoms[D].Position.y-=gamma*Vd.y;
          Cations[CurrentSystem][m].Atoms[D].Position.z-=gamma*Vd.z;

          error=fabs(Phi-Phi0);
          max_error=MAX2(error,max_error);

/*
          switch(DihedralConstraintType)
          {
            case COS_PHI_SQUARED:
              fprintf(stderr, "Gamma %g dihedral angle: %g %g\n",gamma,acos(sqrt(Phi))*RAD2DEG,acos(sqrt(Phi0))*RAD2DEG);
              break;
            case COS_PHI:
              fprintf(stderr, "Gamma %g dihedral angle: %g %g\n",gamma,acos(Phi)*RAD2DEG,acos(Phi0)*RAD2DEG);
              break;
            case PHI:
            default:
              fprintf(stderr, "Gamma %g dihedral angle: %g %g\n",gamma,(Phi)*RAD2DEG,(Phi0)*RAD2DEG);
              break;
          }
*/
        }
      }
    }


    // improper dihedral-angle constraints in adsorbate molecules
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfImproperTorsions;j++)
      {
        if(Components[Type].ImproperTorsionType[j]==FIXED_IMPROPER_DIHEDRAL)
        {
          A=Components[Type].ImproperTorsions[j].A;
          B=Components[Type].ImproperTorsions[j].B;
          C=Components[Type].ImproperTorsions[j].C;
          D=Components[Type].ImproperTorsions[j].D;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[D].Type].Mass;

          switch(DihedralConstraintType)
          {
            case COS_PHI_SQUARED:
              Phi0=SQR(cos(Components[Type].ImproperTorsionArguments[j][0]));
              break;
            case COS_PHI:
              Phi0=cos(Components[Type].ImproperTorsionArguments[j][0]);
              break;
            case PHI:
            default:
              Phi0=Components[Type].ImproperTorsionArguments[j][0];
              break;
          }

          posA=Adsorbates[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].RattleReferencePosition;
          ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&Va,&Vb,&Vc,&Vd);
          Va.x*=SQR(DeltaT)*InverseMassA; Vb.x*=SQR(DeltaT)*InverseMassB; Vc.x*=SQR(DeltaT)*InverseMassC; Vd.x*=SQR(DeltaT)*InverseMassD;
          Va.y*=SQR(DeltaT)*InverseMassA; Vb.y*=SQR(DeltaT)*InverseMassB; Vc.y*=SQR(DeltaT)*InverseMassC; Vd.y*=SQR(DeltaT)*InverseMassD;
          Va.z*=SQR(DeltaT)*InverseMassA; Vb.z*=SQR(DeltaT)*InverseMassB; Vc.z*=SQR(DeltaT)*InverseMassC; Vd.z*=SQR(DeltaT)*InverseMassD;
          Vab.x=Vb.x-Va.x; Vbc.x=Vc.x-Vb.x; Vcd.x=Vd.x-Vc.x;
          Vab.y=Vb.y-Va.y; Vbc.y=Vc.y-Vb.y; Vcd.y=Vd.y-Vc.y;
          Vab.z=Vb.z-Va.z; Vbc.z=Vc.z-Vb.z; Vcd.z=Vd.z-Vc.z;

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;
          Phi=ReturnConstraintDihedralAngle(posA,posB,posC,posD);

          Rab.x=posB.x-posA.x;
          Rab.y=posB.y-posA.y;
          Rab.z=posB.z-posA.z;

          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;

          Rcd.x=posD.x-posC.x;
          Rcd.y=posD.y-posC.y;
          Rcd.z=posD.z-posC.z;

          gamma=(Phi-Phi0)/ReturnTorsionConstrainDerivative(Rab,Rbc,Rcd,Vab,Vbc,Vcd);

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x-=gamma*Va.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y-=gamma*Va.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z-=gamma*Va.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Position.x-=gamma*Vb.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.y-=gamma*Vb.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.z-=gamma*Vb.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Position.x-=gamma*Vc.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.y-=gamma*Vc.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.z-=gamma*Vc.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Position.x-=gamma*Vd.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.y-=gamma*Vd.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.z-=gamma*Vd.z;

          error=fabs(Phi-Phi0);
          max_error=MAX2(error,max_error);

/*
          switch(DihedralConstraintType)
          {
            case COS_PHI_SQUARED:
              fprintf(stderr, "Gamma %g improper dihedral angle: %g %g\n",gamma,acos(sqrt(Phi))*RAD2DEG,acos(sqrt(Phi0))*RAD2DEG);
              break;
            case COS_PHI:
              fprintf(stderr, "Gamma %g improper dihedral angle: %g %g\n",gamma,acos(Phi)*RAD2DEG,acos(Phi0)*RAD2DEG);
              break;
            case PHI:
            default:
              fprintf(stderr, "Gamma %g improper dihedral angle: %g %g\n",gamma,(Phi)*RAD2DEG,(Phi0)*RAD2DEG);
              break;
          }
*/
        }
      }
    }

    // improper dihedral-angle constraints in cation molecules
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      Type=Cations[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfImproperTorsions;j++)
      {
        if(Components[Type].ImproperTorsionType[j]==FIXED_IMPROPER_DIHEDRAL)
        {
          A=Components[Type].ImproperTorsions[j].A;
          B=Components[Type].ImproperTorsions[j].B;
          C=Components[Type].ImproperTorsions[j].C;
          D=Components[Type].ImproperTorsions[j].D;

          InverseMassA=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[D].Type].Mass;

          switch(DihedralConstraintType)
          {
            case COS_PHI_SQUARED:
              Phi0=SQR(cos(Components[Type].ImproperTorsionArguments[j][0]));
              break;
            case COS_PHI:
              Phi0=cos(Components[Type].ImproperTorsionArguments[j][0]);
              break;
            case PHI:
            default:
              Phi0=Components[Type].ImproperTorsionArguments[j][0];
              break;
          }

          posA=Cations[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Cations[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Cations[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          posD=Cations[CurrentSystem][m].Atoms[D].RattleReferencePosition;
          ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&Va,&Vb,&Vc,&Vd);
          Va.x*=SQR(DeltaT)*InverseMassA; Vb.x*=SQR(DeltaT)*InverseMassB; Vc.x*=SQR(DeltaT)*InverseMassC; Vd.x*=SQR(DeltaT)*InverseMassD;
          Va.y*=SQR(DeltaT)*InverseMassA; Vb.y*=SQR(DeltaT)*InverseMassB; Vc.y*=SQR(DeltaT)*InverseMassC; Vd.y*=SQR(DeltaT)*InverseMassD;
          Va.z*=SQR(DeltaT)*InverseMassA; Vb.z*=SQR(DeltaT)*InverseMassB; Vc.z*=SQR(DeltaT)*InverseMassC; Vd.z*=SQR(DeltaT)*InverseMassD;
          Vab.x=Vb.x-Va.x; Vbc.x=Vc.x-Vb.x; Vcd.x=Vd.x-Vc.x;
          Vab.y=Vb.y-Va.y; Vbc.y=Vc.y-Vb.y; Vcd.y=Vd.y-Vc.y;
          Vab.z=Vb.z-Va.z; Vbc.z=Vc.z-Vb.z; Vcd.z=Vd.z-Vc.z;

          posA=Cations[CurrentSystem][m].Atoms[A].Position;
          posB=Cations[CurrentSystem][m].Atoms[B].Position;
          posC=Cations[CurrentSystem][m].Atoms[C].Position;
          posD=Cations[CurrentSystem][m].Atoms[D].Position;
          Phi=ReturnConstraintDihedralAngle(posA,posB,posC,posD);

          Rab.x=posB.x-posA.x;
          Rab.y=posB.y-posA.y;
          Rab.z=posB.z-posA.z;

          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;

          Rcd.x=posD.x-posC.x;
          Rcd.y=posD.y-posC.y;
          Rcd.z=posD.z-posC.z;

          gamma=(Phi-Phi0)/ReturnTorsionConstrainDerivative(Rab,Rbc,Rcd,Vab,Vbc,Vcd);

          Cations[CurrentSystem][m].Atoms[A].Position.x-=gamma*Va.x;
          Cations[CurrentSystem][m].Atoms[A].Position.y-=gamma*Va.y;
          Cations[CurrentSystem][m].Atoms[A].Position.z-=gamma*Va.z;

          Cations[CurrentSystem][m].Atoms[B].Position.x-=gamma*Vb.x;
          Cations[CurrentSystem][m].Atoms[B].Position.y-=gamma*Vb.y;
          Cations[CurrentSystem][m].Atoms[B].Position.z-=gamma*Vb.z;

          Cations[CurrentSystem][m].Atoms[C].Position.x-=gamma*Vc.x;
          Cations[CurrentSystem][m].Atoms[C].Position.y-=gamma*Vc.y;
          Cations[CurrentSystem][m].Atoms[C].Position.z-=gamma*Vc.z;

          Cations[CurrentSystem][m].Atoms[D].Position.x-=gamma*Vd.x;
          Cations[CurrentSystem][m].Atoms[D].Position.y-=gamma*Vd.y;
          Cations[CurrentSystem][m].Atoms[D].Position.z-=gamma*Vd.z;

          error=fabs(Phi-Phi0);
          max_error=MAX2(error,max_error);

/*
          switch(DihedralConstraintType)
          {
            case COS_PHI_SQUARED:
              fprintf(stderr, "Gamma %g improper dihedral angle: %g %g\n",gamma,acos(sqrt(Phi))*RAD2DEG,acos(sqrt(Phi0))*RAD2DEG);
              break;
            case COS_PHI:
              fprintf(stderr, "Gamma %g improper dihedral angle: %g %g\n",gamma,acos(Phi)*RAD2DEG,acos(Phi0)*RAD2DEG);
              break;
            case PHI:
            default:
              fprintf(stderr, "Gamma %g improper dihedral angle: %g %g\n",gamma,(Phi)*RAD2DEG,(Phi0)*RAD2DEG);
              break;
          }
*/
        }
      }
    }


    // inversion-bend-angle constraints in adsorbate molecules
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfInversionBends;j++)
      {
        if(Components[Type].InversionBendType[j]==FIXED_INVERSION_BEND)
        {
          A=Components[Type].InversionBends[j].A;
          B=Components[Type].InversionBends[j].B;
          C=Components[Type].InversionBends[j].C;
          D=Components[Type].InversionBends[j].D;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[D].Type].Mass;

          switch(InversionBendConstraintType)
          {
            case SIN_CHI_SQUARED:
              Chi0=SQR(sin(Components[Type].InversionBendArguments[j][0]));
              break;
            case SIN_CHI:
              Chi0=sin(Components[Type].InversionBendArguments[j][0]);
              break;
            case CHI:
            default:
              Chi0=Components[Type].InversionBendArguments[j][0];
              break;
          }

          posA=Adsorbates[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].RattleReferencePosition;
          ReturnWilsonVectorsInversionBendRATTLE(posA,posB,posC,posD,&Va,&Vb,&Vc,&Vd);
          Va.x*=SQR(DeltaT)*InverseMassA; Vb.x*=SQR(DeltaT)*InverseMassB; Vc.x*=SQR(DeltaT)*InverseMassC; Vd.x*=SQR(DeltaT)*InverseMassD;
          Va.y*=SQR(DeltaT)*InverseMassA; Vb.y*=SQR(DeltaT)*InverseMassB; Vc.y*=SQR(DeltaT)*InverseMassC; Vd.y*=SQR(DeltaT)*InverseMassD;
          Va.z*=SQR(DeltaT)*InverseMassA; Vb.z*=SQR(DeltaT)*InverseMassB; Vc.z*=SQR(DeltaT)*InverseMassC; Vd.z*=SQR(DeltaT)*InverseMassD;
          Vba.x=(Va.x-Vb.x); Vbc.x=(Vc.x-Vb.x); Vbd.x=(Vd.x-Vb.x);
          Vba.y=(Va.y-Vb.y); Vbc.y=(Vc.y-Vb.y); Vbd.y=(Vd.y-Vb.y);
          Vba.z=(Va.z-Vb.z); Vbc.z=(Vc.z-Vb.z); Vbd.z=(Vd.z-Vb.z);

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;
          Chi=ReturnConstraintInversionBendAngle(posA,posB,posC,posD);

          Rba.x=(posA.x-posB.x);
          Rba.y=(posA.y-posB.y);
          Rba.z=(posA.z-posB.z);

          Rbc.x=(posC.x-posB.x);
          Rbc.y=(posC.y-posB.y);
          Rbc.z=(posC.z-posB.z);

          Rbd.x=(posD.x-posB.x);
          Rbd.y=(posD.y-posB.y);
          Rbd.z=(posD.z-posB.z);

          gamma=(Chi-Chi0)/(ReturnInversionBendConstrainDerivative(Rba,Rbc,Rbd,Vba,Vbc,Vbc));

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x-=gamma*Va.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y-=gamma*Va.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z-=gamma*Va.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Position.x-=gamma*Vb.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.y-=gamma*Vb.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.z-=gamma*Vb.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Position.x-=gamma*Vc.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.y-=gamma*Vc.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.z-=gamma*Vc.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Position.x-=gamma*Vd.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.y-=gamma*Vd.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.z-=gamma*Vd.z;


          error=fabs(Chi-Chi0);
          max_error=MAX2(error,max_error);

          switch(InversionBendConstraintType)
          {
            case SIN_CHI_SQUARED:
              fprintf(stderr, "Sin(Chi)^2: Gamma %g inversion-bend angle: %g %g\n",gamma,asin(sqrt(Chi))*RAD2DEG,asin(sqrt(Chi0))*RAD2DEG);
              break;
            case SIN_CHI:
              fprintf(stderr, "Sin(Chi): Gamma %g inversion-bend angle: %g %g\n",gamma,asin(Chi)*RAD2DEG,asin(Chi0)*RAD2DEG);
              break;
            case CHI:
            default:
              fprintf(stderr, "Chi: Gamma %g inversion-bend angle: %g %g\n",gamma,(Chi)*RAD2DEG,(Chi0)*RAD2DEG);
              break;
          }
        }
      }
    }

    // inversion-bend-angle constraints in cations molecules
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      Type=Cations[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfInversionBends;j++)
      {
        if(Components[Type].InversionBendType[j]==FIXED_INVERSION_BEND)
        {
          A=Components[Type].InversionBends[j].A;
          B=Components[Type].InversionBends[j].B;
          C=Components[Type].InversionBends[j].C;
          D=Components[Type].InversionBends[j].D;

          InverseMassA=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Cations[CurrentSystem][m].Atoms[D].Type].Mass;

          switch(InversionBendConstraintType)
          {
            case SIN_CHI_SQUARED:
              Chi0=SQR(cos(Components[Type].InversionBendArguments[j][0]));
              break;
            case SIN_CHI:
              Chi0=cos(Components[Type].InversionBendArguments[j][0]);
              break;
            case CHI:
            default:
              Chi0=Components[Type].InversionBendArguments[j][0];
              break;
          }

          posA=Cations[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Cations[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Cations[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          posD=Cations[CurrentSystem][m].Atoms[D].RattleReferencePosition;
          ReturnWilsonVectorsInversionBendRATTLE(posA,posB,posC,posD,&Va,&Vb,&Vc,&Vd);
          Va.x*=SQR(DeltaT)*InverseMassA; Vb.x*=SQR(DeltaT)*InverseMassB; Vc.x*=SQR(DeltaT)*InverseMassC; Vd.x*=SQR(DeltaT)*InverseMassD;
          Va.y*=SQR(DeltaT)*InverseMassA; Vb.y*=SQR(DeltaT)*InverseMassB; Vc.y*=SQR(DeltaT)*InverseMassC; Vd.y*=SQR(DeltaT)*InverseMassD;
          Va.z*=SQR(DeltaT)*InverseMassA; Vb.z*=SQR(DeltaT)*InverseMassB; Vc.z*=SQR(DeltaT)*InverseMassC; Vd.z*=SQR(DeltaT)*InverseMassD;
          Vba.x=Va.x-Vb.x; Vbc.x=Vc.x-Vb.x; Vbd.x=Vd.x-Vb.x;
          Vba.y=Va.y-Vb.y; Vbc.y=Vc.y-Vb.y; Vbd.y=Vd.y-Vb.y;
          Vba.z=Va.z-Vb.z; Vbc.z=Vc.z-Vb.z; Vbd.z=Vd.z-Vb.z;

          posA=Cations[CurrentSystem][m].Atoms[A].Position;
          posB=Cations[CurrentSystem][m].Atoms[B].Position;
          posC=Cations[CurrentSystem][m].Atoms[C].Position;
          posD=Cations[CurrentSystem][m].Atoms[D].Position;
          Chi=ReturnConstraintInversionBendAngle(posA,posB,posC,posD);

          Rba.x=posA.x-posB.x;
          Rba.y=posA.y-posB.y;
          Rba.z=posA.z-posB.z;

          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;

          Rbd.x=posD.x-posB.x;
          Rbd.y=posD.y-posB.y;
          Rbd.z=posD.z-posB.z;

          gamma=(Chi-Chi0)/(ReturnInversionBendConstrainDerivative(Rba,Rbc,Rbd,Vba,Vbc,Vbd));

          Cations[CurrentSystem][m].Atoms[A].Position.x-=gamma*Va.x;
          Cations[CurrentSystem][m].Atoms[A].Position.y-=gamma*Va.y;
          Cations[CurrentSystem][m].Atoms[A].Position.z-=gamma*Va.z;

          Cations[CurrentSystem][m].Atoms[B].Position.x-=gamma*Vb.x;
          Cations[CurrentSystem][m].Atoms[B].Position.y-=gamma*Vb.y;
          Cations[CurrentSystem][m].Atoms[B].Position.z-=gamma*Vb.z;

          Cations[CurrentSystem][m].Atoms[C].Position.x-=gamma*Vc.x;
          Cations[CurrentSystem][m].Atoms[C].Position.y-=gamma*Vc.y;
          Cations[CurrentSystem][m].Atoms[C].Position.z-=gamma*Vc.z;

          Cations[CurrentSystem][m].Atoms[D].Position.x-=gamma*Vd.x;
          Cations[CurrentSystem][m].Atoms[D].Position.y-=gamma*Vd.y;
          Cations[CurrentSystem][m].Atoms[D].Position.z-=gamma*Vd.z;

          error=fabs(Chi-Chi0);
          max_error=MAX2(error,max_error);

/*
          switch(InversionBendConstraintType)
          {
            case SIN_CHI_SQUARED:
              fprintf(stderr, "Sin(Chi)^2: Gamma %g inversion-bend angle: %g %g\n",gamma,asin(sqrt(Chi))*RAD2DEG,asin(sqrt(Chi0))*RAD2DEG);
              break;
            case SIN_CHI:
              fprintf(stderr, "Sin(Chi): Gamma %g inversion-bend angle: %g %g\n",gamma,asin(Chi)*RAD2DEG,asin(Chi0)*RAD2DEG);
              break;
            case CHI:
            default:
              fprintf(stderr, "Chi: Gamma %g inversion-bend angle: %g %g\n",gamma,(Chi)*RAD2DEG,(Chi0)*RAD2DEG);
              break;
          }
*/
        }
      }
    }



    for(j=0;j<NumberOfDistanceConstraints[CurrentSystem];j++)
    {
      switch(DistanceConstraintType)
      {
        case DISTANCE_R:
          r0=DistanceConstraintParameter[CurrentSystem][j];
          break;
        default:
        case DISTANCE_R_SQUARED:
          r0=SQR(DistanceConstraintParameter[CurrentSystem][j]);
          break;
      }

      Vba.x=DistanceConstraintsDerivatives[CurrentSystem][j][0].x-DistanceConstraintsDerivatives[CurrentSystem][j][1].x;
      Vba.y=DistanceConstraintsDerivatives[CurrentSystem][j][0].y-DistanceConstraintsDerivatives[CurrentSystem][j][1].y;
      Vba.z=DistanceConstraintsDerivatives[CurrentSystem][j][0].z-DistanceConstraintsDerivatives[CurrentSystem][j][1].z;

      posA=DistanceConstraints[CurrentSystem][j][0]->Position;
      posB=DistanceConstraints[CurrentSystem][j][1]->Position;

      Rba.x=posA.x-posB.x;
      Rba.y=posA.y-posB.y;
      Rba.z=posA.z-posB.z;
      Rba=ApplyBoundaryCondition(Rba);

      r=ReturnConstraintDistance(posA,posB);

      gamma=(r-r0)/ReturnDistanceConstrainDerivative(Rba,Vba);

      DistanceConstraints[CurrentSystem][j][0]->Position.x-=gamma*DistanceConstraintsDerivatives[CurrentSystem][j][0].x;
      DistanceConstraints[CurrentSystem][j][0]->Position.y-=gamma*DistanceConstraintsDerivatives[CurrentSystem][j][0].y;
      DistanceConstraints[CurrentSystem][j][0]->Position.z-=gamma*DistanceConstraintsDerivatives[CurrentSystem][j][0].z;

      DistanceConstraints[CurrentSystem][j][1]->Position.x-=gamma*DistanceConstraintsDerivatives[CurrentSystem][j][1].x;
      DistanceConstraints[CurrentSystem][j][1]->Position.y-=gamma*DistanceConstraintsDerivatives[CurrentSystem][j][1].y;
      DistanceConstraints[CurrentSystem][j][1]->Position.z-=gamma*DistanceConstraintsDerivatives[CurrentSystem][j][1].z;

      error=fabs(r-r0)/r0;
      max_error=MAX2(error,max_error);

/*
      switch(DistanceConstraintType)
      {
        case DISTANCE_R:
          fprintf(stderr, "Gamma %g distance: %g %g\n",gamma,r,r0);
          break;
        default:
        case DISTANCE_R_SQUARED:
          fprintf(stderr, "Gamma %g distance: %g %g\n",gamma,sqrt(r),sqrt(r0));
          break;
      }
*/
    }
    for(j=0;j<NumberOfAngleConstraints[CurrentSystem];j++)
    {
      switch(BendConstraintType)
      {
        case COS_THETA_SQUARED:
          Theta0=SQR(cos(AngleConstraintParameter[CurrentSystem][j]));
          break;
        case COS_THETA:
          Theta0=cos(AngleConstraintParameter[CurrentSystem][j]);
          break;
        case THETA:
        default:
          Theta0=AngleConstraintParameter[CurrentSystem][j];
          break;
      }

      posA=AngleConstraints[CurrentSystem][j][0]->Position;
      posB=AngleConstraints[CurrentSystem][j][1]->Position;
      posC=AngleConstraints[CurrentSystem][j][2]->Position;

      Theta=ReturnConstraintBendAngle(posA,posB,posC);

      Vba.x=AngleConstraintsDerivatives[CurrentSystem][j][0].x-AngleConstraintsDerivatives[CurrentSystem][j][1].x;
      Vba.y=AngleConstraintsDerivatives[CurrentSystem][j][0].y-AngleConstraintsDerivatives[CurrentSystem][j][1].y;
      Vba.z=AngleConstraintsDerivatives[CurrentSystem][j][0].z-AngleConstraintsDerivatives[CurrentSystem][j][1].z;

      Vbc.x=AngleConstraintsDerivatives[CurrentSystem][j][2].x-AngleConstraintsDerivatives[CurrentSystem][j][1].x;
      Vbc.y=AngleConstraintsDerivatives[CurrentSystem][j][2].y-AngleConstraintsDerivatives[CurrentSystem][j][1].y;
      Vbc.z=AngleConstraintsDerivatives[CurrentSystem][j][2].z-AngleConstraintsDerivatives[CurrentSystem][j][1].z;

      Rba.x=posA.x-posB.x;
      Rba.y=posA.y-posB.y;
      Rba.z=posA.z-posB.z;
      Rba=ApplyBoundaryCondition(Rba);

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      Rbc=ApplyBoundaryCondition(Rbc);

      gamma=(Theta-Theta0)/ReturnAngleConstrainDerivative(Rba,Rbc,Vba,Vbc);

      AngleConstraints[CurrentSystem][j][0]->Position.x-=gamma*AngleConstraintsDerivatives[CurrentSystem][j][0].x;
      AngleConstraints[CurrentSystem][j][0]->Position.y-=gamma*AngleConstraintsDerivatives[CurrentSystem][j][0].y;
      AngleConstraints[CurrentSystem][j][0]->Position.z-=gamma*AngleConstraintsDerivatives[CurrentSystem][j][0].z;

      AngleConstraints[CurrentSystem][j][1]->Position.x-=gamma*AngleConstraintsDerivatives[CurrentSystem][j][1].x;
      AngleConstraints[CurrentSystem][j][1]->Position.y-=gamma*AngleConstraintsDerivatives[CurrentSystem][j][1].y;
      AngleConstraints[CurrentSystem][j][1]->Position.z-=gamma*AngleConstraintsDerivatives[CurrentSystem][j][1].z;

      AngleConstraints[CurrentSystem][j][2]->Position.x-=gamma*AngleConstraintsDerivatives[CurrentSystem][j][2].x;
      AngleConstraints[CurrentSystem][j][2]->Position.y-=gamma*AngleConstraintsDerivatives[CurrentSystem][j][2].y;
      AngleConstraints[CurrentSystem][j][2]->Position.z-=gamma*AngleConstraintsDerivatives[CurrentSystem][j][2].z;

      error=fabs(Theta-Theta0);
      max_error=MAX2(error,max_error);

/*
      switch(BendConstraintType)
      {
        case COS_THETA_SQUARED:
          fprintf(stderr, "Gamma %g bend angle: %g %g\n",gamma,acos(sqrt(Theta))*RAD2DEG,acos(sqrt(Theta0))*RAD2DEG);
          break;
        case COS_THETA:
          fprintf(stderr, "Gamma %g bend angle: %g %g\n",gamma,acos(Theta)*RAD2DEG,acos(Theta0)*RAD2DEG);
          break;
        case THETA:
        default:
          fprintf(stderr, "Gamma %g bend angle: %g %g\n",gamma,(Theta)*RAD2DEG,(Theta0)*RAD2DEG);
          break;
      }
*/
    }

    for(j=0;j<NumberOfInversionBendConstraints[CurrentSystem];j++)
    {
      switch(InversionBendConstraintType)
      {
        case SIN_CHI_SQUARED:
          Chi0=SQR(cos(InversionBendConstraintParameter[CurrentSystem][j]));
          break;
        case SIN_CHI:
          Chi0=cos(InversionBendConstraintParameter[CurrentSystem][j]);
          break;
        case CHI:
        default:
          Chi0=InversionBendConstraintParameter[CurrentSystem][j];
          break;
      }

      posA=InversionBendConstraints[CurrentSystem][j][0]->Position;
      posB=InversionBendConstraints[CurrentSystem][j][1]->Position;
      posC=InversionBendConstraints[CurrentSystem][j][2]->Position;
      posD=InversionBendConstraints[CurrentSystem][j][3]->Position;

      Vba.x=InversionBendConstraintsDerivatives[CurrentSystem][j][0].x-InversionBendConstraintsDerivatives[CurrentSystem][j][1].x;
      Vba.y=InversionBendConstraintsDerivatives[CurrentSystem][j][0].y-InversionBendConstraintsDerivatives[CurrentSystem][j][1].y;
      Vba.z=InversionBendConstraintsDerivatives[CurrentSystem][j][0].z-InversionBendConstraintsDerivatives[CurrentSystem][j][1].z;

      Vbc.x=InversionBendConstraintsDerivatives[CurrentSystem][j][2].x-InversionBendConstraintsDerivatives[CurrentSystem][j][1].x;
      Vbc.y=InversionBendConstraintsDerivatives[CurrentSystem][j][2].y-InversionBendConstraintsDerivatives[CurrentSystem][j][1].y;
      Vbc.z=InversionBendConstraintsDerivatives[CurrentSystem][j][2].z-InversionBendConstraintsDerivatives[CurrentSystem][j][1].z;

      Vbd.x=InversionBendConstraintsDerivatives[CurrentSystem][j][3].x-InversionBendConstraintsDerivatives[CurrentSystem][j][1].x;
      Vbd.y=InversionBendConstraintsDerivatives[CurrentSystem][j][3].y-InversionBendConstraintsDerivatives[CurrentSystem][j][1].y;
      Vbd.z=InversionBendConstraintsDerivatives[CurrentSystem][j][3].z-InversionBendConstraintsDerivatives[CurrentSystem][j][1].z;

      Chi=ReturnConstraintInversionBendAngle(posA,posB,posC,posD);

      Rba.x=posA.x-posB.x;
      Rba.y=posA.y-posB.y;
      Rba.z=posA.z-posB.z;
      Rba=ApplyBoundaryCondition(Rba);

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      Rbc=ApplyBoundaryCondition(Rbc);

      Rbd.x=posD.x-posB.x;
      Rbd.y=posD.y-posB.y;
      Rbd.z=posD.z-posB.z;
      Rbd=ApplyBoundaryCondition(Rbd);

      gamma=(Chi-Chi0)/(ReturnInversionBendConstrainDerivative(Rba,Rbc,Rbd,Vba,Vbc,Vbd));

      InversionBendConstraints[CurrentSystem][j][0]->Position.x-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][0].x;
      InversionBendConstraints[CurrentSystem][j][0]->Position.y-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][0].y;
      InversionBendConstraints[CurrentSystem][j][0]->Position.z-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][0].z;

      InversionBendConstraints[CurrentSystem][j][1]->Position.x-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][1].x;
      InversionBendConstraints[CurrentSystem][j][1]->Position.y-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][1].y;
      InversionBendConstraints[CurrentSystem][j][1]->Position.z-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][1].z;

      InversionBendConstraints[CurrentSystem][j][2]->Position.x-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][2].x;
      InversionBendConstraints[CurrentSystem][j][2]->Position.y-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][2].y;
      InversionBendConstraints[CurrentSystem][j][2]->Position.z-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][2].z;

      InversionBendConstraints[CurrentSystem][j][3]->Position.x-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][3].x;
      InversionBendConstraints[CurrentSystem][j][3]->Position.y-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][3].y;
      InversionBendConstraints[CurrentSystem][j][3]->Position.z-=gamma*InversionBendConstraintsDerivatives[CurrentSystem][j][3].z;

      error=fabs(Chi-Chi0);
      max_error=MAX2(error,max_error);

/*
      switch(InversionBendConstraintType)
      {
        case COS_CHI_SQUARED:
          fprintf(stderr, "Gamma %g inversion-bend angle: %g %g\n",gamma,acos(sqrt(Chi))*RAD2DEG,acos(sqrt(Chi0))*RAD2DEG);
          break;
        case COS_CHI:
          fprintf(stderr, "Gamma %g inversion-bend angle: %g %g\n",gamma,acos(Chi)*RAD2DEG,acos(Chi0)*RAD2DEG);
          break;
        case CHI:
        default:
          fprintf(stderr, "Gamma %g inversion-bend angle: %g %g\n",gamma,(Chi)*RAD2DEG,(Chi0)*RAD2DEG);
          break;
      }
*/
    }

    for(j=0;j<NumberOfDihedralConstraints[CurrentSystem];j++)
    {
      switch(DihedralConstraintType)
      {
        case COS_PHI_SQUARED:
          Phi0=SQR(cos(DihedralConstraintParameter[CurrentSystem][j]));
          break;
        case COS_PHI:
          Phi0=cos(DihedralConstraintParameter[CurrentSystem][j]);
          break;
        case PHI:
        default:
          Phi0=DihedralConstraintParameter[CurrentSystem][j];
          break;
      }

      posA=DihedralConstraints[CurrentSystem][j][0]->Position;
      posB=DihedralConstraints[CurrentSystem][j][1]->Position;
      posC=DihedralConstraints[CurrentSystem][j][2]->Position;
      posD=DihedralConstraints[CurrentSystem][j][3]->Position;

      Vab.x=DihedralConstraintsDerivatives[CurrentSystem][j][1].x-DihedralConstraintsDerivatives[CurrentSystem][j][0].x;
      Vab.y=DihedralConstraintsDerivatives[CurrentSystem][j][1].y-DihedralConstraintsDerivatives[CurrentSystem][j][0].y;
      Vab.z=DihedralConstraintsDerivatives[CurrentSystem][j][1].z-DihedralConstraintsDerivatives[CurrentSystem][j][0].z;

      Vbc.x=DihedralConstraintsDerivatives[CurrentSystem][j][2].x-DihedralConstraintsDerivatives[CurrentSystem][j][1].x;
      Vbc.y=DihedralConstraintsDerivatives[CurrentSystem][j][2].y-DihedralConstraintsDerivatives[CurrentSystem][j][1].y;
      Vbc.z=DihedralConstraintsDerivatives[CurrentSystem][j][2].z-DihedralConstraintsDerivatives[CurrentSystem][j][1].z;

      Vcd.x=DihedralConstraintsDerivatives[CurrentSystem][j][3].x-DihedralConstraintsDerivatives[CurrentSystem][j][2].x;
      Vcd.y=DihedralConstraintsDerivatives[CurrentSystem][j][3].y-DihedralConstraintsDerivatives[CurrentSystem][j][2].y;
      Vcd.z=DihedralConstraintsDerivatives[CurrentSystem][j][3].z-DihedralConstraintsDerivatives[CurrentSystem][j][2].z;

      Phi=ReturnConstraintDihedralAngle(posA,posB,posC,posD);

      Rab.x=posB.x-posA.x;
      Rab.y=posB.y-posA.y;
      Rab.z=posB.z-posA.z;
      Rab=ApplyBoundaryCondition(Rab);

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      Rbc=ApplyBoundaryCondition(Rbc);

      Rcd.x=posD.x-posC.x;
      Rcd.y=posD.y-posC.y;
      Rcd.z=posD.z-posC.z;
      Rcd=ApplyBoundaryCondition(Rcd);

      gamma=(Phi-Phi0)/ReturnTorsionConstrainDerivative(Rab,Rbc,Rcd,Vab,Vbc,Vcd);

      DihedralConstraints[CurrentSystem][j][0]->Position.x-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][0].x;
      DihedralConstraints[CurrentSystem][j][0]->Position.y-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][0].y;
      DihedralConstraints[CurrentSystem][j][0]->Position.z-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][0].z;

      DihedralConstraints[CurrentSystem][j][1]->Position.x-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][1].x;
      DihedralConstraints[CurrentSystem][j][1]->Position.y-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][1].y;
      DihedralConstraints[CurrentSystem][j][1]->Position.z-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][1].z;

      DihedralConstraints[CurrentSystem][j][2]->Position.x-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][2].x;
      DihedralConstraints[CurrentSystem][j][2]->Position.y-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][2].y;
      DihedralConstraints[CurrentSystem][j][2]->Position.z-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][2].z;

      DihedralConstraints[CurrentSystem][j][3]->Position.x-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][3].x;
      DihedralConstraints[CurrentSystem][j][3]->Position.y-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][3].y;
      DihedralConstraints[CurrentSystem][j][3]->Position.z-=gamma*DihedralConstraintsDerivatives[CurrentSystem][j][3].z;

      error=fabs(Phi-Phi0);
      max_error=MAX2(error,max_error);
    }

    for(j=0;j<NumberOfImproperDihedralConstraints[CurrentSystem];j++)
    {
      switch(ImproperDihedralConstraintType)
      {
        case COS_PHI_SQUARED:
          Phi0=SQR(cos(ImproperDihedralConstraintParameter[CurrentSystem][j]));
          break;
        case COS_PHI:
          Phi0=cos(ImproperDihedralConstraintParameter[CurrentSystem][j]);
          break;
        case PHI:
        default:
          Phi0=ImproperDihedralConstraintParameter[CurrentSystem][j];
          break;
      }

      posA=ImproperDihedralConstraints[CurrentSystem][j][0]->Position;
      posB=ImproperDihedralConstraints[CurrentSystem][j][1]->Position;
      posC=ImproperDihedralConstraints[CurrentSystem][j][2]->Position;
      posD=ImproperDihedralConstraints[CurrentSystem][j][3]->Position;

      Vab.x=ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].x-ImproperDihedralConstraintsDerivatives[CurrentSystem][j][0].x;
      Vab.y=ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].y-ImproperDihedralConstraintsDerivatives[CurrentSystem][j][0].y;
      Vab.z=ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].z-ImproperDihedralConstraintsDerivatives[CurrentSystem][j][0].z;

      Vbc.x=ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].x-ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].x;
      Vbc.y=ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].y-ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].y;
      Vbc.z=ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].z-ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].z;

      Vcd.x=ImproperDihedralConstraintsDerivatives[CurrentSystem][j][3].x-ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].x;
      Vcd.y=ImproperDihedralConstraintsDerivatives[CurrentSystem][j][3].y-ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].y;
      Vcd.z=ImproperDihedralConstraintsDerivatives[CurrentSystem][j][3].z-ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].z;

      Phi=ReturnConstraintDihedralAngle(posA,posB,posC,posD);

      Rab.x=posB.x-posA.x;
      Rab.y=posB.y-posA.y;
      Rab.z=posB.z-posA.z;
      Rab=ApplyBoundaryCondition(Rab);

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      Rbc=ApplyBoundaryCondition(Rbc);

      Rcd.x=posD.x-posC.x;
      Rcd.y=posD.y-posC.y;
      Rcd.z=posD.z-posC.z;
      Rcd=ApplyBoundaryCondition(Rcd);

      gamma=(Phi-Phi0)/ReturnTorsionConstrainDerivative(Rab,Rbc,Rcd,Vab,Vbc,Vcd);

      ImproperDihedralConstraints[CurrentSystem][j][0]->Position.x-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][0].x;
      ImproperDihedralConstraints[CurrentSystem][j][0]->Position.y-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][0].y;
      ImproperDihedralConstraints[CurrentSystem][j][0]->Position.z-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][0].z;

      ImproperDihedralConstraints[CurrentSystem][j][1]->Position.x-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].x;
      ImproperDihedralConstraints[CurrentSystem][j][1]->Position.y-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].y;
      ImproperDihedralConstraints[CurrentSystem][j][1]->Position.z-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][1].z;

      ImproperDihedralConstraints[CurrentSystem][j][2]->Position.x-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].x;
      ImproperDihedralConstraints[CurrentSystem][j][2]->Position.y-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].y;
      ImproperDihedralConstraints[CurrentSystem][j][2]->Position.z-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][2].z;

      ImproperDihedralConstraints[CurrentSystem][j][3]->Position.x-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][3].x;
      ImproperDihedralConstraints[CurrentSystem][j][3]->Position.y-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][3].y;
      ImproperDihedralConstraints[CurrentSystem][j][3]->Position.z-=gamma*ImproperDihedralConstraintsDerivatives[CurrentSystem][j][3].z;

      error=fabs(Phi-Phi0);
      max_error=MAX2(error,max_error);

/*
      switch(ImproperDihedralConstraintType)
      {
        case COS_PHI_SQUARED:
          fprintf(stderr, "Gamma %g dihedral angle: %g %g\n",gamma,acos(sqrt(Phi))*RAD2DEG,acos(sqrt(Phi0))*RAD2DEG);
          break;
        case COS_PHI:
          fprintf(stderr, "Gamma %g dihedral angle: %g %g\n",gamma,acos(Phi)*RAD2DEG,acos(Phi0)*RAD2DEG);
          break;
        case PHI:
        default:
          fprintf(stderr, "Gamma %g dihedral angle: %g %g\n",gamma,(Phi)*RAD2DEG,(Phi0)*RAD2DEG);
          break;
      }
*/
    }
  }while(max_error>1e-8);

  for(j=0;j<NumberOfDistanceConstraints[CurrentSystem];j++)
  {
    posA=DistanceConstraints[CurrentSystem][j][0]->Position;
    posB=DistanceConstraints[CurrentSystem][j][1]->Position;
    fprintf(stderr, "Distance: %18.10f\n",ReturnBondDistance(posA,posB));
  }
  for(j=0;j<NumberOfAngleConstraints[CurrentSystem];j++)
  {
    posA=AngleConstraints[CurrentSystem][j][0]->Position;
    posB=AngleConstraints[CurrentSystem][j][1]->Position;
    posC=AngleConstraints[CurrentSystem][j][2]->Position;
    fprintf(stderr, "Bend Angle: %18.10f\n",ReturnBendAngle(posA,posB,posC)*RAD2DEG);
  }
  for(j=0;j<NumberOfDihedralConstraints[CurrentSystem];j++)
  {
    posA=DihedralConstraints[CurrentSystem][j][0]->Position;
    posB=DihedralConstraints[CurrentSystem][j][1]->Position;
    posC=DihedralConstraints[CurrentSystem][j][2]->Position;
    posD=DihedralConstraints[CurrentSystem][j][3]->Position;
    fprintf(stderr, "Dihedral Angle: %18.10f\n",ReturnDihedralAngle(posA,posB,posC,posD)*RAD2DEG);
  }
  for(j=0;j<NumberOfImproperDihedralConstraints[CurrentSystem];j++)
  {
    posA=ImproperDihedralConstraints[CurrentSystem][j][0]->Position;
    posB=ImproperDihedralConstraints[CurrentSystem][j][1]->Position;
    posC=ImproperDihedralConstraints[CurrentSystem][j][2]->Position;
    posD=ImproperDihedralConstraints[CurrentSystem][j][3]->Position;
    fprintf(stderr, "Improper dihedral Angle: %18.10f\n",ReturnDihedralAngle(posA,posB,posC,posD)*RAD2DEG);
  }
  for(j=0;j<NumberOfInversionBendConstraints[CurrentSystem];j++)
  {
    posA=InversionBendConstraints[CurrentSystem][j][0]->Position;
    posB=InversionBendConstraints[CurrentSystem][j][1]->Position;
    posC=InversionBendConstraints[CurrentSystem][j][2]->Position;
    posD=InversionBendConstraints[CurrentSystem][j][3]->Position;
    Chi=ReturnInversionBendAngle(posA,posB,posC,posD);
    fprintf(stderr, "Inversionbend Angle: %18.10f\n",Chi*RAD2DEG);
  }
}

void RattleStageOne(void)
{
  int j,m;
  int A,B,C,D;
  int Type;
  int NumberOfRattleSteps;
  REAL gamma;
  REAL max_error,error;
  REAL sigma,ConstraintValue,ConstraintTarget;
  REAL InverseMassA,InverseMassB,InverseMassC,InverseMassD;
  VECTOR posA,posB,posC,posD;
  VECTOR ReferenceGradientA,ReferenceGradientB,ReferenceGradientC,ReferenceGradientD;
  VECTOR CurrentGradientA,CurrentGradientB,CurrentGradientC,CurrentGradientD;
  REAL dot_productA,dot_productB,dot_productC,dot_productD;
  VECTOR dr;

  NumberOfRattleSteps=0;
  do
  {
    max_error=0.0;
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(j=0;j<Components[Type].NumberOfBonds;j++)
      {
        if(Components[Type].BondType[j]==FIXED_BOND)
        {
          A=Components[Type].Bonds[j].A;
          B=Components[Type].Bonds[j].B;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;

          switch(DistanceConstraintType)
          {
            case DISTANCE_R:
              ConstraintTarget=Components[Type].BondArguments[j][0];
              break;
            default:
            case DISTANCE_R_SQUARED:
              ConstraintTarget=SQR(Components[Type].BondArguments[j][0]);
              break;
          }

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          ConstraintValue=ReturnWilsonVectorsDistanceRATTLE(posA,posB,&CurrentGradientA,&CurrentGradientB);

          posA=Adsorbates[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          ReturnWilsonVectorsDistanceRATTLE(posA,posB,&ReferenceGradientA,&ReferenceGradientB);
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;

          sigma=ConstraintValue-ConstraintTarget;

          error=0.5*fabs(sigma)/ConstraintTarget;
          max_error=MAX2(error,max_error);

          dot_productA=CurrentGradientA.x*ReferenceGradientA.x+CurrentGradientA.y*ReferenceGradientA.y+CurrentGradientA.z*ReferenceGradientA.z;
          dot_productB=CurrentGradientB.x*ReferenceGradientB.x+CurrentGradientB.y*ReferenceGradientB.y+CurrentGradientB.z*ReferenceGradientB.z;

          gamma=(2.0/SQR(DeltaT))*sigma/(InverseMassA*dot_productA+InverseMassB*dot_productB);

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Position.x-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.y-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.z-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.z;

          // add contribution to the stress tensor
          StrainDerivativeTensor[CurrentSystem].ax-=gamma*dr.x*dr.x;
          StrainDerivativeTensor[CurrentSystem].bx-=gamma*dr.y*dr.x;
          StrainDerivativeTensor[CurrentSystem].cx-=gamma*dr.z*dr.x;

          StrainDerivativeTensor[CurrentSystem].ay-=gamma*dr.x*dr.y;
          StrainDerivativeTensor[CurrentSystem].by-=gamma*dr.y*dr.y;
          StrainDerivativeTensor[CurrentSystem].cy-=gamma*dr.z*dr.y;

          StrainDerivativeTensor[CurrentSystem].az-=gamma*dr.x*dr.z;
          StrainDerivativeTensor[CurrentSystem].bz-=gamma*dr.y*dr.z;
          StrainDerivativeTensor[CurrentSystem].cz-=gamma*dr.z*dr.z;
        }
      }
      for(j=0;j<Components[Type].NumberOfBends;j++)
      {
        if(Components[Type].BendType[j]==FIXED_BEND)
        {
          A=Components[Type].Bends[j].A;
          B=Components[Type].Bends[j].B;
          C=Components[Type].Bends[j].C;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;

          switch(BendConstraintType)
          {
            case COS_THETA_SQUARED:
              ConstraintTarget=SQR(cos(Components[Type].BendArguments[j][0]));
              break;
            case COS_THETA:
              ConstraintTarget=cos(Components[Type].BendArguments[j][0]);
              break;
            default:
            case THETA:
              ConstraintTarget=Components[Type].BendArguments[j][0];
              break;
          }

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          ConstraintValue=ReturnWilsonVectorsBendRATTLE(posA,posB,posC,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC);
          posA=Adsorbates[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          ReturnWilsonVectorsBendRATTLE(posA,posB,posC,&ReferenceGradientA,&ReferenceGradientB,&ReferenceGradientC);

          sigma=ConstraintValue-ConstraintTarget;

          error=0.5*fabs(sigma)/ConstraintTarget;
          max_error=MAX2(error,max_error);

          dot_productA=CurrentGradientA.x*ReferenceGradientA.x+CurrentGradientA.y*ReferenceGradientA.y+CurrentGradientA.z*ReferenceGradientA.z;
          dot_productB=CurrentGradientB.x*ReferenceGradientB.x+CurrentGradientB.y*ReferenceGradientB.y+CurrentGradientB.z*ReferenceGradientB.z;
          dot_productC=CurrentGradientC.x*ReferenceGradientC.x+CurrentGradientC.y*ReferenceGradientC.y+CurrentGradientC.z*ReferenceGradientC.z;

          gamma=(2.0/SQR(DeltaT))*sigma/(InverseMassA*dot_productA+InverseMassB*dot_productB+InverseMassC*dot_productC);

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Position.x-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.y-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.z-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Position.x-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.y-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.z-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.z;

          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.x-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.y-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.z-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.z;
        }
      }
      for(j=0;j<Components[Type].NumberOfTorsions;j++)
      {
        if(Components[Type].TorsionType[j]==FIXED_DIHEDRAL)
        {
          A=Components[Type].Torsions[j].A;
          B=Components[Type].Torsions[j].B;
          C=Components[Type].Torsions[j].C;
          D=Components[Type].Torsions[j].D;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[D].Type].Mass;

          switch(DihedralConstraintType)
          {
            case COS_PHI_SQUARED:
              ConstraintTarget=SQR(cos(Components[Type].TorsionArguments[j][0]));
              break;
            case COS_PHI:
              ConstraintTarget=cos(Components[Type].TorsionArguments[j][0]);
              break;
            default:
            case PHI:
              ConstraintTarget=Components[Type].TorsionArguments[j][0];
              break;
          }

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;
          ConstraintValue=ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC,&CurrentGradientD);
          posA=Adsorbates[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].RattleReferencePosition;
          ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&ReferenceGradientA,&ReferenceGradientB,&ReferenceGradientC,&ReferenceGradientD);

          sigma=ConstraintValue-ConstraintTarget;

          error=0.5*fabs(sigma)/ConstraintTarget;
          max_error=MAX2(error,max_error);

          dot_productA=CurrentGradientA.x*ReferenceGradientA.x+CurrentGradientA.y*ReferenceGradientA.y+CurrentGradientA.z*ReferenceGradientA.z;
          dot_productB=CurrentGradientB.x*ReferenceGradientB.x+CurrentGradientB.y*ReferenceGradientB.y+CurrentGradientB.z*ReferenceGradientB.z;
          dot_productC=CurrentGradientC.x*ReferenceGradientC.x+CurrentGradientC.y*ReferenceGradientC.y+CurrentGradientC.z*ReferenceGradientC.z;
          dot_productD=CurrentGradientD.x*ReferenceGradientD.x+CurrentGradientD.y*ReferenceGradientD.y+CurrentGradientD.z*ReferenceGradientD.z;

          gamma=(2.0/SQR(DeltaT))*sigma/(InverseMassA*dot_productA+InverseMassB*dot_productB+InverseMassC*dot_productC+InverseMassD*dot_productD);

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Position.x-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.y-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.z-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Position.x-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.y-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.z-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Position.x-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.y-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.z-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.z;

          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.x-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.y-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.z-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.x-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.y-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.z-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.z;
        }
      }
      for(j=0;j<Components[Type].NumberOfImproperTorsions;j++)
      {
        if(Components[Type].ImproperTorsionType[j]==FIXED_IMPROPER_DIHEDRAL)
        {
          A=Components[Type].ImproperTorsions[j].A;
          B=Components[Type].ImproperTorsions[j].B;
          C=Components[Type].ImproperTorsions[j].C;
          D=Components[Type].ImproperTorsions[j].D;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[D].Type].Mass;

          switch(DihedralConstraintType)
          {
            case COS_PHI_SQUARED:
              ConstraintTarget=SQR(cos(Components[Type].TorsionArguments[j][0]));
              break;
            case COS_PHI:
              ConstraintTarget=cos(Components[Type].TorsionArguments[j][0]);
              break;
            default:
            case PHI:
              ConstraintTarget=Components[Type].TorsionArguments[j][0];
              break;
          }

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;
          ConstraintValue=ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC,&CurrentGradientD);
          posA=Adsorbates[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].RattleReferencePosition;
          ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&ReferenceGradientA,&ReferenceGradientB,&ReferenceGradientC,&ReferenceGradientD);

          sigma=ConstraintValue-ConstraintTarget;

          error=0.5*fabs(sigma)/ConstraintTarget;
          max_error=MAX2(error,max_error);

          dot_productA=CurrentGradientA.x*ReferenceGradientA.x+CurrentGradientA.y*ReferenceGradientA.y+CurrentGradientA.z*ReferenceGradientA.z;
          dot_productB=CurrentGradientB.x*ReferenceGradientB.x+CurrentGradientB.y*ReferenceGradientB.y+CurrentGradientB.z*ReferenceGradientB.z;
          dot_productC=CurrentGradientC.x*ReferenceGradientC.x+CurrentGradientC.y*ReferenceGradientC.y+CurrentGradientC.z*ReferenceGradientC.z;
          dot_productD=CurrentGradientD.x*ReferenceGradientD.x+CurrentGradientD.y*ReferenceGradientD.y+CurrentGradientD.z*ReferenceGradientD.z;

          gamma=(2.0/SQR(DeltaT))*sigma/(InverseMassA*dot_productA+InverseMassB*dot_productB+InverseMassC*dot_productC+InverseMassD*dot_productD);

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Position.x-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.y-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.z-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Position.x-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.y-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.z-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Position.x-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.y-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.z-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.z;

          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.x-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.y-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.z-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.x-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.y-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.z-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.z;
        }
      }

      for(j=0;j<Components[Type].NumberOfInversionBends;j++)
      {
        if(Components[Type].InversionBendType[j]==FIXED_INVERSION_BEND)
        {
          A=Components[Type].InversionBends[j].A;
          B=Components[Type].InversionBends[j].B;
          C=Components[Type].InversionBends[j].C;
          D=Components[Type].InversionBends[j].D;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[D].Type].Mass;

          switch(InversionBendConstraintType)
          {
            case SIN_CHI_SQUARED:
              ConstraintTarget=SQR(cos(Components[Type].InversionBendArguments[j][0]));
              break;
            case SIN_CHI:
              ConstraintTarget=cos(Components[Type].InversionBendArguments[j][0]);
              break;
            default:
            case CHI:
              ConstraintTarget=Components[Type].InversionBendArguments[j][0];
              break;
          }

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;
          ConstraintValue=ReturnWilsonVectorsInversionBendRATTLE(posA,posB,posC,posD,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC,&CurrentGradientD);
          posA=Adsorbates[CurrentSystem][m].Atoms[A].RattleReferencePosition;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].RattleReferencePosition;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].RattleReferencePosition;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].RattleReferencePosition;
          ReturnWilsonVectorsInversionBendRATTLE(posA,posB,posC,posD,&ReferenceGradientA,&ReferenceGradientB,&ReferenceGradientC,&ReferenceGradientD);

          sigma=ConstraintValue-ConstraintTarget;

          error=0.5*fabs(sigma);
          max_error=MAX2(error,max_error);

          dot_productA=CurrentGradientA.x*ReferenceGradientA.x+CurrentGradientA.y*ReferenceGradientA.y+CurrentGradientA.z*ReferenceGradientA.z;
          dot_productB=CurrentGradientB.x*ReferenceGradientB.x+CurrentGradientB.y*ReferenceGradientB.y+CurrentGradientB.z*ReferenceGradientB.z;
          dot_productC=CurrentGradientC.x*ReferenceGradientC.x+CurrentGradientC.y*ReferenceGradientC.y+CurrentGradientC.z*ReferenceGradientC.z;
          dot_productD=CurrentGradientD.x*ReferenceGradientD.x+CurrentGradientD.y*ReferenceGradientD.y+CurrentGradientD.z*ReferenceGradientD.z;

          gamma=(2.0/SQR(DeltaT))*sigma/(InverseMassA*dot_productA+InverseMassB*dot_productB+InverseMassC*dot_productC+InverseMassD*dot_productD);

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Position.x-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.y-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Position.z-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Position.x-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.y-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Position.z-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Position.x-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.y-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Position.z-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.z;

          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.x-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.y-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.z-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.x-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.y-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.z-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.z;
        }
      }
    }

    for(m=0;m<NumberOfDistanceConstraints[CurrentSystem];m++)
    {
      InverseMassA=1.0/PseudoAtoms[DistanceConstraints[CurrentSystem][m][0]->Type].Mass;
      InverseMassB=1.0/PseudoAtoms[DistanceConstraints[CurrentSystem][m][1]->Type].Mass;

      switch(DistanceConstraintType)
      {
        case DISTANCE_R:
          ConstraintTarget=DistanceConstraintParameter[CurrentSystem][m];
          break;
        default:
        case DISTANCE_R_SQUARED:
          ConstraintTarget=SQR(DistanceConstraintParameter[CurrentSystem][m]);
          break;
      }

      posA=DistanceConstraints[CurrentSystem][m][0]->Position;
      posB=DistanceConstraints[CurrentSystem][m][1]->Position;
      ConstraintValue=ReturnWilsonVectorsDistanceRATTLE(posA,posB,&CurrentGradientA,&CurrentGradientB);

      posA=DistanceConstraints[CurrentSystem][m][0]->RattleReferencePosition;
      posB=DistanceConstraints[CurrentSystem][m][1]->RattleReferencePosition;
      ReturnWilsonVectorsDistanceRATTLE(posA,posB,&ReferenceGradientA,&ReferenceGradientB);

      sigma=ConstraintValue-ConstraintTarget;

      error=0.5*fabs(sigma)/ConstraintTarget;
      max_error=MAX2(error,max_error);

      dot_productA=CurrentGradientA.x*ReferenceGradientA.x+CurrentGradientA.y*ReferenceGradientA.y+CurrentGradientA.z*ReferenceGradientA.z;
      dot_productB=CurrentGradientB.x*ReferenceGradientB.x+CurrentGradientB.y*ReferenceGradientB.y+CurrentGradientB.z*ReferenceGradientB.z;

      gamma=(2.0/SQR(DeltaT))*sigma/(InverseMassA*dot_productA+InverseMassB*dot_productB);

      DistanceConstraints[CurrentSystem][m][0]->Position.x-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.x;
      DistanceConstraints[CurrentSystem][m][0]->Position.y-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.y;
      DistanceConstraints[CurrentSystem][m][0]->Position.z-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.z;

      DistanceConstraints[CurrentSystem][m][1]->Position.x-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.x;
      DistanceConstraints[CurrentSystem][m][1]->Position.y-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.y;
      DistanceConstraints[CurrentSystem][m][1]->Position.z-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.z;

      DistanceConstraints[CurrentSystem][m][0]->Velocity.x-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.x;
      DistanceConstraints[CurrentSystem][m][0]->Velocity.y-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.y;
      DistanceConstraints[CurrentSystem][m][0]->Velocity.z-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.z;

      DistanceConstraints[CurrentSystem][m][1]->Velocity.x-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.x;
      DistanceConstraints[CurrentSystem][m][1]->Velocity.y-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.y;
      DistanceConstraints[CurrentSystem][m][1]->Velocity.z-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.z;
    }

    for(m=0;m<NumberOfAngleConstraints[CurrentSystem];m++)
    {
      InverseMassA=1.0/PseudoAtoms[AngleConstraints[CurrentSystem][m][0]->Type].Mass;
      InverseMassB=1.0/PseudoAtoms[AngleConstraints[CurrentSystem][m][1]->Type].Mass;
      InverseMassC=1.0/PseudoAtoms[AngleConstraints[CurrentSystem][m][2]->Type].Mass;

      switch(BendConstraintType)
      {
        case COS_THETA_SQUARED:
          ConstraintTarget=SQR(cos(AngleConstraintParameter[CurrentSystem][m]));
          break;
        case COS_THETA:
          ConstraintTarget=cos(AngleConstraintParameter[CurrentSystem][m]);
          break;
        default:
        case THETA:
          ConstraintTarget=AngleConstraintParameter[CurrentSystem][m];
          break;
      }

      posA=AngleConstraints[CurrentSystem][m][0]->Position;
      posB=AngleConstraints[CurrentSystem][m][1]->Position;
      posC=AngleConstraints[CurrentSystem][m][2]->Position;
      ConstraintValue=ReturnWilsonVectorsBendRATTLE(posA,posB,posC,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC);
      posA=AngleConstraints[CurrentSystem][m][0]->RattleReferencePosition;
      posB=AngleConstraints[CurrentSystem][m][1]->RattleReferencePosition;
      posC=AngleConstraints[CurrentSystem][m][2]->RattleReferencePosition;
      ReturnWilsonVectorsBendRATTLE(posA,posB,posC,&ReferenceGradientA,&ReferenceGradientB,&ReferenceGradientC);

      sigma=ConstraintValue-ConstraintTarget;

      error=0.5*fabs(sigma)/ConstraintTarget;
      max_error=MAX2(error,max_error);

      dot_productA=CurrentGradientA.x*ReferenceGradientA.x+CurrentGradientA.y*ReferenceGradientA.y+CurrentGradientA.z*ReferenceGradientA.z;
      dot_productB=CurrentGradientB.x*ReferenceGradientB.x+CurrentGradientB.y*ReferenceGradientB.y+CurrentGradientB.z*ReferenceGradientB.z;
      dot_productC=CurrentGradientC.x*ReferenceGradientC.x+CurrentGradientC.y*ReferenceGradientC.y+CurrentGradientC.z*ReferenceGradientC.z;

      gamma=(2.0/SQR(DeltaT))*sigma/(InverseMassA*dot_productA+InverseMassB*dot_productB+InverseMassC*dot_productC);

      AngleConstraints[CurrentSystem][m][0]->Position.x-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.x;
      AngleConstraints[CurrentSystem][m][0]->Position.y-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.y;
      AngleConstraints[CurrentSystem][m][0]->Position.z-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.z;

      AngleConstraints[CurrentSystem][m][1]->Position.x-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.x;
      AngleConstraints[CurrentSystem][m][1]->Position.y-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.y;
      AngleConstraints[CurrentSystem][m][1]->Position.z-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.z;

      AngleConstraints[CurrentSystem][m][2]->Position.x-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.x;
      AngleConstraints[CurrentSystem][m][2]->Position.y-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.y;
      AngleConstraints[CurrentSystem][m][2]->Position.z-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.z;

      AngleConstraints[CurrentSystem][m][0]->Velocity.x-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.x;
      AngleConstraints[CurrentSystem][m][0]->Velocity.y-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.y;
      AngleConstraints[CurrentSystem][m][0]->Velocity.z-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.z;

      AngleConstraints[CurrentSystem][m][1]->Velocity.x-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.x;
      AngleConstraints[CurrentSystem][m][1]->Velocity.y-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.y;
      AngleConstraints[CurrentSystem][m][1]->Velocity.z-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.z;

      AngleConstraints[CurrentSystem][m][2]->Velocity.x-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.x;
      AngleConstraints[CurrentSystem][m][2]->Velocity.y-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.y;
      AngleConstraints[CurrentSystem][m][2]->Velocity.z-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.z;
    }
    for(m=0;m<NumberOfDihedralConstraints[CurrentSystem];m++)
    {
      InverseMassA=1.0/PseudoAtoms[DihedralConstraints[CurrentSystem][m][0]->Type].Mass;
      InverseMassB=1.0/PseudoAtoms[DihedralConstraints[CurrentSystem][m][1]->Type].Mass;
      InverseMassC=1.0/PseudoAtoms[DihedralConstraints[CurrentSystem][m][2]->Type].Mass;
      InverseMassD=1.0/PseudoAtoms[DihedralConstraints[CurrentSystem][m][3]->Type].Mass;

      switch(DihedralConstraintType)
      {
        case COS_PHI_SQUARED:
          ConstraintTarget=SQR(cos(DihedralConstraintParameter[CurrentSystem][m]));
          break;
        case COS_PHI:
          ConstraintTarget=cos(DihedralConstraintParameter[CurrentSystem][m]);
          break;
        default:
        case PHI:
          ConstraintTarget=DihedralConstraintParameter[CurrentSystem][m];
          break;
      }

      posA=DihedralConstraints[CurrentSystem][m][0]->Position;
      posB=DihedralConstraints[CurrentSystem][m][1]->Position;
      posC=DihedralConstraints[CurrentSystem][m][2]->Position;
      posD=DihedralConstraints[CurrentSystem][m][3]->Position;
      ConstraintValue=ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC,&CurrentGradientD);
      posA=DihedralConstraints[CurrentSystem][m][0]->RattleReferencePosition;
      posB=DihedralConstraints[CurrentSystem][m][1]->RattleReferencePosition;
      posC=DihedralConstraints[CurrentSystem][m][2]->RattleReferencePosition;
      posD=DihedralConstraints[CurrentSystem][m][3]->RattleReferencePosition;
      ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&ReferenceGradientA,&ReferenceGradientB,&ReferenceGradientC,&ReferenceGradientD);

      sigma=ConstraintValue-ConstraintTarget;

      error=0.5*fabs(sigma)/ConstraintTarget;
      max_error=MAX2(error,max_error);

      dot_productA=CurrentGradientA.x*ReferenceGradientA.x+CurrentGradientA.y*ReferenceGradientA.y+CurrentGradientA.z*ReferenceGradientA.z;
      dot_productB=CurrentGradientB.x*ReferenceGradientB.x+CurrentGradientB.y*ReferenceGradientB.y+CurrentGradientB.z*ReferenceGradientB.z;
      dot_productC=CurrentGradientC.x*ReferenceGradientC.x+CurrentGradientC.y*ReferenceGradientC.y+CurrentGradientC.z*ReferenceGradientC.z;
      dot_productD=CurrentGradientD.x*ReferenceGradientD.x+CurrentGradientD.y*ReferenceGradientD.y+CurrentGradientD.z*ReferenceGradientD.z;

      gamma=(2.0/SQR(DeltaT))*sigma/(InverseMassA*dot_productA+InverseMassB*dot_productB+InverseMassC*dot_productC+InverseMassD*dot_productD);

      DihedralConstraints[CurrentSystem][m][0]->Position.x-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.x;
      DihedralConstraints[CurrentSystem][m][0]->Position.y-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.y;
      DihedralConstraints[CurrentSystem][m][0]->Position.z-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.z;

      DihedralConstraints[CurrentSystem][m][1]->Position.x-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.x;
      DihedralConstraints[CurrentSystem][m][1]->Position.y-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.y;
      DihedralConstraints[CurrentSystem][m][1]->Position.z-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.z;

      DihedralConstraints[CurrentSystem][m][2]->Position.x-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.x;
      DihedralConstraints[CurrentSystem][m][2]->Position.y-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.y;
      DihedralConstraints[CurrentSystem][m][2]->Position.z-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.z;

      DihedralConstraints[CurrentSystem][m][3]->Position.x-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.x;
      DihedralConstraints[CurrentSystem][m][3]->Position.y-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.y;
      DihedralConstraints[CurrentSystem][m][3]->Position.z-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.z;

      DihedralConstraints[CurrentSystem][m][0]->Velocity.x-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.x;
      DihedralConstraints[CurrentSystem][m][0]->Velocity.y-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.y;
      DihedralConstraints[CurrentSystem][m][0]->Velocity.z-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.z;

      DihedralConstraints[CurrentSystem][m][1]->Velocity.x-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.x;
      DihedralConstraints[CurrentSystem][m][1]->Velocity.y-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.y;
      DihedralConstraints[CurrentSystem][m][1]->Velocity.z-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.z;

      DihedralConstraints[CurrentSystem][m][2]->Velocity.x-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.x;
      DihedralConstraints[CurrentSystem][m][2]->Velocity.y-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.y;
      DihedralConstraints[CurrentSystem][m][2]->Velocity.z-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.z;

      DihedralConstraints[CurrentSystem][m][3]->Velocity.x-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.x;
      DihedralConstraints[CurrentSystem][m][3]->Velocity.y-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.y;
      DihedralConstraints[CurrentSystem][m][3]->Velocity.z-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.z;
    }

    for(m=0;m<NumberOfImproperDihedralConstraints[CurrentSystem];m++)
    {
      InverseMassA=1.0/PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][m][0]->Type].Mass;
      InverseMassB=1.0/PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][m][1]->Type].Mass;
      InverseMassC=1.0/PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][m][2]->Type].Mass;
      InverseMassD=1.0/PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][m][3]->Type].Mass;

      switch(ImproperDihedralConstraintType)
      {
        case COS_PHI_SQUARED:
          ConstraintTarget=SQR(cos(ImproperDihedralConstraintParameter[CurrentSystem][m]));
          break;
        case COS_PHI:
          ConstraintTarget=cos(ImproperDihedralConstraintParameter[CurrentSystem][m]);
          break;
        default:
        case PHI:
          ConstraintTarget=ImproperDihedralConstraintParameter[CurrentSystem][m];
          break;
      }

      posA=ImproperDihedralConstraints[CurrentSystem][m][0]->Position;
      posB=ImproperDihedralConstraints[CurrentSystem][m][1]->Position;
      posC=ImproperDihedralConstraints[CurrentSystem][m][2]->Position;
      posD=ImproperDihedralConstraints[CurrentSystem][m][3]->Position;
      ConstraintValue=ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC,&CurrentGradientD);
      posA=ImproperDihedralConstraints[CurrentSystem][m][0]->RattleReferencePosition;
      posB=ImproperDihedralConstraints[CurrentSystem][m][1]->RattleReferencePosition;
      posC=ImproperDihedralConstraints[CurrentSystem][m][2]->RattleReferencePosition;
      posD=ImproperDihedralConstraints[CurrentSystem][m][3]->RattleReferencePosition;
      ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&ReferenceGradientA,&ReferenceGradientB,&ReferenceGradientC,&ReferenceGradientD);

      sigma=ConstraintValue-ConstraintTarget;

      error=0.5*fabs(sigma)/ConstraintTarget;
      max_error=MAX2(error,max_error);

      dot_productA=CurrentGradientA.x*ReferenceGradientA.x+CurrentGradientA.y*ReferenceGradientA.y+CurrentGradientA.z*ReferenceGradientA.z;
      dot_productB=CurrentGradientB.x*ReferenceGradientB.x+CurrentGradientB.y*ReferenceGradientB.y+CurrentGradientB.z*ReferenceGradientB.z;
      dot_productC=CurrentGradientC.x*ReferenceGradientC.x+CurrentGradientC.y*ReferenceGradientC.y+CurrentGradientC.z*ReferenceGradientC.z;
      dot_productD=CurrentGradientD.x*ReferenceGradientD.x+CurrentGradientD.y*ReferenceGradientD.y+CurrentGradientD.z*ReferenceGradientD.z;

      gamma=(2.0/SQR(DeltaT))*sigma/(InverseMassA*dot_productA+InverseMassB*dot_productB+InverseMassC*dot_productC+InverseMassD*dot_productD);

      ImproperDihedralConstraints[CurrentSystem][m][0]->Position.x-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.x;
      ImproperDihedralConstraints[CurrentSystem][m][0]->Position.y-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.y;
      ImproperDihedralConstraints[CurrentSystem][m][0]->Position.z-=0.5*SQR(DeltaT)*InverseMassA*gamma*ReferenceGradientA.z;

      ImproperDihedralConstraints[CurrentSystem][m][1]->Position.x-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.x;
      ImproperDihedralConstraints[CurrentSystem][m][1]->Position.y-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.y;
      ImproperDihedralConstraints[CurrentSystem][m][1]->Position.z-=0.5*SQR(DeltaT)*InverseMassB*gamma*ReferenceGradientB.z;

      ImproperDihedralConstraints[CurrentSystem][m][2]->Position.x-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.x;
      ImproperDihedralConstraints[CurrentSystem][m][2]->Position.y-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.y;
      ImproperDihedralConstraints[CurrentSystem][m][2]->Position.z-=0.5*SQR(DeltaT)*InverseMassC*gamma*ReferenceGradientC.z;

      ImproperDihedralConstraints[CurrentSystem][m][3]->Position.x-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.x;
      ImproperDihedralConstraints[CurrentSystem][m][3]->Position.y-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.y;
      ImproperDihedralConstraints[CurrentSystem][m][3]->Position.z-=0.5*SQR(DeltaT)*InverseMassD*gamma*ReferenceGradientD.z;

      ImproperDihedralConstraints[CurrentSystem][m][0]->Velocity.x-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.x;
      ImproperDihedralConstraints[CurrentSystem][m][0]->Velocity.y-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.y;
      ImproperDihedralConstraints[CurrentSystem][m][0]->Velocity.z-=0.5*DeltaT*InverseMassA*gamma*ReferenceGradientA.z;

      ImproperDihedralConstraints[CurrentSystem][m][1]->Velocity.x-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.x;
      ImproperDihedralConstraints[CurrentSystem][m][1]->Velocity.y-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.y;
      ImproperDihedralConstraints[CurrentSystem][m][1]->Velocity.z-=0.5*DeltaT*InverseMassB*gamma*ReferenceGradientB.z;

      ImproperDihedralConstraints[CurrentSystem][m][2]->Velocity.x-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.x;
      ImproperDihedralConstraints[CurrentSystem][m][2]->Velocity.y-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.y;
      ImproperDihedralConstraints[CurrentSystem][m][2]->Velocity.z-=0.5*DeltaT*InverseMassC*gamma*ReferenceGradientC.z;

      ImproperDihedralConstraints[CurrentSystem][m][3]->Velocity.x-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.x;
      ImproperDihedralConstraints[CurrentSystem][m][3]->Velocity.y-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.y;
      ImproperDihedralConstraints[CurrentSystem][m][3]->Velocity.z-=0.5*DeltaT*InverseMassD*gamma*ReferenceGradientD.z;
    }


    NumberOfRattleSteps++;
  }while(max_error>1e-8);

  NumberOfRattleCyclesStage1[CurrentSystem]+=NumberOfRattleSteps;
  if(NumberOfRattleSteps>MaximumNumberOfRattleCyclesStage1[CurrentSystem])
    MaximumNumberOfRattleCyclesStage1[CurrentSystem]=NumberOfRattleSteps;
}

void RattleStageTwo(void)
{
  int j,m;
  int A,B,C,D;
  int Type;
  int NumberOfRattleSteps;
  VECTOR velA,velB,velC,velD;
  REAL gamma;
  REAL max_error;
  REAL ConstraintValue;
  REAL InverseMassA,InverseMassB,InverseMassC,InverseMassD;
  VECTOR posA,posB,posC,posD;
  VECTOR CurrentGradientA,CurrentGradientB,CurrentGradientC,CurrentGradientD;
  REAL dot_productA,dot_productB,dot_productC,dot_productD;

  NumberOfRattleSteps=0;
  do
  {
    max_error=0.0;

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;

      for(j=0;j<Components[Type].NumberOfBonds;j++)
      {
        if(Components[Type].BondType[j]==FIXED_BOND)
        {
          A=Components[Type].Bonds[j].A;
          B=Components[Type].Bonds[j].B;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;

          velA=Adsorbates[CurrentSystem][m].Atoms[A].Velocity;
          velB=Adsorbates[CurrentSystem][m].Atoms[B].Velocity;

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

          ConstraintValue=ReturnWilsonVectorsDistanceRATTLE(posA,posB,&CurrentGradientA,&CurrentGradientB);

          dot_productA=CurrentGradientA.x*velA.x+CurrentGradientA.y*velA.y+CurrentGradientA.z*velA.z;
          dot_productB=CurrentGradientB.x*velB.x+CurrentGradientB.y*velB.y+CurrentGradientB.z*velB.z;

          gamma=(2.0/DeltaT)*(dot_productA+dot_productB)/
                 (InverseMassA*(SQR(CurrentGradientA.x)+SQR(CurrentGradientA.y)+SQR(CurrentGradientA.z))+
                  InverseMassB*(SQR(CurrentGradientB.x)+SQR(CurrentGradientB.y)+SQR(CurrentGradientB.z)));

          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.z;

          max_error=MAX2(fabs(gamma),max_error);
        }
      }
      for(j=0;j<Components[Type].NumberOfBends;j++)
      {
        if(Components[Type].BendType[j]==FIXED_BEND)
        {
          A=Components[Type].Bends[j].A;
          B=Components[Type].Bends[j].B;
          C=Components[Type].Bends[j].C;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;

          velA=Adsorbates[CurrentSystem][m].Atoms[A].Velocity;
          velB=Adsorbates[CurrentSystem][m].Atoms[B].Velocity;
          velC=Adsorbates[CurrentSystem][m].Atoms[C].Velocity;

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          ConstraintValue=ReturnWilsonVectorsBendRATTLE(posA,posB,posC,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC);

          dot_productA=CurrentGradientA.x*velA.x+CurrentGradientA.y*velA.y+CurrentGradientA.z*velA.z;
          dot_productB=CurrentGradientB.x*velB.x+CurrentGradientB.y*velB.y+CurrentGradientB.z*velB.z;
          dot_productC=CurrentGradientC.x*velC.x+CurrentGradientC.y*velC.y+CurrentGradientC.z*velC.z;

          gamma=(2.0/DeltaT)*(dot_productA+dot_productB+dot_productC)/
                 (InverseMassA*(SQR(CurrentGradientA.x)+SQR(CurrentGradientA.y)+SQR(CurrentGradientA.z))+
                  InverseMassB*(SQR(CurrentGradientB.x)+SQR(CurrentGradientB.y)+SQR(CurrentGradientB.z))+
                  InverseMassC*(SQR(CurrentGradientC.x)+SQR(CurrentGradientC.y)+SQR(CurrentGradientC.z)));

          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.x-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.y-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.z-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.z;

          max_error=MAX2(fabs(gamma),max_error);
        }
      }
      for(j=0;j<Components[Type].NumberOfTorsions;j++)
      {
        if(Components[Type].TorsionType[j]==FIXED_DIHEDRAL)
        {
          A=Components[Type].Torsions[j].A;
          B=Components[Type].Torsions[j].B;
          C=Components[Type].Torsions[j].C;
          D=Components[Type].Torsions[j].D;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[D].Type].Mass;

          velA=Adsorbates[CurrentSystem][m].Atoms[A].Velocity;
          velB=Adsorbates[CurrentSystem][m].Atoms[B].Velocity;
          velC=Adsorbates[CurrentSystem][m].Atoms[C].Velocity;
          velD=Adsorbates[CurrentSystem][m].Atoms[D].Velocity;

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;
          ConstraintValue=ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC,&CurrentGradientD);

          dot_productA=CurrentGradientA.x*velA.x+CurrentGradientA.y*velA.y+CurrentGradientA.z*velA.z;
          dot_productB=CurrentGradientB.x*velB.x+CurrentGradientB.y*velB.y+CurrentGradientB.z*velB.z;
          dot_productC=CurrentGradientC.x*velC.x+CurrentGradientC.y*velC.y+CurrentGradientC.z*velC.z;
          dot_productD=CurrentGradientD.x*velD.x+CurrentGradientD.y*velD.y+CurrentGradientD.z*velD.z;

          gamma=(2.0/DeltaT)*(dot_productA+dot_productB+dot_productC+dot_productD)/
                 (InverseMassA*(SQR(CurrentGradientA.x)+SQR(CurrentGradientA.y)+SQR(CurrentGradientA.z))+
                  InverseMassB*(SQR(CurrentGradientB.x)+SQR(CurrentGradientB.y)+SQR(CurrentGradientB.z))+
                  InverseMassC*(SQR(CurrentGradientC.x)+SQR(CurrentGradientC.y)+SQR(CurrentGradientC.z))+
                  InverseMassD*(SQR(CurrentGradientD.x)+SQR(CurrentGradientD.y)+SQR(CurrentGradientD.z)));

          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.x-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.y-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.z-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.x-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.y-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.z-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.z;

          max_error=MAX2(fabs(gamma),max_error);
        }
      }
      for(j=0;j<Components[Type].NumberOfImproperTorsions;j++)
      {
        if(Components[Type].ImproperTorsionType[j]==FIXED_IMPROPER_DIHEDRAL)
        {
          A=Components[Type].ImproperTorsions[j].A;
          B=Components[Type].ImproperTorsions[j].B;
          C=Components[Type].ImproperTorsions[j].C;
          D=Components[Type].ImproperTorsions[j].D;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[D].Type].Mass;

          velA=Adsorbates[CurrentSystem][m].Atoms[A].Velocity;
          velB=Adsorbates[CurrentSystem][m].Atoms[B].Velocity;
          velC=Adsorbates[CurrentSystem][m].Atoms[C].Velocity;
          velD=Adsorbates[CurrentSystem][m].Atoms[D].Velocity;

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;
          ConstraintValue=ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC,&CurrentGradientD);

          dot_productA=CurrentGradientA.x*velA.x+CurrentGradientA.y*velA.y+CurrentGradientA.z*velA.z;
          dot_productB=CurrentGradientB.x*velB.x+CurrentGradientB.y*velB.y+CurrentGradientB.z*velB.z;
          dot_productC=CurrentGradientC.x*velC.x+CurrentGradientC.y*velC.y+CurrentGradientC.z*velC.z;
          dot_productD=CurrentGradientD.x*velD.x+CurrentGradientD.y*velD.y+CurrentGradientD.z*velD.z;

          gamma=(2.0/DeltaT)*(dot_productA+dot_productB+dot_productC+dot_productD)/
                 (InverseMassA*(SQR(CurrentGradientA.x)+SQR(CurrentGradientA.y)+SQR(CurrentGradientA.z))+
                  InverseMassB*(SQR(CurrentGradientB.x)+SQR(CurrentGradientB.y)+SQR(CurrentGradientB.z))+
                  InverseMassC*(SQR(CurrentGradientC.x)+SQR(CurrentGradientC.y)+SQR(CurrentGradientC.z))+
                  InverseMassD*(SQR(CurrentGradientD.x)+SQR(CurrentGradientD.y)+SQR(CurrentGradientD.z)));

          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.x-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.y-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.z-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.x-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.y-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.z-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.z;

          max_error=MAX2(fabs(gamma),max_error);
        }
      }

      for(j=0;j<Components[Type].NumberOfInversionBends;j++)
      {
        if(Components[Type].InversionBendType[j]==FIXED_INVERSION_BEND)
        {
          A=Components[Type].InversionBends[j].A;
          B=Components[Type].InversionBends[j].B;
          C=Components[Type].InversionBends[j].C;
          D=Components[Type].InversionBends[j].D;

          InverseMassA=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          InverseMassB=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[B].Type].Mass;
          InverseMassC=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[C].Type].Mass;
          InverseMassD=1.0/PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[D].Type].Mass;

          velA=Adsorbates[CurrentSystem][m].Atoms[A].Velocity;
          velB=Adsorbates[CurrentSystem][m].Atoms[B].Velocity;
          velC=Adsorbates[CurrentSystem][m].Atoms[C].Velocity;
          velD=Adsorbates[CurrentSystem][m].Atoms[D].Velocity;

          posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
          posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
          posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;
          ConstraintValue=ReturnWilsonVectorsInversionBendRATTLE(posA,posB,posC,posD,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC,&CurrentGradientD);

          dot_productA=CurrentGradientA.x*velA.x+CurrentGradientA.y*velA.y+CurrentGradientA.z*velA.z;
          dot_productB=CurrentGradientB.x*velB.x+CurrentGradientB.y*velB.y+CurrentGradientB.z*velB.z;
          dot_productC=CurrentGradientC.x*velC.x+CurrentGradientC.y*velC.y+CurrentGradientC.z*velC.z;
          dot_productD=CurrentGradientD.x*velD.x+CurrentGradientD.y*velD.y+CurrentGradientD.z*velD.z;

          gamma=(2.0/DeltaT)*(dot_productA+dot_productB+dot_productC+dot_productD)/
                 (InverseMassA*(SQR(CurrentGradientA.x)+SQR(CurrentGradientA.y)+SQR(CurrentGradientA.z))+
                  InverseMassB*(SQR(CurrentGradientB.x)+SQR(CurrentGradientB.y)+SQR(CurrentGradientB.z))+
                  InverseMassC*(SQR(CurrentGradientC.x)+SQR(CurrentGradientC.y)+SQR(CurrentGradientC.z))+
                  InverseMassD*(SQR(CurrentGradientD.x)+SQR(CurrentGradientD.y)+SQR(CurrentGradientD.z)));

          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.x;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.y;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.z;

          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.x;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.y;
          Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.z;

          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.x-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.x;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.y-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.y;
          Adsorbates[CurrentSystem][m].Atoms[C].Velocity.z-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.z;

          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.x-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.x;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.y-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.y;
          Adsorbates[CurrentSystem][m].Atoms[D].Velocity.z-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.z;

          max_error=MAX2(fabs(gamma),max_error);
        }
      }
    }

    for(m=0;m<NumberOfDistanceConstraints[CurrentSystem];m++)
    {
      InverseMassA=1.0/PseudoAtoms[DistanceConstraints[CurrentSystem][m][0]->Type].Mass;
      InverseMassB=1.0/PseudoAtoms[DistanceConstraints[CurrentSystem][m][1]->Type].Mass;

      velA=DistanceConstraints[CurrentSystem][m][0]->Velocity;
      velB=DistanceConstraints[CurrentSystem][m][1]->Velocity;

      posA=DistanceConstraints[CurrentSystem][m][0]->Position;
      posB=DistanceConstraints[CurrentSystem][m][1]->Position;
      ConstraintValue=ReturnWilsonVectorsDistanceRATTLE(posA,posB,&CurrentGradientA,&CurrentGradientB);

      dot_productA=CurrentGradientA.x*velA.x+CurrentGradientA.y*velA.y+CurrentGradientA.z*velA.z;
      dot_productB=CurrentGradientB.x*velB.x+CurrentGradientB.y*velB.y+CurrentGradientB.z*velB.z;

      gamma=(2.0/DeltaT)*(dot_productA+dot_productB)/
            (InverseMassA*(SQR(CurrentGradientA.x)+SQR(CurrentGradientA.y)+SQR(CurrentGradientA.z))+
             InverseMassB*(SQR(CurrentGradientB.x)+SQR(CurrentGradientB.y)+SQR(CurrentGradientB.z)));

      DistanceConstraints[CurrentSystem][m][0]->Velocity.x-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.x;
      DistanceConstraints[CurrentSystem][m][0]->Velocity.y-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.y;
      DistanceConstraints[CurrentSystem][m][0]->Velocity.z-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.z;

      DistanceConstraints[CurrentSystem][m][1]->Velocity.x-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.x;
      DistanceConstraints[CurrentSystem][m][1]->Velocity.y-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.y;
      DistanceConstraints[CurrentSystem][m][1]->Velocity.z-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.z;

      max_error=MAX2(fabs(gamma),max_error);
    }

    for(m=0;m<NumberOfAngleConstraints[CurrentSystem];m++)
    {
      InverseMassA=1.0/PseudoAtoms[AngleConstraints[CurrentSystem][m][0]->Type].Mass;
      InverseMassB=1.0/PseudoAtoms[AngleConstraints[CurrentSystem][m][1]->Type].Mass;
      InverseMassC=1.0/PseudoAtoms[AngleConstraints[CurrentSystem][m][2]->Type].Mass;

      velA=AngleConstraints[CurrentSystem][m][0]->Velocity;
      velB=AngleConstraints[CurrentSystem][m][1]->Velocity;
      velC=AngleConstraints[CurrentSystem][m][2]->Velocity;

      posA=AngleConstraints[CurrentSystem][m][0]->Position;
      posB=AngleConstraints[CurrentSystem][m][1]->Position;
      posC=AngleConstraints[CurrentSystem][m][2]->Position;
      ConstraintValue=ReturnWilsonVectorsBendRATTLE(posA,posB,posC,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC);

      dot_productA=CurrentGradientA.x*velA.x+CurrentGradientA.y*velA.y+CurrentGradientA.z*velA.z;
      dot_productB=CurrentGradientB.x*velB.x+CurrentGradientB.y*velB.y+CurrentGradientB.z*velB.z;
      dot_productC=CurrentGradientC.x*velC.x+CurrentGradientC.y*velC.y+CurrentGradientC.z*velC.z;

      gamma=(2.0/DeltaT)*(dot_productA+dot_productB+dot_productC)/
            (InverseMassA*(SQR(CurrentGradientA.x)+SQR(CurrentGradientA.y)+SQR(CurrentGradientA.z))+
             InverseMassB*(SQR(CurrentGradientB.x)+SQR(CurrentGradientB.y)+SQR(CurrentGradientB.z))+
             InverseMassC*(SQR(CurrentGradientC.x)+SQR(CurrentGradientC.y)+SQR(CurrentGradientC.z)));

      AngleConstraints[CurrentSystem][m][0]->Velocity.x-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.x;
      AngleConstraints[CurrentSystem][m][0]->Velocity.y-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.y;
      AngleConstraints[CurrentSystem][m][0]->Velocity.z-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.z;

      AngleConstraints[CurrentSystem][m][1]->Velocity.x-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.x;
      AngleConstraints[CurrentSystem][m][1]->Velocity.y-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.y;
      AngleConstraints[CurrentSystem][m][1]->Velocity.z-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.z;

      AngleConstraints[CurrentSystem][m][2]->Velocity.x-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.x;
      AngleConstraints[CurrentSystem][m][2]->Velocity.y-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.y;
      AngleConstraints[CurrentSystem][m][2]->Velocity.z-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.z;

      max_error=MAX2(fabs(gamma),max_error);
    }

    for(m=0;m<NumberOfDihedralConstraints[CurrentSystem];m++)
    {
      InverseMassA=1.0/PseudoAtoms[DihedralConstraints[CurrentSystem][m][0]->Type].Mass;
      InverseMassB=1.0/PseudoAtoms[DihedralConstraints[CurrentSystem][m][1]->Type].Mass;
      InverseMassC=1.0/PseudoAtoms[DihedralConstraints[CurrentSystem][m][2]->Type].Mass;
      InverseMassD=1.0/PseudoAtoms[DihedralConstraints[CurrentSystem][m][3]->Type].Mass;

      velA=DihedralConstraints[CurrentSystem][m][0]->Velocity;
      velB=DihedralConstraints[CurrentSystem][m][1]->Velocity;
      velC=DihedralConstraints[CurrentSystem][m][2]->Velocity;
      velD=DihedralConstraints[CurrentSystem][m][3]->Velocity;

      posA=DihedralConstraints[CurrentSystem][m][0]->Position;
      posB=DihedralConstraints[CurrentSystem][m][1]->Position;
      posC=DihedralConstraints[CurrentSystem][m][2]->Position;
      posD=DihedralConstraints[CurrentSystem][m][3]->Position;
      ConstraintValue=ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC,&CurrentGradientD);

      dot_productA=CurrentGradientA.x*velA.x+CurrentGradientA.y*velA.y+CurrentGradientA.z*velA.z;
      dot_productB=CurrentGradientB.x*velB.x+CurrentGradientB.y*velB.y+CurrentGradientB.z*velB.z;
      dot_productC=CurrentGradientC.x*velC.x+CurrentGradientC.y*velC.y+CurrentGradientC.z*velC.z;
      dot_productD=CurrentGradientD.x*velD.x+CurrentGradientD.y*velD.y+CurrentGradientD.z*velD.z;

      gamma=(2.0/DeltaT)*(dot_productA+dot_productB+dot_productC+dot_productD)/
            (InverseMassA*(SQR(CurrentGradientA.x)+SQR(CurrentGradientA.y)+SQR(CurrentGradientA.z))+
             InverseMassB*(SQR(CurrentGradientB.x)+SQR(CurrentGradientB.y)+SQR(CurrentGradientB.z))+
             InverseMassC*(SQR(CurrentGradientC.x)+SQR(CurrentGradientC.y)+SQR(CurrentGradientC.z))+
             InverseMassD*(SQR(CurrentGradientD.x)+SQR(CurrentGradientD.y)+SQR(CurrentGradientD.z)));

      DihedralConstraints[CurrentSystem][m][0]->Velocity.x-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.x;
      DihedralConstraints[CurrentSystem][m][0]->Velocity.y-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.y;
      DihedralConstraints[CurrentSystem][m][0]->Velocity.z-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.z;

      DihedralConstraints[CurrentSystem][m][1]->Velocity.x-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.x;
      DihedralConstraints[CurrentSystem][m][1]->Velocity.y-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.y;
      DihedralConstraints[CurrentSystem][m][1]->Velocity.z-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.z;

      DihedralConstraints[CurrentSystem][m][2]->Velocity.x-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.x;
      DihedralConstraints[CurrentSystem][m][2]->Velocity.y-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.y;
      DihedralConstraints[CurrentSystem][m][2]->Velocity.z-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.z;

      DihedralConstraints[CurrentSystem][m][3]->Velocity.x-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.x;
      DihedralConstraints[CurrentSystem][m][3]->Velocity.y-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.y;
      DihedralConstraints[CurrentSystem][m][3]->Velocity.z-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.z;

      max_error=MAX2(fabs(gamma),max_error);
    }

    for(m=0;m<NumberOfImproperDihedralConstraints[CurrentSystem];m++)
    {
      InverseMassA=1.0/PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][m][0]->Type].Mass;
      InverseMassB=1.0/PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][m][1]->Type].Mass;
      InverseMassC=1.0/PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][m][2]->Type].Mass;
      InverseMassD=1.0/PseudoAtoms[ImproperDihedralConstraints[CurrentSystem][m][3]->Type].Mass;

      velA=ImproperDihedralConstraints[CurrentSystem][m][0]->Velocity;
      velB=ImproperDihedralConstraints[CurrentSystem][m][1]->Velocity;
      velC=ImproperDihedralConstraints[CurrentSystem][m][2]->Velocity;
      velD=ImproperDihedralConstraints[CurrentSystem][m][3]->Velocity;

      posA=ImproperDihedralConstraints[CurrentSystem][m][0]->Position;
      posB=ImproperDihedralConstraints[CurrentSystem][m][1]->Position;
      posC=ImproperDihedralConstraints[CurrentSystem][m][2]->Position;
      posD=ImproperDihedralConstraints[CurrentSystem][m][3]->Position;
      ConstraintValue=ReturnWilsonVectorsTorsionRATTLE(posA,posB,posC,posD,&CurrentGradientA,&CurrentGradientB,&CurrentGradientC,&CurrentGradientD);

      dot_productA=CurrentGradientA.x*velA.x+CurrentGradientA.y*velA.y+CurrentGradientA.z*velA.z;
      dot_productB=CurrentGradientB.x*velB.x+CurrentGradientB.y*velB.y+CurrentGradientB.z*velB.z;
      dot_productC=CurrentGradientC.x*velC.x+CurrentGradientC.y*velC.y+CurrentGradientC.z*velC.z;
      dot_productD=CurrentGradientD.x*velD.x+CurrentGradientD.y*velD.y+CurrentGradientD.z*velD.z;

      gamma=(2.0/DeltaT)*(dot_productA+dot_productB+dot_productC+dot_productD)/
            (InverseMassA*(SQR(CurrentGradientA.x)+SQR(CurrentGradientA.y)+SQR(CurrentGradientA.z))+
             InverseMassB*(SQR(CurrentGradientB.x)+SQR(CurrentGradientB.y)+SQR(CurrentGradientB.z))+
             InverseMassC*(SQR(CurrentGradientC.x)+SQR(CurrentGradientC.y)+SQR(CurrentGradientC.z))+
             InverseMassD*(SQR(CurrentGradientD.x)+SQR(CurrentGradientD.y)+SQR(CurrentGradientD.z)));

      ImproperDihedralConstraints[CurrentSystem][m][0]->Velocity.x-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.x;
      ImproperDihedralConstraints[CurrentSystem][m][0]->Velocity.y-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.y;
      ImproperDihedralConstraints[CurrentSystem][m][0]->Velocity.z-=0.5*DeltaT*InverseMassA*gamma*CurrentGradientA.z;

      ImproperDihedralConstraints[CurrentSystem][m][1]->Velocity.x-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.x;
      ImproperDihedralConstraints[CurrentSystem][m][1]->Velocity.y-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.y;
      ImproperDihedralConstraints[CurrentSystem][m][1]->Velocity.z-=0.5*DeltaT*InverseMassB*gamma*CurrentGradientB.z;

      ImproperDihedralConstraints[CurrentSystem][m][2]->Velocity.x-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.x;
      ImproperDihedralConstraints[CurrentSystem][m][2]->Velocity.y-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.y;
      ImproperDihedralConstraints[CurrentSystem][m][2]->Velocity.z-=0.5*DeltaT*InverseMassC*gamma*CurrentGradientC.z;

      ImproperDihedralConstraints[CurrentSystem][m][3]->Velocity.x-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.x;
      ImproperDihedralConstraints[CurrentSystem][m][3]->Velocity.y-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.y;
      ImproperDihedralConstraints[CurrentSystem][m][3]->Velocity.z-=0.5*DeltaT*InverseMassD*gamma*CurrentGradientD.z;

      max_error=MAX2(fabs(gamma),max_error);
    }

    NumberOfRattleSteps++;
  }while(max_error>(1e-8/DeltaT));
  //exit(0);

  NumberOfRattleCyclesStage2[CurrentSystem]+=NumberOfRattleSteps;
  if(NumberOfRattleSteps>MaximumNumberOfRattleCyclesStage2[CurrentSystem])
    MaximumNumberOfRattleCyclesStage2[CurrentSystem]=NumberOfRattleSteps;
}


//
REAL_MATRIX3x3 ComputeRotationMatrix(VECTOR p)
{
  REAL EulerAngle;
  REAL_MATRIX3x3 RotationMatrix;
  REAL CosTheta,SinTheta,fac;

  EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));
  CosTheta=cos(EulerAngle);
  SinTheta=sin(EulerAngle);


  if(EulerAngle<1e-8)
  {
    RotationMatrix.ax=1.0; RotationMatrix.bx=0.0; RotationMatrix.cx=0.0;
    RotationMatrix.ay=0.0; RotationMatrix.by=1.0; RotationMatrix.cy=0.0;
    RotationMatrix.az=0.0; RotationMatrix.bz=0.0; RotationMatrix.cz=1.0;
  }
  else
  {
    fac=1.0/SQR(EulerAngle);
    RotationMatrix.ax=1.0+fac*(SQR(p.y)+SQR(p.z))*(CosTheta-1.0);
    RotationMatrix.ay=fac*(p.x*p.y-p.x*p.y*CosTheta+p.z*EulerAngle*SinTheta);
    RotationMatrix.az=fac*(p.x*p.z-p.x*p.z*CosTheta-p.y*EulerAngle*SinTheta);

    RotationMatrix.bx=fac*(p.x*p.y-p.x*p.y*CosTheta-p.z*EulerAngle*SinTheta);
    RotationMatrix.by=1.0+fac*((SQR(p.x)+SQR(p.z))*(CosTheta-1.0));
    RotationMatrix.bz=fac*(p.y*p.z-p.y*p.z*CosTheta+p.x*EulerAngle*SinTheta);

    RotationMatrix.cx=fac*(p.x*p.z-p.x*p.z*CosTheta+p.y*EulerAngle*SinTheta);
    RotationMatrix.cy=fac*(p.y*p.z-p.y*p.z*CosTheta-p.x*EulerAngle*SinTheta);
    RotationMatrix.cz=1.0+fac*((SQR(p.x)+SQR(p.y))*(CosTheta-1.0));
  }
  return RotationMatrix;
}

REAL_MATRIX3x3 ComputeRotationMatrixDerivativeX(VECTOR p)
{
  REAL EulerAngle;
  REAL_MATRIX3x3 DRotationMatrix;
  REAL CosTheta,SinTheta,fac;

  EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));
  CosTheta=cos(EulerAngle);
  SinTheta=sin(EulerAngle);

  if(EulerAngle<1e-8)
  {
    DRotationMatrix.ax=0.0; DRotationMatrix.bx=0.0; DRotationMatrix.cx=0.0;
    DRotationMatrix.ay=0.0; DRotationMatrix.by=0.0; DRotationMatrix.cy=-1.0;
    DRotationMatrix.az=0.0; DRotationMatrix.bz=1.0; DRotationMatrix.cz=0.0;
  }
  else
  {
    fac=1.0/(SQR(EulerAngle)*SQR(EulerAngle));
    DRotationMatrix.ax=-fac*(p.x*(SQR(p.y)+SQR(p.z))*(-2.0+2.0*CosTheta+EulerAngle*SinTheta));
    DRotationMatrix.ay=fac*(p.y*(-SQR(p.x)+SQR(p.y)+SQR(p.z))+(SQR(p.x)*p.y+SQR(p.x)*p.x*p.z-p.y*(SQR(p.y)+SQR(p.z))+p.x*p.z*(SQR(p.y)+SQR(p.z)))*CosTheta+
                         p.x*(p.x*p.y-p.z)*EulerAngle*SinTheta);
    DRotationMatrix.az=fac*(p.z*(-SQR(p.x)+SQR(p.y)+SQR(p.z))-(SQR(p.x)*p.x*p.y-SQR(p.x)*p.z+p.x*p.y*(SQR(p.y)+SQR(p.z))+p.z*(SQR(p.y)+SQR(p.z)))*CosTheta+
                         p.x*(p.y+p.x*p.z)*EulerAngle*SinTheta);

    DRotationMatrix.bx=fac*(p.y*(-SQR(p.x)+SQR(p.y)+SQR(p.z))-(-SQR(p.x)*p.y+SQR(p.x)*p.x*p.z+p.y*(SQR(p.y)+SQR(p.z))+p.x*p.z*(SQR(p.y)+SQR(p.z)))*CosTheta+
                         p.x*(p.x*p.y+p.z)*EulerAngle*SinTheta);
    DRotationMatrix.by=-fac*(p.x*(2.0*SQR(p.y)-2.0*SQR(p.y)*CosTheta+(SQR(p.x)+SQR(p.z))*EulerAngle*SinTheta));
    DRotationMatrix.bz=fac*(-2.0*p.x*p.y*p.z+p.x*(SQR(p.x)*p.x+2.0*p.y*p.z+p.x*(SQR(p.y)+SQR(p.z)))*CosTheta+EulerAngle*(SQR(p.y)+p.x*p.y*p.z+SQR(p.z))*SinTheta);

    DRotationMatrix.cx=fac*(p.z*(-SQR(p.x)+SQR(p.y)+SQR(p.z))+(SQR(p.x)*p.x*p.y+SQR(p.x)*p.z+p.x*p.y*(SQR(p.y)+SQR(p.z))-p.z*(SQR(p.y)+SQR(p.z)))*CosTheta+
                         p.x*(-p.y+p.x*p.z)*EulerAngle*SinTheta);
    DRotationMatrix.cy=fac*(-2.0*p.x*p.y*p.z-p.x*(SQR(p.x)*p.x-2.0*p.y*p.z+p.x*(SQR(p.y)+SQR(p.z)))*CosTheta-EulerAngle*(SQR(p.y)-p.x*p.y*p.z+SQR(p.z))*SinTheta);
    DRotationMatrix.cz=-fac*(p.x*(2*SQR(p.z) - 2*SQR(p.z)*CosTheta + (SQR(p.x)+SQR(p.y))*EulerAngle*SinTheta));
  }
  return DRotationMatrix;
}

REAL_MATRIX3x3 ComputeRotationMatrixDerivativeY(VECTOR p)
{
  REAL EulerAngle;
  REAL_MATRIX3x3 DRotationMatrix;
  REAL CosTheta,SinTheta,fac;

  EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));
  CosTheta=cos(EulerAngle);
  SinTheta=sin(EulerAngle);

  if(EulerAngle<1e-8)
  {
    DRotationMatrix.ax=0.0;  DRotationMatrix.bx=0.0; DRotationMatrix.cx=1.0;
    DRotationMatrix.ay=0.0;  DRotationMatrix.by=0.0; DRotationMatrix.cy=0.0;
    DRotationMatrix.az=-1.0; DRotationMatrix.bz=0.0; DRotationMatrix.cz=0.0;
  }
  else
  {
    fac=1.0/(SQR(EulerAngle)*SQR(EulerAngle));
    DRotationMatrix.ax=fac*(2.0*SQR(p.x)*p.y*(-1.0+CosTheta)-EulerAngle*p.y*(SQR(p.y)+SQR(p.z))*SinTheta);
    DRotationMatrix.ay=fac*(p.x*(SQR(p.x)-SQR(p.y)+SQR(p.z))+(-CUBE(p.x)+SQR(EulerAngle)*p.y*p.z+p.x*(SQR(p.y)-SQR(p.z)))*CosTheta+EulerAngle*p.y*(p.x*p.y-p.z)*SinTheta);
    DRotationMatrix.az=-fac*(2.0*p.x*p.y*p.z + p.y*(SQR(EulerAngle)*p.y-2.0*p.x*p.z)*CosTheta+EulerAngle*(SQR(p.x)-p.x*p.y*p.z+SQR(p.z))*SinTheta);

    DRotationMatrix.bx=fac*(p.x*(SQR(p.x)-SQR(p.y)+SQR(p.z))-(CUBE(p.x)+SQR(EulerAngle)*p.y*p.z+p.x*(-SQR(p.y)+SQR(p.z)))*CosTheta+EulerAngle*p.y*(p.x*p.y+p.z)*SinTheta);
    DRotationMatrix.by=-fac*(p.y*(SQR(p.x)+SQR(p.z))*(-2.0 + 2.0*CosTheta + EulerAngle*SinTheta));
    DRotationMatrix.bz=fac*(p.z*(SQR(p.x)-SQR(p.y)+SQR(p.z))+(SQR(EulerAngle)*p.x*p.y-p.z*(SQR(p.x)-SQR(p.y)+SQR(p.z)))*CosTheta+EulerAngle*p.y*(-p.x+p.y*p.z)*SinTheta);

    DRotationMatrix.cx=fac*(-2.0*p.x*p.y*p.z+p.y*(SQR(EulerAngle)*p.y+2*p.x*p.z)*CosTheta+EulerAngle*(SQR(p.x)+p.x*p.y*p.z+SQR(p.z))*SinTheta);
    DRotationMatrix.cy=fac*(p.z*(SQR(p.x)-SQR(p.y)+SQR(p.z))-(SQR(EulerAngle)*p.x*p.y+p.z*(SQR(p.x)-SQR(p.y)+SQR(p.z)))*CosTheta+EulerAngle*p.y*(p.x+p.y*p.z)*SinTheta);
    DRotationMatrix.cz=fac*(2.0*p.y*SQR(p.z)*(-1.0+CosTheta)-EulerAngle*p.y*(SQR(p.x)+SQR(p.y))*SinTheta);
  }
  return DRotationMatrix;
}

REAL_MATRIX3x3 ComputeRotationMatrixDerivativeZ(VECTOR p)
{
  REAL EulerAngle;
  REAL_MATRIX3x3 DRotationMatrix;
  REAL CosTheta,SinTheta,fac;

  EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));
  CosTheta=cos(EulerAngle);
  SinTheta=sin(EulerAngle);

  if(EulerAngle<1e-8)
  {
    DRotationMatrix.ax=0.0; DRotationMatrix.bx=-1.0; DRotationMatrix.cx=0.0;
    DRotationMatrix.ay=1.0; DRotationMatrix.by=0.0;  DRotationMatrix.cy=0.0;
    DRotationMatrix.az=0.0; DRotationMatrix.bz=0.0;  DRotationMatrix.cz=0.0;
  }
  else
  {
    fac=1.0/(SQR(EulerAngle)*SQR(EulerAngle));
    DRotationMatrix.ax=fac*(2.0*SQR(p.x)*p.z*(-1.0+CosTheta)-EulerAngle*p.z*(SQR(p.y)+SQR(p.z))*SinTheta);
    DRotationMatrix.ay=fac*(-2.0*p.x*p.y*p.z+p.z*(2.0*p.x*p.y+SQR(EulerAngle)*p.z)*CosTheta+EulerAngle*(SQR(p.x)+SQR(p.y)+p.x*p.y*p.z)*SinTheta);
    DRotationMatrix.az=fac*(p.x*(SQR(p.x)+SQR(p.y)-SQR(p.z))-(CUBE(p.x)+SQR(EulerAngle)*p.y*p.z+p.x*(SQR(p.y)-SQR(p.z)))*CosTheta+EulerAngle*p.z*(p.y+p.x*p.z)*SinTheta);

    DRotationMatrix.bx=fac*(-2.0*p.x*p.y*p.z+p.z*(2.0*p.x*p.y-SQR(EulerAngle)*p.z)*CosTheta-EulerAngle*(SQR(p.x)+SQR(p.y)-p.x*p.y*p.z)*SinTheta);
    DRotationMatrix.by=fac*(2.0*SQR(p.y)*p.z*(-1.0+CosTheta)-EulerAngle*p.z*(SQR(p.x)+SQR(p.z))*SinTheta);
    DRotationMatrix.bz=fac*(p.y*(SQR(p.x)+SQR(p.y)-SQR(p.z))+(-(p.y*(SQR(p.x)+SQR(p.y)))+SQR(EulerAngle)*p.x*p.z+p.y*SQR(p.z))*CosTheta+EulerAngle*p.z*(-p.x+p.y*p.z)*SinTheta);

    DRotationMatrix.cx=fac*(p.x*(SQR(p.x)+SQR(p.y)-SQR(p.z))+(-(p.x*(SQR(p.x)+SQR(p.y)))+SQR(EulerAngle)*p.y*p.z+p.x*SQR(p.z))*CosTheta+EulerAngle*p.z*(-p.y+p.x*p.z)*SinTheta);
    DRotationMatrix.cy=fac*(p.y*(SQR(p.x)+SQR(p.y)-SQR(p.z))-(SQR(p.x)*p.y+CUBE(p.y)+SQR(EulerAngle)*p.x*p.z-p.y*SQR(p.z))*CosTheta+EulerAngle*p.z*(p.x+p.y*p.z)*SinTheta);
    DRotationMatrix.cz=-fac*((SQR(p.x)+SQR(p.y))*p.z*(-2.0+2.0*CosTheta+EulerAngle*SinTheta));
  }
  return DRotationMatrix;
}

REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeAX(VECTOR p)
{
  REAL EulerAngle;
  REAL_MATRIX3x3 DDRotationMatrix;
  REAL CosTheta,SinTheta,fac;

  EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));
  CosTheta=cos(EulerAngle);
  SinTheta=sin(EulerAngle);

  if(EulerAngle<1e-8)
  {
    DDRotationMatrix.ax=0.0; DDRotationMatrix.bx=0.0;  DDRotationMatrix.cx=0.0;
    DDRotationMatrix.ay=0.0; DDRotationMatrix.by=-1.0; DDRotationMatrix.cy=0.0;
    DDRotationMatrix.az=0.0; DDRotationMatrix.bz=0.0;  DDRotationMatrix.cz=-1.0;
  }
  else
  {
    fac=1.0/pow(EulerAngle,6.0);
    DDRotationMatrix.ax=-fac*((SQR(p.y)+SQR(p.z))*(-2.0*SQR(EulerAngle)+8.0*SQR(p.x)+(-8.0*SQR(p.x)+SQR(EulerAngle)*(2.0+SQR(p.x)))*CosTheta+EulerAngle*(SQR(EulerAngle)-5.0*SQR(p.x))*SinTheta));
    DDRotationMatrix.ay=fac*(-6.0*SQR(EulerAngle)*p.x*p.y+8.0*CUBE(p.x)*p.y+(-8.0*CUBE(p.x)*p.y+SQR(EulerAngle)*SQR(EulerAngle)*p.z+SQR(EulerAngle)*p.x*((6.0+SQR(p.x))*p.y-3.0*p.x*p.z))*CosTheta-
                        EulerAngle*(SQR(p.x)*(5.0*p.x*p.y-3.0*p.z)+SQR(EulerAngle)*(-3.0*p.x*p.y+p.z+SQR(p.x)*p.z))*SinTheta);
    DDRotationMatrix.az=fac*(-6.0*SQR(EulerAngle)*p.x*p.z+8.0*CUBE(p.x)*p.z+(-(SQR(EulerAngle)*SQR(EulerAngle)*p.y)-8.0*CUBE(p.x)*p.z+SQR(EulerAngle)*p.x*(3.0*p.x*p.y+(6.0+SQR(p.x))*p.z))*CosTheta+
                        EulerAngle*(SQR(EulerAngle)*(p.y+SQR(p.x)*p.y+3.0*p.x*p.z)-SQR(p.x)*(3.0*p.y+5.0*p.x*p.z))*SinTheta);

    DDRotationMatrix.bx=fac*(-6.0*SQR(EulerAngle)*p.x*p.y+8.0*CUBE(p.x)*p.y+(-8.0*CUBE(p.x)*p.y-SQR(EulerAngle)*SQR(EulerAngle)*p.z+SQR(EulerAngle)*p.x*((6.0+SQR(p.x))*p.y+3.0*p.x*p.z))*CosTheta+
                        EulerAngle*(-(SQR(p.x)*(5.0*p.x*p.y+3.0*p.z))+SQR(EulerAngle)*(3.0*p.x*p.y+p.z+SQR(p.x)*p.z))*SinTheta);
    DDRotationMatrix.by=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.x))*(SQR(EulerAngle)-SQR(p.x)-SQR(p.z))+
                        (2.0*SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.x)*(SQR(p.x)+SQR(p.z))-SQR(EulerAngle)*(SQR(p.x)*SQR(p.x)+2.0*SQR(p.z)+SQR(p.x)*(10.0+SQR(p.z))))*CosTheta+
                        EulerAngle*(5.0*SQR(p.x)*(SQR(p.x)+SQR(p.z))-SQR(EulerAngle)*(5.0*SQR(p.x)+SQR(p.z)))*SinTheta);
    DDRotationMatrix.bz=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.x))*p.y*p.z+(3.0*SQR(EulerAngle)*SQR(EulerAngle)*p.x-8.0*SQR(p.x)*p.y*p.z+SQR(EulerAngle)*(-3.0*CUBE(p.x)+(2.0+SQR(p.x))*p.y*p.z))*CosTheta+
                        EulerAngle*(SQR(p.x)*(3.0*p.x-5.0*p.y*p.z)-SQR(EulerAngle)*(3.0*p.x+CUBE(p.x)-p.y*p.z))*SinTheta);

    DDRotationMatrix.cx=fac*(-6.0*SQR(EulerAngle)*p.x*p.z+8.0*CUBE(p.x)*p.z+(SQR(EulerAngle)*SQR(EulerAngle)*p.y-8.0*CUBE(p.x)*p.z+SQR(EulerAngle)*p.x*(-3.0*p.x*p.y+(6.0+SQR(p.x))*p.z))*CosTheta-
                         EulerAngle*(SQR(EulerAngle)*(p.y+SQR(p.x)*p.y-3.0*p.x*p.z)+SQR(p.x)*(-3.0*p.y+5.0*p.x*p.z))*SinTheta);
    DDRotationMatrix.cy=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.x))*p.y*p.z+(-3.0*SQR(EulerAngle)*SQR(EulerAngle)*p.x-8.0*SQR(p.x)*p.y*p.z+SQR(EulerAngle)*(3.0*CUBE(p.x)+(2.0+SQR(p.x))*p.y*p.z))*CosTheta+
                        EulerAngle*(SQR(EulerAngle)*(3.0*p.x+CUBE(p.x)+p.y*p.z)-SQR(p.x)*(3.0*p.x+5.0*p.y*p.z))*SinTheta);
    DDRotationMatrix.cz=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.x))*(SQR(EulerAngle)-SQR(p.x)-SQR(p.y))+
                        (2.0*SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.x)*(SQR(p.x)+SQR(p.y))-SQR(EulerAngle)*(SQR(p.x)*SQR(p.x)+2.0*SQR(p.y)+SQR(p.x)*(10.0+SQR(p.y))))*CosTheta+
                        EulerAngle*(5.0*SQR(p.x)*(SQR(p.x)+SQR(p.y))-SQR(EulerAngle)*(5.0*SQR(p.x)+SQR(p.y)))*SinTheta);
  }

  return DDRotationMatrix;
}

REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeBY(VECTOR p)
{
  REAL EulerAngle;
  REAL_MATRIX3x3 DDRotationMatrix;
  REAL CosTheta,SinTheta,fac;

  EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));
  CosTheta=cos(EulerAngle);
  SinTheta=sin(EulerAngle);

  if(EulerAngle<1e-8)
  {
    DDRotationMatrix.ax=-1.0; DDRotationMatrix.bx=0.0; DDRotationMatrix.cx=0.0;
    DDRotationMatrix.ay=0.0;  DDRotationMatrix.by=0.0; DDRotationMatrix.cy=0.0;
    DDRotationMatrix.az=0.0;  DDRotationMatrix.bz=0.0; DDRotationMatrix.cz=-1.0;
  }
  else
  {
    fac=1.0/pow(EulerAngle,6.0);
    DDRotationMatrix.ax=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.y))*(SQR(EulerAngle)-SQR(p.y)-SQR(p.z))+
                        (2*SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.y)*(SQR(p.y)+SQR(p.z))-SQR(EulerAngle)*(SQR(p.y)*SQR(p.y)+2.0*SQR(p.z)+SQR(p.y)*(10.0+SQR(p.z))))*CosTheta+
                         EulerAngle*(5*SQR(p.y)*(SQR(p.y)+SQR(p.z))-SQR(EulerAngle)*(5.0*SQR(p.y)+SQR(p.z)))*SinTheta);
    DDRotationMatrix.ay=fac*(-6.0*SQR(EulerAngle)*p.x*p.y+8.0*p.x*CUBE(p.y)+(-8.0*p.x*CUBE(p.y)+SQR(EulerAngle)*SQR(EulerAngle)*p.z+SQR(EulerAngle)*p.y*(p.x*(6.0+SQR(p.y))-3.0*p.y*p.z))*CosTheta+
                        EulerAngle*(SQR(p.y)*(-5.0*p.x*p.y+3.0*p.z)+SQR(EulerAngle)*(3.0*p.x*p.y-(1.0+SQR(p.y))*p.z))*SinTheta);
    DDRotationMatrix.az=fac*(-2.0*p.x*(SQR(EulerAngle)-4.0*SQR(p.y))*p.z+(-3.0*SQR(EulerAngle)*SQR(EulerAngle)*p.y-8.0*p.x*SQR(p.y)*p.z+SQR(EulerAngle)*(3.0*CUBE(p.y)+p.x*(2.0+SQR(p.y))*p.z))*CosTheta+
                        EulerAngle*(SQR(EulerAngle)*(3.0*p.y+CUBE(p.y)+p.x*p.z)-SQR(p.y)*(3.0*p.y+5.0*p.x*p.z))*SinTheta);

    DDRotationMatrix.bx=fac*(-6.0*SQR(EulerAngle)*p.x*p.y+8.0*p.x*CUBE(p.y)+(-8.0*p.x*CUBE(p.y)-SQR(EulerAngle)*SQR(EulerAngle)*p.z+SQR(EulerAngle)*p.y*(p.x*(6.0+SQR(p.y))+3.0*p.y*p.z))*CosTheta+
                         EulerAngle*(-(SQR(p.y)*(5.0*p.x*p.y+3.0*p.z))+SQR(EulerAngle)*(3.0*p.x*p.y+p.z+SQR(p.y)*p.z))*SinTheta);
    DDRotationMatrix.by=-fac*((SQR(p.x)+SQR(p.z))*(-2.0*SQR(EulerAngle)+8.0*SQR(p.y)+(-8.0*SQR(p.y)+SQR(EulerAngle)*(2.0+SQR(p.y)))*CosTheta+EulerAngle*(SQR(EulerAngle)-5.0*SQR(p.y))*SinTheta));
    DDRotationMatrix.bz=fac*(-6.0*SQR(EulerAngle)*p.y*p.z+8.0*CUBE(p.y)*p.z+(SQR(EulerAngle)*SQR(EulerAngle)*p.x-8.0*CUBE(p.y)*p.z+SQR(EulerAngle)*p.y*(-3.0*p.x*p.y+(6.0+SQR(p.y))*p.z))*CosTheta-
                         EulerAngle*(SQR(EulerAngle)*(p.x+p.x*SQR(p.y)-3.0*p.y*p.z)+SQR(p.y)*(-3.0*p.x+5.0*p.y*p.z))*SinTheta);

    DDRotationMatrix.cx=fac*(-2.0*p.x*(SQR(EulerAngle)-4.0*SQR(p.y))*p.z+(3.0*SQR(EulerAngle)*SQR(EulerAngle)*p.y-8.0*p.x*SQR(p.y)*p.z+SQR(EulerAngle)*(-3.0*CUBE(p.y)+p.x*(2.0+SQR(p.y))*p.z))*CosTheta+
                         EulerAngle*(SQR(p.y)*(3.0*p.y-5.0*p.x*p.z)-SQR(EulerAngle)*(3.0*p.y+CUBE(p.y)-p.x*p.z))*SinTheta);
    DDRotationMatrix.cy=fac*(-6.0*SQR(EulerAngle)*p.y*p.z+8.0*CUBE(p.y)*p.z+(-(SQR(EulerAngle)*SQR(EulerAngle)*p.x)-8.0*CUBE(p.y)*p.z+SQR(EulerAngle)*p.y*(3.0*p.x*p.y+(6.0+SQR(p.y))*p.z))*CosTheta+
                         EulerAngle*(SQR(EulerAngle)*(p.x+p.x*SQR(p.y)+3.0*p.y*p.z)-SQR(p.y)*(3.0*p.x+5.0*p.y*p.z))*SinTheta);
    DDRotationMatrix.cz=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.y))*(SQR(EulerAngle)-SQR(p.x)-SQR(p.y))+
                        (2.0*SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.y)*(SQR(p.x)+SQR(p.y))-SQR(EulerAngle)*(2.0*SQR(p.x)+(10.0+SQR(p.x))*SQR(p.y)+SQR(p.y)*SQR(p.y)))*CosTheta+
                        EulerAngle*(5*SQR(p.y)*(SQR(p.x)+SQR(p.y))-SQR(EulerAngle)*(SQR(p.x)+5.0*SQR(p.y)))*SinTheta);
  }

  return DDRotationMatrix;
}

REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeCZ(VECTOR p)
{
  REAL EulerAngle;
  REAL_MATRIX3x3 DDRotationMatrix;
  REAL CosTheta,SinTheta,fac;

  EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));
  CosTheta=cos(EulerAngle);
  SinTheta=sin(EulerAngle);

  if(EulerAngle<1e-8)
  {
    DDRotationMatrix.ax=-1.0; DDRotationMatrix.bx=0.0;  DDRotationMatrix.cx=0.0;
    DDRotationMatrix.ay=0.0;  DDRotationMatrix.by=-1.0; DDRotationMatrix.cy=0.0;
    DDRotationMatrix.az=0.0;  DDRotationMatrix.bz=0.0;  DDRotationMatrix.cz=0.0;
  }
  else
  {
    fac=1.0/pow(EulerAngle,6.0);
    DDRotationMatrix.ax=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.z))*(SQR(EulerAngle)-SQR(p.y)-SQR(p.z))+
                        (2.0*SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.z)*(SQR(p.y)+SQR(p.z))-SQR(EulerAngle)*(2.0*SQR(p.y)+(10.0+SQR(p.y))*SQR(p.z)+SQR(p.z)*SQR(p.z)))*CosTheta +
                        EulerAngle*(5.0*SQR(p.z)*(SQR(p.y)+SQR(p.z))-SQR(EulerAngle)*(SQR(p.y)+5.0*SQR(p.z)))*SinTheta);
    DDRotationMatrix.ay=fac*(-2.0*p.x*p.y*(SQR(EulerAngle)-4.0*SQR(p.z))+(3.0*SQR(EulerAngle)*SQR(EulerAngle)*p.z-8.0*p.x*p.y*SQR(p.z)+SQR(EulerAngle)*(-3.0*CUBE(p.z)+p.x*p.y*(2.0+SQR(p.z))))*CosTheta +
                        EulerAngle*(SQR(p.z)*(-5.0*p.x*p.y+3.0*p.z)+SQR(EulerAngle)*(p.x*p.y-p.z*(3.0+SQR(p.z))))*SinTheta);
    DDRotationMatrix.az=fac*(-6.0*SQR(EulerAngle)*p.x*p.z+8.0*p.x*CUBE(p.z)+(-(SQR(EulerAngle)*SQR(EulerAngle)*p.y)-8.0*p.x*CUBE(p.z)+SQR(EulerAngle)*p.z*(3.0*p.y*p.z+p.x*(6.0+SQR(p.z))))*CosTheta +
                        EulerAngle*(-(SQR(p.z)*(3.0*p.y+5.0*p.x*p.z))+SQR(EulerAngle)*(p.y+3.0*p.x*p.z+p.y*SQR(p.z)))*SinTheta);

    DDRotationMatrix.bx=fac*(-2.0*p.x*p.y*(SQR(EulerAngle)-4.0*SQR(p.z))+(-3.0*SQR(EulerAngle)*SQR(EulerAngle)*p.z-8.0*p.x*p.y*SQR(p.z)+SQR(EulerAngle)*(3.0*CUBE(p.z)+p.x*p.y*(2.0+SQR(p.z))))*CosTheta +
                        EulerAngle*(-(SQR(p.z)*(5.0*p.x*p.y+3.0*p.z))+SQR(EulerAngle)*(p.x*p.y+3.0*p.z+CUBE(p.z)))*SinTheta);
    DDRotationMatrix.by=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.z))*(SQR(EulerAngle)-SQR(p.x)-SQR(p.z))+
                        (2.0*SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.z)*(SQR(p.x)+SQR(p.z))-SQR(EulerAngle)*(2.0*SQR(p.x)+(10.0+SQR(p.x))*SQR(p.z)+SQR(p.z)*SQR(p.z)))*CosTheta +
                        EulerAngle*(5.0*SQR(p.z)*(SQR(p.x)+SQR(p.z))-SQR(EulerAngle)*(SQR(p.x)+5.0*SQR(p.z)))*SinTheta);
    DDRotationMatrix.bz=fac*(-6.0*SQR(EulerAngle)*p.y*p.z+8.0*p.y*CUBE(p.z)+(SQR(EulerAngle)*SQR(EulerAngle)*p.x-8.0*p.y*CUBE(p.z)+SQR(EulerAngle)*p.z*(-3.0*p.x*p.z+p.y*(6.0+SQR(p.z))))*CosTheta -
                        EulerAngle*(SQR(p.z)*(-3.0*p.x+5.0*p.y*p.z)+SQR(EulerAngle)*(p.x-3.0*p.y*p.z+p.x*SQR(p.z)))*SinTheta);

    DDRotationMatrix.cx=fac*(-6.0*SQR(EulerAngle)*p.x*p.z+8.0*p.x*CUBE(p.z)+(SQR(EulerAngle)*SQR(EulerAngle)*p.y-8.0*p.x*CUBE(p.z)+SQR(EulerAngle)*p.z*(-3.0*p.y*p.z+p.x*(6.0+SQR(p.z))))*CosTheta -
                        EulerAngle*(SQR(p.z)*(-3.0*p.y+5.0*p.x*p.z)+SQR(EulerAngle)*(p.y-3.0*p.x*p.z+p.y*SQR(p.z)))*SinTheta);
    DDRotationMatrix.cy=fac*(-6.0*SQR(EulerAngle)*p.y*p.z+8.0*p.y*CUBE(p.z)+(-(SQR(EulerAngle)*SQR(EulerAngle)*p.x)-8.0*p.y*CUBE(p.z)+SQR(EulerAngle)*p.z*(3.0*p.x*p.z+p.y*(6.0+SQR(p.z))))*CosTheta +
                        EulerAngle*(-(SQR(p.z)*(3.0*p.x+5.0*p.y*p.z))+SQR(EulerAngle)*(p.x+3.0*p.y*p.z+p.x*SQR(p.z)))*SinTheta);
    DDRotationMatrix.cz=-fac*((SQR(p.x)+SQR(p.y))*(-2.0*SQR(EulerAngle)+8.0*SQR(p.z)+(-8.0*SQR(p.z)+SQR(EulerAngle)*(2.0+SQR(p.z)))*CosTheta+EulerAngle*(SQR(EulerAngle)-5.0*SQR(p.z))*SinTheta));
  }

  return DDRotationMatrix;
}

REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeAY(VECTOR p)
{
  REAL EulerAngle;
  REAL_MATRIX3x3 DDRotationMatrix;
  REAL CosTheta,SinTheta,fac;

  EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));
  CosTheta=cos(EulerAngle);
  SinTheta=sin(EulerAngle);

  if(EulerAngle<1e-8)
  {
    DDRotationMatrix.ax=0.0; DDRotationMatrix.bx=0.5; DDRotationMatrix.cx=0.0;
    DDRotationMatrix.ay=0.5; DDRotationMatrix.by=0.0; DDRotationMatrix.cy=0.0;
    DDRotationMatrix.az=0.0; DDRotationMatrix.bz=0.0; DDRotationMatrix.cz=0.0;
  }
  else
  {
    fac=1.0/pow(EulerAngle,6.0);
    DDRotationMatrix.ax=-fac*(p.x*p.y*(-4.0*SQR(EulerAngle)+8.0*(SQR(p.y)+SQR(p.z))+(-8.0*(SQR(p.y)+SQR(p.z))+SQR(EulerAngle)*(4.0+SQR(p.y)+SQR(p.z)))*CosTheta+
                        EulerAngle*(2.0*SQR(EulerAngle)-5.0*(SQR(p.y)+SQR(p.z)))*SinTheta));
    DDRotationMatrix.ay=fac*(SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.x)*SQR(p.y)-2.0*SQR(EulerAngle)*(SQR(p.x)+SQR(p.y))+
                        (-SQR(EulerAngle)*SQR(EulerAngle)-8.0*SQR(p.x)*SQR(p.y)+SQR(EulerAngle)*(2.0*SQR(p.y)+SQR(p.x)*(2.0+SQR(p.y))-3.0*p.x*p.y*p.z))*CosTheta +
                        EulerAngle*(p.x*p.y*(-5.0*p.x*p.y+3.0*p.z)+SQR(EulerAngle)*(SQR(p.x)+SQR(p.y)-p.x*p.y*p.z))*SinTheta);
    DDRotationMatrix.az=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.x))*p.y*p.z+(-(SQR(EulerAngle)*SQR(EulerAngle)*p.x)-8.0*SQR(p.x)*p.y*p.z+SQR(EulerAngle)*p.y*(3.0*p.x*p.y+(2.0+SQR(p.x))*p.z))*CosTheta +
                        EulerAngle*(-(p.x*p.y*(3.0*p.y+5.0*p.x*p.z))+SQR(EulerAngle)*(p.x+p.x*SQR(p.y)+p.y*p.z))*SinTheta);

    DDRotationMatrix.bx=fac*(SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.x)*SQR(p.y)-2.0*SQR(EulerAngle)*(SQR(p.x)+SQR(p.y))+
                        (-SQR(EulerAngle)*SQR(EulerAngle)-8.0*SQR(p.x)*SQR(p.y)+SQR(EulerAngle)*(2.0*SQR(p.y)+SQR(p.x)*(2.0+SQR(p.y))+3.0*p.x*p.y*p.z))*CosTheta +
                        EulerAngle*(-(p.x*p.y*(5.0*p.x*p.y+3.0*p.z))+SQR(EulerAngle)*(SQR(p.x)+SQR(p.y)+p.x*p.y*p.z))*SinTheta);
    DDRotationMatrix.by=-fac*(p.x*p.y*(-4.0*SQR(EulerAngle)+8.0*(SQR(p.x)+SQR(p.z))+(-8.0*(SQR(p.x)+SQR(p.z))+SQR(EulerAngle)*(4.0+SQR(p.x)+SQR(p.z)))*CosTheta +
                        EulerAngle*(2.0*SQR(EulerAngle)-5.0*(SQR(p.x)+SQR(p.z)))*SinTheta));
    DDRotationMatrix.bz=fac*(-2.0*p.x*(SQR(EulerAngle)-4.0*SQR(p.y))*p.z+(SQR(EulerAngle)*SQR(EulerAngle)*p.y-8.0*p.x*SQR(p.y)*p.z+SQR(EulerAngle)*p.x*(-3.0*p.x*p.y+(2.0+SQR(p.y))*p.z))*CosTheta -
                        EulerAngle*(SQR(EulerAngle)*(p.y+SQR(p.x)*p.y-p.x*p.z)+p.x*p.y*(-3.0*p.x+5.0*p.y*p.z))*SinTheta);

    DDRotationMatrix.cx=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.x))*p.y*p.z+(SQR(EulerAngle)*SQR(EulerAngle)*p.x-8.0*SQR(p.x)*p.y*p.z+SQR(EulerAngle)*p.y*(-3.0*p.x*p.y+(2.0+SQR(p.x))*p.z))*CosTheta-
                        EulerAngle*(p.x*p.y*(-3.0*p.y+5.0*p.x*p.z)+SQR(EulerAngle)*(p.x+p.x*SQR(p.y)-p.y*p.z))*SinTheta);
    DDRotationMatrix.cy=fac*(-2.0*p.x*(SQR(EulerAngle)-4.0*SQR(p.y))*p.z+(-(SQR(EulerAngle)*SQR(EulerAngle)*p.y)-8.0*p.x*SQR(p.y)*p.z+SQR(EulerAngle)*p.x*(3.0*p.x*p.y+(2.0+SQR(p.y))*p.z))*CosTheta+
                        EulerAngle*(SQR(EulerAngle)*(p.y+SQR(p.x)*p.y+p.x*p.z)-p.x*p.y*(3.0*p.x+5.0*p.y*p.z))*SinTheta);
    DDRotationMatrix.cz=-fac*(p.x*p.y*(8.0*(-SQR(EulerAngle)+SQR(p.x)+SQR(p.y))+(-8.0*(SQR(p.x)+SQR(p.y))+SQR(EulerAngle)*(8.0+SQR(p.x)+SQR(p.y)))*CosTheta +
                        EulerAngle*(4.0*SQR(EulerAngle)-5.0*(SQR(p.x)+SQR(p.y)))*SinTheta));
  }

  return DDRotationMatrix;
}

REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeAZ(VECTOR p)
{
  REAL EulerAngle;
  REAL_MATRIX3x3 DDRotationMatrix;
  REAL CosTheta,SinTheta,fac;

  EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));
  CosTheta=cos(EulerAngle);
  SinTheta=sin(EulerAngle);

  if(EulerAngle<1e-8)
  {
    DDRotationMatrix.ax=0.0; DDRotationMatrix.bx=0.0; DDRotationMatrix.cx=0.5;
    DDRotationMatrix.ay=0.0; DDRotationMatrix.by=0.0; DDRotationMatrix.cy=0.0;
    DDRotationMatrix.az=0.5; DDRotationMatrix.bz=0.0; DDRotationMatrix.cz=0.0;
  }
  else
  {
    fac=1.0/pow(EulerAngle,6.0);
    DDRotationMatrix.ax=-fac*(p.x*p.z*(-4.0*SQR(EulerAngle)+8.0*(SQR(p.y)+SQR(p.z))+(-8.0*(SQR(p.y)+SQR(p.z))+SQR(EulerAngle)*(4.0+SQR(p.y)+SQR(p.z)))*CosTheta+
                        EulerAngle*(2.0*SQR(EulerAngle)-5.0*(SQR(p.y)+SQR(p.z)))*SinTheta));
    DDRotationMatrix.ay=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.x))*p.y*p.z+(SQR(EulerAngle)*SQR(EulerAngle)*p.x-8.0*SQR(p.x)*p.y*p.z+SQR(EulerAngle)*p.z*((2.0+SQR(p.x))*p.y-3.0*p.x*p.z))*CosTheta-
                        EulerAngle*(p.x*(5.0*p.x*p.y-3.0*p.z)*p.z+SQR(EulerAngle)*(p.x-p.y*p.z+p.x*SQR(p.z)))*SinTheta);
    DDRotationMatrix.az=fac*(SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.x)*SQR(p.z)-2.0*SQR(EulerAngle)*(SQR(p.x)+SQR(p.z))+
                        (-SQR(EulerAngle)*SQR(EulerAngle)-8.0*SQR(p.x)*SQR(p.z)+SQR(EulerAngle)*(3.0*p.x*p.y*p.z+2.0*SQR(p.z)+SQR(p.x)*(2.0+SQR(p.z))))*CosTheta+
                        EulerAngle*(-(p.x*p.z*(3.0*p.y+5.0*p.x*p.z))+SQR(EulerAngle)*(SQR(p.x)+p.x*p.y*p.z+SQR(p.z)))*SinTheta);

    DDRotationMatrix.bx=fac*(-2.0*(SQR(EulerAngle)-4.0*SQR(p.x))*p.y*p.z+(-(SQR(EulerAngle)*SQR(EulerAngle)*p.x)-8.0*SQR(p.x)*p.y*p.z+SQR(EulerAngle)*p.z*((2.0+SQR(p.x))*p.y+3.0*p.x*p.z))*CosTheta +
                        EulerAngle*(-(p.x*p.z*(5.0*p.x*p.y+3.0*p.z))+SQR(EulerAngle)*(p.x+p.y*p.z+p.x*SQR(p.z)))*SinTheta);
    DDRotationMatrix.by=-fac*(p.x*p.z*(8.0*(-SQR(EulerAngle)+SQR(p.x)+SQR(p.z))+(-8.0*(SQR(p.x)+SQR(p.z))+SQR(EulerAngle)*(8.0+SQR(p.x)+SQR(p.z)))*CosTheta+
                        EulerAngle*(4.0*SQR(EulerAngle)-5.0*(SQR(p.x)+SQR(p.z)))*SinTheta));
    DDRotationMatrix.bz=fac*(-2.0*p.x*p.y*(SQR(EulerAngle)-4.0*SQR(p.z))+(SQR(EulerAngle)*SQR(EulerAngle)*p.z-8.0*p.x*p.y*SQR(p.z)+SQR(EulerAngle)*p.x*(-3.0*p.x*p.z+p.y*(2.0+SQR(p.z))))*CosTheta-
                        EulerAngle*(SQR(EulerAngle)*(-(p.x*p.y)+p.z+SQR(p.x)*p.z)+p.x*p.z*(-3.0*p.x+5.0*p.y*p.z))*SinTheta);

    DDRotationMatrix.cx=fac*(SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.x)*SQR(p.z)-2.0*SQR(EulerAngle)*(SQR(p.x)+SQR(p.z))+
                        (-SQR(EulerAngle)*SQR(EulerAngle)-8.0*SQR(p.x)*SQR(p.z)+SQR(EulerAngle)*(-3.0*p.x*p.y*p.z+2*SQR(p.z)+SQR(p.x)*(2.0+SQR(p.z))))*CosTheta +
                        EulerAngle*(p.x*p.z*(3.0*p.y-5.0*p.x*p.z)+SQR(EulerAngle)*(SQR(p.x)-p.x*p.y*p.z+SQR(p.z)))*SinTheta);
    DDRotationMatrix.cy=fac*(-2.0*p.x*p.y*(SQR(EulerAngle)-4.0*SQR(p.z))+(-(SQR(EulerAngle)*SQR(EulerAngle)*p.z)-8.0*p.x*p.y*SQR(p.z)+SQR(EulerAngle)*p.x*(3.0*p.x*p.z+p.y*(2.0+SQR(p.z))))*CosTheta +
                        EulerAngle*(-(p.x*p.z*(3.0*p.x+5.0*p.y*p.z))+SQR(EulerAngle)*(p.z+p.x*(p.y+p.x*p.z)))*SinTheta);
    DDRotationMatrix.cz=-fac*(p.x*p.z*(-4.0*SQR(EulerAngle)+8.0*(SQR(p.x)+SQR(p.y))+(-8.0*(SQR(p.x)+SQR(p.y))+SQR(EulerAngle)*(4.0+SQR(p.x)+SQR(p.y)))*CosTheta +
                        EulerAngle*(2.0*SQR(EulerAngle)-5.0*(SQR(p.x)+SQR(p.y)))*SinTheta));
  }

  return DDRotationMatrix;
}

REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeBZ(VECTOR p)
{
  REAL EulerAngle;
  REAL_MATRIX3x3 DDRotationMatrix;
  REAL CosTheta,SinTheta,fac;

  EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));
  CosTheta=cos(EulerAngle);
  SinTheta=sin(EulerAngle);

  if(EulerAngle<1e-8)
  {
    DDRotationMatrix.ax=0.0; DDRotationMatrix.bx=0.0; DDRotationMatrix.cx=0.0;
    DDRotationMatrix.ay=0.0; DDRotationMatrix.by=0.0; DDRotationMatrix.cy=0.5;
    DDRotationMatrix.az=0.0; DDRotationMatrix.bz=0.5; DDRotationMatrix.cz=0.0;
  }
  else
  {
    fac=1.0/pow(EulerAngle,6.0);
    DDRotationMatrix.ax=-fac*(p.y*p.z*(8.0*(-SQR(EulerAngle)+SQR(p.y)+SQR(p.z))+(-8.0*(SQR(p.y)+SQR(p.z))+SQR(EulerAngle)*(8.0+SQR(p.y)+SQR(p.z)))*CosTheta+
                         EulerAngle*(4.0*SQR(EulerAngle)-5.0*(SQR(p.y)+SQR(p.z)))*SinTheta));
    DDRotationMatrix.ay=fac*(-2.0*p.x*(SQR(EulerAngle)-4.0*SQR(p.y))*p.z+(SQR(EulerAngle)*SQR(EulerAngle)*p.y-8.0*p.x*SQR(p.y)*p.z+SQR(EulerAngle)*p.z*(p.x*(2.0+SQR(p.y))-3.0*p.y*p.z))*CosTheta-
                         EulerAngle*(p.y*(5.0*p.x*p.y-3.0*p.z)*p.z+SQR(EulerAngle)*(p.y-p.x*p.z+p.y*SQR(p.z)))*SinTheta);
    DDRotationMatrix.az=fac*(-2.0*p.x*p.y*(SQR(EulerAngle)-4.0*SQR(p.z))+(-(SQR(EulerAngle)*SQR(EulerAngle)*p.z)-8*p.x*p.y*SQR(p.z)+SQR(EulerAngle)*p.y*(3.0*p.y*p.z+p.x*(2.0+SQR(p.z))))*CosTheta+
                        EulerAngle*(-(p.y*p.z*(3.0*p.y+5.0*p.x*p.z))+SQR(EulerAngle)*(p.z+p.y*(p.x+p.y*p.z)))*SinTheta);

    DDRotationMatrix.bx=fac*(-2.0*p.x*(SQR(EulerAngle)-4.0*SQR(p.y))*p.z+(-(SQR(EulerAngle)*SQR(EulerAngle)*p.y)-8.0*p.x*SQR(p.y)*p.z+SQR(EulerAngle)*p.z*(p.x*(2.0+SQR(p.y))+3.0*p.y*p.z))*CosTheta+
                         EulerAngle*(-(p.y*p.z*(5.0*p.x*p.y+3.0*p.z))+SQR(EulerAngle)*(p.y+p.x*p.z+p.y*SQR(p.z)))*SinTheta);
    DDRotationMatrix.by=-fac*(p.y*p.z*(-4.0*SQR(EulerAngle)+8.0*(SQR(p.x)+SQR(p.z))+(-8.0*(SQR(p.x)+SQR(p.z))+SQR(EulerAngle)*(4.0+SQR(p.x)+SQR(p.z)))*CosTheta+
                        EulerAngle*(2.0*SQR(EulerAngle)-5.0*(SQR(p.x)+SQR(p.z)))*SinTheta));
    DDRotationMatrix.bz=fac*(SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.y)*SQR(p.z)-2.0*SQR(EulerAngle)*(SQR(p.y)+SQR(p.z))+
                        (-SQR(EulerAngle)*SQR(EulerAngle)-8.0*SQR(p.y)*SQR(p.z)+SQR(EulerAngle)*(-3.0*p.x*p.y*p.z+2.0*SQR(p.z)+SQR(p.y)*(2.0+SQR(p.z))))*CosTheta +
                        EulerAngle*(p.y*p.z*(3.0*p.x-5.0*p.y*p.z)+SQR(EulerAngle)*(SQR(p.y)-p.x*p.y*p.z+SQR(p.z)))*SinTheta);

    DDRotationMatrix.cx=fac*(-2.0*p.x*p.y*(SQR(EulerAngle)-4.0*SQR(p.z))+(SQR(EulerAngle)*SQR(EulerAngle)*p.z-8.0*p.x*p.y*SQR(p.z)+SQR(EulerAngle)*p.y*(-3.0*p.y*p.z+p.x*(2.0+SQR(p.z))))*CosTheta+
                        EulerAngle*(p.y*p.z*(3.0*p.y-5.0*p.x*p.z)+SQR(EulerAngle)*(p.x*p.y-(1.0+SQR(p.y))*p.z))*SinTheta);
    DDRotationMatrix.cy=fac*(SQR(EulerAngle)*SQR(EulerAngle)+8.0*SQR(p.y)*SQR(p.z)-2.0*SQR(EulerAngle)*(SQR(p.y)+SQR(p.z))+
                        (-SQR(EulerAngle)*SQR(EulerAngle)-8.0*SQR(p.y)*SQR(p.z)+SQR(EulerAngle)*(3.0*p.x*p.y*p.z+2.0*SQR(p.z)+SQR(p.y)*(2.0+SQR(p.z))))*CosTheta+
                        EulerAngle*(-(p.y*p.z*(3.0*p.x+5.0*p.y*p.z))+SQR(EulerAngle)*(SQR(p.y)+p.x*p.y*p.z+SQR(p.z)))*SinTheta);
    DDRotationMatrix.cz=-fac*(p.y*p.z*(-4.0*SQR(EulerAngle)+8.0*(SQR(p.x)+SQR(p.y))+(-8.0*(SQR(p.x)+SQR(p.y))+SQR(EulerAngle)*(4.0+SQR(p.x)+SQR(p.y)))*CosTheta+
                        EulerAngle*(2.0*SQR(EulerAngle)-5.0*(SQR(p.x)+SQR(p.y)))*SinTheta));
  }

  return DDRotationMatrix;
}

void CalculateConstraintsExclusionEnergy(void)
{
  int m;
  int typeA,typeB;
  REAL rr,r;
  REAL energy,force_factor;
  REAL chargeA,chargeB;
  VECTOR posA,posB,dr;

  for(m=0;m<NumberOfDistanceConstraints[CurrentSystem];m++)
  {
    posA=DistanceConstraints[CurrentSystem][m][0]->Position;
    posB=DistanceConstraints[CurrentSystem][m][1]->Position;

    typeA=DistanceConstraints[CurrentSystem][m][0]->Type;
    typeB=DistanceConstraints[CurrentSystem][m][1]->Type;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    dr=ApplyBoundaryCondition(dr);
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

    if(rr<CutOffVDWSquared)
    {
      PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);

      UExclusionConstraints[CurrentSystem]-=energy;

    }
  }

  for(m=0;m<NumberOfDistanceConstraints[CurrentSystem];m++)
  {
    posA=DistanceConstraints[CurrentSystem][m][0]->Position;
    posB=DistanceConstraints[CurrentSystem][m][1]->Position;

    typeA=DistanceConstraints[CurrentSystem][m][0]->Type;
    typeB=DistanceConstraints[CurrentSystem][m][1]->Type;

    chargeA=PseudoAtoms[typeA].Charge1;
    chargeB=PseudoAtoms[typeB].Charge1;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    dr=ApplyBoundaryCondition(dr);
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    if(rr<CutOffChargeChargeSquared[CurrentSystem])
    {
      switch(ChargeMethod)
      {
         case NONE:
           energy=0.0;
           break;
         case SHIFTED_COULOMB:
           energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge[CurrentSystem]);
           break;
         case TRUNCATED_COULOMB:
           energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
           break;
        case SMOOTHED_COULOMB:
           break;
        case EWALD:
        default:
          energy=COULOMBIC_CONVERSION_FACTOR*erfc(Alpha[CurrentSystem]*r)*
                                              chargeA*chargeB/r;

          if(!OmitEwaldFourier)
          {
            energy+=COULOMBIC_CONVERSION_FACTOR*erf(Alpha[CurrentSystem]*r)*
                             chargeA*chargeB/r;
          }
          break;
      }

      UExclusionConstraints[CurrentSystem]-=energy;
    }
  }
}



void CalculateConstraintsExclusionForce(void)
{
  int m;
  int typeA,typeB;
  REAL rr,r;
  REAL energy,force_factor;
  REAL chargeA,chargeB;
  VECTOR posA,posB,f,dr;


  for(m=0;m<NumberOfDistanceConstraints[CurrentSystem];m++)
  {
    posA=DistanceConstraints[CurrentSystem][m][0]->Position;
    posB=DistanceConstraints[CurrentSystem][m][1]->Position;

    typeA=DistanceConstraints[CurrentSystem][m][0]->Type;
    typeB=DistanceConstraints[CurrentSystem][m][1]->Type;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    dr=ApplyBoundaryCondition(dr);
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

    if(rr<CutOffVDWSquared)
    {
      PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);

      UExclusionConstraints[CurrentSystem]-=energy;

      f.x=force_factor*dr.x;
      f.y=force_factor*dr.y;
      f.z=force_factor*dr.z;

      DistanceConstraints[CurrentSystem][m][0]->Force.x+=f.x;
      DistanceConstraints[CurrentSystem][m][0]->Force.y+=f.y;
      DistanceConstraints[CurrentSystem][m][0]->Force.z+=f.z;

      DistanceConstraints[CurrentSystem][m][1]->Force.x-=f.x;
      DistanceConstraints[CurrentSystem][m][1]->Force.y-=f.y;
      DistanceConstraints[CurrentSystem][m][1]->Force.z-=f.z;

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


  for(m=0;m<NumberOfDistanceConstraints[CurrentSystem];m++)
  {
    posA=DistanceConstraints[CurrentSystem][m][0]->Position;
    posB=DistanceConstraints[CurrentSystem][m][1]->Position;

    typeA=DistanceConstraints[CurrentSystem][m][0]->Type;
    typeB=DistanceConstraints[CurrentSystem][m][1]->Type;

    chargeA=PseudoAtoms[typeA].Charge1;
    chargeB=PseudoAtoms[typeB].Charge1;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    dr=ApplyBoundaryCondition(dr);
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    if(rr<CutOffChargeChargeSquared[CurrentSystem])
    {
      switch(ChargeMethod)
      {
        case NONE:
          energy=0.0;
          force_factor=0.0;
          break;
        case SHIFTED_COULOMB:
          energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge[CurrentSystem]);
          force_factor=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
          break;
        case TRUNCATED_COULOMB:
          energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
          force_factor=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
          break;
        case SMOOTHED_COULOMB:
          break;
        case EWALD:
        default:
          energy=COULOMBIC_CONVERSION_FACTOR*erfc(Alpha[CurrentSystem]*r)*
                                              chargeA*chargeB/r;

          force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
              (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
              (r*rr);

          if(!OmitEwaldFourier)
          {
            energy+=COULOMBIC_CONVERSION_FACTOR*erf(Alpha[CurrentSystem]*r)*
                             chargeA*chargeB/r;

            force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                (erf(Alpha[CurrentSystem]*r)-2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
               (r*rr);
          }
          break;
      }

      UExclusionConstraints[CurrentSystem]-=energy;

      f.x=-force_factor*dr.x;
      f.y=-force_factor*dr.y;
      f.z=-force_factor*dr.z;

      DistanceConstraints[CurrentSystem][m][0]->Force.x+=f.x;
      DistanceConstraints[CurrentSystem][m][0]->Force.y+=f.y;
      DistanceConstraints[CurrentSystem][m][0]->Force.z+=f.z;

      DistanceConstraints[CurrentSystem][m][1]->Force.x-=f.x;
      DistanceConstraints[CurrentSystem][m][1]->Force.y-=f.y;
      DistanceConstraints[CurrentSystem][m][1]->Force.z-=f.z;

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


