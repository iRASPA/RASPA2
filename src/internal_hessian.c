/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'internal_hessian.h' is part of RASPA-2.0

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
#include "utils.h"
#include "internal_hessian.h"
#include "inter_hessian.h"
#include "internal_force.h"
#include "ewald.h"
#include "matrix.h"
#include "spectra.h"
#include "minimization.h"

// The first and second derivative of the potential with respect to separation distance are required:
// DF=D[U[r],r]/r
// DDF=D[D[U[r],r]/r,r]/r

static inline void HessianAtomicPositionStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,REAL f1,REAL f2,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;
      if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;
      if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;

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
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.y*dr.y*dr.x;
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.z*dr.z*dr.x;

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.z*dr.z*dr.y;

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.y*dr.y*dr.z;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;

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
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.x*dr.z*dr.x+f1*dr.z;
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.y*dr.y*dr.x;
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=f2*dr.y*dr.z*dr.x;
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=f2*dr.z*dr.z*dr.x;

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.x*dr.z*dr.y;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=f2*dr.y*dr.z*dr.y+f1*dr.z;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=f2*dr.z*dr.z*dr.y;

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.y*dr.z;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.x*dr.z*dr.z+f1*dr.x;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.y*dr.y*dr.z;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=f2*dr.y*dr.z*dr.z+f1*dr.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.x*dr.z*dr.x+f1*dr.z;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=f2*dr.y*dr.y*dr.x;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]-=f2*dr.y*dr.z*dr.x;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]-=f2*dr.z*dr.z*dr.x;

          if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.x*dr.z*dr.y;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]-=f2*dr.y*dr.z*dr.y+f1*dr.z;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]-=f2*dr.z*dr.z*dr.y;

          if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.y*dr.z;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.x*dr.z*dr.z+f1*dr.x;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.y*dr.y*dr.z;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]-=f2*dr.y*dr.z*dr.z+f1*dr.y;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
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
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.z*dr.x+f1*dr.z;
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.y*dr.x;
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.z*dr.y;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.z*dr.z+f1*dr.x;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.y*dr.z;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;

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
              break;
            case MONOCLINIC_GAMMA_ANGLE:
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

static inline void HessianAtomicStrainStrainLocal(REAL_MATRIX HessianMatrix,REAL f1,REAL f2,VECTOR dr)
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
          HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
          HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
          HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz
          HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;   // yyyy
          HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.z*dr.z;                    // yyzz
          HessianMatrix.element[n+2][n+2]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;   // zzzz
          break;
        case REGULAR:
          HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
          HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;    // xxxy
          HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;    // xxxz
          HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
          HessianMatrix.element[n][n+4]+=f2*dr.x*dr.x*dr.y*dr.z;                 // xxyz
          HessianMatrix.element[n][n+5]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

          HessianMatrix.element[n+1][n+1]+=0.5*f2*(dr.x*dr.y*dr.x*dr.y+dr.y*dr.x*dr.x*dr.y)+0.5*f1*(dr.x*dr.x+dr.y*dr.y);  // xyxy
          HessianMatrix.element[n+1][n+2]+=0.5*f2*(dr.x*dr.y*dr.x*dr.z+dr.y*dr.x*dr.x*dr.z)+0.5*f1*dr.y*dr.z;              // xyxz
          HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.x*dr.y*dr.y*dr.y+dr.y*dr.x*dr.y*dr.y)+f1*dr.x*dr.y;                  // xyyy
          HessianMatrix.element[n+1][n+4]+=0.5*f2*(dr.x*dr.y*dr.y*dr.z+dr.y*dr.x*dr.y*dr.z)+0.5*f1*dr.x*dr.z;              // xyyz
          HessianMatrix.element[n+1][n+5]+=0.5*f2*(dr.x*dr.y*dr.z*dr.z+dr.y*dr.x*dr.z*dr.z);                               // xyzz


          HessianMatrix.element[n+2][n+2]+=0.5*f2*(dr.x*dr.z*dr.x*dr.z+dr.z*dr.x*dr.x*dr.z)+0.5*f1*(dr.x*dr.x+dr.z*dr.z);  // xzxz
          HessianMatrix.element[n+2][n+3]+=0.5*f2*(dr.x*dr.z*dr.y*dr.y+dr.x*dr.z*dr.y*dr.y);                               // xzyy
          HessianMatrix.element[n+2][n+4]+=0.5*f2*(dr.x*dr.z*dr.y*dr.z+dr.x*dr.z*dr.y*dr.z)+0.5*f1*dr.x*dr.y;              // xzyz
          HessianMatrix.element[n+2][n+5]+=0.5*f2*(dr.x*dr.z*dr.z*dr.z+dr.x*dr.z*dr.z*dr.z)+f1*dr.x*dr.z;                  // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;   // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;       // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dr.y*dr.z*dr.z;                    // yyzz

          HessianMatrix.element[n+4][n+4]+=0.5*f2*(dr.y*dr.z*dr.y*dr.z+dr.y*dr.z*dr.y*dr.z)+0.5*f1*(dr.y*dr.y+dr.z*dr.z); // yzyz
          HessianMatrix.element[n+4][n+5]+=0.5*f2*(dr.y*dr.z*dr.z*dr.z+dr.y*dr.z*dr.z*dr.z)+f1*dr.y*dr.z;                 // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;   // zzzz
          break;
        case REGULAR_UPPER_TRIANGLE:
          HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
          HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;       // xxxy
          HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;       // xxxz
          HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
          HessianMatrix.element[n][n+4]+=f2*dr.x*dr.x*dr.y*dr.z;                    // xxyz
          HessianMatrix.element[n][n+5]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz


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
              HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.z;                 // xxyz
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;   // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;       // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+2][n+2]+=0.5*f2*(dr.y*dr.z*dr.y*dr.z+dr.y*dr.z*dr.y*dr.z)+0.5*f1*(dr.y*dr.y+dr.z*dr.z); // yzyz
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(dr.y*dr.z*dr.z*dr.z+dr.y*dr.z*dr.z*dr.z)+f1*dr.y*dr.z;                 // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;   // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;    // xxxz
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

              HessianMatrix.element[n+1][n+1]+=0.5*f2*(dr.x*dr.z*dr.x*dr.z+dr.z*dr.x*dr.x*dr.z)+0.5*f1*(dr.x*dr.x+dr.z*dr.z);  // xzxz
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(dr.x*dr.z*dr.y*dr.y+dr.x*dr.z*dr.y*dr.y);                               // xzyy
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.x*dr.z*dr.z*dr.z+dr.x*dr.z*dr.z*dr.z)+f1*dr.x*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;   // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;   // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;    // xxxy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

              HessianMatrix.element[n+1][n+1]+=0.5*f2*(dr.x*dr.y*dr.x*dr.y+dr.y*dr.x*dr.x*dr.y)+0.5*f1*(dr.x*dr.x+dr.y*dr.y);  // xyxy
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(dr.x*dr.y*dr.y*dr.y+dr.y*dr.x*dr.y*dr.y)+f1*dr.x*dr.y;                  // xyyy
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.x*dr.y*dr.z*dr.z+dr.y*dr.x*dr.z*dr.z);                               // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;   // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;   // zzzz
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.z;                    // xxyz
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;     // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.z*dr.y*dr.z+f1*dr.z*dr.z;     // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.z*dr.z*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;       // xxxz
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.z*dr.x*dr.z+f1*dr.z*dr.z;     // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.z*dr.y*dr.y;                  // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.z*dr.z*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;       // xxxy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz


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

static inline void GradientStrain(REAL *Gradient,REAL f1,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]+=f1*dr.x*dr.x+f1*dr.y*dr.y+f1*dr.z*dr.z;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n]+=f1*dr.x*dr.x;
          Gradient[n+1]+=f1*(dr.y)*(dr.y);
          Gradient[n+2]+=f1*(dr.z)*(dr.z);
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n]+=f1*dr.x*dr.x;
          Gradient[n+1]+=f1*dr.x*dr.y;
          Gradient[n+2]+=f1*dr.x*dr.z;
          Gradient[n+3]+=f1*dr.y*dr.y;
          Gradient[n+4]+=f1*dr.y*dr.z;
          Gradient[n+5]+=f1*dr.z*dr.z;
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
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]+=f1*dr.x*dr.x;
              Gradient[n+1]+=f1*dr.x*dr.z;
              Gradient[n+2]+=f1*dr.y*dr.y;
              Gradient[n+3]+=f1*dr.z*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]+=f1*dr.x*dr.x;
              Gradient[n+1]+=f1*dr.x*dr.y;
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




void HessianBendStrainPosition(INT_VECTOR3 index_i,INT_VECTOR3 index_j,INT_VECTOR3 index_k,REAL_MATRIX HessianMatrix,VECTOR vec_u,VECTOR vec_v,REAL u,REAL v,
           REAL rab,REAL rbc,VECTOR Rab,VECTOR Rbc,VECTOR dtA,VECTOR dtC,REAL DF,REAL DDF,REAL_MATRIX3x3 S,REAL CosTheta)
{
  int n;
  VECTOR veci,veck;
  VECTOR vec_ddtax_i,vec_ddtax_k;
  VECTOR vec_ddtay_i,vec_ddtay_k;
  VECTOR vec_ddtaz_i,vec_ddtaz_k;
  VECTOR vec_ddtcx_i,vec_ddtcx_k;
  VECTOR vec_ddtcy_i,vec_ddtcy_k;
  VECTOR vec_ddtcz_i,vec_ddtcz_k;

  n=NumberOfCoordinatesMinimizationVariables;

  veci.x=(1.0/CUBE(u))*(vec_u.x*vec_u.x*Rbc.x+vec_u.y*vec_u.x*Rbc.y+vec_u.z*vec_u.x*Rbc.z)-Rbc.x/u;
  veci.y=(1.0/CUBE(u))*(vec_u.x*vec_u.y*Rbc.x+vec_u.y*vec_u.y*Rbc.y+vec_u.z*vec_u.y*Rbc.z)-Rbc.y/u;
  veci.z=(1.0/CUBE(u))*(vec_u.x*vec_u.z*Rbc.x+vec_u.y*vec_u.z*Rbc.y+vec_u.z*vec_u.z*Rbc.z)-Rbc.z/u;

  veck.x=(1.0/CUBE(v))*(vec_v.x*vec_v.x*Rab.x+vec_v.y*vec_v.x*Rab.y+vec_v.z*vec_v.x*Rab.z)-Rab.x/v;
  veck.y=(1.0/CUBE(v))*(vec_v.x*vec_v.y*Rab.x+vec_v.y*vec_v.y*Rab.y+vec_v.z*vec_v.y*Rab.z)-Rab.y/v;
  veck.z=(1.0/CUBE(v))*(vec_v.x*vec_v.z*Rab.x+vec_v.y*vec_v.z*Rab.y+vec_v.z*vec_v.z*Rab.z)-Rab.z/v;


  // FINAL PART 1: derivative dtA.x=(Rbc.x-CosTheta*Rab.x)/rab
  vec_ddtax_i.x=(veci.x*Rab.x+CosTheta*((1.0/CUBE(u))*(vec_u.x*vec_u.x)-1/u))/rab+CosTheta*Rab.x*(1.0/CUBE(u))*vec_u.x-(1.0/CUBE(u))*vec_u.x*Rbc.x;
  vec_ddtax_i.y=(veci.y*Rab.x+CosTheta*(1.0/CUBE(u))*(vec_u.x*vec_u.y))/rab+CosTheta*Rab.x*(1.0/CUBE(u))*vec_u.y-(1.0/CUBE(u))*vec_u.y*Rbc.x;
  vec_ddtax_i.z=(veci.z*Rab.x+CosTheta*(1.0/CUBE(u))*(vec_u.x*vec_u.z))/rab+CosTheta*Rab.x*(1.0/CUBE(u))*vec_u.z-(1.0/CUBE(u))*vec_u.z*Rbc.x;

  vec_ddtax_k.x=veck.x*Rab.x/rab-((1.0/CUBE(v))*(vec_v.x*vec_v.x)-1/v)/rab;
  vec_ddtax_k.y=veck.y*Rab.x/rab-(1.0/CUBE(v))*(vec_v.x*vec_v.y)/rab;
  vec_ddtax_k.z=veck.z*Rab.x/rab-(1.0/CUBE(v))*(vec_v.x*vec_v.z)/rab;

  // FINAL PART 1: derivative dtA.y=(Rbc.y-CosTheta*Rab.y)/rab
  vec_ddtay_i.x=(veci.x*Rab.y+CosTheta*(1.0/CUBE(u))*(vec_u.y*vec_u.x))/rab+CosTheta*Rab.y*(1.0/CUBE(u))*vec_u.x-(1.0/CUBE(u))*vec_u.x*Rbc.y;
  vec_ddtay_i.y=(veci.y*Rab.y+CosTheta*((1.0/CUBE(u))*(vec_u.y*vec_u.y)-1.0/u))/rab+CosTheta*Rab.y*(1.0/CUBE(u))*vec_u.y-(1.0/CUBE(u))*vec_u.y*Rbc.y;
  vec_ddtay_i.z=(veci.z*Rab.y+CosTheta*(1.0/CUBE(u))*(vec_u.y*vec_u.z))/rab+CosTheta*Rab.y*(1.0/CUBE(u))*vec_u.z-(1.0/CUBE(u))*vec_u.z*Rbc.y;

  vec_ddtay_k.x=veck.x*Rab.y/rab-(1.0/CUBE(v))*(vec_v.y*vec_v.x)/rab;
  vec_ddtay_k.y=veck.y*Rab.y/rab-((1.0/CUBE(v))*(vec_v.y*vec_v.y)-1.0/v)/rab;
  vec_ddtay_k.z=veck.z*Rab.y/rab-(1.0/CUBE(v))*(vec_v.y*vec_v.z)/rab;

  // FINAL PART 1: derivative dtA.z=(Rbc.z-CosTheta*Rab.z)/rab
  vec_ddtaz_i.x=(veci.x*Rab.z+CosTheta*(1.0/CUBE(u))*(vec_u.z*vec_u.x))/rab+CosTheta*Rab.z*(1.0/CUBE(u))*vec_u.x-(1.0/CUBE(u))*vec_u.x*Rbc.z;
  vec_ddtaz_i.y=(veci.y*Rab.z+CosTheta*(1.0/CUBE(u))*(vec_u.z*vec_u.y))/rab+CosTheta*Rab.z*(1.0/CUBE(u))*vec_u.y-(1.0/CUBE(u))*vec_u.y*Rbc.z;
  vec_ddtaz_i.z=(veci.z*Rab.z+CosTheta*((1.0/CUBE(u))*(vec_u.z*vec_u.z)-1.0/u))/rab+CosTheta*Rab.z*(1.0/CUBE(u))*vec_u.z-(1.0/CUBE(u))*vec_u.z*Rbc.z;

  vec_ddtaz_k.x=veck.x*Rab.z/rab-(1.0/CUBE(v))*(vec_v.z*vec_v.x)/rab;
  vec_ddtaz_k.y=veck.y*Rab.z/rab-(1.0/CUBE(v))*(vec_v.z*vec_v.y)/rab;
  vec_ddtaz_k.z=veck.z*Rab.z/rab-((1.0/CUBE(v))*(vec_v.z*vec_v.z)-1.0/v)/rab;

  // FINAL PART 2: derivative dtC.x=(Rab.x-CosTheta*Rbc.x)/rbc
  vec_ddtcx_i.x=veci.x*Rbc.x/rbc-((1.0/CUBE(u))*(vec_u.x*vec_u.x)-1/u)/rbc;
  vec_ddtcx_i.y=veci.y*Rbc.x/rbc-(1.0/CUBE(u))*(vec_u.x*vec_u.y)/rbc;
  vec_ddtcx_i.z=veci.z*Rbc.x/rbc-(1.0/CUBE(u))*(vec_u.x*vec_u.z)/rbc;

  vec_ddtcx_k.x=(veck.x*Rbc.x+CosTheta*((1.0/CUBE(v))*(vec_v.x*vec_v.x)-1/v))/rbc+CosTheta*Rbc.x*(1.0/CUBE(v))*vec_v.x-(1.0/CUBE(v))*vec_v.x*Rab.x;
  vec_ddtcx_k.y=(veck.y*Rbc.x+CosTheta*(1.0/CUBE(v))*(vec_v.x*vec_v.y))/rbc+CosTheta*Rbc.x*(1.0/CUBE(v))*vec_v.y-(1.0/CUBE(v))*vec_v.y*Rab.x;
  vec_ddtcx_k.z=(veck.z*Rbc.x+CosTheta*(1.0/CUBE(v))*(vec_v.x*vec_v.z))/rbc+CosTheta*Rbc.x*(1.0/CUBE(v))*vec_v.z-(1.0/CUBE(v))*vec_v.z*Rab.x;

  // FINAL PART 2: derivative dtC.y=(Rab.y-CosTheta*Rbc.y)/rbc
  vec_ddtcy_i.x=veci.x*Rbc.y/rbc-(1.0/CUBE(u))*(vec_u.y*vec_u.x)/rbc;
  vec_ddtcy_i.y=veci.y*Rbc.y/rbc-((1.0/CUBE(u))*(vec_u.y*vec_u.y)-1.0/u)/rbc;
  vec_ddtcy_i.z=veci.z*Rbc.y/rbc-(1.0/CUBE(u))*(vec_u.y*vec_u.z)/rbc;

  vec_ddtcy_k.x=(veck.x*Rbc.y+CosTheta*(1.0/CUBE(v))*(vec_v.y*vec_v.x))/rbc+CosTheta*Rbc.y*(1.0/CUBE(v))*vec_v.x-(1.0/CUBE(v))*vec_v.x*Rab.y;
  vec_ddtcy_k.y=(veck.y*Rbc.y+CosTheta*((1.0/CUBE(v))*(vec_v.y*vec_v.y)-1.0/v))/rbc+CosTheta*Rbc.y*(1.0/CUBE(v))*vec_v.y-(1.0/CUBE(v))*vec_v.y*Rab.y;
  vec_ddtcy_k.z=(veck.z*Rbc.y+CosTheta*(1.0/CUBE(v))*(vec_v.y*vec_v.z))/rbc+CosTheta*Rbc.y*(1.0/CUBE(v))*vec_v.z-(1.0/CUBE(v))*vec_v.z*Rab.y;

  // FINAL PART 2: derivative dtC.z=(Rab.z-CosTheta*Rbc.z)/rbc
  vec_ddtcz_i.x=veci.x*Rbc.z/rbc-(1.0/CUBE(u))*(vec_u.z*vec_u.x)/rbc;
  vec_ddtcz_i.y=veci.y*Rbc.z/rbc-(1.0/CUBE(u))*(vec_u.z*vec_u.y)/rbc;
  vec_ddtcz_i.z=veci.z*Rbc.z/rbc-((1.0/CUBE(u))*(vec_u.z*vec_u.z)-1.0/u)/rbc;

  vec_ddtcz_k.x=(veck.x*Rbc.z+CosTheta*(1.0/CUBE(v))*(vec_v.z*vec_v.x))/rbc+CosTheta*Rbc.z*(1.0/CUBE(v))*vec_v.x-(1.0/CUBE(v))*vec_v.x*Rab.z;
  vec_ddtcz_k.y=(veck.y*Rbc.z+CosTheta*(1.0/CUBE(v))*(vec_v.z*vec_v.y))/rbc+CosTheta*Rbc.z*(1.0/CUBE(v))*vec_v.y-(1.0/CUBE(v))*vec_v.y*Rab.z;
  vec_ddtcz_k.z=(veck.z*Rbc.z+CosTheta*((1.0/CUBE(v))*(vec_v.z*vec_v.z)-1.0/v))/rbc+CosTheta*Rbc.z*(1.0/CUBE(v))*vec_v.z-(1.0/CUBE(v))*vec_v.z*Rab.z;


  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // contribution is zero (i.e. cancels out)
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          // S.ax=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [0]
          veci.x=vec_ddtax_i.x*DF*vec_u.x+dtA.x*(DDF*dtA.x*vec_u.x+DF) +vec_ddtcx_i.x*DF*vec_v.x+dtC.x*DDF*dtA.x*vec_v.x;
          veci.y=vec_ddtax_i.y*DF*vec_u.x+dtA.x*DDF*dtA.y*vec_u.x      +vec_ddtcx_i.y*DF*vec_v.x+dtC.x*DDF*dtA.y*vec_v.x;
          veci.z=vec_ddtax_i.z*DF*vec_u.x+dtA.x*DDF*dtA.z*vec_u.x      +vec_ddtcx_i.z*DF*vec_v.x+dtC.x*DDF*dtA.z*vec_v.x;

          veck.x=vec_ddtax_k.x*DF*vec_u.x+dtA.x*DDF*dtC.x*vec_u.x      +vec_ddtcx_k.x*DF*vec_v.x+dtC.x*(DDF*dtC.x*vec_v.x+DF);
          veck.y=vec_ddtax_k.y*DF*vec_u.x+dtA.x*DDF*dtC.y*vec_u.x      +vec_ddtcx_k.y*DF*vec_v.x+dtC.x*DDF*dtC.y*vec_v.x;
          veck.z=vec_ddtax_k.z*DF*vec_u.x+dtA.x*DDF*dtC.z*vec_u.x      +vec_ddtcx_k.z*DF*vec_v.x+dtC.x*DDF*dtC.z*vec_v.x;

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=veci.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(veci.x+veck.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(veci.y+veck.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(veci.z+veck.z);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][n]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n]+=veck.z;

          veci.x=vec_ddtay_i.x*DF*vec_u.y+dtA.y*DDF*dtA.x*vec_u.y      +vec_ddtcy_i.x*DF*vec_v.y+dtC.y*DDF*dtA.x*vec_v.y;
          veci.y=vec_ddtay_i.y*DF*vec_u.y+dtA.y*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcy_i.y*DF*vec_v.y+dtC.y*DDF*dtA.y*vec_v.y;
          veci.z=vec_ddtay_i.z*DF*vec_u.y+dtA.y*DDF*dtA.z*vec_u.y      +vec_ddtcy_i.z*DF*vec_v.y+dtC.y*DDF*dtA.z*vec_v.y;

          veck.x=vec_ddtay_k.x*DF*vec_u.y+dtA.y*DDF*dtC.x*vec_u.y      +vec_ddtcy_k.x*DF*vec_v.y+dtC.y*DDF*dtC.x*vec_v.y;
          veck.y=vec_ddtay_k.y*DF*vec_u.y+dtA.y*DDF*dtC.y*vec_u.y      +vec_ddtcy_k.y*DF*vec_v.y+dtC.y*(DDF*dtC.y*vec_v.y+DF);
          veck.z=vec_ddtay_k.z*DF*vec_u.y+dtA.y*DDF*dtC.z*vec_u.y      +vec_ddtcy_k.z*DF*vec_v.y+dtC.y*DDF*dtC.z*vec_v.y;

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=veci.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=(veci.x+veck.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=(veci.y+veck.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=(veci.z+veck.z);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+1]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+1]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+1]+=veck.z;

          veci.x=vec_ddtaz_i.x*DF*vec_u.z+dtA.z*DDF*dtA.x*vec_u.z      +vec_ddtcz_i.x*DF*vec_v.z+dtC.z*DDF*dtA.x*vec_v.z;
          veci.y=vec_ddtaz_i.y*DF*vec_u.z+dtA.z*DDF*dtA.y*vec_u.z      +vec_ddtcz_i.y*DF*vec_v.z+dtC.z*DDF*dtA.y*vec_v.z;
          veci.z=vec_ddtaz_i.z*DF*vec_u.z+dtA.z*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcz_i.z*DF*vec_v.z+dtC.z*DDF*dtA.z*vec_v.z;

          veck.x=vec_ddtaz_k.x*DF*vec_u.z+dtA.z*DDF*dtC.x*vec_u.z      +vec_ddtcz_k.x*DF*vec_v.z+dtC.z*DDF*dtC.x*vec_v.z;
          veck.y=vec_ddtaz_k.y*DF*vec_u.z+dtA.z*DDF*dtC.y*vec_u.z      +vec_ddtcz_k.y*DF*vec_v.z+dtC.z*DDF*dtC.y*vec_v.z;
          veck.z=vec_ddtaz_k.z*DF*vec_u.z+dtA.z*DDF*dtC.z*vec_u.z      +vec_ddtcz_k.z*DF*vec_v.z+dtC.z*(DDF*dtC.z*vec_v.z+DF);

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=veci.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=(veci.x+veck.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=(veci.y+veck.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=(veci.z+veck.z);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+2]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+2]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+2]+=veck.z;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          // S.ax=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [0]
          veci.x=vec_ddtax_i.x*DF*vec_u.x+dtA.x*(DDF*dtA.x*vec_u.x+DF) +vec_ddtcx_i.x*DF*vec_v.x+dtC.x*DDF*dtA.x*vec_v.x;
          veci.y=vec_ddtax_i.y*DF*vec_u.x+dtA.x*DDF*dtA.y*vec_u.x      +vec_ddtcx_i.y*DF*vec_v.x+dtC.x*DDF*dtA.y*vec_v.x;
          veci.z=vec_ddtax_i.z*DF*vec_u.x+dtA.x*DDF*dtA.z*vec_u.x      +vec_ddtcx_i.z*DF*vec_v.x+dtC.x*DDF*dtA.z*vec_v.x;

          veck.x=vec_ddtax_k.x*DF*vec_u.x+dtA.x*DDF*dtC.x*vec_u.x      +vec_ddtcx_k.x*DF*vec_v.x+dtC.x*(DDF*dtC.x*vec_v.x+DF);
          veck.y=vec_ddtax_k.y*DF*vec_u.x+dtA.x*DDF*dtC.y*vec_u.x      +vec_ddtcx_k.y*DF*vec_v.x+dtC.x*DDF*dtC.y*vec_v.x;
          veck.z=vec_ddtax_k.z*DF*vec_u.x+dtA.x*DDF*dtC.z*vec_u.x      +vec_ddtcx_k.z*DF*vec_v.x+dtC.x*DDF*dtC.z*vec_v.x;

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=veci.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(veci.x+veck.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(veci.y+veck.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(veci.z+veck.z);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][n]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n]+=veck.z;

          // S.ay=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [1] and [3]
          veci.x=vec_ddtax_i.x*DF*vec_u.y+dtA.x*DDF*dtA.x*vec_u.y      +vec_ddtcx_i.x*DF*vec_v.y+dtC.x*DDF*dtA.x*vec_v.y;
          veci.y=vec_ddtax_i.y*DF*vec_u.y+dtA.x*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcx_i.y*DF*vec_v.y+dtC.x*DDF*dtA.y*vec_v.y;
          veci.z=vec_ddtax_i.z*DF*vec_u.y+dtA.x*DDF*dtA.z*vec_u.y      +vec_ddtcx_i.z*DF*vec_v.y+dtC.x*DDF*dtA.z*vec_v.y;

          veck.x=vec_ddtax_k.x*DF*vec_u.y+dtA.x*DDF*dtC.x*vec_u.y      +vec_ddtcx_k.x*DF*vec_v.y+dtC.x*DDF*dtC.x*vec_v.y;
          veck.y=vec_ddtax_k.y*DF*vec_u.y+dtA.x*DDF*dtC.y*vec_u.y      +vec_ddtcx_k.y*DF*vec_v.y+dtC.x*(DDF*dtC.y*vec_v.y+DF);
          veck.z=vec_ddtax_k.z*DF*vec_u.y+dtA.x*DDF*dtC.z*vec_u.y      +vec_ddtcx_k.z*DF*vec_v.y+dtC.x*DDF*dtC.z*vec_v.y;

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=veci.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=(veci.x+veck.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=(veci.y+veck.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=(veci.z+veck.z);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+1]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+1]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+1]+=veck.z;

          // S.cx=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;   [2] and [6]
          veci.x=vec_ddtax_i.x*DF*vec_u.z+dtA.x*DDF*dtA.x*vec_u.z      +vec_ddtcx_i.x*DF*vec_v.z+dtC.x*DDF*dtA.x*vec_v.z;
          veci.y=vec_ddtax_i.y*DF*vec_u.z+dtA.x*DDF*dtA.y*vec_u.z      +vec_ddtcx_i.y*DF*vec_v.z+dtC.x*DDF*dtA.y*vec_v.z;
          veci.z=vec_ddtax_i.z*DF*vec_u.z+dtA.x*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcx_i.z*DF*vec_v.z+dtC.x*DDF*dtA.z*vec_v.z;

          veck.x=vec_ddtax_k.x*DF*vec_u.z+dtA.x*DDF*dtC.x*vec_u.z      +vec_ddtcx_k.x*DF*vec_v.z+dtC.x*DDF*dtC.x*vec_v.z;
          veck.y=vec_ddtax_k.y*DF*vec_u.z+dtA.x*DDF*dtC.y*vec_u.z      +vec_ddtcx_k.y*DF*vec_v.z+dtC.x*DDF*dtC.y*vec_v.z;
          veck.z=vec_ddtax_k.z*DF*vec_u.z+dtA.x*DDF*dtC.z*vec_u.z      +vec_ddtcx_k.z*DF*vec_v.z+dtC.x*(DDF*dtC.z*vec_v.z+DF);

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=veci.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=(veci.x+veck.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=(veci.y+veck.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=(veci.z+veck.z);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+2]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+2]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+2]+=veck.z;

          // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [4]
          veci.x=vec_ddtay_i.x*DF*vec_u.y+dtA.y*DDF*dtA.x*vec_u.y      +vec_ddtcy_i.x*DF*vec_v.y+dtC.y*DDF*dtA.x*vec_v.y;
          veci.y=vec_ddtay_i.y*DF*vec_u.y+dtA.y*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcy_i.y*DF*vec_v.y+dtC.y*DDF*dtA.y*vec_v.y;
          veci.z=vec_ddtay_i.z*DF*vec_u.y+dtA.y*DDF*dtA.z*vec_u.y      +vec_ddtcy_i.z*DF*vec_v.y+dtC.y*DDF*dtA.z*vec_v.y;

          veck.x=vec_ddtay_k.x*DF*vec_u.y+dtA.y*DDF*dtC.x*vec_u.y      +vec_ddtcy_k.x*DF*vec_v.y+dtC.y*DDF*dtC.x*vec_v.y;
          veck.y=vec_ddtay_k.y*DF*vec_u.y+dtA.y*DDF*dtC.y*vec_u.y      +vec_ddtcy_k.y*DF*vec_v.y+dtC.y*(DDF*dtC.y*vec_v.y+DF);
          veck.z=vec_ddtay_k.z*DF*vec_u.y+dtA.y*DDF*dtC.z*vec_u.y      +vec_ddtcy_k.z*DF*vec_v.y+dtC.y*DDF*dtC.z*vec_v.y;

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=veci.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=(veci.x+veck.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=(veci.y+veck.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=(veci.z+veck.z);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+3]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+3]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+3]+=veck.z;

          // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
          veci.x=vec_ddtay_i.x*DF*vec_u.z+dtA.y*DDF*dtA.x*vec_u.z      +vec_ddtcy_i.x*DF*vec_v.z+dtC.y*DDF*dtA.x*vec_v.z;
          veci.y=vec_ddtay_i.y*DF*vec_u.z+dtA.y*DDF*dtA.y*vec_u.z      +vec_ddtcy_i.y*DF*vec_v.z+dtC.y*DDF*dtA.y*vec_v.z;
          veci.z=vec_ddtay_i.z*DF*vec_u.z+dtA.y*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcy_i.z*DF*vec_v.z+dtC.y*DDF*dtA.z*vec_v.z;

          veck.x=vec_ddtay_k.x*DF*vec_u.z+dtA.y*DDF*dtC.x*vec_u.z      +vec_ddtcy_k.x*DF*vec_v.z+dtC.y*DDF*dtC.x*vec_v.z;
          veck.y=vec_ddtay_k.y*DF*vec_u.z+dtA.y*DDF*dtC.y*vec_u.z      +vec_ddtcy_k.y*DF*vec_v.z+dtC.y*DDF*dtC.y*vec_v.z;
          veck.z=vec_ddtay_k.z*DF*vec_u.z+dtA.y*DDF*dtC.z*vec_u.z      +vec_ddtcy_k.z*DF*vec_v.z+dtC.y*(DDF*dtC.z*vec_v.z+DF);

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=veci.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]-=(veci.x+veck.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]-=(veci.y+veck.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]-=(veci.z+veck.z);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+4]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+4]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+4]+=veck.z;

          // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
          veci.x=vec_ddtaz_i.x*DF*vec_u.z+dtA.z*DDF*dtA.x*vec_u.z      +vec_ddtcz_i.x*DF*vec_v.z+dtC.z*DDF*dtA.x*vec_v.z;
          veci.y=vec_ddtaz_i.y*DF*vec_u.z+dtA.z*DDF*dtA.y*vec_u.z      +vec_ddtcz_i.y*DF*vec_v.z+dtC.z*DDF*dtA.y*vec_v.z;
          veci.z=vec_ddtaz_i.z*DF*vec_u.z+dtA.z*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcz_i.z*DF*vec_v.z+dtC.z*DDF*dtA.z*vec_v.z;

          veck.x=vec_ddtaz_k.x*DF*vec_u.z+dtA.z*DDF*dtC.x*vec_u.z      +vec_ddtcz_k.x*DF*vec_v.z+dtC.z*DDF*dtC.x*vec_v.z;
          veck.y=vec_ddtaz_k.y*DF*vec_u.z+dtA.z*DDF*dtC.y*vec_u.z      +vec_ddtcz_k.y*DF*vec_v.z+dtC.z*DDF*dtC.y*vec_v.z;
          veck.z=vec_ddtaz_k.z*DF*vec_u.z+dtA.z*DDF*dtC.z*vec_u.z      +vec_ddtcz_k.z*DF*vec_v.z+dtC.z*(DDF*dtC.z*vec_v.z+DF);

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=veci.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]-=(veci.x+veck.x);
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]-=(veci.y+veck.y);
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]-=(veci.z+veck.z);

          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+5]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+5]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+5]+=veck.z;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // S.ax=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [0]
              veci.x=vec_ddtax_i.x*DF*vec_u.x+dtA.x*(DDF*dtA.x*vec_u.x+DF) +vec_ddtcx_i.x*DF*vec_v.x+dtC.x*DDF*dtA.x*vec_v.x;
              veci.y=vec_ddtax_i.y*DF*vec_u.x+dtA.x*DDF*dtA.y*vec_u.x      +vec_ddtcx_i.y*DF*vec_v.x+dtC.x*DDF*dtA.y*vec_v.x;
              veci.z=vec_ddtax_i.z*DF*vec_u.x+dtA.x*DDF*dtA.z*vec_u.x      +vec_ddtcx_i.z*DF*vec_v.x+dtC.x*DDF*dtA.z*vec_v.x;

              veck.x=vec_ddtax_k.x*DF*vec_u.x+dtA.x*DDF*dtC.x*vec_u.x      +vec_ddtcx_k.x*DF*vec_v.x+dtC.x*(DDF*dtC.x*vec_v.x+DF);
              veck.y=vec_ddtax_k.y*DF*vec_u.x+dtA.x*DDF*dtC.y*vec_u.x      +vec_ddtcx_k.y*DF*vec_v.x+dtC.x*DDF*dtC.y*vec_v.x;
              veck.z=vec_ddtax_k.z*DF*vec_u.x+dtA.x*DDF*dtC.z*vec_u.x      +vec_ddtcx_k.z*DF*vec_v.x+dtC.x*DDF*dtC.z*vec_v.x;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [4]
              veci.x=vec_ddtay_i.x*DF*vec_u.y+dtA.y*DDF*dtA.x*vec_u.y      +vec_ddtcy_i.x*DF*vec_v.y+dtC.y*DDF*dtA.x*vec_v.y;
              veci.y=vec_ddtay_i.y*DF*vec_u.y+dtA.y*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcy_i.y*DF*vec_v.y+dtC.y*DDF*dtA.y*vec_v.y;
              veci.z=vec_ddtay_i.z*DF*vec_u.y+dtA.y*DDF*dtA.z*vec_u.y      +vec_ddtcy_i.z*DF*vec_v.y+dtC.y*DDF*dtA.z*vec_v.y;

              veck.x=vec_ddtay_k.x*DF*vec_u.y+dtA.y*DDF*dtC.x*vec_u.y      +vec_ddtcy_k.x*DF*vec_v.y+dtC.y*DDF*dtC.x*vec_v.y;
              veck.y=vec_ddtay_k.y*DF*vec_u.y+dtA.y*DDF*dtC.y*vec_u.y      +vec_ddtcy_k.y*DF*vec_v.y+dtC.y*(DDF*dtC.y*vec_v.y+DF);
              veck.z=vec_ddtay_k.z*DF*vec_u.y+dtA.y*DDF*dtC.z*vec_u.y      +vec_ddtcy_k.z*DF*vec_v.y+dtC.y*DDF*dtC.z*vec_v.y;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+1]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+1]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+1]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
              veci.x=vec_ddtay_i.x*DF*vec_u.z+dtA.y*DDF*dtA.x*vec_u.z      +vec_ddtcy_i.x*DF*vec_v.z+dtC.y*DDF*dtA.x*vec_v.z;
              veci.y=vec_ddtay_i.y*DF*vec_u.z+dtA.y*DDF*dtA.y*vec_u.z      +vec_ddtcy_i.y*DF*vec_v.z+dtC.y*DDF*dtA.y*vec_v.z;
              veci.z=vec_ddtay_i.z*DF*vec_u.z+dtA.y*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcy_i.z*DF*vec_v.z+dtC.y*DDF*dtA.z*vec_v.z;

              veck.x=vec_ddtay_k.x*DF*vec_u.z+dtA.y*DDF*dtC.x*vec_u.z      +vec_ddtcy_k.x*DF*vec_v.z+dtC.y*DDF*dtC.x*vec_v.z;
              veck.y=vec_ddtay_k.y*DF*vec_u.z+dtA.y*DDF*dtC.y*vec_u.z      +vec_ddtcy_k.y*DF*vec_v.z+dtC.y*DDF*dtC.y*vec_v.z;
              veck.z=vec_ddtay_k.z*DF*vec_u.z+dtA.y*DDF*dtC.z*vec_u.z      +vec_ddtcy_k.z*DF*vec_v.z+dtC.y*(DDF*dtC.z*vec_v.z+DF);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+2]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+2]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+2]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
              veci.x=vec_ddtaz_i.x*DF*vec_u.z+dtA.z*DDF*dtA.x*vec_u.z      +vec_ddtcz_i.x*DF*vec_v.z+dtC.z*DDF*dtA.x*vec_v.z;
              veci.y=vec_ddtaz_i.y*DF*vec_u.z+dtA.z*DDF*dtA.y*vec_u.z      +vec_ddtcz_i.y*DF*vec_v.z+dtC.z*DDF*dtA.y*vec_v.z;
              veci.z=vec_ddtaz_i.z*DF*vec_u.z+dtA.z*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcz_i.z*DF*vec_v.z+dtC.z*DDF*dtA.z*vec_v.z;

              veck.x=vec_ddtaz_k.x*DF*vec_u.z+dtA.z*DDF*dtC.x*vec_u.z      +vec_ddtcz_k.x*DF*vec_v.z+dtC.z*DDF*dtC.x*vec_v.z;
              veck.y=vec_ddtaz_k.y*DF*vec_u.z+dtA.z*DDF*dtC.y*vec_u.z      +vec_ddtcz_k.y*DF*vec_v.z+dtC.z*DDF*dtC.y*vec_v.z;
              veck.z=vec_ddtaz_k.z*DF*vec_u.z+dtA.z*DDF*dtC.z*vec_u.z      +vec_ddtcz_k.z*DF*vec_v.z+dtC.z*(DDF*dtC.z*vec_v.z+DF);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+3]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+3]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+3]+=veck.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              // S.ax=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [0]
              veci.x=vec_ddtax_i.x*DF*vec_u.x+dtA.x*(DDF*dtA.x*vec_u.x+DF) +vec_ddtcx_i.x*DF*vec_v.x+dtC.x*DDF*dtA.x*vec_v.x;
              veci.y=vec_ddtax_i.y*DF*vec_u.x+dtA.x*DDF*dtA.y*vec_u.x      +vec_ddtcx_i.y*DF*vec_v.x+dtC.x*DDF*dtA.y*vec_v.x;
              veci.z=vec_ddtax_i.z*DF*vec_u.x+dtA.x*DDF*dtA.z*vec_u.x      +vec_ddtcx_i.z*DF*vec_v.x+dtC.x*DDF*dtA.z*vec_v.x;

              veck.x=vec_ddtax_k.x*DF*vec_u.x+dtA.x*DDF*dtC.x*vec_u.x      +vec_ddtcx_k.x*DF*vec_v.x+dtC.x*(DDF*dtC.x*vec_v.x+DF);
              veck.y=vec_ddtax_k.y*DF*vec_u.x+dtA.x*DDF*dtC.y*vec_u.x      +vec_ddtcx_k.y*DF*vec_v.x+dtC.x*DDF*dtC.y*vec_v.x;
              veck.z=vec_ddtax_k.z*DF*vec_u.x+dtA.x*DDF*dtC.z*vec_u.x      +vec_ddtcx_k.z*DF*vec_v.x+dtC.x*DDF*dtC.z*vec_v.x;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n]+=veck.z;

              // S.cx=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;   [2] and [6]
              veci.x=vec_ddtax_i.x*DF*vec_u.z+dtA.x*DDF*dtA.x*vec_u.z      +vec_ddtcx_i.x*DF*vec_v.z+dtC.x*DDF*dtA.x*vec_v.z;
              veci.y=vec_ddtax_i.y*DF*vec_u.z+dtA.x*DDF*dtA.y*vec_u.z      +vec_ddtcx_i.y*DF*vec_v.z+dtC.x*DDF*dtA.y*vec_v.z;
              veci.z=vec_ddtax_i.z*DF*vec_u.z+dtA.x*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcx_i.z*DF*vec_v.z+dtC.x*DDF*dtA.z*vec_v.z;

              veck.x=vec_ddtax_k.x*DF*vec_u.z+dtA.x*DDF*dtC.x*vec_u.z      +vec_ddtcx_k.x*DF*vec_v.z+dtC.x*DDF*dtC.x*vec_v.z;
              veck.y=vec_ddtax_k.y*DF*vec_u.z+dtA.x*DDF*dtC.y*vec_u.z      +vec_ddtcx_k.y*DF*vec_v.z+dtC.x*DDF*dtC.y*vec_v.z;
              veck.z=vec_ddtax_k.z*DF*vec_u.z+dtA.x*DDF*dtC.z*vec_u.z      +vec_ddtcx_k.z*DF*vec_v.z+dtC.x*(DDF*dtC.z*vec_v.z+DF);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+1]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+1]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+1]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [4]
              veci.x=vec_ddtay_i.x*DF*vec_u.y+dtA.y*DDF*dtA.x*vec_u.y      +vec_ddtcy_i.x*DF*vec_v.y+dtC.y*DDF*dtA.x*vec_v.y;
              veci.y=vec_ddtay_i.y*DF*vec_u.y+dtA.y*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcy_i.y*DF*vec_v.y+dtC.y*DDF*dtA.y*vec_v.y;
              veci.z=vec_ddtay_i.z*DF*vec_u.y+dtA.y*DDF*dtA.z*vec_u.y      +vec_ddtcy_i.z*DF*vec_v.y+dtC.y*DDF*dtA.z*vec_v.y;

              veck.x=vec_ddtay_k.x*DF*vec_u.y+dtA.y*DDF*dtC.x*vec_u.y      +vec_ddtcy_k.x*DF*vec_v.y+dtC.y*DDF*dtC.x*vec_v.y;
              veck.y=vec_ddtay_k.y*DF*vec_u.y+dtA.y*DDF*dtC.y*vec_u.y      +vec_ddtcy_k.y*DF*vec_v.y+dtC.y*(DDF*dtC.y*vec_v.y+DF);
              veck.z=vec_ddtay_k.z*DF*vec_u.y+dtA.y*DDF*dtC.z*vec_u.y      +vec_ddtcy_k.z*DF*vec_v.y+dtC.y*DDF*dtC.z*vec_v.y;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+2]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+2]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+2]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
              veci.x=vec_ddtaz_i.x*DF*vec_u.z+dtA.z*DDF*dtA.x*vec_u.z      +vec_ddtcz_i.x*DF*vec_v.z+dtC.z*DDF*dtA.x*vec_v.z;
              veci.y=vec_ddtaz_i.y*DF*vec_u.z+dtA.z*DDF*dtA.y*vec_u.z      +vec_ddtcz_i.y*DF*vec_v.z+dtC.z*DDF*dtA.y*vec_v.z;
              veci.z=vec_ddtaz_i.z*DF*vec_u.z+dtA.z*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcz_i.z*DF*vec_v.z+dtC.z*DDF*dtA.z*vec_v.z;

              veck.x=vec_ddtaz_k.x*DF*vec_u.z+dtA.z*DDF*dtC.x*vec_u.z      +vec_ddtcz_k.x*DF*vec_v.z+dtC.z*DDF*dtC.x*vec_v.z;
              veck.y=vec_ddtaz_k.y*DF*vec_u.z+dtA.z*DDF*dtC.y*vec_u.z      +vec_ddtcz_k.y*DF*vec_v.z+dtC.z*DDF*dtC.y*vec_v.z;
              veck.z=vec_ddtaz_k.z*DF*vec_u.z+dtA.z*DDF*dtC.z*vec_u.z      +vec_ddtcz_k.z*DF*vec_v.z+dtC.z*(DDF*dtC.z*vec_v.z+DF);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+3]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+3]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+3]+=veck.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // S.ax=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [0]
              veci.x=vec_ddtax_i.x*DF*vec_u.x+dtA.x*(DDF*dtA.x*vec_u.x+DF) +vec_ddtcx_i.x*DF*vec_v.x+dtC.x*DDF*dtA.x*vec_v.x;
              veci.y=vec_ddtax_i.y*DF*vec_u.x+dtA.x*DDF*dtA.y*vec_u.x      +vec_ddtcx_i.y*DF*vec_v.x+dtC.x*DDF*dtA.y*vec_v.x;
              veci.z=vec_ddtax_i.z*DF*vec_u.x+dtA.x*DDF*dtA.z*vec_u.x      +vec_ddtcx_i.z*DF*vec_v.x+dtC.x*DDF*dtA.z*vec_v.x;

              veck.x=vec_ddtax_k.x*DF*vec_u.x+dtA.x*DDF*dtC.x*vec_u.x      +vec_ddtcx_k.x*DF*vec_v.x+dtC.x*(DDF*dtC.x*vec_v.x+DF);
              veck.y=vec_ddtax_k.y*DF*vec_u.x+dtA.x*DDF*dtC.y*vec_u.x      +vec_ddtcx_k.y*DF*vec_v.x+dtC.x*DDF*dtC.y*vec_v.x;
              veck.z=vec_ddtax_k.z*DF*vec_u.x+dtA.x*DDF*dtC.z*vec_u.x      +vec_ddtcx_k.z*DF*vec_v.x+dtC.x*DDF*dtC.z*vec_v.x;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n]+=veck.z;

              // S.ay=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [1] and [3]
              veci.x=vec_ddtax_i.x*DF*vec_u.y+dtA.x*DDF*dtA.x*vec_u.y      +vec_ddtcx_i.x*DF*vec_v.y+dtC.x*DDF*dtA.x*vec_v.y;
              veci.y=vec_ddtax_i.y*DF*vec_u.y+dtA.x*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcx_i.y*DF*vec_v.y+dtC.x*DDF*dtA.y*vec_v.y;
              veci.z=vec_ddtax_i.z*DF*vec_u.y+dtA.x*DDF*dtA.z*vec_u.y      +vec_ddtcx_i.z*DF*vec_v.y+dtC.x*DDF*dtA.z*vec_v.y;

              veck.x=vec_ddtax_k.x*DF*vec_u.y+dtA.x*DDF*dtC.x*vec_u.y      +vec_ddtcx_k.x*DF*vec_v.y+dtC.x*DDF*dtC.x*vec_v.y;
              veck.y=vec_ddtax_k.y*DF*vec_u.y+dtA.x*DDF*dtC.y*vec_u.y      +vec_ddtcx_k.y*DF*vec_v.y+dtC.x*(DDF*dtC.y*vec_v.y+DF);
              veck.z=vec_ddtax_k.z*DF*vec_u.y+dtA.x*DDF*dtC.z*vec_u.y      +vec_ddtcx_k.z*DF*vec_v.y+dtC.x*DDF*dtC.z*vec_v.y;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+1]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+1]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+1]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [4]
              veci.x=vec_ddtay_i.x*DF*vec_u.y+dtA.y*DDF*dtA.x*vec_u.y      +vec_ddtcy_i.x*DF*vec_v.y+dtC.y*DDF*dtA.x*vec_v.y;
              veci.y=vec_ddtay_i.y*DF*vec_u.y+dtA.y*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcy_i.y*DF*vec_v.y+dtC.y*DDF*dtA.y*vec_v.y;
              veci.z=vec_ddtay_i.z*DF*vec_u.y+dtA.y*DDF*dtA.z*vec_u.y      +vec_ddtcy_i.z*DF*vec_v.y+dtC.y*DDF*dtA.z*vec_v.y;

              veck.x=vec_ddtay_k.x*DF*vec_u.y+dtA.y*DDF*dtC.x*vec_u.y      +vec_ddtcy_k.x*DF*vec_v.y+dtC.y*DDF*dtC.x*vec_v.y;
              veck.y=vec_ddtay_k.y*DF*vec_u.y+dtA.y*DDF*dtC.y*vec_u.y      +vec_ddtcy_k.y*DF*vec_v.y+dtC.y*(DDF*dtC.y*vec_v.y+DF);
              veck.z=vec_ddtay_k.z*DF*vec_u.y+dtA.y*DDF*dtC.z*vec_u.y      +vec_ddtcy_k.z*DF*vec_v.y+dtC.y*DDF*dtC.z*vec_v.y;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+2]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+2]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+2]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
              veci.x=vec_ddtaz_i.x*DF*vec_u.z+dtA.z*DDF*dtA.x*vec_u.z      +vec_ddtcz_i.x*DF*vec_v.z+dtC.z*DDF*dtA.x*vec_v.z;
              veci.y=vec_ddtaz_i.y*DF*vec_u.z+dtA.z*DDF*dtA.y*vec_u.z      +vec_ddtcz_i.y*DF*vec_v.z+dtC.z*DDF*dtA.y*vec_v.z;
              veci.z=vec_ddtaz_i.z*DF*vec_u.z+dtA.z*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcz_i.z*DF*vec_v.z+dtC.z*DDF*dtA.z*vec_v.z;

              veck.x=vec_ddtaz_k.x*DF*vec_u.z+dtA.z*DDF*dtC.x*vec_u.z      +vec_ddtcz_k.x*DF*vec_v.z+dtC.z*DDF*dtC.x*vec_v.z;
              veck.y=vec_ddtaz_k.y*DF*vec_u.z+dtA.z*DDF*dtC.y*vec_u.z      +vec_ddtcz_k.y*DF*vec_v.z+dtC.z*DDF*dtC.y*vec_v.z;
              veck.z=vec_ddtaz_k.z*DF*vec_u.z+dtA.z*DDF*dtC.z*vec_u.z      +vec_ddtcz_k.z*DF*vec_v.z+dtC.z*(DDF*dtC.z*vec_v.z+DF);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=veci.z;

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=(veci.x+veck.x);
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=(veci.y+veck.y);
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=(veci.z+veck.z);

              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+3]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+3]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+3]+=veck.z;
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

// strain-strain

void HessianBendStrainStrain(REAL_MATRIX HessianMatrix,VECTOR vec_u,VECTOR vec_v,REAL u,REAL v,REAL DF,REAL DDF,REAL_MATRIX3x3 S,REAL CosTheta)
{
  int n;
  REAL_MATRIX3x3 T;

  n=NumberOfCoordinatesMinimizationVariables;

  T.ax=(vec_u.x*vec_v.x+vec_v.x*vec_u.x)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.x*vec_u.x)/SQR(u)+(vec_v.x*vec_v.x)/SQR(v));
  T.ay=(vec_u.x*vec_v.y+vec_v.x*vec_u.y)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.x*vec_u.y)/SQR(u)+(vec_v.x*vec_v.y)/SQR(v));
  T.az=(vec_u.x*vec_v.z+vec_v.x*vec_u.z)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.x*vec_u.z)/SQR(u)+(vec_v.x*vec_v.z)/SQR(v));

  T.bx=(vec_u.y*vec_v.x+vec_v.y*vec_u.x)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.y*vec_u.x)/SQR(u)+(vec_v.y*vec_v.x)/SQR(v));
  T.by=(vec_u.y*vec_v.y+vec_v.y*vec_u.y)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.y*vec_u.y)/SQR(u)+(vec_v.y*vec_v.y)/SQR(v));
  T.bz=(vec_u.y*vec_v.z+vec_v.y*vec_u.z)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.y*vec_u.z)/SQR(u)+(vec_v.y*vec_v.z)/SQR(v));

  T.cx=(vec_u.z*vec_v.x+vec_v.z*vec_u.x)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.z*vec_u.x)/SQR(u)+(vec_v.z*vec_v.x)/SQR(v));
  T.cy=(vec_u.z*vec_v.y+vec_v.z*vec_u.y)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.z*vec_u.y)/SQR(u)+(vec_v.z*vec_v.y)/SQR(v));
  T.cz=(vec_u.z*vec_v.z+vec_v.z*vec_u.z)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.z*vec_u.z)/SQR(u)+(vec_v.z*vec_v.z)/SQR(v));

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                    +T.ax*T.ax)+2.0*S.ax;
      HessianMatrix.element[n][n]+=2.0*(DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                    +T.ax*T.by));
      HessianMatrix.element[n][n]+=2.0*(DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                    +T.ax*T.cz));
      HessianMatrix.element[n][n]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                    -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                    +T.by*T.by)+2.0*S.by;
      HessianMatrix.element[n][n]+=2.0*(DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                    -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                    +T.by*T.cz));
      HessianMatrix.element[n][n]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                    -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                    +T.cz*T.cz)+2.0*S.cz;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                        +T.ax*T.ax)+2.0*S.ax;
          HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.ax*T.by);
          HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.cz);

          HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.by*T.by)+2.0*S.by;
          HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.by*T.cz);

          HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.cz*T.cz)+2.0*S.cz;
          break;
        case REGULAR:
         HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                        +T.ax*T.ax)+2.0*S.ax;
          HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.ay+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.y)/SQR(SQR(v)))
                        +T.ax*T.ay)+S.ay;
          HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.az)+S.az;
          HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.ax*T.by);
          HessianMatrix.element[n][n+4]+=DDF*CosTheta*T.ax*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.bz);
          HessianMatrix.element[n][n+5]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.cz);

          HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.ay*CosTheta*T.ay+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.y)/SQR(SQR(v)))
                        +T.ay*T.ay)+0.5*(S.ax+S.by);
          HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.ay*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.az)+0.5*(S.bz);
          HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.ay*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.ay*T.by)+S.ay;
          HessianMatrix.element[n+1][n+4]+=DDF*CosTheta*T.ay*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.bz)+0.5*S.az;
          HessianMatrix.element[n+1][n+5]+=DDF*CosTheta*T.ay*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.cz);

          HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.az*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.az)+0.5*(S.ax+S.cz);
          HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.az*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.az*T.by);
          HessianMatrix.element[n+2][n+4]+=DDF*CosTheta*T.az*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.bz)+0.5*S.ay;
          HessianMatrix.element[n+2][n+5]+=DDF*CosTheta*T.az*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.cz)+S.az;

          HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.by*T.by)+2.0*S.by;
          HessianMatrix.element[n+3][n+4]+=DDF*CosTheta*T.by*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.by*T.bz)+S.bz;
          HessianMatrix.element[n+3][n+5]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.by*T.cz);

          HessianMatrix.element[n+4][n+4]+=DDF*CosTheta*T.bz*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.bz*T.bz)+0.5*(S.by+S.cz);
          HessianMatrix.element[n+4][n+5]+=DDF*CosTheta*T.bz*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.bz*T.cz)+S.bz;

          HessianMatrix.element[n+5][n+5]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.cz*T.cz)+2.0*S.cz;
          break;
        case REGULAR_UPPER_TRIANGLE:
          HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                        +T.ax*T.ax)+2.0*S.ax;
          HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.ay+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.y)/SQR(SQR(v)))
                        +T.ax*T.ay)+S.ay;
          HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.az)+S.az;
          HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.ax*T.by);
          HessianMatrix.element[n][n+4]+=DDF*CosTheta*T.ax*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.bz);
          HessianMatrix.element[n][n+5]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.cz);

          HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.ay*CosTheta*T.ay+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.y)/SQR(SQR(v)))
                        +T.ay*T.ay)+S.by;
          HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.ay*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.az)+S.bz;
          HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.ay*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.ay*T.by);
          HessianMatrix.element[n+1][n+4]+=DDF*CosTheta*T.ay*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.bz);
          HessianMatrix.element[n+1][n+5]+=DDF*CosTheta*T.ay*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.cz);

          HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.az*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.az)+S.cz;
          HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.az*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.az*T.by);
          HessianMatrix.element[n+2][n+4]+=DDF*CosTheta*T.az*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.bz);
          HessianMatrix.element[n+2][n+5]+=DDF*CosTheta*T.az*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.cz);

          HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.by*T.by)+2.0*S.by;
          HessianMatrix.element[n+3][n+4]+=DDF*CosTheta*T.by*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.by*T.bz)+S.bz;
          HessianMatrix.element[n+3][n+5]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.by*T.cz);

          HessianMatrix.element[n+4][n+4]+=DDF*CosTheta*T.bz*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.bz*T.bz)+S.cz;
          HessianMatrix.element[n+4][n+5]+=DDF*CosTheta*T.bz*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.bz*T.cz);

          HessianMatrix.element[n+5][n+5]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.cz*T.cz)+2.0*S.cz;
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.bz);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.by*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.bz)+S.bz;
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.bz*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.bz*T.bz)+0.5*(S.by+S.cz);
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.bz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.bz*T.cz)+S.bz;

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.az+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.az)+S.az;
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.az*CosTheta*T.az+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.x*vec_v.z)/SQR(SQR(v)))
                            +T.az*T.az)+0.5*(S.ax+S.cz);
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.az*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.az*T.by);
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.az*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.az*T.cz)+S.az;

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.ay+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.ay)+S.ay;
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.ay*CosTheta*T.ay+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.y)/SQR(SQR(v)))
                            +T.ay*T.ay)+0.5*(S.ax+S.by);
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.ay*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ay*T.by)+S.ay;
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.ay*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ay*T.cz);

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.bz);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.by*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.bz)+S.bz;
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.bz*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.bz*T.bz)+S.cz;
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.bz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.bz*T.cz);

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.az+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.az)+S.az;
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.az*CosTheta*T.az+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.x*vec_v.z)/SQR(SQR(v)))
                            +T.az*T.az)+S.cz;
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.az*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.az*T.by);
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.az*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.az*T.cz);

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.ay+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.ay)+S.ay;
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.ay*CosTheta*T.ay+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.y)/SQR(SQR(v)))
                            +T.ay*T.ay)+S.by;
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.ay*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ay*T.by);
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.ay*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ay*T.cz);

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
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






void HessianTorsionStrainPosition(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,INT_VECTOR3 index_k,INT_VECTOR3 index_l,
             VECTOR vec_u,VECTOR vec_v,VECTOR vec_w,REAL u,REAL v,REAL w,VECTOR fa,VECTOR fc,VECTOR fd,
             VECTOR Dab,VECTOR Dcb,VECTOR Ddc,REAL rbc,REAL_MATRIX3x3 D2I,REAL_MATRIX3x3 D2K,REAL_MATRIX3x3 D2L,
             REAL_MATRIX3x3 D2IJ,REAL_MATRIX3x3 D2IK,REAL_MATRIX3x3 D2IL,REAL_MATRIX3x3 D2JK,REAL_MATRIX3x3 D2JL,REAL_MATRIX3x3 D2KL,
             VECTOR dtA,VECTOR dtB,VECTOR dtC,VECTOR dtD,REAL DDF,REAL_MATRIX3x3 S,REAL CosPhi)
{
  int n;
  VECTOR veci,vecj,veck,vecl;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // contribution is zero (i.e. cancels out)
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          veci.x=Dab.x*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.x*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.x*(DDF*dtD.x*dtA.x+D2IL.ax)+fa.x;
          veci.y=Dab.x*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.x*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.x*(DDF*dtD.x*dtA.y+D2IL.bx);
          veci.z=Dab.x*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.x*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.x*(DDF*dtD.x*dtA.z+D2IL.cx);

          vecj.x=Dab.x*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.x*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.x*(DDF*dtD.x*dtB.x+D2JL.ax)-fa.x-fc.x-fd.x;
          vecj.y=Dab.x*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.x*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.x*(DDF*dtD.x*dtB.y+D2JL.bx);
          vecj.z=Dab.x*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.x*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.x*(DDF*dtD.x*dtB.z+D2JL.cx);

          veck.x=Dab.x*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.x*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.x*(DDF*dtD.x*dtC.x+D2KL.ax)+fc.x+fd.x-fd.x;
          veck.y=Dab.x*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.x*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.x*(DDF*dtD.x*dtC.y+D2KL.bx);
          veck.z=Dab.x*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.x*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.x*(DDF*dtD.x*dtC.z+D2KL.cx);

          vecl.x=Dab.x*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.x*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.x*(DDF*dtD.x*dtD.x+D2L.ax)+fd.x;
          vecl.y=Dab.x*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.x*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.x*(DDF*dtD.x*dtD.y+D2L.ay);
          vecl.z=Dab.x*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.x*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.x*(DDF*dtD.x*dtD.z+D2L.az);

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=veci.z;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=vecj.x;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=vecj.y;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=vecj.z;
          if(index_k.x>=0) HessianMatrix.element[index_k.x][n]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n]+=veck.z;
          if(index_l.x>=0) HessianMatrix.element[index_l.x][n]+=vecl.x;
          if(index_l.y>=0) HessianMatrix.element[index_l.y][n]+=vecl.y;
          if(index_l.z>=0) HessianMatrix.element[index_l.z][n]+=vecl.z;

          veci.x=Dab.y*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.y*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.y*(DDF*dtD.y*dtA.x+D2IL.ay);
          veci.y=Dab.y*(DDF*dtA.y*dtA.y+D2I.by)+fa.y+Dcb.y*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.y*(DDF*dtD.y*dtA.y+D2IL.by);
          veci.z=Dab.y*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.y*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.y*(DDF*dtD.y*dtA.z+D2IL.cy);

          vecj.x=Dab.y*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.y*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.y*(DDF*dtD.y*dtB.x+D2JL.ay);
          vecj.y=Dab.y*(DDF*dtA.y*dtB.y+D2IJ.by)-fa.y+Dcb.y*rbc*(DDF*dtC.y*dtB.y+D2JK.by)-fc.y+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.y+D2JL.by)-fd.y+Ddc.y*(DDF*dtD.y*dtB.y+D2JL.by);
          vecj.z=Dab.y*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.y*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.y*(DDF*dtD.y*dtB.z+D2JL.cy);

          veck.x=Dab.y*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.y*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.y*(DDF*dtD.y*dtC.x+D2KL.ay);
          veck.y=Dab.y*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.y*rbc*(DDF*dtC.y*dtC.y+D2K.by)+fc.y+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+fd.y+Ddc.y*(DDF*dtD.y*dtC.y+D2KL.by)-fd.y;
          veck.z=Dab.y*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.y*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.y*(DDF*dtD.y*dtC.z+D2KL.cy);

          vecl.x=Dab.y*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.y*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.y*(DDF*dtD.y*dtD.x+D2L.ay);
          vecl.y=Dab.y*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.y*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.y*(DDF*dtD.y*dtD.y+D2L.by)+fd.y;
          vecl.z=Dab.y*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.y*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.y*(DDF*dtD.y*dtD.z+D2L.bz);

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=veci.z;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=vecj.x;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=vecj.y;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=vecj.z;
          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+1]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+1]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+1]+=veck.z;
          if(index_l.x>=0) HessianMatrix.element[index_l.x][n+1]+=vecl.x;
          if(index_l.y>=0) HessianMatrix.element[index_l.y][n+1]+=vecl.y;
          if(index_l.z>=0) HessianMatrix.element[index_l.z][n+1]+=vecl.z;

          veci.x=Dab.z*(DDF*dtA.z*dtA.x+D2I.az)+Dcb.z*rbc*(DDF*dtC.z*dtA.x+D2IK.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.x+D2IL.az)+Ddc.z*(DDF*dtD.z*dtA.x+D2IL.az);
          veci.y=Dab.z*(DDF*dtA.z*dtA.y+D2I.bz)+Dcb.z*rbc*(DDF*dtC.z*dtA.y+D2IK.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.y+D2IL.bz)+Ddc.z*(DDF*dtD.z*dtA.y+D2IL.bz);
          veci.z=Dab.z*(DDF*dtA.z*dtA.z+D2I.cz)+fa.z+Dcb.z*rbc*(DDF*dtC.z*dtA.z+D2IK.cz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.z+D2IL.cz)+Ddc.z*(DDF*dtD.z*dtA.z+D2IL.cz);

          vecj.x=Dab.z*(DDF*dtA.z*dtB.x+D2IJ.cx)+Dcb.z*rbc*(DDF*dtC.z*dtB.x+D2JK.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.x+D2JL.az)+Ddc.z*(DDF*dtD.z*dtB.x+D2JL.az);
          vecj.y=Dab.z*(DDF*dtA.z*dtB.y+D2IJ.cy)+Dcb.z*rbc*(DDF*dtC.z*dtB.y+D2JK.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.y+D2JL.bz)+Ddc.z*(DDF*dtD.z*dtB.y+D2JL.bz);
          vecj.z=Dab.z*(DDF*dtA.z*dtB.z+D2IJ.cz)-fa.z+Dcb.z*rbc*(DDF*dtC.z*dtB.z+D2JK.cz)-fc.z+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.z+D2JL.cz)-fd.z+Ddc.z*(DDF*dtD.z*dtB.z+D2JL.cz);

          veck.x=Dab.z*(DDF*dtA.z*dtC.x+D2IK.cx)+Dcb.z*rbc*(DDF*dtC.z*dtC.x+D2K.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.x+D2KL.az)+Ddc.z*(DDF*dtD.z*dtC.x+D2KL.az);
          veck.y=Dab.z*(DDF*dtA.z*dtC.y+D2IK.cy)+Dcb.z*rbc*(DDF*dtC.z*dtC.y+D2K.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.y+D2KL.bz)+Ddc.z*(DDF*dtD.z*dtC.y+D2KL.bz);
          veck.z=Dab.z*(DDF*dtA.z*dtC.z+D2IK.cz)+Dcb.z*rbc*(DDF*dtC.z*dtC.z+D2K.cz)+fc.z+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.z+D2KL.cz)+fd.z+Ddc.z*(DDF*dtD.z*dtC.z+D2KL.cz)-fd.z;

          vecl.x=Dab.z*(DDF*dtA.z*dtD.x+D2IL.cx)+Dcb.z*rbc*(DDF*dtC.z*dtD.x+D2KL.cx)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.x+D2L.az)+Ddc.z*(DDF*dtD.z*dtD.x+D2L.az);
          vecl.y=Dab.z*(DDF*dtA.z*dtD.y+D2IL.cy)+Dcb.z*rbc*(DDF*dtC.z*dtD.y+D2KL.cy)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.y+D2L.bz)+Ddc.z*(DDF*dtD.z*dtD.y+D2L.bz);
          vecl.z=Dab.z*(DDF*dtA.z*dtD.z+D2IL.cz)+Dcb.z*rbc*(DDF*dtC.z*dtD.z+D2KL.cz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.z+D2L.cz)+Ddc.z*(DDF*dtD.z*dtD.z+D2L.cz)+fd.z;

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=veci.z;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=vecj.x;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=vecj.y;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=vecj.z;
          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+2]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+2]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+2]+=veck.z;
          if(index_l.x>=0) HessianMatrix.element[index_l.x][n+2]+=vecl.x;
          if(index_l.y>=0) HessianMatrix.element[index_l.y][n+2]+=vecl.y;
          if(index_l.z>=0) HessianMatrix.element[index_l.z][n+2]+=vecl.z;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          //S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
          veci.x=Dab.x*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.x*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.x*(DDF*dtD.x*dtA.x+D2IL.ax)+fa.x;
          veci.y=Dab.x*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.x*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.x*(DDF*dtD.x*dtA.y+D2IL.bx);
          veci.z=Dab.x*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.x*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.x*(DDF*dtD.x*dtA.z+D2IL.cx);

          vecj.x=Dab.x*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.x*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.x*(DDF*dtD.x*dtB.x+D2JL.ax)-fa.x-fc.x-fd.x;
          vecj.y=Dab.x*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.x*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.x*(DDF*dtD.x*dtB.y+D2JL.bx);
          vecj.z=Dab.x*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.x*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.x*(DDF*dtD.x*dtB.z+D2JL.cx);

          veck.x=Dab.x*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.x*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.x*(DDF*dtD.x*dtC.x+D2KL.ax)+fc.x+fd.x-fd.x;
          veck.y=Dab.x*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.x*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.x*(DDF*dtD.x*dtC.y+D2KL.bx);
          veck.z=Dab.x*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.x*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.x*(DDF*dtD.x*dtC.z+D2KL.cx);

          vecl.x=Dab.x*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.x*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.x*(DDF*dtD.x*dtD.x+D2L.ax)+fd.x;
          vecl.y=Dab.x*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.x*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.x*(DDF*dtD.x*dtD.y+D2L.ay);
          vecl.z=Dab.x*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.x*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.x*(DDF*dtD.x*dtD.z+D2L.az);

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=veci.z;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=vecj.x;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=vecj.y;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=vecj.z;
          if(index_k.x>=0) HessianMatrix.element[index_k.x][n]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n]+=veck.z;
          if(index_l.x>=0) HessianMatrix.element[index_l.x][n]+=vecl.x;
          if(index_l.y>=0) HessianMatrix.element[index_l.y][n]+=vecl.y;
          if(index_l.z>=0) HessianMatrix.element[index_l.z][n]+=vecl.z;

          //S.bx=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
          veci.x=Dab.y*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.y*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                 Dcb.y*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.y*(DDF*dtD.x*dtA.x+D2IL.ax);
          veci.y=Dab.y*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.y*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                 Dcb.y*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.y*(DDF*dtD.x*dtA.y+D2IL.bx)+fa.x;
          veci.z=Dab.y*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.y*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                 Dcb.y*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.y*(DDF*dtD.x*dtA.z+D2IL.cx);

          vecj.x=Dab.y*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.y*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                 Dcb.y*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.y*(DDF*dtD.x*dtB.x+D2JL.ax);
          vecj.y=Dab.y*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.y*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                 Dcb.y*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.y*(DDF*dtD.x*dtB.y+D2JL.bx)-fa.x-fc.x-fd.x;
          vecj.z=Dab.y*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.y*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                 Dcb.y*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.y*(DDF*dtD.x*dtB.z+D2JL.cx);

          veck.x=Dab.y*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.y*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                 Dcb.y*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.y*(DDF*dtD.x*dtC.x+D2KL.ax);
          veck.y=Dab.y*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.y*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                 Dcb.y*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.y*(DDF*dtD.x*dtC.y+D2KL.bx)+fc.x+fd.x-fd.x;
          veck.z=Dab.y*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.y*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                 Dcb.y*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.y*(DDF*dtD.x*dtC.z+D2KL.cx);

          vecl.x=Dab.y*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.y*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                 Dcb.y*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.y*(DDF*dtD.x*dtD.x+D2L.ax);
          vecl.y=Dab.y*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.y*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                 Dcb.y*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.y*(DDF*dtD.x*dtD.y+D2L.ay)+fd.x;
          vecl.z=Dab.y*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.y*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                 Dcb.y*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.y*(DDF*dtD.x*dtD.z+D2L.az);

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=veci.z;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=vecj.x;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=vecj.y;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=vecj.z;
          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+1]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+1]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+1]+=veck.z;
          if(index_l.x>=0) HessianMatrix.element[index_l.x][n+1]+=vecl.x;
          if(index_l.y>=0) HessianMatrix.element[index_l.y][n+1]+=vecl.y;
          if(index_l.z>=0) HessianMatrix.element[index_l.z][n+1]+=vecl.z;


          veci.x=Dab.z*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.z*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                Dcb.z*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.z*(DDF*dtD.x*dtA.x+D2IL.ax);
          veci.y=Dab.z*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.z*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                Dcb.z*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.z*(DDF*dtD.x*dtA.y+D2IL.bx);
          veci.z=Dab.z*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.z*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                 Dcb.z*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.z*(DDF*dtD.x*dtA.z+D2IL.cx)+fa.x;

          vecj.x=Dab.z*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.z*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                 Dcb.z*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.z*(DDF*dtD.x*dtB.x+D2JL.ax);
          vecj.y=Dab.z*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.z*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                 Dcb.z*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.z*(DDF*dtD.x*dtB.y+D2JL.bx);
          vecj.z=Dab.z*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.z*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                 Dcb.z*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.z*(DDF*dtD.x*dtB.z+D2JL.cx)-fa.x-fc.x-fd.x;

          veck.x=Dab.z*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.z*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                 Dcb.z*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.z*(DDF*dtD.x*dtC.x+D2KL.ax);
          veck.y=Dab.z*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.z*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                 Dcb.z*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.z*(DDF*dtD.x*dtC.y+D2KL.bx);
          veck.z=Dab.z*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.z*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                 Dcb.z*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.z*(DDF*dtD.x*dtC.z+D2KL.cx)+fc.x+fd.x-fd.x;

          vecl.x=Dab.z*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.z*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                 Dcb.z*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.z*(DDF*dtD.x*dtD.x+D2L.ax);
          vecl.y=Dab.z*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.z*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                 Dcb.z*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.z*(DDF*dtD.x*dtD.y+D2L.ay);
          vecl.z=Dab.z*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.z*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                 Dcb.z*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.z*(DDF*dtD.x*dtD.z+D2L.az)+fd.x;

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=veci.z;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=vecj.x;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=vecj.y;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=vecj.z;
          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+2]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+2]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+2]+=veck.z;
          if(index_l.x>=0) HessianMatrix.element[index_l.x][n+2]+=vecl.x;
          if(index_l.y>=0) HessianMatrix.element[index_l.y][n+2]+=vecl.y;
          if(index_l.z>=0) HessianMatrix.element[index_l.z][n+2]+=vecl.z;


          veci.x=Dab.y*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.y*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.y*(DDF*dtD.y*dtA.x+D2IL.ay);
          veci.y=Dab.y*(DDF*dtA.y*dtA.y+D2I.by)+fa.y+Dcb.y*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.y*(DDF*dtD.y*dtA.y+D2IL.by);
          veci.z=Dab.y*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.y*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.y*(DDF*dtD.y*dtA.z+D2IL.cy);

          vecj.x=Dab.y*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.y*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.y*(DDF*dtD.y*dtB.x+D2JL.ay);
          vecj.y=Dab.y*(DDF*dtA.y*dtB.y+D2IJ.by)-fa.y+Dcb.y*rbc*(DDF*dtC.y*dtB.y+D2JK.by)-fc.y+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.y+D2JL.by)-fd.y+Ddc.y*(DDF*dtD.y*dtB.y+D2JL.by);
          vecj.z=Dab.y*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.y*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.y*(DDF*dtD.y*dtB.z+D2JL.cy);

          veck.x=Dab.y*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.y*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.y*(DDF*dtD.y*dtC.x+D2KL.ay);
          veck.y=Dab.y*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.y*rbc*(DDF*dtC.y*dtC.y+D2K.by)+fc.y+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+fd.y+Ddc.y*(DDF*dtD.y*dtC.y+D2KL.by)-fd.y;
          veck.z=Dab.y*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.y*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.y*(DDF*dtD.y*dtC.z+D2KL.cy);

          vecl.x=Dab.y*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.y*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.y*(DDF*dtD.y*dtD.x+D2L.ay);
          vecl.y=Dab.y*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.y*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.y*(DDF*dtD.y*dtD.y+D2L.by)+fd.y;
          vecl.z=Dab.y*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.y*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.y*(DDF*dtD.y*dtD.z+D2L.bz);

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=veci.z;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=vecj.x;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=vecj.y;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=vecj.z;
          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+3]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+3]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+3]+=veck.z;
          if(index_l.x>=0) HessianMatrix.element[index_l.x][n+3]+=vecl.x;
          if(index_l.y>=0) HessianMatrix.element[index_l.y][n+3]+=vecl.y;
          if(index_l.z>=0) HessianMatrix.element[index_l.z][n+3]+=vecl.z;


          veci.x=Dab.z*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.z*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                 Dcb.z*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.z*(DDF*dtD.y*dtA.x+D2IL.ay);
          veci.y=Dab.z*(DDF*dtA.y*dtA.y+D2I.by)+Dcb.z*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                 Dcb.z*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.z*(DDF*dtD.y*dtA.y+D2IL.by);
          veci.z=Dab.z*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.z*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                 Dcb.z*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.z*(DDF*dtD.y*dtA.z+D2IL.cy)+fa.y;

          vecj.x=Dab.z*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.z*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                 Dcb.z*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.z*(DDF*dtD.y*dtB.x+D2JL.ay);
          vecj.y=Dab.z*(DDF*dtA.y*dtB.y+D2IJ.by)+Dcb.z*rbc*(DDF*dtC.y*dtB.y+D2JK.by)+
                 Dcb.z*rbc*(DDF*dtD.y*dtB.y+D2JL.by)+Ddc.z*(DDF*dtD.y*dtB.y+D2JL.by);
          vecj.z=Dab.z*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.z*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                 Dcb.z*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.z*(DDF*dtD.y*dtB.z+D2JL.cy)-fa.y-fc.y-fd.y;

          veck.x=Dab.z*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.z*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                 Dcb.z*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.z*(DDF*dtD.y*dtC.x+D2KL.ay);
          veck.y=Dab.z*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.z*rbc*(DDF*dtC.y*dtC.y+D2K.by)+
                 Dcb.z*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+Ddc.z*(DDF*dtD.y*dtC.y+D2KL.by);
          veck.z=Dab.z*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.z*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                 Dcb.z*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.z*(DDF*dtD.y*dtC.z+D2KL.cy)+fc.y+fd.y-fd.y;

          vecl.x=Dab.z*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.z*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                 Dcb.z*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.z*(DDF*dtD.y*dtD.x+D2L.ay);
          vecl.y=Dab.z*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.z*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                 Dcb.z*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.z*(DDF*dtD.y*dtD.y+D2L.by);
          vecl.z=Dab.z*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.z*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                 Dcb.z*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.z*(DDF*dtD.y*dtD.z+D2L.bz)+fd.y;

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=veci.z;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]+=vecj.x;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]+=vecj.y;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]+=vecj.z;
          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+4]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+4]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+4]+=veck.z;
          if(index_l.x>=0) HessianMatrix.element[index_l.x][n+4]+=vecl.x;
          if(index_l.y>=0) HessianMatrix.element[index_l.y][n+4]+=vecl.y;
          if(index_l.z>=0) HessianMatrix.element[index_l.z][n+4]+=vecl.z;

          veci.x=Dab.z*(DDF*dtA.z*dtA.x+D2I.az)+Dcb.z*rbc*(DDF*dtC.z*dtA.x+D2IK.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.x+D2IL.az)+Ddc.z*(DDF*dtD.z*dtA.x+D2IL.az);
          veci.y=Dab.z*(DDF*dtA.z*dtA.y+D2I.bz)+Dcb.z*rbc*(DDF*dtC.z*dtA.y+D2IK.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.y+D2IL.bz)+Ddc.z*(DDF*dtD.z*dtA.y+D2IL.bz);
          veci.z=Dab.z*(DDF*dtA.z*dtA.z+D2I.cz)+fa.z+Dcb.z*rbc*(DDF*dtC.z*dtA.z+D2IK.cz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.z+D2IL.cz)+Ddc.z*(DDF*dtD.z*dtA.z+D2IL.cz);

          vecj.x=Dab.z*(DDF*dtA.z*dtB.x+D2IJ.cx)+Dcb.z*rbc*(DDF*dtC.z*dtB.x+D2JK.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.x+D2JL.az)+Ddc.z*(DDF*dtD.z*dtB.x+D2JL.az);
          vecj.y=Dab.z*(DDF*dtA.z*dtB.y+D2IJ.cy)+Dcb.z*rbc*(DDF*dtC.z*dtB.y+D2JK.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.y+D2JL.bz)+Ddc.z*(DDF*dtD.z*dtB.y+D2JL.bz);
          vecj.z=Dab.z*(DDF*dtA.z*dtB.z+D2IJ.cz)-fa.z+Dcb.z*rbc*(DDF*dtC.z*dtB.z+D2JK.cz)-fc.z+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.z+D2JL.cz)-fd.z+Ddc.z*(DDF*dtD.z*dtB.z+D2JL.cz);

          veck.x=Dab.z*(DDF*dtA.z*dtC.x+D2IK.cx)+Dcb.z*rbc*(DDF*dtC.z*dtC.x+D2K.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.x+D2KL.az)+Ddc.z*(DDF*dtD.z*dtC.x+D2KL.az);
          veck.y=Dab.z*(DDF*dtA.z*dtC.y+D2IK.cy)+Dcb.z*rbc*(DDF*dtC.z*dtC.y+D2K.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.y+D2KL.bz)+Ddc.z*(DDF*dtD.z*dtC.y+D2KL.bz);
          veck.z=Dab.z*(DDF*dtA.z*dtC.z+D2IK.cz)+Dcb.z*rbc*(DDF*dtC.z*dtC.z+D2K.cz)+fc.z+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.z+D2KL.cz)+fd.z+Ddc.z*(DDF*dtD.z*dtC.z+D2KL.cz)-fd.z;

          vecl.x=Dab.z*(DDF*dtA.z*dtD.x+D2IL.cx)+Dcb.z*rbc*(DDF*dtC.z*dtD.x+D2KL.cx)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.x+D2L.az)+Ddc.z*(DDF*dtD.z*dtD.x+D2L.az);
          vecl.y=Dab.z*(DDF*dtA.z*dtD.y+D2IL.cy)+Dcb.z*rbc*(DDF*dtC.z*dtD.y+D2KL.cy)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.y+D2L.bz)+Ddc.z*(DDF*dtD.z*dtD.y+D2L.bz);
          vecl.z=Dab.z*(DDF*dtA.z*dtD.z+D2IL.cz)+Dcb.z*rbc*(DDF*dtC.z*dtD.z+D2KL.cz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.z+D2L.cz)+Ddc.z*(DDF*dtD.z*dtD.z+D2L.cz)+fd.z;

          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=veci.x;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=veci.y;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=veci.z;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]+=vecj.x;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]+=vecj.y;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]+=vecj.z;
          if(index_k.x>=0) HessianMatrix.element[index_k.x][n+5]+=veck.x;
          if(index_k.y>=0) HessianMatrix.element[index_k.y][n+5]+=veck.y;
          if(index_k.z>=0) HessianMatrix.element[index_k.z][n+5]+=veck.z;
          if(index_l.x>=0) HessianMatrix.element[index_l.x][n+5]+=vecl.x;
          if(index_l.y>=0) HessianMatrix.element[index_l.y][n+5]+=vecl.y;
          if(index_l.z>=0) HessianMatrix.element[index_l.z][n+5]+=vecl.z;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              //S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
              veci.x=Dab.x*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.x*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.x*(DDF*dtD.x*dtA.x+D2IL.ax)+fa.x;
              veci.y=Dab.x*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.x*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.x*(DDF*dtD.x*dtA.y+D2IL.bx);
              veci.z=Dab.x*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.x*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.x*(DDF*dtD.x*dtA.z+D2IL.cx);

              vecj.x=Dab.x*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.x*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.x*(DDF*dtD.x*dtB.x+D2JL.ax)-fa.x-fc.x-fd.x;
              vecj.y=Dab.x*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.x*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.x*(DDF*dtD.x*dtB.y+D2JL.bx);
              vecj.z=Dab.x*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.x*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.x*(DDF*dtD.x*dtB.z+D2JL.cx);

              veck.x=Dab.x*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.x*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.x*(DDF*dtD.x*dtC.x+D2KL.ax)+fc.x+fd.x-fd.x;
              veck.y=Dab.x*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.x*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.x*(DDF*dtD.x*dtC.y+D2KL.bx);
              veck.z=Dab.x*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.x*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.x*(DDF*dtD.x*dtC.z+D2KL.cx);

              vecl.x=Dab.x*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.x*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.x*(DDF*dtD.x*dtD.x+D2L.ax)+fd.x;
              vecl.y=Dab.x*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.x*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.x*(DDF*dtD.x*dtD.y+D2L.ay);
              vecl.z=Dab.x*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.x*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.x*(DDF*dtD.x*dtD.z+D2L.az);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n]+=vecl.z;

              veci.x=Dab.y*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.y*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.y*(DDF*dtD.y*dtA.x+D2IL.ay);
              veci.y=Dab.y*(DDF*dtA.y*dtA.y+D2I.by)+fa.y+Dcb.y*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.y*(DDF*dtD.y*dtA.y+D2IL.by);
              veci.z=Dab.y*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.y*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.y*(DDF*dtD.y*dtA.z+D2IL.cy);

              vecj.x=Dab.y*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.y*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.y*(DDF*dtD.y*dtB.x+D2JL.ay);
              vecj.y=Dab.y*(DDF*dtA.y*dtB.y+D2IJ.by)-fa.y+Dcb.y*rbc*(DDF*dtC.y*dtB.y+D2JK.by)-fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.y+D2JL.by)-fd.y+Ddc.y*(DDF*dtD.y*dtB.y+D2JL.by);
              vecj.z=Dab.y*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.y*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.y*(DDF*dtD.y*dtB.z+D2JL.cy);

              veck.x=Dab.y*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.y*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.y*(DDF*dtD.y*dtC.x+D2KL.ay);
              veck.y=Dab.y*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.y*rbc*(DDF*dtC.y*dtC.y+D2K.by)+fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+fd.y+Ddc.y*(DDF*dtD.y*dtC.y+D2KL.by)-fd.y;
              veck.z=Dab.y*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.y*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.y*(DDF*dtD.y*dtC.z+D2KL.cy);

              vecl.x=Dab.y*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.y*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.y*(DDF*dtD.y*dtD.x+D2L.ay);
              vecl.y=Dab.y*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.y*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.y*(DDF*dtD.y*dtD.y+D2L.by)+fd.y;
              vecl.z=Dab.y*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.y*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.y*(DDF*dtD.y*dtD.z+D2L.bz);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+1]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+1]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+1]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n+1]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n+1]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n+1]+=vecl.z;


              veci.x=Dab.z*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.z*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                     Dcb.z*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.z*(DDF*dtD.y*dtA.x+D2IL.ay);
              veci.y=Dab.z*(DDF*dtA.y*dtA.y+D2I.by)+Dcb.z*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                     Dcb.z*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.z*(DDF*dtD.y*dtA.y+D2IL.by);
              veci.z=Dab.z*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.z*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                     Dcb.z*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.z*(DDF*dtD.y*dtA.z+D2IL.cy)+fa.y;

              vecj.x=Dab.z*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.z*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                     Dcb.z*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.z*(DDF*dtD.y*dtB.x+D2JL.ay);
              vecj.y=Dab.z*(DDF*dtA.y*dtB.y+D2IJ.by)+Dcb.z*rbc*(DDF*dtC.y*dtB.y+D2JK.by)+
                     Dcb.z*rbc*(DDF*dtD.y*dtB.y+D2JL.by)+Ddc.z*(DDF*dtD.y*dtB.y+D2JL.by);
              vecj.z=Dab.z*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.z*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                     Dcb.z*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.z*(DDF*dtD.y*dtB.z+D2JL.cy)-fa.y-fc.y-fd.y;

              veck.x=Dab.z*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.z*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                     Dcb.z*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.z*(DDF*dtD.y*dtC.x+D2KL.ay);
              veck.y=Dab.z*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.z*rbc*(DDF*dtC.y*dtC.y+D2K.by)+
                     Dcb.z*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+Ddc.z*(DDF*dtD.y*dtC.y+D2KL.by);
              veck.z=Dab.z*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.z*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                     Dcb.z*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.z*(DDF*dtD.y*dtC.z+D2KL.cy)+fc.y+fd.y-fd.y;

              vecl.x=Dab.z*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.z*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                     Dcb.z*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.z*(DDF*dtD.y*dtD.x+D2L.ay);
              vecl.y=Dab.z*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.z*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                     Dcb.z*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.z*(DDF*dtD.y*dtD.y+D2L.by);
              vecl.z=Dab.z*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.z*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                     Dcb.z*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.z*(DDF*dtD.y*dtD.z+D2L.bz)+fd.y;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+2]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+2]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+2]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n+2]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n+2]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n+2]+=vecl.z;

              veci.x=Dab.z*(DDF*dtA.z*dtA.x+D2I.az)+Dcb.z*rbc*(DDF*dtC.z*dtA.x+D2IK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.x+D2IL.az)+Ddc.z*(DDF*dtD.z*dtA.x+D2IL.az);
              veci.y=Dab.z*(DDF*dtA.z*dtA.y+D2I.bz)+Dcb.z*rbc*(DDF*dtC.z*dtA.y+D2IK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.y+D2IL.bz)+Ddc.z*(DDF*dtD.z*dtA.y+D2IL.bz);
              veci.z=Dab.z*(DDF*dtA.z*dtA.z+D2I.cz)+fa.z+Dcb.z*rbc*(DDF*dtC.z*dtA.z+D2IK.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.z+D2IL.cz)+Ddc.z*(DDF*dtD.z*dtA.z+D2IL.cz);

              vecj.x=Dab.z*(DDF*dtA.z*dtB.x+D2IJ.cx)+Dcb.z*rbc*(DDF*dtC.z*dtB.x+D2JK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.x+D2JL.az)+Ddc.z*(DDF*dtD.z*dtB.x+D2JL.az);
              vecj.y=Dab.z*(DDF*dtA.z*dtB.y+D2IJ.cy)+Dcb.z*rbc*(DDF*dtC.z*dtB.y+D2JK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.y+D2JL.bz)+Ddc.z*(DDF*dtD.z*dtB.y+D2JL.bz);
              vecj.z=Dab.z*(DDF*dtA.z*dtB.z+D2IJ.cz)-fa.z+Dcb.z*rbc*(DDF*dtC.z*dtB.z+D2JK.cz)-fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.z+D2JL.cz)-fd.z+Ddc.z*(DDF*dtD.z*dtB.z+D2JL.cz);

              veck.x=Dab.z*(DDF*dtA.z*dtC.x+D2IK.cx)+Dcb.z*rbc*(DDF*dtC.z*dtC.x+D2K.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.x+D2KL.az)+Ddc.z*(DDF*dtD.z*dtC.x+D2KL.az);
              veck.y=Dab.z*(DDF*dtA.z*dtC.y+D2IK.cy)+Dcb.z*rbc*(DDF*dtC.z*dtC.y+D2K.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.y+D2KL.bz)+Ddc.z*(DDF*dtD.z*dtC.y+D2KL.bz);
              veck.z=Dab.z*(DDF*dtA.z*dtC.z+D2IK.cz)+Dcb.z*rbc*(DDF*dtC.z*dtC.z+D2K.cz)+fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.z+D2KL.cz)+fd.z+Ddc.z*(DDF*dtD.z*dtC.z+D2KL.cz)-fd.z;

              vecl.x=Dab.z*(DDF*dtA.z*dtD.x+D2IL.cx)+Dcb.z*rbc*(DDF*dtC.z*dtD.x+D2KL.cx)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.x+D2L.az)+Ddc.z*(DDF*dtD.z*dtD.x+D2L.az);
              vecl.y=Dab.z*(DDF*dtA.z*dtD.y+D2IL.cy)+Dcb.z*rbc*(DDF*dtC.z*dtD.y+D2KL.cy)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.y+D2L.bz)+Ddc.z*(DDF*dtD.z*dtD.y+D2L.bz);
              vecl.z=Dab.z*(DDF*dtA.z*dtD.z+D2IL.cz)+Dcb.z*rbc*(DDF*dtC.z*dtD.z+D2KL.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.z+D2L.cz)+Ddc.z*(DDF*dtD.z*dtD.z+D2L.cz)+fd.z;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+3]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+3]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+3]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n+3]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n+3]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n+3]+=vecl.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              //S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
              veci.x=Dab.x*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.x*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.x*(DDF*dtD.x*dtA.x+D2IL.ax)+fa.x;
              veci.y=Dab.x*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.x*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.x*(DDF*dtD.x*dtA.y+D2IL.bx);
              veci.z=Dab.x*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.x*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.x*(DDF*dtD.x*dtA.z+D2IL.cx);

              vecj.x=Dab.x*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.x*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.x*(DDF*dtD.x*dtB.x+D2JL.ax)-fa.x-fc.x-fd.x;
              vecj.y=Dab.x*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.x*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.x*(DDF*dtD.x*dtB.y+D2JL.bx);
              vecj.z=Dab.x*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.x*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.x*(DDF*dtD.x*dtB.z+D2JL.cx);

              veck.x=Dab.x*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.x*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.x*(DDF*dtD.x*dtC.x+D2KL.ax)+fc.x+fd.x-fd.x;
              veck.y=Dab.x*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.x*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.x*(DDF*dtD.x*dtC.y+D2KL.bx);
              veck.z=Dab.x*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.x*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.x*(DDF*dtD.x*dtC.z+D2KL.cx);

              vecl.x=Dab.x*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.x*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.x*(DDF*dtD.x*dtD.x+D2L.ax)+fd.x;
              vecl.y=Dab.x*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.x*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.x*(DDF*dtD.x*dtD.y+D2L.ay);
              vecl.z=Dab.x*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.x*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.x*(DDF*dtD.x*dtD.z+D2L.az);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n]+=vecl.z;

              veci.x=Dab.z*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.z*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                    Dcb.z*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.z*(DDF*dtD.x*dtA.x+D2IL.ax);
              veci.y=Dab.z*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.z*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                    Dcb.z*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.z*(DDF*dtD.x*dtA.y+D2IL.bx);
              veci.z=Dab.z*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.z*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                     Dcb.z*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.z*(DDF*dtD.x*dtA.z+D2IL.cx)+fa.x;

              vecj.x=Dab.z*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.z*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                     Dcb.z*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.z*(DDF*dtD.x*dtB.x+D2JL.ax);
              vecj.y=Dab.z*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.z*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                     Dcb.z*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.z*(DDF*dtD.x*dtB.y+D2JL.bx);
              vecj.z=Dab.z*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.z*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                     Dcb.z*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.z*(DDF*dtD.x*dtB.z+D2JL.cx)-fa.x-fc.x-fd.x;

              veck.x=Dab.z*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.z*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                     Dcb.z*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.z*(DDF*dtD.x*dtC.x+D2KL.ax);
              veck.y=Dab.z*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.z*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                     Dcb.z*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.z*(DDF*dtD.x*dtC.y+D2KL.bx);
              veck.z=Dab.z*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.z*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                     Dcb.z*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.z*(DDF*dtD.x*dtC.z+D2KL.cx)+fc.x+fd.x-fd.x;

              vecl.x=Dab.z*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.z*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                     Dcb.z*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.z*(DDF*dtD.x*dtD.x+D2L.ax);
              vecl.y=Dab.z*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.z*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                     Dcb.z*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.z*(DDF*dtD.x*dtD.y+D2L.ay);
              vecl.z=Dab.z*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.z*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                     Dcb.z*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.z*(DDF*dtD.x*dtD.z+D2L.az)+fd.x;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+1]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+1]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+1]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n+1]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n+1]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n+1]+=vecl.z;

              veci.x=Dab.y*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.y*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.y*(DDF*dtD.y*dtA.x+D2IL.ay);
              veci.y=Dab.y*(DDF*dtA.y*dtA.y+D2I.by)+fa.y+Dcb.y*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.y*(DDF*dtD.y*dtA.y+D2IL.by);
              veci.z=Dab.y*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.y*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.y*(DDF*dtD.y*dtA.z+D2IL.cy);

              vecj.x=Dab.y*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.y*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.y*(DDF*dtD.y*dtB.x+D2JL.ay);
              vecj.y=Dab.y*(DDF*dtA.y*dtB.y+D2IJ.by)-fa.y+Dcb.y*rbc*(DDF*dtC.y*dtB.y+D2JK.by)-fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.y+D2JL.by)-fd.y+Ddc.y*(DDF*dtD.y*dtB.y+D2JL.by);
              vecj.z=Dab.y*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.y*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.y*(DDF*dtD.y*dtB.z+D2JL.cy);

              veck.x=Dab.y*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.y*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.y*(DDF*dtD.y*dtC.x+D2KL.ay);
              veck.y=Dab.y*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.y*rbc*(DDF*dtC.y*dtC.y+D2K.by)+fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+fd.y+Ddc.y*(DDF*dtD.y*dtC.y+D2KL.by)-fd.y;
              veck.z=Dab.y*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.y*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.y*(DDF*dtD.y*dtC.z+D2KL.cy);

              vecl.x=Dab.y*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.y*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.y*(DDF*dtD.y*dtD.x+D2L.ay);
              vecl.y=Dab.y*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.y*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.y*(DDF*dtD.y*dtD.y+D2L.by)+fd.y;
              vecl.z=Dab.y*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.y*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.y*(DDF*dtD.y*dtD.z+D2L.bz);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+2]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+2]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+2]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n+2]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n+2]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n+2]+=vecl.z;


              veci.x=Dab.z*(DDF*dtA.z*dtA.x+D2I.az)+Dcb.z*rbc*(DDF*dtC.z*dtA.x+D2IK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.x+D2IL.az)+Ddc.z*(DDF*dtD.z*dtA.x+D2IL.az);
              veci.y=Dab.z*(DDF*dtA.z*dtA.y+D2I.bz)+Dcb.z*rbc*(DDF*dtC.z*dtA.y+D2IK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.y+D2IL.bz)+Ddc.z*(DDF*dtD.z*dtA.y+D2IL.bz);
              veci.z=Dab.z*(DDF*dtA.z*dtA.z+D2I.cz)+fa.z+Dcb.z*rbc*(DDF*dtC.z*dtA.z+D2IK.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.z+D2IL.cz)+Ddc.z*(DDF*dtD.z*dtA.z+D2IL.cz);

              vecj.x=Dab.z*(DDF*dtA.z*dtB.x+D2IJ.cx)+Dcb.z*rbc*(DDF*dtC.z*dtB.x+D2JK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.x+D2JL.az)+Ddc.z*(DDF*dtD.z*dtB.x+D2JL.az);
              vecj.y=Dab.z*(DDF*dtA.z*dtB.y+D2IJ.cy)+Dcb.z*rbc*(DDF*dtC.z*dtB.y+D2JK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.y+D2JL.bz)+Ddc.z*(DDF*dtD.z*dtB.y+D2JL.bz);
              vecj.z=Dab.z*(DDF*dtA.z*dtB.z+D2IJ.cz)-fa.z+Dcb.z*rbc*(DDF*dtC.z*dtB.z+D2JK.cz)-fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.z+D2JL.cz)-fd.z+Ddc.z*(DDF*dtD.z*dtB.z+D2JL.cz);

              veck.x=Dab.z*(DDF*dtA.z*dtC.x+D2IK.cx)+Dcb.z*rbc*(DDF*dtC.z*dtC.x+D2K.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.x+D2KL.az)+Ddc.z*(DDF*dtD.z*dtC.x+D2KL.az);
              veck.y=Dab.z*(DDF*dtA.z*dtC.y+D2IK.cy)+Dcb.z*rbc*(DDF*dtC.z*dtC.y+D2K.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.y+D2KL.bz)+Ddc.z*(DDF*dtD.z*dtC.y+D2KL.bz);
              veck.z=Dab.z*(DDF*dtA.z*dtC.z+D2IK.cz)+Dcb.z*rbc*(DDF*dtC.z*dtC.z+D2K.cz)+fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.z+D2KL.cz)+fd.z+Ddc.z*(DDF*dtD.z*dtC.z+D2KL.cz)-fd.z;

              vecl.x=Dab.z*(DDF*dtA.z*dtD.x+D2IL.cx)+Dcb.z*rbc*(DDF*dtC.z*dtD.x+D2KL.cx)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.x+D2L.az)+Ddc.z*(DDF*dtD.z*dtD.x+D2L.az);
              vecl.y=Dab.z*(DDF*dtA.z*dtD.y+D2IL.cy)+Dcb.z*rbc*(DDF*dtC.z*dtD.y+D2KL.cy)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.y+D2L.bz)+Ddc.z*(DDF*dtD.z*dtD.y+D2L.bz);
              vecl.z=Dab.z*(DDF*dtA.z*dtD.z+D2IL.cz)+Dcb.z*rbc*(DDF*dtC.z*dtD.z+D2KL.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.z+D2L.cz)+Ddc.z*(DDF*dtD.z*dtD.z+D2L.cz)+fd.z;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+3]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+3]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+3]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n+3]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n+3]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n+3]+=vecl.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              //S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
              veci.x=Dab.x*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.x*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.x*(DDF*dtD.x*dtA.x+D2IL.ax)+fa.x;
              veci.y=Dab.x*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.x*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.x*(DDF*dtD.x*dtA.y+D2IL.bx);
              veci.z=Dab.x*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.x*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.x*(DDF*dtD.x*dtA.z+D2IL.cx);

              vecj.x=Dab.x*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.x*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.x*(DDF*dtD.x*dtB.x+D2JL.ax)-fa.x-fc.x-fd.x;
              vecj.y=Dab.x*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.x*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.x*(DDF*dtD.x*dtB.y+D2JL.bx);
              vecj.z=Dab.x*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.x*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.x*(DDF*dtD.x*dtB.z+D2JL.cx);

              veck.x=Dab.x*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.x*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.x*(DDF*dtD.x*dtC.x+D2KL.ax)+fc.x+fd.x-fd.x;
              veck.y=Dab.x*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.x*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.x*(DDF*dtD.x*dtC.y+D2KL.bx);
              veck.z=Dab.x*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.x*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.x*(DDF*dtD.x*dtC.z+D2KL.cx);

              vecl.x=Dab.x*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.x*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.x*(DDF*dtD.x*dtD.x+D2L.ax)+fd.x;
              vecl.y=Dab.x*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.x*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.x*(DDF*dtD.x*dtD.y+D2L.ay);
              vecl.z=Dab.x*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.x*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.x*(DDF*dtD.x*dtD.z+D2L.az);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n]+=vecl.z;

              //S.bx=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
              veci.x=Dab.y*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.y*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                     Dcb.y*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.y*(DDF*dtD.x*dtA.x+D2IL.ax);
              veci.y=Dab.y*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.y*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                     Dcb.y*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.y*(DDF*dtD.x*dtA.y+D2IL.bx)+fa.x;
              veci.z=Dab.y*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.y*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                     Dcb.y*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.y*(DDF*dtD.x*dtA.z+D2IL.cx);

              vecj.x=Dab.y*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.y*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                     Dcb.y*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.y*(DDF*dtD.x*dtB.x+D2JL.ax);
              vecj.y=Dab.y*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.y*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                     Dcb.y*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.y*(DDF*dtD.x*dtB.y+D2JL.bx)-fa.x-fc.x-fd.x;
              vecj.z=Dab.y*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.y*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                     Dcb.y*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.y*(DDF*dtD.x*dtB.z+D2JL.cx);

              veck.x=Dab.y*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.y*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                     Dcb.y*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.y*(DDF*dtD.x*dtC.x+D2KL.ax);
              veck.y=Dab.y*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.y*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                     Dcb.y*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.y*(DDF*dtD.x*dtC.y+D2KL.bx)+fc.x+fd.x-fd.x;
              veck.z=Dab.y*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.y*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                     Dcb.y*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.y*(DDF*dtD.x*dtC.z+D2KL.cx);

              vecl.x=Dab.y*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.y*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                     Dcb.y*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.y*(DDF*dtD.x*dtD.x+D2L.ax);
              vecl.y=Dab.y*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.y*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                     Dcb.y*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.y*(DDF*dtD.x*dtD.y+D2L.ay)+fd.x;
              vecl.z=Dab.y*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.y*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                     Dcb.y*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.y*(DDF*dtD.x*dtD.z+D2L.az);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+1]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+1]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+1]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n+1]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n+1]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n+1]+=vecl.z;

              veci.x=Dab.y*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.y*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.y*(DDF*dtD.y*dtA.x+D2IL.ay);
              veci.y=Dab.y*(DDF*dtA.y*dtA.y+D2I.by)+fa.y+Dcb.y*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.y*(DDF*dtD.y*dtA.y+D2IL.by);
              veci.z=Dab.y*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.y*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.y*(DDF*dtD.y*dtA.z+D2IL.cy);

              vecj.x=Dab.y*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.y*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.y*(DDF*dtD.y*dtB.x+D2JL.ay);
              vecj.y=Dab.y*(DDF*dtA.y*dtB.y+D2IJ.by)-fa.y+Dcb.y*rbc*(DDF*dtC.y*dtB.y+D2JK.by)-fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.y+D2JL.by)-fd.y+Ddc.y*(DDF*dtD.y*dtB.y+D2JL.by);
              vecj.z=Dab.y*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.y*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.y*(DDF*dtD.y*dtB.z+D2JL.cy);

              veck.x=Dab.y*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.y*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.y*(DDF*dtD.y*dtC.x+D2KL.ay);
              veck.y=Dab.y*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.y*rbc*(DDF*dtC.y*dtC.y+D2K.by)+fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+fd.y+Ddc.y*(DDF*dtD.y*dtC.y+D2KL.by)-fd.y;
              veck.z=Dab.y*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.y*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.y*(DDF*dtD.y*dtC.z+D2KL.cy);

              vecl.x=Dab.y*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.y*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.y*(DDF*dtD.y*dtD.x+D2L.ay);
              vecl.y=Dab.y*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.y*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.y*(DDF*dtD.y*dtD.y+D2L.by)+fd.y;
              vecl.z=Dab.y*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.y*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.y*(DDF*dtD.y*dtD.z+D2L.bz);

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+2]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+2]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+2]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n+2]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n+2]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n+2]+=vecl.z;


              veci.x=Dab.z*(DDF*dtA.z*dtA.x+D2I.az)+Dcb.z*rbc*(DDF*dtC.z*dtA.x+D2IK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.x+D2IL.az)+Ddc.z*(DDF*dtD.z*dtA.x+D2IL.az);
              veci.y=Dab.z*(DDF*dtA.z*dtA.y+D2I.bz)+Dcb.z*rbc*(DDF*dtC.z*dtA.y+D2IK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.y+D2IL.bz)+Ddc.z*(DDF*dtD.z*dtA.y+D2IL.bz);
              veci.z=Dab.z*(DDF*dtA.z*dtA.z+D2I.cz)+fa.z+Dcb.z*rbc*(DDF*dtC.z*dtA.z+D2IK.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.z+D2IL.cz)+Ddc.z*(DDF*dtD.z*dtA.z+D2IL.cz);

              vecj.x=Dab.z*(DDF*dtA.z*dtB.x+D2IJ.cx)+Dcb.z*rbc*(DDF*dtC.z*dtB.x+D2JK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.x+D2JL.az)+Ddc.z*(DDF*dtD.z*dtB.x+D2JL.az);
              vecj.y=Dab.z*(DDF*dtA.z*dtB.y+D2IJ.cy)+Dcb.z*rbc*(DDF*dtC.z*dtB.y+D2JK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.y+D2JL.bz)+Ddc.z*(DDF*dtD.z*dtB.y+D2JL.bz);
              vecj.z=Dab.z*(DDF*dtA.z*dtB.z+D2IJ.cz)-fa.z+Dcb.z*rbc*(DDF*dtC.z*dtB.z+D2JK.cz)-fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.z+D2JL.cz)-fd.z+Ddc.z*(DDF*dtD.z*dtB.z+D2JL.cz);

              veck.x=Dab.z*(DDF*dtA.z*dtC.x+D2IK.cx)+Dcb.z*rbc*(DDF*dtC.z*dtC.x+D2K.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.x+D2KL.az)+Ddc.z*(DDF*dtD.z*dtC.x+D2KL.az);
              veck.y=Dab.z*(DDF*dtA.z*dtC.y+D2IK.cy)+Dcb.z*rbc*(DDF*dtC.z*dtC.y+D2K.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.y+D2KL.bz)+Ddc.z*(DDF*dtD.z*dtC.y+D2KL.bz);
              veck.z=Dab.z*(DDF*dtA.z*dtC.z+D2IK.cz)+Dcb.z*rbc*(DDF*dtC.z*dtC.z+D2K.cz)+fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.z+D2KL.cz)+fd.z+Ddc.z*(DDF*dtD.z*dtC.z+D2KL.cz)-fd.z;

              vecl.x=Dab.z*(DDF*dtA.z*dtD.x+D2IL.cx)+Dcb.z*rbc*(DDF*dtC.z*dtD.x+D2KL.cx)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.x+D2L.az)+Ddc.z*(DDF*dtD.z*dtD.x+D2L.az);
              vecl.y=Dab.z*(DDF*dtA.z*dtD.y+D2IL.cy)+Dcb.z*rbc*(DDF*dtC.z*dtD.y+D2KL.cy)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.y+D2L.bz)+Ddc.z*(DDF*dtD.z*dtD.y+D2L.bz);
              vecl.z=Dab.z*(DDF*dtA.z*dtD.z+D2IL.cz)+Dcb.z*rbc*(DDF*dtC.z*dtD.z+D2KL.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.z+D2L.cz)+Ddc.z*(DDF*dtD.z*dtD.z+D2L.cz)+fd.z;

              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=veci.x;
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=veci.y;
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=veci.z;
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]+=vecj.x;
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]+=vecj.y;
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]+=vecj.z;
              if(index_k.x>=0) HessianMatrix.element[index_k.x][n+3]+=veck.x;
              if(index_k.y>=0) HessianMatrix.element[index_k.y][n+3]+=veck.y;
              if(index_k.z>=0) HessianMatrix.element[index_k.z][n+3]+=veck.z;
              if(index_l.x>=0) HessianMatrix.element[index_l.x][n+3]+=vecl.x;
              if(index_l.y>=0) HessianMatrix.element[index_l.y][n+3]+=vecl.y;
              if(index_l.z>=0) HessianMatrix.element[index_l.z][n+3]+=vecl.z;
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


void HessianTorsionStrainStrain(REAL_MATRIX HessianMatrix,VECTOR vec_u,VECTOR vec_v,VECTOR vec_w,REAL u,REAL v,REAL w,REAL DF,REAL DDF,REAL_MATRIX3x3 S,REAL CosPhi)
{
  int n;
  REAL_MATRIX3x3 dme,dne,dmn,T;
  REAL m2,n2,dot_mn;
  VECTOR vec_m,vec_n;

  n=NumberOfCoordinatesMinimizationVariables;

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

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
          HessianMatrix.element[n][n]+=
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
          HessianMatrix.element[n][n]+=
             2.0*(DDF*T.ax*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ax*T.by);
          HessianMatrix.element[n][n]+=
             2.0*(DDF*T.ax*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.cz);
          HessianMatrix.element[n][n]+=
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
          HessianMatrix.element[n][n]+=
             2.0*(DDF*T.by*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.by*T.cz);
          HessianMatrix.element[n][n]+=
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
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          HessianMatrix.element[n][n]+=
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
          HessianMatrix.element[n][n+1]+=
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
          HessianMatrix.element[n][n+2]+=
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

          HessianMatrix.element[n+1][n+1]+=
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
          HessianMatrix.element[n+1][n+2]+=
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
          HessianMatrix.element[n+2][n+2]+=
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
          break;
        case REGULAR:
          HessianMatrix.element[n][n]+=
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
          HessianMatrix.element[n][n+1]+=
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
          HessianMatrix.element[n][n+2]+=
             DDF*T.ax*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.ax*dme.az/SQR(m2)+dne.ax*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.az+S.az;
          HessianMatrix.element[n][n+3]+=
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
          HessianMatrix.element[n][n+4]+=
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
          HessianMatrix.element[n][n+5]+=
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

          HessianMatrix.element[n+1][n+1]+=
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
          HessianMatrix.element[n+1][n+2]+=
             DDF*T.ay*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.ay*dme.az/SQR(m2)+dne.ay*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.z+vec_u.x*vec_u.y*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.z+vec_u.x*vec_u.y*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.az+0.5*S.bz;
          HessianMatrix.element[n+1][n+3]+=
             DDF*T.ay*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.ay*dme.by/SQR(m2)+dne.ay*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ay*T.by+S.ay;
          HessianMatrix.element[n+1][n+4]+=
             DDF*T.ay*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.bz)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.ay*dme.bz/SQR(m2)+dne.ay*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.z+vec_u.x*vec_u.y*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.z+vec_u.x*vec_u.y*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.bz+0.5*S.az;
          HessianMatrix.element[n+1][n+5]+=
             DDF*T.ay*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.ay*dme.cz/SQR(m2)+dne.ay*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.cz;

          HessianMatrix.element[n+2][n+2]+=DDF*T.az*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.az*dme.az/SQR(m2)+dne.az*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.az+0.5*(S.ax+S.cz);
          HessianMatrix.element[n+2][n+3]+=
             DDF*T.az*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.az*dme.by/SQR(m2)+dne.az*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.az*T.by;
          HessianMatrix.element[n+2][n+4]+=
             DDF*T.az*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.bz)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.az*dme.bz/SQR(m2)+dne.az*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.z+vec_u.x*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.z+vec_u.x*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.bz+0.5*S.ay;
          HessianMatrix.element[n+2][n+5]+=
             DDF*T.az*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.az*dme.cz/SQR(m2)+dne.az*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.cz+S.az;

          HessianMatrix.element[n+3][n+3]+=
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
          HessianMatrix.element[n+3][n+4]+=
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
          HessianMatrix.element[n+3][n+5]+=
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

          HessianMatrix.element[n+4][n+4]+=
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
          HessianMatrix.element[n+4][n+5]+=
             DDF*T.bz*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.bz*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.cz)+2.0*dot_mn*
             ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.bz*dme.cz/SQR(m2)+dne.bz*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.bz*T.cz+S.bz;

          HessianMatrix.element[n+5][n+5]+=
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
          break;
        case REGULAR_UPPER_TRIANGLE:
          HessianMatrix.element[n][n]+=
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
          HessianMatrix.element[n][n+1]+=
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
          HessianMatrix.element[n][n+2]+=
             DDF*T.ax*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.ax*dme.az/SQR(m2)+dne.ax*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.az+S.az;
          HessianMatrix.element[n][n+3]+=
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
          HessianMatrix.element[n][n+4]+=
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
          HessianMatrix.element[n][n+5]+=
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

          HessianMatrix.element[n+1][n+1]+=
             DDF*T.ay*T.ay+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.ay/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.ay)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
             +(dme.ay*dme.ay/SQR(m2)+dne.ay*dne.ay/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ay*T.ay+S.by;
          HessianMatrix.element[n+1][n+2]+=
             DDF*T.ay*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.ay*dme.az/SQR(m2)+dne.ay*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.z+vec_u.x*vec_u.y*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.z+vec_u.x*vec_u.y*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.az+S.bz;
          HessianMatrix.element[n+1][n+3]+=
             DDF*T.ay*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.ay*dme.by/SQR(m2)+dne.ay*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ay*T.by;
          HessianMatrix.element[n+1][n+4]+=
             DDF*T.ay*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.bz)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.ay*dme.bz/SQR(m2)+dne.ay*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.z+vec_u.x*vec_u.y*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.z+vec_u.x*vec_u.y*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.bz;
          HessianMatrix.element[n+1][n+5]+=
             DDF*T.ay*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.ay*dme.cz/SQR(m2)+dne.ay*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.cz;

          HessianMatrix.element[n+2][n+2]+=DDF*T.az*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.az*dme.az/SQR(m2)+dne.az*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.az+S.cz;
          HessianMatrix.element[n+2][n+3]+=
             DDF*T.az*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.az*dme.by/SQR(m2)+dne.az*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.az*T.by;
          HessianMatrix.element[n+2][n+4]+=
             DDF*T.az*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.bz)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.az*dme.bz/SQR(m2)+dne.az*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.z+vec_u.x*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.z+vec_u.x*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.bz;
          HessianMatrix.element[n+2][n+5]+=
             DDF*T.az*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.az*dme.cz/SQR(m2)+dne.az*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.cz;

          HessianMatrix.element[n+3][n+3]+=
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
          HessianMatrix.element[n+3][n+4]+=
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
          HessianMatrix.element[n+3][n+5]+=
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

          HessianMatrix.element[n+4][n+4]+=
             DDF*T.bz*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.bz*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.bz)+2.0*dot_mn*
             ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.bz*dme.bz/SQR(m2)+dne.bz*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.bz*T.bz+S.cz;
          HessianMatrix.element[n+4][n+5]+=
             DDF*T.bz*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.bz*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.cz)+2.0*dot_mn*
             ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.bz*dme.cz/SQR(m2)+dne.bz*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.bz*T.cz;

          HessianMatrix.element[n+5][n+5]+=
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
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n]+=
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
              HessianMatrix.element[n][n+1]+=
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
              HessianMatrix.element[n][n+2]+=
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
              HessianMatrix.element[n][n+3]+=
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

              HessianMatrix.element[n+1][n+1]+=
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
              HessianMatrix.element[n+1][n+2]+=
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
              HessianMatrix.element[n+1][n+3]+=
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

              HessianMatrix.element[n+2][n+2]+=
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
              HessianMatrix.element[n+2][n+3]+=
                 DDF*T.bz*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.bz*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.cz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.bz*dme.cz/SQR(m2)+dne.bz*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.bz*T.cz+S.bz;

              HessianMatrix.element[n+3][n+3]+=
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
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=
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
              HessianMatrix.element[n][n+1]+=
                 DDF*T.ax*T.az+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.az/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.az)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
                 +(dme.ax*dme.az/SQR(m2)+dne.ax*dne.az/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.az+S.az;
              HessianMatrix.element[n][n+2]+=
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
              HessianMatrix.element[n][n+3]+=
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

              HessianMatrix.element[n+1][n+1]+=DDF*T.az*T.az+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.az/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.az)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
                 +(dme.az*dme.az/SQR(m2)+dne.az*dne.az/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.az*T.az+0.5*(S.ax+S.cz);
              HessianMatrix.element[n+1][n+2]+=
                 DDF*T.az*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.az*dme.by/SQR(m2)+dne.az*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.az*T.by;
              HessianMatrix.element[n+1][n+3]+=
                 DDF*T.az*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.az*dme.cz/SQR(m2)+dne.az*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.az*T.cz+S.az;

              HessianMatrix.element[n+2][n+2]+=
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
              HessianMatrix.element[n+2][n+3]+=
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

              HessianMatrix.element[n+3][n+3]+=
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
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=
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
              HessianMatrix.element[n][n+1]+=
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
              HessianMatrix.element[n][n+2]+=
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
              HessianMatrix.element[n][n+3]+=
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

              HessianMatrix.element[n+1][n+1]+=
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
              HessianMatrix.element[n+1][n+2]+=
                 DDF*T.ay*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.ay*dme.by/SQR(m2)+dne.ay*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ay*T.by+S.ay;
              HessianMatrix.element[n+1][n+3]+=
                 DDF*T.ay*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.ay*dme.cz/SQR(m2)+dne.ay*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ay*T.cz;

              HessianMatrix.element[n+2][n+2]+=
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
              HessianMatrix.element[n+2][n+3]+=
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

              HessianMatrix.element[n+3][n+3]+=
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
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n]+=
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
              HessianMatrix.element[n][n+1]+=
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
              HessianMatrix.element[n][n+2]+=
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
              HessianMatrix.element[n][n+3]+=
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

              HessianMatrix.element[n+1][n+1]+=
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
              HessianMatrix.element[n+1][n+2]+=
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
              HessianMatrix.element[n+1][n+3]+=
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

              HessianMatrix.element[n+2][n+2]+=
                 DDF*T.bz*T.bz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.bz*dmn.bz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.bz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
                 -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
                 +(dme.bz*dme.bz/SQR(m2)+dne.bz*dne.bz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.bz*T.bz+S.cz;
              HessianMatrix.element[n+2][n+3]+=
                 DDF*T.bz*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.bz*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.cz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.bz*dme.cz/SQR(m2)+dne.bz*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.bz*T.cz;

              HessianMatrix.element[n+3][n+3]+=
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
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=
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
              HessianMatrix.element[n][n+1]+=
                 DDF*T.ax*T.az+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.az/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.az)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
                 +(dme.ax*dme.az/SQR(m2)+dne.ax*dne.az/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.az+S.az;
              HessianMatrix.element[n][n+2]+=
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
              HessianMatrix.element[n][n+3]+=
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

              HessianMatrix.element[n+1][n+1]+=DDF*T.az*T.az+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.az/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.az)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
                 +(dme.az*dme.az/SQR(m2)+dne.az*dne.az/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.az*T.az+S.cz;
              HessianMatrix.element[n+1][n+2]+=
                 DDF*T.az*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.az*dme.by/SQR(m2)+dne.az*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.az*T.by;
              HessianMatrix.element[n+1][n+3]+=
                 DDF*T.az*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.az*dme.cz/SQR(m2)+dne.az*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.az*T.cz;

              HessianMatrix.element[n+2][n+2]+=
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
              HessianMatrix.element[n+2][n+3]+=
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

              HessianMatrix.element[n+3][n+3]+=
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
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=
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
              HessianMatrix.element[n][n+1]+=
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
              HessianMatrix.element[n][n+2]+=
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
              HessianMatrix.element[n][n+3]+=
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

              HessianMatrix.element[n+1][n+1]+=
                 DDF*T.ay*T.ay+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.ay/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.ay)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
                 +(dme.ay*dme.ay/SQR(m2)+dne.ay*dne.ay/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ay*T.ay+S.by;
              HessianMatrix.element[n+1][n+2]+=
                 DDF*T.ay*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.ay*dme.by/SQR(m2)+dne.ay*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ay*T.by;
              HessianMatrix.element[n+1][n+3]+=
                 DDF*T.ay*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.ay*dme.cz/SQR(m2)+dne.ay*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ay*T.cz;

              HessianMatrix.element[n+2][n+2]+=
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
              HessianMatrix.element[n+2][n+3]+=
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

              HessianMatrix.element[n+3][n+3]+=
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



void CalculateAdsorbateBondHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,Type,NumberOfBonds,A,B;
  REAL U,DF,DDF,r,rr,temp,temp2,exp_term;
  REAL *parms;
  POINT posA,posB,comA,comB;
  VECTOR dr;
  int grpA,grpB,RigidI,RigidJ;
  INT_VECTOR3 index_i,index_j;
  INT_VECTOR3 index_i2,index_j2;
  int index1,index2;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfBonds=Components[Type].NumberOfBonds;
    for(i=0;i<NumberOfBonds;i++)
    {
      A=MIN2(Components[Type].Bonds[i].A,Components[Type].Bonds[i].B);
      B=MAX2(Components[Type].Bonds[i].A,Components[Type].Bonds[i].B);

      comA=posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      comB=posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

      index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
      index_i2=UNDEFINED_INT_VECTOR3;
      index_j2=UNDEFINED_INT_VECTOR3;
      index1=-1;
      index2=-1;

      grpA=Components[Type].group[A];
      RigidI=Components[Type].Groups[grpA].Rigid;
      if(RigidI)
      {
        comA=Adsorbates[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
        index_i=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndex;
        index_i2=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
        index1=Adsorbates[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      }

      grpB=Components[Type].group[B];
      RigidJ=Components[Type].Groups[grpB].Rigid;
      if(RigidJ)
      {
        comB=Adsorbates[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
        index_j=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndex;
        index_j2=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
        index2=Adsorbates[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      }

      parms=(REAL*)&Components[Type].BondArguments[i];

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      switch(Components[Type].BondType[i])
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
          temp=CUBE(SQR(parms[1])/rr);
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
        default:
          U=DF=DDF=0.0;
          fprintf(stderr, "Undefined Bond potential (Bond Hessian adsorbate)\n");
          exit(0);
          break;
      }

      // add contribution to the energy
      *Energy+=U;

      //if((index_i<0)&&(index_j<0)) continue;

      // add contribution to the first derivatives

      if(ComputeGradient)
      {
        if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
        if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
        if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

        if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
        if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
        if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

        if(RigidI)
        {
          GradientStrainI(Gradient,DF,dr,posA,comA);

          if(index_i2.x>=0) Gradient[index_i2.x]+=DF*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
          if(index_i2.y>=0) Gradient[index_i2.y]+=DF*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
          if(index_i2.z>=0) Gradient[index_i2.z]+=DF*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
        }

        if(RigidJ)
        {
          GradientStrainJ(Gradient,DF,dr,posB,comB);

          if(index_j2.x>=0) Gradient[index_j2.x]-=DF*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
          if(index_j2.y>=0) Gradient[index_j2.y]-=DF*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
          if(index_j2.z>=0) Gradient[index_j2.z]-=DF*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
        }

        GradientStrain(Gradient,DF,dr);
      }

      // add contribution to the strain derivative tensor
      StrainDerivativeTensor->ax+=dr.x*DF*dr.x;
      StrainDerivativeTensor->bx+=dr.y*DF*dr.x;
      StrainDerivativeTensor->cx+=dr.z*DF*dr.x;

      StrainDerivativeTensor->ay+=dr.x*DF*dr.y;
      StrainDerivativeTensor->by+=dr.y*DF*dr.y;
      StrainDerivativeTensor->cy+=dr.z*DF*dr.y;

      StrainDerivativeTensor->az+=dr.x*DF*dr.z;
      StrainDerivativeTensor->bz+=dr.y*DF*dr.z;
      StrainDerivativeTensor->cz+=dr.z*DF*dr.z;

      if(ComputeHessian)
      {
        HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);
        HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);

        HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB,RigidI,RigidJ);
        HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,DF,DDF,posA,comA,posB,comB,dr);

        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
        HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB);
      }
    }
  }
}

void CalculateCationBondHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,Type,NumberOfBonds,A,B;
  REAL U,DF,DDF,r,rr,temp,temp2,exp_term;
  REAL *parms;
  POINT posA,posB,comA,comB;
  VECTOR dr;
  int grpA,grpB,RigidI,RigidJ;
  INT_VECTOR3 index_i,index_j;
  INT_VECTOR3 index_i2,index_j2;
  int index1,index2;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    NumberOfBonds=Components[Type].NumberOfBonds;
    for(i=0;i<NumberOfBonds;i++)
    {
      A=MIN2(Components[Type].Bonds[i].A,Components[Type].Bonds[i].B);
      B=MAX2(Components[Type].Bonds[i].A,Components[Type].Bonds[i].B);

      comA=posA=Cations[CurrentSystem][m].Atoms[A].Position;
      comB=posB=Cations[CurrentSystem][m].Atoms[B].Position;

      index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Cations[CurrentSystem][m].Atoms[B].HessianIndex;
      index_i2=UNDEFINED_INT_VECTOR3;
      index_j2=UNDEFINED_INT_VECTOR3;
      index1=-1;
      index2=-1;

      grpA=Components[Type].group[A];
      RigidI=Components[Type].Groups[grpA].Rigid;
      if(RigidI)
      {
        comA=Cations[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
        index_i=Cations[CurrentSystem][m].Groups[grpA].HessianIndex;
        index_i2=Cations[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
        index1=Cations[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      }

      grpB=Components[Type].group[B];
      RigidJ=Components[Type].Groups[grpB].Rigid;
      if(RigidJ)
      {
        comB=Cations[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
        index_j=Cations[CurrentSystem][m].Groups[grpB].HessianIndex;
        index_j2=Cations[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
        index2=Cations[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      }

      parms=(REAL*)&Components[Type].BondArguments[i];

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      switch(Components[Type].BondType[i])
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
        default:
          U=DF=DDF=0.0;
          fprintf(stderr, "Undefined Bond potential (Bond Hessian adsorbate)\n");
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

        if(RigidI)
        {
          GradientStrainI(Gradient,DF,dr,posA,comA);

          if(index_i2.x>=0) Gradient[index_i2.x]+=DF*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
          if(index_i2.y>=0) Gradient[index_i2.y]+=DF*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
          if(index_i2.z>=0) Gradient[index_i2.z]+=DF*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
        }

        if(RigidJ)
        {
          GradientStrainJ(Gradient,DF,dr,posB,comB);

          if(index_j2.x>=0) Gradient[index_j2.x]-=DF*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
          if(index_j2.y>=0) Gradient[index_j2.y]-=DF*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
          if(index_j2.z>=0) Gradient[index_j2.z]-=DF*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
        }

        GradientStrain(Gradient,DF,dr);
      }


      // add contribution to the strain derivative tensor
      StrainDerivativeTensor->ax+=dr.x*DF*dr.x;
      StrainDerivativeTensor->bx+=dr.y*DF*dr.x;
      StrainDerivativeTensor->cx+=dr.z*DF*dr.x;

      StrainDerivativeTensor->ay+=dr.x*DF*dr.y;
      StrainDerivativeTensor->by+=dr.y*DF*dr.y;
      StrainDerivativeTensor->cy+=dr.z*DF*dr.y;

      StrainDerivativeTensor->az+=dr.x*DF*dr.z;
      StrainDerivativeTensor->bz+=dr.y*DF*dr.z;
      StrainDerivativeTensor->cz+=dr.z*DF*dr.z;

      if(ComputeHessian)
      {
        HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);
        HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);

        HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB,RigidI,RigidJ);
        HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,DF,DDF,posA,comA,posB,comB,dr);

        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
        HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB);
      }
    }
  }
}

// Calculate Urey-Bradley Hessian
void CalculateAdsorbateUreyBradleyHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
}


// Calculate Urey-Bradley Hessian
void CalculateCationUreyBradleyHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
}

// The first and second derivative of the potential with respect to the angle theta are required:
// DF=D[U[theta],theta]/sin(x)
// DDF=D[D[U[theta]]/sin(theta),theta]/sin(x)

void CalculateAdsorbateBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,Type,NumberOfBends,A,B,C;
  REAL *parms,U;
  REAL CosTheta,Theta,SinTheta,temp,temp2;
  REAL rab,rbc,rac;
  POINT posA,posB,posC;
  POINT comA,comB,comC;
  VECTOR Rab,Rbc,Rac;
  int grpA,grpB,grpC;
  INT_VECTOR3 index_i,index_j,index_k;
  INT_VECTOR3 index_i2,index_j2,index_k2;
  int index1,index2,index3;
  int RigidI,RigidJ,RigidK;
  REAL DTDX,DF,DDF;
  REAL_MATRIX3x3 D2I,D2K,D2IK,S;
  VECTOR dtA,dtB,dtC;
  VECTOR vec_u,vec_v;
  REAL u,v;

  S.ax=S.bx=S.cx=0.0;
  S.ay=S.by=S.cy=0.0;
  S.az=S.bz=S.cz=0.0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfBends=Components[Type].NumberOfBends;

    for(i=0;i<NumberOfBends;i++)
    {
      A=Components[Type].Bends[i].A;
      B=Components[Type].Bends[i].B;
      C=Components[Type].Bends[i].C;
      parms=(REAL*)&Components[Type].BendArguments[i];

      index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
      index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;
      index_i2=UNDEFINED_INT_VECTOR3;
      index_j2=UNDEFINED_INT_VECTOR3;
      index_k2=UNDEFINED_INT_VECTOR3;
      index1=-1;
      index2=-1;
      index3=-1;

      posA=comA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=comB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=comC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

      grpA=Components[Type].group[A];
      RigidI=Components[Type].Groups[grpA].Rigid;
      index1=Adsorbates[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      if(RigidI)
      {
        comA=Adsorbates[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
        index_i=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndex;
        index_i2=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
      }

      grpB=Components[Type].group[B];
      RigidJ=Components[Type].Groups[grpB].Rigid;
      index2=Adsorbates[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      if(RigidJ)
      {
        comB=Adsorbates[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
        index_j=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndex;
        index_j2=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
      }

      grpC=Components[Type].group[C];
      RigidK=Components[Type].Groups[grpC].Rigid;
      index3=Adsorbates[CurrentSystem][m].Atoms[C].HessianAtomIndex;
      if(RigidK)
      {
        comC=Adsorbates[CurrentSystem][m].Groups[grpC].CenterOfMassPosition;
        index_k=Adsorbates[CurrentSystem][m].Groups[grpC].HessianIndex;
        index_k2=Adsorbates[CurrentSystem][m].Groups[grpC].HessianIndexOrientation;
      }

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      vec_u=Rab;
      rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
      u=rab;
      Rab.x/=rab;
      Rab.y/=rab;
      Rab.z/=rab;

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      vec_v=Rbc;
      rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
      v=rbc;
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
      SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));
      DTDX=-1.0/SinTheta;

      switch(Components[Type].BendType[i])
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

      if(ComputeGradient)
      {
        // add contribution to the first derivatives
        if(index_i.x>=0) Gradient[index_i.x]+=DF*dtA.x;
        if(index_i.y>=0) Gradient[index_i.y]+=DF*dtA.y;
        if(index_i.z>=0) Gradient[index_i.z]+=DF*dtA.z;

        if(index_j.x>=0) Gradient[index_j.x]+=DF*dtB.x;
        if(index_j.y>=0) Gradient[index_j.y]+=DF*dtB.y;
        if(index_j.z>=0) Gradient[index_j.z]+=DF*dtB.z;

        if(index_k.x>=0) Gradient[index_k.x]+=DF*dtC.x;
        if(index_k.y>=0) Gradient[index_k.y]+=DF*dtC.y;
        if(index_k.z>=0) Gradient[index_k.z]+=DF*dtC.z;

        if(RigidI)
        {
          GradientStrainI(Gradient,DF,dtA,posA,comA);

          if(index_i2.x>=0) Gradient[index_i2.x]+=DF*(dtA.x*DVecX[index1].x+dtA.y*DVecX[index1].y+dtA.z*DVecX[index1].z);
          if(index_i2.y>=0) Gradient[index_i2.y]+=DF*(dtA.x*DVecY[index1].x+dtA.y*DVecY[index1].y+dtA.z*DVecY[index1].z);
          if(index_i2.z>=0) Gradient[index_i2.z]+=DF*(dtA.x*DVecZ[index1].x+dtA.y*DVecZ[index1].y+dtA.z*DVecZ[index1].z);
        }

        if(RigidJ)
        {
          GradientStrainI(Gradient,DF,dtB,posB,comB);

          if(index_j2.x>=0) Gradient[index_j2.x]+=DF*(dtB.x*DVecX[index2].x+dtB.y*DVecX[index2].y+dtB.z*DVecX[index2].z);
          if(index_j2.y>=0) Gradient[index_j2.y]+=DF*(dtB.x*DVecY[index2].x+dtB.y*DVecY[index2].y+dtB.z*DVecY[index2].z);
          if(index_j2.z>=0) Gradient[index_j2.z]+=DF*(dtB.x*DVecZ[index2].x+dtB.y*DVecZ[index2].y+dtB.z*DVecZ[index2].z);
        }

        if(RigidK)
        {
          GradientStrainI(Gradient,DF,dtC,posC,comC);

          if(index_k2.x>=0) Gradient[index_k2.x]+=DF*(dtC.x*DVecX[index3].x+dtC.y*DVecX[index3].y+dtC.z*DVecX[index3].z);
          if(index_k2.y>=0) Gradient[index_k2.y]+=DF*(dtC.x*DVecY[index3].x+dtC.y*DVecY[index3].y+dtC.z*DVecY[index3].z);
          if(index_k2.z>=0) Gradient[index_k2.z]+=DF*(dtC.x*DVecZ[index3].x+dtC.y*DVecZ[index3].y+dtC.z*DVecZ[index3].z);
        }

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

void CalculateCationBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,Type,NumberOfBends,A,B,C;
  REAL *parms,U;
  REAL CosTheta,Theta,SinTheta,temp,temp2;
  REAL rab,rbc,rac;
  POINT posA,posB,posC;
  POINT comA,comB,comC;
  VECTOR Rab,Rbc,Rac;
  int grpA,grpB,grpC;
  INT_VECTOR3 index_i,index_j,index_k;
  INT_VECTOR3 index_i2,index_j2,index_k2;
  int index1,index2,index3;
  int RigidI,RigidJ,RigidK;
  REAL DTDX,DF,DDF;
  REAL_MATRIX3x3 D2I,D2K,D2IK,S;
  VECTOR dtA,dtB,dtC;
  VECTOR vec_u,vec_v;
  REAL u,v;

  S.ax=S.bx=S.cx=0.0;
  S.ay=S.by=S.cy=0.0;
  S.az=S.bz=S.cz=0.0;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    NumberOfBends=Components[Type].NumberOfBends;

    for(i=0;i<NumberOfBends;i++)
    {
      A=Components[Type].Bends[i].A;
      B=Components[Type].Bends[i].B;
      C=Components[Type].Bends[i].C;
      parms=(REAL*)&Components[Type].BendArguments[i];

      index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Cations[CurrentSystem][m].Atoms[B].HessianIndex;
      index_k=Cations[CurrentSystem][m].Atoms[C].HessianIndex;
      index_i2=UNDEFINED_INT_VECTOR3;
      index_j2=UNDEFINED_INT_VECTOR3;
      index_k2=UNDEFINED_INT_VECTOR3;
      index1=-1;
      index2=-1;
      index3=-1;

      posA=comA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=comB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=comC=Cations[CurrentSystem][m].Atoms[C].Position;

      grpA=Components[Type].group[A];
      RigidI=Components[Type].Groups[grpA].Rigid;
      index1=Cations[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      if(RigidI)
      {
        comA=Cations[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
        index_i=Cations[CurrentSystem][m].Groups[grpA].HessianIndex;
        index_i2=Cations[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
      }

      grpB=Components[Type].group[B];
      RigidJ=Components[Type].Groups[grpB].Rigid;
      index2=Cations[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      if(RigidJ)
      {
        comB=Cations[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
        index_j=Cations[CurrentSystem][m].Groups[grpB].HessianIndex;
        index_j2=Cations[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
      }

      grpC=Components[Type].group[C];
      RigidK=Components[Type].Groups[grpC].Rigid;
      index3=Cations[CurrentSystem][m].Atoms[C].HessianAtomIndex;
      if(RigidK)
      {
        comC=Cations[CurrentSystem][m].Groups[grpC].CenterOfMassPosition;
        index_k=Cations[CurrentSystem][m].Groups[grpC].HessianIndex;
        index_k2=Cations[CurrentSystem][m].Groups[grpC].HessianIndexOrientation;
      }

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      vec_u=Rab;
      rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
      u=rab;
      Rab.x/=rab;
      Rab.y/=rab;
      Rab.z/=rab;

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      vec_v=Rbc;
      rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
      v=rbc;
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
      SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));
      DTDX=-1.0/SinTheta;

      switch(Components[Type].BendType[i])
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

      if(ComputeGradient)
      {
        // add contribution to the first derivatives
        if(index_i.x>=0) Gradient[index_i.x]+=DF*dtA.x;
        if(index_i.y>=0) Gradient[index_i.y]+=DF*dtA.y;
        if(index_i.z>=0) Gradient[index_i.z]+=DF*dtA.z;

        if(index_j.x>=0) Gradient[index_j.x]+=DF*dtB.x;
        if(index_j.y>=0) Gradient[index_j.y]+=DF*dtB.y;
        if(index_j.z>=0) Gradient[index_j.z]+=DF*dtB.z;

        if(index_k.x>=0) Gradient[index_k.x]+=DF*dtC.x;
        if(index_k.y>=0) Gradient[index_k.y]+=DF*dtC.y;
        if(index_k.z>=0) Gradient[index_k.z]+=DF*dtC.z;

        if(RigidI)
        {
          GradientStrainI(Gradient,DF,dtA,posA,comA);

          if(index_i2.x>=0) Gradient[index_i2.x]+=DF*(dtA.x*DVecX[index1].x+dtA.y*DVecX[index1].y+dtA.z*DVecX[index1].z);
          if(index_i2.y>=0) Gradient[index_i2.y]+=DF*(dtA.x*DVecY[index1].x+dtA.y*DVecY[index1].y+dtA.z*DVecY[index1].z);
          if(index_i2.z>=0) Gradient[index_i2.z]+=DF*(dtA.x*DVecZ[index1].x+dtA.y*DVecZ[index1].y+dtA.z*DVecZ[index1].z);
        }

        if(RigidJ)
        {
          GradientStrainI(Gradient,DF,dtB,posB,comB);

          if(index_j2.x>=0) Gradient[index_j2.x]+=DF*(dtB.x*DVecX[index2].x+dtB.y*DVecX[index2].y+dtB.z*DVecX[index2].z);
          if(index_j2.y>=0) Gradient[index_j2.y]+=DF*(dtB.x*DVecY[index2].x+dtB.y*DVecY[index2].y+dtB.z*DVecY[index2].z);
          if(index_j2.z>=0) Gradient[index_j2.z]+=DF*(dtB.x*DVecZ[index2].x+dtB.y*DVecZ[index2].y+dtB.z*DVecZ[index2].z);
        }

        if(RigidK)
        {
          GradientStrainI(Gradient,DF,dtC,posC,comC);

          if(index_k2.x>=0) Gradient[index_k2.x]+=DF*(dtC.x*DVecX[index3].x+dtC.y*DVecX[index3].y+dtC.z*DVecX[index3].z);
          if(index_k2.y>=0) Gradient[index_k2.y]+=DF*(dtC.x*DVecY[index3].x+dtC.y*DVecY[index3].y+dtC.z*DVecY[index3].z);
          if(index_k2.z>=0) Gradient[index_k2.z]+=DF*(dtC.x*DVecZ[index3].x+dtC.y*DVecZ[index3].y+dtC.z*DVecZ[index3].z);
        }

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

void CalculateAdsorbateTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,A,B,C,D,Type;
  POINT posA,posB,posC,posD;
  POINT comA,comB,comC,comD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  REAL ShiftedCosPhi,ShiftedCosPhi2,ShiftedSinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL *parms;
  int RigidI,RigidJ,RigidK,RigidL;
  int grpA,grpB,grpC,grpD;
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  INT_VECTOR3 index_i2,index_j2,index_k2,index_l2;
  int index1,index2,index3,index4;
  REAL_MATRIX3x3 D2I,D2L,D2IL,D2IJ,D2IK,D2JL;
  REAL_MATRIX3x3 D2KL,D2J,D2K,D2JK,S;
  VECTOR DIL,DDJ,DDK,DEJ,DEK;
  VECTOR fa,fb,fc,fd;
  REAL DF,DDF;
  VECTOR vec_u,vec_v,vec_w;
  REAL u,w,v;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfTorsions;i++)
    {
      A=Components[Type].Torsions[i].A;
      B=Components[Type].Torsions[i].B;
      C=Components[Type].Torsions[i].C;
      D=Components[Type].Torsions[i].D;
      parms=(REAL*)&Components[Type].TorsionArguments[i];

      posA=comA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=comB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=comC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
      posD=comD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

      index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
      index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;
      index_l=Adsorbates[CurrentSystem][m].Atoms[D].HessianIndex;
      index_i2=UNDEFINED_INT_VECTOR3;
      index_j2=UNDEFINED_INT_VECTOR3;
      index_k2=UNDEFINED_INT_VECTOR3;
      index_l2=UNDEFINED_INT_VECTOR3;
      index1=-1;
      index2=-1;
      index3=-1;
      index4=-1;

      grpA=Components[Type].group[A];
      RigidI=Components[Type].Groups[grpA].Rigid;
      index1=Adsorbates[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      if(RigidI)
      {
        comA=Adsorbates[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
        index_i=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndex;
        index_i2=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
      }

      grpB=Components[Type].group[B];
      RigidJ=Components[Type].Groups[grpB].Rigid;
      index2=Adsorbates[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      if(RigidJ)
      {
        comB=Adsorbates[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
        index_j=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndex;
        index_j2=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
      }

      grpC=Components[Type].group[C];
      RigidK=Components[Type].Groups[grpC].Rigid;
      index3=Adsorbates[CurrentSystem][m].Atoms[C].HessianAtomIndex;
      if(RigidK)
      {
        comC=Adsorbates[CurrentSystem][m].Groups[grpC].CenterOfMassPosition;
        index_k=Adsorbates[CurrentSystem][m].Groups[grpC].HessianIndex;
        index_k2=Adsorbates[CurrentSystem][m].Groups[grpC].HessianIndexOrientation;
      }

      grpD=Components[Type].group[D];
      RigidL=Components[Type].Groups[grpD].Rigid;
      index4=Adsorbates[CurrentSystem][m].Atoms[D].HessianAtomIndex;
      if(RigidL)
      {
        comD=Adsorbates[CurrentSystem][m].Groups[grpD].CenterOfMassPosition;
        index_l=Adsorbates[CurrentSystem][m].Groups[grpD].HessianIndex;
        index_l2=Adsorbates[CurrentSystem][m].Groups[grpD].HessianIndexOrientation;
      }

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
          SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
          U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
          DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
          DDF=-(SinPhi*parms[0]*SQR(parms[1])*cos(parms[1]*Phi-parms[2])-parms[0]*parms[1]*sin(parms[1]*Phi-parms[2])*CosPhi)/CUBE(SinPhi);
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
        case FIXED_DIHEDRAL:
          U=DF=DDF=0.0;
          break;
        default:
          fprintf(stderr, "Undefined Torsion potential\n");
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

        if(RigidI)
        {
          GradientStrainI(Gradient,1.0,fa,posA,comA);

          if(index_i2.x>=0) Gradient[index_i2.x]+=(fa.x*DVecX[index1].x+fa.y*DVecX[index1].y+fa.z*DVecX[index1].z);
          if(index_i2.y>=0) Gradient[index_i2.y]+=(fa.x*DVecY[index1].x+fa.y*DVecY[index1].y+fa.z*DVecY[index1].z);
          if(index_i2.z>=0) Gradient[index_i2.z]+=(fa.x*DVecZ[index1].x+fa.y*DVecZ[index1].y+fa.z*DVecZ[index1].z);
        }

        if(RigidJ)
        {
          GradientStrainI(Gradient,1.0,fb,posB,comB);

          if(index_j2.x>=0) Gradient[index_j2.x]+=(fb.x*DVecX[index2].x+fb.y*DVecX[index2].y+fb.z*DVecX[index2].z);
          if(index_j2.y>=0) Gradient[index_j2.y]+=(fb.x*DVecY[index2].x+fb.y*DVecY[index2].y+fb.z*DVecY[index2].z);
          if(index_j2.z>=0) Gradient[index_j2.z]+=(fb.x*DVecZ[index2].x+fb.y*DVecZ[index2].y+fb.z*DVecZ[index2].z);
        }

        if(RigidK)
        {
          GradientStrainI(Gradient,1.0,fc,posC,comC);

          if(index_k2.x>=0) Gradient[index_k2.x]+=(fc.x*DVecX[index3].x+fc.y*DVecX[index3].y+fc.z*DVecX[index3].z);
          if(index_k2.y>=0) Gradient[index_k2.y]+=(fc.x*DVecY[index3].x+fc.y*DVecY[index3].y+fc.z*DVecY[index3].z);
          if(index_k2.z>=0) Gradient[index_k2.z]+=(fc.x*DVecZ[index3].x+fc.y*DVecZ[index3].y+fc.z*DVecZ[index3].z);
        }

        if(RigidL)
        {
          GradientStrainI(Gradient,1.0,fd,posD,comD);

          if(index_l2.x>=0) Gradient[index_l2.x]+=(fd.x*DVecX[index4].x+fd.y*DVecX[index4].y+fd.z*DVecX[index4].z);
          if(index_l2.y>=0) Gradient[index_l2.y]+=(fd.x*DVecY[index4].x+fd.y*DVecY[index4].y+fd.z*DVecY[index4].z);
          if(index_l2.z>=0) Gradient[index_l2.z]+=(fd.x*DVecZ[index4].x+fd.y*DVecZ[index4].y+fd.z*DVecZ[index4].z);
        }

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

void CalculateCationTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,A,B,C,D,Type;
  POINT posA,posB,posC,posD;
  POINT comA,comB,comC,comD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  REAL ShiftedCosPhi,ShiftedCosPhi2,ShiftedSinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL *parms;
  int RigidI,RigidJ,RigidK,RigidL;
  int grpA,grpB,grpC,grpD;
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  INT_VECTOR3 index_i2,index_j2,index_k2,index_l2;
  int index1,index2,index3,index4;
  REAL_MATRIX3x3 D2I,D2L,D2IL,D2IJ,D2IK,D2JL;
  REAL_MATRIX3x3 D2KL,D2J,D2K,D2JK,S;
  VECTOR DIL,DDJ,DDK,DEJ,DEK;
  VECTOR fa,fb,fc,fd;
  REAL DF,DDF;
  VECTOR vec_u,vec_v,vec_w;
  REAL u,w,v;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfTorsions;i++)
    {
      A=Components[Type].Torsions[i].A;
      B=Components[Type].Torsions[i].B;
      C=Components[Type].Torsions[i].C;
      D=Components[Type].Torsions[i].D;
      parms=(REAL*)&Components[Type].TorsionArguments[i];

      posA=comA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=comB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=comC=Cations[CurrentSystem][m].Atoms[C].Position;
      posD=comD=Cations[CurrentSystem][m].Atoms[D].Position;

      index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Cations[CurrentSystem][m].Atoms[B].HessianIndex;
      index_k=Cations[CurrentSystem][m].Atoms[C].HessianIndex;
      index_l=Cations[CurrentSystem][m].Atoms[D].HessianIndex;
      index_i2=UNDEFINED_INT_VECTOR3;
      index_j2=UNDEFINED_INT_VECTOR3;
      index_k2=UNDEFINED_INT_VECTOR3;
      index_l2=UNDEFINED_INT_VECTOR3;
      index1=-1;
      index2=-1;
      index3=-1;
      index4=-1;

      grpA=Components[Type].group[A];
      RigidI=Components[Type].Groups[grpA].Rigid;
      index1=Cations[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      if(RigidI)
      {
        comA=Cations[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
        index_i=Cations[CurrentSystem][m].Groups[grpA].HessianIndex;
        index_i2=Cations[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
      }

      grpB=Components[Type].group[B];
      RigidJ=Components[Type].Groups[grpB].Rigid;
      index2=Cations[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      if(RigidJ)
      {
        comB=Cations[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
        index_j=Cations[CurrentSystem][m].Groups[grpB].HessianIndex;
        index_j2=Cations[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
      }

      grpC=Components[Type].group[C];
      RigidK=Components[Type].Groups[grpC].Rigid;
      index3=Cations[CurrentSystem][m].Atoms[C].HessianAtomIndex;
      if(RigidK)
      {
        comC=Cations[CurrentSystem][m].Groups[grpC].CenterOfMassPosition;
        index_k=Cations[CurrentSystem][m].Groups[grpC].HessianIndex;
        index_k2=Cations[CurrentSystem][m].Groups[grpC].HessianIndexOrientation;
      }

      grpD=Components[Type].group[D];
      RigidL=Components[Type].Groups[grpD].Rigid;
      index4=Cations[CurrentSystem][m].Atoms[D].HessianAtomIndex;
      if(RigidL)
      {
        comD=Cations[CurrentSystem][m].Groups[grpD].CenterOfMassPosition;
        index_l=Cations[CurrentSystem][m].Groups[grpD].HessianIndex;
        index_l2=Cations[CurrentSystem][m].Groups[grpD].HessianIndexOrientation;
      }

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
          SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
          U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
          DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
          DDF=-(SinPhi*parms[0]*SQR(parms[1])*cos(parms[1]*Phi-parms[2])-parms[0]*parms[1]*sin(parms[1]*Phi-parms[2])*CosPhi)/CUBE(SinPhi);
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
        case FIXED_DIHEDRAL:
          U=DF=DDF=0.0;
          break;
        default:
          fprintf(stderr, "Undefined Torsion potential\n");
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

        if(RigidI)
        {
          GradientStrainI(Gradient,1.0,fa,posA,comA);

          if(index_i2.x>=0) Gradient[index_i2.x]+=(fa.x*DVecX[index1].x+fa.y*DVecX[index1].y+fa.z*DVecX[index1].z);
          if(index_i2.y>=0) Gradient[index_i2.y]+=(fa.x*DVecY[index1].x+fa.y*DVecY[index1].y+fa.z*DVecY[index1].z);
          if(index_i2.z>=0) Gradient[index_i2.z]+=(fa.x*DVecZ[index1].x+fa.y*DVecZ[index1].y+fa.z*DVecZ[index1].z);
        }

        if(RigidJ)
        {
          GradientStrainI(Gradient,1.0,fb,posB,comB);

          if(index_j2.x>=0) Gradient[index_j2.x]+=(fb.x*DVecX[index2].x+fb.y*DVecX[index2].y+fb.z*DVecX[index2].z);
          if(index_j2.y>=0) Gradient[index_j2.y]+=(fb.x*DVecY[index2].x+fb.y*DVecY[index2].y+fb.z*DVecY[index2].z);
          if(index_j2.z>=0) Gradient[index_j2.z]+=(fb.x*DVecZ[index2].x+fb.y*DVecZ[index2].y+fb.z*DVecZ[index2].z);
        }

        if(RigidK)
        {
          GradientStrainI(Gradient,1.0,fc,posC,comC);

          if(index_k2.x>=0) Gradient[index_k2.x]+=(fc.x*DVecX[index3].x+fc.y*DVecX[index3].y+fc.z*DVecX[index3].z);
          if(index_k2.y>=0) Gradient[index_k2.y]+=(fc.x*DVecY[index3].x+fc.y*DVecY[index3].y+fc.z*DVecY[index3].z);
          if(index_k2.z>=0) Gradient[index_k2.z]+=(fc.x*DVecZ[index3].x+fc.y*DVecZ[index3].y+fc.z*DVecZ[index3].z);
        }

        if(RigidL)
        {
          GradientStrainI(Gradient,1.0,fd,posD,comD);

          if(index_l2.x>=0) Gradient[index_l2.x]+=(fd.x*DVecX[index4].x+fd.y*DVecX[index4].y+fd.z*DVecX[index4].z);
          if(index_l2.y>=0) Gradient[index_l2.y]+=(fd.x*DVecY[index4].x+fd.y*DVecY[index4].y+fd.z*DVecY[index4].z);
          if(index_l2.z>=0) Gradient[index_l2.z]+=(fd.x*DVecZ[index4].x+fd.y*DVecZ[index4].y+fd.z*DVecZ[index4].z);
        }

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

void CalculateAdsorbateImproperTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,A,B,C,D,Type;
  POINT posA,posB,posC,posD;
  POINT comA,comB,comC,comD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL *parms;
  int RigidI,RigidJ,RigidK,RigidL;
  int grpA,grpB,grpC,grpD;
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  INT_VECTOR3 index_i2,index_j2,index_k2,index_l2;
  int index1,index2,index3,index4;
  REAL_MATRIX3x3 D2I,D2L,D2IL,D2IJ,D2IK,D2JL;
  REAL_MATRIX3x3 D2KL,D2J,D2K,D2JK,S;
  VECTOR DIL,DDJ,DDK,DEJ,DEK;
  VECTOR fa,fb,fc,fd;
  REAL DF,DDF;
  VECTOR vec_u,vec_v,vec_w;
  REAL u,w,v;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfImproperTorsions;i++)
    {
      A=Components[Type].ImproperTorsions[i].A;
      B=Components[Type].ImproperTorsions[i].B;
      C=Components[Type].ImproperTorsions[i].C;
      D=Components[Type].ImproperTorsions[i].D;
      parms=(REAL*)&Components[Type].ImproperTorsionArguments[i];

      posA=comA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=comB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC=comC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
      posD=comD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

      index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
      index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;
      index_l=Adsorbates[CurrentSystem][m].Atoms[D].HessianIndex;
      index_i2=UNDEFINED_INT_VECTOR3;
      index_j2=UNDEFINED_INT_VECTOR3;
      index_k2=UNDEFINED_INT_VECTOR3;
      index_l2=UNDEFINED_INT_VECTOR3;
      index1=-1;
      index2=-1;
      index3=-1;
      index4=-1;

      grpA=Components[Type].group[A];
      RigidI=Components[Type].Groups[grpA].Rigid;
      index1=Adsorbates[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      if(RigidI)
      {
        comA=Adsorbates[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
        index_i=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndex;
        index_i2=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
      }

      grpB=Components[Type].group[B];
      RigidJ=Components[Type].Groups[grpB].Rigid;
      index2=Adsorbates[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      if(RigidJ)
      {
        comB=Adsorbates[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
        index_j=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndex;
        index_j2=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
      }

      grpC=Components[Type].group[C];
      RigidK=Components[Type].Groups[grpC].Rigid;
      index3=Adsorbates[CurrentSystem][m].Atoms[C].HessianAtomIndex;
      if(RigidK)
      {
        comC=Adsorbates[CurrentSystem][m].Groups[grpC].CenterOfMassPosition;
        index_k=Adsorbates[CurrentSystem][m].Groups[grpC].HessianIndex;
        index_k2=Adsorbates[CurrentSystem][m].Groups[grpC].HessianIndexOrientation;
      }

      grpD=Components[Type].group[D];
      RigidL=Components[Type].Groups[grpD].Rigid;
      index4=Adsorbates[CurrentSystem][m].Atoms[D].HessianAtomIndex;
      if(RigidL)
      {
        comD=Adsorbates[CurrentSystem][m].Groups[grpD].CenterOfMassPosition;
        index_l=Adsorbates[CurrentSystem][m].Groups[grpD].HessianIndex;
        index_l2=Adsorbates[CurrentSystem][m].Groups[grpD].HessianIndexOrientation;
      }

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

      switch(Components[Type].ImproperTorsionType[i])
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
        case FIXED_IMPROPER_DIHEDRAL:
          U=DF=DDF=0.0;
          break;
        default:
          fprintf(stderr, "Undefined Improper Torsion potential\n");
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

        if(RigidI)
        {
          GradientStrainI(Gradient,1.0,fa,posA,comA);

          if(index_i2.x>=0) Gradient[index_i2.x]+=(fa.x*DVecX[index1].x+fa.y*DVecX[index1].y+fa.z*DVecX[index1].z);
          if(index_i2.y>=0) Gradient[index_i2.y]+=(fa.x*DVecY[index1].x+fa.y*DVecY[index1].y+fa.z*DVecY[index1].z);
          if(index_i2.z>=0) Gradient[index_i2.z]+=(fa.x*DVecZ[index1].x+fa.y*DVecZ[index1].y+fa.z*DVecZ[index1].z);
        }

        if(RigidJ)
        {
          GradientStrainI(Gradient,1.0,fb,posB,comB);

          if(index_j2.x>=0) Gradient[index_j2.x]+=(fb.x*DVecX[index2].x+fb.y*DVecX[index2].y+fb.z*DVecX[index2].z);
          if(index_j2.y>=0) Gradient[index_j2.y]+=(fb.x*DVecY[index2].x+fb.y*DVecY[index2].y+fb.z*DVecY[index2].z);
          if(index_j2.z>=0) Gradient[index_j2.z]+=(fb.x*DVecZ[index2].x+fb.y*DVecZ[index2].y+fb.z*DVecZ[index2].z);
        }

        if(RigidK)
        {
          GradientStrainI(Gradient,1.0,fc,posC,comC);

          if(index_k2.x>=0) Gradient[index_k2.x]+=(fc.x*DVecX[index3].x+fc.y*DVecX[index3].y+fc.z*DVecX[index3].z);
          if(index_k2.y>=0) Gradient[index_k2.y]+=(fc.x*DVecY[index3].x+fc.y*DVecY[index3].y+fc.z*DVecY[index3].z);
          if(index_k2.z>=0) Gradient[index_k2.z]+=(fc.x*DVecZ[index3].x+fc.y*DVecZ[index3].y+fc.z*DVecZ[index3].z);
        }

        if(RigidL)
        {
          GradientStrainI(Gradient,1.0,fd,posD,comD);

          if(index_l2.x>=0) Gradient[index_l2.x]+=(fd.x*DVecX[index4].x+fd.y*DVecX[index4].y+fd.z*DVecX[index4].z);
          if(index_l2.y>=0) Gradient[index_l2.y]+=(fd.x*DVecY[index4].x+fd.y*DVecY[index4].y+fd.z*DVecY[index4].z);
          if(index_l2.z>=0) Gradient[index_l2.z]+=(fd.x*DVecZ[index4].x+fd.y*DVecZ[index4].y+fd.z*DVecZ[index4].z);
        }

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

void CalculateCationImproperTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,A,B,C,D,Type;
  POINT posA,posB,posC,posD;
  POINT comA,comB,comC,comD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL *parms;
  int RigidI,RigidJ,RigidK,RigidL;
  int grpA,grpB,grpC,grpD;
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  INT_VECTOR3 index_i2,index_j2,index_k2,index_l2;
  int index1,index2,index3,index4;
  REAL_MATRIX3x3 D2I,D2L,D2IL,D2IJ,D2IK,D2JL;
  REAL_MATRIX3x3 D2KL,D2J,D2K,D2JK,S;
  VECTOR DIL,DDJ,DDK,DEJ,DEK;
  VECTOR fa,fb,fc,fd;
  REAL DF,DDF;
  VECTOR vec_u,vec_v,vec_w;
  REAL u,w,v;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfImproperTorsions;i++)
    {
      A=Components[Type].ImproperTorsions[i].A;
      B=Components[Type].ImproperTorsions[i].B;
      C=Components[Type].ImproperTorsions[i].C;
      D=Components[Type].ImproperTorsions[i].D;
      parms=(REAL*)&Components[Type].ImproperTorsionArguments[i];

      posA=comA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=comB=Cations[CurrentSystem][m].Atoms[B].Position;
      posC=comC=Cations[CurrentSystem][m].Atoms[C].Position;
      posD=comD=Cations[CurrentSystem][m].Atoms[D].Position;

      index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Cations[CurrentSystem][m].Atoms[B].HessianIndex;
      index_k=Cations[CurrentSystem][m].Atoms[C].HessianIndex;
      index_l=Cations[CurrentSystem][m].Atoms[D].HessianIndex;
      index_i2=UNDEFINED_INT_VECTOR3;
      index_j2=UNDEFINED_INT_VECTOR3;
      index_k2=UNDEFINED_INT_VECTOR3;
      index_l2=UNDEFINED_INT_VECTOR3;
      index1=-1;
      index2=-1;
      index3=-1;
      index4=-1;

      grpA=Components[Type].group[A];
      RigidI=Components[Type].Groups[grpA].Rigid;
      index1=Cations[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      if(RigidI)
      {
        comA=Cations[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
        index_i=Cations[CurrentSystem][m].Groups[grpA].HessianIndex;
        index_i2=Cations[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
      }

      grpB=Components[Type].group[B];
      RigidJ=Components[Type].Groups[grpB].Rigid;
      index2=Cations[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      if(RigidJ)
      {
        comB=Cations[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
        index_j=Cations[CurrentSystem][m].Groups[grpB].HessianIndex;
        index_j2=Cations[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
      }

      grpC=Components[Type].group[C];
      RigidK=Components[Type].Groups[grpC].Rigid;
      index3=Cations[CurrentSystem][m].Atoms[C].HessianAtomIndex;
      if(RigidK)
      {
        comC=Cations[CurrentSystem][m].Groups[grpC].CenterOfMassPosition;
        index_k=Cations[CurrentSystem][m].Groups[grpC].HessianIndex;
        index_k2=Cations[CurrentSystem][m].Groups[grpC].HessianIndexOrientation;
      }

      grpD=Components[Type].group[D];
      RigidL=Components[Type].Groups[grpD].Rigid;
      index4=Cations[CurrentSystem][m].Atoms[D].HessianAtomIndex;
      if(RigidL)
      {
        comD=Cations[CurrentSystem][m].Groups[grpD].CenterOfMassPosition;
        index_l=Cations[CurrentSystem][m].Groups[grpD].HessianIndex;
        index_l2=Cations[CurrentSystem][m].Groups[grpD].HessianIndexOrientation;
      }

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

      switch(Components[Type].ImproperTorsionType[i])
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
        case FIXED_IMPROPER_DIHEDRAL:
          U=DF=DDF=0.0;
          break;
        default:
          fprintf(stderr, "Undefined Improper Torsion potential\n");
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

        if(RigidI)
        {
          GradientStrainI(Gradient,1.0,fa,posA,comA);

          if(index_i2.x>=0) Gradient[index_i2.x]+=(fa.x*DVecX[index1].x+fa.y*DVecX[index1].y+fa.z*DVecX[index1].z);
          if(index_i2.y>=0) Gradient[index_i2.y]+=(fa.x*DVecY[index1].x+fa.y*DVecY[index1].y+fa.z*DVecY[index1].z);
          if(index_i2.z>=0) Gradient[index_i2.z]+=(fa.x*DVecZ[index1].x+fa.y*DVecZ[index1].y+fa.z*DVecZ[index1].z);
        }

        if(RigidJ)
        {
          GradientStrainI(Gradient,1.0,fb,posB,comB);

          if(index_j2.x>=0) Gradient[index_j2.x]+=(fb.x*DVecX[index2].x+fb.y*DVecX[index2].y+fb.z*DVecX[index2].z);
          if(index_j2.y>=0) Gradient[index_j2.y]+=(fb.x*DVecY[index2].x+fb.y*DVecY[index2].y+fb.z*DVecY[index2].z);
          if(index_j2.z>=0) Gradient[index_j2.z]+=(fb.x*DVecZ[index2].x+fb.y*DVecZ[index2].y+fb.z*DVecZ[index2].z);
        }

        if(RigidK)
        {
          GradientStrainI(Gradient,1.0,fc,posC,comC);

          if(index_k2.x>=0) Gradient[index_k2.x]+=(fc.x*DVecX[index3].x+fc.y*DVecX[index3].y+fc.z*DVecX[index3].z);
          if(index_k2.y>=0) Gradient[index_k2.y]+=(fc.x*DVecY[index3].x+fc.y*DVecY[index3].y+fc.z*DVecY[index3].z);
          if(index_k2.z>=0) Gradient[index_k2.z]+=(fc.x*DVecZ[index3].x+fc.y*DVecZ[index3].y+fc.z*DVecZ[index3].z);
        }

        if(RigidL)
        {
          GradientStrainI(Gradient,1.0,fd,posD,comD);

          if(index_l2.x>=0) Gradient[index_l2.x]+=(fd.x*DVecX[index4].x+fd.y*DVecX[index4].y+fd.z*DVecX[index4].z);
          if(index_l2.y>=0) Gradient[index_l2.y]+=(fd.x*DVecY[index4].x+fd.y*DVecY[index4].y+fd.z*DVecY[index4].z);
          if(index_l2.z>=0) Gradient[index_l2.z]+=(fd.x*DVecZ[index4].x+fd.y*DVecZ[index4].y+fd.z*DVecZ[index4].z);
        }

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

// NEW

void CalculateAdsorbateBondBondHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,Type,NumberOfBondBonds,A,B,C;
  REAL U,DF,DDF,r,rr,temp,temp2,exp_term,energy;
  REAL DFA,DFC,DDFA,DDFC;
  REAL *parms;
  POINT posA,posB,posC;
  VECTOR dr;
  REAL_MATRIX3x3 Hessian,Hessian1,Hessian2,Hessian3;
  VECTOR Rab,Rbc;
  REAL rab,rbc;
  INT_VECTOR3 index_i,index_j,index_k;
  int index1,index2,index3;
  REAL F,DA,DC;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfBondBonds=Components[Type].NumberOfBondBonds;
    for(i=0;i<NumberOfBondBonds;i++)
    {
      A=Components[Type].BondBonds[i].A;
      B=Components[Type].BondBonds[i].B;
      C=Components[Type].BondBonds[i].C;
      parms=(REAL*)&Components[Type].BondBondArguments[i];

      index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
      index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;

      index1=Adsorbates[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      index2=Adsorbates[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      index3=Adsorbates[CurrentSystem][m].Atoms[C].HessianAtomIndex;

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
          // =============================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2     [A]
          energy=parms[0]*(rab-parms[1])*(rbc-parms[2]);
          F=parms[0];
          DA=1.0;
          DC=1.0;
          DFA=parms[0]*(rbc-parms[2])/rab;
          DFC=parms[0]*(rab-parms[1])/rbc;
          DDFA=-parms[0]*(rbc-parms[2])/CUBE(rab);
          DDFC=-parms[0]*(rab-parms[1])/CUBE(rbc);
          break;
        case HARMONIC_BOND_BOND_CROSS:
          // p_0*SQR(rab-p_1)*SQR(rbc-p_2)
          // =============================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2     [A]
          energy=0.5*parms[0]*SQR(rab-parms[1])*SQR(rbc-parms[2]);
          F=0.5*parms[0];
          DA=2.0*(rab-parms[1]);
          DC=2.0*(rbc-parms[2]);
          DFA=2.0*F*SQR(rbc-parms[2])*(rab-parms[1])/rab;
          DFC=2.0*F*SQR(rab-parms[1])*(rbc-parms[2])/rbc;
          DDFA=2.0*F*parms[1]*SQR(rbc-parms[2])/CUBE(rab);
          DDFC=2.0*F*parms[2]*SQR(rab-parms[1])/CUBE(rbc);
          break;
        default:
          fprintf(stderr, "Undefined Bond-Bond potential in routine 'CalculateAdsorbateBondBondForce' ('internal_force.c')\n");
          exit(0);
          break;
      }

      // add contribution to the energy
      *Energy+=energy;

      // add contribution to the first derivatives
      if(ComputeGradient)
      {
        Gradient[index_i.x]+=DFA*Rab.x;
        Gradient[index_i.y]+=DFA*Rab.y;
        Gradient[index_i.z]+=DFA*Rab.z;

        Gradient[index_j.x]-=DFA*Rab.x+DFC*Rbc.x;
        Gradient[index_j.y]-=DFA*Rab.y+DFC*Rbc.y;
        Gradient[index_j.z]-=DFA*Rab.z+DFC*Rbc.z;

        Gradient[index_k.x]+=DFC*Rbc.x;
        Gradient[index_k.y]+=DFC*Rbc.y;
        Gradient[index_k.z]+=DFC*Rbc.z;

        GradientStrain(Gradient,DFA,Rab);
        GradientStrain(Gradient,DFC,Rbc);
      }

      // add contribution to the strain derivative tensor
      StrainDerivativeTensor->ax-=DFA*Rab.x*Rab.x+DFC*Rbc.x*Rbc.x;
      StrainDerivativeTensor->bx-=DFA*Rab.y*Rab.x+DFC*Rbc.y*Rbc.x;
      StrainDerivativeTensor->cx-=DFA*Rab.z*Rab.x+DFC*Rbc.z*Rbc.x;

      StrainDerivativeTensor->ay-=DFA*Rab.x*Rab.y+DFC*Rbc.x*Rbc.y;
      StrainDerivativeTensor->by-=DFA*Rab.y*Rab.y+DFC*Rbc.y*Rbc.y;
      StrainDerivativeTensor->cy-=DFA*Rab.z*Rab.y+DFC*Rbc.z*Rbc.y;

      StrainDerivativeTensor->az-=DFA*Rab.x*Rab.z+DFC*Rbc.x*Rbc.z;
      StrainDerivativeTensor->bz-=DFA*Rab.y*Rab.z+DFC*Rbc.y*Rbc.z;
      StrainDerivativeTensor->cz-=DFA*Rab.z*Rab.z+DFC*Rbc.z*Rbc.z;


      if(ComputeHessian)
      {
        // add contribution to the second derivatives (Hessian matrix)
        Hessian1.ax=(DDFA*Rab.x*Rab.x+DFA);
        Hessian1.ay=(DDFA*Rab.x*Rab.y);
        Hessian1.az=(DDFA*Rab.x*Rab.z);
        Hessian1.by=(DDFA*Rab.y*Rab.y+DFA);
        Hessian1.bz=(DDFA*Rab.y*Rab.z);
        Hessian1.cz=(DDFA*Rab.z*Rab.z+DFA);

        Hessian2.ax=(DDFC*Rbc.x*Rbc.x+DFC);
        Hessian2.ay=(DDFC*Rbc.x*Rbc.y);
        Hessian2.az=(DDFC*Rbc.x*Rbc.z);
        Hessian2.by=(DDFC*Rbc.y*Rbc.y+DFC);
        Hessian2.bz=(DDFC*Rbc.y*Rbc.z);
        Hessian2.cz=(DDFC*Rbc.z*Rbc.z+DFC);

        Hessian3.ax=F*DA*DC*Rab.x*Rbc.x/(rab*rbc);
        Hessian3.ay=F*DA*DC*Rab.x*Rbc.y/(rab*rbc);
        Hessian3.az=F*DA*DC*Rab.x*Rbc.z/(rab*rbc);
        Hessian3.bx=F*DA*DC*Rab.y*Rbc.x/(rab*rbc);
        Hessian3.by=F*DA*DC*Rab.y*Rbc.y/(rab*rbc);
        Hessian3.bz=F*DA*DC*Rab.y*Rbc.z/(rab*rbc);
        Hessian3.cx=F*DA*DC*Rab.z*Rbc.x/(rab*rbc);
        Hessian3.cy=F*DA*DC*Rab.z*Rbc.y/(rab*rbc);
        Hessian3.cz=F*DA*DC*Rab.z*Rbc.z/(rab*rbc);

        if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x]+=Hessian1.ax;
        if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y]+=Hessian1.ay;
        if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z]+=Hessian1.az;
        if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y]+=Hessian1.by;
        if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z]+=Hessian1.bz;
        if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z]+=Hessian1.cz;

        if((index_j.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_j.x][index_j.x]+=Hessian1.ax+2.0*Hessian3.ax+Hessian2.ax;
        if((index_j.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.x][index_j.y]+=Hessian1.ay+Hessian3.ay+Hessian3.bx+Hessian2.ay;
        if((index_j.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.x][index_j.z]+=Hessian1.az+Hessian3.az+Hessian3.cx+Hessian2.az;
        if((index_j.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.y][index_j.y]+=Hessian1.by+2.0*Hessian3.by+Hessian2.by;
        if((index_j.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.y][index_j.z]+=Hessian1.bz+Hessian3.bz+Hessian3.cy+Hessian2.bz;
        if((index_j.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.z][index_j.z]+=Hessian1.cz+2.0*Hessian3.cz+Hessian2.cz;

        if((index_k.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_k.x][index_k.x]+=Hessian2.ax;
        if((index_k.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_k.x][index_k.y]+=Hessian2.ay;
        if((index_k.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.x][index_k.z]+=Hessian2.az;
        if((index_k.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_k.y][index_k.y]+=Hessian2.by;
        if((index_k.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.y][index_k.z]+=Hessian2.bz;
        if((index_k.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.z][index_k.z]+=Hessian2.cz;

        if(index1<index2)
        {
          if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]-=Hessian1.ax+Hessian3.ax;
          if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]-=Hessian1.ay+Hessian3.ay;
          if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]-=Hessian1.az+Hessian3.az;
          if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]-=Hessian1.ay+Hessian3.bx;
          if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]-=Hessian1.by+Hessian3.by;
          if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]-=Hessian1.bz+Hessian3.bz;
          if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]-=Hessian1.az+Hessian3.cx;
          if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]-=Hessian1.bz+Hessian3.cy;
          if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]-=Hessian1.cz+Hessian3.cz;
        }
        else
        {
          if((index_j.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.x][index_i.x]-=Hessian1.ax+Hessian3.ax;
          if((index_j.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.x][index_i.y]-=Hessian1.ay+Hessian3.bx;
          if((index_j.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.x][index_i.z]-=Hessian1.az+Hessian3.cx;
          if((index_j.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.y][index_i.x]-=Hessian1.ay+Hessian3.ay;
          if((index_j.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.y][index_i.y]-=Hessian1.by+Hessian3.by;
          if((index_j.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.y][index_i.z]-=Hessian1.bz+Hessian3.cy;
          if((index_j.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.z][index_i.x]-=Hessian1.az+Hessian3.az;
          if((index_j.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.z][index_i.y]-=Hessian1.bz+Hessian3.bz;
          if((index_j.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.z][index_i.z]-=Hessian1.cz+Hessian3.cz;
        }

        if(index1<index3)
        {
          if((index_i.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.x][index_k.x]+=Hessian3.ax;
          if((index_i.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.x][index_k.y]+=Hessian3.ay;
          if((index_i.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.x][index_k.z]+=Hessian3.az;
          if((index_i.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.y][index_k.x]+=Hessian3.bx;
          if((index_i.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.y][index_k.y]+=Hessian3.by;
          if((index_i.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.y][index_k.z]+=Hessian3.bz;
          if((index_i.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.z][index_k.x]+=Hessian3.cx;
          if((index_i.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.z][index_k.y]+=Hessian3.cy;
          if((index_i.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.z][index_k.z]+=Hessian3.cz;
        }
        else
        {
          if((index_k.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.x][index_i.x]+=Hessian3.ax;
          if((index_k.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.x][index_i.y]+=Hessian3.bx;
          if((index_k.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.x][index_i.z]+=Hessian3.cx;
          if((index_k.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.y][index_i.x]+=Hessian3.ay;
          if((index_k.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.y][index_i.y]+=Hessian3.by;
          if((index_k.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.y][index_i.z]+=Hessian3.cy;
          if((index_k.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.z][index_i.x]+=Hessian3.az;
          if((index_k.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.z][index_i.y]+=Hessian3.bz;
          if((index_k.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.z][index_i.z]+=Hessian3.cz;
        }

        if(index2<index3)
        {
          if((index_j.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.x][index_k.x]-=Hessian2.ax+Hessian3.ax;
          if((index_j.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.x][index_k.y]-=Hessian2.ay+Hessian3.ay;
          if((index_j.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.x][index_k.z]-=Hessian2.az+Hessian3.az;
          if((index_j.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.y][index_k.x]-=Hessian2.ay+Hessian3.bx;
          if((index_j.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.y][index_k.y]-=Hessian2.by+Hessian3.by;
          if((index_j.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.y][index_k.z]-=Hessian2.bz+Hessian3.bz;
          if((index_j.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.z][index_k.x]-=Hessian2.az+Hessian3.cx;
          if((index_j.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.z][index_k.y]-=Hessian2.bz+Hessian3.cy;
          if((index_j.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.z][index_k.z]-=Hessian2.cz+Hessian3.cz;
        }
        else
        {
          if((index_k.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.x][index_j.x]-=Hessian2.ax+Hessian3.ax;
          if((index_k.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.x][index_j.y]-=Hessian2.ay+Hessian3.bx;
          if((index_k.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.x][index_j.z]-=Hessian2.az+Hessian3.cx;
          if((index_k.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.y][index_j.x]-=Hessian2.ay+Hessian3.ay;
          if((index_k.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.y][index_j.y]-=Hessian2.by+Hessian3.by;
          if((index_k.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.y][index_j.z]-=Hessian2.bz+Hessian3.cy;
          if((index_k.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.z][index_j.x]-=Hessian2.az+Hessian3.az;
          if((index_k.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.z][index_j.y]-=Hessian2.bz+Hessian3.bz;
          if((index_k.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.z][index_j.z]-=Hessian2.cz+Hessian3.cz;
        }

        // add contribution to the cross term
        //HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,DF,DDF,dr);
        //HessianAtomicStrainStrainLocal(HessianMatrix,index_i,index_j,DF,DDF,dr);
      }
    }
  }
}




void CalculateAdsorbateIntraVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,NumberOfIntraVDW,Type,typeA,typeB,A,B;
  REAL rr,energy,DF,DDF,Scaling;
  VECTOR dr;
  POINT posA,posB,comA,comB;
  int grpA,grpB;
  INT_VECTOR3 index_i,index_j;
  INT_VECTOR3 index_i2,index_j2;
  int index1,index2;
  int RigidI,RigidJ;
  REAL scalingA,scalingB;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfIntraVDW=Components[Type].NumberOfIntraVDW;
    for(i=0;i<NumberOfIntraVDW;i++)
    {
      A=MIN2(Components[Type].IntraVDW[i].A,Components[Type].IntraVDW[i].B);
      B=MAX2(Components[Type].IntraVDW[i].A,Components[Type].IntraVDW[i].B);
      Scaling=Components[Type].IntraVDWScaling[i];
      posA=comA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=comB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      // note: a cutoff is also used for intra-van der Waals forces
      if(rr<CutOffVDWSquared)
      {
        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
        index_i2=UNDEFINED_INT_VECTOR3;
        index_j2=UNDEFINED_INT_VECTOR3;
        index1=-1;
        index2=-1;

        grpA=Components[Type].group[A];
        RigidI=Components[Type].Groups[grpA].Rigid;
        if(RigidI)
        {
          comA=Adsorbates[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
          index_i=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndex;
          index_i2=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
          index1=Adsorbates[CurrentSystem][m].Atoms[A].HessianAtomIndex;
        }

        grpB=Components[Type].group[B];
        RigidJ=Components[Type].Groups[grpB].Rigid;
        if(RigidJ)
        {
          comB=Adsorbates[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
          index_j=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndex;
          index_j2=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
          index2=Adsorbates[CurrentSystem][m].Atoms[B].HessianAtomIndex;
        }

        typeA=Adsorbates[CurrentSystem][m].Atoms[A].Type;
        typeB=Adsorbates[CurrentSystem][m].Atoms[B].Type;

        scalingA=Adsorbates[CurrentSystem][m].Atoms[A].CFVDWScalingParameter;
        scalingB=Adsorbates[CurrentSystem][m].Atoms[B].CFVDWScalingParameter;

        PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,scalingA*scalingB);
        energy*=Scaling;
        DF*=Scaling;
        DDF*=Scaling;

        // add contribution to the energy
        *Energy+=energy;

        //if((index_i<0)&&(index_j<0)) continue;

        // add contribution to the first derivatives
        if(ComputeGradient)
        {
          if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
          if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
          if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

          if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
          if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
          if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

          if(RigidI)
          {
            GradientStrainI(Gradient,DF,dr,posA,comA);

            if(index_i2.x>=0) Gradient[index_i2.x]+=DF*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
            if(index_i2.y>=0) Gradient[index_i2.y]+=DF*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
            if(index_i2.z>=0) Gradient[index_i2.z]+=DF*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
          }

          if(RigidJ)
          {
            GradientStrainJ(Gradient,DF,dr,posB,comB);

            if(index_j2.x>=0) Gradient[index_j2.x]-=DF*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
            if(index_j2.y>=0) Gradient[index_j2.y]-=DF*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
            if(index_j2.z>=0) Gradient[index_j2.z]-=DF*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
          }

          GradientStrain(Gradient,DF,dr);
        }

        // add contribution to the strain derivative tensor
        StrainDerivativeTensor->ax+=DF*dr.x*dr.x;
        StrainDerivativeTensor->bx+=DF*dr.y*dr.x;
        StrainDerivativeTensor->cx+=DF*dr.z*dr.x;

        StrainDerivativeTensor->ay+=DF*dr.x*dr.y;
        StrainDerivativeTensor->by+=DF*dr.y*dr.y;
        StrainDerivativeTensor->cy+=DF*dr.z*dr.y;

        StrainDerivativeTensor->az+=DF*dr.x*dr.z;
        StrainDerivativeTensor->bz+=DF*dr.y*dr.z;
        StrainDerivativeTensor->cz+=DF*dr.z*dr.z;


        if(ComputeHessian)
        {
          HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);
          HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);

          HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB,RigidI,RigidJ);
          HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,DF,DDF,posA,comA,posB,comB,dr);

          HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
          HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB);
        }
      }
    }
  }
}

void CalculateCationIntraVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,NumberOfIntraVDW,Type,typeA,typeB,A,B;
  REAL rr,energy,DF,DDF,Scaling;
  VECTOR dr;
  POINT posA,posB,comA,comB;
  int grpA,grpB;
  INT_VECTOR3 index_i,index_j;
  INT_VECTOR3 index_i2,index_j2;
  int index1,index2;
  int RigidI,RigidJ;
  REAL scalingA,scalingB;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    NumberOfIntraVDW=Components[Type].NumberOfIntraVDW;
    for(i=0;i<NumberOfIntraVDW;i++)
    {
      A=MIN2(Components[Type].IntraVDW[i].A,Components[Type].IntraVDW[i].B);
      B=MAX2(Components[Type].IntraVDW[i].A,Components[Type].IntraVDW[i].B);
      Scaling=Components[Type].IntraVDWScaling[i];
      posA=comA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=comB=Cations[CurrentSystem][m].Atoms[B].Position;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      // note: a cutoff is also used for intra-van der Waals forces
      if(rr<CutOffVDWSquared)
      {
        index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Cations[CurrentSystem][m].Atoms[B].HessianIndex;
        index_i2=UNDEFINED_INT_VECTOR3;
        index_j2=UNDEFINED_INT_VECTOR3;
        index1=-1;
        index2=-1;

        grpA=Components[Type].group[A];
        RigidI=Components[Type].Groups[grpA].Rigid;
        if(RigidI)
        {
          comA=Cations[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
          index_i=Cations[CurrentSystem][m].Groups[grpA].HessianIndex;
          index_i2=Cations[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
          index1=Cations[CurrentSystem][m].Atoms[A].HessianAtomIndex;
        }

        grpB=Components[Type].group[B];
        RigidJ=Components[Type].Groups[grpB].Rigid;
        if(RigidJ)
        {
          comB=Cations[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
          index_j=Cations[CurrentSystem][m].Groups[grpB].HessianIndex;
          index_j2=Cations[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
          index2=Cations[CurrentSystem][m].Atoms[B].HessianAtomIndex;
        }

        typeA=Cations[CurrentSystem][m].Atoms[A].Type;
        typeB=Cations[CurrentSystem][m].Atoms[B].Type;

        scalingA=Cations[CurrentSystem][m].Atoms[A].CFVDWScalingParameter;
        scalingB=Cations[CurrentSystem][m].Atoms[B].CFVDWScalingParameter;

        PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,scalingA*scalingB);
        energy*=Scaling;
        DF*=Scaling;
        DDF*=Scaling;

        // add contribution to the energy
        *Energy+=energy;

        //if((index_i<0)&&(index_j<0)) continue;

        // add contribution to the first derivatives
        if(ComputeGradient)
        {
          if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
          if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
          if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

          if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
          if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
          if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

          if(RigidI)
          {
            GradientStrainI(Gradient,DF,dr,posA,comA);

            if(index_i2.x>=0) Gradient[index_i2.x]+=DF*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
            if(index_i2.y>=0) Gradient[index_i2.y]+=DF*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
            if(index_i2.z>=0) Gradient[index_i2.z]+=DF*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
          }

          if(RigidJ)
          {
            GradientStrainJ(Gradient,DF,dr,posB,comB);

            if(index_j2.x>=0) Gradient[index_j2.x]-=DF*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
            if(index_j2.y>=0) Gradient[index_j2.y]-=DF*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
            if(index_j2.z>=0) Gradient[index_j2.z]-=DF*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
          }

          GradientStrain(Gradient,DF,dr);
        }

        // add contribution to the strain derivative tensor
        StrainDerivativeTensor->ax+=DF*dr.x*dr.x;
        StrainDerivativeTensor->bx+=DF*dr.y*dr.x;
        StrainDerivativeTensor->cx+=DF*dr.z*dr.x;

        StrainDerivativeTensor->ay+=DF*dr.x*dr.y;
        StrainDerivativeTensor->by+=DF*dr.y*dr.y;
        StrainDerivativeTensor->cy+=DF*dr.z*dr.y;

        StrainDerivativeTensor->az+=DF*dr.x*dr.z;
        StrainDerivativeTensor->bz+=DF*dr.y*dr.z;
        StrainDerivativeTensor->cz+=DF*dr.z*dr.z;

        if(ComputeHessian)
        {
          HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);
          HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);

          HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB,RigidI,RigidJ);
          HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,DF,DDF,posA,comA,posB,comB,dr);

          HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
          HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB);
        }
      }
    }
  }
}

void CalculateAdsorbateIntraCoulombHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,NumberOfIntraCoulomb,Type,typeA,typeB,A,B;
  REAL r,rr;
  VECTOR dr;
  POINT posA,posB,comA,comB;
  REAL ChargeA,ChargeB;
  REAL energy,DF,DDF,Scaling;
  int grpA,grpB;
  INT_VECTOR3 index_i,index_j;
  INT_VECTOR3 index_i2,index_j2;
  int index1,index2;
  int RigidI,RigidJ;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfIntraCoulomb=Components[Type].NumberOfIntraChargeCharge;
    for(i=0;i<NumberOfIntraCoulomb;i++)
    {
      A=MIN2(Components[Type].IntraChargeCharge[i].A,Components[Type].IntraChargeCharge[i].B);
      B=MAX2(Components[Type].IntraChargeCharge[i].A,Components[Type].IntraChargeCharge[i].B);
      Scaling=Components[Type].IntraChargeChargeScaling[i];

      posA=comA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB=comB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

      typeA=Adsorbates[CurrentSystem][m].Atoms[A].Type;
      ChargeA=Adsorbates[CurrentSystem][m].Atoms[A].Charge;
      typeB=Adsorbates[CurrentSystem][m].Atoms[B].Type;
      ChargeB=Adsorbates[CurrentSystem][m].Atoms[B].Charge;

      index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
      index_i2=UNDEFINED_INT_VECTOR3;
      index_j2=UNDEFINED_INT_VECTOR3;
      index1=-1;
      index2=-1;

      grpA=Components[Type].group[A];
      RigidI=Components[Type].Groups[grpA].Rigid;
      if(RigidI)
      {
        comA=Adsorbates[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
        index_i=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndex;
        index_i2=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
        index1=Adsorbates[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      }

      grpB=Components[Type].group[B];
      RigidJ=Components[Type].Groups[grpB].Rigid;
      if(RigidJ)
      {
        comB=Adsorbates[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
        index_j=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndex;
        index_j2=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
        index2=Adsorbates[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      }

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      // note: no cutoff used here
      switch(ChargeMethod)
      {
        case NONE:
          energy=DF=DDF=0.0;
          break;
        case SHIFTED_COULOMB:
        case TRUNCATED_COULOMB:
        case SMOOTHED_COULOMB:
        case EWALD:
        default:
          energy=Scaling*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/r;
          DF=-Scaling*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
          DDF=Scaling*3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
          break;
      }

      // add contribution to the energy
      *Energy+=energy;

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

        if(RigidI)
        {
          GradientStrainI(Gradient,DF,dr,posA,comA);

          if(index_i2.x>=0) Gradient[index_i2.x]+=DF*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
          if(index_i2.y>=0) Gradient[index_i2.y]+=DF*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
          if(index_i2.z>=0) Gradient[index_i2.z]+=DF*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
        }

        if(RigidJ)
        {
          GradientStrainJ(Gradient,DF,dr,posB,comB);

          if(index_j2.x>=0) Gradient[index_j2.x]-=DF*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
          if(index_j2.y>=0) Gradient[index_j2.y]-=DF*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
          if(index_j2.z>=0) Gradient[index_j2.z]-=DF*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
        }

        GradientStrain(Gradient,DF,dr);
      }

      // add contribution to the strain derivative tensor
      StrainDerivativeTensor->ax+=DF*dr.x*dr.x;
      StrainDerivativeTensor->bx+=DF*dr.y*dr.x;
      StrainDerivativeTensor->cx+=DF*dr.z*dr.x;

      StrainDerivativeTensor->ay+=DF*dr.x*dr.y;
      StrainDerivativeTensor->by+=DF*dr.y*dr.y;
      StrainDerivativeTensor->cy+=DF*dr.z*dr.y;

      StrainDerivativeTensor->az+=DF*dr.x*dr.z;
      StrainDerivativeTensor->bz+=DF*dr.y*dr.z;
      StrainDerivativeTensor->cz+=DF*dr.z*dr.z;

      if(ComputeHessian)
      {
        HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);
        HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);

        HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB,RigidI,RigidJ);
        HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,DF,DDF,posA,comA,posB,comB,dr);

        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
        HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB);
      }
    }
  }
}

void CalculateCationIntraCoulombHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,m,NumberOfIntraCoulomb,Type,typeA,typeB,A,B;
  REAL r,rr;
  VECTOR dr;
  POINT posA,posB,comA,comB;
  REAL ChargeA,ChargeB;
  REAL energy,DF,DDF,Scaling;
  int grpA,grpB;
  INT_VECTOR3 index_i,index_j;
  INT_VECTOR3 index_i2,index_j2;
  int index1,index2;
  int RigidI,RigidJ;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    NumberOfIntraCoulomb=Components[Type].NumberOfIntraChargeCharge;
    for(i=0;i<NumberOfIntraCoulomb;i++)
    {
      A=MIN2(Components[Type].IntraChargeCharge[i].A,Components[Type].IntraChargeCharge[i].B);
      B=MAX2(Components[Type].IntraChargeCharge[i].A,Components[Type].IntraChargeCharge[i].B);
      Scaling=Components[Type].IntraChargeChargeScaling[i];

      typeA=Cations[CurrentSystem][m].Atoms[A].Type;
      ChargeA=Cations[CurrentSystem][m].Atoms[A].Charge;
      typeB=Cations[CurrentSystem][m].Atoms[B].Type;
      ChargeB=Cations[CurrentSystem][m].Atoms[B].Charge;

      posA=comA=Cations[CurrentSystem][m].Atoms[A].Position;
      posB=comB=Cations[CurrentSystem][m].Atoms[B].Position;

      index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Cations[CurrentSystem][m].Atoms[B].HessianIndex;
      index_i2=UNDEFINED_INT_VECTOR3;
      index_j2=UNDEFINED_INT_VECTOR3;
      index1=-1;
      index2=-1;

      grpA=Components[Type].group[A];
      RigidI=Components[Type].Groups[grpA].Rigid;
      if(RigidI)
      {
        comA=Cations[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
        index_i=Cations[CurrentSystem][m].Groups[grpA].HessianIndex;
        index_i2=Cations[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
        index1=Cations[CurrentSystem][m].Atoms[A].HessianAtomIndex;
      }

      grpB=Components[Type].group[B];
      RigidJ=Components[Type].Groups[grpB].Rigid;
      if(RigidJ)
      {
        comB=Cations[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
        index_j=Cations[CurrentSystem][m].Groups[grpB].HessianIndex;
        index_j2=Cations[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
        index2=Cations[CurrentSystem][m].Atoms[B].HessianAtomIndex;
      }

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      // note: no cutoff used here
      switch(ChargeMethod)
      {
        case NONE:
          energy=DF=DDF=0.0;
          break;
        case SHIFTED_COULOMB:
        case TRUNCATED_COULOMB:
        case SMOOTHED_COULOMB:
        case EWALD:
        default:
          energy=Scaling*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/r;
          DF=-Scaling*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
          DDF=Scaling*3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
          break;
      }

      // add contribution to the energy
      *Energy+=energy;

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

        if(RigidI)
        {
          GradientStrainI(Gradient,DF,dr,posA,comA);

          if(index_i2.x>=0) Gradient[index_i2.x]+=DF*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
          if(index_i2.y>=0) Gradient[index_i2.y]+=DF*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
          if(index_i2.z>=0) Gradient[index_i2.z]+=DF*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
        }

        if(RigidJ)
        {
          GradientStrainJ(Gradient,DF,dr,posB,comB);

          if(index_j2.x>=0) Gradient[index_j2.x]-=DF*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
          if(index_j2.y>=0) Gradient[index_j2.y]-=DF*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
          if(index_j2.z>=0) Gradient[index_j2.z]-=DF*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
        }

        GradientStrain(Gradient,DF,dr);
      }

      // add contribution to the strain derivative tensor
      StrainDerivativeTensor->ax+=DF*dr.x*dr.x;
      StrainDerivativeTensor->bx+=DF*dr.y*dr.x;
      StrainDerivativeTensor->cx+=DF*dr.z*dr.x;

      StrainDerivativeTensor->ay+=DF*dr.x*dr.y;
      StrainDerivativeTensor->by+=DF*dr.y*dr.y;
      StrainDerivativeTensor->cy+=DF*dr.z*dr.y;

      StrainDerivativeTensor->az+=DF*dr.x*dr.z;
      StrainDerivativeTensor->bz+=DF*dr.y*dr.z;
      StrainDerivativeTensor->cz+=DF*dr.z*dr.z;

      if(ComputeHessian)
      {
        HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);
        HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);

        HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB,RigidI,RigidJ);
        HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,DF,DDF,posA,comA,posB,comB,dr);

        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
        HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB);
      }
    }
  }
}

// Numerical routines
// TODO: convert them to analytical expressions

void CalculateAdsorbateInversionBendForces(int m,int i,VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,REAL *Energy,VECTOR *fa,VECTOR *fb,VECTOR *fc,VECTOR *fd,REAL_MATRIX3x3 *strain_derivative)
{
  int Type;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2,rrbc,rbc2,rbd2,rad2,rac2,dot;
  REAL CosChi,Chi;
  REAL temp,temp2;
  VECTOR Rab,Rbc,Rbd,Rad,Rac;
  REAL term,dedcos;
  VECTOR dccdia,dccdic,dccdid;
  VECTOR deedia,deedic,deedid;
  REAL energy;

  Type=Adsorbates[CurrentSystem][m].Type;
  parms=Components[Type].InversionBendArguments[i];

    Rab.x=posA.x-posB.x;
    Rab.y=posA.y-posB.y;
    Rab.z=posA.z-posB.z;
    rab2=Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z;
    rrab=sqrt(rab2);

    Rbc.x=posC.x-posB.x;
    Rbc.y=posC.y-posB.y;
    Rbc.z=posC.z-posB.z;
    rbc2=Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z;
    rrbc=sqrt(rbc2);

    Rbd.x=posD.x-posB.x;
    Rbd.y=posD.y-posB.y;
    Rbd.z=posD.z-posB.z;
    rbd2=Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z;

    Rac.x=posC.x-posA.x;
    Rac.y=posC.y-posA.y;
    Rac.z=posC.z-posA.z;
    rac2=Rac.x*Rac.x+Rac.y*Rac.y+Rac.z*Rac.z;

    Rad.x=posD.x-posA.x;
    Rad.y=posD.y-posA.y;
    Rad.z=posD.z-posA.z;
    rad2=Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z;

    switch(Components[Type].InversionBendType[i])
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

    switch(Components[Type].InversionBendType[i])
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

    switch(Components[Type].InversionBendType[i])
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

  // add contribution to the stress tensor
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


void CalculateAdsorbateInversionBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivativeTensor,int ComputeGradient,int ComputeHessian)
{
  int i,m;
  int A,B,C,D,Type;
  const REAL Delta=1e-7;
  VECTOR posA0,posB0,posC0,posD0;
  VECTOR posA,posB,posC,posD;
  VECTOR fa0,fb0,fc0,fd0,fa[4],fb[4],fc[4],fd[4];
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  REAL U;
  REAL_MATRIX3x3 strain_derivative;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfInversionBends;i++)
    {
      A=Components[Type].InversionBends[i].A;
      B=Components[Type].InversionBends[i].B;
      C=Components[Type].InversionBends[i].C;
      D=Components[Type].InversionBends[i].D;

      index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
      index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;
      index_l=Adsorbates[CurrentSystem][m].Atoms[D].HessianIndex;

      posA0=Adsorbates[CurrentSystem][m].Atoms[A].Position;
      posB0=Adsorbates[CurrentSystem][m].Atoms[B].Position;
      posC0=Adsorbates[CurrentSystem][m].Atoms[C].Position;
      posD0=Adsorbates[CurrentSystem][m].Atoms[D].Position;


      CalculateAdsorbateInversionBendForces(m,i,posA0,posB0,posC0,posD0,&U,&fa0,&fb0,&fc0,&fd0,&strain_derivative);

      *Energy+=U;

      //if((index_i<0)&&(index_j<0)&&(index_k<0)&&(index_l<0)) continue;

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
        if(index_k.x>=0) Gradient[index_k.y]-=fc0.y;
        if(index_k.x>=0) Gradient[index_k.z]-=fc0.z;

        if(index_l.x>=0) Gradient[index_l.x]-=fd0.x;
        if(index_l.y>=0) Gradient[index_l.y]-=fd0.y;
        if(index_l.z>=0) Gradient[index_l.z]-=fd0.z;
      }

      if(ComputeHessian)
      {
        // Atom A

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x+=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z+=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z-=0.5*Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z-=Delta;
        CalculateAdsorbateInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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


/// CATION
void CalculateCationInversionBendForces(int m,int i,VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,REAL *Energy,VECTOR *fa,VECTOR *fb,VECTOR *fc,VECTOR *fd,REAL_MATRIX3x3 *strain_derivative)
{
  int Type;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2,rrbc,rbc2,rbd2,rad2,rac2,dot;
  REAL CosChi,Chi;
  REAL temp,temp2;
  VECTOR Rab,Rbc,Rbd,Rad,Rac;
  REAL term,dedcos;
  VECTOR dccdia,dccdic,dccdid;
  VECTOR deedia,deedic,deedid;

  Type=Cations[CurrentSystem][m].Type;
  parms=Components[Type].InversionBendArguments[i];

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

  switch(Components[Type].InversionBendType[i])
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
      fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculatecationInversionBendForces' ('internal_hessian.c')\n");
      exit(0);
      break;
  }

  e=Rab.x*(Rbc.y*Rbd.z-Rbc.z*Rbd.y)+Rab.y*(Rbc.z*Rbd.x-Rbc.x*Rbd.z)+Rab.z*(Rbc.x*Rbd.y-Rbc.y*Rbd.x);
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
      *Energy=0.5*parms[0]*SQR(Chi-parms[1]);
      dedcos=-SIGN(1.0,e)*(parms[0]*(Chi-parms[1])/sqrt(c*(rab2-e*e/c)));
      break;
    case HARMONIC_COSINE_INVERSION:
    case HARMONIC_COSINE_INVERSION2:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      Chi=acos(CosChi);
      *Energy=0.5*parms[0]*SQR(CosChi-parms[1]);
      dedcos=SIGN(1.0,e)*parms[0]*(CosChi-parms[1])*sin(Chi)/sqrt(c*(rab2-e*e/c));
      break;
    case PLANAR_INVERSION:
    case PLANAR_INVERSION2:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      Chi=acos(CosChi);
      *Energy=parms[0]*(1.0-CosChi);
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
      *Energy=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
      dedcos=-SIGN(1.0,e)*parms[0]*temp*RAD2DEG*(2.0-3.0*0.014*temp+4.0*5.6e-5*temp2-5.0*7.0e-7*temp*temp2+6.0*2.2e-8*SQR(temp2))/sqrt(c*(rab2-e*e/c));
      break;
    default:
      fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateCationInversionBendForces' ('internal_hessian.c')\n");
      exit(0);
      break;
  }

  switch(Components[Type].InversionBendType[i])
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
      fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateCationInversionBendForces' ('internal_hessian.c')\n");
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

  fa->x=dedcos*(dccdia.x+deedia.x);
  fa->y=dedcos*(dccdia.y+deedia.y);
  fa->z=dedcos*(dccdia.z+deedia.z);
  fc->x=dedcos*(dccdic.x+deedic.x);
  fc->y=dedcos*(dccdic.y+deedic.y);
  fc->z=dedcos*(dccdic.z+deedic.z);
  fd->x=dedcos*(dccdid.x+deedid.x);
  fd->y=dedcos*(dccdid.y+deedid.y);
  fd->z=dedcos*(dccdid.z+deedid.z);

  fb->x=-(fa->x+fc->x+fd->x);
  fb->y=-(fa->y+fc->y+fd->y);
  fb->z=-(fa->z+fc->z+fd->z);

  // add contribution to the stress tensor
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


void CalculateCationInversionBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivativeTensor,int ComputeGradient,int ComputeHessian)
{
  int i,m;
  int A,B,C,D,Type;
  const REAL Delta=1e-7;
  VECTOR posA0,posB0,posC0,posD0;
  VECTOR posA,posB,posC,posD;
  VECTOR fa0,fb0,fc0,fd0,fa[4],fb[4],fc[4],fd[4];
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  REAL U;
  REAL_MATRIX3x3 strain_derivative;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfInversionBends;i++)
    {
      A=Components[Type].InversionBends[i].A;
      B=Components[Type].InversionBends[i].B;
      C=Components[Type].InversionBends[i].C;
      D=Components[Type].InversionBends[i].D;

      index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
      index_j=Cations[CurrentSystem][m].Atoms[B].HessianIndex;
      index_k=Cations[CurrentSystem][m].Atoms[C].HessianIndex;
      index_l=Cations[CurrentSystem][m].Atoms[D].HessianIndex;

      posA0=Cations[CurrentSystem][m].Atoms[A].Position;
      posB0=Cations[CurrentSystem][m].Atoms[B].Position;
      posC0=Cations[CurrentSystem][m].Atoms[C].Position;
      posD0=Cations[CurrentSystem][m].Atoms[D].Position;


      CalculateCationInversionBendForces(m,i,posA0,posB0,posC0,posD0,&U,&fa0,&fb0,&fc0,&fd0,&strain_derivative);

      *Energy+=U;

      //if((index_i<0)&&(index_j<0)&&(index_k<0)&&(index_l<0)) continue;

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
        if(index_k.x>=0) Gradient[index_k.y]-=fc0.y;
        if(index_k.x>=0) Gradient[index_k.z]-=fc0.z;

        if(index_l.x>=0) Gradient[index_l.x]-=fd0.x;
        if(index_l.y>=0) Gradient[index_l.y]-=fd0.y;
        if(index_l.z>=0) Gradient[index_l.z]-=fd0.z;
      }

      if(ComputeHessian)
      {
        // Atom A

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x+=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.x-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.y-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posA.z-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.x-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.y-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posB.z-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.x-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.y-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posC.z-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.x-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.y-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z+=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z-=0.5*Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

        posA=posA0; posB=posB0; posC=posC0; posD=posD0;
        posD.z-=Delta;
        CalculateCationInversionBendForces(m,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

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



void CalculateHarmonicBondConstraintHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int m;
  REAL U,DF,DDF,r,rr;
  REAL parms0,parms1;
  POINT posA,posB;
  VECTOR dr;
  REAL_MATRIX3x3 Hessian;
  INT_VECTOR3 index_i,index_j;

  for(m=0;m<NumberOfHarmonicDistanceConstraints[CurrentSystem];m++)
  {
    index_i=HarmonicDistanceConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=HarmonicDistanceConstraints[CurrentSystem][m][1]->HessianIndex;

    posA=HarmonicDistanceConstraints[CurrentSystem][m][0]->Position;
    posB=HarmonicDistanceConstraints[CurrentSystem][m][1]->Position;

    parms0=HarmonicDistanceConstraintParameters[CurrentSystem][m][0];
    parms1=HarmonicDistanceConstraintParameters[CurrentSystem][m][1];

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    fprintf(stderr, "Distance: %g\n",r);

    U=0.5*parms0*SQR(r-parms1);
    DF=parms0*(r-parms1)/r;
    DDF=(parms0*parms1)/(r*rr);

    // add contribution to the energy
    *Energy+=U;

    //if((index_i<0)&&(index_j<0)) continue;

    // add contribution to the first derivatives
    if(ComputeGradient)
    {
      if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
      if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
      if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

      if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
      if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
      if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

      GradientStrain(Gradient,DF,dr);
    }

    // add contribution to the strain derivative tensor
    StrainDerivativeTensor->ax+=dr.x*DF*dr.x;
    StrainDerivativeTensor->bx+=dr.y*DF*dr.x;
    StrainDerivativeTensor->cx+=dr.z*DF*dr.x;

    StrainDerivativeTensor->ay+=dr.x*DF*dr.y;
    StrainDerivativeTensor->by+=dr.y*DF*dr.y;
    StrainDerivativeTensor->cy+=dr.z*DF*dr.y;

    StrainDerivativeTensor->az+=dr.x*DF*dr.z;
    StrainDerivativeTensor->bz+=dr.y*DF*dr.z;
    StrainDerivativeTensor->cz+=dr.z*DF*dr.z;

    if(ComputeHessian)
    {
      // add contribution to the second derivatives (Hessian matrix)
      Hessian.ax=DDF*dr.x*dr.x+DF;
      Hessian.ay=DDF*dr.x*dr.y;
      Hessian.az=DDF*dr.x*dr.z;
      Hessian.by=DDF*dr.y*dr.y+DF;
      Hessian.bz=DDF*dr.y*dr.z;
      Hessian.cz=DDF*dr.z*dr.z+DF;

      if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x]+=Hessian.ax;
      if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y]+=Hessian.ay;
      if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z]+=Hessian.az;
      if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y]+=Hessian.by;
      if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z]+=Hessian.bz;
      if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z]+=Hessian.cz;

      if((index_j.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_j.x][index_j.x]+=Hessian.ax;
      if((index_j.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.x][index_j.y]+=Hessian.ay;
      if((index_j.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.x][index_j.z]+=Hessian.az;
      if((index_j.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.y][index_j.y]+=Hessian.by;
      if((index_j.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.y][index_j.z]+=Hessian.bz;
      if((index_j.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.z][index_j.z]+=Hessian.cz;

      if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.x)][MAX2(index_i.x,index_j.x)]-=Hessian.ax;
      if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.y)][MAX2(index_i.x,index_j.y)]-=Hessian.ay;
      if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.z)][MAX2(index_i.x,index_j.z)]-=Hessian.az;
      if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.x)][MAX2(index_i.y,index_j.x)]-=Hessian.ay;
      if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.y)][MAX2(index_i.y,index_j.y)]-=Hessian.by;
      if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.z)][MAX2(index_i.y,index_j.z)]-=Hessian.bz;
      if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.x)][MAX2(index_i.z,index_j.x)]-=Hessian.az;
      if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.y)][MAX2(index_i.z,index_j.y)]-=Hessian.bz;
      if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.z)][MAX2(index_i.z,index_j.z)]-=Hessian.cz;

      // add contribution to the cross term
      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,DF,DDF,dr);
      HessianAtomicStrainStrainLocal(HessianMatrix,DF,DDF,dr);
    }
  }
}

void CalculateHarmonicBendConstraintHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int m;
  REAL parms0,parms1,U;
  REAL CosTheta,Theta,SinTheta,temp;
  REAL rab,rbc,rac;
  POINT posA,posB,posC;
  VECTOR Rab,Rbc,Rac;
  INT_VECTOR3 index_i,index_j,index_k;
  REAL DTDX,DF,DDF;
  REAL_MATRIX3x3 D2I,D2K,D2IK,S;
  VECTOR dtA,dtB,dtC;
  VECTOR vec_u,vec_v;
  REAL u,v;

  S.ax=S.bx=S.cx=0.0;
  S.ay=S.by=S.cy=0.0;
  S.az=S.bz=S.cz=0.0;

  for(m=0;m<NumberOfHarmonicAngleConstraints[CurrentSystem];m++)
  {
    index_i=HarmonicAngleConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=HarmonicAngleConstraints[CurrentSystem][m][1]->HessianIndex;
    index_k=HarmonicAngleConstraints[CurrentSystem][m][2]->HessianIndex;

    posA=HarmonicAngleConstraints[CurrentSystem][m][0]->Position;
    posB=HarmonicAngleConstraints[CurrentSystem][m][1]->Position;
    posC=HarmonicAngleConstraints[CurrentSystem][m][2]->Position;

    parms0=HarmonicAngleConstraintParameters[CurrentSystem][m][0];
    parms1=HarmonicAngleConstraintParameters[CurrentSystem][m][1];

    Rab.x=posA.x-posB.x;
    Rab.y=posA.y-posB.y;
    Rab.z=posA.z-posB.z;
    vec_u=Rab;
    rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
    u=rab;
    Rab.x/=rab;
    Rab.y/=rab;
    Rab.z/=rab;

    Rbc.x=posC.x-posB.x;
    Rbc.y=posC.y-posB.y;
    Rbc.z=posC.z-posB.z;
    vec_v=Rbc;
    rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
    v=rbc;
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
    SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));
    DTDX=-1.0/SinTheta;

    temp=Theta-parms1;
    U=0.5*parms0*SQR(temp);
    DF=parms0*temp*DTDX;
    DDF=parms0*SQR(DTDX)+parms0*temp*CosTheta*CUBE(DTDX);

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

    if(ComputeGradient)
    {
      // add contribution to the first derivatives
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
      if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.x)][MAX2(index_i.x,index_j.x)]+=DDF*dtA.x*dtB.x-D2I.ax-D2IK.ax;
      if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.y)][MAX2(index_i.x,index_j.y)]+=DDF*dtA.x*dtB.y-D2I.ay-D2IK.ay;
      if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.z)][MAX2(index_i.x,index_j.z)]+=DDF*dtA.x*dtB.z-D2I.az-D2IK.az;
      if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.x)][MAX2(index_i.y,index_j.x)]+=DDF*dtA.y*dtB.x-D2I.ay-D2IK.bx;
      if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.y)][MAX2(index_i.y,index_j.y)]+=DDF*dtA.y*dtB.y-D2I.by-D2IK.by;
      if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.z)][MAX2(index_i.y,index_j.z)]+=DDF*dtA.y*dtB.z-D2I.bz-D2IK.bz;
      if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.x)][MAX2(index_i.z,index_j.x)]+=DDF*dtA.z*dtB.x-D2I.az-D2IK.cx;
      if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.y)][MAX2(index_i.z,index_j.y)]+=DDF*dtA.z*dtB.y-D2I.bz-D2IK.cy;
      if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.z)][MAX2(index_i.z,index_j.z)]+=DDF*dtA.z*dtB.z-D2I.cz-D2IK.cz;

      // Calculate the AC-block of the Hessian
      if((index_i.x>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_i.x,index_k.x)][MAX2(index_i.x,index_k.x)]+=DDF*dtA.x*dtC.x+D2IK.ax;
      if((index_i.x>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_i.x,index_k.y)][MAX2(index_i.x,index_k.y)]+=DDF*dtA.x*dtC.y+D2IK.ay;
      if((index_i.x>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_i.x,index_k.z)][MAX2(index_i.x,index_k.z)]+=DDF*dtA.x*dtC.z+D2IK.az;
      if((index_i.y>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_i.y,index_k.x)][MAX2(index_i.y,index_k.x)]+=DDF*dtA.y*dtC.x+D2IK.bx;
      if((index_i.y>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_i.y,index_k.y)][MAX2(index_i.y,index_k.y)]+=DDF*dtA.y*dtC.y+D2IK.by;
      if((index_i.y>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_i.y,index_k.z)][MAX2(index_i.y,index_k.z)]+=DDF*dtA.y*dtC.z+D2IK.bz;
      if((index_i.z>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_i.z,index_k.x)][MAX2(index_i.z,index_k.x)]+=DDF*dtA.z*dtC.x+D2IK.cx;
      if((index_i.z>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_i.z,index_k.y)][MAX2(index_i.z,index_k.y)]+=DDF*dtA.z*dtC.y+D2IK.cy;
      if((index_i.z>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_i.z,index_k.z)][MAX2(index_i.z,index_k.z)]+=DDF*dtA.z*dtC.z+D2IK.cz;

      // Calculate the BC-block of the Hessian
      if((index_j.x>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_j.x,index_k.x)][MAX2(index_j.x,index_k.x)]+=DDF*dtC.x*dtB.x-D2K.ax-D2IK.ax;
      if((index_j.x>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_j.x,index_k.y)][MAX2(index_j.x,index_k.y)]+=DDF*dtC.x*dtB.y-D2K.ay-D2IK.bx;
      if((index_j.x>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_j.x,index_k.z)][MAX2(index_j.x,index_k.z)]+=DDF*dtC.x*dtB.z-D2K.az-D2IK.cx;
      if((index_j.y>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_j.y,index_k.x)][MAX2(index_j.y,index_k.x)]+=DDF*dtC.y*dtB.x-D2K.ay-D2IK.ay;
      if((index_j.y>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_j.y,index_k.y)][MAX2(index_j.y,index_k.y)]+=DDF*dtC.y*dtB.y-D2K.by-D2IK.by;
      if((index_j.y>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_j.y,index_k.z)][MAX2(index_j.y,index_k.z)]+=DDF*dtC.y*dtB.z-D2K.bz-D2IK.cy;
      if((index_j.z>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_j.z,index_k.x)][MAX2(index_j.z,index_k.x)]+=DDF*dtC.z*dtB.x-D2K.az-D2IK.az;
      if((index_j.z>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_j.z,index_k.y)][MAX2(index_j.z,index_k.y)]+=DDF*dtC.z*dtB.y-D2K.bz-D2IK.bz;
      if((index_j.z>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_j.z,index_k.z)][MAX2(index_j.z,index_k.z)]+=DDF*dtC.z*dtB.z-D2K.cz-D2IK.cz;

      HessianBendStrainPosition(index_i,index_j,index_k,HessianMatrix,vec_u,vec_v,u,v,
                                rab,rbc,Rab,Rbc,dtA,dtC,DF,DDF,S,CosTheta);

      HessianBendStrainStrain(HessianMatrix,vec_u,vec_v,u,v,DF,DDF,S,CosTheta);
    }
  }
}

void CalculateHarmonicDihedralConstraintHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int m;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL parms0,parms1;
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  REAL_MATRIX3x3 D2I,D2L,D2IL,D2IJ,D2IK,D2JL;
  REAL_MATRIX3x3 D2KL,D2J,D2K,D2JK,S;
  VECTOR DIL,DDJ,DDK,DEJ,DEK;
  VECTOR fa,fb,fc,fd;
  REAL DF,DDF;
  VECTOR vec_u,vec_v,vec_w;
  REAL u,w,v;

  for(m=0;m<NumberOfHarmonicDihedralConstraints[CurrentSystem];m++)
  {
    index_i=HarmonicDihedralConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=HarmonicDihedralConstraints[CurrentSystem][m][1]->HessianIndex;
    index_k=HarmonicDihedralConstraints[CurrentSystem][m][2]->HessianIndex;
    index_l=HarmonicDihedralConstraints[CurrentSystem][m][3]->HessianIndex;

    posA=HarmonicDihedralConstraints[CurrentSystem][m][0]->Position;
    posB=HarmonicDihedralConstraints[CurrentSystem][m][1]->Position;
    posC=HarmonicDihedralConstraints[CurrentSystem][m][2]->Position;
    posD=HarmonicDihedralConstraints[CurrentSystem][m][3]->Position;

    parms0=HarmonicDihedralConstraintParameters[CurrentSystem][m][0];
    parms1=HarmonicDihedralConstraintParameters[CurrentSystem][m][1];

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
    fprintf(stderr, "Dihedral angle: %g %g %g\n",Phi*RAD2DEG,parms0,parms1*RAD2DEG);
    SinPhi=sin(Phi);
    Phi-=parms1;
    Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
    SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
    U=0.5*parms0*SQR(Phi);
    DF=-parms0*Phi/SinPhi;
    DDF=-parms0*(Phi*CosPhi-SinPhi)/CUBE(SinPhi);

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


      if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.x)][MAX2(index_i.x,index_j.x)]+=DDF*dtA.x*dtB.x+D2IJ.ax;
      if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.y)][MAX2(index_i.x,index_j.y)]+=DDF*dtA.x*dtB.y+D2IJ.ay;
      if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.z)][MAX2(index_i.x,index_j.z)]+=DDF*dtA.x*dtB.z+D2IJ.az;
      if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.x)][MAX2(index_i.y,index_j.x)]+=DDF*dtA.y*dtB.x+D2IJ.bx;
      if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.y)][MAX2(index_i.y,index_j.y)]+=DDF*dtA.y*dtB.y+D2IJ.by;
      if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.z)][MAX2(index_i.y,index_j.z)]+=DDF*dtA.y*dtB.z+D2IJ.bz;
      if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.x)][MAX2(index_i.z,index_j.x)]+=DDF*dtA.z*dtB.x+D2IJ.cx;
      if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.y)][MAX2(index_i.z,index_j.y)]+=DDF*dtA.z*dtB.y+D2IJ.cy;
      if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.z)][MAX2(index_i.z,index_j.z)]+=DDF*dtA.z*dtB.z+D2IJ.cz;

      if((index_i.x>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_i.x,index_k.x)][MAX2(index_i.x,index_k.x)]+=DDF*dtA.x*dtC.x+D2IK.ax;
      if((index_i.x>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_i.x,index_k.y)][MAX2(index_i.x,index_k.y)]+=DDF*dtA.x*dtC.y+D2IK.ay;
      if((index_i.x>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_i.x,index_k.z)][MAX2(index_i.x,index_k.z)]+=DDF*dtA.x*dtC.z+D2IK.az;
      if((index_i.y>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_i.y,index_k.x)][MAX2(index_i.y,index_k.x)]+=DDF*dtA.y*dtC.x+D2IK.bx;
      if((index_i.y>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_i.y,index_k.y)][MAX2(index_i.y,index_k.y)]+=DDF*dtA.y*dtC.y+D2IK.by;
      if((index_i.y>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_i.y,index_k.z)][MAX2(index_i.y,index_k.z)]+=DDF*dtA.y*dtC.z+D2IK.bz;
      if((index_i.z>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_i.z,index_k.x)][MAX2(index_i.z,index_k.x)]+=DDF*dtA.z*dtC.x+D2IK.cx;
      if((index_i.z>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_i.z,index_k.y)][MAX2(index_i.z,index_k.y)]+=DDF*dtA.z*dtC.y+D2IK.cy;
      if((index_i.z>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_i.z,index_k.z)][MAX2(index_i.z,index_k.z)]+=DDF*dtA.z*dtC.z+D2IK.cz;

      if((index_i.x>=0)&&(index_l.x>=0)) HessianMatrix.element[MIN2(index_i.x,index_l.x)][MAX2(index_i.x,index_l.x)]+=DDF*dtA.x*dtD.x+D2IL.ax;
      if((index_i.x>=0)&&(index_l.y>=0)) HessianMatrix.element[MIN2(index_i.x,index_l.y)][MAX2(index_i.x,index_l.y)]+=DDF*dtA.x*dtD.y+D2IL.ay;
      if((index_i.x>=0)&&(index_l.z>=0)) HessianMatrix.element[MIN2(index_i.x,index_l.z)][MAX2(index_i.x,index_l.z)]+=DDF*dtA.x*dtD.z+D2IL.az;
      if((index_i.y>=0)&&(index_l.x>=0)) HessianMatrix.element[MIN2(index_i.y,index_l.x)][MAX2(index_i.y,index_l.x)]+=DDF*dtA.y*dtD.x+D2IL.bx;
      if((index_i.y>=0)&&(index_l.y>=0)) HessianMatrix.element[MIN2(index_i.y,index_l.y)][MAX2(index_i.y,index_l.y)]+=DDF*dtA.y*dtD.y+D2IL.by;
      if((index_i.y>=0)&&(index_l.z>=0)) HessianMatrix.element[MIN2(index_i.y,index_l.z)][MAX2(index_i.y,index_l.z)]+=DDF*dtA.y*dtD.z+D2IL.bz;
      if((index_i.z>=0)&&(index_l.x>=0)) HessianMatrix.element[MIN2(index_i.z,index_l.x)][MAX2(index_i.z,index_l.x)]+=DDF*dtA.z*dtD.x+D2IL.cx;
      if((index_i.z>=0)&&(index_l.y>=0)) HessianMatrix.element[MIN2(index_i.z,index_l.y)][MAX2(index_i.z,index_l.y)]+=DDF*dtA.z*dtD.y+D2IL.cy;
      if((index_i.z>=0)&&(index_l.z>=0)) HessianMatrix.element[MIN2(index_i.z,index_l.z)][MAX2(index_i.z,index_l.z)]+=DDF*dtA.z*dtD.z+D2IL.cz;

      if((index_j.x>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_j.x,index_k.x)][MAX2(index_j.x,index_k.x)]+=DDF*dtB.x*dtC.x+D2JK.ax;
      if((index_j.x>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_j.x,index_k.y)][MAX2(index_j.x,index_k.y)]+=DDF*dtB.x*dtC.y+D2JK.ay;
      if((index_j.x>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_j.x,index_k.z)][MAX2(index_j.x,index_k.z)]+=DDF*dtB.x*dtC.z+D2JK.az;
      if((index_j.y>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_j.y,index_k.x)][MAX2(index_j.y,index_k.x)]+=DDF*dtB.y*dtC.x+D2JK.bx;
      if((index_j.y>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_j.y,index_k.y)][MAX2(index_j.y,index_k.y)]+=DDF*dtB.y*dtC.y+D2JK.by;
      if((index_j.y>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_j.y,index_k.z)][MAX2(index_j.y,index_k.z)]+=DDF*dtB.y*dtC.z+D2JK.bz;
      if((index_j.z>=0)&&(index_k.x>=0)) HessianMatrix.element[MIN2(index_j.z,index_k.x)][MAX2(index_j.z,index_k.x)]+=DDF*dtB.z*dtC.x+D2JK.cx;
      if((index_j.z>=0)&&(index_k.y>=0)) HessianMatrix.element[MIN2(index_j.z,index_k.y)][MAX2(index_j.z,index_k.y)]+=DDF*dtB.z*dtC.y+D2JK.cy;
      if((index_j.z>=0)&&(index_k.z>=0)) HessianMatrix.element[MIN2(index_j.z,index_k.z)][MAX2(index_j.z,index_k.z)]+=DDF*dtB.z*dtC.z+D2JK.cz;

      if((index_j.x>=0)&&(index_l.x>=0)) HessianMatrix.element[MIN2(index_j.x,index_l.x)][MAX2(index_j.x,index_l.x)]+=DDF*dtB.x*dtD.x+D2JL.ax;
      if((index_j.x>=0)&&(index_l.y>=0)) HessianMatrix.element[MIN2(index_j.x,index_l.y)][MAX2(index_j.x,index_l.y)]+=DDF*dtB.x*dtD.y+D2JL.ay;
      if((index_j.x>=0)&&(index_l.z>=0)) HessianMatrix.element[MIN2(index_j.x,index_l.z)][MAX2(index_j.x,index_l.z)]+=DDF*dtB.x*dtD.z+D2JL.az;
      if((index_j.y>=0)&&(index_l.x>=0)) HessianMatrix.element[MIN2(index_j.y,index_l.x)][MAX2(index_j.y,index_l.x)]+=DDF*dtB.y*dtD.x+D2JL.bx;
      if((index_j.y>=0)&&(index_l.y>=0)) HessianMatrix.element[MIN2(index_j.y,index_l.y)][MAX2(index_j.y,index_l.y)]+=DDF*dtB.y*dtD.y+D2JL.by;
      if((index_j.y>=0)&&(index_l.z>=0)) HessianMatrix.element[MIN2(index_j.y,index_l.z)][MAX2(index_j.y,index_l.z)]+=DDF*dtB.y*dtD.z+D2JL.bz;
      if((index_j.z>=0)&&(index_l.x>=0)) HessianMatrix.element[MIN2(index_j.z,index_l.x)][MAX2(index_j.z,index_l.x)]+=DDF*dtB.z*dtD.x+D2JL.cx;
      if((index_j.z>=0)&&(index_l.y>=0)) HessianMatrix.element[MIN2(index_j.z,index_l.y)][MAX2(index_j.z,index_l.y)]+=DDF*dtB.z*dtD.y+D2JL.cy;
      if((index_j.z>=0)&&(index_l.z>=0)) HessianMatrix.element[MIN2(index_j.z,index_l.z)][MAX2(index_j.z,index_l.z)]+=DDF*dtB.z*dtD.z+D2JL.cz;

      if((index_k.x>=0)&&(index_l.x>=0)) HessianMatrix.element[MIN2(index_k.x,index_l.x)][MAX2(index_k.x,index_l.x)]+=DDF*dtC.x*dtD.x+D2KL.ax;
      if((index_k.x>=0)&&(index_l.y>=0)) HessianMatrix.element[MIN2(index_k.x,index_l.y)][MAX2(index_k.x,index_l.y)]+=DDF*dtC.x*dtD.y+D2KL.ay;
      if((index_k.x>=0)&&(index_l.z>=0)) HessianMatrix.element[MIN2(index_k.x,index_l.z)][MAX2(index_k.x,index_l.z)]+=DDF*dtC.x*dtD.z+D2KL.az;
      if((index_k.y>=0)&&(index_l.x>=0)) HessianMatrix.element[MIN2(index_k.y,index_l.x)][MAX2(index_k.y,index_l.x)]+=DDF*dtC.y*dtD.x+D2KL.bx;
      if((index_k.y>=0)&&(index_l.y>=0)) HessianMatrix.element[MIN2(index_k.y,index_l.y)][MAX2(index_k.y,index_l.y)]+=DDF*dtC.y*dtD.y+D2KL.by;
      if((index_k.y>=0)&&(index_l.z>=0)) HessianMatrix.element[MIN2(index_k.y,index_l.z)][MAX2(index_k.y,index_l.z)]+=DDF*dtC.y*dtD.z+D2KL.bz;
      if((index_k.z>=0)&&(index_l.x>=0)) HessianMatrix.element[MIN2(index_k.z,index_l.x)][MAX2(index_k.z,index_l.x)]+=DDF*dtC.z*dtD.x+D2KL.cx;
      if((index_k.z>=0)&&(index_l.y>=0)) HessianMatrix.element[MIN2(index_k.z,index_l.y)][MAX2(index_k.z,index_l.y)]+=DDF*dtC.z*dtD.y+D2KL.cy;
      if((index_k.z>=0)&&(index_l.z>=0)) HessianMatrix.element[MIN2(index_k.z,index_l.z)][MAX2(index_k.z,index_l.z)]+=DDF*dtC.z*dtD.z+D2KL.cz;

      HessianTorsionStrainPosition(HessianMatrix,index_i,index_j,index_k,index_l,
             vec_u,vec_v,vec_w,u,v,w,fa,fc,fd,Dab,Dcb,Ddc,rbc,D2I,D2K,D2L,
             D2IJ,D2IK,D2IL,D2JK,D2JL,D2KL,dtA,dtB,dtC,dtD,DDF,S,CosPhi);

      HessianTorsionStrainStrain(HessianMatrix,vec_u,vec_v,vec_w,u,v,w,DF,DDF,S,CosPhi);
    }
  }
}

