/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'vector.c' is part of RASPA-2.0

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
#include "vector.h"
#include "utils.h"

INT_VECTOR3 UNDEFINED_INT_VECTOR3={-1,-1,-1};
INT_VECTOR3 ZERO_INT_VECTOR3={0,0,0};
INT_VECTOR3 UNIT_INT_VECTOR3={1,1,1};

INT_VECTOR CreateIntVector(int m)
{
  INT_VECTOR c;

  c.m=m;
  c.element=(int*)calloc(m,sizeof(int));
  return c;
}

void DeleteIntVector(INT_VECTOR m)
{
  free(m.element);
}

REAL_VECTOR CreateRealVector(int m)
{
  REAL_VECTOR c;

  c.m=m;
  c.element=(REAL*)calloc(m,sizeof(REAL));
  return c;
}

void DeleteRealVector(REAL_VECTOR m)
{
  free(m.element);
}

COMPLEX_VECTOR CreateComplexVector(int m)
{
  COMPLEX_VECTOR c;

  c.m=m;
  c.element=(Complex*)calloc(m,sizeof(Complex));
  return c;
}

void DeleteComplexVector(COMPLEX_VECTOR m)
{
  free(m.element);
}


POINT_VECTOR CreatePointVector(int m)
{
  POINT_VECTOR c;

  c.m=m;
  c.element=(POINT*)calloc(m,sizeof(POINT));
  return c;
}

void DeletePointVector(POINT_VECTOR m)
{
  free(m.element);
}

REAL DotProduct(VECTOR a,VECTOR b)
{
  return (a.x*b.x+a.y*b.y+a.z*b.z);
}

VECTOR CrossProduct(VECTOR a,VECTOR b)
{
  VECTOR c;

  c.x=a.y*b.z-a.z*b.y;
  c.y=a.z*b.x-a.x*b.z;
  c.z=a.x*b.y-a.y*b.x;
  return c;
}

VECTOR NormalizeVector(VECTOR a)
{
  REAL length;
  VECTOR c;

  length=sqrt(SQR(a.x)+SQR(a.y)+SQR(a.z));
  c.x=a.x/length;
  c.y=a.y/length;
  c.z=a.z/length;
  return c;
}

VECTOR Negative(VECTOR a)
{
  VECTOR c;

  c.x=-a.x;
  c.y=-a.y;
  c.z=-a.z;
  return c;
}
