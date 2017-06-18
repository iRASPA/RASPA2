/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'vector.h' is part of RASPA-2.0

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

#ifndef VECTOR_H
#define VECTOR_H

#include "complex.h"
#include "cubic_spline_1d.h"
#include "constants.h"

typedef struct int_vector3
{
  int x;
  int y;
  int z;
} INT_VECTOR3;

extern INT_VECTOR3 UNDEFINED_INT_VECTOR3;
extern INT_VECTOR3 ZERO_INT_VECTOR3;
extern INT_VECTOR3 UNIT_INT_VECTOR3;

typedef struct int_vector
{
  int m;
  int *element;
} INT_VECTOR;

typedef struct real_vector
{
  int m;
  REAL *element;
} REAL_VECTOR;

typedef struct complex_vector
{
  int m;
  Complex *element;
} COMPLEX_VECTOR;

typedef struct point_vector
{
  int m;
  POINT *element;
} POINT_VECTOR;

INT_VECTOR CreateIntVector(int m);
void DeleteIntVector(INT_VECTOR c);

REAL_VECTOR CreateRealVector(int m);
void DeleteRealVector(REAL_VECTOR c);

COMPLEX_VECTOR CreateComplexVector(int m);
void DeleteComplexVector(COMPLEX_VECTOR c);

POINT_VECTOR CreatePointVector(int m);
void DeletePointVector(POINT_VECTOR c);

REAL DotProduct(VECTOR a,VECTOR b);
VECTOR CrossProduct(VECTOR a,VECTOR b);
VECTOR NormalizeVector(VECTOR a);
VECTOR Negative(VECTOR a);
#endif
