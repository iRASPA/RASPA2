/*************************************************************************************************************
 The MIT License

 Copyright (C) 2006-2019 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

     D.Dubbeldam@uva.nl            http://www.uva.nl/profiel/d/u/d.dubbeldam/d.dubbeldam.html
     scaldia@upo.es                http://www.upo.es/raspa/
     t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
     don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
     snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************************/

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
