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
