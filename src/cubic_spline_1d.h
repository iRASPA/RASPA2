/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'cubic_spline_1d.h' is part of RASPA-2.0

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

#ifndef CUBIC_SPLINE_1D_H
#define CUBIC_SPLINE_1D_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "constants.h"
#include "complex.h"

typedef short DEGREE;
typedef int INDEX;

// a general point in three-dimensional space
typedef struct point
{
  REAL x;
  REAL y;
  REAL z;
} POINT,VECTOR;

typedef struct complex_point
{
  COMPLEX x;
  COMPLEX y;
  COMPLEX z;
} CVECTOR;


// a control point
typedef struct cpoint
{
  REAL x;
  REAL y;
  REAL z;
} CPOINT;

// a polygon
typedef struct polygon
{
  int n;
  POINT *Pw;
} POLYGON;

// a control polygon
typedef struct cpolygon
{
  int n;
  CPOINT *Pw;
} CPOLYGON;

// a knotvector
typedef struct knotvector
{
  int m;
  REAL *U;
} KNOTVECTOR;

// a curve
typedef struct curve
{
  CPOLYGON pol;
  DEGREE p;
  KNOTVECTOR knt;
} BSPLINE_CURVE;

// a control net
typedef struct cnet
{
  INDEX m;
  INDEX n;
  CPOINT **Pw;
} CNET;

// a surface
typedef struct surface
{
  CNET net;
  DEGREE p;
  DEGREE q;
  KNOTVECTOR knu;
  KNOTVECTOR knv;
} BSPLINE_SURFACE;

// a line
typedef struct line
{
  POINT p1;
  POINT p2;
  VECTOR v;
} LINE;

// a plane
typedef struct plane
{
  POINT p;
  VECTOR n;
} PLANE;

// a circle
typedef struct circle
{
  POINT c;
  REAL r;
  VECTOR n;
} CIRCLE;


enum{NOT_A_NODE,FIRST_DERIVATIVE,SECOND_DERIVATIVE,THIRD_DERIVATIVE,PERIODIC_SPLINE};

typedef struct cubic_spline
{
  int n;

  REAL *x;
  REAL *y;

  REAL *a;
  REAL *b;
  REAL *c;
  REAL *d;
} CUBIC_SPLINE;

int Interval(int n,REAL xwert,REAL *x);
REAL EvaluateCubicSplineFast(CUBIC_SPLINE spline,REAL r);
REAL EvaluateCubicSpline(CUBIC_SPLINE spline,int n,REAL r,
              REAL *x,REAL *derivative);
void DeleteCubicSpline(CUBIC_SPLINE spline);
REAL IntegrateCubicSpline(CUBIC_SPLINE spline,int n,
                 REAL r1,REAL r2,REAL *x);
REAL QuadratureSimpson(CUBIC_SPLINE spline,REAL f(REAL),REAL a,REAL b);
REAL QuadratureSimpsonWeight(CUBIC_SPLINE spline,REAL f(REAL),
                REAL g(REAL),REAL a,REAL b);

CUBIC_SPLINE CreateCubicSpline(int m,REAL *x,REAL *y,int boundary_condition,
                      REAL left_boundary,REAL right_boundary);

CUBIC_SPLINE CreateCubicFittingSpline(int n,REAL *xn,REAL *fn,REAL *w,
              int marg_cond,REAL marg_0,REAL marg_n);

KNOTVECTOR CreateKnotVector(int m);
void DeleteKnotVector(KNOTVECTOR c);
void PrintKnotVector(KNOTVECTOR c);


CNET CreateControlNet(int m,int n);
void DeleteControlNet(CNET m);

REAL Distance3D(POINT a,POINT b);

REAL ShortestDistancePointToLine(POINT point,LINE line);
POINT PointToLine(POINT point,LINE line);

#endif
