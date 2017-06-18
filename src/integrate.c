/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'integrate.c' is part of RASPA-2.0

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
#include "vector.h"
#include "integrate.h"
#include "cubic_spline_1d.h"

/*
#define EPS 1.0e-8
#define JMAX 20

#define FUNC(x) ((*func)(x))
REAL Trapezoidal(CUBIC_SPLINE spline,REAL f(REAL),REAL a,REAL b,int n)
{
  REAL x,tnm,sum,del,r;
  static REAL s;
  int it,j;
  static int i=0;
  REAL valuea,valueb;
  REAL tmp1,tmp2,tmp3;

  fprintf(stderr, "entering\n");
  if(n==1)
  {
    r=a;
    if (r<spline.x[i]||r>=spline.x[i+1])
      i=Interval(spline.n,r,spline.x);
    r-=spline.x[i];
    tmp1=3.0*spline.d[i];
    tmp2=2.0*spline.c[i];
    tmp3=2.0*tmp1;
    valuea=((spline.d[i]*r+spline.c[i])*r+spline.b[i])*r+spline.a[i];

    r=b;
    if (r<spline.x[i]||r>=spline.x[i+1])
      i=Interval(spline.n,r,spline.x);
    r-=spline.x[i];
    tmp1=3.0*spline.d[i];
    tmp2=2.0*spline.c[i];
    tmp3=2.0*tmp1;
    valueb=((spline.d[i]*r+spline.c[i])*r+spline.b[i])*r+spline.a[i];
    fprintf(stderr, "leaving\n");
    return (s=0.5*(b-a)*(f(valuea)+f(valueb)));
  }
  else
  {
    for (it=1,j=1;j<n-1;j++)
      it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del)
    {
      r=x;
      if (r<spline.x[i]||r>=spline.x[i+1])
        i=Interval(spline.n,r,spline.x);
      r-=spline.x[i];
      tmp1=3.0*spline.d[i];
      tmp2=2.0*spline.c[i];
      tmp3=2.0*tmp1;
      valuea=((spline.d[i]*r+spline.c[i])*r+spline.b[i])*r+spline.a[i];
      sum+=f(valuea);
    }
    s=0.5*(s+(b-a)*sum/tnm);
    fprintf(stderr, "leaving\n");
    return s;
  }
}

REAL QuadratureSimpson(CUBIC_SPLINE spline,REAL f(REAL),REAL a,REAL b)
{
  int j;
  REAL s,st,ost,os;

  ost = os = -1.0e30;
  for (j=1;j<=JMAX;j++)
  {
    st=Trapezoidal(spline,f,a,b,j);
    s=(4.0*st-ost)/3.0;
    if(fabs(s-os) < EPS*fabs(os)) return s;
    if((s==0)&&(os==0)&&(j>6)) return s;
    os=s;
    ost=st;
  }
  fprintf(stderr, "Too many steps in routine qsimp");
  return 0.0;
}
#undef EPS
#undef JMAX

*/


// Simpson's rule for abitrary number of intervals (even or odd)
// Jack W. Hollingsworth and Henry F. Hunter
// Rensselaer Polytechnic Institute, Troy, N.Y.
// A generalization of the Simpson's rule in the sense that it is exact for cubic polynomials and is valid
// for an odd as well as even number of intervals. The weights are equal except for the three points at the end.
// For an even number of intervaks of length 'h', the error coefficient is 23/3 h^5, corresponding to 20/3 h^5 for Simpson's rule.

REAL IntegrateGeneralizedSimpsonOrderN(REAL *data,REAL *count,REAL delta,int n)
{
  int i;
  REAL value;

  switch(n)
  {
    case 0:
      return 0.0;
    case 1:
      return delta*data[0]/count[0];
    case 2:
      return (delta/24.0)*(12.0*data[0]/count[0]+12.0*data[1]/count[1]);
    case 3:
      return (delta/24.0)*(8.0*data[0]/count[0]+32.0*data[1]/count[1]+8.0*data[2]/count[2]);
    case 4:
      return (delta/24.0)*(9.0*data[0]/count[0]+27.0*data[1]/count[1]+27.0*data[2]/count[2]+9.0*data[3]/count[3]);
    case 5:
      return (delta/24.0)*(9.0*data[0]/count[0]+28.0*data[1]/count[1]+22.0*data[2]/count[2]+28.0*data[3]/count[3]+9.0*data[4]/count[4]);
    case 6:
      return (delta/24.0)*(9.0*data[0]/count[0]+28.0*data[1]/count[1]+23.0*data[2]/count[2]+23.0*data[3]/count[3]+28.0*data[4]/count[4]+9.0*data[5]/count[5]);
   default:
      value=9.0*data[0]/count[0]+28.0*data[1]/count[1]+23.0*data[2]/count[2]+23.0*data[n-3]/count[n-3]+28.0*data[n-2]/count[n-2]+9.0*data[n-1]/count[n-1];
      for(i=3;i<n-3;i++)
        value+=24.0*data[i]/count[i];
      return (delta/24.0)*value;
  }
  return 0.0;
}

REAL IntegrateGeneralizedSimpsonConventional(REAL *data,REAL delta,int n)
{
  int i;
  REAL value;

  switch(n)
  {
    case 0:
      return 0.0;
    case 1:
      return delta*data[0];
    case 2:
      return (delta/24.0)*(12.0*data[0]+12.0*data[1]);
    case 3:
      return (delta/24.0)*(8.0*data[0]+32.0*data[1]+8.0*data[2]);
    case 4:
      return (delta/24.0)*(9.0*data[0]+27.0*data[1]+27.0*data[2]+9.0*data[3]);
    case 5:
      return (delta/24.0)*(9.0*data[0]+28.0*data[1]+22.0*data[2]+28.0*data[3]+9.0*data[4]);
    case 6:
      return (delta/24.0)*(9.0*data[0]+28.0*data[1]+23.0*data[2]+23.0*data[3]+28.0*data[4]+9.0*data[5]);
   default:
      value=9.0*data[0]+28.0*data[1]+23.0*data[2]+23.0*data[n-3]+28.0*data[n-2]+9.0*data[n-1];
      for(i=3;i<n-3;i++)
        value+=24.0*data[i];
      return (delta/24.0)*value;
  }
  return 0.0;
}

REAL IntegrateGeneralizedSimpsonConventionalVectorX(VECTOR *data,REAL delta,int n)
{
  int i;
  REAL value;

  switch(n)
  {
    case 0:
      return 0.0;
    case 1:
      return delta*data[0].x;
    case 2:
      return (delta/24.0)*(12.0*data[0].x+12.0*data[1].x);
    case 3:
      return (delta/24.0)*(8.0*data[0].x+32.0*data[1].x+8.0*data[2].x);
    case 4:
      return (delta/24.0)*(9.0*data[0].x+27.0*data[1].x+27.0*data[2].x+9.0*data[3].x);
    case 5:
      return (delta/24.0)*(9.0*data[0].x+28.0*data[1].x+22.0*data[2].x+28.0*data[3].x+9.0*data[4].x);
    case 6:
      return (delta/24.0)*(9.0*data[0].x+28.0*data[1].x+23.0*data[2].x+23.0*data[3].x+28.0*data[4].x+9.0*data[5].x);
   default:
      value=9.0*data[0].x+28.0*data[1].x+23.0*data[2].x+23.0*data[n-3].x+28.0*data[n-2].x+9.0*data[n-1].x;
      for(i=3;i<n-3;i++)
        value+=24.0*data[i].x;
      return (delta/24.0)*value;
  }
  return 0.0;
}

REAL IntegrateGeneralizedSimpsonConventionalVectorY(VECTOR *data,REAL delta,int n)
{
  int i;
  REAL value;

  switch(n)
  {
    case 0:
      return 0.0;
    case 1:
      return delta*data[0].y;
    case 2:
      return (delta/24.0)*(12.0*data[0].y+12.0*data[1].y);
    case 3:
      return (delta/24.0)*(8.0*data[0].y+32.0*data[1].y+8.0*data[2].y);
    case 4:
      return (delta/24.0)*(9.0*data[0].y+27.0*data[1].y+27.0*data[2].y+9.0*data[3].y);
    case 5:
      return (delta/24.0)*(9.0*data[0].y+28.0*data[1].y+22.0*data[2].y+28.0*data[3].y+9.0*data[4].y);
    case 6:
      return (delta/24.0)*(9.0*data[0].y+28.0*data[1].y+23.0*data[2].y+23.0*data[3].y+28.0*data[4].y+9.0*data[5].y);
   default:
      value=9.0*data[0].y+28.0*data[1].y+23.0*data[2].y+23.0*data[n-3].y+28.0*data[n-2].y+9.0*data[n-1].y;
      for(i=3;i<n-3;i++)
        value+=24.0*data[i].y;
      return (delta/24.0)*value;
  }
  return 0.0;
}

REAL IntegrateGeneralizedSimpsonConventionalVectorZ(VECTOR *data,REAL delta,int n)
{
  int i;
  REAL value;

  switch(n)
  {
    case 0:
      return 0.0;
    case 1:
      return delta*data[0].z;
    case 2:
      return (delta/24.0)*(12.0*data[0].z+12.0*data[1].z);
    case 3:
      return (delta/24.0)*(8.0*data[0].z+32.0*data[1].z+8.0*data[2].z);
    case 4:
      return (delta/24.0)*(9.0*data[0].z+27.0*data[1].z+27.0*data[2].z+9.0*data[3].z);
    case 5:
      return (delta/24.0)*(9.0*data[0].z+28.0*data[1].z+22.0*data[2].z+28.0*data[3].z+9.0*data[4].z);
    case 6:
      return (delta/24.0)*(9.0*data[0].z+28.0*data[1].z+23.0*data[2].z+23.0*data[3].z+28.0*data[4].z+9.0*data[5].z);
   default:
      value=9.0*data[0].z+28.0*data[1].z+23.0*data[2].z+23.0*data[n-3].z+28.0*data[n-2].z+9.0*data[n-1].z;
      for(i=3;i<n-3;i++)
        value+=24.0*data[i].z;
      return (delta/24.0)*value;
  }
  return 0.0;
}
