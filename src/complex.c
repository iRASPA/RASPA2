/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'complex.c' is part of RASPA-2.0

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
#include "complex.h"

REAL Re(Complex a)
{
  return a.re;
}

REAL Im(Complex a)
{
  return a.im;
}

// creates a Complex number from the real and imaginary part
Complex MakeComplex(REAL re,REAL im)
{
  Complex c;
  c.re=re;
  c.im=im;
  return c;
}

// returns the Complex addition of Complex a and Complex b
Complex ComplexAdd(Complex a,Complex b)
{
  Complex c;

  c.re=a.re+b.re;
  c.im=a.im+b.im;
  return c;
}

// returns the Complex subtraction of Complex a and Complex b
Complex ComplexSubtraction(Complex a,Complex b)
{
  Complex c;
  c.re=a.re-b.re;
  c.im=a.im-b.im;
  return c;
}

// returns the Complex multiplication of Complex a and Complex b
Complex ComplexMultiplication(Complex a,Complex b)
{
  Complex c;
  c.re=a.re*b.re-a.im*b.im;
  c.im=a.im*b.re+a.re*b.im;
  return c;
}

// returns the Complex division of Complex a and Complex b (without under/overflow)
Complex ComplexDivision(Complex a,Complex b)
{
  Complex c;
  REAL r,den;
  if (fabs(b.re)>=fabs(b.im))
  {
    r=b.im/b.re;
    den=b.re+r*b.im;
    c.re=(a.re+r*a.im)/den;
    c.im=(a.im-r*a.re)/den;
  }
  else
  {
    r=b.re/b.im;
    den=b.im+r*b.re;
    c.re=(a.re*r+a.im)/den;
    c.im=(a.im*r-a.re)/den;
  }
  return c;
}

// returns the multiplication of real x with Complex a
Complex ComplexRealMuliplication(REAL x,Complex a)
{
  Complex c;
  c.re=x*a.re;
  c.im=x*a.im;
  return c;
}
// returns the conjugate of Complex z
Complex Conjugate(Complex z)
{
  Complex c;
  c.re=z.re;
  c.im=-z.im;
  return c;
}

// returns the absolute value of Complex z (without under/overflow)
REAL ComplexAbs(Complex z)
{
  REAL x,y,ans,temp;
  x=z.re>=0.0?z.re:-z.re;
  y=z.im>=0.0?z.im:-z.im;
  if(x==0.0) ans=y;
  else if(y==0.0) ans=x;
  else if(x>y)
  {
    temp=y/x;
    ans=x*sqrt(1.0+temp*temp);
  }
  else
  {
    temp=x/y;
    ans=y*sqrt(1.0+temp*temp);
  }
  return ans;
}

// returns the squareroot of Complex z (without under/overflow)
Complex ComplexSqrt(Complex z)
{
  Complex c;
  REAL x,y,w,r;
  if ((z.re==0.0)&&(z.im==0.0))
  {
    c.re=0.0;
    c.im=0.0;
    return c;
  }
  else
  {
    x=z.re>=0.0?z.re:-z.re;
    y=z.im>=0.0?z.im:-z.im;
    if (x >= y)
    {
      r=y/x;
      w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
    }
    else
    {
      r=x/y;
      w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
    }
    if (z.re>=0.0)
    {
      c.re=w;
      c.im=z.im/(2.0*w);
    }
    else
    {
      c.im=(z.im>=0.0)?w:-w;
      c.re=z.im/(2.0*c.im);
    }
    return c;
  }
}

Complex ComplexLog(Complex a)
{
  Complex c;

  c.re=log(sqrt(a.re*a.re+a.im*a.im));
  c.im=atan2(a.im,a.re);
  return c;
}
