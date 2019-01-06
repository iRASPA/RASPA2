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

#ifndef COMPLEX_H
#define COMPLEX_H

#include "constants.h"

typedef struct Complex
{
  REAL re;
  REAL im;
} COMPLEX;
typedef struct Complex Complex;

REAL Re(Complex a);
REAL Im(Complex a);

Complex MakeComplex(REAL re,REAL im);         // create complex number
Complex ComplexAdd(Complex a,Complex b);                    // complex addition
Complex ComplexSubtraction(Complex a,Complex b);            // complex subtraction
Complex ComplexMultiplication(Complex a,Complex b);         // complex multiplication
Complex ComplexDivision(Complex a,Complex b);               // complex division
Complex ComplexRealMuliplication(REAL x,Complex a);  // complex-real multiplication
Complex Conjugate(Complex z);                               // complex conjugate
REAL ComplexAbs(Complex z);                          // absolute value of a complex number
Complex ComplexSqrt(Complex z);                             // squareroot of a complex number
Complex ComplexLog(Complex a);


#endif
