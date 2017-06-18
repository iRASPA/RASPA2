/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'complex.h' is part of RASPA-2.0

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
