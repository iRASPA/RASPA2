/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'utils.h' is part of RASPA-2.0

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

#ifndef UTILS_H
#define UTILS_H

#include "constants.h"
#include "complex.h"
#include "matrix.h"

enum {WINDOW,BIG_CAGE,SMALL_CAGE};
enum {ABSOLUTE,RELATIVE};
enum {ANISOTROPIC_BISECTION,ANISOTROPIC_MID_POINT};

enum{BOX,CYLINDER};

typedef struct quaternion
{
  REAL r;
  REAL i;
  REAL j;
  REAL k;
} QUATERNION,VECTOR4;

#define TRUE 1
#define FALSE 0

#define UNDEFINED -99999

#define NINT(x) ((int)((x)>=0.0?((x)+0.5):((x)-0.5)) )
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MIN2(x,y) (((x)<=(y))?(x):(y))                 // the minimum of two numbers
#define MAX2(x,y) (((x)>(y))?(x):(y))                 // the maximum of two numbers
#define SIGN(a,b) ((b)>=0.0?fabs(a):-fabs(a))    // a with sign of b
#define RANGE(x,a,b,eps) (((x)>=(a-eps))&&((x)<=(b+eps)))
#define RANGE2LR(x,a,b,eps) (((x)>=(a-eps))&&((x)<=(b+eps)))
#define RANGE2L(x,a,b,eps) (((x)>=(a-eps))&&((x)<(b+eps)))
#define RANGE2R(x,a,b,eps) (((x)>(a-eps))&&((x)<=(b+eps)))
#define RANGE2(a,x,b) (((a)<(x))&&((x)<(b)))
#define MIN3(x,y,z) MIN2((x),MIN2((y),(z)))
#define MAX3(x,y,z) MAX2((x),MAX2((y),(z)))
#define MIN4(u,v,w,x) MIN2(MIN2((u),(v)),MIN2((w),(x)))
#define MAX4(u,v,w,x) MAX2(MAX2((u),(v)),MAX2((w),(x)))
#define SWAP(x,y,z) {z=(x);x=(y);y=(z);}

#define ZERO_TEST (1e-8)

#ifndef M_1_SQRTPI
#define M_1_SQRTPI 0.56418958354775628694807945156077259
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.7071067811865475244008443621048490
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#ifndef M_SQRT1_3
#  define M_SQRT1_3   0.57735026918962576451      /* 1/sqrt(3) */
#endif

#define BIT(x) (1 << (x))
#define SETBITS(x,y) ((x) |= (y))
#define CLEARBITS(x,y) ((x) &= (~(y)))
#define SETBIT(x,y) SETBITS((x), (BIT((y))))
#define CLEARBIT(x,y) CLEARBITS((x), (BIT((y))))
#define BITSET(x,y) ((x) & (BIT(y)))
#define BITCLEAR(x,y) !BITSET((x), (y))
#define BITSSET(x,y) (((x) & (y)) == (y))
#define BITSCLEAR(x,y) (((x) & (y)) == 0)
#define BITVAL(x,y) (((x)>>(y)) & 1)

long Factorial(int n);
int isInArrayOfSize(int value,int n,int *array);
void BubbleSort(int list[], int n);
double get_wall_time(void);
double get_cpu_time(void);

REAL Smoothing(REAL theta);
REAL SmoothingDerivative(REAL theta);

REAL Exp(REAL x);
REAL MExp(REAL x);
REAL Exp2(REAL x);
REAL MExp2(REAL x);
REAL Identity(REAL x);
REAL ExpRM(REAL x);

VECTOR RotateAboutArbitraryLine(VECTOR vec, VECTOR origin, VECTOR dir,REAL theta);
VECTOR RotateAboutBondVector(VECTOR p, VECTOR first_bead, VECTOR second_bead, REAL theta);

void ConvertStringToLowercase(char *buffer);
void ConvertStringToUppercase(char *buffer);

#if defined (__LP64__) || defined (__64BIT__) || defined (_LP64) || (__WORDSIZE == 64)
  void InitializeRandomNumberGenerator(unsigned long long seed);
#else
  void InitializeRandomNumberGenerator(unsigned long seed);
#endif
REAL RandomNumber(void);                                // returns a random number
REAL RandomGaussianNumber(void);                        // returns a Gaussian number
REAL RandomGaussianNumberMeanVariance(REAL mean,REAL sigma);

void RetrieveRandomNumberStatus(void);
void StoreRandomNumberStatus(void);

VECTOR RandomNumberOnUnitSphere(void);
VECTOR RandomNumberOnCone(VECTOR v,REAL theta);
VECTOR Perpendicular(VECTOR a,VECTOR b);
void CalculateRotationMatrix(VECTOR a,VECTOR b,REAL_MATRIX3x3 *rot);
void RotationAroundXYZAxis(VECTOR v,VECTOR *Cord,int n,REAL theta);
void RotationAroundXAxis(VECTOR *Cord,int n,REAL theta);
void RotationAroundYAxis(VECTOR *Cord,int n,REAL theta);
void RotationAroundZAxis(VECTOR *Cord,int n,REAL theta);
void RandomArrayRotationMatrix(VECTOR *Cord,int n);
VECTOR RotationAroundXYZAxisAtPoint(POINT origin,VECTOR v,POINT Cord,REAL theta);
void RotationArrayAroundXYZAxisAtPoint(POINT origin,VECTOR v,POINT *Cord,int n,REAL theta);

REAL DotProdukt(VECTOR a,VECTOR b);
VECTOR RotateVectorAboutX(VECTOR vec,REAL angle);
VECTOR RotateVectorAboutY(VECTOR vec,REAL angle);
VECTOR RotateVectorAboutZ(VECTOR vec,REAL angle);

VECTOR GenerateRandomCylinderPoint(int dir,VECTOR origin,VECTOR rotation,REAL radius);
VECTOR TransformMapping(REAL_MATRIX3x3 m,VECTOR t);

QUATERNION MultiplyQuarternions(QUATERNION q1,QUATERNION q2);
QUATERNION RandomQuarternion(void);
void EulerToQuaternions(QUATERNION *q,VECTOR Ang);
QUATERNION ScaleQuarternion(QUATERNION q1,REAL s);

int quintic(REAL [], REAL [], REAL [], int*, REAL);
int quartic(REAL[], REAL[], REAL[], int* );
int cubic(REAL[], REAL[], int*);
int signR(REAL);
REAL CBRT(REAL);

void LineSearch(int np,int nb,int n,REAL *xold,REAL fold,REAL *g,REAL *p,REAL *x,
                REAL *f,REAL stpmax,int *check,REAL (*func)(int,int,REAL []));

REAL mnbrak(int np,int nb,REAL *ax, REAL *bx, REAL *cx, REAL *fa, REAL *fb, REAL *fc,
        REAL (*func)(int,int,REAL));

REAL dbrent(int np,int nb,REAL ax, REAL bx, REAL cx, REAL (*f)(int,int,REAL),
        REAL (*df)(int,int,REAL), REAL tol, REAL *xmin);

void dlinmin(int np,int nb,REAL *p, REAL *xi, int n, REAL *fret, REAL (*func)(int,int,REAL []),
        void (*dfunc)(int,int,REAL [], REAL []));

void linmin(int np,int nb,REAL *p, REAL *xi, int n, REAL *fret, REAL (*func)(int,int,REAL []));

void powell(int np,int nb,REAL *p, REAL **xi, int n, REAL ftol, int *iter, REAL *fret,
        REAL (*func)(int,int,REAL []));

void FastFourierTransform(REAL(* data)[2], unsigned long nn, int isign);

REAL_MATRIX3x3 ConvertToVoigt2D(REAL_MATRIX9x9 b);
REAL_MATRIX6x6 ConvertToVoigt3D(REAL_MATRIX9x9 b);

int CheckEnantioFace(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR posE);

void WriteRestartUtils(FILE *FilePtr);
void ReadRestartUtils(FILE *FilePtr);

void TrimStringInPlace(char *s);
void CompressSpacesInString(char *str);
void StripLeadingAndTrailingQuotesInPlace(char *str);
void StripWhiteSpacesInPlace(char *str);
void ReplaceCharacterInString(char * string,char searchchar,char replacechar);
void RemoveWhiteSpacesFromString(char *string);
char* StringReverse(char *s);

#endif
