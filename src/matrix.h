/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'matrix.h' is part of RASPA-2.0

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

#ifndef REAL_MATRIX_H
#define REAL_MATRIX_H

#include <stdio.h>
#include "complex.h"
#include "cubic_spline_1d.h"
#include "constants.h"

typedef struct real_matrix3x3
{
  REAL ax;
  REAL ay;
  REAL az;

  REAL bx;
  REAL by;
  REAL bz;

  REAL cx;
  REAL cy;
  REAL cz;

} REAL_MATRIX3x3;

typedef struct real_matrix5x5
{
  REAL C11; REAL C21; REAL C31; REAL C41; REAL C51;
  REAL C12; REAL C22; REAL C32; REAL C42; REAL C52;
  REAL C13; REAL C23; REAL C33; REAL C43; REAL C53;
  REAL C14; REAL C24; REAL C34; REAL C44; REAL C54;
  REAL C15; REAL C25; REAL C35; REAL C45; REAL C55;

} REAL_MATRIX5x5;


typedef struct real_matrix6x6
{
  REAL C11; REAL C21; REAL C31; REAL C41; REAL C51; REAL C61;
  REAL C12; REAL C22; REAL C32; REAL C42; REAL C52; REAL C62;
  REAL C13; REAL C23; REAL C33; REAL C43; REAL C53; REAL C63;
  REAL C14; REAL C24; REAL C34; REAL C44; REAL C54; REAL C64;
  REAL C15; REAL C25; REAL C35; REAL C45; REAL C55; REAL C65;
  REAL C16; REAL C26; REAL C36; REAL C46; REAL C56; REAL C66;

} REAL_MATRIX6x6;

typedef struct real_matrix6x6x6
{
  REAL C111; REAL C211; REAL C311; REAL C411; REAL C511; REAL C611;
  REAL C121; REAL C221; REAL C321; REAL C421; REAL C521; REAL C621;
  REAL C131; REAL C231; REAL C331; REAL C431; REAL C531; REAL C631;
  REAL C141; REAL C241; REAL C341; REAL C441; REAL C541; REAL C641;
  REAL C151; REAL C251; REAL C351; REAL C451; REAL C551; REAL C651;
  REAL C161; REAL C261; REAL C361; REAL C461; REAL C561; REAL C661;

  REAL C112; REAL C212; REAL C312; REAL C412; REAL C512; REAL C612;
  REAL C122; REAL C222; REAL C322; REAL C422; REAL C522; REAL C622;
  REAL C132; REAL C232; REAL C332; REAL C432; REAL C532; REAL C632;
  REAL C142; REAL C242; REAL C342; REAL C442; REAL C542; REAL C642;
  REAL C152; REAL C252; REAL C352; REAL C452; REAL C552; REAL C652;
  REAL C162; REAL C262; REAL C362; REAL C462; REAL C562; REAL C662;

  REAL C113; REAL C213; REAL C313; REAL C413; REAL C513; REAL C613;
  REAL C123; REAL C223; REAL C323; REAL C423; REAL C523; REAL C623;
  REAL C133; REAL C233; REAL C333; REAL C433; REAL C533; REAL C633;
  REAL C143; REAL C243; REAL C343; REAL C443; REAL C543; REAL C643;
  REAL C153; REAL C253; REAL C353; REAL C453; REAL C553; REAL C653;
  REAL C163; REAL C263; REAL C363; REAL C463; REAL C563; REAL C663;

  REAL C114; REAL C214; REAL C314; REAL C414; REAL C514; REAL C614;
  REAL C124; REAL C224; REAL C324; REAL C424; REAL C524; REAL C624;
  REAL C134; REAL C234; REAL C334; REAL C434; REAL C534; REAL C634;
  REAL C144; REAL C244; REAL C344; REAL C444; REAL C544; REAL C644;
  REAL C154; REAL C254; REAL C354; REAL C454; REAL C554; REAL C654;
  REAL C164; REAL C264; REAL C364; REAL C464; REAL C564; REAL C664;

  REAL C115; REAL C215; REAL C315; REAL C415; REAL C515; REAL C615;
  REAL C125; REAL C225; REAL C325; REAL C425; REAL C525; REAL C625;
  REAL C135; REAL C235; REAL C335; REAL C435; REAL C535; REAL C635;
  REAL C145; REAL C245; REAL C345; REAL C445; REAL C545; REAL C645;
  REAL C155; REAL C255; REAL C355; REAL C455; REAL C555; REAL C655;
  REAL C165; REAL C265; REAL C365; REAL C465; REAL C565; REAL C665;

  REAL C116; REAL C216; REAL C316; REAL C416; REAL C516; REAL C616;
  REAL C126; REAL C226; REAL C326; REAL C426; REAL C526; REAL C626;
  REAL C136; REAL C236; REAL C336; REAL C436; REAL C536; REAL C636;
  REAL C146; REAL C246; REAL C346; REAL C446; REAL C546; REAL C646;
  REAL C156; REAL C256; REAL C356; REAL C456; REAL C556; REAL C656;
  REAL C166; REAL C266; REAL C366; REAL C466; REAL C566; REAL C666;
} REAL_MATRIX6x6x6;


typedef struct real_matrix9x9
{
  REAL xxxx; REAL yxxx; REAL zxxx;  REAL xxyx; REAL yxyx; REAL zxyx;  REAL xxzx; REAL yxzx; REAL zxzx;
  REAL xyxx; REAL yyxx; REAL zyxx;  REAL xyyx; REAL yyyx; REAL zyyx;  REAL xyzx; REAL yyzx; REAL zyzx;
  REAL xzxx; REAL yzxx; REAL zzxx;  REAL xzyx; REAL yzyx; REAL zzyx;  REAL xzzx; REAL yzzx; REAL zzzx;
  REAL xxxy; REAL yxxy; REAL zxxy;  REAL xxyy; REAL yxyy; REAL zxyy;  REAL xxzy; REAL yxzy; REAL zxzy;
  REAL xyxy; REAL yyxy; REAL zyxy;  REAL xyyy; REAL yyyy; REAL zyyy;  REAL xyzy; REAL yyzy; REAL zyzy;
  REAL xzxy; REAL yzxy; REAL zzxy;  REAL xzyy; REAL yzyy; REAL zzyy;  REAL xzzy; REAL yzzy; REAL zzzy;
  REAL xxxz; REAL yxxz; REAL zxxz;  REAL xxyz; REAL yxyz; REAL zxyz;  REAL xxzz; REAL yxzz; REAL zxzz;
  REAL xyxz; REAL yyxz; REAL zyxz;  REAL xyyz; REAL yyyz; REAL zyyz;  REAL xyzz; REAL yyzz; REAL zyzz;
  REAL xzxz; REAL yzxz; REAL zzxz;  REAL xzyz; REAL yzyz; REAL zzyz;  REAL xzzz; REAL yzzz; REAL zzzz;

} REAL_MATRIX9x9;


typedef struct real_matrix4x4
{
  REAL ar;
  REAL ai;
  REAL aj;
  REAL ak;

  REAL br;
  REAL bi;
  REAL bj;
  REAL bk;

  REAL cr;
  REAL ci;
  REAL cj;
  REAL ck;

  REAL dr;
  REAL di;
  REAL dj;
  REAL dk;
} REAL_MATRIX4x4;


typedef struct int_matrix
{
  int m;
  int n;
  int **element;
} INT_MATRIX;

typedef struct real_matrix
{
  int m;
  int n;
  REAL **element;
} REAL_MATRIX;

typedef struct real_fortran_matrix
{
  int m;
  int n;
  REAL *element;
} REAL_FORTRAN_MATRIX;


typedef struct complex_matrix
{
  int m;
  int n;
  Complex **element;
} COMPLEX_MATRIX;

typedef struct point_matrix
{
  int m;
  int n;
  POINT **element;
} POINT_MATRIX;

void PrintRealMatrix3x3(REAL_MATRIX3x3 *m);
REAL Trace3x3Matrix(REAL_MATRIX3x3 *in);

INT_MATRIX CreateIntMatrix(int m,int n);
void DeleteIntMatrix(INT_MATRIX c);

void InitializeMatrix6x6(REAL_MATRIX6x6 *m);
void InitializeMatrix6x6x6(REAL_MATRIX6x6x6 *m);
void InitializeMatrix9x9(REAL_MATRIX9x9 *m);
void PrintRealMatrix9x9(REAL_MATRIX9x9 *m);

REAL_MATRIX CreateRealMatrix(int m,int n);
void DeleteRealMatrix(REAL_MATRIX c);
void PrintRealMatrix(REAL_MATRIX *c);
void PrintComplexMatrix(COMPLEX_MATRIX *c);

REAL_FORTRAN_MATRIX CreateRealFortranMatrix(int m,int n);
void DeleteRealFortranMatrix(REAL_FORTRAN_MATRIX c);
void PrintRealFortranMatrix(REAL_FORTRAN_MATRIX *c);

COMPLEX_MATRIX CreateComplexMatrix(int m,int n);
void DeleteComplexMatrix(COMPLEX_MATRIX c);

POINT_MATRIX CreatePointMatrix(int m,int n);
void DeletePointMatrix(POINT_MATRIX c);

void AddRealMatrix3x3(REAL_MATRIX3x3 *c,REAL_MATRIX3x3 a,REAL_MATRIX3x3 b);
void SubtractRealMatrix3x3(REAL_MATRIX3x3 *c,REAL_MATRIX3x3 a,REAL_MATRIX3x3 b);
void DivideRealMatrix3x3ByReal(REAL_MATRIX3x3 *c,REAL_MATRIX3x3 a,REAL b);

void AddRealMatrix6x6(REAL_MATRIX6x6 *c,REAL_MATRIX6x6 a,REAL_MATRIX6x6 b);
void SubtractRealMatrix6x6(REAL_MATRIX6x6 *c,REAL_MATRIX6x6 a,REAL_MATRIX6x6 b);
void DivideRealMatrix6x6ByReal(REAL_MATRIX6x6 *c,REAL_MATRIX6x6 a,REAL b);

void AddRealMatrix9x9(REAL_MATRIX9x9 *c,REAL_MATRIX9x9 a,REAL_MATRIX9x9 b);
void SubtractRealMatrix9x9(REAL_MATRIX9x9 *c,REAL_MATRIX9x9 a,REAL_MATRIX9x9 b);
void DivideRealMatrix9x9ByReal(REAL_MATRIX9x9 *c,REAL_MATRIX9x9 a,REAL b);

void EigenSystemVoigt3x3(REAL_MATRIX3x3 in,REAL *eigenvalues);
void InverseMatrix3x3(REAL_MATRIX3x3 in,REAL_MATRIX3x3 *out);
void EigenSystemVoigt6x6(REAL_MATRIX6x6 in,REAL *eigenvalues);
void InverseMatrix6x6(REAL_MATRIX6x6 in,REAL_MATRIX6x6 *out);

REAL_MATRIX3x3 MatrixMatrixMultiplication3x3(REAL_MATRIX3x3 a,REAL_MATRIX3x3 b);
VECTOR MatrixVectorMultiplication3x3(REAL_MATRIX3x3 a,VECTOR b);
void RescaleMatrix3x3x(REAL_MATRIX3x3 *a,REAL c);
void Invert3x3Matrix(REAL_MATRIX3x3 *in,REAL_MATRIX3x3 *out,REAL *determinant);
REAL GetLargestElementMatrix3x3(REAL_MATRIX3x3 *a);
REAL RescaleMatrix3x3(REAL_MATRIX3x3 *a,REAL c);
void EigenSystem3x3(REAL_MATRIX3x3 in,REAL_MATRIX3x3 *eigenvectors,VECTOR *eigenvalues);
void TransposeMatrix3x3(REAL_MATRIX3x3 *in);

void tqli(REAL d[], REAL e[], int n, REAL *z);
void tred2(REAL *a, int n, REAL d[], REAL e[]);
void eigsrt(REAL d[], REAL *v, int n);
void eigsrt3(REAL d[], REAL *v, int n);
void jacobi(REAL *a, int n, REAL d[], REAL *v, int *nrot);
void ModifiedGramSchmidt(REAL *A,int m,int n);

void CheckMatrixInversion(void);
void InverseRealMatrix(REAL_MATRIX a);
void InverseComplexMatrix(COMPLEX_MATRIX a);
void SingularValueDecomposition(REAL_MATRIX a,int m,int n,REAL w[],REAL_MATRIX v);
void SingularValueDecompositionMatrixInversion(REAL_MATRIX a);
void PrintRealMatrixmMathematica(REAL_MATRIX *c);
void PrintRealMatrix3x3ToFile(REAL_MATRIX3x3 *m,FILE *FilePtr,REAL a);
void PrintRealMatrix6x6ToFile(REAL_MATRIX6x6 *m,FILE *FilePtr,REAL a);
void PrintRealMatrix6x6x6ToFile(REAL_MATRIX6x6x6 *m,FILE *FilePtr,REAL a);
void PrintRealMatrix9x9ToFile(REAL_MATRIX9x9 *m,FILE *FilePtr,REAL a);
void MultiplyRealMatrix(REAL_MATRIX c,REAL_MATRIX a,REAL_MATRIX b);
void MultiplyComplexMatrix(COMPLEX_MATRIX c,COMPLEX_MATRIX a,COMPLEX_MATRIX b);
void MultiplyRealMatrixVector(REAL *c,REAL_MATRIX a,REAL *b);
void TransposeRealMatrix(REAL_MATRIX c,REAL_MATRIX a);
void ScaleRealMatrix(REAL_MATRIX a,REAL fac);
void PrintRealMatrixToFile(REAL_MATRIX *a,char *name);
void Convert9x9ToRealMatrix(REAL_MATRIX9x9 *m,REAL_MATRIX c);
void TestEigenStystem(REAL_MATRIX H);

void PolynomialFit(double *x,double *y,double *std_dev,unsigned int ndata,int num_terms);

#endif
