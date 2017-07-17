/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'matrix.c' is part of RASPA-2.0

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
#include "matrix.h"
#include "utils.h"

// Matrices are stored in "1D Fortran-style", so:
//   A_ij, i=0,M-1   j=0,N-1
//   v_k=m+n*N
//
// a 3x2 (m=3,n=2) matrix is stored in memory as
//  0:  0,0
//  1:  1,0
//  2:  2,0
//  3:  0,1
//  4:  0,2
//  5:  0,3
//
// index using the [][] where the second index changes fastest:
//  0:   [0][0]
//  1:   [0][1]
//  2:   [0][2]
//  3:   [1][0]
//  4:   [1][1]
//  5:   [1][2]
//
// the matrix will be printed out as:
//   [0][0]  [1][0]
//   [0][1]  [1][1]
//   [0][2]  [1][2]
//
//  Matrices are stored as the transpose of the usual C-style in order to
//  be compatible with Lapack/Fortran-style.
//    0) allocate the transpose of the matrix
//    1) store the initial elements transposed:
//       Matrix[column][row]=Matrix[n][m]
//    2) use C-style code if indexed as [m][n]
//    3) use fortran-style code if indexed as v_k=m+n*N
//    4) print the transpose of a matrix

void InitializeMatrix6x6(REAL_MATRIX6x6 *m)
{
  m->C11=0.0; m->C21=0.0; m->C31=0.0; m->C41=0.0; m->C51=0.0; m->C61=0.0;
  m->C12=0.0; m->C22=0.0; m->C32=0.0; m->C42=0.0; m->C52=0.0; m->C62=0.0;
  m->C13=0.0; m->C23=0.0; m->C33=0.0; m->C43=0.0; m->C53=0.0; m->C63=0.0;
  m->C14=0.0; m->C24=0.0; m->C34=0.0; m->C44=0.0; m->C54=0.0; m->C64=0.0;
  m->C15=0.0; m->C25=0.0; m->C35=0.0; m->C45=0.0; m->C55=0.0; m->C65=0.0;
  m->C16=0.0; m->C26=0.0; m->C36=0.0; m->C46=0.0; m->C56=0.0; m->C66=0.0;
}

void InitializeMatrix6x6x6(REAL_MATRIX6x6x6 *m)
{
  m->C111=0.0; m->C211=0.0; m->C311=0.0; m->C411=0.0; m->C511=0.0; m->C611=0.0;
  m->C121=0.0; m->C221=0.0; m->C321=0.0; m->C421=0.0; m->C521=0.0; m->C621=0.0;
  m->C131=0.0; m->C231=0.0; m->C331=0.0; m->C431=0.0; m->C531=0.0; m->C631=0.0;
  m->C141=0.0; m->C241=0.0; m->C341=0.0; m->C441=0.0; m->C541=0.0; m->C641=0.0;
  m->C151=0.0; m->C251=0.0; m->C351=0.0; m->C451=0.0; m->C551=0.0; m->C651=0.0;
  m->C161=0.0; m->C261=0.0; m->C361=0.0; m->C461=0.0; m->C561=0.0; m->C661=0.0;

  m->C112=0.0; m->C212=0.0; m->C312=0.0; m->C412=0.0; m->C512=0.0; m->C612=0.0;
  m->C122=0.0; m->C222=0.0; m->C322=0.0; m->C422=0.0; m->C522=0.0; m->C622=0.0;
  m->C132=0.0; m->C232=0.0; m->C332=0.0; m->C432=0.0; m->C532=0.0; m->C632=0.0;
  m->C142=0.0; m->C242=0.0; m->C342=0.0; m->C442=0.0; m->C542=0.0; m->C642=0.0;
  m->C152=0.0; m->C252=0.0; m->C352=0.0; m->C452=0.0; m->C552=0.0; m->C652=0.0;
  m->C162=0.0; m->C262=0.0; m->C362=0.0; m->C462=0.0; m->C562=0.0; m->C662=0.0;

  m->C113=0.0; m->C213=0.0; m->C313=0.0; m->C413=0.0; m->C513=0.0; m->C613=0.0;
  m->C123=0.0; m->C223=0.0; m->C323=0.0; m->C423=0.0; m->C523=0.0; m->C623=0.0;
  m->C133=0.0; m->C233=0.0; m->C333=0.0; m->C433=0.0; m->C533=0.0; m->C633=0.0;
  m->C143=0.0; m->C243=0.0; m->C343=0.0; m->C443=0.0; m->C543=0.0; m->C643=0.0;
  m->C153=0.0; m->C253=0.0; m->C353=0.0; m->C453=0.0; m->C553=0.0; m->C653=0.0;
  m->C163=0.0; m->C263=0.0; m->C363=0.0; m->C463=0.0; m->C563=0.0; m->C663=0.0;

  m->C114=0.0; m->C214=0.0; m->C314=0.0; m->C414=0.0; m->C514=0.0; m->C614=0.0;
  m->C124=0.0; m->C224=0.0; m->C324=0.0; m->C424=0.0; m->C524=0.0; m->C624=0.0;
  m->C134=0.0; m->C234=0.0; m->C334=0.0; m->C434=0.0; m->C534=0.0; m->C634=0.0;
  m->C144=0.0; m->C244=0.0; m->C344=0.0; m->C444=0.0; m->C544=0.0; m->C644=0.0;
  m->C154=0.0; m->C254=0.0; m->C354=0.0; m->C454=0.0; m->C554=0.0; m->C654=0.0;
  m->C164=0.0; m->C264=0.0; m->C364=0.0; m->C464=0.0; m->C564=0.0; m->C664=0.0;

  m->C115=0.0; m->C215=0.0; m->C315=0.0; m->C415=0.0; m->C515=0.0; m->C615=0.0;
  m->C125=0.0; m->C225=0.0; m->C325=0.0; m->C425=0.0; m->C525=0.0; m->C625=0.0;
  m->C135=0.0; m->C235=0.0; m->C335=0.0; m->C435=0.0; m->C535=0.0; m->C635=0.0;
  m->C145=0.0; m->C245=0.0; m->C345=0.0; m->C445=0.0; m->C545=0.0; m->C645=0.0;
  m->C155=0.0; m->C255=0.0; m->C355=0.0; m->C455=0.0; m->C555=0.0; m->C655=0.0;
  m->C165=0.0; m->C265=0.0; m->C365=0.0; m->C465=0.0; m->C565=0.0; m->C665=0.0;

  m->C116=0.0; m->C216=0.0; m->C316=0.0; m->C416=0.0; m->C516=0.0; m->C616=0.0;
  m->C126=0.0; m->C226=0.0; m->C326=0.0; m->C426=0.0; m->C526=0.0; m->C626=0.0;
  m->C136=0.0; m->C236=0.0; m->C336=0.0; m->C436=0.0; m->C536=0.0; m->C636=0.0;
  m->C146=0.0; m->C246=0.0; m->C346=0.0; m->C446=0.0; m->C546=0.0; m->C646=0.0;
  m->C156=0.0; m->C256=0.0; m->C356=0.0; m->C456=0.0; m->C556=0.0; m->C656=0.0;
  m->C166=0.0; m->C266=0.0; m->C366=0.0; m->C466=0.0; m->C566=0.0; m->C666=0.0;
}

void AddRealMatrix3x3(REAL_MATRIX3x3 *c,REAL_MATRIX3x3 a,REAL_MATRIX3x3 b)
{
  c->ax=a.ax+b.ax; c->bx=a.bx+b.bx; c->cx=a.cx+b.cx;
  c->ay=a.ay+b.ay; c->by=a.by+b.by; c->cy=a.cy+b.cy;
  c->az=a.az+b.az; c->bz=a.bz+b.bz; c->cz=a.cz+b.cz;
}

void SubtractRealMatrix3x3(REAL_MATRIX3x3 *c,REAL_MATRIX3x3 a,REAL_MATRIX3x3 b)
{
  c->ax=a.ax-b.ax; c->bx=a.bx-b.bx; c->cx=a.cx-b.cx;
  c->ay=a.ay-b.ay; c->by=a.by-b.by; c->cy=a.cy-b.cy;
  c->az=a.az-b.az; c->bz=a.bz-b.bz; c->cz=a.cz-b.cz;
}

void DivideRealMatrix3x3ByReal(REAL_MATRIX3x3 *c,REAL_MATRIX3x3 a,REAL b)
{
  c->ax=a.ax/b; c->bx=a.bx/b; c->cx=a.cx/b;
  c->ay=a.ay/b; c->by=a.by/b; c->cy=a.cy/b;
  c->az=a.az/b; c->bz=a.bz/b; c->cz=a.cz/b;
}



void AddRealMatrix6x6(REAL_MATRIX6x6 *c,REAL_MATRIX6x6 a,REAL_MATRIX6x6 b)
{
  c->C11=a.C11+b.C11; c->C21=a.C21+b.C21; c->C31=a.C31+b.C31; c->C41=a.C41+b.C41; c->C51=a.C51+b.C51; c->C61=a.C61+b.C61;
  c->C12=a.C12+b.C12; c->C22=a.C22+b.C22; c->C32=a.C32+b.C32; c->C42=a.C42+b.C42; c->C52=a.C52+b.C52; c->C62=a.C62+b.C62;
  c->C13=a.C13+b.C13; c->C23=a.C23+b.C23; c->C33=a.C33+b.C33; c->C43=a.C43+b.C43; c->C53=a.C53+b.C53; c->C63=a.C63+b.C63;
  c->C14=a.C14+b.C14; c->C24=a.C24+b.C24; c->C34=a.C34+b.C34; c->C44=a.C44+b.C44; c->C54=a.C54+b.C54; c->C64=a.C64+b.C64;
  c->C15=a.C15+b.C15; c->C25=a.C25+b.C25; c->C35=a.C35+b.C35; c->C45=a.C45+b.C45; c->C55=a.C55+b.C55; c->C65=a.C65+b.C65;
  c->C16=a.C16+b.C16; c->C26=a.C26+b.C26; c->C36=a.C36+b.C36; c->C46=a.C46+b.C46; c->C56=a.C56+b.C56; c->C66=a.C66+b.C66;
}

void SubtractRealMatrix6x6(REAL_MATRIX6x6 *c,REAL_MATRIX6x6 a,REAL_MATRIX6x6 b)
{
  c->C11=a.C11-b.C11; c->C21=a.C21-b.C21; c->C31=a.C31-b.C31; c->C41=a.C41-b.C41; c->C51=a.C51-b.C51; c->C61=a.C61-b.C61;
  c->C12=a.C12-b.C12; c->C22=a.C22-b.C22; c->C32=a.C32-b.C32; c->C42=a.C42-b.C42; c->C52=a.C52-b.C52; c->C62=a.C62-b.C62;
  c->C13=a.C13-b.C13; c->C23=a.C23-b.C23; c->C33=a.C33-b.C33; c->C43=a.C43-b.C43; c->C53=a.C53-b.C53; c->C63=a.C63-b.C63;
  c->C14=a.C14-b.C14; c->C24=a.C24-b.C24; c->C34=a.C34-b.C34; c->C44=a.C44-b.C44; c->C54=a.C54-b.C54; c->C64=a.C64-b.C64;
  c->C15=a.C15-b.C15; c->C25=a.C25-b.C25; c->C35=a.C35-b.C35; c->C45=a.C45-b.C45; c->C55=a.C55-b.C55; c->C65=a.C65-b.C65;
  c->C16=a.C16-b.C16; c->C26=a.C26-b.C26; c->C36=a.C36-b.C36; c->C46=a.C46-b.C46; c->C56=a.C56-b.C56; c->C66=a.C66-b.C66;
}

void DivideRealMatrix6x6ByReal(REAL_MATRIX6x6 *c,REAL_MATRIX6x6 a,REAL b)
{
  c->C11=a.C11/b; c->C12=a.C12/b; c->C13=a.C13/b; c->C14=a.C14/b; c->C15=a.C15/b; c->C16=a.C16/b;
  c->C21=a.C21/b; c->C22=a.C22/b; c->C23=a.C23/b; c->C24=a.C24/b; c->C25=a.C25/b; c->C26=a.C26/b;
  c->C31=a.C31/b; c->C32=a.C32/b; c->C33=a.C33/b; c->C34=a.C34/b; c->C35=a.C35/b; c->C36=a.C36/b;
  c->C41=a.C41/b; c->C42=a.C42/b; c->C43=a.C43/b; c->C44=a.C44/b; c->C45=a.C45/b; c->C46=a.C46/b;
  c->C51=a.C51/b; c->C52=a.C52/b; c->C53=a.C53/b; c->C54=a.C54/b; c->C55=a.C55/b; c->C56=a.C56/b;
  c->C61=a.C61/b; c->C62=a.C62/b; c->C63=a.C63/b; c->C64=a.C64/b; c->C65=a.C65/b; c->C66=a.C66/b;
}

void AddRealMatrix9x9(REAL_MATRIX9x9 *c,REAL_MATRIX9x9 a,REAL_MATRIX9x9 b)
{
  c->xxxx=a.xxxx+b.xxxx;
  c->yxxx=a.yxxx+b.yxxx;
  c->zxxx=a.zxxx+b.zxxx;
  c->xxyx=a.xxyx+b.xxyx;
  c->yxyx=a.yxyx+b.yxyx;
  c->zxyx=a.zxyx+b.zxyx;
  c->xxzx=a.xxzx+b.xxzx;
  c->yxzx=a.yxzx+b.yxzx;
  c->zxzx=a.zxzx+b.zxzx;

  c->xyxx=a.xyxx+b.xyxx;
  c->yyxx=a.yyxx+b.yyxx;
  c->zyxx=a.zyxx+b.zyxx;
  c->xyyx=a.xyyx+b.xyyx;
  c->yyyx=a.yyyx+b.yyyx;
  c->zyyx=a.zyyx+b.zyyx;
  c->xyzx=a.xyzx+b.xyzx;
  c->yyzx=a.yyzx+b.yyzx;
  c->zyzx=a.zyzx+b.zyzx;

  c->xzxx=a.xzxx+b.xzxx;
  c->yzxx=a.yzxx+b.yzxx;
  c->zzxx=a.zzxx+b.zzxx;
  c->xzyx=a.xzyx+b.xzyx;
  c->yzyx=a.yzyx+b.yzyx;
  c->zzyx=a.zzyx+b.zzyx;
  c->xzzx=a.xzzx+b.xzzx;
  c->yzzx=a.yzzx+b.yzzx;
  c->zzzx=a.zzzx+b.zzzx;

  c->xxxy=a.xxxy+b.xxxy;
  c->yxxy=a.yxxy+b.yxxy;
  c->zxxy=a.zxxy+b.zxxy;
  c->xxyy=a.xxyy+b.xxyy;
  c->yxyy=a.yxyy+b.yxyy;
  c->zxyy=a.zxyy+b.zxyy;
  c->xxzy=a.xxzy+b.xxzy;
  c->yxzy=a.yxzy+b.yxzy;
  c->zxzy=a.zxzy+b.zxzy;

  c->xyxy=a.xyxy+b.xyxy;
  c->yyxy=a.yyxy+b.yyxy;
  c->zyxy=a.zyxy+b.zyxy;
  c->xyyy=a.xyyy+b.xyyy;
  c->yyyy=a.yyyy+b.yyyy;
  c->zyyy=a.zyyy+b.zyyy;
  c->xyzy=a.xyzy+b.xyzx;
  c->yyzy=a.yyzy+b.yyzy;
  c->zyzy=a.zyzy+b.zyzy;

  c->xzxy=a.xzxy+b.xzxy;
  c->yzxy=a.yzxy+b.yzxy;
  c->zzxy=a.zzxy+b.zzxy;
  c->xzyy=a.xzyy+b.xzyy;
  c->yzyy=a.yzyy+b.yzyy;
  c->zzyy=a.zzyy+b.zzyy;
  c->xzzy=a.xzzy+b.xzzy;
  c->yzzy=a.yzzy+b.yzzy;
  c->zzzy=a.zzzy+b.zzzy;

  c->xxxz=a.xxxz+b.xxxz;
  c->yxxz=a.yxxz+b.yxxz;
  c->zxxz=a.zxxz+b.zxxz;
  c->xxyz=a.xxyz+b.xxyz;
  c->yxyz=a.yxyz+b.yxyz;
  c->zxyz=a.zxyz+b.zxyz;
  c->xxzz=a.xxzz+b.xxzz;
  c->yxzz=a.yxzz+b.yxzz;
  c->zxzz=a.zxzz+b.zxzz;

  c->xyxz=a.xyxz+b.xyxz;
  c->yyxz=a.yyxz+b.yyxz;
  c->zyxz=a.zyxz+b.zyxz;
  c->xyyz=a.xyyz+b.xyyz;
  c->yyyz=a.yyyz+b.yyyz;
  c->zyyz=a.zyyz+b.zyyz;
  c->xyzz=a.xyzz+b.xyzz;
  c->yyzz=a.yyzz+b.yyzz;
  c->zyzz=a.zyzz+b.zyzz;

  c->xzxz=a.xzxz+b.xzxz;
  c->yzxz=a.yzxz+b.yzxz;
  c->zzxz=a.zzxz+b.zzxz;
  c->xzyz=a.xzyz+b.xzyz;
  c->yzyz=a.yzyz+b.yzyz;
  c->zzyz=a.zzyz+b.zzyz;
  c->xzzz=a.xzzz+b.xzzz;
  c->yzzz=a.yzzz+b.yzzz;
  c->zzzz=a.zzzz+b.zzzz;

}

void SubtractRealMatrix9x9(REAL_MATRIX9x9 *c,REAL_MATRIX9x9 a,REAL_MATRIX9x9 b)
{
  c->xxxx=a.xxxx-b.xxxx;
  c->yxxx=a.yxxx-b.yxxx;
  c->zxxx=a.zxxx-b.zxxx;
  c->xxyx=a.xxyx-b.xxyx;
  c->yxyx=a.yxyx-b.yxyx;
  c->zxyx=a.zxyx-b.zxyx;
  c->xxzx=a.xxzx-b.xxzx;
  c->yxzx=a.yxzx-b.yxzx;
  c->zxzx=a.zxzx-b.zxzx;

  c->xyxx=a.xyxx-b.xyxx;
  c->yyxx=a.yyxx-b.yyxx;
  c->zyxx=a.zyxx-b.zyxx;
  c->xyyx=a.xyyx-b.xyyx;
  c->yyyx=a.yyyx-b.yyyx;
  c->zyyx=a.zyyx-b.zyyx;
  c->xyzx=a.xyzx-b.xyzx;
  c->yyzx=a.yyzx-b.yyzx;
  c->zyzx=a.zyzx-b.zyzx;

  c->xzxx=a.xzxx-b.xzxx;
  c->yzxx=a.yzxx-b.yzxx;
  c->zzxx=a.zzxx-b.zzxx;
  c->xzyx=a.xzyx-b.xzyx;
  c->yzyx=a.yzyx-b.yzyx;
  c->zzyx=a.zzyx-b.zzyx;
  c->xzzx=a.xzzx-b.xzzx;
  c->yzzx=a.yzzx-b.yzzx;
  c->zzzx=a.zzzx-b.zzzx;

  c->xxxy=a.xxxy-b.xxxy;
  c->yxxy=a.yxxy-b.yxxy;
  c->zxxy=a.zxxy-b.zxxy;
  c->xxyy=a.xxyy-b.xxyy;
  c->yxyy=a.yxyy-b.yxyy;
  c->zxyy=a.zxyy-b.zxyy;
  c->xxzy=a.xxzy-b.xxzy;
  c->yxzy=a.yxzy-b.yxzy;
  c->zxzy=a.zxzy-b.zxzy;

  c->xyxy=a.xyxy-b.xyxy;
  c->yyxy=a.yyxy-b.yyxy;
  c->zyxy=a.zyxy-b.zyxy;
  c->xyyy=a.xyyy-b.xyyy;
  c->yyyy=a.yyyy-b.yyyy;
  c->zyyy=a.zyyy-b.zyyy;
  c->xyzy=a.xyzy-b.xyzx;
  c->yyzy=a.yyzy-b.yyzy;
  c->zyzy=a.zyzy-b.zyzy;

  c->xzxy=a.xzxy-b.xzxy;
  c->yzxy=a.yzxy-b.yzxy;
  c->zzxy=a.zzxy-b.zzxy;
  c->xzyy=a.xzyy-b.xzyy;
  c->yzyy=a.yzyy-b.yzyy;
  c->zzyy=a.zzyy-b.zzyy;
  c->xzzy=a.xzzy-b.xzzy;
  c->yzzy=a.yzzy-b.yzzy;
  c->zzzy=a.zzzy-b.zzzy;

  c->xxxz=a.xxxz-b.xxxz;
  c->yxxz=a.yxxz-b.yxxz;
  c->zxxz=a.zxxz-b.zxxz;
  c->xxyz=a.xxyz-b.xxyz;
  c->yxyz=a.yxyz-b.yxyz;
  c->zxyz=a.zxyz-b.zxyz;
  c->xxzz=a.xxzz-b.xxzz;
  c->yxzz=a.yxzz-b.yxzz;
  c->zxzz=a.zxzz-b.zxzz;

  c->xyxz=a.xyxz-b.xyxz;
  c->yyxz=a.yyxz-b.yyxz;
  c->zyxz=a.zyxz-b.zyxz;
  c->xyyz=a.xyyz-b.xyyz;
  c->yyyz=a.yyyz-b.yyyz;
  c->zyyz=a.zyyz-b.zyyz;
  c->xyzz=a.xyzz-b.xyzz;
  c->yyzz=a.yyzz-b.yyzz;
  c->zyzz=a.zyzz-b.zyzz;

  c->xzxz=a.xzxz-b.xzxz;
  c->yzxz=a.yzxz-b.yzxz;
  c->zzxz=a.zzxz-b.zzxz;
  c->xzyz=a.xzyz-b.xzyz;
  c->yzyz=a.yzyz-b.yzyz;
  c->zzyz=a.zzyz-b.zzyz;
  c->xzzz=a.xzzz-b.xzzz;
  c->yzzz=a.yzzz-b.yzzz;
  c->zzzz=a.zzzz-b.zzzz;

}

void DivideRealMatrix9x9ByReal(REAL_MATRIX9x9 *c,REAL_MATRIX9x9 a,REAL b)
{
  c->xxxx=a.xxxx/b;
  c->yxxx=a.yxxx/b;
  c->zxxx=a.zxxx/b;
  c->xxyx=a.xxyx/b;
  c->yxyx=a.yxyx/b;
  c->zxyx=a.zxyx/b;
  c->xxzx=a.xxzx/b;
  c->yxzx=a.yxzx/b;
  c->zxzx=a.zxzx/b;

  c->xyxx=a.xyxx/b;
  c->yyxx=a.yyxx/b;
  c->zyxx=a.zyxx/b;
  c->xyyx=a.xyyx/b;
  c->yyyx=a.yyyx/b;
  c->zyyx=a.zyyx/b;
  c->xyzx=a.xyzx/b;
  c->yyzx=a.yyzx/b;
  c->zyzx=a.zyzx/b;

  c->xzxx=a.xzxx/b;
  c->yzxx=a.yzxx/b;
  c->zzxx=a.zzxx/b;
  c->xzyx=a.xzyx/b;
  c->yzyx=a.yzyx/b;
  c->zzyx=a.zzyx/b;
  c->xzzx=a.xzzx/b;
  c->yzzx=a.yzzx/b;
  c->zzzx=a.zzzx/b;

  c->xxxy=a.xxxy/b;
  c->yxxy=a.yxxy/b;
  c->zxxy=a.zxxy/b;
  c->xxyy=a.xxyy/b;
  c->yxyy=a.yxyy/b;
  c->zxyy=a.zxyy/b;
  c->xxzy=a.xxzy/b;
  c->yxzy=a.yxzy/b;
  c->zxzy=a.zxzy/b;

  c->xyxy=a.xyxy/b;
  c->yyxy=a.yyxy/b;
  c->zyxy=a.zyxy/b;
  c->xyyy=a.xyyy/b;
  c->yyyy=a.yyyy/b;
  c->zyyy=a.zyyy/b;
  c->xyzy=a.xyzy/b;
  c->yyzy=a.yyzy/b;
  c->zyzy=a.zyzy/b;

  c->xzxy=a.xzxy/b;
  c->yzxy=a.yzxy/b;
  c->zzxy=a.zzxy/b;
  c->xzyy=a.xzyy/b;
  c->yzyy=a.yzyy/b;
  c->zzyy=a.zzyy/b;
  c->xzzy=a.xzzy/b;
  c->yzzy=a.yzzy/b;
  c->zzzy=a.zzzy/b;

  c->xxxz=a.xxxz/b;
  c->yxxz=a.yxxz/b;
  c->zxxz=a.zxxz/b;
  c->xxyz=a.xxyz/b;
  c->yxyz=a.yxyz/b;
  c->zxyz=a.zxyz/b;
  c->xxzz=a.xxzz/b;
  c->yxzz=a.yxzz/b;
  c->zxzz=a.zxzz/b;

  c->xyxz=a.xyxz/b;
  c->yyxz=a.yyxz/b;
  c->zyxz=a.zyxz/b;
  c->xyyz=a.xyyz/b;
  c->yyyz=a.yyyz/b;
  c->zyyz=a.zyyz/b;
  c->xyzz=a.xyzz/b;
  c->yyzz=a.yyzz/b;
  c->zyzz=a.zyzz/b;

  c->xzxz=a.xzxz/b;
  c->yzxz=a.yzxz/b;
  c->zzxz=a.zzxz/b;
  c->xzyz=a.xzyz/b;
  c->yzyz=a.yzyz/b;
  c->zzyz=a.zzyz/b;
  c->xzzz=a.xzzz/b;
  c->yzzz=a.yzzz/b;
  c->zzzz=a.zzzz/b;

}

void PrintRealMatrix9x9(REAL_MATRIX9x9 *m)
{
  fprintf(stderr, "% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xxxx,(double)m->yxxx,(double)m->zxxx,(double)m->xxyx,(double)m->yxyx,(double)m->zxyx,(double)m->xxzx,(double)m->yxzx,(double)m->zxzx);
  fprintf(stderr, "% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xyxx,(double)m->yyxx,(double)m->zyxx,(double)m->xyyx,(double)m->yyyx,(double)m->zyyx,(double)m->xyzx,(double)m->yyzx,(double)m->zyzx);
  fprintf(stderr, "% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xzxx,(double)m->yzxx,(double)m->zzxx,(double)m->xzyx,(double)m->yzyx,(double)m->zzyx,(double)m->xzzx,(double)m->yzzx,(double)m->zzzx);
  fprintf(stderr, "% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xxxy,(double)m->yxxy,(double)m->zxxy,(double)m->xxyy,(double)m->yxyy,(double)m->zxyy,(double)m->xxzy,(double)m->yxzy,(double)m->zxzy);
  fprintf(stderr, "% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xyxy,(double)m->yyxy,(double)m->zyxy,(double)m->xyyy,(double)m->yyyy,(double)m->zyyy,(double)m->xyzy,(double)m->yyzy,(double)m->zyzy);
  fprintf(stderr, "% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xzxy,(double)m->yzxy,(double)m->zzxy,(double)m->xzyy,(double)m->yzyy,(double)m->zzyy,(double)m->xzzy,(double)m->yzzy,(double)m->zzzy);
  fprintf(stderr, "% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xxxz,(double)m->yxxz,(double)m->zxxz,(double)m->xxyz,(double)m->yxyz,(double)m->zxyz,(double)m->xxzz,(double)m->yxzz,(double)m->zxzz);
  fprintf(stderr, "% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xyxz,(double)m->yyxz,(double)m->zyxz,(double)m->xyyz,(double)m->yyyz,(double)m->zyyz,(double)m->xyzz,(double)m->yyzz,(double)m->zyzz);
  fprintf(stderr, "% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xzxz,(double)m->yzxz,(double)m->zzxz,(double)m->xzyz,(double)m->yzyz,(double)m->zzyz,(double)m->xzzz,(double)m->yzzz,(double)m->zzzz);
  fprintf(stderr, "\n");
}

void PrintRealMatrix3x3ToFile(REAL_MATRIX3x3 *m,FILE *FilePtr,REAL a)
{
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f\n",
    (double)m->ax/a,(double)m->bx/a,(double)m->cx/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f\n",
    (double)m->ay/a,(double)m->by/a,(double)m->cy/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f\n",
    (double)m->az/a,(double)m->bz/a,(double)m->cz/a);
}

void PrintRealMatrix6x6ToFile(REAL_MATRIX6x6 *m,FILE *FilePtr,REAL a)
{
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C11/a,(double)m->C21/a,(double)m->C31/a,(double)m->C41/a,(double)m->C51/a,(double)m->C61/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C12/a,(double)m->C22/a,(double)m->C32/a,(double)m->C42/a,(double)m->C52/a,(double)m->C62/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C13/a,(double)m->C23/a,(double)m->C33/a,(double)m->C43/a,(double)m->C53/a,(double)m->C63/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C14/a,(double)m->C24/a,(double)m->C34/a,(double)m->C44/a,(double)m->C54/a,(double)m->C64/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C15/a,(double)m->C25/a,(double)m->C35/a,(double)m->C45/a,(double)m->C55/a,(double)m->C65/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C16/a,(double)m->C26/a,(double)m->C36/a,(double)m->C46/a,(double)m->C56/a,(double)m->C66/a);
}

void PrintRealMatrix6x6x6ToFile(REAL_MATRIX6x6x6 *m,FILE *FilePtr,REAL a)
{
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C111/a,(double)m->C211/a,(double)m->C311/a,(double)m->C411/a,(double)m->C511/a,(double)m->C611/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C121/a,(double)m->C221/a,(double)m->C321/a,(double)m->C421/a,(double)m->C521/a,(double)m->C621/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C131/a,(double)m->C231/a,(double)m->C331/a,(double)m->C431/a,(double)m->C531/a,(double)m->C631/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C141/a,(double)m->C241/a,(double)m->C341/a,(double)m->C441/a,(double)m->C541/a,(double)m->C641/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C151/a,(double)m->C251/a,(double)m->C351/a,(double)m->C451/a,(double)m->C551/a,(double)m->C651/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n\n",
    (double)m->C161/a,(double)m->C261/a,(double)m->C361/a,(double)m->C461/a,(double)m->C561/a,(double)m->C661/a);

  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C112/a,(double)m->C212/a,(double)m->C312/a,(double)m->C412/a,(double)m->C512/a,(double)m->C612/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C122/a,(double)m->C222/a,(double)m->C322/a,(double)m->C422/a,(double)m->C522/a,(double)m->C622/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C132/a,(double)m->C232/a,(double)m->C332/a,(double)m->C432/a,(double)m->C532/a,(double)m->C632/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C142/a,(double)m->C242/a,(double)m->C342/a,(double)m->C442/a,(double)m->C542/a,(double)m->C642/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C152/a,(double)m->C252/a,(double)m->C352/a,(double)m->C452/a,(double)m->C552/a,(double)m->C652/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n\n",
    (double)m->C162/a,(double)m->C262/a,(double)m->C362/a,(double)m->C462/a,(double)m->C562/a,(double)m->C662/a);

  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C113/a,(double)m->C213/a,(double)m->C313/a,(double)m->C413/a,(double)m->C513/a,(double)m->C613/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C123/a,(double)m->C223/a,(double)m->C323/a,(double)m->C423/a,(double)m->C523/a,(double)m->C623/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C133/a,(double)m->C233/a,(double)m->C333/a,(double)m->C433/a,(double)m->C533/a,(double)m->C633/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C143/a,(double)m->C243/a,(double)m->C343/a,(double)m->C443/a,(double)m->C543/a,(double)m->C643/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C153/a,(double)m->C253/a,(double)m->C353/a,(double)m->C453/a,(double)m->C553/a,(double)m->C653/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n\n",
    (double)m->C163/a,(double)m->C263/a,(double)m->C363/a,(double)m->C463/a,(double)m->C563/a,(double)m->C663/a);

  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C114/a,(double)m->C214/a,(double)m->C314/a,(double)m->C414/a,(double)m->C514/a,(double)m->C614/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C124/a,(double)m->C224/a,(double)m->C324/a,(double)m->C424/a,(double)m->C524/a,(double)m->C624/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C134/a,(double)m->C234/a,(double)m->C334/a,(double)m->C434/a,(double)m->C534/a,(double)m->C634/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C144/a,(double)m->C244/a,(double)m->C344/a,(double)m->C444/a,(double)m->C544/a,(double)m->C644/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C154/a,(double)m->C254/a,(double)m->C354/a,(double)m->C454/a,(double)m->C554/a,(double)m->C654/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n\n",
    (double)m->C164/a,(double)m->C264/a,(double)m->C364/a,(double)m->C464/a,(double)m->C564/a,(double)m->C664/a);

  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C115/a,(double)m->C215/a,(double)m->C315/a,(double)m->C415/a,(double)m->C515/a,(double)m->C615/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C125/a,(double)m->C225/a,(double)m->C325/a,(double)m->C425/a,(double)m->C525/a,(double)m->C625/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C135/a,(double)m->C235/a,(double)m->C335/a,(double)m->C435/a,(double)m->C535/a,(double)m->C635/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C145/a,(double)m->C245/a,(double)m->C345/a,(double)m->C445/a,(double)m->C545/a,(double)m->C645/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C155/a,(double)m->C255/a,(double)m->C355/a,(double)m->C455/a,(double)m->C555/a,(double)m->C655/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n\n",
    (double)m->C165/a,(double)m->C265/a,(double)m->C365/a,(double)m->C465/a,(double)m->C565/a,(double)m->C665/a);

  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C116/a,(double)m->C216/a,(double)m->C316/a,(double)m->C416/a,(double)m->C516/a,(double)m->C616/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C126/a,(double)m->C226/a,(double)m->C326/a,(double)m->C426/a,(double)m->C526/a,(double)m->C626/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C136/a,(double)m->C236/a,(double)m->C336/a,(double)m->C436/a,(double)m->C536/a,(double)m->C636/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C146/a,(double)m->C246/a,(double)m->C346/a,(double)m->C446/a,(double)m->C546/a,(double)m->C646/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->C156/a,(double)m->C256/a,(double)m->C356/a,(double)m->C456/a,(double)m->C556/a,(double)m->C656/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n\n",
    (double)m->C166/a,(double)m->C266/a,(double)m->C366/a,(double)m->C466/a,(double)m->C566/a,(double)m->C666/a);
}


void PrintRealMatrix9x9ToFile(REAL_MATRIX9x9 *m,FILE *FilePtr,REAL a)
{
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xxxx/a,(double)m->yxxx/a,(double)m->zxxx/a,(double)m->xxyx/a,(double)m->yxyx/a,(double)m->zxyx/a,(double)m->xxzx/a,(double)m->yxzx/a,(double)m->zxzx/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xyxx/a,(double)m->yyxx/a,(double)m->zyxx/a,(double)m->xyyx/a,(double)m->yyyx/a,(double)m->zyyx/a,(double)m->xyzx/a,(double)m->yyzx/a,(double)m->zyzx/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xzxx/a,(double)m->yzxx/a,(double)m->zzxx/a,(double)m->xzyx/a,(double)m->yzyx/a,(double)m->zzyx/a,(double)m->xzzx/a,(double)m->yzzx/a,(double)m->zzzx/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xxxy/a,(double)m->yxxy/a,(double)m->zxxy/a,(double)m->xxyy/a,(double)m->yxyy/a,(double)m->zxyy/a,(double)m->xxzy/a,(double)m->yxzy/a,(double)m->zxzy/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xyxy/a,(double)m->yyxy/a,(double)m->zyxy/a,(double)m->xyyy/a,(double)m->yyyy/a,(double)m->zyyy/a,(double)m->xyzy/a,(double)m->yyzy/a,(double)m->zyzy/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xzxy/a,(double)m->yzxy/a,(double)m->zzxy/a,(double)m->xzyy/a,(double)m->yzyy/a,(double)m->zzyy/a,(double)m->xzzy/a,(double)m->yzzy/a,(double)m->zzzy/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xxxz/a,(double)m->yxxz/a,(double)m->zxxz/a,(double)m->xxyz/a,(double)m->yxyz/a,(double)m->zxyz/a,(double)m->xxzz/a,(double)m->yxzz/a,(double)m->zxzz/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xyxz/a,(double)m->yyxz/a,(double)m->zyxz/a,(double)m->xyyz/a,(double)m->yyyz/a,(double)m->zyyz/a,(double)m->xyzz/a,(double)m->yyzz/a,(double)m->zyzz/a);
  fprintf(FilePtr,"% 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f % 12.5f\n",
    (double)m->xzxz/a,(double)m->yzxz/a,(double)m->zzxz/a,(double)m->xzyz/a,(double)m->yzyz/a,(double)m->zzyz/a,(double)m->xzzz/a,(double)m->yzzz/a,(double)m->zzzz/a);
  fprintf(FilePtr,"\n");
}

void PrintRealMatrixToFile(REAL_MATRIX *a,char *name)
{
  int i,j;
  FILE *FilePtr;

  FilePtr=fopen(name,"w");
  for(i=0;i<a->m;i++)
  {
    for(j=0;j<a->n;j++)
      fprintf(FilePtr,"%12.6f ",a->element[i][j]);
    fprintf(FilePtr,"\n");
  }
  fclose(FilePtr);
}

void InitializeMatrix9x9(REAL_MATRIX9x9 *m)
{
  m->xxxx=0.0; m->yxxx=0.0; m->zxxx=0.0;  m->xxyx=0.0; m->yxyx=0.0; m->zxyx=0.0;  m->xxzx=0.0; m->yxzx=0.0; m->zxzx=0.0;
  m->xyxx=0.0; m->yyxx=0.0; m->zyxx=0.0;  m->xyyx=0.0; m->yyyx=0.0; m->zyyx=0.0;  m->xyzx=0.0; m->yyzx=0.0; m->zyzx=0.0;
  m->xzxx=0.0; m->yzxx=0.0; m->zzxx=0.0;  m->xzyx=0.0; m->yzyx=0.0; m->zzyx=0.0;  m->xzzx=0.0; m->yzzx=0.0; m->zzzx=0.0;
  m->xxxy=0.0; m->yxxy=0.0; m->zxxy=0.0;  m->xxyy=0.0; m->yxyy=0.0; m->zxyy=0.0;  m->xxzy=0.0; m->yxzy=0.0; m->zxzy=0.0;
  m->xyxy=0.0; m->yyxy=0.0; m->zyxy=0.0;  m->xyyy=0.0; m->yyyy=0.0; m->zyyy=0.0;  m->xyzy=0.0; m->yyzy=0.0; m->zyzy=0.0;
  m->xzxy=0.0; m->yzxy=0.0; m->zzxy=0.0;  m->xzyy=0.0; m->yzyy=0.0; m->zzyy=0.0;  m->xzzy=0.0; m->yzzy=0.0; m->zzzy=0.0;
  m->xxxz=0.0; m->yxxz=0.0; m->zxxz=0.0;  m->xxyz=0.0; m->yxyz=0.0; m->zxyz=0.0;  m->xxzz=0.0; m->yxzz=0.0; m->zxzz=0.0;
  m->xyxz=0.0; m->yyxz=0.0; m->zyxz=0.0;  m->xyyz=0.0; m->yyyz=0.0; m->zyyz=0.0;  m->xyzz=0.0; m->yyzz=0.0; m->zyzz=0.0;
  m->xzxz=0.0; m->yzxz=0.0; m->zzxz=0.0;  m->xzyz=0.0; m->yzyz=0.0; m->zzyz=0.0;  m->xzzz=0.0; m->yzzz=0.0; m->zzzz=0.0;
}

void Convert9x9ToRealMatrix(REAL_MATRIX9x9 *m,REAL_MATRIX c)
{
  c.element[0][0]=m->xxxx; c.element[1][0]=m->yxxx; c.element[2][0]=m->zxxx; c.element[3][0]=m->xxyx; c.element[4][0]=m->yxyx;
                           c.element[5][0]=m->zxyx; c.element[6][0]=m->xxzx; c.element[7][0]=m->yxzx; c.element[8][0]=m->zxzx;
  c.element[0][1]=m->xyxx; c.element[1][1]=m->yyxx; c.element[2][1]=m->zyxx; c.element[3][1]=m->xyyx; c.element[4][1]=m->yyyx;
                           c.element[5][1]=m->zyyx; c.element[6][1]=m->xyzx; c.element[7][1]=m->yyzx; c.element[8][1]=m->zyzx;
  c.element[0][2]=m->xzxx; c.element[1][2]=m->yzxx; c.element[2][2]=m->zzxx; c.element[3][2]=m->xzyx; c.element[4][2]=m->yzyx;
                           c.element[5][2]=m->zzyx; c.element[6][2]=m->xzzx; c.element[7][2]=m->yzzx; c.element[8][2]=m->zzzx;
  c.element[0][3]=m->xxxy; c.element[1][3]=m->yxxy; c.element[2][3]=m->zxxy; c.element[3][3]=m->xxyy; c.element[4][3]=m->yxyy;
                           c.element[5][3]=m->zxyy; c.element[6][3]=m->xxzy; c.element[7][3]=m->yxzy; c.element[8][3]=m->zxzy;
  c.element[0][4]=m->xyxy; c.element[1][4]=m->yyxy; c.element[2][4]=m->zyxy; c.element[3][4]=m->xyyy; c.element[4][4]=m->yyyy;
                           c.element[5][4]=m->zyyy; c.element[6][4]=m->xyzy; c.element[7][4]=m->yyzy; c.element[8][4]=m->zyzy;
  c.element[0][5]=m->xzxy; c.element[1][5]=m->yzxy; c.element[2][5]=m->zzxy; c.element[3][5]=m->xzyy; c.element[4][5]=m->yzyy;
                           c.element[5][5]=m->zzyy; c.element[6][5]=m->xzzy; c.element[7][5]=m->yzzy; c.element[8][5]=m->zzzy;
  c.element[0][6]=m->xxxz; c.element[1][6]=m->yxxz; c.element[2][6]=m->zxxz; c.element[3][6]=m->xxyz; c.element[4][6]=m->yxyz;
                           c.element[5][6]=m->zxyz; c.element[6][6]=m->xxzz; c.element[7][6]=m->yxzz; c.element[8][6]=m->zxzz;
  c.element[0][7]=m->xyxz; c.element[1][7]=m->yyxz; c.element[2][7]=m->zyxz; c.element[3][7]=m->xyyz; c.element[4][7]=m->yyyz;
                           c.element[5][7]=m->zyyz; c.element[6][7]=m->xyzz; c.element[7][7]=m->yyzz; c.element[8][7]=m->zyzz;
  c.element[0][8]=m->xzxz; c.element[1][8]=m->yzxz; c.element[2][8]=m->zzxz; c.element[3][8]=m->xzyz; c.element[4][8]=m->yzyz;
                           c.element[5][8]=m->zzyz; c.element[6][8]=m->xzzz; c.element[7][8]=m->yzzz; c.element[8][8]=m->zzzz;

}


void PrintRealMatrix3x3(REAL_MATRIX3x3 *m)
{
  fprintf(stderr, "%16.10f %16.10f %16.10f\n",(double)(m->ax),(double)(m->bx),(double)(m->cx));
  fprintf(stderr, "%16.10f %16.10f %16.10f\n",(double)(m->ay),(double)(m->by),(double)(m->cy));
  fprintf(stderr, "%16.10f %16.10f %16.10f\n",(double)(m->az),(double)(m->bz),(double)(m->cz));
}


INT_MATRIX CreateIntMatrix(int m,int n)
{
  int i;
  INT_MATRIX c;

  c.m=m;
  c.n=n;
  c.element=(int**)calloc(m,sizeof(int*));
  c.element[0]=(int*)calloc(m*n,sizeof(int));

  for(i=1;i<m;i++)
    c.element[i]=c.element[i-1]+m;
  return c;
}

void DeleteIntMatrix(INT_MATRIX m)
{
  free(m.element[0]);
  free(m.element);
}

// creates a matrix that is accessible in both
// the fortran- and c-convention
REAL_MATRIX CreateRealMatrix(int m,int n)
{
  int i;
  REAL_MATRIX c;

  c.m=m;
  c.n=n;
  c.element=(REAL**)calloc(m,sizeof(REAL*));
  c.element[0]=(REAL*)calloc(m*n,sizeof(REAL));

  for(i=1;i<m;i++)
    c.element[i]=c.element[i-1]+n;
  return c;
}

void DeleteRealMatrix(REAL_MATRIX m)
{
  free(m.element[0]);
  free(m.element);
}

REAL_FORTRAN_MATRIX CreateRealFortranMatrix(int m,int n)
{
  REAL_FORTRAN_MATRIX c;

  c.m=m;
  c.n=n;
  c.element=(REAL*)calloc(m*n,sizeof(REAL));
  return c;
}

void DeleteRealFortranMatrix(REAL_FORTRAN_MATRIX m)
{
  free(m.element);
}


COMPLEX_MATRIX CreateComplexMatrix(int m,int n)
{
  int i;
  COMPLEX_MATRIX c;

  c.m=m;
  c.n=n;
  c.element=(Complex**)calloc(m,sizeof(Complex*));
  c.element[0]=(Complex*)calloc(m*n,sizeof(Complex));

  for(i=1;i<m;i++)
    c.element[i]=c.element[i-1]+n;
  return c;
}

void DeleteComplexMatrix(COMPLEX_MATRIX m)
{
  free(m.element[0]);
  free(m.element);
}


POINT_MATRIX CreatePointMatrix(int m,int n)
{
  int i;
  POINT_MATRIX c;

  c.m=m;
  c.n=n;
  c.element=(POINT**)calloc(m,sizeof(POINT*));
  c.element[0]=(POINT*)calloc(m*n,sizeof(POINT));

  for(i=1;i<m;i++)
    c.element[i]=c.element[i-1]+n;
  return c;
}

void DeletePointMatrix(POINT_MATRIX m)
{
  free(m.element[0]);
  free(m.element);
}

void PrintRealMatrix(REAL_MATRIX *c)
{
  int i,j;
  int n,m;

  m=c->m;
  n=c->n;

  for(i=0;i<m;i++)
  {
    for(j=0;j<n;j++)
      fprintf(stderr, "%18.10lf ",(double)c->element[i][j]);
    fprintf(stderr, "\n");
  }
}

void PrintComplexMatrix(COMPLEX_MATRIX *c)
{
  int i,j;
  int n,m;

  m=c->m;
  n=c->n;

  for(i=0;i<m;i++)
  {
    for(j=0;j<n;j++)
      fprintf(stderr, "(%10.5lf,%10.5lf) ",(double)c->element[i][j].re,(double)c->element[i][j].im);
    fprintf(stderr, "\n");
  }
}


void PrintRealMatrixmMathematica(REAL_MATRIX *c)
{
  int i,j;
  int n,m;

  m=c->m;
  n=c->n;

  fprintf(stderr, "{");
  for(i=0;i<m;i++)
  {
    fprintf(stderr, "{");
    for(j=0;j<n;j++)
      fprintf(stderr, "%lf,",(double)c->element[i][j]);
    fprintf(stderr, "},");
  }
  fprintf(stderr, "}\n");
}



void PrintRealFortranMatrix(REAL_FORTRAN_MATRIX *c)
{
  int i,j,size;

  size=c->m*c->n;
  for(i=0;i<c->m;i++)
  {
    for(j=0;j<c->n;j++)
     fprintf(stderr, "%lg ",(double)c->element[i+j*(c->m)]);
    fprintf(stderr, "\n");
  }
}

REAL Trace3x3Matrix(REAL_MATRIX3x3 *in)
{
  return (in->ax+in->by+in->cz);
}

REAL_MATRIX3x3 MatrixMatrixMultiplication3x3(REAL_MATRIX3x3 a,REAL_MATRIX3x3 b)
{
  REAL_MATRIX3x3 c;

  c.ax=a.ax*b.ax+a.bx*b.ay+a.cx*b.az;
  c.ay=a.ay*b.ax+a.by*b.ay+a.cy*b.az;
  c.az=a.az*b.ax+a.bz*b.ay+a.cz*b.az;

  c.bx=a.ax*b.bx+a.bx*b.by+a.cx*b.bz;
  c.by=a.ay*b.bx+a.by*b.by+a.cy*b.bz;
  c.bz=a.az*b.bx+a.bz*b.by+a.cz*b.bz;

  c.cx=a.ax*b.cx+a.bx*b.cy+a.cx*b.cz;
  c.cy=a.ay*b.cx+a.by*b.cy+a.cy*b.cz;
  c.cz=a.az*b.cx+a.bz*b.cy+a.cz*b.cz;

/*
  c.ax=b.ax*a.ax+b.bx*a.ay+b.cx*a.az;
  c.ay=b.ay*a.ax+b.by*a.ay+b.cy*a.az;
  c.az=b.az*a.ax+b.bz*a.ay+b.cz*a.az;
  c.bx=b.ax*a.bx+b.bx*a.by+b.cx*a.bz;
  c.by=b.ay*a.bx+b.by*a.by+b.cy*a.bz;
  c.bz=b.az*a.bx+b.bz*a.by+b.cz*a.bz;
  c.cx=b.ax*a.cx+b.bx*a.cy+b.cx*a.cz;
  c.cy=b.ay*a.cx+b.by*a.cy+b.cy*a.cz;
  c.cz=b.az*a.cx+b.bz*a.cy+b.cz*a.cz;
*/


  return c;
}



VECTOR MatrixVectorMultiplication3x3(REAL_MATRIX3x3 a,VECTOR b)
{
  VECTOR c;

  c.x=a.ax*b.x+a.bx*b.y+a.cx*b.z;
  c.y=a.ay*b.x+a.by*b.y+a.cy*b.z;
  c.z=a.az*b.x+a.bz*b.y+a.cz*b.z;
  return c;
}


//  subroutine to invert a 3 * 3 matrix using cofactors
void Invert3x3Matrix(REAL_MATRIX3x3 *in,REAL_MATRIX3x3 *out,REAL *determinant)
{
  REAL r,d;

  // calculate adjoint matrix
  out->ax=in->by*in->cz-in->bz*in->cy;
  out->ay=in->az*in->cy-in->ay*in->cz;
  out->az=in->ay*in->bz-in->az*in->by;
  out->bx=in->bz*in->cx-in->bx*in->cz;
  out->by=in->ax*in->cz-in->az*in->cx;
  out->bz=in->az*in->bx-in->ax*in->bz;
  out->cx=in->bx*in->cy-in->by*in->cx;
  out->cy=in->ay*in->cx-in->ax*in->cy;
  out->cz=in->ax*in->by-in->ay*in->bx;

  // calculate determinant
  d=in->ax*out->ax+in->bx*out->ay+in->cx*out->az;
  r=0.0;
  if(fabs(d)>0.0) r=1.0/d;

  // complete inverse matrix
  out->ax*=r; out->bx*=r; out->cx*=r;
  out->ay*=r; out->by*=r; out->cy*=r;
  out->az*=r; out->bz*=r; out->cz*=r;
  *determinant=d;
}

REAL GetLargestElementMatrix3x3(REAL_MATRIX3x3 *a)
{
  VECTOR largest;

  largest.x=MAX3(a->ax,a->bx,a->cx);
  largest.y=MAX3(a->ay,a->by,a->cy);
  largest.z=MAX3(a->az,a->bz,a->cz);
  return MAX3(largest.x,largest.y,largest.z);
}

void RescaleMatrix3x3x(REAL_MATRIX3x3 *a,REAL c)
{
  a->ax*=c; a->bx*=c; a->cx*=c;
  a->ay*=c; a->by*=c; a->cy*=c;
  a->az*=c; a->bz*=c; a->cz*=c;
}

void TransposeMatrix3x3(REAL_MATRIX3x3 *in)
{
  REAL temp;

  temp=in->bx;
  in->bx=in->ay;
  in->ay=temp;

  temp=in->cx;
  in->cx=in->az;
  in->az=temp;

  temp=in->cy;
  in->cy=in->bz;
  in->bz=temp;
}

void EigenSystem3x3(REAL_MATRIX3x3 in,REAL_MATRIX3x3 *eigenvectors,VECTOR *eigenvalues)
{
  int i,j,k,n;
  const REAL rho=1.0e-16;
  REAL tes,scl;
  int pass;
  REAL a[3][3],v[3][3];
  REAL v1,v2,v3,tem,c,omg,s,u;

  // initialize eigenvectors
  eigenvectors->ax=eigenvectors->bx=eigenvectors->cx=0.0;
  eigenvectors->ay=eigenvectors->by=eigenvectors->cy=0.0;
  eigenvectors->az=eigenvectors->bz=eigenvectors->cz=0.0;

  eigenvalues->x=0.0;
  eigenvalues->y=0.0;
  eigenvalues->z=0.0;

  v[0][0]=1.0; v[1][0]=0.0; v[2][0]=0.0;
  v[0][1]=0.0; v[1][1]=1.0; v[2][1]=0.0;
  v[0][2]=0.0; v[1][2]=0.0; v[2][2]=1.0;

  n=3;
  a[0][0]=in.ax; a[1][0]=in.bx; a[2][0]=in.cx;
  a[0][1]=in.ay; a[1][1]=in.by; a[2][1]=in.cy;
  a[0][2]=in.az; a[1][2]=in.bz; a[2][2]=in.cz;

  // rescale matrix for optimal accuracy
  scl=1.0;
  for(i=1;i<=n;i++)
    if(fabs(a[i-1][i-1])>scl) scl=fabs(a[i-1][i-1]);
  for(i=1;i<=n;i++)
    for(j=1;j<=i;j++)
      a[i-1][j-1]/=scl;

  // set initial value of moving tolerance
  tes=0.0;
  for(i=2;i<=n;i++)
    for(j=1;j<=i-1;j++)
      tes+=2.0*a[i-1][j-1]*a[i-1][j-1];
  tes=sqrt(tes);

  do
  {
    tes/=(REAL)n;
    if(tes<rho) tes=rho;

    // jacobi diagonalisation
    do
    {
      pass=FALSE;
      for(i=2;i<=n;i++)
      {
        for(j=1;j<=i-1;j++)
        {
          if(fabs(a[i-1][j-1])>=tes)
          {
            pass=TRUE;
            v1=a[j-1][j-1];
            v2=a[i-1][j-1];
            v3=a[i-1][i-1];
            u=0.5*(v1-v3);
            if(fabs(u)<rho)
              omg=-1.0;
            else
            {
              omg=-v2/sqrt(v2*v2+u*u);
              if(u<0.0) omg=-omg;
            }
            s=omg/sqrt(2.0*(1.0+sqrt(1.0-omg*omg)));
            c=sqrt(1.0-s*s);
            for(k=1;k<=n;k++)
            {
              if(k>=i)
              {
                tem=a[k-1][j-1]*c-a[k-1][i-1]*s;
                a[k-1][i-1]=a[k-1][j-1]*s+a[k-1][i-1]*c;
                a[k-1][j-1]=tem;
              }
              else if(k<j)
              {
                tem=a[j-1][k-1]*c-a[i-1][k-1]*s;
                a[i-1][k-1]=a[j-1][k-1]*s+a[i-1][k-1]*c;
                a[j-1][k-1]=tem;
              }
              else
              {
                tem=a[k-1][j-1]*c-a[i-1][k-1]*s;
                a[i-1][k-1]=a[k-1][j-1]*s+a[i-1][k-1]*c;
                a[k-1][j-1]=tem;
              }
              tem=v[k-1][j-1]*c-v[k-1][i-1]*s;
              v[k-1][i-1]=v[k-1][j-1]*s+v[k-1][i-1]*c;
              v[k-1][j-1]=tem;
            }
            a[j-1][j-1]=v1*c*c+v3*s*s-2.0*v2*s*c;
            a[i-1][i-1]=v1*s*s+v3*c*c+2.0*v2*s*c;
            a[i-1][j-1]=(v1-v3)*s*c+v2*(c*c-s*s);
          }
        }
      }
    //recycle until moving tolerance satisfied
    }
    while(pass);
    // recycle until absolute tolerance satisfied
  }
  while(tes>rho);

  // rescale matrix
  for(i=1;i<=n;i++)
    for(j=1;j<=i;j++)
      a[i-1][j-1]*=scl;

  eigenvectors->ax=v[0][0]; eigenvectors->bx=v[1][0]; eigenvectors->cx=v[2][0];
  eigenvectors->ay=v[0][1]; eigenvectors->by=v[1][1]; eigenvectors->cy=v[2][1];
  eigenvectors->az=v[0][2]; eigenvectors->bz=v[1][2]; eigenvectors->cz=v[2][2];

  eigenvalues->x=a[0][0];
  eigenvalues->y=a[1][1];
  eigenvalues->z=a[2][2];
}

void EigenSystemVoigt3x3(REAL_MATRIX3x3 in,REAL *eigenvalues)
{
  int i;
  REAL d[3],e[3];
  REAL a[9];

  a[0]=in.ax; a[1]=in.bx; a[2]=in.cx;
  a[3]=in.ay; a[4]=in.by; a[5]=in.cy;
  a[6]=in.az; a[7]=in.bz; a[8]=in.cz;

  tred2(a,3,d,e);
  tqli(d,e,3,a);
  //eigsrt(d,a,6);

  for(i=0;i<3;i++)
    eigenvalues[i]=d[i];
}


void EigenSystemVoigt6x6(REAL_MATRIX6x6 in,REAL *eigenvalues)
{
  int i;
  REAL d[6],e[6];
  REAL a[36];

  a[0]=in.C11;  a[1]=in.C12;  a[2]=in.C13;  a[3]=in.C14;  a[4]=in.C15;  a[5]=in.C16;
  a[6]=in.C12;  a[7]=in.C22;  a[8]=in.C23;  a[9]=in.C24;  a[10]=in.C25; a[11]=in.C26;
  a[12]=in.C13; a[13]=in.C23; a[14]=in.C33; a[15]=in.C34; a[16]=in.C35; a[17]=in.C36;
  a[18]=in.C14; a[19]=in.C24; a[20]=in.C34; a[21]=in.C44; a[22]=in.C45; a[23]=in.C46;
  a[24]=in.C15; a[25]=in.C25; a[26]=in.C35; a[27]=in.C45; a[28]=in.C55; a[29]=in.C56;
  a[30]=in.C16; a[31]=in.C26; a[32]=in.C36; a[33]=in.C46; a[34]=in.C56; a[35]=in.C66;

  tred2(a,6,d,e);
  tqli(d,e,6,a);
  //eigsrt(d,a,6);

  for(i=0;i<6;i++)
    eigenvalues[i]=d[i];
}


REAL pythag(REAL a, REAL b)
{
  REAL absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return ((absb == 0.0)? (REAL)0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void eigsrt(REAL d[], REAL *v, int n)
{
  int k,j,i;
  REAL p;

  for(i=1;i<n;i++)
  {
    k=i;
    p=d[k-1];
    for (j=i+1;j<=n;j++)
      if (d[j-1] < p)
      {
        k=j;
        p=d[k-1];
      }
    if(k!=i)
    {
      d[k-1]=d[i-1];
      d[i-1]=p;
      for (j=1;j<=n;j++)
      {
        p=v[j-1+(i-1)*n];
        v[j-1+(i-1)*n]=v[j-1+(k-1)*n];
        v[j-1+(k-1)*n]=p;
      }
    }
  }
}

void eigsrt3(REAL d[], REAL *v, int n)
{
  int k,j,i;
  REAL p;

  for(i=4;i<n;i++)
  {
    k=i;
    p=d[k-1];
    for (j=i+1;j<=n;j++)
      if (d[j-1] < p)
      {
        k=j;
        p=d[k-1];
      }
    if(k!=i)
    {
      d[k-1]=d[i-1];
      d[i-1]=p;
      for (j=4;j<=n;j++)
      {
        p=v[j-1+(i-1)*n];
        v[j-1+(i-1)*n]=v[j-1+(k-1)*n];
        v[j-1+(k-1)*n]=p;
      }
    }
  }
}


void tred2(REAL *a, int n, REAL d[], REAL e[])
{
  int l,k,j,i;
  REAL scale,hh,h,g,f;

  for (i=n;i>=2;i--)
  {
    l=i-1;
    h=scale=0.0;
    if(l>1)
    {
      for(k=1;k<=l;k++)
        scale += fabs(a[i-1+(k-1)*n]);
      if (scale == 0.0)
        e[i-1]=a[i-1+(l-1)*n];
      else
      {
        for (k=1;k<=l;k++)
        {
          a[i-1+(k-1)*n] /= scale;
          h += a[i-1+(k-1)*n]*a[i-1+(k-1)*n];
        }
        f=a[i-1+(l-1)*n];
        g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i-1]=scale*g;
        h -= f*g;
        a[i-1+(l-1)*n]=f-g;
        f=0.0;
        for (j=1;j<=l;j++)
        {
          a[j-1+(i-1)*n]=a[i-1+(j-1)*n]/h;
          g=0.0;
          for (k=1;k<=j;k++)
            g += a[j-1+(k-1)*n]*a[i-1+(k-1)*n];
          for (k=j+1;k<=l;k++)
            g += a[k-1+(j-1)*n]*a[i-1+(k-1)*n];
          e[j-1]=g/h;
          f += e[j-1]*a[i-1+(j-1)*n];
        }
        hh=f/(h+h);
        for (j=1;j<=l;j++)
        {
          f=a[i-1+(j-1)*n];
          e[j-1]=g=e[j-1]-hh*f;
          for (k=1;k<=j;k++)
          a[j-1+(k-1)*n] -= (f*e[k-1]+g*a[i-1+(k-1)*n]);
        }
      }
    }
    else
      e[i-1]=a[i-1+(l-1)*n];
    d[i-1]=h;
  }
  d[1-1]=0.0;
  e[1-1]=0.0;
  for (i=1;i<=n;i++)
  {
    l=i-1;
    if(d[i-1]!=0.0)
    {
      for (j=1;j<=l;j++)
      {
        g=0.0;
        for (k=1;k<=l;k++)
          g += a[i-1+(k-1)*n]*a[k-1+(j-1)*n];
        for (k=1;k<=l;k++)
          a[k-1+(j-1)*n] -= g*a[k-1+(i-1)*n];
      }
    }
    d[i-1]=a[i-1+(i-1)*n];
    a[i-1+(i-1)*n]=1.0;
    for (j=1;j<=l;j++)
      a[j-1+(i-1)*n]=a[i-1+(j-1)*n]=0.0;
  }
}

void tqli(REAL d[], REAL e[], int n, REAL *z)
{
  int m,l,iter,i,k;
  REAL s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1-1]=e[i-1];
    e[n-1]=0.0;
  for (l=1;l<=n;l++)
  {
    iter=0;
    do
    {
      for (m=l;m<=n-1;m++)
      {
        dd=fabs(d[m-1])+fabs(d[m+1-1]);
        if ((REAL)(fabs(e[m-1])+dd) == dd) break;
      }
      if(m!=l)
      {
        if (iter++ == 100) fprintf(stderr, "Too many iterations in tqli\n");
        g=(d[l+1-1]-d[l-1])/(2.0*e[l-1]);
        r=pythag(g,1.0);
        g=d[m-1]-d[l-1]+e[l-1]/(g+(g>=0.0?fabs(r):-fabs(r)));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--)
        {
          f=s*e[i-1];
          b=c*e[i-1];
          e[i+1-1]=(r=pythag(f,g));
          if (r == 0.0)
          {
            d[i+1-1] -= p;
            e[m-1]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1-1]-p;
          r=(d[i-1]-g)*s+2.0*c*b;
          d[i+1-1]=g+(p=s*r);
          g=c*r-b;
          for (k=1;k<=n;k++)
          {
            f=z[k-1+(i+1-1)*n];
            z[k-1+(i+1-1)*n]=s*z[k-1+(i-1)*n]+c*f;
            z[k-1+(i-1)*n]=c*z[k-1+(i-1)*n]-s*f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l-1] -= p;
        e[l-1]=g;
        e[m-1]=0.0;
      }
    } while (m != l);
  }
}

#define ROTATE(a,i,j,k,l,n) g=a[(i-1)+n*(j-1)];h=a[(k-1)+n*(l-1)];a[(i-1)+n*(j-1)]=g-s*(h+g*tau);\
  a[(k-1)+n*(l-1)]=h+s*(g-h*tau);

void jacobi(REAL *a,int n,REAL d[],REAL *v,int *nrot)
{
  int j,iq,ip,i;
  REAL tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  b=(REAL *)calloc(n,sizeof(REAL));
  z=(REAL *)calloc(n,sizeof(REAL));
  for (ip=1;ip<=n;ip++)
  {
    for (iq=1;iq<=n;iq++) v[(ip-1)+n*(iq-1)]=0.0;
      v[(ip-1)+n*(ip-1)]=1.0;
  }
  for (ip=1;ip<=n;ip++)
  {
    b[ip-1]=d[ip-1]=a[(ip-1)+n*(ip-1)];
    z[ip-1]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++)
  {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++)
    {
      for (iq=ip+1;iq<=n;iq++)
        sm += fabs(a[(ip-1)+n*(iq-1)]);
    }
    if (sm == 0.0)
    {
      free(z);
      free(b);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(REAL)(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++)
    {
      for (iq=ip+1;iq<=n;iq++)
      {
        g=100.0*fabs(a[(ip-1)+n*(iq-1)]);
        if (i > 4 && (REAL)(fabs(d[ip-1])+g) == (REAL)fabs(d[ip-1])
            && (REAL)(fabs(d[iq-1])+g) == (REAL)fabs(d[iq-1]))
          a[(ip-1)+n*(iq-1)]=0.0;
        else if (fabs(a[(ip-1)+n*(iq-1)]) > tresh)
        {
          h=d[iq-1]-d[ip-1];
          if ((REAL)(fabs(h)+g) == (REAL)fabs(h))
            t=(a[(ip-1)+n*(iq-1)])/h;
          else
          {
            theta=0.5*h/(a[(ip-1)+n*(iq-1)]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1.0+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[(ip-1)+n*(iq-1)];
          z[ip-1] -= h;
          z[iq-1] += h;
          d[ip-1] -= h;
          d[iq-1] += h;
          a[(ip-1)+n*(iq-1)]=0.0;
          for (j=1;j<=ip-1;j++)
          {
            ROTATE(a,j,ip,j,iq,n)
          }
          for (j=ip+1;j<=iq-1;j++)
          {
            ROTATE(a,ip,j,j,iq,n)
          }
          for (j=iq+1;j<=n;j++)
          {
            ROTATE(a,ip,j,iq,j,n)
          }
          for (j=1;j<=n;j++)
          {
            ROTATE(v,j,ip,j,iq,n)
          }
          ++(*nrot);
        }
      }
    }
    for (ip=1;ip<=n;ip++)
    {
      b[ip-1] += z[ip-1];
      d[ip-1]=b[ip-1];
      z[ip-1]=0.0;
    }
  }
  fprintf(stderr, "Too many iterations in routine jacobi\n");
}
#undef ROTATE

// Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
// is the input matrix. b[1..n][1..m] is input containing the m right-hand side vectors. On
// output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
// vectors.

void GaussJordan(REAL_MATRIX a,int n,REAL_MATRIX b,int m)
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  REAL big,dum,pivinv,temp;

  indxc=(int*)calloc(n,sizeof(int)); // The integer arrays ipiv, indxr, and indxc are used
  indxr=(int*)calloc(n,sizeof(int)); // for bookkeeping on the pivoting.
  ipiv=(int*)calloc(n,sizeof(int));

  icol=irow=0;

  for (j=0;j<n;j++) ipiv[j]=0;
  for (i=0;i<n;i++)
  {
    //This is the main loop over the columns to be reduced.
    big=0.0;
    for (j=0;j<n;j++)  //This is the outer loop of the search for a pivot element.
    {
      if (ipiv[j] != 1)
      {
        for (k=0;k<n;k++)
        {
          if (ipiv[k] == 0)
          {
            if (fabs(a.element[j][k]) >= big)
            {
              big=fabs(a.element[j][k]);
              irow=j;
              icol=k;
            }
          }
        }
      }
    }
    ++(ipiv[icol]);
    // We now have the pivot element, so we interchange rows, if needed, to put the pivot
    // element on the diagonal. The columns are not physically interchanged, only relabeled:
    // indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
    // indxr[i] is the row in which that pivot element was originally located. If indxr[i] =
    // indxc[i] there is an implied column interchange. With this form of bookkeeping, the
    // solution b’s will end up in the correct order, and the inverse matrix will be scrambled
    // by columns.
    if(irow != icol)
    {
      for (l=0;l<n;l++) SWAP(a.element[irow][l],a.element[icol][l],temp)
      for (l=0;l<m;l++) SWAP(b.element[irow][l],b.element[icol][l],temp)
    }
    indxr[i]=irow; // We are now ready to divide the pivot row by the
    indxc[i]=icol; // pivot element, located at irow and icol.
    if(a.element[icol][icol]==0.0)
    {
      fprintf(stderr, "Matrix Inversion, Gauss-Jordan: Singular Matrix\n");
    }
    pivinv=1.0/a.element[icol][icol];
    a.element[icol][icol]=1.0;
    for (l=0;l<n;l++) a.element[icol][l] *= pivinv;
    for (l=0;l<m;l++) b.element[icol][l] *= pivinv;

    for (ll=0;ll<n;ll++) // Next, we reduce the rows...
    {
      if (ll != icol) // ...except for the pivot one, of course.
      {
        dum=a.element[ll][icol];
        a.element[ll][icol]=0.0;
        for (l=0;l<n;l++) a.element[ll][l] -= a.element[icol][l]*dum;
        for (l=0;l<m;l++) b.element[ll][l] -= b.element[icol][l]*dum;
      }
    }
  }
  // This is the end of the main loop over columns of the reduction. It only remains to unscramble
  // the solution in view of the column interchanges. We do this by interchanging pairs of
  // columns in the reverse order that the permutation was built up.
  for (l=n-1;l>=0;l--)
  {
    if (indxr[l] != indxc[l])
      for (k=0;k<n;k++)
        SWAP(a.element[k][indxr[l]],a.element[k][indxc[l]],temp);
  }
  free(ipiv);
  free(indxr);
  free(indxc);
}

void CheckMatrixInversion(void)
{
  int i,j,k,l,n;
  REAL_MATRIX a,as,b,check;

  n=5;
  a=CreateRealMatrix(5,5);
  as=CreateRealMatrix(5,5);
  b=CreateRealMatrix(5,5);
  check=CreateRealMatrix(5,5);

  a.element[0][0]=4;  a.element[1][0]=3;  a.element[2][0]=1; a.element[3][0]=2;  a.element[4][0]=2;
  a.element[0][1]=0;  a.element[1][1]=-3; a.element[2][1]=2; a.element[3][1]=-2; a.element[4][1]=-2;
  a.element[0][2]=14; a.element[1][2]=3;  a.element[2][2]=3; a.element[3][2]=-2; a.element[4][2]=-2;
  a.element[0][3]=-4; a.element[1][3]=-3; a.element[2][3]=4; a.element[3][3]=0;  a.element[4][3]=0;
  a.element[0][4]=-3; a.element[1][4]=-3; a.element[2][4]=4; a.element[3][4]=0;  a.element[4][4]=1;
  PrintRealMatrix(&a);
  as.element[0][0]=4;  as.element[1][0]=3;  as.element[2][0]=1; as.element[3][0]=2;  as.element[4][0]=2;
  as.element[0][1]=0;  as.element[1][1]=-3; as.element[2][1]=2; as.element[3][1]=-2; as.element[4][1]=-2;
  as.element[0][2]=14; as.element[1][2]=3;  as.element[2][2]=3; as.element[3][2]=-2; as.element[4][2]=-2;
  as.element[0][3]=-4; as.element[1][3]=-3; as.element[2][3]=4; as.element[3][3]=0;  as.element[4][3]=0;
  as.element[0][4]=-3; as.element[1][4]=-3; as.element[2][4]=4; as.element[3][4]=0;  as.element[4][4]=1;

  for(i=0;i<5;i++)
   b.element[i][i]=1.0;

  GaussJordan(a,n,b,n);

  fprintf(stderr, "\n");
  PrintRealMatrix(&a);

  for(k=0;k<n;k++)
  {
    for(l=0;l<n;l++)
    {
      check.element[k][l]=0.0;
      for (j=0;j<n;j++)
        check.element[k][l]+=(as.element[k][j]*a.element[j][l]);
    }
    for(l=0;l<n;l++) fprintf(stderr, "%14.6f",check.element[k][l]);
    fprintf(stderr, "\n");
  }

  a.element[0][0]=4;  a.element[1][0]=3;  a.element[2][0]=1; a.element[3][0]=2;  a.element[4][0]=2;
  a.element[0][1]=0;  a.element[1][1]=-3; a.element[2][1]=2; a.element[3][1]=-2; a.element[4][1]=-2;
  a.element[0][2]=14; a.element[1][2]=3;  a.element[2][2]=3; a.element[3][2]=-2; a.element[4][2]=-2;
  a.element[0][3]=-4; a.element[1][3]=-3; a.element[2][3]=4; a.element[3][3]=0;  a.element[4][3]=0;
  a.element[0][4]=-3; a.element[1][4]=-3; a.element[2][4]=4; a.element[3][4]=0;  a.element[4][4]=1;
  SingularValueDecompositionMatrixInversion(a);
  PrintRealMatrix(&a);

  DeleteRealMatrix(a);
}

void InverseRealMatrix(REAL_MATRIX a)
{
  int n;
  REAL_MATRIX b;

  n=a.n;

  b=CreateRealMatrix(n,n);
  GaussJordan(a,n,b,n);
  DeleteRealMatrix(b);
}

#ifdef HAVE_LAPACK

extern void zgetri_(int *,COMPLEX *,int *,int *,COMPLEX *,int *,int *);

extern void zgetrf_(int *m,int *n,COMPLEX *a,int *lda,int *ipiv,int *info);

void InverseComplexMatrix(COMPLEX_MATRIX a)
{
  int n,m;
  int *ipiv;
  int info,lwork;
  COMPLEX *work;

  m=a.m;
  n=a.n;

  // LU decomposition
  ipiv=(int*)calloc(n,sizeof(int));
  zgetrf_(&n,&n,&a.element[0][0],&n,ipiv,&info);

  // ask for the amount of required memory
  lwork=-1;
  work=calloc(1,sizeof(COMPLEX));
  zgetri_(&n,&a.element[0][0],&n,ipiv,work,&lwork,&info);
  lwork=(int)work[0].re;
  free(work);

  // do the real call
  work=calloc(lwork,sizeof(COMPLEX));
  zgetri_(&n,&a.element[0][0],&n,ipiv,work,&lwork,&info);
}

#else

// TODO: internal complex inverse routine
void InverseComplexMatrix(COMPLEX_MATRIX a)
{
}

#endif


/*
void InverseComplexMatrix(Matrix RealA, Matrix ImagA,Matrix& RealAinv, Matrix& ImagAinv )
// Inputs
//   RealA  -    Real part of matrix A (N by N)
//   ImagA  -    Imaginary part of matrix A (N by N)
// Outputs
//   RealAinv  -    Real part of inverse of matrix A (N by N)
//   ImagAinv  -    Imaginary part of A inverse (N by N)
{

  int N = RealA.nRow();
  assert( N == RealA.nCol() && N == ImagA.nRow()
                            && N == ImagA.nCol());
    RealAinv = RealA; // Copy matrices to ensure they are same size
  ImagAinv = ImagA;

  int i, j, k;
  Matrix scale(N);   // Scale factor
  int *index;  index = new int [N+1];

  // Matrix B is initialized to the identity matrix
  Matrix RealB(N,N), ImagB(N,N);
  RealB.set(0.0);  ImagB.set(0.0);
  for( i=1; i<=N; i++ )
    RealB(i,i) = 1.0;

  // Set scale factor, scale(i) = max( |a(i,j)| ), for each row
  for( i=1; i<=N; i++ ) {
    index[i] = i;       // Initialize row index list
    double scaleMax = 0.;
    for( j=1; j<=N; j++ ) {
      double MagA = RealA(i,j)*RealA(i,j) + ImagA(i,j)*ImagA(i,j);
      scaleMax = (scaleMax > MagA) ? scaleMax : MagA;
    }
    scale(i) = scaleMax;
  }

  // Loop over rows k = 1, ..., (N-1)
  for( k=1; k<=N-1; k++ ) {
    // Select pivot row from max( |a(j,k)/s(j)| )
    double ratiomax = 0.0;
    int jPivot = k;
    for( i=k; i<=N; i++ ) {
      double MagA = RealA(index[i],k)*RealA(index[i],k) +
                    ImagA(index[i],k)*ImagA(index[i],k);
      double ratio = MagA/scale(index[i]);
      if( ratio > ratiomax ) {
        jPivot=i;
        ratiomax = ratio;
      }
    }
    // Perform pivoting using row index list
    int indexJ = index[k];
    if( jPivot != k ) {           // Pivot
      indexJ = index[jPivot];
      index[jPivot] = index[k];   // Swap index jPivot and k
      index[k] = indexJ;
    }
    // Perform forward elimination
    for( i=k+1; i<=N; i++ ) {
      double denom = RealA(indexJ,k)*RealA(indexJ,k)
                   + ImagA(indexJ,k)*ImagA(indexJ,k);
      double RealCoeff = (RealA(index[i],k)*RealA(indexJ,k)
                       + ImagA(index[i],k)*ImagA(indexJ,k))/denom;
      double ImagCoeff = (ImagA(index[i],k)*RealA(indexJ,k)
                       - RealA(index[i],k)*ImagA(indexJ,k))/denom;
      for( j=k+1; j<=N; j++ ) {
        RealA(index[i],j) -= RealCoeff*RealA(indexJ,j)
                           - ImagCoeff*ImagA(indexJ,j);
        ImagA(index[i],j) -= RealCoeff*ImagA(indexJ,j)
                           + ImagCoeff*RealA(indexJ,j);
      }
      RealA(index[i],k) = RealCoeff;
      ImagA(index[i],k) = ImagCoeff;
      for( j=1; j<=N; j++ ) {
        RealB(index[i],j) -= RealA(index[i],k)*RealB(indexJ,j)
                           - ImagA(index[i],k)*ImagB(indexJ,j);
        ImagB(index[i],j) -= RealA(index[i],k)*ImagB(indexJ,j)
                           + ImagA(index[i],k)*RealB(indexJ,j);
      }
    }
  }

  // Perform backsubstitution
  for( k=1; k<=N; k++ ) {
    double denom = RealA(index[N],N)*RealA(index[N],N)
                 + ImagA(index[N],N)*ImagA(index[N],N);
    RealAinv(N,k) = (RealB(index[N],k)*RealA(index[N],N)
                  + ImagB(index[N],k)*ImagA(index[N],N))/denom;
    ImagAinv(N,k) = (ImagB(index[N],k)*RealA(index[N],N)
                  - RealB(index[N],k)*ImagA(index[N],N))/denom;
    for( i=N-1; i>=1; i--) {
      double RealSum = RealB(index[i],k);
      double ImagSum = ImagB(index[i],k);
      for( j=i+1; j<=N; j++ ) {
        RealSum -= RealA(index[i],j)*RealAinv(j,k)
                 - ImagA(index[i],j)*ImagAinv(j,k);
        ImagSum -= RealA(index[i],j)*ImagAinv(j,k)
                 + ImagA(index[i],j)*RealAinv(j,k);
      }
      double denom = RealA(index[i],i)*RealA(index[i],i)
                   + ImagA(index[i],i)*ImagA(index[i],i);
      RealAinv(i,k) = (RealSum*RealA(index[i],i)
                    + ImagSum*ImagA(index[i],i))/denom;
      ImagAinv(i,k) = (ImagSum*RealA(index[i],i)
                    - RealSum*ImagA(index[i],i))/denom;
    }
  }

  delete [] index;  // Release allocated memory
}
*/

/*
  int ndata;
  int num_terms;
  double *x,*y,*std_dev;

  ndata=5;
  num_terms=2;
  x=(double*)calloc(ndata,sizeof(double));
  y=(double*)calloc(ndata,sizeof(double));
  std_dev=(double*)calloc(ndata,sizeof(double));
  x[0]=0.0; y[0]=1.0; std_dev[0]=0.05;
  x[1]=1.0; y[1]=2.5; std_dev[1]=0.05;
  x[2]=2.0; y[2]=3.1; std_dev[2]=0.05;
  x[3]=3.0; y[3]=3.8; std_dev[3]=0.05;
  x[4]=4.0; y[4]=5.5; std_dev[4]=0.05;

  PolynomialFit(x,y,std_dev,ndata,num_terms);
*/


void funcs( double x, double *afunc, unsigned int ma );
void SingularValueDecompositionFit( double *X, double *Y, double *Sig, unsigned int NData,
    double *A, unsigned int MA,
    double **U, double **V, double *W, unsigned int MP, unsigned int NP,
    double *ChiSq, void funcs(double x, double *afunc, unsigned int ma) );
void svdcmp( double **A, unsigned int M, unsigned int N,double *W, double **V );
void SingularValueDecompositionBackSubstitution( double **U, double *W, double **V, unsigned int M,
    unsigned int N,double *B, double *X );
void SingularValueDecompositionCovarianceMatrix( double **V, unsigned int MA,double *W, double **CVM);

# define true ((int)1)
# define false ((int)0)
# define nmax ((int)1000)
# define mmax ((int)50)
# define tol ((double)1.0e-5)

void fpoly(double x, double *p,unsigned int np)
{
  int j;

  p[0]=1.0;
  for(j=1;j<np;j++)
    p[j]=p[j-1]*x;
}

void SingularValueDecompositionFit( double *X, double *Y, double *Sig, unsigned int NData,
    double *A, unsigned int MA,
    double **U, double **V, double *W, unsigned int MP, unsigned int NP,
    double *ChiSq, void funcs(double x, double *afunc, unsigned int ma) )
{
    /*
       Given a set of NData points X[], Y[] with individual standard
       deviations of Sig[], use chi-square minimization to determine the
       MA coefficients, A[], of the fitting function
       y = sum over i A_i * func_i(x).
       Here we solve the fitting equation using singular value decomposition
       of the NData by MA matrix. The arrays U, V and W provide workspace
       on input. On output they define the singular value decomposition and
       can be used to obtaint he covariance matrix. MP and NP are the
       physical dimensions of the matrices U, V, and W as indicated below.
       It is necessary that MP be greater than or equal to NData and that
       NP be greather than or equal to MP. The program returns values for
       the MA fit parameters A[] and the chi-square, ChiSq. The user
       supplies a subroutine, funcs(), that returns the MA basis functions
       evaluated at x in the array afunc[].
    */

    int i;
    int j;
    double sum;
    double thresh;
    double tmp;
    double wmax;
    double wmin;

    double beta[nmax];
    double afunc[mmax];

    /* Accumulate coefficients of the fitting matrix. */
    for( i = 0; i < NData; ++i ) {
        funcs( X[i], afunc, MA );
        tmp = 1.0 / Sig[i];
        for( j = 0; j < MA; ++j ) {
            U[i][j] = afunc[j] * tmp;
        }
        beta[i] = Y[i] * tmp;
    }

    /* Singular value decomposition. */
    svdcmp( U, NData, MA, W, V );

    /* Edit the singular values, given tol from the parameter statement,
       between here ... */
    wmax = 0.0;
    wmin = 1.0e99;
    for( j = 0; j < MA; ++j ) {
        if( W[j] > wmax )
            wmax = W[j];
        if( W[j] < wmin )
            wmin = W[j];
    }

    thresh = tol * wmax;
    for( j = 0; j < MA; ++j ) {
        if( W[j] < thresh ) {
            W[j] = 0.0;
        }
    }
    /* ... and here. */

    SingularValueDecompositionBackSubstitution( U, W, V, NData, MA, beta, A );

    /* Evaluate chi-square. */
    *ChiSq = 0.0;
    for( i = 0; i < NData; ++i ) {
        funcs( X[i], afunc, MA );
        sum = 0.0;
        for( j = 0; j < MA; ++j ) {
            sum = sum + A[j] * afunc[j];
        }
        tmp = ((Y[i] - sum) / Sig[i]);
        *ChiSq = *ChiSq + tmp*tmp;
    }

    return;
}

void SingularValueDecompositionCovarianceMatrix( double **V, unsigned int MA,double *W, double **CVM)
{
    /*
       To evaluate the covariance matrix CVM of the fit for MA paramaters
       obtained by svdfit, call this routine with matrix V and W as returned
       from svdfit. NP, NCVM give the physical dimensions of V, W and CVM as
       indicated below.
    */

    int i;
    int j;
    int k;
    double sum;

    double wti[mmax];

    for( i = 0; i < MA; ++i ) {
        wti[i] = 0.0;
        if( W[i] != 0.0 )
            wti[i] = 1.0 / (W[i] * W[i]);
    }

    for( i = 0; i < MA; ++i ) {
        for( j = 0; j <= i; ++j ) {
            sum = 0.0;
            for( k = 0; k < MA; ++k ) {
                sum+=V[i][k] * V[j][k] * wti[k];
            }
            CVM[i][j] = sum;
            CVM[j][i] = sum;
        }
    }

    return;
}

void SingularValueDecompositionBackSubstitution(double **U,double *W,double **V,unsigned int M,
                                                unsigned int N,double *B,double *X)
{
  int i;
  int j;
  double S;

  double tmp[nmax];

  /* Calculate transpose U * B */
  for(j=0;j<N;j++)
  {
    S=0.0;
    /* Nonzero result only if W[j] is nonzero. */
    if(W[j]!=0.0)
    {
      for(i=0;i<M;i++)
        S+=U[i][j]*B[i];
      S/=W[j];
    }
    tmp[j]=S;
  }

  /* Multiply by V to get answer. */
  for(j=0;j<N;j++)
  {
    S=0.0;
    for(i=0;i<N;i++)
      S+=V[j][i]*tmp[i];
    X[j]=S;
  }

  return;
}

void svdcmp( double **A, unsigned int M, unsigned int N,double *W, double **V )
{
    /*
       Give a matrix A, with logical dimensions M by N and physical
       dimensions MP by NP, this routine computes its singular value
       decomposition, A = U * W * transpose V. The matrix U replaces
       A on output. The diagonal matrix of singular values, W, is output
       as a vector W. The matrix V (not the transpose of V) is output as
       V. M must be greater or equal to N. If it is smaller then A should
       be filled up to square with zero rows.
    */

    double rv1[nmax];

    /* Householder reduction to bidiagonal form. */
    int NM;
    double C;
    double F;
    double G = 0.0;
    double H;
    double S;
    double X;
    double Y;
    double Z;
    double Scale = 0.0;
    double ANorm = 0.0;
    double tmp;
    int flag;
    int i;
    int its;
    int j;
    int jj;
    int k;
    int l;

    if( M < N ) {
        fprintf( stderr, "You must augment A with extra zero rows.\n" );
        return;
    }

    for( i = 0; i < N; ++i ) {
        l = i + 1;
        rv1[i] = Scale * G;
        G = 0.0;
        S = 0.0;
        Scale = 0.0;
        if( i < M ) {
            for( k = i; k < M; ++k ) {
                Scale = Scale + fabs( A[k][i] );
            }
            if( Scale != 0.0 ) {
                for( k = i; k < M; ++k ) {
                    A[k][i] = A[k][i] / Scale;
                    S = S + A[k][i] * A[k][i];
                }
                F = A[i][i];
                G = sqrt(S);
                if( F > 0.0 ) {
                    G = -G;
                }
                H = F * G - S;
                A[i][i] = F - G;
                if( i != (N-1) ) {
                    for( j = l; j < N; ++j ) {
                        S = 0.0;
                        for( k = i; k < M; ++k ) {
                            S = S + A[k][i] * A[k][j];
                        }
                        F = S / H;
                        for( k = i; k < M; ++k ) {
                            A[k][j] = A[k][j] + F * A[k][i];
                        }
                    }
                }
                for( k = i; k < M; ++k ) {
                    A[k][i] = Scale * A[k][i];
                }
            }
        }

        W[i] = Scale * G;
        G = 0.0;
        S = 0.0;
        Scale = 0.0;
        if( (i < M) && (i != (N-1)) ) {
            for( k = l; k < N; ++k ) {
                Scale = Scale + fabs( A[i][k] );
            }
            if( Scale != 0.0 ) {
                for( k = l; k < N; ++k ) {
                    A[i][k] = A[i][k] / Scale;
                    S = S + A[i][k] * A[i][k];
                }
                F = A[i][l];
                G = sqrt(S);
                if( F > 0.0 ) {
                    G = -G;
                }
                H = F * G - S;
                A[i][l] = F - G;
                for( k = l; k < N; ++k ) {
                    rv1[k] = A[i][k] / H;
                }
                if( i != (M-1) ) {
                    for( j = l; j < M; ++j ) {
                        S = 0.0;
                        for( k = l; k < N; ++k ) {
                            S = S + A[j][k] * A[i][k];
                        }
                        for( k = l; k < N; ++k ) {
                            A[j][k] = A[j][k] + S * rv1[k];
                        }
                    }
                }
                for( k = l; k < N; ++k ) {
                    A[i][k] = Scale * A[i][k];
                }
            }
        }
        tmp = fabs( W[i] ) + fabs( rv1[i] );
        if( tmp > ANorm )
            ANorm = tmp;
    }

    /* Accumulation of right-hand transformations. */
    for( i = N-1; i >= 0; --i ) {
        if( i < (N-1) ) {
            if( G != 0.0 ) {
                for( j = l; j < N; ++j ) {
                    V[j][i] = (A[i][j] / A[i][l]) / G;
                }
                for( j = l; j < N; ++j ) {
                    S = 0.0;
                    for( k = l; k < N; ++k ) {
                        S = S + A[i][k] * V[k][j];
                    }
                    for( k = l; k < N; ++k ) {
                        V[k][j] = V[k][j] + S * V[k][i];
                    }
                }
            }
            for( j = l; j < N; ++j ) {
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        }
        V[i][i] = 1.0;
        G = rv1[i];
        l = i;
    }

    /* Accumulation of left-hand transformations. */
    for( i = N-1; i >= 0; --i ) {
        l = i + 1;
        G = W[i];
        if( i < (N-1) ) {
            for( j = l; j < N; ++j ) {
                A[i][j] = 0.0;
            }
        }
        if( G != 0.0 ) {
            G = 1.0 / G;
            if( i != (N-1) ) {
                for( j = l; j < N; ++j ) {
                    S = 0.0;
                    for( k = l; k < M; ++k ) {
                        S = S + A[k][i] * A[k][j];
                    }
                    F = (S / A[i][i]) * G;
                    for( k = i; k < M; ++k ) {
                        A[k][j] = A[k][j] + F * A[k][i];
                    }
                }
            }
            for( j = i; j < M; ++j ) {
                A[j][i] = A[j][i] * G;
            }
        } else {
            for( j = i; j < M; ++j ) {
                A[j][i] = 0.0;
            }
        }
        A[i][i] = A[i][i] + 1.0;
    }

    /* Diagonalization of the bidiagonal form.
       Loop over singular values. */
    for( k = (N-1); k >= 0; --k ) {
        /* Loop over allowed iterations. */
        for( its = 1; its <= 30; ++its ) {
            /* Test for splitting.
               Note that rv1[0] is always zero. */
            flag = true;
            for( l = k; l >= 0; --l ) {
                NM = l - 1;
                if( (fabs(rv1[l]) + ANorm) == ANorm ) {
                    flag = false;
                    break;
                } else if( (fabs(W[NM]) + ANorm) == ANorm ) {
                    break;
                }
            }

            /* Cancellation of rv1[l], if l > 0; */
            if( flag ) {
                C = 0.0;
                S = 1.0;
                for( i = l; i <= k; ++i ) {
                    F = S * rv1[i];
                    if( (fabs(F) + ANorm) != ANorm ) {
                        G = W[i];
                        H = sqrt( F * F + G * G );
                        W[i] = H;
                        H = 1.0 / H;
                        C = ( G * H );
                        S = -( F * H );
                        for( j = 0; j < M; ++j ) {
                            Y = A[j][NM];
                            Z = A[j][i];
                            A[j][NM] = (Y * C) + (Z * S);
                            A[j][i] = -(Y * S) + (Z * C);
                        }
                    }
                }
            }
            Z = W[k];
            /* Convergence. */
            if( l == k ) {
                /* Singular value is made nonnegative. */
                if( Z < 0.0 ) {
                    W[k] = -Z;
                    for( j = 0; j < N; ++j ) {
                        V[j][k] = -V[j][k];
                    }
                }
                break;
            }

            if( its >= 30 ) {
                fprintf( stderr, "No convergence in 30 iterations.\n" );
                return;
            }

            X = W[l];
            NM = k - 1;
            Y = W[NM];
            G = rv1[NM];
            H = rv1[k];
            F = ((Y-Z)*(Y+Z) + (G-H)*(G+H)) / (2.0*H*Y);
            G = sqrt( F * F + 1.0 );
            tmp = G;
            if( F < 0.0 )
                tmp = -tmp;
            F = ((X-Z)*(X+Z) + H*((Y/(F+tmp))-H)) / X;

            /* Next QR transformation. */
            C = 1.0;
            S = 1.0;
            for( j = l; j <= NM; ++j ) {
                i = j + 1;
                G = rv1[i];
                Y = W[i];
                H = S * G;
                G = C * G;
                Z = sqrt( F * F + H * H );
                rv1[j] = Z;
                C = F / Z;
                S = H / Z;
                F = (X * C) + (G * S);
                G = -(X * S) + (G * C);
                H = Y * S;
                Y = Y * C;
                for( jj = 0; jj < N; ++jj ) {
                    X = V[jj][j];
                    Z = V[jj][i];
                    V[jj][j] = (X * C) + (Z * S);
                    V[jj][i] = -(X * S) + (Z * C);
                }
                Z = sqrt( F * F + H * H );
                W[j] = Z;

                /* Rotation can be arbitrary if Z = 0. */
                if( Z != 0.0 ) {
                    Z = 1.0 / Z;
                    C = F * Z;
                    S = H * Z;
                }
                F = (C * G) + (S * Y);
                X = -(S * G) + (C * Y);
                for( jj = 0; jj < M; ++jj ) {
                    Y = A[jj][j];
                    Z = A[jj][i];
                    A[jj][j] = (Y * C) + (Z * S);
                    A[jj][i] = -(Y * S) + (Z * C);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = F;
            W[k] = X;
        }
    }

    return;
}


void PolynomialFit(double *x,double *y,double *std_dev,unsigned int ndata,int num_terms)
{
  int i,j;
/*
 *  Singular value decomposition storage.
 *    OrbFunc        = routine which computes basis function values
 *    afunc          = evaluated basis function values
 *    x_coef         = X interpolation coefficients
 *    y_coef         = Y interpolation coefficients
 *    z_coef         = Z interpolation coefficients
 *    u              = singular value decomposition workspace
 *    v              = singular value decomposition workspace
 *    w              = singular value decomposition workspace
 *    cvm            = singular value decomposition covariance matrix
 *    mp             = number of data points
 *    np             = number of coefficients
 */
    double chisq;
    double *x_coef = NULL;
    double *afunc = NULL;
    double **cvm = NULL;
    double *vm = NULL;
    double **u = NULL;
    double **v = NULL;
    double *w = NULL;

    unsigned int mp;
    unsigned int np;

/*
 *  Command line argument storage
 *    n_data       = the number of data points to use in when computing
 *                     interpolation terms
 *    num_terms      = the number of interpolation terms
 */
        if(ndata >= num_terms)
          mp = ndata;
        else
          mp = num_terms;
        np = num_terms;

        x_coef = (double *)malloc( num_terms * sizeof( double ) );
        for( i = 0; i < num_terms; ++i )
            x_coef[i] = 1.0;

        afunc = (double *)malloc( num_terms * sizeof( double ) );
        for( i = 0; i < num_terms; ++i )
          afunc[i] = 0.0;

        u = (double **)malloc( mp * sizeof( double * ) );
        for( i = 0; i < mp; ++i )
        {
            u[i] = NULL;
            u[i] = (double *)malloc( np * sizeof( double ) );
            for( j = 0; j < np; ++j )
                u[i][j] = 0.0;
        }

        v = (double **)malloc( np * sizeof( double * ) );
        for( i = 0; i < np; ++i )
        {
            v[i] = NULL;
            v[i] = (double *)malloc( np * sizeof( double ) );
            for( j = 0; j < np; ++j )
              v[i][j] = 0.0;
        }

        w = (double *)malloc( np * sizeof( double ) );
        for( i = 0; i < np; ++i )
            w[i] = 0.0;

        cvm = (double **)malloc( np * sizeof( double * ) );
        vm = (double *)malloc( np * sizeof( double) );
        for( i = 0; i < np; ++i )
        {
            cvm[i] = (double *)malloc( np * sizeof( double ) );
            for( j = 0; j < np; ++j )
                cvm[i][j] = 0.0;
        }

        chisq = 0.0;
        SingularValueDecompositionFit(x,y,std_dev,ndata,x_coef,num_terms,u, v,w,mp,np,&chisq,fpoly);
        SingularValueDecompositionCovarianceMatrix(v,num_terms,w,cvm);

        for( i = 0; i < num_terms; ++i )
          fprintf(stderr, "FIT-term %d: %g +/- %g\n",i,x_coef[i],sqrt(cvm[i][i])*sqrt(chisq/(ndata-num_terms)));
        for( i = 0; i < num_terms; ++i )
        {
          for( j = 0; j < num_terms; ++j )
            fprintf(stderr, "%g ",cvm[i][j]);
          fprintf(stderr, "\n");
        }

        fprintf(stderr, "Chi-sq: %g\n",chisq);

        exit(0);
}

void SingularValueDecomposition(REAL_MATRIX a,int m,int n,REAL *w,REAL_MATRIX v)
{
  int flag,i,its,j,jj,k,l,nm;
  REAL anorm,c,f,g,h,s,scale,x,y,z,*rv1;

  rv1=(REAL*)calloc(n,sizeof(REAL));
  g=scale=anorm=0.0;
  for (i=1;i<=n;i++)
  {
    l=i+1;
    rv1[i-1]=scale*g;
    g=s=scale=0.0;
    if (i <= m)
    {
      for (k=i;k<=m;k++) scale += fabs(a.element[k-1][i-1]);
      if (scale)
      {
        for (k=i;k<=m;k++)
        {
          a.element[k-1][i-1] /= scale;
          s += a.element[k-1][i-1]*a.element[k-1][i-1];
        }
        f=a.element[i-1][i-1];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a.element[i-1][i-1]=f-g;
        for (j=l;j<=n;j++)
        {
          for (s=0.0,k=i;k<=m;k++) s += a.element[k-1][i-1]*a.element[k-1][j-1];
          f=s/h;
          for (k=i;k<=m;k++) a.element[k-1][j-1] += f*a.element[k-1][i-1];
        }
        for (k=i;k<=m;k++) a.element[k-1][i-1] *= scale;
      }
    }
    w[i-1]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n)
    {
      for (k=l;k<=n;k++) scale += fabs(a.element[i-1][k-1]);
      if (scale)
      {
        for (k=l;k<=n;k++)
        {
          a.element[i-1][k-1] /= scale;
          s += a.element[i-1][k-1]*a.element[i-1][k-1];
        }
        f=a.element[i-1][l-1];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a.element[i-1][l-1]=f-g;
        for (k=l;k<=n;k++) rv1[k-1]=a.element[i-1][k-1]/h;
        for (j=l;j<=m;j++)
        {
          for (s=0.0,k=l;k<=n;k++) s += a.element[j-1][k-1]*a.element[i-1][k-1];
          for (k=l;k<=n;k++) a.element[j-1][k-1] += s*rv1[k-1];
        }
        for (k=l;k<=n;k++) a.element[i-1][k-1] *= scale;
      }
    }
    anorm=MAX2(anorm,(fabs(w[i-1])+fabs(rv1[i-1])));
  }

  for (i=n;i>=1;i--)
  {
    if (i < n)
    {
      if (g)
      {
        for (j=l;j<=n;j++)
          v.element[j-1][i-1]=(a.element[i-1][j-1]/a.element[i-1][l-1])/g;
        for (j=l;j<=n;j++)
        {
          for (s=0.0,k=l;k<=n;k++) s += a.element[i-1][k-1]*v.element[k-1][j-1];
          for (k=l;k<=n;k++) v.element[k-1][j-1] += s*v.element[k-1][i-1];
        }
      }
      for (j=l;j<=n;j++) v.element[i-1][j-1]=v.element[j-1][i-1]=0.0;
    }
    v.element[i-1][i-1]=1.0;
    g=rv1[i-1];
    l=i;
  }

  for (i=MIN2(m,n);i>=1;i--)
  {
    l=i+1;
    g=w[i-1];
    for (j=l;j<=n;j++) a.element[i-1][j-1]=0.0;
    if (g)
    {
      g=1.0/g;
      for (j=l;j<=n;j++)
      {
        for (s=0.0,k=l;k<=m;k++) s += a.element[k-1][i-1]*a.element[k-1][j-1];
        f=(s/a.element[i-1][i-1])*g;
        for (k=i;k<=m;k++) a.element[k-1][j-1] += f*a.element[k-1][i-1];
      }
      for (j=i;j<=m;j++) a.element[j-1][i-1] *= g;
    }
    else for (j=i;j<=m;j++) a.element[j-1][i-1]=0.0;
    ++a.element[i-1][i-1];
  }

  for (k=n;k>=1;k--)
  {
    for (its=1;its<=30;its++)
    {
      flag=1;
      for (l=k;l>=1;l--)
      {
        nm=l-1;
        if ((REAL)(fabs(rv1[l-1])+anorm) == anorm)
        {
          flag=0;
          break;
        }
        if ((REAL)(fabs(w[nm-1])+anorm) == anorm) break;
      }
      if (flag)
      {
        c=0.0;
        s=1.0;
        for (i=l;i<=k;i++)
        {
          f=s*rv1[i-1];
          rv1[i-1]=c*rv1[i-1];
          if ((REAL)(fabs(f)+anorm) == anorm) break;
          g=w[i-1];
          h=pythag(f,g);
          w[i-1]=h;
          h=1.0/h;
          c=g*h;
          s = -f*h;
          for (j=1;j<=m;j++)
          {
            y=a.element[j-1][nm-1];
            z=a.element[j-1][i-1];
            a.element[j-1][nm-1]=y*c+z*s;
            a.element[j-1][i-1]=z*c-y*s;
          }
        }
      }
      z=w[k-1];
      if (l == k)
      {
        if (z < 0.0)
        {
          w[k-1] = -z;
          for (j=1;j<=n;j++) v.element[j-1][k-1] = -v.element[j-1][k-1];
        }
        break;
      }
      if (its == 30)
      {
        fprintf(stderr, "no convergence in 30 svdcmp iterations");
        exit(0);
      }
      x=w[l-1];
      nm=k-1;
      y=w[nm-1];
      g=rv1[nm-1];
      h=rv1[k-1];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++)
      {
        i=j+1;
        g=rv1[i-1];
        y=w[i-1];
        h=s*g;
        g=c*g;
        z=pythag(f,h);
        rv1[j-1]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g = g*c-x*s;
        h=y*s;
        y *= c;
        for (jj=1;jj<=n;jj++)
        {
          x=v.element[jj-1][j-1];
          z=v.element[jj-1][i-1];
          v.element[jj-1][j-1]=x*c+z*s;
          v.element[jj-1][i-1]=z*c-x*s;
        }
        z=pythag(f,h);
        w[j-1]=z;
        if (z)
        {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for (jj=1;jj<=m;jj++)
        {
          y=a.element[jj-1][j-1];
          z=a.element[jj-1][i-1];
          a.element[jj-1][j-1]=y*c+z*s;
          a.element[jj-1][i-1]=z*c-y*s;
        }
      }
      rv1[l-1]=0.0;
      rv1[k-1]=f;
      w[k-1]=x;
    }
  }
  free(rv1);
}

void SingularValueDecompositionMatrixInversion(REAL_MATRIX a)
{
  int i,j,k,l,n;
  REAL *w;
  REAL_MATRIX v,wn,check1,check2;

  n=a.n;
  w=(REAL*)calloc(n,sizeof(REAL));
  v=CreateRealMatrix(n,n);
  wn=CreateRealMatrix(n,n);
  check1=CreateRealMatrix(n,n);
  check2=CreateRealMatrix(n,n);

  SingularValueDecomposition(a,n,n,w,v);

  for(i=0;i<n;i++)
  {
    if(fabs(w[i])<1e-8)
      wn.element[i][i]=0.0;
    else wn.element[i][i]=1.0/w[i];
  }

  for(k=0;k<n;k++)
  {
    for(l=0;l<n;l++)
    {
      check1.element[k][l]=0.0;
      for (j=0;j<n;j++)
        check1.element[k][l]+=(v.element[k][j]*wn.element[j][l]);
    }
  }

  for(k=0;k<n;k++)
  {
    for(l=0;l<n;l++)
    {
      check2.element[k][l]=0.0;
      for (j=0;j<n;j++)
        check2.element[k][l]+=(check1.element[k][j]*a.element[l][j]);
    }
  }
/*
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
    {
      if(i==j)
        wn.element[i][j]=w[i];
      else
        wn.element[i][j]=0.0;
    }

  for(k=0;k<n;k++)
  {
    for(l=0;l<n;l++)
    {
      check1.element[k][l]=0.0;
      for (j=0;j<n;j++)
        check1.element[k][l]+=(a.element[k][j]*wn.element[j][l]);
    }
  }

  for(k=0;k<n;k++)
  {
    for(l=0;l<n;l++)
    {
      check2.element[k][l]=0.0;
      for (j=0;j<n;j++)
        check2.element[k][l]+=(check1.element[k][j]*v.element[l][j]);
    }
  }

*/
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      a.element[i][j]=check2.element[i][j];

  DeleteRealMatrix(v);
  DeleteRealMatrix(check1);
  DeleteRealMatrix(check2);
  free(w);
}

void ScaleRealMatrix(REAL_MATRIX a,REAL fac)
{
  int i,j;

  for(i=0;i<a.m;i++)
    for(j=0;j<a.n;j++)
      a.element[i][j]*=fac;
}

void TransposeRealMatrix(REAL_MATRIX c,REAL_MATRIX a)
{
  int i,j;

  for(i=0;i<a.m;i++)
    for(j=0;j<a.n;j++)
      c.element[j][i]=a.element[i][j];

}

void MultiplyRealMatrix(REAL_MATRIX c,REAL_MATRIX a,REAL_MATRIX b)
{
  int i,j,k;

  if(a.n!=b.m) fprintf(stderr, "Matrix Multiplication error: %dx%d %dx%d\n",a.m,a.n,b.m,b.n);
  if(c.m!=a.m) fprintf(stderr, "Matrix Multiplication error: %dx%d %dx%d\n",a.m,a.n,b.m,b.n);
  if(c.n!=b.n) fprintf(stderr, "Matrix Multiplication error: %dx%d %dx%d\n",a.m,a.n,b.m,b.n);

  for(i=0;i<a.m;i++)
  {
    for(j=0;j<b.n;j++)
    {
      c.element[i][j]=0.0;
      for (k=0;k<a.n;k++)
        c.element[i][j]+=a.element[i][k]*b.element[k][j];
    }
  }
}

void MultiplyComplexMatrix(COMPLEX_MATRIX c,COMPLEX_MATRIX a,COMPLEX_MATRIX b)
{
  int i,j,k;

  if(a.n!=b.m) fprintf(stderr, "Matrix Multiplication error: %dx%d %dx%d\n",a.m,a.n,b.m,b.n);
  if(c.m!=a.m) fprintf(stderr, "Matrix Multiplication error: %dx%d %dx%d\n",a.m,a.n,b.m,b.n);
  if(c.n!=b.n) fprintf(stderr, "Matrix Multiplication error: %dx%d %dx%d\n",a.m,a.n,b.m,b.n);

  for(i=0;i<a.m;i++)
  {
    for(j=0;j<b.n;j++)
    {
      c.element[i][j].re=0.0;
      c.element[i][j].im=0.0;
      for (k=0;k<a.n;k++)
      {
        c.element[i][j].re+=a.element[i][k].re*b.element[k][j].re-a.element[i][k].im*b.element[k][j].im;
        c.element[i][j].im+=a.element[i][k].im*b.element[k][j].re+a.element[i][k].re*b.element[k][j].im;
      }
    }
  }
}


void MultiplyRealMatrixVector(REAL *c,REAL_MATRIX a,REAL *b)
{
  int i,k;

  for(i=0;i<a.m;i++)
  {
    c[i]=0.0;
    for (k=0;k<a.n;k++)
      c[i]+=a.element[i][k]*b[k];
  }
}


void TestEigenStystem(REAL_MATRIX H)
{
  int i,j,n;
  REAL_MATRIX Eigenvectors;
  REAL *Eigenvalues,*e;
  REAL *first_vector,*second_vector;

  fprintf(stderr, "testing eigensystem\n");

  n=H.n;

  Eigenvectors=CreateRealMatrix(n,n);
  Eigenvalues=(REAL*)calloc(n,sizeof(REAL));
  e=(REAL*)calloc(n,sizeof(REAL));
  first_vector=(REAL*)calloc(n,sizeof(REAL));
  second_vector=(REAL*)calloc(n,sizeof(REAL));

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      Eigenvectors.element[i][j]=H.element[i][j];

  tred2(Eigenvectors.element[0],n,Eigenvalues,e);
  tqli(Eigenvalues,e,n,Eigenvectors.element[0]);
  eigsrt(Eigenvalues,Eigenvectors.element[0],n);

  for(i=0;i<n;i++)
  {
    first_vector[i]=0.0;
    for(j=0;j<n;j++)
      first_vector[i]+=H.element[i][j]*Eigenvectors.element[0][j];
  }

  for(i=0;i<n;i++)
    second_vector[i]=Eigenvalues[0]*Eigenvectors.element[0][i];

  for(i=0;i<n;i++)
    fprintf(stderr, "%g %g <-> %g\n",first_vector[i],second_vector[i],fabs(first_vector[i]-second_vector[i]));

  free(Eigenvalues);
  free(e);
  free(first_vector);
  free(second_vector);
  DeleteRealMatrix(Eigenvectors);
}

void InverseMatrix3x3(REAL_MATRIX3x3 in,REAL_MATRIX3x3 *out)
{
  REAL_MATRIX a;

  a=CreateRealMatrix(3,3);

  a.element[0][0]=in.ax;  a.element[0][1]=in.bx;  a.element[0][2]=in.cx;
  a.element[1][0]=in.ay;  a.element[1][1]=in.by;  a.element[1][2]=in.cy;
  a.element[2][0]=in.az;  a.element[2][1]=in.bz;  a.element[2][2]=in.cz;

  SingularValueDecompositionMatrixInversion(a);

  out->ax=a.element[0][0]; out->bx=a.element[0][1]; out->cx=a.element[0][2];
  out->ay=a.element[1][0]; out->by=a.element[1][1]; out->cy=a.element[1][2];
  out->az=a.element[2][0]; out->bz=a.element[2][1]; out->cz=a.element[2][2];

  DeleteRealMatrix(a);
}

void InverseMatrix6x6(REAL_MATRIX6x6 in,REAL_MATRIX6x6 *out)
{
  REAL_MATRIX a;

  a=CreateRealMatrix(6,6);

  a.element[0][0]=in.C11;  a.element[0][1]=in.C12;  a.element[0][2]=in.C13;  a.element[0][3]=in.C14;  a.element[0][4]=in.C15;  a.element[0][5]=in.C16;
  a.element[1][0]=in.C21;  a.element[1][1]=in.C22;  a.element[1][2]=in.C23;  a.element[1][3]=in.C24;  a.element[1][4]=in.C25;  a.element[1][5]=in.C26;
  a.element[2][0]=in.C31;  a.element[2][1]=in.C32;  a.element[2][2]=in.C33;  a.element[2][3]=in.C34;  a.element[2][4]=in.C35;  a.element[2][5]=in.C36;
  a.element[3][0]=in.C41;  a.element[3][1]=in.C42;  a.element[3][2]=in.C43;  a.element[3][3]=in.C44;  a.element[3][4]=in.C45;  a.element[3][5]=in.C46;
  a.element[4][0]=in.C51;  a.element[4][1]=in.C52;  a.element[4][2]=in.C53;  a.element[4][3]=in.C54;  a.element[4][4]=in.C55;  a.element[4][5]=in.C56;
  a.element[5][0]=in.C61;  a.element[5][1]=in.C62;  a.element[5][2]=in.C63;  a.element[5][3]=in.C64;  a.element[5][4]=in.C65;  a.element[5][5]=in.C66;

  SingularValueDecompositionMatrixInversion(a);

  out->C11=a.element[0][0]; out->C12=a.element[0][1]; out->C13=a.element[0][2]; out->C14=a.element[0][3]; out->C15=a.element[0][4]; out->C16=a.element[0][5];
  out->C21=a.element[1][0]; out->C22=a.element[1][1]; out->C23=a.element[1][2]; out->C24=a.element[1][3]; out->C25=a.element[1][4]; out->C26=a.element[1][5];
  out->C31=a.element[2][0]; out->C32=a.element[2][1]; out->C33=a.element[2][2]; out->C34=a.element[2][3]; out->C35=a.element[2][4]; out->C36=a.element[2][5];
  out->C41=a.element[3][0]; out->C42=a.element[3][1]; out->C43=a.element[3][2]; out->C44=a.element[3][3]; out->C45=a.element[3][4]; out->C46=a.element[3][5];
  out->C51=a.element[4][0]; out->C52=a.element[4][1]; out->C53=a.element[4][2]; out->C54=a.element[4][3]; out->C55=a.element[4][4]; out->C56=a.element[4][5];
  out->C61=a.element[5][0]; out->C62=a.element[5][1]; out->C63=a.element[5][2]; out->C64=a.element[5][3]; out->C65=a.element[5][4]; out->C66=a.element[5][5];

  DeleteRealMatrix(a);
}
