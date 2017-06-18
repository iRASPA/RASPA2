/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'linear_equations.c' is part of RASPA-2.0

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
#include "linear_equations.h"

#define SQR(x) ((x)*(x))

int fzyfsy(int n,REAL * md,REAL *ud1, REAL *ud2,
           REAL *rs,REAL *x,REAL *cmd,REAL *cld_1,
           REAL *cld_2,REAL *cld_l2,REAL *cld_l1,
           REAL *bud_1,REAL *bud_2,REAL *brs_2,REAL *brs_1);
int fzyfsz(int n,REAL *md,REAL *ud1,REAL *ud2,REAL *cmd,
           REAL *cld_1,REAL *cld_2,REAL *cld_l2,REAL *cld_l1,
           REAL *bud_1,REAL *bud_2,REAL *brs_2, REAL *brs_1);
int fzyfsl(int n,REAL *rs,REAL *x,REAL *cmd,REAL *cld_1,
           REAL *cld_2,REAL *cld_l2,REAL *cld_l1,REAL *bud_1,
           REAL *bud_2,REAL *brs_2,REAL *brs_1);
int diag5dec(int   n,REAL *ld2, REAL *ld1,REAL *d,
             REAL *ud1,REAL *ud2);
int diag5sol(int n,REAL *ld2,REAL *ld1,REAL *d,
             REAL *ud1,REAL *ud2,REAL *b);
int diag5pddec(int n,REAL *d,REAL *ud1,REAL *ud2);
int diag5pdsol(int n,REAL *d,REAL *ud1,REAL *ud2,REAL *b);
int diag5pd(int mod,int n,REAL *d,REAL *ud1,
            REAL *ud2,REAL *b);
int diag5(int mod,int n,REAL *ld2,REAL *ld1, REAL *d,
          REAL *ud1,REAL *ud2,REAL *b);

// used lapack routines
int dgetrf_(int *m,int *n,REAL *a,int *lda,int *ipiv,int *info);
int dgetri_(int *n,REAL *a,int *lda,int *ipiv,REAL *work,int *lwork,int *info);
int dgesv_(int *n,int *nhrs,REAL *a,int *lda,int *ipiv,REAL *b,int *ldb,int *info);

int Inverse(REAL_MATRIX c)
{
  int arg1,arg2;
  int *ipiv;
  REAL *workspace;

  ipiv=(int*)calloc(c.n,sizeof(int));
  workspace=(REAL*)calloc(8*c.n,sizeof(REAL));
  arg1=c.n;
  arg2=8*c.n;
  //dgetrf_(&arg1,&arg1,c.element[0],&arg1,ipiv,&ok);
  //dgetri_(&arg1,c.element[0],&arg1,ipiv,workspace,&arg2,&ok);
  free(workspace);
  free(ipiv);
  return 0;
}

int SolveLinearSystem(REAL_FORTRAN_MATRIX *a,REAL_FORTRAN_MATRIX *b)
{
  int *ipiv;

  ipiv=(int*)calloc(a->n,sizeof(int));
  //dgesv_(&a->n,&b->n,a->element,&a->m,ipiv,b->element,&b->m,&ok);
  free(ipiv);
  return 0;
}

// Tridiagonal linear systems
int trdiag(int n,REAL *lower,REAL *diag,REAL *upper,REAL *b,int rep)
{
  int i;

  if(n < 2) return (1);                    /*  n at least 2          */

  if (lower == NULL || diag == NULL || upper == NULL ||
      b == NULL) return (1);

  if (rep == 0)                             /*  for rep = 0, determine*/
  {                                         /*  LU factorization      */
    for (i = 1; i < n; i++)
    {
      if (fabs(diag[i-1])<(REAL)DBL_EPSILON)        /*  if one   diag[i] = 0  */
        return (2);                         /*  we have no LU fact.   */
      lower[i]/=diag[i-1];
      diag[i]-=lower[i]*upper[i-1];
    }
  }

  if (fabs(diag[n-1])<(REAL)DBL_EPSILON) return (2);

  for (i = 1; i < n; i++)                   /*  update  b            */
    b[i] -= lower[i] * b[i-1];

  b[n-1] /= diag[n-1];                      /* backsubsitute         */
  for (i = n-2; i >= 0; i--)
    b[i] = ( b[i] - upper[i] * b[i+1] ) / diag[i];

  return (0);
}

int tzdiag(int n,REAL *lower,REAL *diag,REAL *upper,
           REAL *lowrow,REAL *ricol,REAL *b,int rep)
{
  REAL tmp;
  int i;

  if (n<3) return (1);
  if (lower==NULL||diag==NULL||upper==NULL||
      lowrow==NULL||ricol==NULL) return (1);

  if(rep==0)                             /*  if rep = 0,           */
  {                                         /*  decompose matrix      */
    lower[0]=upper[n-1]=0.0;

    if (fabs(diag[0])<(REAL)DBL_EPSILON) return (2);
                                          /* If a diagonal entry is   */
    tmp=1.0/diag[0];                  /* too close to zero, stop  */
    upper[0]*=tmp;
    ricol[0]*=tmp;

    for(i=1;i<n-2;i++)
    {
      diag[i]-=lower[i]*upper[i-1];
      if(fabs(diag[i])<(REAL)DBL_EPSILON) return (2);
      tmp=1.0/diag[i];
      upper[i]*=tmp;
      ricol[i]=-lower[i]*ricol[i-1]*tmp;
    }

    diag[n-2]-=lower[n-2]*upper[n-3];
    if (fabs(diag[n-2])<(REAL)DBL_EPSILON) return (2);

    for (i=1;i<n-2;i++)
      lowrow[i]=-lowrow[i-1]*upper[i-1];

    lower[n-1]-=lowrow[n-3]*upper[n-3];
    upper[n-2]=(upper[n-2]-lower[n-2]*ricol[n-3])/diag[n-2];

    for(tmp=0.0,i=0;i<n-2;i++)
      tmp-=lowrow[i]*ricol[i];
    diag[n-1]+=tmp-lower[n-1]*upper[n-2];

    if (fabs(diag[n-1])<(REAL)DBL_EPSILON) return (2);
  }  /* end if ( rep == 0 ) */

  b[0]/=diag[0];                                     /* update b   */
  for(i=1;i<n-1;i++)
    b[i]=(b[i]-b[i-1]*lower[i])/diag[i];

  for(tmp=0.0,i=0;i<n-2;i++)
    tmp-=lowrow[i]*b[i];

  b[n-1]=(b[n-1]+tmp-lower[n-1]*b[n-2])/diag[n-1];

  b[n-2]-=b[n-1]*upper[n-2];           /* back substitute         */
  for(i=n-3;i>=0;i--)
    b[i]-=upper[i]*b[i+1]+ricol[i]*b[n-1];

  return (0);
}

int fzyfsy(int n,REAL * md,REAL *ud1, REAL *ud2,
           REAL *rs,REAL *x,REAL *cmd,REAL *cld_1,
           REAL *cld_2,REAL *cld_l2,REAL *cld_l1,
           REAL *bud_1,REAL *bud_2,REAL *brs_2,REAL *brs_1)
{
  int error;

  if(n<6) return (2);

  // Factor system matrix into C * B for triangular matrices C, B
  error=fzyfsz(n,md,ud1,ud2,cmd,cld_1,cld_2,cld_l2,cld_l1,
          bud_1,bud_2,brs_2,brs_1);

  if (!error)  /* factorization without error */
    fzyfsl(n,rs,x,cmd,cld_1,cld_2,cld_l2,cld_l1,bud_1,bud_2,brs_2,brs_1);
  return (0);
}

int fzyfsz(int n,REAL *md,REAL *ud1,REAL *ud2,REAL *cmd,
           REAL *cld_1,REAL *cld_2,REAL *cld_l2,REAL *cld_l1,
           REAL *bud_1,REAL *bud_2,REAL *brs_2, REAL *brs_1)
{
  int i,j,k;
  REAL h_var_1,h_var_2,h_var_3;

  //  Factor A = C * B.
  cmd[1]=md[1];
  if(cmd[1]<=4.0*DBL_EPSILON) return (1);
  bud_1[1]=ud1[1]/cmd[1];
  brs_2[1]=ud2[n-1]/cmd[1];
  brs_1[1]=ud1[n]/cmd[1];
  cld_1[2]=ud1[1];
  cmd[2]=md[2]-cld_1[2]*bud_1[1];
  if(cmd[2]<=4.0*DBL_EPSILON) return (1);
  brs_2[2]=-(brs_2[1]*cld_1[2])/cmd[2];
  brs_1[2]=(ud2[n]-cld_1[2]*brs_1[1])/cmd[2];

  for(i=3;i<=n-2;++i)
  {
    j=i-2;
    k=i-1;
    cld_2[i]=ud2[i-2];
    bud_2[j]=ud2[j]/cmd[j];
    bud_1[k]=(ud1[k]-cld_1[k]*bud_2[j])/cmd[k];
    cld_1[i]=ud1[i-1]-cld_2[i]*bud_1[j];
    cmd[i]=md[i]-cld_1[i]*bud_1[k]-cld_2[i]*bud_2[j];
    if(cmd[i]<=4.0*DBL_EPSILON) return (1);
  }
  for(i=3;i<=n-4;++i)
    brs_2[i]=-(cld_2[i]*brs_2[i-2]+cld_1[i]*brs_2[i-1])/cmd[i];

  for(i=3;i<=n-3;++i)
    brs_1[i]=-(cld_2[i]*brs_1[i-2]+cld_1[i]*brs_1[i-1])/cmd[i];

  bud_2[n-3]=(ud2[n-3]-cld_1[n-3]*brs_2[n-4]-cld_2[n-3]*brs_2[n-5])/cmd[n-3];
  bud_2[n-2]=(ud2[n-2]-cld_1[n-2]*brs_1[n-3]-cld_2[n-2]*brs_1[n-4])/cmd[n-2];
  bud_1[n-2]=(ud1[n-2]-cld_1[n-2]*bud_2[n-3]-cld_2[n-2]*brs_2[n-4])/cmd[n-2];

  cld_l2[1]=ud2[n-1];
  cld_l2[2]=-cld_l2[1]*bud_1[1];
  for (i=3;i<=n-4;++i)
    cld_l2[i]=-(cld_l2[i-2]*bud_2[i-2]+cld_l2[i-1]*bud_1[i-1]);

  cld_l1[1]=ud1[n];
  cld_l1[2]=ud2[n]-cld_l1[1]*bud_1[1];
  for (i=3;i<=n-3;++i)
    cld_l1[i]=-(cld_l1[i-2]*bud_2[i-2]+cld_l1[i-1]*bud_1[i-1]);

  cld_2[n-1]=ud2[n-3]-(cld_l2[n-5]*bud_2[n-5]+cld_l2[n-4]*bud_1[n-4]);
  cld_2[n]=ud2[n-2]-(cld_l1[n-4]*bud_2[n-4]+ cld_l1[n-3]*bud_1[n-3]);
  cld_1[n-1]=ud1[n-2]-(cld_l2[n-4]*bud_2[n-4]+cld_2[n-1]*bud_1[n-3]);

  h_var_1=h_var_2=h_var_3=0.0;
  for(i=1;i<=n-4;++i)
  {
    h_var_1+=cld_l1[i]*brs_2[i];
    h_var_2+=cld_l2[i]*brs_2[i];
    h_var_3+=cld_l2[i]*brs_1[i];
  }

  cld_1[n]=ud1[n-1]-h_var_1-cld_l1[n-3]*bud_2[n-3]-cld_2[n]*bud_1[n-2];

  cmd[n-1]=md[n-1]-h_var_2-cld_2[n-1]*bud_2[n-3]-cld_1[n-1]*bud_1[n-2];
  if(cmd[n-1]<=4.0*DBL_EPSILON) return (1);

  bud_1[n-1]=(ud1[n-1]-h_var_3-cld_2[n-1]*brs_1[n-3]
              -cld_1[n-1]*bud_2[n-2])/cmd[n-1];

  for(h_var_1=0.0,i=1;i<=n-3;++i)
    h_var_1+=cld_l1[i]*brs_1[i];
  cmd[n]=md[n]-h_var_1-cld_2[n]*bud_2[n-2]-cld_1[n]*bud_1[n-1];
  if(cmd[n]<=4.0*DBL_EPSILON) return (1);
  return (0);
}

int fzyfsl(int n,REAL *rs,REAL *x,REAL *cmd,REAL *cld_1,
           REAL *cld_2,REAL *cld_l2,REAL *cld_l1,REAL *bud_1,
           REAL *bud_2,REAL *brs_2,REAL *brs_1)
{
  int i;
  REAL h_var_1;

  //  Solve system by updating right hand side and backsubstitution
  x[1]=rs[1]/cmd[1];                   /*  update right hand side */
  x[2]=(rs[2]-x[1]*cld_1[2])/cmd[2];
  for(i=3;i<=n-2;++i)
    x[i]=(rs[i]-x[i-2]*cld_2[i]-x[i-1]*cld_1[i])/cmd[i];

  for(h_var_1=0.0,i=1;i<=n-4;++i)
    h_var_1+=x[i]*cld_l2[i];
  x[n-1]=(rs[n-1]-h_var_1-x[n-3]*cld_2[n-1]
            - x[n-2]*cld_1[n-1])/cmd[n-1];

  for(h_var_1=0.0,i=1;i<=n-3;++i)
    h_var_1+=x[i]*cld_l1[i];
  x[n]=(rs[n]-h_var_1-x[n-2]*cld_2[n]-x[n-1]*cld_1[n])/cmd[n];

  x[n-1]-=bud_1[n-1]*x[n];       /* back substitute  */
  x[n-2]-=(bud_1[n-2]*x[n-1]+bud_2[n-2]*x[n]);
  x[n-3]-=(bud_1[n-3]*x[n-2]+bud_2[n-3]*x[n-1]+brs_1[n-3]*x[n]);

  for(i=n-4;i>=1;--i)
    x[i]-=(bud_1[i]*x[i+1]+bud_2[i]*x[i+2]+brs_2[i]*x[n-1]+brs_1[i]*x[n]);
  return 0;
}

// 5 diagonal linear systems
int diag5(int mod,int n,REAL *ld2,REAL *ld1, REAL *d,
          REAL *ud1,REAL *ud2,REAL *b)
{
  int rc;

  switch (mod)
  {
    case 0:
      rc=diag5dec(n,ld2,ld1,d,ud1,ud2);
      if (rc==0)
        return (diag5sol(n,ld2,ld1,d,ud1,ud2,b));
      else
        return (rc);
    case 1:
      return (diag5dec(n,ld2,ld1,d,ud1,ud2));
    case 2:
      return (diag5sol(n,ld2,ld1,d,ud1,ud2,b));
  }
  return (3);
}


int diag5dec(int   n,REAL *ld2, REAL *ld1,REAL *d,
             REAL *ud1,REAL *ud2)
{
  register int i;
  REAL row,D,ud1i,ud2i;

  if (n<3) return (1);

  if (ld2==NULL||ld1==NULL||d==NULL||ud1==NULL||
      ud2==NULL) return (1);

  row=fabs(d[0])+fabs(ud1[0])+fabs(ud2[0]);
  if (row==0.0) return 2;
  D=1.0/row;
  if (fabs(d[0])*D<=4.0*DBL_EPSILON) return 2;

  ud1[0]/=d[0];
  ud2[0]/=d[0];
  row=fabs(ld1[1])+fabs(d[1])+fabs(ud1[1])+fabs(ud2[1]);
  if (row==0.0) return 2;
  D=1.0/row;

  d[1]-=ld1[1]*ud1[0];

  if (fabs(d[1])*D<=4.0*DBL_EPSILON) return 2;

  ud1[1]=(ud1[1]-ld1[1]*ud2[0])/d[1];
  ud2[1]/=d[1];

  ud1i=ud1[2];
  ud2i=ud2[2];
  for(i=2;i<n;i++)
  {
    row=fabs(ld2[i])+fabs(ld1[i])+fabs(d[i])+fabs(ud1i)+fabs(ud2i);
    if (row==0.0) return 2;
    D=1.0/row;
    ld1[i]-=ld2[i]*ud1[i-2];
    d[i]-=(ld2[i]*ud2[i-2]+ld1[i]*ud1[i-1]);

    if(fabs(d[i])*D<=4.0*DBL_EPSILON) return 2;

    if (i<n-1)
      ud1[i]=(ud1[i]-ld1[i]*ud2[i-1])/d[i];

    if (i<n-2)
      ud2[i]/=d[i],
      ud1i=ud1[i+1];
    else
      ud1i=0.0;

    if (i<n-3)
      ud2i=ud2[i+1];
    else
      ud2i=0.0;
  }
  return (0);
}


int diag5sol(int n,REAL *ld2,REAL *ld1,REAL *d,
             REAL *ud1,REAL *ud2,REAL *b)
{
  register int i;

  if (n<3) return (1);

  if (ld2==NULL||ld1==NULL||d==NULL||ud1==NULL||ud2==NULL||b==NULL) return (1);

  if (fabs(d[0])<(REAL)DBL_EPSILON) return (2);

  b[0]/=d[0];

  if (fabs(d[1])<(REAL)DBL_EPSILON) return (2);
  b[1]=(b[1]-ld1[1]*b[0])/d[1];

  for(i=2;i<n;i++)
  {
    if (fabs(d[i])<(REAL)DBL_EPSILON) return (2);
    b[i]=(b[i]-ld2[i]*b[i-2]-ld1[i]*b[i-1])/d[i];
  }
  b[n-2]-=ud1[n-2]*b[n-1];

  for(i=n-3;i>=0;i--)
    b[i]-=(ud1[i]*b[i+1]+ud2[i]*b[i+2]);

  return (0);
}

int diag5pd(int mod,int n,REAL *d,REAL *ud1,
            REAL *ud2,REAL *b)
{
  int rc;

  if(n<3) return (1);

  switch (mod)
  {
    case 0:
      rc = diag5pddec (n, d, ud1, ud2);
      if (rc == 0)
        return (diag5pdsol (n, d, ud1, ud2, b));
      else
        return (rc);
    case 1:
      return (diag5pddec(n,d,ud1,ud2));

    case 2:
      return (diag5pdsol(n,d,ud1,ud2,b));
  }
  return (3);
}


int diag5pddec(int n,REAL *d,REAL *ud1,REAL *ud2)
{
  register int i;
  REAL e_1,e_2,tmp,sum;

  if (n<3) return (1);

  if (d==NULL||ud1==NULL||ud2==NULL) return (1);

  ud1[n-1]=0.0;
  ud2[n-1]=0.0;
  ud2[n-2]=0.0;

  sum=fabs(d[0])+fabs(ud1[0])+fabs(ud2[0]);
  if(sum==0.0) return (3);

  if(d[0]<sum*(REAL)DBL_EPSILON) return (2);

  tmp=ud1[0];
  ud1[0]/=d[0];
  e_2=ud2[0];
  ud2[0]/=d[0];

  sum=fabs(tmp)+fabs(d[1])+fabs(ud1[1])+fabs(ud2[1]);
  if(sum==0.0) return (2);

  d[1]-=tmp*ud1[0];

  if(fabs(d[1])<sum*(REAL)DBL_EPSILON) return (2);

  tmp=ud1[1];
  ud1[1]=(ud1[1]-e_2*ud1[0])/d[1];
  e_1=ud2[1];
  ud2[1]/=d[1];

  for(i=2;i<n;i++)
  {
    sum=fabs(e_2)+fabs(tmp)+fabs(d[i])+fabs(ud1[i])+fabs(ud2[i]);
    if(sum==0.0) return (2);

    d[i]-=e_2*ud2[i-2]+d[i-1]*SQR(ud1[i-1]);

    if(fabs(d[i])<sum*(REAL)DBL_EPSILON) return (2);

    if (i<n-1)
    {
      tmp=ud1[i];
      ud1[i]=(ud1[i]-e_1*ud1[i-1])/d[i];
    }

    if(i<n-2)
    {
      e_2=e_1;
      e_1=ud2[i];
      ud2[i]/=d[i];
    }
  }
  return (0);
}

int diag5pdsol(int n,REAL *d,REAL *ud1,REAL *ud2,REAL *b)
{
  register int i;

  if(n<3) return (1);

  if(d==NULL||ud1==NULL||ud2==NULL||b==NULL) return (1);

  if((REAL)fabs(d[1])<(REAL)DBL_EPSILON) return (2);
  b[1]-=ud1[0]*b[0];

  for(i=2;i<n;i++)
    b[i]-=ud1[i-1]*b[i-1]+ud2[i-2]*b[i-2];

  for(i=0;i<n;i++)
  {
    if(fabs(d[i])<(REAL)DBL_EPSILON) return (2);
    b[i]/=d[i];
  }

  b[n-2]-=ud1[n-2]*b[n-1];

  for(i=n-3;i>=0;i--)
    b[i]-=(ud1[i]*b[i+1]+ud2[i]*b[i+2]);

  return (0);
}
