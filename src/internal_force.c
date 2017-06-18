/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'internal_force.c' is part of RASPA-2.0

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "potentials.h"
#include "molecule.h"
#include "cbmc.h"
#include "utils.h"
#include "internal_force.h"
#include "ewald.h"
#include "matrix.h"
#include "spectra.h"
#include "minimization.h"
#include "rigid.h"

REAL ReturnBondDistance(VECTOR posA,VECTOR posB)
{
  VECTOR dr;
  REAL rr;

  dr.x=posA.x-posB.x;
  dr.y=posA.y-posB.y;
  dr.z=posA.z-posB.z;
  dr=ApplyBoundaryCondition(dr);
  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  return sqrt(rr);
}

REAL ReturnBendAngle(VECTOR posA,VECTOR posB,VECTOR posC)
{
  REAL CosTheta;
  VECTOR e12,e32;
  REAL r12,r32;

  e12.x=posA.x-posB.x;
  e12.y=posA.y-posB.y;
  e12.z=posA.z-posB.z;
  e12=ApplyBoundaryCondition(e12);
  r12=sqrt(SQR(e12.x)+SQR(e12.y)+SQR(e12.z));
  e12.x/=r12;
  e12.y/=r12;
  e12.z/=r12;

  e32.x=posC.x-posB.x;
  e32.y=posC.y-posB.y;
  e32.z=posC.z-posB.z;
  e32=ApplyBoundaryCondition(e32);
  r32=sqrt(SQR(e32.x)+SQR(e32.y)+SQR(e32.z));
  e32.x/=r32;
  e32.y/=r32;
  e32.z/=r32;

  CosTheta=(e12.x*e32.x+e12.y*e32.y+e12.z*e32.z);
  return acos(CosTheta);
}

REAL ReturnDihedralAngle(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD)
{
  REAL rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi;
  VECTOR Pb,Pc;
  VECTOR vec_u,vec_v,vec_w;
  REAL u,w,v;

  vec_u.x=posB.x-posC.x;
  vec_u.y=posB.y-posC.y;
  vec_u.z=posB.z-posC.z;
  vec_u=ApplyBoundaryCondition(vec_u);
  u=sqrt(SQR(vec_u.x)+SQR(vec_u.y)+SQR(vec_u.z));

  vec_v.x=posD.x-posC.x;
  vec_v.y=posD.y-posC.y;
  vec_v.z=posD.z-posC.z;
  vec_v=ApplyBoundaryCondition(vec_v);
  v=sqrt(SQR(vec_v.x)+SQR(vec_v.y)+SQR(vec_v.z));

  vec_w.x=posA.x-posB.x;
  vec_w.y=posA.y-posB.y;
  vec_w.z=posA.z-posB.z;
  vec_w=ApplyBoundaryCondition(vec_w);
  w=sqrt(SQR(vec_w.x)+SQR(vec_w.y)+SQR(vec_w.z));

  Dab.x=posA.x-posB.x;
  Dab.y=posA.y-posB.y;
  Dab.z=posA.z-posB.z;
  Dab=ApplyBoundaryCondition(Dab);

  Dcb.x=posC.x-posB.x;
  Dcb.y=posC.y-posB.y;
  Dcb.z=posC.z-posB.z;
  Dcb=ApplyBoundaryCondition(Dcb);
  rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
  Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

  Ddc.x=posD.x-posC.x;
  Ddc.y=posD.y-posC.y;
  Ddc.z=posD.z-posC.z;
  Ddc=ApplyBoundaryCondition(Ddc);

  dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
  dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;


  dr.x=Dab.x-dot_ab*Dcb.x;
  dr.y=Dab.y-dot_ab*Dcb.y;
  dr.z=Dab.z-dot_ab*Dcb.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  dr.x/=r; dr.y/=r; dr.z/=r;

  ds.x=Ddc.x-dot_cd*Dcb.x;
  ds.y=Ddc.y-dot_cd*Dcb.y;
  ds.z=Ddc.z-dot_cd*Dcb.z;
  s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
  ds.x/=s; ds.y/=s; ds.z/=s;

  // compute Cos(Phi)
  // Phi is defined in protein convention Phi(trans)=Pi
  CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

  // Ensure CosPhi is between -1 and 1.
  CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);

  // potential defined in terms of 'phi' and therefore contains a singularity
  // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
  // same direction as Rbc, and negative otherwise
  Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
  Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
  Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
  Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
  Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
  Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
  sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
      +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
  Phi=SIGN(acos(CosPhi),sign);
  return Phi;
}


REAL ReturnInversionBendAngle(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD)
{
  REAL r;
  VECTOR cross_product_dc,cross_product_ad,cross_product_ca;
  REAL dot_product;
  VECTOR Dab,Dcb,Ddb;
  REAL ChiA,ChiC,ChiD;

  Dab.x=posA.x-posB.x;
  Dab.y=posA.y-posB.y;
  Dab.z=posA.z-posB.z;
  Dab=ApplyBoundaryCondition(Dab);
  r=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));
  Dab.x/=r;
  Dab.y/=r;
  Dab.z/=r;

  Dcb.x=posC.x-posB.x;
  Dcb.y=posC.y-posB.y;
  Dcb.z=posC.z-posB.z;
  Dcb=ApplyBoundaryCondition(Dcb);
  r=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
  Dcb.x/=r;
  Dcb.y/=r;
  Dcb.z/=r;

  Ddb.x=posD.x-posB.x;
  Ddb.y=posD.y-posB.y;
  Ddb.z=posD.z-posB.z;
  Ddb=ApplyBoundaryCondition(Ddb);
  r=sqrt(SQR(Ddb.x)+SQR(Ddb.y)+SQR(Ddb.z));
  Ddb.x/=r;
  Ddb.y/=r;
  Ddb.z/=r;

  cross_product_dc.x=Ddb.y*Dcb.z-Ddb.z*Dcb.y;
  cross_product_dc.y=Ddb.z*Dcb.x-Ddb.x*Dcb.z;
  cross_product_dc.z=Ddb.x*Dcb.y-Ddb.y*Dcb.x;
  dot_product=cross_product_dc.x*Dab.x+cross_product_dc.y*Dab.y+cross_product_dc.z*Dab.z;
  ChiA=asin(dot_product/sin(ReturnBendAngle(posC,posB,posD)));

  cross_product_ad.x=Dab.y*Ddb.z-Dab.z*Ddb.y;
  cross_product_ad.y=Dab.z*Ddb.x-Dab.x*Ddb.z;
  cross_product_ad.z=Dab.x*Ddb.y-Dab.y*Ddb.x;
  dot_product=cross_product_ad.x*Dcb.x+cross_product_ad.y*Dcb.y+cross_product_ad.z*Dcb.z;
  ChiC=asin(dot_product/sin(ReturnBendAngle(posA,posB,posD)));

  cross_product_ca.x=Dcb.y*Dab.z-Dcb.z*Dab.y;
  cross_product_ca.y=Dcb.z*Dab.x-Dcb.x*Dab.z;
  cross_product_ca.z=Dcb.x*Dab.y-Dcb.y*Dab.x;
  dot_product=cross_product_ca.x*Ddb.x+cross_product_ca.y*Ddb.y+cross_product_ca.z*Ddb.z;
  ChiD=asin(dot_product/sin(ReturnBendAngle(posA,posB,posC)));

  //return (ChiA+ChiC+ChiD)/3.0;
  return ChiA;
}

REAL ReturnOutOfPlaneDistance(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD)
{
  return 0.0;
}

// the Wilson vectors are the conversion from internal bond-distance derivative to Cartesian derivatives: d\phi\dr_{1,2}
void ReturnWilsonVectorsBond(VECTOR posA,VECTOR posB,VECTOR *wa,VECTOR *wb)
{
  VECTOR e12;
  REAL r12;

  e12.x=posA.x-posB.x;
  e12.y=posA.y-posB.y;
  e12.z=posA.z-posB.z;
  e12=ApplyBoundaryCondition(e12);
  r12=sqrt(SQR(e12.x)+SQR(e12.y)+SQR(e12.z));
  e12.x/=r12;
  e12.y/=r12;
  e12.z/=r12;

  wa->x=e12.x;
  wa->y=e12.y;
  wa->z=e12.z;
  wb->x=-e12.x;
  wb->y=-e12.y;
  wb->z=-e12.z;
}


// the Wilson vectors are the conversion from internal theta-angle derivative to Cartesian derivatives: d\phi\dr_{1,2,3}
void ReturnWilsonVectorsBend(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR *wa,VECTOR *wb,VECTOR *wc)
{
  VECTOR dtA,dtB,dtC;
  VECTOR e12,e32;
  REAL r12,r32;
  REAL CosTheta,SinTheta;

  e12.x=posA.x-posB.x;
  e12.y=posA.y-posB.y;
  e12.z=posA.z-posB.z;
  e12=ApplyBoundaryCondition(e12);
  r12=sqrt(SQR(e12.x)+SQR(e12.y)+SQR(e12.z));
  e12.x/=r12;
  e12.y/=r12;
  e12.z/=r12;

  e32.x=posC.x-posB.x;
  e32.y=posC.y-posB.y;
  e32.z=posC.z-posB.z;
  e32=ApplyBoundaryCondition(e32);
  r32=sqrt(SQR(e32.x)+SQR(e32.y)+SQR(e32.z));
  e32.x/=r32;
  e32.y/=r32;
  e32.z/=r32;

  CosTheta=(e12.x*e32.x+e12.y*e32.y+e12.z*e32.z);
  SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));

  // Calculate the components of the derivatives.
  dtA.x=(CosTheta*e12.x-e32.x)/(r12*SinTheta);
  dtA.y=(CosTheta*e12.y-e32.y)/(r12*SinTheta);
  dtA.z=(CosTheta*e12.z-e32.z)/(r12*SinTheta);

  dtC.x=(CosTheta*e32.x-e12.x)/(r32*SinTheta);
  dtC.y=(CosTheta*e32.y-e12.y)/(r32*SinTheta);
  dtC.z=(CosTheta*e32.z-e12.z)/(r32*SinTheta);

  dtB.x=-(dtA.x+dtC.x);
  dtB.y=-(dtA.y+dtC.y);
  dtB.z=-(dtA.z+dtC.z);

  wa->x=dtA.x;
  wa->y=dtA.y;
  wa->z=dtA.z;

  wb->x=dtB.x;
  wb->y=dtB.y;
  wb->z=dtB.z;

  wc->x=dtC.x;
  wc->y=dtC.y;
  wc->z=dtC.z;
}

// the Wilson vectors are the conversion from internal phi-angle derivative to Cartesian derivatives: d\phi\dr_{1,2,3,4}
void ReturnWilsonVectorsTorsion(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR *wa,VECTOR *wb,VECTOR *wc,VECTOR *wd)
{
  REAL CosTheta1,CosTheta2;
  REAL SinTheta1,SinTheta2;
  VECTOR e12,e23,e32,e34,e43;
  REAL r12,r23,r34,temp1,temp2;

  e12.x=posB.x-posA.x;
  e12.y=posB.y-posA.y;
  e12.z=posB.z-posA.z;
  e12=ApplyBoundaryCondition(e12);
  r12=sqrt(SQR(e12.x)+SQR(e12.y)+SQR(e12.z));
  e12.x/=r12;
  e12.y/=r12;
  e12.z/=r12;

  e23.x=posC.x-posB.x;
  e23.y=posC.y-posB.y;
  e23.z=posC.z-posB.z;
  e23=ApplyBoundaryCondition(e23);
  r23=sqrt(SQR(e23.x)+SQR(e23.y)+SQR(e23.z));
  e23.x/=r23;
  e23.y/=r23;
  e23.z/=r23;
  e32.x=-e23.x;
  e32.y=-e23.y;
  e32.z=-e23.z;

  e34.x=posD.x-posC.x;
  e34.y=posD.y-posC.y;
  e34.z=posD.z-posC.z;
  e34=ApplyBoundaryCondition(e34);
  r34=sqrt(SQR(e34.x)+SQR(e34.y)+SQR(e34.z));
  e34.x/=r34;
  e34.y/=r34;
  e34.z/=r34;
  e43.x=-e34.x;
  e43.y=-e34.y;
  e43.z=-e34.z;

  CosTheta1=(e12.x*e32.x+e12.y*e32.y+e12.z*e32.z);
  SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

  CosTheta2=(e23.x*e43.x+e23.y*e43.y+e23.z*e43.z);
  SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

  temp1=1.0/(r12*SQR(SinTheta1));
  wa->x=-temp1*(e12.y*e23.z-e12.z*e23.y);
  wa->y=-temp1*(e12.z*e23.x-e12.x*e23.z);
  wa->z=-temp1*(e12.x*e23.y-e12.y*e23.x);

  temp1=(r23-r12*CosTheta1)/(r23*r12*SQR(SinTheta1));
  temp2=CosTheta2/(r23*SQR(SinTheta2));
  wb->x=temp1*(e12.y*e23.z-e12.z*e23.y)+temp2*(e43.y*e32.z-e43.z*e32.y);
  wb->y=temp1*(e12.z*e23.x-e12.x*e23.z)+temp2*(e43.z*e32.x-e43.x*e32.z);
  wb->z=temp1*(e12.x*e23.y-e12.y*e23.x)+temp2*(e43.x*e32.y-e43.y*e32.x);

  temp1=(r23-r34*CosTheta2)/(r23*r34*SQR(SinTheta2));
  temp2=CosTheta1/(r23*SQR(SinTheta1));
  wc->x=temp1*(e43.y*e32.z-e43.z*e32.y)+temp2*(e12.y*e23.z-e12.z*e23.y);
  wc->y=temp1*(e43.z*e32.x-e43.x*e32.z)+temp2*(e12.z*e23.x-e12.x*e23.z);
  wc->z=temp1*(e43.x*e32.y-e43.y*e32.x)+temp2*(e12.x*e23.y-e12.y*e23.x);

  temp1=1.0/(r34*SQR(SinTheta2));
  wd->x=-temp1*(e43.y*e32.z-e43.z*e32.y);
  wd->y=-temp1*(e43.z*e32.x-e43.x*e32.z);
  wd->z=-temp1*(e43.x*e32.y-e43.y*e32.x);
}

void ReturnWilsonVectorsInversionBend(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR *wa,VECTOR *wb,VECTOR *wc,VECTOR *wd)
{
  REAL rab,rbc,rbd;
  VECTOR cross_product_dc,cross_product_ad,cross_product_ca;
  REAL dot_product;
  REAL ThetaA,SinThetaA,CosThetaA,ChiA,CosChiA,TanChiA;
  REAL ThetaC,SinThetaC,CosThetaC,ChiC,CosChiC,TanChiC;
  REAL ThetaD,SinThetaD,CosThetaD,ChiD,CosChiD,TanChiD;
  VECTOR Dab,Dcb,Ddb;

  Dab.x=posA.x-posB.x;
  Dab.y=posA.y-posB.y;
  Dab.z=posA.z-posB.z;
  Dab=ApplyBoundaryCondition(Dab);
  rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));
  Dab.x/=rab;
  Dab.y/=rab;
  Dab.z/=rab;

  Dcb.x=posC.x-posB.x;
  Dcb.y=posC.y-posB.y;
  Dcb.z=posC.z-posB.z;
  Dcb=ApplyBoundaryCondition(Dcb);
  rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
  Dcb.x/=rbc;
  Dcb.y/=rbc;
  Dcb.z/=rbc;

  Ddb.x=posD.x-posB.x;
  Ddb.y=posD.y-posB.y;
  Ddb.z=posD.z-posB.z;
  Ddb=ApplyBoundaryCondition(Ddb);
  rbd=sqrt(SQR(Ddb.x)+SQR(Ddb.y)+SQR(Ddb.z));
  Ddb.x/=rbd;
  Ddb.y/=rbd;
  Ddb.z/=rbd;

  ThetaA=ReturnBendAngle(posC,posB,posD);
  CosThetaA=cos(ReturnBendAngle(posC,posB,posD));
  SinThetaA=sin(ReturnBendAngle(posC,posB,posD));

  ThetaC=ReturnBendAngle(posA,posB,posD);
  CosThetaC=cos(ReturnBendAngle(posA,posB,posD));
  SinThetaC=sin(ReturnBendAngle(posA,posB,posD));

  ThetaD=ReturnBendAngle(posA,posB,posC);
  CosThetaD=cos(ReturnBendAngle(posA,posB,posC));
  SinThetaD=sin(ReturnBendAngle(posA,posB,posC));

  cross_product_dc.x=Ddb.y*Dcb.z-Ddb.z*Dcb.y;
  cross_product_dc.y=Ddb.z*Dcb.x-Ddb.x*Dcb.z;
  cross_product_dc.z=Ddb.x*Dcb.y-Ddb.y*Dcb.x;
  dot_product=cross_product_dc.x*Dab.x+cross_product_dc.y*Dab.y+cross_product_dc.z*Dab.z;
  ChiA=asin(dot_product/SinThetaA);
  CosChiA=cos(ChiA);
  TanChiA=tan(ChiA);

  cross_product_ad.x=Dab.y*Ddb.z-Dab.z*Ddb.y;
  cross_product_ad.y=Dab.z*Ddb.x-Dab.x*Ddb.z;
  cross_product_ad.z=Dab.x*Ddb.y-Dab.y*Ddb.x;
  dot_product=cross_product_ad.x*Dcb.x+cross_product_ad.y*Dcb.y+cross_product_ad.z*Dcb.z;
  ChiC=asin(dot_product/sin(ReturnBendAngle(posA,posB,posD)));
  CosChiC=cos(ChiC);
  TanChiC=tan(ChiC);

  cross_product_ca.x=Dcb.y*Dab.z-Dcb.z*Dab.y;
  cross_product_ca.y=Dcb.z*Dab.x-Dcb.x*Dab.z;
  cross_product_ca.z=Dcb.x*Dab.y-Dcb.y*Dab.x;
  dot_product=cross_product_ca.x*Ddb.x+cross_product_ca.y*Ddb.y+cross_product_ca.z*Ddb.z;
  ChiD=asin(dot_product/sin(ReturnBendAngle(posA,posB,posC)));
  CosChiD=cos(ChiD);
  TanChiD=tan(ChiD);

  wa->x=(cross_product_dc.x/(CosChiA*SinThetaA)-Dab.x*TanChiA)/(3.0*rab);
  wa->y=(cross_product_dc.y/(CosChiA*SinThetaA)-Dab.y*TanChiA)/(3.0*rab);
  wa->z=(cross_product_dc.z/(CosChiA*SinThetaA)-Dab.z*TanChiA)/(3.0*rab);

  wc->x=(cross_product_ad.x/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Dcb.x-Ddb.x*CosThetaA))/(3.0*rbc);
  wc->y=(cross_product_ad.y/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Dcb.y-Ddb.y*CosThetaA))/(3.0*rbc);
  wc->z=(cross_product_ad.z/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Dcb.z-Ddb.z*CosThetaA))/(3.0*rbc);

  wd->x=(cross_product_ca.x/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Ddb.x-Dcb.x*CosThetaA))/(3.0*rbd);
  wd->y=(cross_product_ca.y/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Ddb.y-Dcb.y*CosThetaA))/(3.0*rbd);
  wd->z=(cross_product_ca.z/(CosChiA*SinThetaA)-(TanChiA/(SQR(SinThetaA)))*(Ddb.z-Dcb.z*CosThetaA))/(3.0*rbd);


  wc->x+=(cross_product_ad.x/(CosChiC*SinThetaC)-Dcb.x*TanChiC)/(3.0*rbc);
  wc->y+=(cross_product_ad.y/(CosChiC*SinThetaC)-Dcb.y*TanChiC)/(3.0*rbc);
  wc->z+=(cross_product_ad.z/(CosChiC*SinThetaC)-Dcb.z*TanChiC)/(3.0*rbc);

  wd->x+=(cross_product_ca.x/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Ddb.x-Dab.x*CosThetaC))/(3.0*rbd);
  wd->y+=(cross_product_ca.y/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Ddb.y-Dab.y*CosThetaC))/(3.0*rbd);
  wd->z+=(cross_product_ca.z/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Ddb.z-Dab.z*CosThetaC))/(3.0*rbd);

  wa->x+=(cross_product_dc.x/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Dab.x-Ddb.x*CosThetaC))/(3.0*rab);
  wa->y+=(cross_product_dc.y/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Dab.y-Ddb.y*CosThetaC))/(3.0*rab);
  wa->z+=(cross_product_dc.z/(CosChiC*SinThetaC)-(TanChiC/(SQR(SinThetaC)))*(Dab.z-Ddb.z*CosThetaC))/(3.0*rab);


  wd->x+=(cross_product_ca.x/(CosChiD*SinThetaD)-Ddb.x*TanChiD)/(3.0*rbd);
  wd->y+=(cross_product_ca.y/(CosChiD*SinThetaD)-Ddb.y*TanChiD)/(3.0*rbd);
  wd->z+=(cross_product_ca.z/(CosChiD*SinThetaD)-Ddb.z*TanChiD)/(3.0*rbd);

  wa->x+=(cross_product_dc.x/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Dab.x-Dcb.x*CosThetaD))/(3.0*rab);
  wa->y+=(cross_product_dc.y/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Dab.y-Dcb.y*CosThetaD))/(3.0*rab);
  wa->z+=(cross_product_dc.z/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Dab.z-Dcb.z*CosThetaD))/(3.0*rab);

  wc->x+=(cross_product_ad.x/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Dcb.x-Dab.x*CosThetaD))/(3.0*rbc);
  wc->y+=(cross_product_ad.y/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Dcb.y-Dab.y*CosThetaD))/(3.0*rbc);
  wc->z+=(cross_product_ad.z/(CosChiD*SinThetaD)-(TanChiD/(SQR(SinThetaD)))*(Dcb.z-Dab.z*CosThetaD))/(3.0*rbc);




  wb->x=-(wa->x+wc->x+wd->x);
  wb->y=-(wa->y+wc->y+wd->y);
  wb->z=-(wa->z+wc->z+wd->z);
}

// Calculate Bondforce
void CalculateAdsorbateBondForce(int m)
{
  int i,Type,NumberOfBonds,A,B;
  REAL U,DF,r,rr,temp,temp2,exp_term,r1;
  REAL *parms;
  POINT posA,posB;
  VECTOR dr,f;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBonds=Components[Type].NumberOfBonds;
  for(i=0;i<NumberOfBonds;i++)
  {
    A=Components[Type].Bonds[i].A;
    B=Components[Type].Bonds[i].B;

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    parms=(REAL*)&Components[Type].BondArguments[i];

    switch(Components[Type].BondType[i])
    {
      case HARMONIC_BOND:
        // 0.5*p0*SQR(r-p1);
        // ===============================================
        // p_0/k_B [K/A^2]   force constant
        // p_1     [A]       reference bond distance
        U=0.5*parms[0]*SQR(r-parms[1]);
        DF=parms[0]*(r-parms[1])/r;
        break;
      case CORE_SHELL_SPRING:
        U=0.5*parms[0]*SQR(r);
        DF=parms[0];
        break;
      case MORSE_BOND:
        // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
        // ===============================================
        // p_0/k_B [K]       force constant
        // p_1     [A^-1]    parameter
        // p_2     [A]       reference bond distance
        temp=exp(parms[1]*(parms[2]-r));
        U=parms[0]*(SQR(1.0-temp)-1.0);
        DF=2.0*parms[0]*parms[1]*(1.0-temp)*temp/r;
        break;
      case LJ_12_6_BOND:
        // A/r_ij^12-B/r_ij^6
        // ===============================================
        // p_0/k_B [K A^12]
        // p_1/k_B [K A^6]
        temp=CUBE(1.0/rr);
        U=parms[0]*SQR(temp)-parms[1]*temp;
        DF=6.0*(parms[1]*temp-2.0*parms[0]*SQR(temp))/rr;
        break;
      case LENNARD_JONES_BOND:
        // 4*p_0*((p_1/r)^12-(p_1/r)^6)
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A]
        temp=CUBE(parms[1]/rr);
        U=4.0*parms[0]*(temp*(temp-1.0));
        DF=24.0*parms[0]*(temp*(1.0-2.0*temp))/rr;
        break;
      case BUCKINGHAM_BOND:
        // p_0*exp(-p_1 r)-p_2/r^6
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A^-1]
        // p_2/k_B [K A^6]
        temp=parms[2]*CUBE(1.0/rr);
        exp_term=parms[0]*exp(-parms[1]*r);
        U=-temp+exp_term;
        DF=(6.0/rr)*temp-parms[1]*exp_term/r;
        break;
      case RESTRAINED_HARMONIC_BOND:
        // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
        // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        r1=r-parms[1];
        U=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
              +parms[0]*parms[2]*MAX2(fabs(r1)-parms[2],(REAL)0.0);
        DF=-parms[0]*(SIGN(MIN2(fabs(r1),parms[2]),r1))/r;
        break;
      case QUARTIC_BOND:
        // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
        // ===========================================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2/k_B [K/A^3]
        // p_3/k_B [K/A^4]
        temp=r-parms[1];
        temp2=SQR(r-parms[1]);
        U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
        DF=temp*(parms[0]+parms[2]*temp+parms[3]*temp2)/r;
        break;
      case CFF_QUARTIC_BOND:
        // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2/k_B [K/A^3]
        // p_3/k_B [K/A^4]
        temp=r-parms[1];
        temp2=SQR(r-parms[1]);
        U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
        DF=temp*(2.0*parms[0]+3.0*parms[2]*temp+4.0*parms[3]*temp2)/r;
        break;
      case MM3_BOND:
        // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
        // =================================================================
        // p_0     [mdyne/A molecule]
        // p_1     [A]
        temp=r-parms[1];
        temp2=SQR(temp);
        U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
        DF=parms[0]*(2.0+2.55*(4.0*2.55*(7.0/12.0)*temp-3.0)*temp)*temp/r;
        break;
      case RIGID_BOND:
        U=DF=0.0;
        break;
      case FIXED_BOND:
        U=DF=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Bond potential in routine 'CalculateAdsorbateBondForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // add contribution to the Adsorbate stretch energy
    UAdsorbateBond[CurrentSystem]+=U;

    // forces are oppositely directed to the gradient
    f.x=-DF*dr.x;
    f.y=-DF*dr.y;
    f.z=-DF*dr.z;

    // add contribution to the forces
    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=f.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=f.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=f.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x-=f.x;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y-=f.y;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z-=f.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax-=dr.x*f.x;
    StrainDerivativeTensor[CurrentSystem].bx-=dr.y*f.x;
    StrainDerivativeTensor[CurrentSystem].cx-=dr.z*f.x;

    StrainDerivativeTensor[CurrentSystem].ay-=dr.x*f.y;
    StrainDerivativeTensor[CurrentSystem].by-=dr.y*f.y;
    StrainDerivativeTensor[CurrentSystem].cy-=dr.z*f.y;

    StrainDerivativeTensor[CurrentSystem].az-=dr.x*f.z;
    StrainDerivativeTensor[CurrentSystem].bz-=dr.y*f.z;
    StrainDerivativeTensor[CurrentSystem].cz-=dr.z*f.z;
  }
}

// Calculate Bondforce
void CalculateCationBondForce(int m)
{
  int i,Type,NumberOfBonds,A,B;
  REAL U,DF,r,rr,temp,temp2,exp_term,r1;
  REAL *parms;
  POINT posA,posB;
  VECTOR dr,f;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfBonds=Components[Type].NumberOfBonds;
  for(i=0;i<NumberOfBonds;i++)
  {
    A=Components[Type].Bonds[i].A;
    B=Components[Type].Bonds[i].B;

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    parms=(REAL*)&Components[Type].BondArguments[i];

    switch(Components[Type].BondType[i])
    {
      case HARMONIC_BOND:
        // 0.5*p0*SQR(r-p1);
        // ===============================================
        // p_0/k_B [K/A^2]   force constant
        // p_1     [A]       reference bond distance
        U=0.5*parms[0]*SQR(r-parms[1]);
        DF=parms[0]*(r-parms[1])/r;
        break;
      case CORE_SHELL_SPRING:
        U=0.5*parms[0]*SQR(r);
        DF=parms[0];
        break;
      case MORSE_BOND:
        // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
        // ===============================================
        // p_0/k_B [K]       force constant
        // p_1     [A^-1]    parameter
        // p_2     [A]       reference bond distance
        temp=exp(parms[1]*(parms[2]-r));
        U=parms[0]*(SQR(1.0-temp)-1.0);
        DF=2.0*parms[0]*parms[1]*(1.0-temp)*temp/r;
        break;
      case LJ_12_6_BOND:
        // A/r_ij^12-B/r_ij^6
        // ===============================================
        // p_0/k_B [K A^12]
        // p_1/k_B [K A^6]
        temp=CUBE(1.0/rr);
        U=parms[0]*SQR(temp)-parms[1]*temp;
        DF=6.0*(parms[1]*temp-2.0*parms[0]*SQR(temp))/rr;
        break;
      case LENNARD_JONES_BOND:
        // 4*p_0*((p_1/r)^12-(p_1/r)^6)
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A]
        temp=CUBE(parms[1]/rr);
        U=4.0*parms[0]*(temp*(temp-1.0));
        DF=24.0*parms[0]*(temp*(1.0-2.0*temp))/rr;
        break;
      case BUCKINGHAM_BOND:
        // p_0*exp(-p_1 r)-p_2/r^6
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A^-1]
        // p_2/k_B [K A^6]
        temp=parms[2]*CUBE(1.0/rr);
        exp_term=parms[0]*exp(-parms[1]*r);
        U=-temp+exp_term;
        DF=(6.0/rr)*temp-parms[1]*exp_term/r;
        break;
      case RESTRAINED_HARMONIC_BOND:
        // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
        // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        r1=r-parms[1];
        U=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
              +parms[0]*parms[2]*MAX2(fabs(r1)-parms[2],(REAL)0.0);
        DF=-parms[0]*(SIGN(MIN2(fabs(r1),parms[2]),r1))/r;
        break;
      case QUARTIC_BOND:
        // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
        // ===========================================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2/k_B [K/A^3]
        // p_3/k_B [K/A^4]
        temp=r-parms[1];
        temp2=SQR(r-parms[1]);
        U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
        DF=temp*(parms[0]+parms[2]*temp+parms[3]*temp2)/r;
        break;
      case CFF_QUARTIC_BOND:
        // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2/k_B [K/A^3]
        // p_3/k_B [K/A^4]
        temp=r-parms[1];
        temp2=SQR(r-parms[1]);
        U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
        DF=temp*(2.0*parms[0]+3.0*parms[2]*temp+4.0*parms[3]*temp2)/r;
        break;
      case MM3_BOND:
        // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
        // =================================================================
        // p_0     [mdyne/A molecule]
        // p_1     [A]
        temp=r-parms[1];
        temp2=SQR(temp);
        U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
        DF=parms[0]*(2.0+2.55*(4.0*2.55*(7.0/12.0)*temp-3.0)*temp)*temp/r;
        break;
      case RIGID_BOND:
        U=DF=0.0;
        break;
      case FIXED_BOND:
        U=DF=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Bond potential in routine 'CalculateCationBondForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // add contribution to the Cation stretch energy
    UCationBond[CurrentSystem]+=U;

    // forces are oppositely directed to the gradient
    f.x=-DF*dr.x;
    f.y=-DF*dr.y;
    f.z=-DF*dr.z;

    // add contribution to the forces
    Cations[CurrentSystem][m].Atoms[A].Force.x+=f.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=f.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=f.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x-=f.x;
    Cations[CurrentSystem][m].Atoms[B].Force.y-=f.y;
    Cations[CurrentSystem][m].Atoms[B].Force.z-=f.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax-=dr.x*f.x;
    StrainDerivativeTensor[CurrentSystem].bx-=dr.y*f.x;
    StrainDerivativeTensor[CurrentSystem].cx-=dr.z*f.x;

    StrainDerivativeTensor[CurrentSystem].ay-=dr.x*f.y;
    StrainDerivativeTensor[CurrentSystem].by-=dr.y*f.y;
    StrainDerivativeTensor[CurrentSystem].cy-=dr.z*f.y;

    StrainDerivativeTensor[CurrentSystem].az-=dr.x*f.z;
    StrainDerivativeTensor[CurrentSystem].bz-=dr.y*f.z;
    StrainDerivativeTensor[CurrentSystem].cz-=dr.z*f.z;
  }
}

// Calculate Urey-Bradley force
void CalculateAdsorbateUreyBradleyForce(int m)
{
  int i,Type,NumberOfUreyBradleys,A,C;
  REAL U,DF,r,rr,temp,temp2,exp_term,r1;
  REAL *parms;
  POINT posA,posC;
  VECTOR dr,f;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfUreyBradleys=Components[Type].NumberOfUreyBradleys;
  for(i=0;i<NumberOfUreyBradleys;i++)
  {
    A=Components[Type].UreyBradleys[i].A;
    C=Components[Type].UreyBradleys[i].C;

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

    dr.x=posA.x-posC.x;
    dr.y=posA.y-posC.y;
    dr.z=posA.z-posC.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    parms=(REAL*)&Components[Type].UreyBradleyArguments[i];

    switch(Components[Type].UreyBradleyType[i])
    {
      case HARMONIC_UREYBRADLEY:
        // 0.5*p0*SQR(r-p1);
        // ===============================================
        // p_0/k_B [K/A^2]   force constant
        // p_1     [A]       reference bond distance
        U=0.5*parms[0]*SQR(r-parms[1]);
        DF=parms[0]*(r-parms[1])/r;
        break;
      case MORSE_UREYBRADLEY:
        // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
        // ===============================================
        // p_0/k_B [K]       force constant
        // p_1     [A^-1]    parameter
        // p_2     [A]       reference bond distance
        temp=exp(parms[1]*(parms[2]-r));
        U=parms[0]*(SQR(1.0-temp)-1.0);
        DF=2.0*parms[0]*parms[1]*(1.0-temp)*temp/r;
        break;
      case LJ_12_6_UREYBRADLEY:
        // A/r_ij^12-B/r_ij^6
        // ===============================================
        // p_0/k_B [K A^12]
        // p_1/k_B [K A^6]
        temp=CUBE(1.0/rr);
        U=parms[0]*SQR(temp)-parms[1]*temp;
        DF=6.0*(parms[1]*temp-2.0*parms[0]*SQR(temp))/rr;
        break;
      case LENNARD_JONES_UREYBRADLEY:
        // 4*p_0*((p_1/r)^12-(p_1/r)^6)
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A]
        temp=CUBE(parms[1]/rr);
        U=4.0*parms[0]*(temp*(temp-1.0));
        DF=24.0*parms[0]*(temp*(1.0-2.0*temp))/rr;
        break;
      case BUCKINGHAM_UREYBRADLEY:
        // p_0*exp(-p_1 r)-p_2/r^6
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A^-1]
        // p_2/k_B [K A^6]
        temp=parms[2]*CUBE(1.0/rr);
        exp_term=parms[0]*exp(-parms[1]*r);
        U=-temp+exp_term;
        DF=(6.0/rr)*temp-parms[1]*exp_term/r;
        break;
      case RESTRAINED_HARMONIC_UREYBRADLEY:
        // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
        // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        r1=r-parms[1];
        U=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
              +parms[0]*parms[2]*MAX2(fabs(r1)-parms[2],(REAL)0.0);
        DF=-parms[0]*(SIGN(MIN2(fabs(r1),parms[2]),r1))/r;
        break;
      case QUARTIC_UREYBRADLEY:
        // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
        // ===========================================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2/k_B [K/A^3]
        // p_3/k_B [K/A^4]
        temp=r-parms[1];
        temp2=SQR(r-parms[1]);
        U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
        DF=temp*(parms[0]+parms[2]*temp+parms[3]*temp2)/r;
        break;
      case CFF_QUARTIC_UREYBRADLEY:
        // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2/k_B [K/A^3]
        // p_3/k_B [K/A^4]
        temp=r-parms[1];
        temp2=SQR(r-parms[1]);
        U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
        DF=temp*(2.0*parms[0]+3.0*parms[2]*temp+4.0*parms[3]*temp2)/r;
        break;
      case MM3_UREYBRADLEY:
        // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
        // =================================================================
        // p_0     [mdyne/A molecule]
        // p_1     [A]
        temp=r-parms[1];
        temp2=SQR(temp);
        U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
        DF=parms[0]*(2.0+2.55*(4.0*2.55*(7.0/12.0)*temp-3.0)*temp)*temp/r;
        break;
      case RIGID_UREYBRADLEY:
        U=DF=0.0;
        break;
      case FIXED_UREYBRADLEY:
        U=DF=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Urey-Bradley potential in routine 'CalculateAdsorbateUreyBradleyForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // add contribution to the Adsorbate Urey-Bradley energy
    UAdsorbateUreyBradley[CurrentSystem]+=U;

    // forces are oppositely directed to the gradient
    f.x=-DF*dr.x;
    f.y=-DF*dr.y;
    f.z=-DF*dr.z;

    // add contribution to the forces
    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=f.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=f.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=f.z;

    Adsorbates[CurrentSystem][m].Atoms[C].Force.x-=f.x;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.y-=f.y;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.z-=f.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax-=dr.x*f.x;
    StrainDerivativeTensor[CurrentSystem].bx-=dr.y*f.x;
    StrainDerivativeTensor[CurrentSystem].cx-=dr.z*f.x;

    StrainDerivativeTensor[CurrentSystem].ay-=dr.x*f.y;
    StrainDerivativeTensor[CurrentSystem].by-=dr.y*f.y;
    StrainDerivativeTensor[CurrentSystem].cy-=dr.z*f.y;

    StrainDerivativeTensor[CurrentSystem].az-=dr.x*f.z;
    StrainDerivativeTensor[CurrentSystem].bz-=dr.y*f.z;
    StrainDerivativeTensor[CurrentSystem].cz-=dr.z*f.z;
  }
}

// Calculate Urey-Bradley force
void CalculateCationUreyBradleyForce(int m)
{
  int i,Type,NumberOfUreyBradleys,A,C;
  REAL U,DF,r,rr,temp,temp2,exp_term,r1;
  REAL *parms;
  POINT posA,posC;
  VECTOR dr,f;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfUreyBradleys=Components[Type].NumberOfUreyBradleys;
  for(i=0;i<NumberOfUreyBradleys;i++)
  {
    A=Components[Type].UreyBradleys[i].A;
    C=Components[Type].UreyBradleys[i].C;

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;

    dr.x=posA.x-posC.x;
    dr.y=posA.y-posC.y;
    dr.z=posA.z-posC.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    parms=(REAL*)&Components[Type].UreyBradleyArguments[i];

    switch(Components[Type].UreyBradleyType[i])
    {
      case HARMONIC_UREYBRADLEY:
        // 0.5*p0*SQR(r-p1);
        // ===============================================
        // p_0/k_B [K/A^2]   force constant
        // p_1     [A]       reference bond distance
        U=0.5*parms[0]*SQR(r-parms[1]);
        DF=parms[0]*(r-parms[1])/r;
        break;
      case MORSE_UREYBRADLEY:
        // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
        // ===============================================
        // p_0/k_B [K]       force constant
        // p_1     [A^-1]    parameter
        // p_2     [A]       reference bond distance
        temp=exp(parms[1]*(parms[2]-r));
        U=parms[0]*(SQR(1.0-temp)-1.0);
        DF=2.0*parms[0]*parms[1]*(1.0-temp)*temp/r;
        break;
      case LJ_12_6_UREYBRADLEY:
        // A/r_ij^12-B/r_ij^6
        // ===============================================
        // p_0/k_B [K A^12]
        // p_1/k_B [K A^6]
        temp=CUBE(1.0/rr);
        U=parms[0]*SQR(temp)-parms[1]*temp;
        DF=6.0*(parms[1]*temp-2.0*parms[0]*SQR(temp))/rr;
        break;
      case LENNARD_JONES_UREYBRADLEY:
        // 4*p_0*((p_1/r)^12-(p_1/r)^6)
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A]
        temp=CUBE(parms[1]/rr);
        U=4.0*parms[0]*(temp*(temp-1.0));
        DF=24.0*parms[0]*(temp*(1.0-2.0*temp))/rr;
        break;
      case BUCKINGHAM_UREYBRADLEY:
        // p_0*exp(-p_1 r)-p_2/r^6
        // ===============================================
        // p_0/k_B [K]
        // p_1     [A^-1]
        // p_2/k_B [K A^6]
        temp=parms[2]*CUBE(1.0/rr);
        exp_term=parms[0]*exp(-parms[1]*r);
        U=-temp+exp_term;
        DF=(6.0/rr)*temp-parms[1]*exp_term/r;
        break;
      case RESTRAINED_HARMONIC_UREYBRADLEY:
        // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
        // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        r1=r-parms[1];
        U=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
              +parms[0]*parms[2]*MAX2(fabs(r1)-parms[2],(REAL)0.0);
        DF=-parms[0]*(SIGN(MIN2(fabs(r1),parms[2]),r1))/r;
        break;
      case QUARTIC_UREYBRADLEY:
        // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
        // ===========================================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2/k_B [K/A^3]
        // p_3/k_B [K/A^4]
        temp=r-parms[1];
        temp2=SQR(r-parms[1]);
        U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
        DF=temp*(parms[0]+parms[2]*temp+parms[3]*temp2)/r;
        break;
      case CFF_QUARTIC_UREYBRADLEY:
        // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
        // ===============================================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2/k_B [K/A^3]
        // p_3/k_B [K/A^4]
        temp=r-parms[1];
        temp2=SQR(r-parms[1]);
        U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
        DF=temp*(2.0*parms[0]+3.0*parms[2]*temp+4.0*parms[3]*temp2)/r;
        break;
      case MM3_UREYBRADLEY:
        // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
        // =================================================================
        // p_0     [mdyne/A molecule]
        // p_1     [A]
        temp=r-parms[1];
        temp2=SQR(temp);
        U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
        DF=parms[0]*(2.0+2.55*(4.0*2.55*(7.0/12.0)*temp-3.0)*temp)*temp/r;
        break;
      case RIGID_UREYBRADLEY:
        U=DF=0.0;
        break;
      case FIXED_UREYBRADLEY:
        U=DF=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Urey-Bradley potential in routine 'CalculateCationUreyBradleyForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // add contribution to the Cation Urey-Bradley energy
    UCationUreyBradley[CurrentSystem]+=U;

    // forces are oppositely directed to the gradient
    f.x=-DF*dr.x;
    f.y=-DF*dr.y;
    f.z=-DF*dr.z;

    // add contribution to the forces
    Cations[CurrentSystem][m].Atoms[A].Force.x+=f.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=f.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=f.z;

    Cations[CurrentSystem][m].Atoms[C].Force.x-=f.x;
    Cations[CurrentSystem][m].Atoms[C].Force.y-=f.y;
    Cations[CurrentSystem][m].Atoms[C].Force.z-=f.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax-=dr.x*f.x;
    StrainDerivativeTensor[CurrentSystem].bx-=dr.y*f.x;
    StrainDerivativeTensor[CurrentSystem].cx-=dr.z*f.x;

    StrainDerivativeTensor[CurrentSystem].ay-=dr.x*f.y;
    StrainDerivativeTensor[CurrentSystem].by-=dr.y*f.y;
    StrainDerivativeTensor[CurrentSystem].cy-=dr.z*f.y;

    StrainDerivativeTensor[CurrentSystem].az-=dr.x*f.z;
    StrainDerivativeTensor[CurrentSystem].bz-=dr.y*f.z;
    StrainDerivativeTensor[CurrentSystem].cz-=dr.z*f.z;
  }
}


void CalculateAdsorbateBendForce(int m)
{
  int Type,NumberOfBends;
  int i,A,B,C,D;
  REAL *parms,DF,U,temp,temp2;
  REAL CosTheta,Theta,SinTheta;
  REAL rab,rbc,rac,DTDX;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac,fa,fb,fc,dtA,dtB,dtC;;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL rt2,rap2,rcp2;
  REAL terma,termc,delta,delta2,rm,ptrt2,term;
  VECTOR dpdia,dpdic,dedip,fd,n;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBends=Components[Type].NumberOfBends;
  for(i=0;i<NumberOfBends;i++)
  {
    A=Components[Type].Bends[i].A;
    B=Components[Type].Bends[i].B;
    C=Components[Type].Bends[i].C;
    D=Components[Type].Bends[i].D;
    parms=Components[Type].BendArguments[i];

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
    posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

    Rab.x=Rab.y=Rab.z=0.0;
    Rcd.x=Rcd.y=Rcd.z=0.0;
    Rbd.x=Rbd.y=Rbd.z=0.0;
    Rad.x=Rad.y=Rad.z=0.0;
    Rbc.x=Rbc.y=Rbc.z=0.0;
    ap.x=ap.y=ap.z=0.0;
    cp.x=cp.y=cp.z=0.0;
    t.x=t.y=t.z=0.0;
    delta=rcp2=rap2=rt2=rab=rbc=temp=0.0;
    switch(Components[Type].BendType[i])
    {
      case MM3_IN_PLANE_BEND:
        Rad.x=posA.x-posD.x;
        Rad.y=posA.y-posD.y;
        Rad.z=posA.z-posD.z;

        Rbd.x=posB.x-posD.x;
        Rbd.y=posB.y-posD.y;
        Rbd.z=posB.z-posD.z;

        Rcd.x=posC.x-posD.x;
        Rcd.y=posC.y-posD.y;
        Rcd.z=posC.z-posD.z;

        t.x=Rad.y*Rcd.z-Rad.z*Rcd.y;
        t.y=Rad.z*Rcd.x-Rad.x*Rcd.z;
        t.z=Rad.x*Rcd.y-Rad.y*Rcd.x;
        rt2=t.x*t.x+t.y*t.y+t.z*t.z;
        delta=-(t.x*Rbd.x+t.y*Rbd.y+t.z*Rbd.z)/rt2;

        ip.x=posB.x+t.x*delta;
        ip.y=posB.y+t.y*delta;
        ip.z=posB.z+t.z*delta;
        ap.x=posA.x-ip.x;
        ap.y=posA.y-ip.y;
        ap.z=posA.z-ip.z;
        cp.x=posC.x-ip.x;
        cp.y=posC.y-ip.y;
        cp.z=posC.z-ip.z;
        ap=ApplyBoundaryCondition(ap);
        cp=ApplyBoundaryCondition(cp);

        rap2=ap.x*ap.x+ap.y*ap.y+ap.z*ap.z;
        rcp2=cp.x*cp.x+cp.y*cp.y+cp.z*cp.z;

        CosTheta=(ap.x*cp.x+ap.y*cp.y+ap.z*cp.z)/sqrt(rap2*rcp2);
        break;
      default:
        Rab.x=posA.x-posB.x;
        Rab.y=posA.y-posB.y;
        Rab.z=posA.z-posB.z;
        rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
        Rab.x/=rab;
        Rab.y/=rab;
        Rab.z/=rab;

        Rbc.x=posC.x-posB.x;
        Rbc.y=posC.y-posB.y;
        Rbc.z=posC.z-posB.z;
        rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
        Rbc.x/=rbc;
        Rbc.y/=rbc;
        Rbc.z/=rbc;

        Rac.x=posC.x-posA.x;
        Rac.y=posC.y-posA.y;
        Rac.z=posC.z-posA.z;
        rac=sqrt(SQR(Rac.x)+SQR(Rac.y)+SQR(Rac.z));
        Rac.x/=rac;
        Rac.y/=rac;
        Rac.z/=rac;

        CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
        break;
    }

    CosTheta=MIN2(1.0,MAX2(-1.0,CosTheta));
    Theta=acos(CosTheta);
    SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));
    DTDX=-1.0/sqrt(1.0-SQR(CosTheta));

    switch(Components[Type].BendType[i])
    {
      case HARMONIC_BEND:
        // (1/2)p_0*(theta-p_1)^2
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        U=0.5*parms[0]*SQR(Theta-parms[1]);
        DF=parms[0]*(Theta-parms[1])*DTDX;
        break;
      case CORE_SHELL_BEND:
        // (1/2)p_0*(theta-p_1)^2
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        U=0.5*parms[0]*SQR(Theta-parms[1]);
        DF=parms[0]*(Theta-parms[1])*DTDX;
        break;
      case QUARTIC_BEND:
        // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
        // ======================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2/k_B [K/rad^3]
        // p_3/k_B [K/rad^4]
        temp=(Theta-parms[1]);
        temp2=SQR(temp);
        U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
        DF=(parms[0]*temp+parms[2]*temp2+parms[3]*temp*temp2)*DTDX;
        break;
      case CFF_QUARTIC_BEND:
        // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
        // =====================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2/k_B [K/rad^3]
        // p_3/k_B [K/rad^4]
        temp=(Theta-parms[1]);
        temp2=SQR(temp);
        U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
        DF=(2.0*parms[0]*temp+3.0*parms[2]*temp2+4.0*parms[3]*temp*temp2)*DTDX;
        break;
      case HARMONIC_COSINE_BEND:
        // (1/2)*p_0*(cos(theta)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        U=0.5*parms[0]*SQR(CosTheta-parms[1]);
        DF=parms[0]*(CosTheta-parms[1]);
        break;
      case COSINE_BEND:
        // p_0*(1+cos(p_1*theta-p_2))
        // ===============================================
        // p_0/k_B [K]
        // p_1     [-]
        // p_2     [degrees]
        temp=parms[1]*Theta-parms[2];
        U=parms[0]*(1.0+cos(temp));
        DF=-(parms[0]*parms[1]*sin(temp))*DTDX;
        break;
      case TAFIPOLSKY_BEND:
        // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
        // ===============================================
        // p_0/k_B [K]
        U=0.5*parms[0]*(1.0+cos(Theta))*(1.0+cos(2.0*Theta));
        DF=parms[0]*CosTheta*(2.0+3.0*CosTheta);
        break;
      case MM3_BEND:
      case MM3_IN_PLANE_BEND:
        // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
        // =================================================================================================
        // p_0/k_B [mdyne A/rad^2]
        // p_1     [degrees]
        U=DF=0;
        temp=RAD2DEG*(Theta-parms[1]);
        temp2=SQR(temp);
        U=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
        DF=parms[0]*RAD2DEG*(2.0-(3.0*0.014-(4.0*5.6e-5-(5.0*7.0e-7-6.0*2.2e-8*temp)*temp)*temp)*temp)*temp*DTDX;
        break;
      case FIXED_BEND:
        U=DF=0.0;
        break;
      case MEASURE_BEND:
        U=DF=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Bend potential in routine 'CalculateFrameworkBendForce' ('framework_force.c')\n");
        exit(0);
        break;
    }

    // add contribution to the energy
    UAdsorbateBend[CurrentSystem]+=U;

    switch(Components[Type].BendType[i])
    {
      case MM3_IN_PLANE_BEND:
        DF=parms[0]*RAD2DEG*(2.0-(3.0*0.014-(4.0*5.6e-5-(5.0*7.0e-7-6.0*2.2e-8*temp)*temp)*temp)*temp)*temp;
        n.x=cp.y*ap.z-cp.z*ap.y;
        n.y=cp.z*ap.x-cp.x*ap.z;
        n.z=cp.x*ap.y-cp.y*ap.x;
        rm=sqrt(n.x*n.x+n.y*n.y+n.z*n.z);
        rm=MAX2(rm,0.000001);

        terma=-DF/(rap2*rm);
        termc=DF/(rcp2*rm);
        fa.x=terma*(ap.y*n.z-ap.z*n.y);
        fa.y=terma*(ap.z*n.x-ap.x*n.z);
        fa.z=terma*(ap.x*n.y-ap.y*n.x);
        fc.x=termc*(cp.y*n.z-cp.z*n.y);
        fc.y=termc*(cp.z*n.x-cp.x*n.z);
        fc.z=termc*(cp.x*n.y-cp.y*n.x);
        dedip.x=-(fa.x+fc.x);
        dedip.y=-(fa.y+fc.y);
        dedip.z=-(fa.z+fc.z);

        delta2=2.0*delta;
        ptrt2=(dedip.x*t.x+dedip.y*t.y+dedip.z*t.z)/rt2;

        term=(Rcd.z*Rbd.y-Rcd.y*Rbd.z)+delta2*(t.y*Rcd.z-t.z*Rcd.y);
        dpdia.x=delta*(Rcd.y*dedip.z-Rcd.z*dedip.y)+term*ptrt2;

        term=(Rcd.x*Rbd.z-Rcd.z*Rbd.x)+delta2*(t.z*Rcd.x-t.x*Rcd.z);
        dpdia.y=delta*(Rcd.z*dedip.x-Rcd.x*dedip.z)+term*ptrt2;

        term=(Rcd.y*Rbd.x-Rcd.x*Rbd.y)+delta2*(t.x*Rcd.y-t.y*Rcd.x);
        dpdia.z=delta*(Rcd.x*dedip.y-Rcd.y*dedip.x)+term*ptrt2;

        term=(Rad.y*Rbd.z-Rad.z*Rbd.y)+delta2*(t.z*Rad.y-t.y*Rad.z);
        dpdic.x=delta*(Rad.z*dedip.y-Rad.y*dedip.z)+term*ptrt2;

        term=(Rad.z*Rbd.x-Rad.x*Rbd.z)+delta2*(t.x*Rad.z-t.z*Rad.x);
        dpdic.y=delta*(Rad.x*dedip.z-Rad.z*dedip.x)+term*ptrt2;

        term=(Rad.x*Rbd.y-Rad.y*Rbd.x)+delta2*(t.y*Rad.x-t.x*Rad.y);
        dpdic.z=delta*(Rad.y*dedip.x-Rad.x*dedip.y)+term*ptrt2;

        fa.x=fa.x+dpdia.x;
        fa.y=fa.y+dpdia.y;
        fa.z=fa.z+dpdia.z;
        fb.x=dedip.x;
        fb.y=dedip.y;
        fb.z=dedip.z;
        fc.x=fc.x+dpdic.x;
        fc.y=fc.y+dpdic.y;
        fc.z=fc.z+dpdic.z;
        fd.x=-(fa.x+fb.x+fc.x);
        fd.y=-(fa.y+fb.y+fc.y);
        fd.z=-(fa.z+fb.z+fc.z);

        // add contribution to the forces
        Adsorbates[CurrentSystem][m].Atoms[A].Force.x-=fa.x;
        Adsorbates[CurrentSystem][m].Atoms[A].Force.y-=fa.y;
        Adsorbates[CurrentSystem][m].Atoms[A].Force.z-=fa.z;

        Adsorbates[CurrentSystem][m].Atoms[B].Force.x-=fb.x;
        Adsorbates[CurrentSystem][m].Atoms[B].Force.y-=fb.y;
        Adsorbates[CurrentSystem][m].Atoms[B].Force.z-=fb.z;

        Adsorbates[CurrentSystem][m].Atoms[C].Force.x-=fc.x;
        Adsorbates[CurrentSystem][m].Atoms[C].Force.y-=fc.y;
        Adsorbates[CurrentSystem][m].Atoms[C].Force.z-=fc.z;

        Adsorbates[CurrentSystem][m].Atoms[D].Force.x-=fd.x;
        Adsorbates[CurrentSystem][m].Atoms[D].Force.y-=fd.y;
        Adsorbates[CurrentSystem][m].Atoms[D].Force.z-=fd.z;

        // add contribution to the stress tensor
        StrainDerivativeTensor[CurrentSystem].ax+=Rad.x*fa.x+Rbd.x*fb.x+Rcd.x*fc.x;
        StrainDerivativeTensor[CurrentSystem].bx+=Rad.y*fa.x+Rbd.y*fb.x+Rcd.y*fc.x;
        StrainDerivativeTensor[CurrentSystem].cx+=Rad.z*fa.x+Rbd.z*fb.x+Rcd.z*fc.x;

        StrainDerivativeTensor[CurrentSystem].ay+=Rad.x*fa.y+Rbd.x*fb.y+Rcd.x*fc.y;
        StrainDerivativeTensor[CurrentSystem].by+=Rad.y*fa.y+Rbd.y*fb.y+Rcd.y*fc.y;
        StrainDerivativeTensor[CurrentSystem].cy+=Rad.z*fa.y+Rbd.z*fb.y+Rcd.z*fc.y;

        StrainDerivativeTensor[CurrentSystem].az+=Rad.x*fa.z+Rbd.x*fb.z+Rcd.x*fc.z;
        StrainDerivativeTensor[CurrentSystem].bz+=Rad.y*fa.z+Rbd.y*fb.z+Rcd.y*fc.z;
        StrainDerivativeTensor[CurrentSystem].cz+=Rad.z*fa.z+Rbd.z*fb.z+Rcd.z*fc.z;
        break;
      default:
        // Calculate the components of the derivatives.
        dtA.x=(Rbc.x-CosTheta*Rab.x)/rab;
        dtA.y=(Rbc.y-CosTheta*Rab.y)/rab;
        dtA.z=(Rbc.z-CosTheta*Rab.z)/rab;

        dtC.x=(Rab.x-CosTheta*Rbc.x)/rbc;
        dtC.y=(Rab.y-CosTheta*Rbc.y)/rbc;
        dtC.z=(Rab.z-CosTheta*Rbc.z)/rbc;

        dtB.x=-(dtA.x+dtC.x);
        dtB.y=-(dtA.y+dtC.y);
        dtB.z=-(dtA.z+dtC.z);

        // forces are oppositely directed to the gradient
        fa.x=-DF*dtA.x;
        fa.y=-DF*dtA.y;
        fa.z=-DF*dtA.z;

        fb.x=-DF*dtB.x;
        fb.y=-DF*dtB.y;
        fb.z=-DF*dtB.z;

        fc.x=-DF*dtC.x;
        fc.y=-DF*dtC.y;
        fc.z=-DF*dtC.z;

          // add contribution to the forces
        Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
        Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
        Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

        Adsorbates[CurrentSystem][m].Atoms[B].Force.x+=fb.x;
        Adsorbates[CurrentSystem][m].Atoms[B].Force.y+=fb.y;
        Adsorbates[CurrentSystem][m].Atoms[B].Force.z+=fb.z;

        Adsorbates[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
        Adsorbates[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
        Adsorbates[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

        // add contribution to the stress tensor
        // Note: rab and rbc are here because the vectors were normalized before
        StrainDerivativeTensor[CurrentSystem].ax-=rab*Rab.x*fa.x+rbc*Rbc.x*fc.x;
        StrainDerivativeTensor[CurrentSystem].bx-=rab*Rab.y*fa.x+rbc*Rbc.y*fc.x;
        StrainDerivativeTensor[CurrentSystem].cx-=rab*Rab.z*fa.x+rbc*Rbc.z*fc.x;

        StrainDerivativeTensor[CurrentSystem].ay-=rab*Rab.x*fa.y+rbc*Rbc.x*fc.y;
        StrainDerivativeTensor[CurrentSystem].by-=rab*Rab.y*fa.y+rbc*Rbc.y*fc.y;
        StrainDerivativeTensor[CurrentSystem].cy-=rab*Rab.z*fa.y+rbc*Rbc.z*fc.y;

        StrainDerivativeTensor[CurrentSystem].az-=rab*Rab.x*fa.z+rbc*Rbc.x*fc.z;
        StrainDerivativeTensor[CurrentSystem].bz-=rab*Rab.y*fa.z+rbc*Rbc.y*fc.z;
        StrainDerivativeTensor[CurrentSystem].cz-=rab*Rab.z*fa.z+rbc*Rbc.z*fc.z;
        break;
    }
  }
}

void CalculateCationBendForce(int m)
{
  int Type,NumberOfBends;
  int i,A,B,C,D;
  REAL *parms,DF,U,temp,temp2;
  REAL CosTheta,Theta,SinTheta;
  REAL rab,rbc,rac,DTDX;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac,fa,fb,fc,dtA,dtB,dtC;;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL rt2,rap2,rcp2;
  REAL terma,termc,delta,delta2,rm,ptrt2,term;
  VECTOR dpdia,dpdic,dedip,fd,n;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfBends=Components[Type].NumberOfBends;
  for(i=0;i<NumberOfBends;i++)
  {
    A=Components[Type].Bends[i].A;
    B=Components[Type].Bends[i].B;
    C=Components[Type].Bends[i].C;
    D=Components[Type].Bends[i].D;
    parms=Components[Type].BendArguments[i];

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;
    posD=Cations[CurrentSystem][m].Atoms[D].Position;

    Rab.x=Rab.y=Rab.z=0.0;
    Rcd.x=Rcd.y=Rcd.z=0.0;
    Rbd.x=Rbd.y=Rbd.z=0.0;
    Rad.x=Rad.y=Rad.z=0.0;
    Rbc.x=Rbc.y=Rbc.z=0.0;
    ap.x=ap.y=ap.z=0.0;
    cp.x=cp.y=cp.z=0.0;
    t.x=t.y=t.z=0.0;
    delta=rcp2=rap2=rt2=rab=rbc=temp=0.0;
    switch(Components[Type].BendType[i])
    {
      case MM3_IN_PLANE_BEND:
        Rad.x=posA.x-posD.x;
        Rad.y=posA.y-posD.y;
        Rad.z=posA.z-posD.z;

        Rbd.x=posB.x-posD.x;
        Rbd.y=posB.y-posD.y;
        Rbd.z=posB.z-posD.z;

        Rcd.x=posC.x-posD.x;
        Rcd.y=posC.y-posD.y;
        Rcd.z=posC.z-posD.z;

        t.x=Rad.y*Rcd.z-Rad.z*Rcd.y;
        t.y=Rad.z*Rcd.x-Rad.x*Rcd.z;
        t.z=Rad.x*Rcd.y-Rad.y*Rcd.x;
        rt2=t.x*t.x+t.y*t.y+t.z*t.z;
        delta=-(t.x*Rbd.x+t.y*Rbd.y+t.z*Rbd.z)/rt2;

        ip.x=posB.x+t.x*delta;
        ip.y=posB.y+t.y*delta;
        ip.z=posB.z+t.z*delta;
        ap.x=posA.x-ip.x;
        ap.y=posA.y-ip.y;
        ap.z=posA.z-ip.z;
        cp.x=posC.x-ip.x;
        cp.y=posC.y-ip.y;
        cp.z=posC.z-ip.z;
        ap=ApplyBoundaryCondition(ap);
        cp=ApplyBoundaryCondition(cp);

        rap2=ap.x*ap.x+ap.y*ap.y+ap.z*ap.z;
        rcp2=cp.x*cp.x+cp.y*cp.y+cp.z*cp.z;

        CosTheta=(ap.x*cp.x+ap.y*cp.y+ap.z*cp.z)/sqrt(rap2*rcp2);
        break;
      default:
        Rab.x=posA.x-posB.x;
        Rab.y=posA.y-posB.y;
        Rab.z=posA.z-posB.z;
        rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
        Rab.x/=rab;
        Rab.y/=rab;
        Rab.z/=rab;

        Rbc.x=posC.x-posB.x;
        Rbc.y=posC.y-posB.y;
        Rbc.z=posC.z-posB.z;
        rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
        Rbc.x/=rbc;
        Rbc.y/=rbc;
        Rbc.z/=rbc;

        Rac.x=posC.x-posA.x;
        Rac.y=posC.y-posA.y;
        Rac.z=posC.z-posA.z;
        rac=sqrt(SQR(Rac.x)+SQR(Rac.y)+SQR(Rac.z));
        Rac.x/=rac;
        Rac.y/=rac;
        Rac.z/=rac;

        CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
        break;
    }

    CosTheta=MIN2(1.0,MAX2(-1.0,CosTheta));
    Theta=acos(CosTheta);
    SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));
    DTDX=-1.0/sqrt(1.0-SQR(CosTheta));

    switch(Components[Type].BendType[i])
    {
      case HARMONIC_BEND:
        // (1/2)p_0*(theta-p_1)^2
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        U=0.5*parms[0]*SQR(Theta-parms[1]);
        DF=parms[0]*(Theta-parms[1])*DTDX;
        break;
      case CORE_SHELL_BEND:
        // (1/2)p_0*(theta-p_1)^2
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        U=0.5*parms[0]*SQR(Theta-parms[1]);
        DF=parms[0]*(Theta-parms[1])*DTDX;
        break;
      case QUARTIC_BEND:
        // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
        // ======================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2/k_B [K/rad^3]
        // p_3/k_B [K/rad^4]
        temp=(Theta-parms[1]);
        temp2=SQR(temp);
        U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
        DF=(parms[0]*temp+parms[2]*temp2+parms[3]*temp*temp2)*DTDX;
        break;
      case CFF_QUARTIC_BEND:
        // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
        // =====================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2/k_B [K/rad^3]
        // p_3/k_B [K/rad^4]
        temp=(Theta-parms[1]);
        temp2=SQR(temp);
        U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
        DF=(2.0*parms[0]*temp+3.0*parms[2]*temp2+4.0*parms[3]*temp*temp2)*DTDX;
        break;
      case HARMONIC_COSINE_BEND:
        // (1/2)*p_0*(cos(theta)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        U=0.5*parms[0]*SQR(CosTheta-parms[1]);
        DF=parms[0]*(CosTheta-parms[1]);
        break;
      case COSINE_BEND:
        // p_0*(1+cos(p_1*theta-p_2))
        // ===============================================
        // p_0/k_B [K]
        // p_1     [-]
        // p_2     [degrees]
        temp=parms[1]*Theta-parms[2];
        U=parms[0]*(1.0+cos(temp));
        DF=-(parms[0]*parms[1]*sin(temp))*DTDX;
        break;
      case TAFIPOLSKY_BEND:
        // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
        // ===============================================
        // p_0/k_B [K]
        U=0.5*parms[0]*(1.0+cos(Theta))*(1.0+cos(2.0*Theta));
        DF=parms[0]*CosTheta*(2.0+3.0*CosTheta);
        break;
      case MM3_BEND:
      case MM3_IN_PLANE_BEND:
        // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
        // =================================================================================================
        // p_0/k_B [mdyne A/rad^2]
        // p_1     [degrees]
        U=DF=0;
        temp=RAD2DEG*(Theta-parms[1]);
        temp2=SQR(temp);
        U=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
        DF=parms[0]*RAD2DEG*(2.0-(3.0*0.014-(4.0*5.6e-5-(5.0*7.0e-7-6.0*2.2e-8*temp)*temp)*temp)*temp)*temp*DTDX;
        break;
      case FIXED_BEND:
        U=DF=0.0;
        break;
      case MEASURE_BEND:
        U=DF=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Bend potential in routine 'CalculateFrameworkBendForce' ('framework_force.c')\n");
        exit(0);
        break;
    }

    // add contribution to the energy
    UCationBend[CurrentSystem]+=U;

    switch(Components[Type].BendType[i])
    {
      case MM3_IN_PLANE_BEND:
        DF=parms[0]*RAD2DEG*(2.0-(3.0*0.014-(4.0*5.6e-5-(5.0*7.0e-7-6.0*2.2e-8*temp)*temp)*temp)*temp)*temp;
        n.x=cp.y*ap.z-cp.z*ap.y;
        n.y=cp.z*ap.x-cp.x*ap.z;
        n.z=cp.x*ap.y-cp.y*ap.x;
        rm=sqrt(n.x*n.x+n.y*n.y+n.z*n.z);
        rm=MAX2(rm,0.000001);

        terma=-DF/(rap2*rm);
        termc=DF/(rcp2*rm);
        fa.x=terma*(ap.y*n.z-ap.z*n.y);
        fa.y=terma*(ap.z*n.x-ap.x*n.z);
        fa.z=terma*(ap.x*n.y-ap.y*n.x);
        fc.x=termc*(cp.y*n.z-cp.z*n.y);
        fc.y=termc*(cp.z*n.x-cp.x*n.z);
        fc.z=termc*(cp.x*n.y-cp.y*n.x);
        dedip.x=-(fa.x+fc.x);
        dedip.y=-(fa.y+fc.y);
        dedip.z=-(fa.z+fc.z);

        delta2=2.0*delta;
        ptrt2=(dedip.x*t.x+dedip.y*t.y+dedip.z*t.z)/rt2;

        term=(Rcd.z*Rbd.y-Rcd.y*Rbd.z)+delta2*(t.y*Rcd.z-t.z*Rcd.y);
        dpdia.x=delta*(Rcd.y*dedip.z-Rcd.z*dedip.y)+term*ptrt2;

        term=(Rcd.x*Rbd.z-Rcd.z*Rbd.x)+delta2*(t.z*Rcd.x-t.x*Rcd.z);
        dpdia.y=delta*(Rcd.z*dedip.x-Rcd.x*dedip.z)+term*ptrt2;

        term=(Rcd.y*Rbd.x-Rcd.x*Rbd.y)+delta2*(t.x*Rcd.y-t.y*Rcd.x);
        dpdia.z=delta*(Rcd.x*dedip.y-Rcd.y*dedip.x)+term*ptrt2;

        term=(Rad.y*Rbd.z-Rad.z*Rbd.y)+delta2*(t.z*Rad.y-t.y*Rad.z);
        dpdic.x=delta*(Rad.z*dedip.y-Rad.y*dedip.z)+term*ptrt2;

        term=(Rad.z*Rbd.x-Rad.x*Rbd.z)+delta2*(t.x*Rad.z-t.z*Rad.x);
        dpdic.y=delta*(Rad.x*dedip.z-Rad.z*dedip.x)+term*ptrt2;

        term=(Rad.x*Rbd.y-Rad.y*Rbd.x)+delta2*(t.y*Rad.x-t.x*Rad.y);
        dpdic.z=delta*(Rad.y*dedip.x-Rad.x*dedip.y)+term*ptrt2;

        fa.x=fa.x+dpdia.x;
        fa.y=fa.y+dpdia.y;
        fa.z=fa.z+dpdia.z;
        fb.x=dedip.x;
        fb.y=dedip.y;
        fb.z=dedip.z;
        fc.x=fc.x+dpdic.x;
        fc.y=fc.y+dpdic.y;
        fc.z=fc.z+dpdic.z;
        fd.x=-(fa.x+fb.x+fc.x);
        fd.y=-(fa.y+fb.y+fc.y);
        fd.z=-(fa.z+fb.z+fc.z);

        // add contribution to the forces
        Cations[CurrentSystem][m].Atoms[A].Force.x-=fa.x;
        Cations[CurrentSystem][m].Atoms[A].Force.y-=fa.y;
        Cations[CurrentSystem][m].Atoms[A].Force.z-=fa.z;

        Cations[CurrentSystem][m].Atoms[B].Force.x-=fb.x;
        Cations[CurrentSystem][m].Atoms[B].Force.y-=fb.y;
        Cations[CurrentSystem][m].Atoms[B].Force.z-=fb.z;

        Cations[CurrentSystem][m].Atoms[C].Force.x-=fc.x;
        Cations[CurrentSystem][m].Atoms[C].Force.y-=fc.y;
        Cations[CurrentSystem][m].Atoms[C].Force.z-=fc.z;

        Cations[CurrentSystem][m].Atoms[D].Force.x-=fd.x;
        Cations[CurrentSystem][m].Atoms[D].Force.y-=fd.y;
        Cations[CurrentSystem][m].Atoms[D].Force.z-=fd.z;

        // add contribution to the stress tensor
        StrainDerivativeTensor[CurrentSystem].ax+=Rad.x*fa.x+Rbd.x*fb.x+Rcd.x*fc.x;
        StrainDerivativeTensor[CurrentSystem].bx+=Rad.y*fa.x+Rbd.y*fb.x+Rcd.y*fc.x;
        StrainDerivativeTensor[CurrentSystem].cx+=Rad.z*fa.x+Rbd.z*fb.x+Rcd.z*fc.x;

        StrainDerivativeTensor[CurrentSystem].ay+=Rad.x*fa.y+Rbd.x*fb.y+Rcd.x*fc.y;
        StrainDerivativeTensor[CurrentSystem].by+=Rad.y*fa.y+Rbd.y*fb.y+Rcd.y*fc.y;
        StrainDerivativeTensor[CurrentSystem].cy+=Rad.z*fa.y+Rbd.z*fb.y+Rcd.z*fc.y;

        StrainDerivativeTensor[CurrentSystem].az+=Rad.x*fa.z+Rbd.x*fb.z+Rcd.x*fc.z;
        StrainDerivativeTensor[CurrentSystem].bz+=Rad.y*fa.z+Rbd.y*fb.z+Rcd.y*fc.z;
        StrainDerivativeTensor[CurrentSystem].cz+=Rad.z*fa.z+Rbd.z*fb.z+Rcd.z*fc.z;
        break;
      default:
        // Calculate the components of the derivatives.
        dtA.x=(Rbc.x-CosTheta*Rab.x)/rab;
        dtA.y=(Rbc.y-CosTheta*Rab.y)/rab;
        dtA.z=(Rbc.z-CosTheta*Rab.z)/rab;

        dtC.x=(Rab.x-CosTheta*Rbc.x)/rbc;
        dtC.y=(Rab.y-CosTheta*Rbc.y)/rbc;
        dtC.z=(Rab.z-CosTheta*Rbc.z)/rbc;

        dtB.x=-(dtA.x+dtC.x);
        dtB.y=-(dtA.y+dtC.y);
        dtB.z=-(dtA.z+dtC.z);

        // forces are oppositely directed to the gradient
        fa.x=-DF*dtA.x;
        fa.y=-DF*dtA.y;
        fa.z=-DF*dtA.z;

        fb.x=-DF*dtB.x;
        fb.y=-DF*dtB.y;
        fb.z=-DF*dtB.z;

        fc.x=-DF*dtC.x;
        fc.y=-DF*dtC.y;
        fc.z=-DF*dtC.z;

          // add contribution to the forces
        Cations[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
        Cations[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
        Cations[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

        Cations[CurrentSystem][m].Atoms[B].Force.x+=fb.x;
        Cations[CurrentSystem][m].Atoms[B].Force.y+=fb.y;
        Cations[CurrentSystem][m].Atoms[B].Force.z+=fb.z;

        Cations[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
        Cations[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
        Cations[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

        // add contribution to the stress tensor
        // Note: rab and rbc are here because the vectors were normalized before
        StrainDerivativeTensor[CurrentSystem].ax-=rab*Rab.x*fa.x+rbc*Rbc.x*fc.x;
        StrainDerivativeTensor[CurrentSystem].bx-=rab*Rab.y*fa.x+rbc*Rbc.y*fc.x;
        StrainDerivativeTensor[CurrentSystem].cx-=rab*Rab.z*fa.x+rbc*Rbc.z*fc.x;

        StrainDerivativeTensor[CurrentSystem].ay-=rab*Rab.x*fa.y+rbc*Rbc.x*fc.y;
        StrainDerivativeTensor[CurrentSystem].by-=rab*Rab.y*fa.y+rbc*Rbc.y*fc.y;
        StrainDerivativeTensor[CurrentSystem].cy-=rab*Rab.z*fa.y+rbc*Rbc.z*fc.y;

        StrainDerivativeTensor[CurrentSystem].az-=rab*Rab.x*fa.z+rbc*Rbc.x*fc.z;
        StrainDerivativeTensor[CurrentSystem].bz-=rab*Rab.y*fa.z+rbc*Rbc.y*fc.z;
        StrainDerivativeTensor[CurrentSystem].cz-=rab*Rab.z*fa.z+rbc*Rbc.z*fc.z;
        break;
    }
  }
}


// derivative results verified by mathematica 4.1
void CalculateAdsorbateInversionBendForce(int m)
{
  int i,A,B,C,D;
  int Type,NumberOfInversionBends;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2,rrbc,rbc2,rbd2,rad2,rac2,dot;
  REAL CosChi,Chi,energy;
  REAL temp,temp2;
  VECTOR Rab,Rbc,Rbd,Rad,Rac;
  POINT posA,posB,posC,posD;
  VECTOR fa,fb,fc,fd;
  REAL term,dedcos;
  VECTOR dccdia,dccdic,dccdid;
  VECTOR deedia,deedic,deedid;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfInversionBends=Components[Type].NumberOfInversionBends;
  for(i=0;i<NumberOfInversionBends;i++)
  {
    A=Components[Type].InversionBends[i].A;
    B=Components[Type].InversionBends[i].B;
    C=Components[Type].InversionBends[i].C;
    D=Components[Type].InversionBends[i].D;
    parms=Components[Type].InversionBendArguments[i];

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
    posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

    Rab.x=posA.x-posB.x;
    Rab.y=posA.y-posB.y;
    Rab.z=posA.z-posB.z;
    rab2=Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z;
    rrab=sqrt(rab2);

    Rbc.x=posC.x-posB.x;
    Rbc.y=posC.y-posB.y;
    Rbc.z=posC.z-posB.z;
    rbc2=Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z;
    rrbc=sqrt(rbc2);

    Rbd.x=posD.x-posB.x;
    Rbd.y=posD.y-posB.y;
    Rbd.z=posD.z-posB.z;
    rbd2=Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z;

    Rac.x=posC.x-posA.x;
    Rac.y=posC.y-posA.y;
    Rac.z=posC.z-posA.z;
    rac2=Rac.x*Rac.x+Rac.y*Rac.y+Rac.z*Rac.z;

    Rad.x=posD.x-posA.x;
    Rad.y=posD.y-posA.y;
    Rad.z=posD.z-posA.z;
    rad2=Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z;

    switch(Components[Type].InversionBendType[i])
    {
      case HARMONIC_INVERSION:
      case HARMONIC_COSINE_INVERSION:
      case PLANAR_INVERSION:
        // w is a vector perpendicular to the B-C-D plane
        // c=w.w=(Rbc x Rbd).(Rbc x Rbd)= r_bc^2 r_bd^2 - (r_cb . r_bd)^2
        dot=Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z;
        c=rbc2*rbd2-SQR(dot);
        //c=(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z)*(Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z)-SQR(Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z);
        break;
      case MM3_INVERSION:
      case HARMONIC_INVERSION2:
      case HARMONIC_COSINE_INVERSION2:
      case PLANAR_INVERSION2:
        // w is a vector perpendicular to the A-C-D plane
        // c=w.w=(Rcd x Rad).(Rcd x Rad)=r_cd^2 r_ad^2 - (r_da . r_cd)^2
        dot=Rac.x*Rad.x+Rac.y*Rad.y+Rac.z*Rad.z;
        c=rac2*rad2-SQR(dot);
        //c=(Rcd.x*Rcd.x+Rcd.y*Rcd.y+Rcd.z*Rcd.z)*(Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z)-SQR(Rad.x*Rcd.x+Rad.y*Rcd.y+Rad.z*Rcd.z);
        break;
      default:
        fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
        exit(0);
        break;
    }

    e=Rab.x*(Rbc.y*Rbd.z-Rbc.z*Rbd.y)+Rab.y*(Rbc.z*Rbd.x-Rbc.x*Rbd.z)+Rab.z*(Rbc.x*Rbd.y-Rbc.y*Rbd.x);
    CosChi=sqrt((Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z)-SQR(e)/c)/rrab;

    // Ensure CosChi is between -1 and 1.
    CosChi=SIGN(MIN2(fabs(CosChi),(REAL)1.0),CosChi);

    switch(Components[Type].InversionBendType[i])
    {
      case HARMONIC_INVERSION:
      case HARMONIC_INVERSION2:
        // (1/2)*p_0*(chi-p_1)^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        Chi=acos(CosChi);
        energy=0.5*parms[0]*SQR(Chi-parms[1]);
        dedcos=-SIGN(1.0,e)*(parms[0]*(Chi-parms[1])/sqrt(c*(rab2-e*e/c)));
        break;
      case HARMONIC_COSINE_INVERSION:
      case HARMONIC_COSINE_INVERSION2:
        // (1/2)*p_0*(cos(phi)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        Chi=acos(CosChi);
        energy=0.5*parms[0]*SQR(CosChi-parms[1]);
        dedcos=SIGN(1.0,e)*parms[0]*(CosChi-parms[1])*sin(Chi)/sqrt(c*(rab2-e*e/c));
        break;
      case PLANAR_INVERSION:
      case PLANAR_INVERSION2:
        // (1/2)*p_0*(1-cos(phi))
        // ===============================================
        // p_0/k_B [K]
        Chi=acos(CosChi);
        energy=parms[0]*(1.0-CosChi);
        dedcos=-SIGN(1.0,e)*parms[0]*sin(Chi)/sqrt(c*(rab2-e*e/c));
        break;
      case MM3_INVERSION:
        // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
        // =================================================================================================
        // p_0/k_B [mdyne A/rad^2]
        // p_1     [degrees]
        Chi=acos(CosChi);
        temp=RAD2DEG*(Chi-parms[1]);
        temp2=SQR(temp);
        energy=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
        dedcos=-SIGN(1.0,e)*parms[0]*temp*RAD2DEG*(2.0-3.0*0.014*temp+4.0*5.6e-5*temp2-5.0*7.0e-7*temp*temp2+6.0*2.2e-8*SQR(temp2))/sqrt(c*(rab2-e*e/c));
        break;
      default:
        fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
        exit(0);
        break;
    }

    // energy
    UAdsorbateInversionBend[CurrentSystem]+=energy;

    switch(Components[Type].InversionBendType[i])
    {
      case HARMONIC_COSINE_INVERSION:
      case PLANAR_INVERSION:
      case HARMONIC_INVERSION:
        term=e/c;
        dccdia.x=0.0;
        dccdia.y=0.0;
        dccdia.z=0.0;
        dccdic.x=(Rbc.x*rbd2-Rbd.x*dot)*term;
        dccdic.y=(Rbc.y*rbd2-Rbd.y*dot)*term;
        dccdic.z=(Rbc.z*rbd2-Rbd.z*dot)*term;
        dccdid.x=(Rbd.x*rbc2-Rbc.x*dot)*term;
        dccdid.y=(Rbd.y*rbc2-Rbc.y*dot)*term;
        dccdid.z=(Rbd.z*rbc2-Rbc.z*dot)*term;
        break;
      case HARMONIC_INVERSION2:
      case HARMONIC_COSINE_INVERSION2:
      case PLANAR_INVERSION2:
      case MM3_INVERSION:
        term=e/c;
        dccdic.x=(Rac.x*rad2-Rad.x*dot)*term;
        dccdic.y=(Rac.y*rad2-Rad.y*dot)*term;
        dccdic.z=(Rac.z*rad2-Rad.z*dot)*term;
        dccdid.x=(Rad.x*rac2-Rac.x*dot)*term;
        dccdid.y=(Rad.y*rac2-Rac.y*dot)*term;
        dccdid.z=(Rad.z*rac2-Rac.z*dot)*term;
        dccdia.x=-(dccdic.x+dccdid.x);
        dccdia.y=-(dccdic.y+dccdid.y);
        dccdia.z=-(dccdic.z+dccdid.z);
        break;
      default:
        fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
        exit(0);
        break;
    }


    term=e/rab2;
    deedia.x=Rbd.y*Rbc.z-Rbd.z*Rbc.y+Rab.x*term;
    deedia.y=Rbd.z*Rbc.x-Rbd.x*Rbc.z+Rab.y*term;
    deedia.z=Rbd.x*Rbc.y-Rbd.y*Rbc.x+Rab.z*term;
    deedic.x=Rab.y*Rbd.z-Rab.z*Rbd.y;
    deedic.y=Rab.z*Rbd.x-Rab.x*Rbd.z;
    deedic.z=Rab.x*Rbd.y-Rab.y*Rbd.x;
    deedid.x=Rbc.y*Rab.z-Rbc.z*Rab.y;
    deedid.y=Rbc.z*Rab.x-Rbc.x*Rab.z;
    deedid.z=Rbc.x*Rab.y-Rbc.y*Rab.x;

    fa.x=dedcos*(dccdia.x+deedia.x);
    fa.y=dedcos*(dccdia.y+deedia.y);
    fa.z=dedcos*(dccdia.z+deedia.z);
    fc.x=dedcos*(dccdic.x+deedic.x);
    fc.y=dedcos*(dccdic.y+deedic.y);
    fc.z=dedcos*(dccdic.z+deedic.z);
    fd.x=dedcos*(dccdid.x+deedid.x);
    fd.y=dedcos*(dccdid.y+deedid.y);
    fd.z=dedcos*(dccdid.z+deedid.z);

    fb.x=-(fa.x+fc.x+fd.x);
    fb.y=-(fa.y+fc.y+fd.y);
    fb.z=-(fa.z+fc.z+fd.z);

    // add contribution to the forces
    Adsorbates[CurrentSystem][m].Atoms[A].Force.x-=fa.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y-=fa.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z-=fa.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x-=fb.x;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y-=fb.y;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z-=fb.z;

    Adsorbates[CurrentSystem][m].Atoms[C].Force.x-=fc.x;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.y-=fc.y;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.z-=fc.z;

    Adsorbates[CurrentSystem][m].Atoms[D].Force.x-=fd.x;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.y-=fd.y;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.z-=fd.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax+=Rab.x*fa.x+Rbc.x*fc.x+Rbd.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].ay+=Rab.x*fa.y+Rbc.x*fc.y+Rbd.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].az+=Rab.x*fa.z+Rbc.x*fc.z+Rbd.x*fd.z;

    StrainDerivativeTensor[CurrentSystem].bx+=Rab.y*fa.x+Rbc.y*fc.x+Rbd.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].by+=Rab.y*fa.y+Rbc.y*fc.y+Rbd.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].bz+=Rab.y*fa.z+Rbc.y*fc.z+Rbd.y*fd.z;

    StrainDerivativeTensor[CurrentSystem].cx+=Rab.z*fa.x+Rbc.z*fc.x+Rbd.z*fd.x;
    StrainDerivativeTensor[CurrentSystem].cy+=Rab.z*fa.y+Rbc.z*fc.y+Rbd.z*fd.y;
    StrainDerivativeTensor[CurrentSystem].cz+=Rab.z*fa.z+Rbc.z*fc.z+Rbd.z*fd.z;
  }
}


// derivative results verified by mathematica 4.1
void CalculateCationInversionBendForce(int m)
{
  int i,A,B,C,D;
  int Type,NumberOfInversionBends;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2,rrbc,rbc2,rbd2,rad2,rac2,dot;
  REAL CosChi,Chi,energy;
  REAL temp,temp2;
  VECTOR Rab,Rbc,Rbd,Rad,Rac;
  POINT posA,posB,posC,posD;
  VECTOR fa,fb,fc,fd;
  REAL term,dedcos;
  VECTOR dccdia,dccdic,dccdid;
  VECTOR deedia,deedic,deedid;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfInversionBends=Components[Type].NumberOfInversionBends;
  for(i=0;i<NumberOfInversionBends;i++)
  {
    A=Components[Type].InversionBends[i].A;
    B=Components[Type].InversionBends[i].B;
    C=Components[Type].InversionBends[i].C;
    D=Components[Type].InversionBends[i].D;
    parms=Components[Type].InversionBendArguments[i];

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;
    posD=Cations[CurrentSystem][m].Atoms[D].Position;

    Rab.x=posA.x-posB.x;
    Rab.y=posA.y-posB.y;
    Rab.z=posA.z-posB.z;
    rab2=Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z;
    rrab=sqrt(rab2);

    Rbc.x=posC.x-posB.x;
    Rbc.y=posC.y-posB.y;
    Rbc.z=posC.z-posB.z;
    rbc2=Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z;
    rrbc=sqrt(rbc2);

    Rbd.x=posD.x-posB.x;
    Rbd.y=posD.y-posB.y;
    Rbd.z=posD.z-posB.z;
    rbd2=Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z;

    Rac.x=posC.x-posA.x;
    Rac.y=posC.y-posA.y;
    Rac.z=posC.z-posA.z;
    rac2=Rac.x*Rac.x+Rac.y*Rac.y+Rac.z*Rac.z;

    Rad.x=posD.x-posA.x;
    Rad.y=posD.y-posA.y;
    Rad.z=posD.z-posA.z;
    rad2=Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z;

    switch(Components[Type].InversionBendType[i])
    {
      case HARMONIC_INVERSION:
      case HARMONIC_COSINE_INVERSION:
      case PLANAR_INVERSION:
        // w is a vector perpendicular to the B-C-D plane
        // c=w.w=(Rbc x Rbd).(Rbc x Rbd)= r_bc^2 r_bd^2 - (r_cb . r_bd)^2
        dot=Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z;
        c=rbc2*rbd2-SQR(dot);
        //c=(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z)*(Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z)-SQR(Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z);
        break;
      case MM3_INVERSION:
      case HARMONIC_INVERSION2:
      case HARMONIC_COSINE_INVERSION2:
      case PLANAR_INVERSION2:
        // w is a vector perpendicular to the A-C-D plane
        // c=w.w=(Rcd x Rad).(Rcd x Rad)=r_cd^2 r_ad^2 - (r_da . r_cd)^2
        dot=Rac.x*Rad.x+Rac.y*Rad.y+Rac.z*Rad.z;
        c=rac2*rad2-SQR(dot);
        //c=(Rcd.x*Rcd.x+Rcd.y*Rcd.y+Rcd.z*Rcd.z)*(Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z)-SQR(Rad.x*Rcd.x+Rad.y*Rcd.y+Rad.z*Rcd.z);
        break;
      default:
        fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
        exit(0);
        break;
    }

    e=Rab.x*(Rbc.y*Rbd.z-Rbc.z*Rbd.y)+Rab.y*(Rbc.z*Rbd.x-Rbc.x*Rbd.z)+Rab.z*(Rbc.x*Rbd.y-Rbc.y*Rbd.x);
    CosChi=sqrt((Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z)-SQR(e)/c)/rrab;

    // Ensure CosChi is between -1 and 1.
    CosChi=SIGN(MIN2(fabs(CosChi),(REAL)1.0),CosChi);

    switch(Components[Type].InversionBendType[i])
    {
      case HARMONIC_INVERSION:
      case HARMONIC_INVERSION2:
        // (1/2)*p_0*(chi-p_1)^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        Chi=acos(CosChi);
        energy=0.5*parms[0]*SQR(Chi-parms[1]);
        dedcos=-SIGN(1.0,e)*(parms[0]*(Chi-parms[1])/sqrt(c*(rab2-e*e/c)));
        break;
      case HARMONIC_COSINE_INVERSION:
      case HARMONIC_COSINE_INVERSION2:
        // (1/2)*p_0*(cos(phi)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        Chi=acos(CosChi);
        energy=0.5*parms[0]*SQR(CosChi-parms[1]);
        dedcos=SIGN(1.0,e)*parms[0]*(CosChi-parms[1])*sin(Chi)/sqrt(c*(rab2-e*e/c));
        break;
      case PLANAR_INVERSION:
      case PLANAR_INVERSION2:
        // (1/2)*p_0*(1-cos(phi))
        // ===============================================
        // p_0/k_B [K]
        Chi=acos(CosChi);
        energy=parms[0]*(1.0-CosChi);
        dedcos=-SIGN(1.0,e)*parms[0]*sin(Chi)/sqrt(c*(rab2-e*e/c));
        break;
      case MM3_INVERSION:
        // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
        // =================================================================================================
        // p_0/k_B [mdyne A/rad^2]
        // p_1     [degrees]
        Chi=acos(CosChi);
        temp=RAD2DEG*(Chi-parms[1]);
        temp2=SQR(temp);
        energy=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
        dedcos=-SIGN(1.0,e)*parms[0]*temp*RAD2DEG*(2.0-3.0*0.014*temp+4.0*5.6e-5*temp2-5.0*7.0e-7*temp*temp2+6.0*2.2e-8*SQR(temp2))/sqrt(c*(rab2-e*e/c));
        break;
      default:
        fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
        exit(0);
        break;
    }

    // energy
    UCationInversionBend[CurrentSystem]+=energy;

    switch(Components[Type].InversionBendType[i])
    {
      case HARMONIC_COSINE_INVERSION:
      case PLANAR_INVERSION:
      case HARMONIC_INVERSION:
        term=e/c;
        dccdia.x=0.0;
        dccdia.y=0.0;
        dccdia.z=0.0;
        dccdic.x=(Rbc.x*rbd2-Rbd.x*dot)*term;
        dccdic.y=(Rbc.y*rbd2-Rbd.y*dot)*term;
        dccdic.z=(Rbc.z*rbd2-Rbd.z*dot)*term;
        dccdid.x=(Rbd.x*rbc2-Rbc.x*dot)*term;
        dccdid.y=(Rbd.y*rbc2-Rbc.y*dot)*term;
        dccdid.z=(Rbd.z*rbc2-Rbc.z*dot)*term;
        break;
      case HARMONIC_INVERSION2:
      case HARMONIC_COSINE_INVERSION2:
      case PLANAR_INVERSION2:
      case MM3_INVERSION:
        term=e/c;
        dccdic.x=(Rac.x*rad2-Rad.x*dot)*term;
        dccdic.y=(Rac.y*rad2-Rad.y*dot)*term;
        dccdic.z=(Rac.z*rad2-Rad.z*dot)*term;
        dccdid.x=(Rad.x*rac2-Rac.x*dot)*term;
        dccdid.y=(Rad.y*rac2-Rac.y*dot)*term;
        dccdid.z=(Rad.z*rac2-Rac.z*dot)*term;
        dccdia.x=-(dccdic.x+dccdid.x);
        dccdia.y=-(dccdic.y+dccdid.y);
        dccdia.z=-(dccdic.z+dccdid.z);
        break;
      default:
        fprintf(stderr, "Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
        exit(0);
        break;
    }


    term=e/rab2;
    deedia.x=Rbd.y*Rbc.z-Rbd.z*Rbc.y+Rab.x*term;
    deedia.y=Rbd.z*Rbc.x-Rbd.x*Rbc.z+Rab.y*term;
    deedia.z=Rbd.x*Rbc.y-Rbd.y*Rbc.x+Rab.z*term;
    deedic.x=Rab.y*Rbd.z-Rab.z*Rbd.y;
    deedic.y=Rab.z*Rbd.x-Rab.x*Rbd.z;
    deedic.z=Rab.x*Rbd.y-Rab.y*Rbd.x;
    deedid.x=Rbc.y*Rab.z-Rbc.z*Rab.y;
    deedid.y=Rbc.z*Rab.x-Rbc.x*Rab.z;
    deedid.z=Rbc.x*Rab.y-Rbc.y*Rab.x;

    fa.x=dedcos*(dccdia.x+deedia.x);
    fa.y=dedcos*(dccdia.y+deedia.y);
    fa.z=dedcos*(dccdia.z+deedia.z);
    fc.x=dedcos*(dccdic.x+deedic.x);
    fc.y=dedcos*(dccdic.y+deedic.y);
    fc.z=dedcos*(dccdic.z+deedic.z);
    fd.x=dedcos*(dccdid.x+deedid.x);
    fd.y=dedcos*(dccdid.y+deedid.y);
    fd.z=dedcos*(dccdid.z+deedid.z);

    fb.x=-(fa.x+fc.x+fd.x);
    fb.y=-(fa.y+fc.y+fd.y);
    fb.z=-(fa.z+fc.z+fd.z);

    // add contribution to the forces
    Cations[CurrentSystem][m].Atoms[A].Force.x-=fa.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y-=fa.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z-=fa.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x-=fb.x;
    Cations[CurrentSystem][m].Atoms[B].Force.y-=fb.y;
    Cations[CurrentSystem][m].Atoms[B].Force.z-=fb.z;

    Cations[CurrentSystem][m].Atoms[C].Force.x-=fc.x;
    Cations[CurrentSystem][m].Atoms[C].Force.y-=fc.y;
    Cations[CurrentSystem][m].Atoms[C].Force.z-=fc.z;

    Cations[CurrentSystem][m].Atoms[D].Force.x-=fd.x;
    Cations[CurrentSystem][m].Atoms[D].Force.y-=fd.y;
    Cations[CurrentSystem][m].Atoms[D].Force.z-=fd.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax+=Rab.x*fa.x+Rbc.x*fc.x+Rbd.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].ay+=Rab.x*fa.y+Rbc.x*fc.y+Rbd.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].az+=Rab.x*fa.z+Rbc.x*fc.z+Rbd.x*fd.z;

    StrainDerivativeTensor[CurrentSystem].bx+=Rab.y*fa.x+Rbc.y*fc.x+Rbd.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].by+=Rab.y*fa.y+Rbc.y*fc.y+Rbd.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].bz+=Rab.y*fa.z+Rbc.y*fc.z+Rbd.y*fd.z;

    StrainDerivativeTensor[CurrentSystem].cx+=Rab.z*fa.x+Rbc.z*fc.x+Rbd.z*fd.x;
    StrainDerivativeTensor[CurrentSystem].cy+=Rab.z*fa.y+Rbc.z*fc.y+Rbd.z*fd.y;
    StrainDerivativeTensor[CurrentSystem].cz+=Rab.z*fa.z+Rbc.z*fc.z+Rbd.z*fd.z;
  }
}


void CalculateAdsorbateTorsionForce(int m)
{
  int i,Type,NumberOfTorsions,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi,U,DF;
  REAL ShiftedCosPhi,ShiftedCosPhi2,ShiftedSinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  VECTOR fa,fb,fc,fd;
  REAL *parms;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfTorsions=Components[Type].NumberOfTorsions;
  for(i=0;i<NumberOfTorsions;i++)
  {
    A=Components[Type].Torsions[i].A;
    B=Components[Type].Torsions[i].B;
    C=Components[Type].Torsions[i].C;
    D=Components[Type].Torsions[i].D;
    parms=(REAL*)&Components[Type].TorsionArguments[i];

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
    posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;

    Dcb.x=posC.x-posB.x;
    Dcb.y=posC.y-posB.y;
    Dcb.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
    Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

    Ddc.x=posD.x-posC.x;
    Ddc.y=posD.y-posC.y;
    Ddc.z=posD.z-posC.z;
    dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
    dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

    dr.x=Dab.x-dot_ab*Dcb.x;
    dr.y=Dab.y-dot_ab*Dcb.y;
    dr.z=Dab.z-dot_ab*Dcb.z;
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    dr.x/=r; dr.y/=r; dr.z/=r;

    ds.x=Ddc.x-dot_cd*Dcb.x;
    ds.y=Ddc.y-dot_cd*Dcb.y;
    ds.z=Ddc.z-dot_cd*Dcb.z;
    s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
    ds.x/=s; ds.y/=s; ds.z/=s;

    // compute Cos(Phi)
    // Phi is defined in protein convention Phi(trans)=Pi
    CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].TorsionType[i])
    {
      case HARMONIC_DIHEDRAL:
        // (1/2)*p_0*(phi-p_1)^2
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
        Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
        Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
        Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
        Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
        Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
        sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        Phi-=parms[1];
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
        U=0.5*parms[0]*SQR(Phi);
        DF=-parms[0]*(Phi)/SinPhi;
        break;
      case HARMONIC_COSINE_DIHEDRAL:
        // (1/2)*p_0*(cos(phi)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        U=0.5*parms[0]*SQR(CosPhi-parms[1]);
        DF=parms[0]*(CosPhi-parms[1]);
        break;
      case THREE_COSINE_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
        break;
      case MM3_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0     [kcal/mol]
        // p_1     [kcal/mol]
        // p_2     [kcal/mol]
        U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
        break;
      case CVFF_BLOCKED_DIHEDRAL:
        //
        // ========================================================================
        // p_0     [rad]
        // p_1     [K]
        // p_2     [-]
        // p_3     [rad]
        // p_4     [rad]
        U=0.0;
        DF=0.0;
        break;
      case CFF_DIHEDRAL:
        // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
        DF=-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2);
        break;
      case CFF_DIHEDRAL2:
        // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
        DF=parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi);
        break;
      case SIX_COSINE_DIHEDRAL:
        // Prod_i=0^5 p_i*cos(phi)^i
        // =========================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
        // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
        // polymer convention Phi'(trans)=0.
        U=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
               parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
        DF=-parms[1]+2.0*parms[2]*CosPhi-3.0*parms[3]*CosPhi2+4.0*parms[4]*CosPhi2*CosPhi-5.0*parms[5]*SQR(CosPhi2);
        break;
      case TRAPPE_DIHEDRAL:
        // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
        // ==========================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
        DF=parms[1]-4.0*parms[2]*CosPhi+3.0*parms[3]*(4.0*CosPhi2-1.0);
        break;
      case TRAPPE_DIHEDRAL_EXTENDED:
        // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
        // =============================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        U=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
        DF=parms[1]-3.0*parms[3]+4.0*(parms[2]-4.0*parms[4])*CosPhi+12.0*parms[3]*CosPhi2+32.0*parms[4]*CosPhi2*CosPhi;
        break;
      case MOD_TRAPPE_DIHEDRAL:
        /* Salvador modification: 16/08/2016
         add phase in cos function:
         p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
        */
        Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
        Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
        Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
        Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
        Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
        Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
        sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
        Phi-=parms[4];           // shift Phi as Phi+parms[4]
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        ShiftedCosPhi=cos(Phi);
        ShiftedSinPhi=sin(Phi);
        ShiftedCosPhi2=SQR(ShiftedCosPhi);
        U=parms[0]+parms[1]+parms[3]+(parms[1]-3.0*parms[3])*ShiftedCosPhi -2.0*parms[2]*ShiftedCosPhi2 +4.0*parms[3]*ShiftedCosPhi*ShiftedCosPhi2;
        DF=((parms[1]-3.0*parms[3])*sin(Phi) -4.0*parms[2]*ShiftedCosPhi*ShiftedSinPhi + 12.0*parms[3]*ShiftedCosPhi2*ShiftedSinPhi)/SinPhi;
        break;
      case CVFF_DIHEDRAL:
        // p_0*(1+cos(p_1*phi-p_2))
        // ========================
        // p_0/k_B [K]
        // p_1     [-]
        // p_2     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
        Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
        Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
        Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
        Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
        Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
        sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
        U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
        DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
        break;
      case OPLS_DIHEDRAL:
        // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
        // =================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
        DF=0.5*parms[1]-2.0*parms[2]*CosPhi+1.5*parms[3]*(4.0*CosPhi2-1.0);
        break;
      case FOURIER_SERIES_DIHEDRAL:
        // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
        // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
        // =======================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        U=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
               2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi2+4.0*(parms[2]-5.0*parms[4])*CosPhi2*CosPhi-
               8.0*(parms[3]-6.0*parms[5])*SQR(CosPhi2)+16.0*parms[4]*SQR(CosPhi2)*CosPhi-32.0*parms[5]*CUBE(CosPhi2));
        DF=0.5*(parms[0]-3.0*parms[2]+5.0*parms[4])-2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi+
           6.0*(parms[2]-5.0*parms[4])*CosPhi2-16.0*(parms[3]-6.0*parms[5])*CosPhi2*CosPhi+40.0*parms[4]*SQR(CosPhi2)-96.0*parms[5]*CosPhi2*CUBE(CosPhi);
        break;
      case FOURIER_SERIES_DIHEDRAL2:
        // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
        // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
        // =======================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        U=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
          2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
          2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
        DF=0.5*parms[0]+parms[2]*(6.0*CosPhi2-1.5)+parms[4]*(2.5-30.0*CosPhi2+40.0*SQR(CosPhi2))+
           CosPhi*(-2.0*parms[1]+parms[3]*(16.0*CosPhi2-8.0)+parms[5]*(18.0-96.0*CosPhi2+96.0*SQR(CosPhi2)));
        break;
      case FIXED_DIHEDRAL:
        U=DF=0;
        break;
      default:
        fprintf(stderr, "Undefined Torsion potential in routine 'CalculateAdsorbateTorsionForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UAdsorbateTorsion[CurrentSystem]+=U;

    // Calculate the first derivative vectors.
    d=dot_ab/rbc;
    e=dot_cd/rbc;

    dtA.x=(ds.x-CosPhi*dr.x)/r;
    dtA.y=(ds.y-CosPhi*dr.y)/r;
    dtA.z=(ds.z-CosPhi*dr.z)/r;

    dtD.x=(dr.x-CosPhi*ds.x)/s;
    dtD.y=(dr.y-CosPhi*ds.y)/s;
    dtD.z=(dr.z-CosPhi*ds.z)/s;

    dtB.x=dtA.x*(d-1.0)+e*dtD.x;
    dtB.y=dtA.y*(d-1.0)+e*dtD.y;
    dtB.z=dtA.z*(d-1.0)+e*dtD.z;

    dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
    dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
    dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

    // forces are oppositely directed to the gradient
    fa.x=-DF*dtA.x;
    fa.y=-DF*dtA.y;
    fa.z=-DF*dtA.z;

    fb.x=-DF*dtB.x;
    fb.y=-DF*dtB.y;
    fb.z=-DF*dtB.z;

    fc.x=-DF*dtC.x;
    fc.y=-DF*dtC.y;
    fc.z=-DF*dtC.z;

    fd.x=-DF*dtD.x;
    fd.y=-DF*dtD.y;
    fd.z=-DF*dtD.z;

    // add contribution to the forces
    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x+=fb.x;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y+=fb.y;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z+=fb.z;

    Adsorbates[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    Adsorbates[CurrentSystem][m].Atoms[D].Force.x+=fd.x;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.y+=fd.y;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.z+=fd.z;

    // add contribution to the stress tensor
    // Note: rbc is here because the vector was normalized before
    StrainDerivativeTensor[CurrentSystem].ax-=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;
  }
}


void CalculateCationTorsionForce(int m)
{
  int i,Type,NumberOfTorsions,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi,U,DF;
  REAL ShiftedCosPhi,ShiftedCosPhi2,ShiftedSinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  VECTOR fa,fb,fc,fd;
  REAL *parms;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfTorsions=Components[Type].NumberOfTorsions;
  for(i=0;i<NumberOfTorsions;i++)
  {
    A=Components[Type].Torsions[i].A;
    B=Components[Type].Torsions[i].B;
    C=Components[Type].Torsions[i].C;
    D=Components[Type].Torsions[i].D;
    parms=(REAL*)&Components[Type].TorsionArguments[i];

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;
    posD=Cations[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;

    Dcb.x=posC.x-posB.x;
    Dcb.y=posC.y-posB.y;
    Dcb.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
    Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

    Ddc.x=posD.x-posC.x;
    Ddc.y=posD.y-posC.y;
    Ddc.z=posD.z-posC.z;
    dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
    dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

    dr.x=Dab.x-dot_ab*Dcb.x;
    dr.y=Dab.y-dot_ab*Dcb.y;
    dr.z=Dab.z-dot_ab*Dcb.z;
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    dr.x/=r; dr.y/=r; dr.z/=r;

    ds.x=Ddc.x-dot_cd*Dcb.x;
    ds.y=Ddc.y-dot_cd*Dcb.y;
    ds.z=Ddc.z-dot_cd*Dcb.z;
    s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
    ds.x/=s; ds.y/=s; ds.z/=s;

    // compute Cos(Phi)
    // Phi is defined in protein convention Phi(trans)=Pi
    CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].TorsionType[i])
    {
      case HARMONIC_DIHEDRAL:
        // (1/2)*p_0*(phi-p_1)^2
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
        Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
        Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
        Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
        Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
        Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
        sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        Phi-=parms[1];
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
        U=0.5*parms[0]*SQR(Phi);
        DF=-parms[0]*(Phi)/SinPhi;
        break;
      case HARMONIC_COSINE_DIHEDRAL:
        // (1/2)*p_0*(cos(phi)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        U=0.5*parms[0]*SQR(CosPhi-parms[1]);
        DF=parms[0]*(CosPhi-parms[1]);
        break;
      case THREE_COSINE_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
        break;
      case MM3_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0     [kcal/mol]
        // p_1     [kcal/mol]
        // p_2     [kcal/mol]
        U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
        break;
      case CVFF_BLOCKED_DIHEDRAL:
        //
        // ========================================================================
        // p_0     [rad]
        // p_1     [K]
        // p_2     [-]
        // p_3     [rad]
        // p_4     [rad]
        U=0.0;
        DF=0.0;
        break;
      case CFF_DIHEDRAL:
        // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
        DF=-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2);
        break;
      case CFF_DIHEDRAL2:
        // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
        DF=parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi);
        break;
      case SIX_COSINE_DIHEDRAL:
        // Prod_i=0^5 p_i*cos(phi)^i
        // =========================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
        // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
        // polymer convention Phi'(trans)=0.
        U=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
               parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
        DF=-parms[1]+2.0*parms[2]*CosPhi-3.0*parms[3]*CosPhi2+4.0*parms[4]*CosPhi2*CosPhi-5.0*parms[5]*SQR(CosPhi2);
        break;
      case TRAPPE_DIHEDRAL:
        // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
        // ==========================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
        DF=parms[1]-4.0*parms[2]*CosPhi+3.0*parms[3]*(4.0*CosPhi2-1.0);
        break;
      case TRAPPE_DIHEDRAL_EXTENDED:
        // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
        // =============================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        U=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
        DF=parms[1]-3.0*parms[3]+4.0*(parms[2]-4.0*parms[4])*CosPhi+12.0*parms[3]*CosPhi2+32.0*parms[4]*CosPhi2*CosPhi;
        break;
      case MOD_TRAPPE_DIHEDRAL:
        /* Salvador modification: 16/08/2016
         add phase in cos function:
         p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
        */
        Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
        Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
        Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
        Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
        Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
        Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
        sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
        Phi-=parms[4];           // shift Phi as Phi+parms[4]
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        ShiftedCosPhi=cos(Phi);
        ShiftedSinPhi=sin(Phi);
        ShiftedCosPhi2=SQR(ShiftedCosPhi);
        U=parms[0]+parms[1]+parms[3]+(parms[1]-3.0*parms[3])*ShiftedCosPhi -2.0*parms[2]*ShiftedCosPhi2 +4.0*parms[3]*ShiftedCosPhi*ShiftedCosPhi2;
        DF=((parms[1]-3.0*parms[3])*sin(Phi) -4.0*parms[2]*ShiftedCosPhi*ShiftedSinPhi + 12.0*parms[3]*ShiftedCosPhi2*ShiftedSinPhi)/SinPhi;
        break;
      case CVFF_DIHEDRAL:
        // p_0*(1+cos(p_1*phi-p_2))
        // ========================
        // p_0/k_B [K]
        // p_1     [-]
        // p_2     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
        Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
        Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
        Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
        Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
        Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
        sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
        U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
        DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
        break;
      case OPLS_DIHEDRAL:
        // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
        // =================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
        DF=0.5*parms[1]-2.0*parms[2]*CosPhi+1.5*parms[3]*(4.0*CosPhi2-1.0);
        break;
      case FOURIER_SERIES_DIHEDRAL:
        // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
        // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
        // =======================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        U=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
               2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi2+4.0*(parms[2]-5.0*parms[4])*CosPhi2*CosPhi-
               8.0*(parms[3]-6.0*parms[5])*SQR(CosPhi2)+16.0*parms[4]*SQR(CosPhi2)*CosPhi-32.0*parms[5]*CUBE(CosPhi2));
        DF=0.5*(parms[0]-3.0*parms[2]+5.0*parms[4])-2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi+
           6.0*(parms[2]-5.0*parms[4])*CosPhi2-16.0*(parms[3]-6.0*parms[5])*CosPhi2*CosPhi+40.0*parms[4]*SQR(CosPhi2)-96.0*parms[5]*CosPhi2*CUBE(CosPhi);
        break;
      case FOURIER_SERIES_DIHEDRAL2:
        // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
        // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
        // =======================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        U=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
          2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
          2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
        DF=0.5*parms[0]+parms[2]*(6.0*CosPhi2-1.5)+parms[4]*(2.5-30.0*CosPhi2+40.0*SQR(CosPhi2))+
           CosPhi*(-2.0*parms[1]+parms[3]*(16.0*CosPhi2-8.0)+parms[5]*(18.0-96.0*CosPhi2+96.0*SQR(CosPhi2)));
        break;
      case FIXED_DIHEDRAL:
        U=DF=0;
        break;
      default:
        fprintf(stderr, "Undefined Torsion potential in routine 'CalculateCationTorsionForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UCationTorsion[CurrentSystem]+=U;

    // virial
    // pure torsion has *no* contributions to the virial

    // Calculate the first derivative vectors.
    d=dot_ab/rbc;
    e=dot_cd/rbc;

    dtA.x=(ds.x-CosPhi*dr.x)/r;
    dtA.y=(ds.y-CosPhi*dr.y)/r;
    dtA.z=(ds.z-CosPhi*dr.z)/r;

    dtD.x=(dr.x-CosPhi*ds.x)/s;
    dtD.y=(dr.y-CosPhi*ds.y)/s;
    dtD.z=(dr.z-CosPhi*ds.z)/s;

    dtB.x=dtA.x*(d-1.0)+e*dtD.x;
    dtB.y=dtA.y*(d-1.0)+e*dtD.y;
    dtB.z=dtA.z*(d-1.0)+e*dtD.z;

    dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
    dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
    dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

    // forces are oppositely directed to the gradient
    fa.x=-DF*dtA.x;
    fa.y=-DF*dtA.y;
    fa.z=-DF*dtA.z;

    fb.x=-DF*dtB.x;
    fb.y=-DF*dtB.y;
    fb.z=-DF*dtB.z;

    fc.x=-DF*dtC.x;
    fc.y=-DF*dtC.y;
    fc.z=-DF*dtC.z;

    fd.x=-DF*dtD.x;
    fd.y=-DF*dtD.y;
    fd.z=-DF*dtD.z;

    // add contribution to the forces
    Cations[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x+=fb.x;
    Cations[CurrentSystem][m].Atoms[B].Force.y+=fb.y;
    Cations[CurrentSystem][m].Atoms[B].Force.z+=fb.z;

    Cations[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Cations[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Cations[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    Cations[CurrentSystem][m].Atoms[D].Force.x+=fd.x;
    Cations[CurrentSystem][m].Atoms[D].Force.y+=fd.y;
    Cations[CurrentSystem][m].Atoms[D].Force.z+=fd.z;

    // add contribution to the stress tensor
    // Note: rbc is here because the vector was normalized before
    StrainDerivativeTensor[CurrentSystem].ax-=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;
  }
}



void CalculateAdsorbateImproperTorsionForce(int m)
{
  int i,Type,NumberOfImproperTorsions,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi,U,DF;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  VECTOR fa,fb,fc,fd;
  REAL *parms;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfImproperTorsions=Components[Type].NumberOfImproperTorsions;
  for(i=0;i<NumberOfImproperTorsions;i++)
  {
    A=Components[Type].ImproperTorsions[i].A;
    B=Components[Type].ImproperTorsions[i].B;
    C=Components[Type].ImproperTorsions[i].C;
    D=Components[Type].ImproperTorsions[i].D;
    parms=(REAL*)&Components[Type].ImproperTorsionArguments[i];

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
    posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;

    Dcb.x=posC.x-posB.x;
    Dcb.y=posC.y-posB.y;
    Dcb.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
    Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

    Ddc.x=posD.x-posC.x;
    Ddc.y=posD.y-posC.y;
    Ddc.z=posD.z-posC.z;
    dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
    dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

    dr.x=Dab.x-dot_ab*Dcb.x;
    dr.y=Dab.y-dot_ab*Dcb.y;
    dr.z=Dab.z-dot_ab*Dcb.z;
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    dr.x/=r; dr.y/=r; dr.z/=r;

    ds.x=Ddc.x-dot_cd*Dcb.x;
    ds.y=Ddc.y-dot_cd*Dcb.y;
    ds.z=Ddc.z-dot_cd*Dcb.z;
    s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
    ds.x/=s; ds.y/=s; ds.z/=s;

    // compute Cos(Phi)
    // Phi is defined in protein convention Phi(trans)=Pi
    CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].ImproperTorsionType[i])
    {
      case HARMONIC_IMPROPER_DIHEDRAL:
        // (1/2)*p_0*(phi-p_1)^2
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
        Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
        Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
        Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
        Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
        Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
        sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        Phi-=parms[1];
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
        U=0.5*parms[0]*SQR(Phi);
        DF=-parms[0]*(Phi)/SinPhi;
        break;
      case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
        // (1/2)*p_0*(cos(phi)-cos(p_1))^2
        // ===============================================
        // p_0/k_B [K]
        // p_1     [degrees]
        U=0.5*parms[0]*SQR(CosPhi-parms[1]);
        DF=parms[0]*(CosPhi-parms[1]);
        break;
      case THREE_COSINE_IMPROPER_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
        break;
      case MM3_IMPROPER_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0     [kcal/mol]
        // p_1     [kcal/mol]
        // p_2     [kcal/mol]
        U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
        break;
      case CFF_IMPROPER_DIHEDRAL:
        // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
        DF=-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2);
        break;
      case CFF_IMPROPER_DIHEDRAL2:
        // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
        DF=parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi);
        break;
      case SIX_COSINE_IMPROPER_DIHEDRAL:
        // Prod_i=0^5 p_i*cos(phi)^i
        // =========================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
        // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
        // polymer convention Phi'(trans)=0.
        U=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
               parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
        DF=-parms[1]+2.0*parms[2]*CosPhi-3.0*parms[3]*CosPhi2+4.0*parms[4]*CosPhi2*CosPhi-5.0*parms[5]*SQR(CosPhi2);
        break;
      case TRAPPE_IMPROPER_DIHEDRAL:
        // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
        // ==========================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
        DF=parms[1]-4.0*parms[2]*CosPhi+3.0*parms[3]*(4.0*CosPhi2-1.0);
        break;
      case TRAPPE_IMPROPER_DIHEDRAL_EXTENDED:
        // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
        // =============================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        U=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
        DF=parms[1]-3.0*parms[3]+4.0*(parms[2]-4.0*parms[4])*CosPhi+12.0*parms[3]*CosPhi2+32.0*parms[4]*CosPhi2*CosPhi;
        break;
      case CVFF_IMPROPER_DIHEDRAL:
        // p_0*(1+cos(p_1*phi-p_2))
        // ========================
        // p_0/k_B [K]
        // p_1     [-]
        // p_2     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
        Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
        Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
        Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
        Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
        Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
        sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
        U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
        DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
        break;
      case OPLS_IMPROPER_DIHEDRAL:
        // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
        // =================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
        DF=0.5*parms[1]-2.0*parms[2]*CosPhi+1.5*parms[3]*(4.0*CosPhi2-1.0);
        break;
      case FOURIER_SERIES_IMPROPER_DIHEDRAL:
        // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
        // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
        // =======================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        U=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
               2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi2+4.0*(parms[2]-5.0*parms[4])*CosPhi2*CosPhi-
               8.0*(parms[3]-6.0*parms[5])*SQR(CosPhi2)+16.0*parms[4]*SQR(CosPhi2)*CosPhi-32.0*parms[5]*CUBE(CosPhi2));
        DF=0.5*(parms[0]-3.0*parms[2]+5.0*parms[4])-2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi+
           6.0*(parms[2]-5.0*parms[4])*CosPhi2-16.0*(parms[3]-6.0*parms[5])*CosPhi2*CosPhi+40.0*parms[4]*SQR(CosPhi2)-96.0*parms[5]*CosPhi2*CUBE(CosPhi);
        break;
      case FOURIER_SERIES_IMPROPER_DIHEDRAL2:
        // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
        // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
        // =======================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        U=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
          2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
          2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
        DF=0.5*parms[0]+parms[2]*(6.0*CosPhi2-1.5)+parms[4]*(2.5-30.0*CosPhi2+40.0*SQR(CosPhi2))+
           CosPhi*(-2.0*parms[1]+parms[3]*(16.0*CosPhi2-8.0)+parms[5]*(18.0-96.0*CosPhi2+96.0*SQR(CosPhi2)));
        break;
      case FIXED_IMPROPER_DIHEDRAL:
        U=DF=0;
        break;
      default:
        fprintf(stderr, "Undefined Improper-Torsion potential in routine 'CalculateAdsorbateImproperTorsionForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UAdsorbateImproperTorsion[CurrentSystem]+=U;

    // virial
    // pure torsion has *no* contributions to the virial

    // Calculate the first derivative vectors.
    d=dot_ab/rbc;
    e=dot_cd/rbc;

    dtA.x=(ds.x-CosPhi*dr.x)/r;
    dtA.y=(ds.y-CosPhi*dr.y)/r;
    dtA.z=(ds.z-CosPhi*dr.z)/r;

    dtD.x=(dr.x-CosPhi*ds.x)/s;
    dtD.y=(dr.y-CosPhi*ds.y)/s;
    dtD.z=(dr.z-CosPhi*ds.z)/s;

    dtB.x=dtA.x*(d-1.0)+e*dtD.x;
    dtB.y=dtA.y*(d-1.0)+e*dtD.y;
    dtB.z=dtA.z*(d-1.0)+e*dtD.z;

    dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
    dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
    dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

    // forces are oppositely directed to the gradient
    fa.x=-DF*dtA.x;
    fa.y=-DF*dtA.y;
    fa.z=-DF*dtA.z;

    fb.x=-DF*dtB.x;
    fb.y=-DF*dtB.y;
    fb.z=-DF*dtB.z;

    fc.x=-DF*dtC.x;
    fc.y=-DF*dtC.y;
    fc.z=-DF*dtC.z;

    fd.x=-DF*dtD.x;
    fd.y=-DF*dtD.y;
    fd.z=-DF*dtD.z;

    // add contribution to the forces
    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x+=fb.x;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y+=fb.y;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z+=fb.z;

    Adsorbates[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    Adsorbates[CurrentSystem][m].Atoms[D].Force.x+=fd.x;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.y+=fd.y;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.z+=fd.z;

    // add contribution to the stress tensor
    // Note: rbc is here because the vector was normalized before
    StrainDerivativeTensor[CurrentSystem].ax-=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;
  }
}


void CalculateCationImproperTorsionForce(int m)
{
  int i,Type,NumberOfImproperTorsions,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi,U,DF;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  VECTOR fa,fb,fc,fd;
  REAL *parms;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfImproperTorsions=Components[Type].NumberOfImproperTorsions;
  for(i=0;i<NumberOfImproperTorsions;i++)
  {
    A=Components[Type].ImproperTorsions[i].A;
    B=Components[Type].ImproperTorsions[i].B;
    C=Components[Type].ImproperTorsions[i].C;
    D=Components[Type].ImproperTorsions[i].D;
    parms=(REAL*)&Components[Type].ImproperTorsionArguments[i];

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;
    posD=Cations[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;

    Dcb.x=posC.x-posB.x;
    Dcb.y=posC.y-posB.y;
    Dcb.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
    Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

    Ddc.x=posD.x-posC.x;
    Ddc.y=posD.y-posC.y;
    Ddc.z=posD.z-posC.z;
    dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
    dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

    dr.x=Dab.x-dot_ab*Dcb.x;
    dr.y=Dab.y-dot_ab*Dcb.y;
    dr.z=Dab.z-dot_ab*Dcb.z;
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    dr.x/=r; dr.y/=r; dr.z/=r;

    ds.x=Ddc.x-dot_cd*Dcb.x;
    ds.y=Ddc.y-dot_cd*Dcb.y;
    ds.z=Ddc.z-dot_cd*Dcb.z;
    s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
    ds.x/=s; ds.y/=s; ds.z/=s;

    // compute Cos(Phi)
    // Phi is defined in protein convention Phi(trans)=Pi
    CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].ImproperTorsionType[i])
    {
      case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
        U=0.5*parms[0]*SQR(CosPhi-cos(parms[1]));
        DF=parms[0]*(CosPhi-cos(parms[1]));
        break;
      case HARMONIC_IMPROPER_DIHEDRAL:
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
        Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
        Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
        Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
        Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
        Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
        sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        Phi-=parms[1];
        Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
        U=0.5*parms[0]*SQR(Phi);
        DF=-parms[0]*(Phi)/SinPhi;
        break;
      case THREE_COSINE_IMPROPER_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
        break;
      case MM3_IMPROPER_DIHEDRAL:
        // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
        // ========================================================================
        // p_0     [kcal/mol]
        // p_1     [kcal/mol]
        // p_2     [kcal/mol]
        U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
        DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
        break;
      case CFF_IMPROPER_DIHEDRAL:
        // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
        DF=-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2);
        break;
      case CFF_IMPROPER_DIHEDRAL2:
        // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
        // ======================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
        DF=parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi);
        break;
      case SIX_COSINE_IMPROPER_DIHEDRAL:
        // Prod_i=0^5 p_i*cos(phi)^i
        // =========================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
        // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
        // polymer convention Phi'(trans)=0.
        U=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
               parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
        DF=-parms[1]+2.0*parms[2]*CosPhi-3.0*parms[3]*CosPhi2+4.0*parms[4]*CosPhi2*CosPhi-5.0*parms[5]*SQR(CosPhi2);
        break;
      case TRAPPE_IMPROPER_DIHEDRAL:
        // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
        // ==========================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
        DF=parms[1]-4.0*parms[2]*CosPhi+3.0*parms[3]*(4.0*CosPhi2-1.0);
        break;
      case TRAPPE_IMPROPER_DIHEDRAL_EXTENDED:
        // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
        // =============================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        U=parms[0]-parms[2]+parms[4]+(parms[1]-3.0*parms[3])*CosPhi+(2.0*parms[2]-8.0*parms[4])*CosPhi2+4.0*parms[3]*CosPhi2*CosPhi+8.0*parms[4]*SQR(CosPhi2);
        DF=parms[1]-3.0*parms[3]+4.0*(parms[2]-4.0*parms[4])*CosPhi+12.0*parms[3]*CosPhi2+32.0*parms[4]*CosPhi2*CosPhi;
        break;
      case CVFF_IMPROPER_DIHEDRAL:
        // p_0*(1+cos(p_1*phi-p_2))
        // ========================
        // p_0/k_B [K]
        // p_1     [-]
        // p_2     [degrees]
        // potential defined in terms of 'phi' and therefore contains a singularity
        // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
        // same direction as Rbc, and negative otherwise
        Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
        Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
        Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
        Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
        Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
        Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
        sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
        U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
        DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
        break;
      case OPLS_IMPROPER_DIHEDRAL:
        // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
        // =================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
        DF=0.5*parms[1]-2.0*parms[2]*CosPhi+1.5*parms[3]*(4.0*CosPhi2-1.0);
        break;
      case FOURIER_SERIES_IMPROPER_DIHEDRAL:
        // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
        // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
        // =======================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        U=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
               2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi2+4.0*(parms[2]-5.0*parms[4])*CosPhi2*CosPhi-
               8.0*(parms[3]-6.0*parms[5])*SQR(CosPhi2)+16.0*parms[4]*SQR(CosPhi2)*CosPhi-32.0*parms[5]*CUBE(CosPhi2));
        DF=0.5*(parms[0]-3.0*parms[2]+5.0*parms[4])-2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi+
           6.0*(parms[2]-5.0*parms[4])*CosPhi2-16.0*(parms[3]-6.0*parms[5])*CosPhi2*CosPhi+40.0*parms[4]*SQR(CosPhi2)-96.0*parms[5]*CosPhi2*CUBE(CosPhi);
        break;
      case FOURIER_SERIES_IMPROPER_DIHEDRAL2:
        // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
        // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
        // =======================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        // p_3/k_B [K]
        // p_4/k_B [K]
        // p_5/k_B [K]
        U=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
          2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
          2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
        DF=0.5*parms[0]+parms[2]*(6.0*CosPhi2-1.5)+parms[4]*(2.5-30.0*CosPhi2+40.0*SQR(CosPhi2))+
           CosPhi*(-2.0*parms[1]+parms[3]*(16.0*CosPhi2-8.0)+parms[5]*(18.0-96.0*CosPhi2+96.0*SQR(CosPhi2)));
        break;
      case FIXED_IMPROPER_DIHEDRAL:
        U=DF=0;
        break;
      default:
        fprintf(stderr, "Undefined Improper-Torsion potential in routine 'CalculateCationImproperTorsionForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UCationImproperTorsion[CurrentSystem]+=U;

    // virial
    // pure torsion has *no* contributions to the virial

    // Calculate the first derivative vectors.
    d=dot_ab/rbc;
    e=dot_cd/rbc;

    dtA.x=(ds.x-CosPhi*dr.x)/r;
    dtA.y=(ds.y-CosPhi*dr.y)/r;
    dtA.z=(ds.z-CosPhi*dr.z)/r;

    dtD.x=(dr.x-CosPhi*ds.x)/s;
    dtD.y=(dr.y-CosPhi*ds.y)/s;
    dtD.z=(dr.z-CosPhi*ds.z)/s;

    dtB.x=dtA.x*(d-1.0)+e*dtD.x;
    dtB.y=dtA.y*(d-1.0)+e*dtD.y;
    dtB.z=dtA.z*(d-1.0)+e*dtD.z;

    dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
    dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
    dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

    // forces are oppositely directed to the gradient
    fa.x=-DF*dtA.x;
    fa.y=-DF*dtA.y;
    fa.z=-DF*dtA.z;

    fb.x=-DF*dtB.x;
    fb.y=-DF*dtB.y;
    fb.z=-DF*dtB.z;

    fc.x=-DF*dtC.x;
    fc.y=-DF*dtC.y;
    fc.z=-DF*dtC.z;

    fd.x=-DF*dtD.x;
    fd.y=-DF*dtD.y;
    fd.z=-DF*dtD.z;

    // add contribution to the forces
    Cations[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x+=fb.x;
    Cations[CurrentSystem][m].Atoms[B].Force.y+=fb.y;
    Cations[CurrentSystem][m].Atoms[B].Force.z+=fb.z;

    Cations[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Cations[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Cations[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    Cations[CurrentSystem][m].Atoms[D].Force.x+=fd.x;
    Cations[CurrentSystem][m].Atoms[D].Force.y+=fd.y;
    Cations[CurrentSystem][m].Atoms[D].Force.z+=fd.z;

    // add contribution to the stress tensor
    // Note: rbc is here because the vector was normalized before
    StrainDerivativeTensor[CurrentSystem].ax-=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;
  }
}


void CalculateAdsorbateBondBondForce(int m)
{
  int i,A,B,C;
  int NumberOfBondBonds,Type;
  REAL *parms,gamma,gammc;
  REAL energy;
  REAL rab,rbc;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc,fa,fc;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBondBonds=Components[Type].NumberOfBondBonds;
  for(i=0;i<NumberOfBondBonds;i++)
  {
    A=Components[Type].BondBonds[i].A;
    B=Components[Type].BondBonds[i].B;
    C=Components[Type].BondBonds[i].C;
    parms=(REAL*)&Components[Type].BondBondArguments[i];

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

    Rab.x=posA.x-posB.x;
    Rab.y=posA.y-posB.y;
    Rab.z=posA.z-posB.z;
    rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));

    Rbc.x=posC.x-posB.x;
    Rbc.y=posC.y-posB.y;
    Rbc.z=posC.z-posB.z;
    rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));

    switch(Components[Type].BondBondType[i])
    {
      case CVFF_BOND_BOND_CROSS:
      case CFF_BOND_BOND_CROSS:
        // p_0*(rab-p_1)*(rbc-p_2)
        // =======================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        energy=parms[0]*(rab-parms[1])*(rbc-parms[2]);
        gamma=-parms[0]*(rbc-parms[2])/rab;
        gammc=-parms[0]*(rab-parms[1])/rbc;
        break;
      default:
        fprintf(stderr, "Undefined Bond-Bond potential in routine 'CalculateAdsorbateBondBondForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    UAdsorbateBondBond[CurrentSystem]+=energy;

    // forces
    fa.x=gamma*Rab.x;
    fa.y=gamma*Rab.y;
    fa.z=gamma*Rab.z;

    fc.x=gammc*Rbc.x;
    fc.y=gammc*Rbc.y;
    fc.z=gammc*Rbc.z;

    // add contribution to the forces
    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x-=(fa.x+fc.x);
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y-=(fa.y+fc.y);
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z-=(fa.z+fc.z);

    Adsorbates[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax-=Rab.x*fa.x+Rbc.x*fc.x;
    StrainDerivativeTensor[CurrentSystem].bx-=Rab.y*fa.x+Rbc.y*fc.x;
    StrainDerivativeTensor[CurrentSystem].cx-=Rab.z*fa.x+Rbc.z*fc.x;

    StrainDerivativeTensor[CurrentSystem].ay-=Rab.x*fa.y+Rbc.x*fc.y;
    StrainDerivativeTensor[CurrentSystem].by-=Rab.y*fa.y+Rbc.y*fc.y;
    StrainDerivativeTensor[CurrentSystem].cy-=Rab.z*fa.y+Rbc.z*fc.y;

    StrainDerivativeTensor[CurrentSystem].az-=Rab.x*fa.z+Rbc.x*fc.z;
    StrainDerivativeTensor[CurrentSystem].bz-=Rab.y*fa.z+Rbc.y*fc.z;
    StrainDerivativeTensor[CurrentSystem].cz-=Rab.z*fa.z+Rbc.z*fc.z;
  }
}

void CalculateCationBondBondForce(int m)
{
  int i,A,B,C;
  int NumberOfBondBonds,Type;
  REAL *parms,gamma,gammc;
  REAL energy;
  REAL rab,rbc;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc,fa,fc;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfBondBonds=Components[Type].NumberOfBondBonds;
  for(i=0;i<NumberOfBondBonds;i++)
  {
    A=Components[Type].BondBonds[i].A;
    B=Components[Type].BondBonds[i].B;
    C=Components[Type].BondBonds[i].C;
    parms=(REAL*)&Components[Type].BondBondArguments[i];

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;

    Rab.x=posA.x-posB.x;
    Rab.y=posA.y-posB.y;
    Rab.z=posA.z-posB.z;
    rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));

    Rbc.x=posC.x-posB.x;
    Rbc.y=posC.y-posB.y;
    Rbc.z=posC.z-posB.z;
    rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));

    switch(Components[Type].BondBondType[i])
    {
      case CVFF_BOND_BOND_CROSS:
      case CFF_BOND_BOND_CROSS:
        // p_0*(rab-p_1)*(rbc-p_2)
        // =======================
        // p_0/k_B [K/A^2]
        // p_1     [A]
        // p_2     [A]
        energy=parms[0]*(rab-parms[1])*(rbc-parms[2]);
        gamma=-parms[0]*(rbc-parms[2])/rab;
        gammc=-parms[0]*(rab-parms[1])/rbc;
        break;
      default:
        fprintf(stderr, "Undefined Bond-Bond potential in routine 'CalculateCationBondBondForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    UCationBondBond[CurrentSystem]+=energy;

    // forces
    fa.x=gamma*Rab.x;
    fa.y=gamma*Rab.y;
    fa.z=gamma*Rab.z;

    fc.x=gammc*Rbc.x;
    fc.y=gammc*Rbc.y;
    fc.z=gammc*Rbc.z;

    // add contribution to the forces
    Cations[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x-=(fa.x+fc.x);
    Cations[CurrentSystem][m].Atoms[B].Force.y-=(fa.y+fc.y);
    Cations[CurrentSystem][m].Atoms[B].Force.z-=(fa.z+fc.z);

    Cations[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Cations[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Cations[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax-=Rab.x*fa.x+Rbc.x*fc.x;
    StrainDerivativeTensor[CurrentSystem].bx-=Rab.y*fa.x+Rbc.y*fc.x;
    StrainDerivativeTensor[CurrentSystem].cx-=Rab.z*fa.x+Rbc.z*fc.x;

    StrainDerivativeTensor[CurrentSystem].ay-=Rab.x*fa.y+Rbc.x*fc.y;
    StrainDerivativeTensor[CurrentSystem].by-=Rab.y*fa.y+Rbc.y*fc.y;
    StrainDerivativeTensor[CurrentSystem].cy-=Rab.z*fa.y+Rbc.z*fc.y;

    StrainDerivativeTensor[CurrentSystem].az-=Rab.x*fa.z+Rbc.x*fc.z;
    StrainDerivativeTensor[CurrentSystem].bz-=Rab.y*fa.z+Rbc.y*fc.z;
    StrainDerivativeTensor[CurrentSystem].cz-=Rab.z*fa.z+Rbc.z*fc.z;
  }
}



void CalculateAdsorbateBondBendForce(int m)
{
  int i,Type,NumberOfBondBends,A,B,C;
  REAL *parms,gamma,gamsa,gamsc,pterm,vterm;
  REAL CosTheta,Theta,SinTheta;
  REAL rab,rbc,rac;
  POINT posA,posB,posC;
  VECTOR Rab,Rbc,Rac,fa,fc;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBondBends=Components[Type].NumberOfBondBends;
  for(i=0;i<NumberOfBondBends;i++)
  {
    A=Components[Type].BondBends[i].A;
    B=Components[Type].BondBends[i].B;
    C=Components[Type].BondBends[i].C;
    parms=(REAL*)&Components[Type].BondBendArguments[i];

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

    Rab.x=posA.x-posB.x;
    Rab.y=posA.y-posB.y;
    Rab.z=posA.z-posB.z;
    rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
    Rab.x/=rab;
    Rab.y/=rab;
    Rab.z/=rab;

    Rbc.x=posC.x-posB.x;
    Rbc.y=posC.y-posB.y;
    Rbc.z=posC.z-posB.z;
    rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
    Rbc.x/=rbc;
    Rbc.y/=rbc;
    Rbc.z/=rbc;

    Rac.x=posC.x-posA.x;
    Rac.y=posC.y-posA.y;
    Rac.z=posC.z-posA.z;
    rac=sqrt(SQR(Rac.x)+SQR(Rac.y)+SQR(Rac.z));
    Rac.x/=rac;
    Rac.y/=rac;
    Rac.z/=rac;

    CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
    CosTheta=SIGN(MIN2(fabs(CosTheta),(REAL)1.0),CosTheta);
    Theta=acos(CosTheta);
    SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));

    switch(Components[Type].BondBendType[i])
    {
      case CVFF_BOND_BEND_CROSS:
      case CFF_BOND_BEND_CROSS:
        // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
        // =========================================
        // p_0     [degrees]
        // p_1/k_B [K/A/rad]
        // p_2     [A]
        // p_3/k_B [K/A/rad]
        // p_4     [A]
        pterm=(Theta-parms[0])*(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]));
        vterm=0.0;
        gamma=(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]))/SinTheta;
        gamsa=-parms[1]*(Theta-parms[0]);
        gamsc=-parms[3]*(Theta-parms[0]);
        break;
      case MM3_BOND_BEND_CROSS:
        // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
        // =====================================
        // p_0     [mdyne/rad]
        // p_1     [A]
        // p_2     [A]
        // p_3     [degrees]
        pterm=parms[0]*((rab-parms[1])+(rbc-parms[2]))*RAD2DEG*(Theta-parms[3]);
        vterm=0.0;
        gamma=parms[0]*RAD2DEG*((rab-parms[1])+(rbc-parms[2]))/SinTheta;
        gamsa=-parms[0]*RAD2DEG*(Theta-parms[3]);
        gamsc=-parms[0]*RAD2DEG*(Theta-parms[3]);
        break;
      case TRUNCATED_HARMONIC:
        // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
        // ================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        pterm=0.5*parms[0]*SQR(Theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8));
        vterm=-8.0*pterm*(pow(rab,8)+pow(rbc,8))/pow(parms[2],8);
        gamma=(parms[0]*(Theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8)))/SinTheta;
        gamsa=(8.0*pterm/pow(parms[2],8))*pow(rab,7);
        gamsc=(8.0*pterm/pow(parms[2],8))*pow(rbc,7);
        break;
      case SCREENED_HARMONIC:
        // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        // p_3     [A]
        pterm=0.5*parms[0]*SQR(Theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3]));
        vterm=-pterm*(rab/parms[2]+rbc/parms[3]);
        gamma=(parms[0]*(Theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3])))/SinTheta;
        gamsa=(pterm/parms[2]);
        gamsc=(pterm/parms[3]);
        break;
      case SCREENED_VESSAL:
        // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
        // ============================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        // p_3     [A]
        pterm=(parms[0]/(8.0*SQR(Theta-M_PI)))*SQR(SQR(parms[1]-M_PI)-SQR(Theta-M_PI))
              *exp(-(rab/parms[2]+rbc/parms[3]));
        vterm=-pterm*(rab/parms[2]+rbc/parms[3]);
        gamma=(parms[0]/(2.0*SQR(Theta-M_PI))*(SQR(parms[1]-M_PI)-SQR(Theta-M_PI))*(Theta-M_PI)
              *exp(-(rab/parms[2]+rbc/parms[3])))/SinTheta;
        gamsa=(pterm/parms[2]);
        gamsc=(pterm/parms[3]);
        break;
      case TRUNCATED_VESSAL:
        // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
        //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
        // ============================================================================
        // p_0/k_B [K/rad^(4+p_2)]
        // p_1     [degrees]
        // p_2     [-]
        // p_3     [A]
        pterm=parms[0]*(pow(Theta,parms[2])*SQR(Theta-parms[1])*SQR(Theta+parms[1]-2.0*M_PI)
              -0.5*parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*SQR(Theta-parms[1])*pow(M_PI-parms[1],(REAL)3.0))
              *exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8));
        vterm=-8.0*pterm*(pow(rab,8)+pow(rbc,8))/pow(parms[3],8);
        gamma=(parms[0]*(pow(Theta,(parms[2]-1.0))*(Theta-parms[1])*(Theta+parms[1]-2.0*M_PI)
               *((parms[2]+4.0)*SQR(Theta)-2.0*M_PI*(parms[2]+2.0)*Theta+parms[2]*parms[1]*
               (2.0*M_PI-parms[1]))-parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*(Theta-parms[1])*
               pow(M_PI-parms[1],3))*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8)))/SinTheta;
        gamsa=(8.0*pterm/pow(parms[3],8))*pow(rab,7);
        gamsc=(8.0*pterm/pow(parms[3],8))*pow(rbc,7);
        break;
      default:
        fprintf(stderr, "Undefined Bond-Bend potential in routine 'CalculateAdsorbateBondBendForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UAdsorbateBondBend[CurrentSystem]+=pterm;

    // forces
    fa.x=gamma*(Rbc.x-Rab.x*CosTheta)/rab+gamsa*Rab.x;
    fa.y=gamma*(Rbc.y-Rab.y*CosTheta)/rab+gamsa*Rab.y;
    fa.z=gamma*(Rbc.z-Rab.z*CosTheta)/rab+gamsa*Rab.z;

    fc.x=gamma*(Rab.x-Rbc.x*CosTheta)/rbc+gamsc*Rbc.x;
    fc.y=gamma*(Rab.y-Rbc.y*CosTheta)/rbc+gamsc*Rbc.y;
    fc.z=gamma*(Rab.z-Rbc.z*CosTheta)/rbc+gamsc*Rbc.z;

    // add contribution to the forces
    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x-=(fa.x+fc.x);
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y-=(fa.y+fc.y);
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z-=(fa.z+fc.z);

    Adsorbates[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax-=(rab*Rab.x*fa.x+rbc*Rbc.x*fc.x);
    StrainDerivativeTensor[CurrentSystem].ay-=(rab*Rab.x*fa.y+rbc*Rbc.x*fc.y);
    StrainDerivativeTensor[CurrentSystem].az-=(rab*Rab.x*fa.z+rbc*Rbc.x*fc.z);

    StrainDerivativeTensor[CurrentSystem].bx-=(rab*Rab.y*fa.x+rbc*Rbc.y*fc.x);
    StrainDerivativeTensor[CurrentSystem].by-=(rab*Rab.y*fa.y+rbc*Rbc.y*fc.y);
    StrainDerivativeTensor[CurrentSystem].bz-=(rab*Rab.y*fa.z+rbc*Rbc.y*fc.z);

    StrainDerivativeTensor[CurrentSystem].cx-=(rab*Rab.z*fa.x+rbc*Rbc.z*fc.x);
    StrainDerivativeTensor[CurrentSystem].cy-=(rab*Rab.z*fa.y+rbc*Rbc.z*fc.y);
    StrainDerivativeTensor[CurrentSystem].cz-=(rab*Rab.z*fa.z+rbc*Rbc.z*fc.z);
  }
}

void CalculateCationBondBendForce(int m)
{
  int i,Type,NumberOfBondBends,A,B,C;
  REAL *parms,gamma,gamsa,gamsc,pterm,vterm;
  REAL CosTheta,Theta,SinTheta;
  REAL rab,rbc,rac;
  POINT posA,posB,posC;
  VECTOR Rab,Rbc,Rac,fa,fc;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfBondBends=Components[Type].NumberOfBondBends;
  for(i=0;i<NumberOfBondBends;i++)
  {
    A=Components[Type].BondBends[i].A;
    B=Components[Type].BondBends[i].B;
    C=Components[Type].BondBends[i].C;
    parms=(REAL*)&Components[Type].BondBendArguments[i];

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;

    Rab.x=posA.x-posB.x;
    Rab.y=posA.y-posB.y;
    Rab.z=posA.z-posB.z;
    rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
    Rab.x/=rab;
    Rab.y/=rab;
    Rab.z/=rab;

    Rbc.x=posC.x-posB.x;
    Rbc.y=posC.y-posB.y;
    Rbc.z=posC.z-posB.z;
    rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
    Rbc.x/=rbc;
    Rbc.y/=rbc;
    Rbc.z/=rbc;

    Rac.x=posC.x-posA.x;
    Rac.y=posC.y-posA.y;
    Rac.z=posC.z-posA.z;
    rac=sqrt(SQR(Rac.x)+SQR(Rac.y)+SQR(Rac.z));
    Rac.x/=rac;
    Rac.y/=rac;
    Rac.z/=rac;

    CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
    CosTheta=SIGN(MIN2(fabs(CosTheta),(REAL)1.0),CosTheta);
    Theta=acos(CosTheta);
    SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));

    switch(Components[Type].BondBendType[i])
    {
      case CVFF_BOND_BEND_CROSS:
      case CFF_BOND_BEND_CROSS:
        // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
        // =========================================
        // p_0     [degrees]
        // p_1/k_B [K/A/rad]
        // p_2     [A]
        // p_3/k_B [K/A/rad]
        // p_4     [A]
        pterm=(Theta-parms[0])*(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]));
        vterm=0.0;
        gamma=(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]))/SinTheta;
        gamsa=-parms[1]*(Theta-parms[0]);
        gamsc=-parms[3]*(Theta-parms[0]);
        break;
      case MM3_BOND_BEND_CROSS:
        // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
        // =====================================
        // p_0     [mdyne/rad]
        // p_1     [A]
        // p_2     [A]
        // p_3     [degrees]
        pterm=parms[0]*((rab-parms[1])+(rbc-parms[2]))*RAD2DEG*(Theta-parms[3]);
        vterm=0.0;
        gamma=parms[0]*RAD2DEG*((rab-parms[1])+(rbc-parms[2]))/SinTheta;
        gamsa=-parms[0]*RAD2DEG*(Theta-parms[3]);
        gamsc=-parms[0]*RAD2DEG*(Theta-parms[3]);
        break;
      case TRUNCATED_HARMONIC:
        // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
        // ================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        pterm=0.5*parms[0]*SQR(Theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8));
        vterm=-8.0*pterm*(pow(rab,8)+pow(rbc,8))/pow(parms[2],8);
        gamma=(parms[0]*(Theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8)))/SinTheta;
        gamsa=(8.0*pterm/pow(parms[2],8))*pow(rab,7);
        gamsc=(8.0*pterm/pow(parms[2],8))*pow(rbc,7);
        break;
      case SCREENED_HARMONIC:
        // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
        // ===============================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        // p_3     [A]
        pterm=0.5*parms[0]*SQR(Theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3]));
        vterm=-pterm*(rab/parms[2]+rbc/parms[3]);
        gamma=(parms[0]*(Theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3])))/SinTheta;
        gamsa=(pterm/parms[2]);
        gamsc=(pterm/parms[3]);
        break;
      case SCREENED_VESSAL:
        // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
        // ============================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [degrees]
        // p_2     [A]
        // p_3     [A]
        pterm=(parms[0]/(8.0*SQR(Theta-M_PI)))*SQR(SQR(parms[1]-M_PI)-SQR(Theta-M_PI))
              *exp(-(rab/parms[2]+rbc/parms[3]));
        vterm=-pterm*(rab/parms[2]+rbc/parms[3]);
        gamma=(parms[0]/(2.0*SQR(Theta-M_PI))*(SQR(parms[1]-M_PI)-SQR(Theta-M_PI))*(Theta-M_PI)
              *exp(-(rab/parms[2]+rbc/parms[3])))/SinTheta;
        gamsa=(pterm/parms[2]);
        gamsc=(pterm/parms[3]);
        break;
      case TRUNCATED_VESSAL:
        // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
        //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
        // ============================================================================
        // p_0/k_B [K/rad^(4+p_2)]
        // p_1     [degrees]
        // p_2     [-]
        // p_3     [A]
        pterm=parms[0]*(pow(Theta,parms[2])*SQR(Theta-parms[1])*SQR(Theta+parms[1]-2.0*M_PI)
              -0.5*parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*SQR(Theta-parms[1])*pow(M_PI-parms[1],3))
              *exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8));
        vterm=-8.0*pterm*(pow(rab,8)+pow(rbc,8))/pow(parms[3],8);
        gamma=(parms[0]*(pow(Theta,(parms[2]-1.0))*(Theta-parms[1])*(Theta+parms[1]-2.0*M_PI)
               *((parms[2]+4.0)*SQR(Theta)-2.0*M_PI*(parms[2]+2.0)*Theta+parms[2]*parms[1]*
               (2.0*M_PI-parms[1]))-parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*(Theta-parms[1])*
               pow(M_PI-parms[1],3))*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8)))/SinTheta;
        gamsa=(8.0*pterm/pow(parms[3],8))*pow(rab,7);
        gamsc=(8.0*pterm/pow(parms[3],8))*pow(rbc,7);
        break;
      default:
        fprintf(stderr, "Undefined Bond-Bend potential in routine 'CalculateCationBondBendForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UCationBondBend[CurrentSystem]+=pterm;

    // forces
    fa.x=gamma*(Rbc.x-Rab.x*CosTheta)/rab+gamsa*Rab.x;
    fa.y=gamma*(Rbc.y-Rab.y*CosTheta)/rab+gamsa*Rab.y;
    fa.z=gamma*(Rbc.z-Rab.z*CosTheta)/rab+gamsa*Rab.z;

    fc.x=gamma*(Rab.x-Rbc.x*CosTheta)/rbc+gamsc*Rbc.x;
    fc.y=gamma*(Rab.y-Rbc.y*CosTheta)/rbc+gamsc*Rbc.y;
    fc.z=gamma*(Rab.z-Rbc.z*CosTheta)/rbc+gamsc*Rbc.z;

    // add contribution to the forces
    Cations[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x-=(fa.x+fc.x);
    Cations[CurrentSystem][m].Atoms[B].Force.y-=(fa.y+fc.y);
    Cations[CurrentSystem][m].Atoms[B].Force.z-=(fa.z+fc.z);

    Cations[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Cations[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Cations[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax-=(rab*Rab.x*fa.x+rbc*Rbc.x*fc.x);
    StrainDerivativeTensor[CurrentSystem].ay-=(rab*Rab.x*fa.y+rbc*Rbc.x*fc.y);
    StrainDerivativeTensor[CurrentSystem].az-=(rab*Rab.x*fa.z+rbc*Rbc.x*fc.z);

    StrainDerivativeTensor[CurrentSystem].bx-=(rab*Rab.y*fa.x+rbc*Rbc.y*fc.x);
    StrainDerivativeTensor[CurrentSystem].by-=(rab*Rab.y*fa.y+rbc*Rbc.y*fc.y);
    StrainDerivativeTensor[CurrentSystem].bz-=(rab*Rab.y*fa.z+rbc*Rbc.y*fc.z);

    StrainDerivativeTensor[CurrentSystem].cx-=(rab*Rab.z*fa.x+rbc*Rbc.z*fc.x);
    StrainDerivativeTensor[CurrentSystem].cy-=(rab*Rab.z*fa.y+rbc*Rbc.z*fc.y);
    StrainDerivativeTensor[CurrentSystem].cz-=(rab*Rab.z*fa.z+rbc*Rbc.z*fc.z);
  }
}

void CalculateAdsorbateBendBendForce(int m)
{
  int i,Type,NumberOfBendBends,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL energy,rab,rbc,rbd;
  VECTOR Dab,Dbc,Dbd;
  REAL dot_abc,dot_abd;
  VECTOR fa,fb,fc,fd;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL DTheta1,DTheta2;
  REAL *parms;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfBendBends=Components[Type].NumberOfBendBends;
  for(i=0;i<NumberOfBendBends;i++)
  {
    A=Components[Type].BendBends[i].A;
    B=Components[Type].BendBends[i].B;
    C=Components[Type].BendBends[i].C;
    D=Components[Type].BendBends[i].D;
    parms=(REAL*)&Components[Type].BendBendArguments[i];

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
    posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;
    rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));
    Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;

    Dbc.x=posC.x-posB.x;
    Dbc.y=posC.y-posB.y;
    Dbc.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
    Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

    Dbd.x=posD.x-posB.x;
    Dbd.y=posD.y-posB.y;
    Dbd.z=posD.z-posB.z;
    rbd=sqrt(SQR(Dbd.x)+SQR(Dbd.y)+SQR(Dbd.z));
    Dbd.x/=rbd; Dbd.y/=rbd; Dbd.z/=rbd;

    dot_abc=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
    CosTheta1=dot_abc;
    CosTheta1=SIGN(MIN2(fabs(CosTheta1),(REAL)1.0),CosTheta1);
    Theta1=acos(CosTheta1);
    SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

    dot_abd=Dab.x*Dbd.x+Dab.y*Dbd.y+Dab.z*Dbd.z;
    CosTheta2=dot_abd;
    CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
    Theta2=acos(CosTheta2);
    SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

    switch(Components[Type].BendBendType[i])
    {
      case CVFF_BEND_BEND_CROSS:
      case CFF_BEND_BEND_CROSS:
        // p_0*(Theta1-p_1)*(Theta2-p_2)
        // ===================================
        // p_0/k_B [K/rad^2)]
        // p_1     [degrees]
        // p_2     [degrees]
        energy=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
        DTheta1=parms[0]*(Theta2-parms[2])/SinTheta1;
        DTheta2=parms[0]*(Theta1-parms[1])/SinTheta2;
        break;
      case MM3_BEND_BEND_CROSS:
        // -p_0*(Theta1-p_1)*(Theta2-p_2)
        // ===================================
        // p_0     [mdyne A/rad^2]
        // p_1     [degrees]
        // p_2     [degrees]
        energy=parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])*(Theta2-parms[2]);
        DTheta1=parms[0]*SQR(RAD2DEG)*(Theta2-parms[2])/SinTheta1;
        DTheta2=parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])/SinTheta2;
        break;
      default:
        fprintf(stderr, "Undefined Bend-Bend potential in routine 'CalculateAdsorbateBendBendForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

      // energy
    UAdsorbateBendBend[CurrentSystem]+=energy;

    // forces bend
    fa.x=DTheta1*(Dbc.x-CosTheta1*Dab.x)/rab;
    fa.y=DTheta1*(Dbc.y-CosTheta1*Dab.y)/rab;
    fa.z=DTheta1*(Dbc.z-CosTheta1*Dab.z)/rab;

    fc.x=DTheta1*(Dab.x-CosTheta1*Dbc.x)/rbc;
    fc.y=DTheta1*(Dab.y-CosTheta1*Dbc.y)/rbc;
    fc.z=DTheta1*(Dab.z-CosTheta1*Dbc.z)/rbc;

    fb.x=-DTheta1*(fa.x+fc.x);
    fb.y=-DTheta1*(fa.y+fc.y);
    fb.z=-DTheta1*(fa.z+fc.z);

    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x-=fa.x+fc.x;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y-=fa.y+fc.y;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z-=fa.z+fc.z;

    Adsorbates[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    StrainDerivativeTensor[CurrentSystem].ax-=rab*Dab.x*fa.x+rbc*Dbc.x*fc.x;
    StrainDerivativeTensor[CurrentSystem].bx-=rab*Dab.y*fa.x+rbc*Dbc.y*fc.x;
    StrainDerivativeTensor[CurrentSystem].cx-=rab*Dab.z*fa.x+rbc*Dbc.z*fc.x;

    StrainDerivativeTensor[CurrentSystem].ay-=rab*Dab.x*fa.y+rbc*Dbc.x*fc.y;
    StrainDerivativeTensor[CurrentSystem].by-=rab*Dab.y*fa.y+rbc*Dbc.y*fc.y;
    StrainDerivativeTensor[CurrentSystem].cy-=rab*Dab.z*fa.y+rbc*Dbc.z*fc.y;

    StrainDerivativeTensor[CurrentSystem].az-=rab*Dab.x*fa.z+rbc*Dbc.x*fc.z;
    StrainDerivativeTensor[CurrentSystem].bz-=rab*Dab.y*fa.z+rbc*Dbc.y*fc.z;
    StrainDerivativeTensor[CurrentSystem].cz-=rab*Dab.z*fa.z+rbc*Dbc.z*fc.z;

    fa.x=DTheta2*(Dbd.x-CosTheta2*Dab.x)/rab;
    fa.y=DTheta2*(Dbd.y-CosTheta2*Dab.y)/rab;
    fa.z=DTheta2*(Dbd.z-CosTheta2*Dab.z)/rab;

    fd.x=DTheta2*(Dab.x-CosTheta2*Dbd.x)/rbd;
    fd.y=DTheta2*(Dab.y-CosTheta2*Dbd.y)/rbd;
    fd.z=DTheta2*(Dab.z-CosTheta2*Dbd.z)/rbd;

    fb.x=-DTheta2*(fa.x+fd.x);
    fb.y=-DTheta2*(fa.y+fd.y);
    fb.z=-DTheta2*(fa.z+fd.z);

    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x-=fa.x+fd.x;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y-=fa.y+fd.y;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z-=fa.z+fd.z;

    Adsorbates[CurrentSystem][m].Atoms[D].Force.x+=fd.x;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.y+=fd.y;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.z+=fd.z;

    StrainDerivativeTensor[CurrentSystem].ax-=rab*Dab.x*fa.x+rbd*Dbd.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=rab*Dab.y*fa.x+rbd*Dbd.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=rab*Dab.z*fa.x+rbd*Dbd.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=rab*Dab.x*fa.y+rbd*Dbd.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=rab*Dab.y*fa.y+rbd*Dbd.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=rab*Dab.z*fa.y+rbd*Dbd.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=rab*Dab.x*fa.z+rbd*Dbd.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=rab*Dab.y*fa.z+rbd*Dbd.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=rab*Dab.z*fa.z+rbd*Dbd.z*fd.z;
  }
}

void CalculateCationBendBendForce(int m)
{
  int i,Type,NumberOfBendBends,A,B,C,D;
  POINT posA,posB,posC,posD;
  REAL energy,rab,rbc,rbd;
  VECTOR Dab,Dbc,Dbd;
  REAL dot_abc,dot_abd;
  VECTOR fa,fb,fc,fd;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL DTheta1,DTheta2;
  REAL *parms;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfBendBends=Components[Type].NumberOfBendBends;
  for(i=0;i<NumberOfBendBends;i++)
  {
    A=Components[Type].BendBends[i].A;
    B=Components[Type].BendBends[i].B;
    C=Components[Type].BendBends[i].C;
    D=Components[Type].BendBends[i].D;
    parms=(REAL*)&Components[Type].BendBendArguments[i];

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;
    posD=Cations[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;
    rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));
    Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;

    Dbc.x=posC.x-posB.x;
    Dbc.y=posC.y-posB.y;
    Dbc.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
    Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

    Dbd.x=posD.x-posB.x;
    Dbd.y=posD.y-posB.y;
    Dbd.z=posD.z-posB.z;
    rbd=sqrt(SQR(Dbd.x)+SQR(Dbd.y)+SQR(Dbd.z));
    Dbd.x/=rbd; Dbd.y/=rbd; Dbd.z/=rbd;

    dot_abc=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
    CosTheta1=dot_abc;
    CosTheta1=SIGN(MIN2(fabs(CosTheta1),(REAL)1.0),CosTheta1);
    Theta1=acos(CosTheta1);
    SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

    dot_abd=Dab.x*Dbd.x+Dab.y*Dbd.y+Dab.z*Dbd.z;
    CosTheta2=dot_abd;
    CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
    Theta2=acos(CosTheta2);
    SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

    switch(Components[Type].BendBendType[i])
    {
      case CVFF_BEND_BEND_CROSS:
      case CFF_BEND_BEND_CROSS:
        // p_0*(Theta1-p_1)*(Theta2-p_2)
        // ===================================
        // p_0/k_B [K/rad^2)]
        // p_1     [degrees]
        // p_2     [degrees]
        energy=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
        DTheta1=parms[0]*(Theta2-parms[2])/SinTheta1;
        DTheta2=parms[0]*(Theta1-parms[1])/SinTheta2;
        break;
      case MM3_BEND_BEND_CROSS:
        // -p_0*(Theta1-p_1)*(Theta2-p_2)
        // ===================================
        // p_0     [mdyne A/rad^2]
        // p_1     [degrees]
        // p_2     [degrees]
        energy=parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])*(Theta2-parms[2]);
        DTheta1=parms[0]*SQR(RAD2DEG)*(Theta2-parms[2])/SinTheta1;
        DTheta2=parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])/SinTheta2;
        break;
      default:
        fprintf(stderr, "Undefined Bend-Bend potential in routine 'CalculateCationBendBendForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UCationBendBend[CurrentSystem]+=energy;

    // forces bend
    fa.x=DTheta1*(Dbc.x-CosTheta1*Dab.x)/rab;
    fa.y=DTheta1*(Dbc.y-CosTheta1*Dab.y)/rab;
    fa.z=DTheta1*(Dbc.z-CosTheta1*Dab.z)/rab;

    fc.x=DTheta1*(Dab.x-CosTheta1*Dbc.x)/rbc;
    fc.y=DTheta1*(Dab.y-CosTheta1*Dbc.y)/rbc;
    fc.z=DTheta1*(Dab.z-CosTheta1*Dbc.z)/rbc;

    fb.x=-DTheta1*(fa.x+fc.x);
    fb.y=-DTheta1*(fa.y+fc.y);
    fb.z=-DTheta1*(fa.z+fc.z);

    Cations[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x-=fa.x+fc.x;
    Cations[CurrentSystem][m].Atoms[B].Force.y-=fa.y+fc.y;
    Cations[CurrentSystem][m].Atoms[B].Force.z-=fa.z+fc.z;

    Cations[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Cations[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Cations[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    StrainDerivativeTensor[CurrentSystem].ax-=rab*Dab.x*fa.x+rbc*Dbc.x*fc.x;
    StrainDerivativeTensor[CurrentSystem].bx-=rab*Dab.y*fa.x+rbc*Dbc.y*fc.x;
    StrainDerivativeTensor[CurrentSystem].cx-=rab*Dab.z*fa.x+rbc*Dbc.z*fc.x;

    StrainDerivativeTensor[CurrentSystem].ay-=rab*Dab.x*fa.y+rbc*Dbc.x*fc.y;
    StrainDerivativeTensor[CurrentSystem].by-=rab*Dab.y*fa.y+rbc*Dbc.y*fc.y;
    StrainDerivativeTensor[CurrentSystem].cy-=rab*Dab.z*fa.y+rbc*Dbc.z*fc.y;

    StrainDerivativeTensor[CurrentSystem].az-=rab*Dab.x*fa.z+rbc*Dbc.x*fc.z;
    StrainDerivativeTensor[CurrentSystem].bz-=rab*Dab.y*fa.z+rbc*Dbc.y*fc.z;
    StrainDerivativeTensor[CurrentSystem].cz-=rab*Dab.z*fa.z+rbc*Dbc.z*fc.z;

    fa.x=DTheta2*(Dbd.x-CosTheta2*Dab.x)/rab;
    fa.y=DTheta2*(Dbd.y-CosTheta2*Dab.y)/rab;
    fa.z=DTheta2*(Dbd.z-CosTheta2*Dab.z)/rab;

    fd.x=DTheta2*(Dab.x-CosTheta2*Dbd.x)/rbd;
    fd.y=DTheta2*(Dab.y-CosTheta2*Dbd.y)/rbd;
    fd.z=DTheta2*(Dab.z-CosTheta2*Dbd.z)/rbd;

    fb.x=-DTheta2*(fa.x+fd.x);
    fb.y=-DTheta2*(fa.y+fd.y);
    fb.z=-DTheta2*(fa.z+fd.z);

    Cations[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x-=fa.x+fd.x;
    Cations[CurrentSystem][m].Atoms[B].Force.y-=fa.y+fd.y;
    Cations[CurrentSystem][m].Atoms[B].Force.z-=fa.z+fd.z;

    Cations[CurrentSystem][m].Atoms[D].Force.x+=fd.x;
    Cations[CurrentSystem][m].Atoms[D].Force.y+=fd.y;
    Cations[CurrentSystem][m].Atoms[D].Force.z+=fd.z;

    StrainDerivativeTensor[CurrentSystem].ax-=rab*Dab.x*fa.x+rbd*Dbd.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=rab*Dab.y*fa.x+rbd*Dbd.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=rab*Dab.z*fa.x+rbd*Dbd.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=rab*Dab.x*fa.y+rbd*Dbd.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=rab*Dab.y*fa.y+rbd*Dbd.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=rab*Dab.z*fa.y+rbd*Dbd.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=rab*Dab.x*fa.z+rbd*Dbd.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=rab*Dab.y*fa.z+rbd*Dbd.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=rab*Dab.z*fa.z+rbd*Dbd.z*fd.z;
  }
}

void CalculateAdsorbateBendTorsionForce(int m)
{
  int i,A,B,C,D,Type;
  POINT posA,posB,posC,posD;
  REAL d,e, energy,rab,rbc,rcd;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2,DCos;
  VECTOR dtA,dtB,dtC,dtD,fa,fb,fc,fd,Pb,Pc;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL DTheta1,DTheta2,sign,Phi,SinPhi;
  REAL *parms;

  Type=Adsorbates[CurrentSystem][m].Type;
  for(i=0;i<Components[Type].NumberOfBendTorsions;i++)
  {
    A=Components[Type].BendTorsions[i].A;
    B=Components[Type].BendTorsions[i].B;
    C=Components[Type].BendTorsions[i].C;
    D=Components[Type].BendTorsions[i].D;
    parms=(REAL*)&Components[Type].BendTorsionArguments[i];

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
    posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;
    rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

    Dbc.x=posC.x-posB.x;
    Dbc.y=posC.y-posB.y;
    Dbc.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
    Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

    Dcd.x=posD.x-posC.x;
    Dcd.y=posD.y-posC.y;
    Dcd.z=posD.z-posC.z;
    rcd=sqrt(SQR(Dcd.x)+SQR(Dcd.y)+SQR(Dcd.z));

    dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
    CosTheta1=dot_ab/rab;
    CosTheta1=SIGN(MIN2(fabs(CosTheta1),(REAL)1.0),CosTheta1);
    Theta1=acos(CosTheta1);
    SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

    dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;
    CosTheta2=-dot_cd/rcd;
    CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
    Theta2=acos(CosTheta2);
    SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

    dr.x=Dab.x-dot_ab*Dbc.x;
    dr.y=Dab.y-dot_ab*Dbc.y;
    dr.z=Dab.z-dot_ab*Dbc.z;
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    dr.x/=r; dr.y/=r; dr.z/=r;

    ds.x=Dcd.x-dot_cd*Dbc.x;
    ds.y=Dcd.y-dot_cd*Dbc.y;
    ds.z=Dcd.z-dot_cd*Dbc.z;
    s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
    ds.x/=s; ds.y/=s; ds.z/=s;

    // compute Cos(Phi)
    // Phi is defined in protein convention Phi(trans)=Pi
    CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].BendTorsionType[i])
    {
      case CVFF_BEND_TORSION_CROSS:
      case CFF_BEND_TORSION_CROSS:
        // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
        // =====================================================================================
        // p_0/k_B [K/rad^3]
        // p_1     [degrees]
        // p_2     [degrees]
        energy=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi;
        DCos=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
        DTheta1=parms[0]*(Theta2-parms[2])*CosPhi/SinTheta1;
        DTheta2=parms[0]*(Theta1-parms[1])*CosPhi/SinTheta2;
        break;
      case SMOOTHED_DIHEDRAL:
        // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2)]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [-]
        // p_2     [degrees]
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
        energy=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2);
        DCos=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2)/SinPhi;
        DTheta1=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
        DTheta2=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*SmoothingDerivative(Theta2)/SinTheta2;
        break;
      case SMOOTHED_THREE_COSINE_DIHEDRAL:
        // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        energy=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
               Smoothing(Theta1)*Smoothing(Theta2);
        DCos=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0)*Smoothing(Theta1)*Smoothing(Theta2);
        DTheta1=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                 SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
        DTheta2=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                 SmoothingDerivative(Theta2)*Smoothing(Theta1)/SinTheta2;
        break;
      case NICHOLAS_DIHEDRAL:
        // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        energy=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
               Smoothing(Theta1);
        DCos=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0)*Smoothing(Theta1);
        DTheta1=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                SmoothingDerivative(Theta1)/SinTheta1;
        DTheta2=0.0;
        break;
      case SMOOTHED_CFF_DIHEDRAL:
        // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        energy=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
        DCos=(-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
        DTheta1=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*
                SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
        DTheta2=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*
                Smoothing(Theta1)*SmoothingDerivative(Theta2)/SinTheta2;
        break;
      case SMOOTHED_CFF_DIHEDRAL2:
        // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        energy=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*Smoothing(Theta1)*Smoothing(Theta2);
        DCos=(parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi))*Smoothing(Theta1)*Smoothing(Theta2);
        DTheta1=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*
                SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
        DTheta2=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*
                Smoothing(Theta1)*SmoothingDerivative(Theta2)/SinTheta2;
        break;
      case SMOOTHED_CFF_BEND_TORSION_CROSS:
        // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K/rad^3]
        // p_1     [degrees]
        // p_2     [degrees]
        energy=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi*Smoothing(Theta1)*Smoothing(Theta2);
        DCos=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*Smoothing(Theta1)*Smoothing(Theta2);
        DTheta1=CosPhi*parms[0]*(Theta2-parms[2])*Smoothing(Theta2)*(Smoothing(Theta1)+(Theta1-parms[1])*SmoothingDerivative(Theta1))/SinTheta1;
        DTheta2=CosPhi*parms[0]*(Theta1-parms[1])*Smoothing(Theta1)*(Smoothing(Theta2)+(Theta2-parms[2])*SmoothingDerivative(Theta2))/SinTheta2;
        break;
      default:
        fprintf(stderr, "Undefined Bend-Torsion potential in routine 'CalculateAdsorbateBendTorsionForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UAdsorbateBendTorsion[CurrentSystem]+=energy;

    // Calculate the first derivative vectors.
    d=dot_ab/rbc;
    e=dot_cd/rbc;

    dtA.x=(ds.x-CosPhi*dr.x)/r;
    dtA.y=(ds.y-CosPhi*dr.y)/r;
    dtA.z=(ds.z-CosPhi*dr.z)/r;

    dtD.x=(dr.x-CosPhi*ds.x)/s;
    dtD.y=(dr.y-CosPhi*ds.y)/s;
    dtD.z=(dr.z-CosPhi*ds.z)/s;

    dtB.x=dtA.x*(d-1.0)+e*dtD.x;
    dtB.y=dtA.y*(d-1.0)+e*dtD.y;
    dtB.z=dtA.z*(d-1.0)+e*dtD.z;

    dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
    dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
    dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

    // forces are oppositely directed to the gradient
    fa.x=-DCos*dtA.x;
    fa.y=-DCos*dtA.y;
    fa.z=-DCos*dtA.z;

    fb.x=-DCos*dtB.x;
    fb.y=-DCos*dtB.y;
    fb.z=-DCos*dtB.z;

    fc.x=-DCos*dtC.x;
    fc.y=-DCos*dtC.y;
    fc.z=-DCos*dtC.z;

    fd.x=-DCos*dtD.x;
    fd.y=-DCos*dtD.y;
    fd.z=-DCos*dtD.z;

    // forces torsion
    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x+=fb.x;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y+=fb.y;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z+=fb.z;

    Adsorbates[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    Adsorbates[CurrentSystem][m].Atoms[D].Force.x+=fd.x;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.y+=fd.y;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.z+=fd.z;

    // add contribution to the stress tensor
    // Note: rbc is here because the vector was normalized before
    StrainDerivativeTensor[CurrentSystem].ax-=Dab.x*fa.x+Dbc.x*rbc*(fc.x+fd.x)+Dcd.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=Dab.y*fa.x+Dbc.y*rbc*(fc.x+fd.x)+Dcd.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=Dab.z*fa.x+Dbc.z*rbc*(fc.x+fd.x)+Dcd.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=Dab.x*fa.y+Dbc.x*rbc*(fc.y+fd.y)+Dcd.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=Dab.y*fa.y+Dbc.y*rbc*(fc.y+fd.y)+Dcd.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=Dab.z*fa.y+Dbc.z*rbc*(fc.y+fd.y)+Dcd.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=Dab.x*fa.z+Dbc.x*rbc*(fc.z+fd.z)+Dcd.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=Dab.y*fa.z+Dbc.y*rbc*(fc.z+fd.z)+Dcd.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=Dab.z*fa.z+Dbc.z*rbc*(fc.z+fd.z)+Dcd.z*fd.z;

    Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;
    Dcd.x/=rcd; Dcd.y/=rcd; Dcd.z/=rcd;

    // forces bends
    fa.x=DTheta1*(Dbc.x-Dab.x*CosTheta1)/rab;
    fa.y=DTheta1*(Dbc.y-Dab.y*CosTheta1)/rab;
    fa.z=DTheta1*(Dbc.z-Dab.z*CosTheta1)/rab;

    fb.x=DTheta2*(Dcd.x+Dbc.x*CosTheta2)/rbc;
    fb.y=DTheta2*(Dcd.y+Dbc.y*CosTheta2)/rbc;
    fb.z=DTheta2*(Dcd.z+Dbc.z*CosTheta2)/rbc;

    fc.x=DTheta1*(Dab.x-Dbc.x*CosTheta1)/rbc;
    fc.y=DTheta1*(Dab.y-Dbc.y*CosTheta1)/rbc;
    fc.z=DTheta1*(Dab.z-Dbc.z*CosTheta1)/rbc;

    fd.x=DTheta2*(-Dbc.x-Dcd.x*CosTheta2)/rcd;
    fd.y=DTheta2*(-Dbc.y-Dcd.y*CosTheta2)/rcd;
    fd.z=DTheta2*(-Dbc.z-Dcd.z*CosTheta2)/rcd;

    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x+=fb.x-fa.x-fc.x;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y+=fb.y-fa.y-fc.y;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z+=fb.z-fa.z-fc.z;

    Adsorbates[CurrentSystem][m].Atoms[C].Force.x+=fc.x-fb.x-fd.x;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.y+=fc.y-fb.y-fd.y;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.z+=fc.z-fb.z-fd.z;

    Adsorbates[CurrentSystem][m].Atoms[D].Force.x+=fd.x;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.y+=fd.y;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.z+=fd.z;

    StrainDerivativeTensor[CurrentSystem].ax-=rab*Dab.x*fa.x+rbc*Dbc.x*fc.x+rbc*Dbc.x*fb.x+rcd*Dcd.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=rab*Dab.y*fa.x+rbc*Dbc.y*fc.x+rbc*Dbc.y*fb.x+rcd*Dcd.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=rab*Dab.z*fa.x+rbc*Dbc.z*fc.x+rbc*Dbc.z*fb.x+rcd*Dcd.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=rab*Dab.x*fa.y+rbc*Dbc.x*fc.y+rbc*Dbc.x*fb.y+rcd*Dcd.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=rab*Dab.y*fa.y+rbc*Dbc.y*fc.y+rbc*Dbc.y*fb.y+rcd*Dcd.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=rab*Dab.z*fa.y+rbc*Dbc.z*fc.y+rbc*Dbc.z*fb.y+rcd*Dcd.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=rab*Dab.x*fa.z+rbc*Dbc.x*fc.z+rbc*Dbc.x*fb.z+rcd*Dcd.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=rab*Dab.y*fa.z+rbc*Dbc.y*fc.z+rbc*Dbc.y*fb.z+rcd*Dcd.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=rab*Dab.z*fa.z+rbc*Dbc.z*fc.z+rbc*Dbc.z*fb.z+rcd*Dcd.z*fd.z;
  }
}

void CalculateCationBendTorsionForce(int m)
{
  int i,A,B,C,D,Type;
  POINT posA,posB,posC,posD;
  REAL d,e, energy,rab,rbc,rcd;
  VECTOR Dab,Dbc,Dcd,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2,DCos;
  VECTOR dtA,dtB,dtC,dtD,fa,fb,fc,fd,Pb,Pc;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL DTheta1,DTheta2,sign,Phi,SinPhi;
  REAL *parms;

  Type=Cations[CurrentSystem][m].Type;
  for(i=0;i<Components[Type].NumberOfBendTorsions;i++)
  {
    A=Components[Type].BendTorsions[i].A;
    B=Components[Type].BendTorsions[i].B;
    C=Components[Type].BendTorsions[i].C;
    D=Components[Type].BendTorsions[i].D;
    parms=(REAL*)&Components[Type].BendTorsionArguments[i];

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;
    posD=Cations[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;
    Dab=ApplyBoundaryCondition(Dab);
    rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

    Dbc.x=posC.x-posB.x;
    Dbc.y=posC.y-posB.y;
    Dbc.z=posC.z-posB.z;
    Dbc=ApplyBoundaryCondition(Dbc);
    rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
    Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;

    Dcd.x=posD.x-posC.x;
    Dcd.y=posD.y-posC.y;
    Dcd.z=posD.z-posC.z;
    Dcd=ApplyBoundaryCondition(Dcd);
    rcd=sqrt(SQR(Dcd.x)+SQR(Dcd.y)+SQR(Dcd.z));

    dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
    CosTheta1=dot_ab/rab;
    CosTheta1=SIGN(MIN2(fabs(CosTheta1),(REAL)1.0),CosTheta1);
    Theta1=acos(CosTheta1);
    SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));

    dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;
    CosTheta2=-dot_cd/rcd;
    CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
    Theta2=acos(CosTheta2);
    SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));

    dr.x=Dab.x-dot_ab*Dbc.x;
    dr.y=Dab.y-dot_ab*Dbc.y;
    dr.z=Dab.z-dot_ab*Dbc.z;
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    dr.x/=r; dr.y/=r; dr.z/=r;

    ds.x=Dcd.x-dot_cd*Dbc.x;
    ds.y=Dcd.y-dot_cd*Dbc.y;
    ds.z=Dcd.z-dot_cd*Dbc.z;
    s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
    ds.x/=s; ds.y/=s; ds.z/=s;

    // compute Cos(Phi)
    // Phi is defined in protein convention Phi(trans)=Pi
    CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].BendTorsionType[i])
    {
      case CVFF_BEND_TORSION_CROSS:
      case CFF_BEND_TORSION_CROSS:
        // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
        // =====================================================================================
        // p_0/k_B [K/rad^3]
        // p_1     [degrees]
        // p_2     [degrees]
        energy=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi;
        DCos=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
        DTheta1=parms[0]*(Theta2-parms[2])*CosPhi/SinTheta1;
        DTheta2=parms[0]*(Theta1-parms[1])*CosPhi/SinTheta2;
        break;
      case SMOOTHED_DIHEDRAL:
        // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2)]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K/rad^2]
        // p_1     [-]
        // p_2     [degrees]
        Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
        Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
        Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
        Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
        Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
        Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
        sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
              +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
        Phi=SIGN(acos(CosPhi),sign);
        SinPhi=sin(Phi);
        SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
        energy=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2);
        DCos=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2)/SinPhi;
        DTheta1=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
        DTheta2=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*SmoothingDerivative(Theta2)/SinTheta2;
        break;
      case SMOOTHED_THREE_COSINE_DIHEDRAL:
        // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        energy=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
               Smoothing(Theta1)*Smoothing(Theta2);
        DCos=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0)*Smoothing(Theta1)*Smoothing(Theta2);
        DTheta1=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                 SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
        DTheta2=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                 SmoothingDerivative(Theta2)*Smoothing(Theta1)/SinTheta2;
        break;
      case NICHOLAS_DIHEDRAL:
        // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        energy=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
               Smoothing(Theta1);
        DCos=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0)*Smoothing(Theta1);
        DTheta1=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                SmoothingDerivative(Theta1)/SinTheta1;
        DTheta2=0.0;
        break;
      case SMOOTHED_CFF_DIHEDRAL:
        // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        energy=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
        DCos=(-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
        DTheta1=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*
                SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
        DTheta2=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*
                Smoothing(Theta1)*SmoothingDerivative(Theta2)/SinTheta2;
        break;
      case SMOOTHED_CFF_DIHEDRAL2:
        // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K]
        // p_1/k_B [K]
        // p_2/k_B [K]
        energy=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*Smoothing(Theta1)*Smoothing(Theta2);
        DCos=(parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi))*Smoothing(Theta1)*Smoothing(Theta2);
        DTheta1=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*
                SmoothingDerivative(Theta1)*Smoothing(Theta2)/SinTheta1;
        DTheta2=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*
                Smoothing(Theta1)*SmoothingDerivative(Theta2)/SinTheta2;
        break;
      case SMOOTHED_CFF_BEND_TORSION_CROSS:
        // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
        // ======================================================================================
        // p_0/k_B [K/rad^3]
        // p_1     [degrees]
        // p_2     [degrees]
        energy=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi*Smoothing(Theta1)*Smoothing(Theta2);
        DCos=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*Smoothing(Theta1)*Smoothing(Theta2);
        DTheta1=CosPhi*parms[0]*(Theta2-parms[2])*Smoothing(Theta2)*(Smoothing(Theta1)+(Theta1-parms[1])*SmoothingDerivative(Theta1))/SinTheta1;
        DTheta2=CosPhi*parms[0]*(Theta1-parms[1])*Smoothing(Theta1)*(Smoothing(Theta2)+(Theta2-parms[2])*SmoothingDerivative(Theta2))/SinTheta2;
        break;
      default:
        fprintf(stderr, "Undefined Bend-Torsion potential in routine 'CalculateCationBendTorsionForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UCationBendTorsion[CurrentSystem]+=energy;

    // Calculate the first derivative vectors.
    d=dot_ab/rbc;
    e=dot_cd/rbc;

    dtA.x=(ds.x-CosPhi*dr.x)/r;
    dtA.y=(ds.y-CosPhi*dr.y)/r;
    dtA.z=(ds.z-CosPhi*dr.z)/r;

    dtD.x=(dr.x-CosPhi*ds.x)/s;
    dtD.y=(dr.y-CosPhi*ds.y)/s;
    dtD.z=(dr.z-CosPhi*ds.z)/s;

    dtB.x=dtA.x*(d-1.0)+e*dtD.x;
    dtB.y=dtA.y*(d-1.0)+e*dtD.y;
    dtB.z=dtA.z*(d-1.0)+e*dtD.z;

    dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
    dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
    dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

    // forces are oppositely directed to the gradient
    fa.x=-DCos*dtA.x;
    fa.y=-DCos*dtA.y;
    fa.z=-DCos*dtA.z;

    fb.x=-DCos*dtB.x;
    fb.y=-DCos*dtB.y;
    fb.z=-DCos*dtB.z;

    fc.x=-DCos*dtC.x;
    fc.y=-DCos*dtC.y;
    fc.z=-DCos*dtC.z;

    fd.x=-DCos*dtD.x;
    fd.y=-DCos*dtD.y;
    fd.z=-DCos*dtD.z;

    // forces torsion
    Cations[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x+=fb.x;
    Cations[CurrentSystem][m].Atoms[B].Force.y+=fb.y;
    Cations[CurrentSystem][m].Atoms[B].Force.z+=fb.z;

    Cations[CurrentSystem][m].Atoms[C].Force.x+=fc.x;
    Cations[CurrentSystem][m].Atoms[C].Force.y+=fc.y;
    Cations[CurrentSystem][m].Atoms[C].Force.z+=fc.z;

    Cations[CurrentSystem][m].Atoms[D].Force.x+=fd.x;
    Cations[CurrentSystem][m].Atoms[D].Force.y+=fd.y;
    Cations[CurrentSystem][m].Atoms[D].Force.z+=fd.z;

    // add contribution to the stress tensor
    // Note: rbc is here because the vector was normalized before
    StrainDerivativeTensor[CurrentSystem].ax-=Dab.x*fa.x+Dbc.x*rbc*(fc.x+fd.x)+Dcd.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=Dab.y*fa.x+Dbc.y*rbc*(fc.x+fd.x)+Dcd.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=Dab.z*fa.x+Dbc.z*rbc*(fc.x+fd.x)+Dcd.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=Dab.x*fa.y+Dbc.x*rbc*(fc.y+fd.y)+Dcd.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=Dab.y*fa.y+Dbc.y*rbc*(fc.y+fd.y)+Dcd.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=Dab.z*fa.y+Dbc.z*rbc*(fc.y+fd.y)+Dcd.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=Dab.x*fa.z+Dbc.x*rbc*(fc.z+fd.z)+Dcd.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=Dab.y*fa.z+Dbc.y*rbc*(fc.z+fd.z)+Dcd.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=Dab.z*fa.z+Dbc.z*rbc*(fc.z+fd.z)+Dcd.z*fd.z;

    Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;
    Dcd.x/=rcd; Dcd.y/=rcd; Dcd.z/=rcd;

    // forces bends
    fa.x=DTheta1*(Dbc.x-Dab.x*CosTheta1)/rab;
    fa.y=DTheta1*(Dbc.y-Dab.y*CosTheta1)/rab;
    fa.z=DTheta1*(Dbc.z-Dab.z*CosTheta1)/rab;

    fb.x=DTheta2*(Dcd.x+Dbc.x*CosTheta2)/rbc;
    fb.y=DTheta2*(Dcd.y+Dbc.y*CosTheta2)/rbc;
    fb.z=DTheta2*(Dcd.z+Dbc.z*CosTheta2)/rbc;

    fc.x=DTheta1*(Dab.x-Dbc.x*CosTheta1)/rbc;
    fc.y=DTheta1*(Dab.y-Dbc.y*CosTheta1)/rbc;
    fc.z=DTheta1*(Dab.z-Dbc.z*CosTheta1)/rbc;

    fd.x=DTheta2*(-Dbc.x-Dcd.x*CosTheta2)/rcd;
    fd.y=DTheta2*(-Dbc.y-Dcd.y*CosTheta2)/rcd;
    fd.z=DTheta2*(-Dbc.z-Dcd.z*CosTheta2)/rcd;

    Cations[CurrentSystem][m].Atoms[A].Force.x+=fa.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=fa.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=fa.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x+=fb.x-fa.x-fc.x;
    Cations[CurrentSystem][m].Atoms[B].Force.y+=fb.y-fa.y-fc.y;
    Cations[CurrentSystem][m].Atoms[B].Force.z+=fb.z-fa.z-fc.z;

    Cations[CurrentSystem][m].Atoms[C].Force.x+=fc.x-fb.x-fd.x;
    Cations[CurrentSystem][m].Atoms[C].Force.y+=fc.y-fb.y-fd.y;
    Cations[CurrentSystem][m].Atoms[C].Force.z+=fc.z-fb.z-fd.z;

    Cations[CurrentSystem][m].Atoms[D].Force.x+=fd.x;
    Cations[CurrentSystem][m].Atoms[D].Force.y+=fd.y;
    Cations[CurrentSystem][m].Atoms[D].Force.z+=fd.z;

    StrainDerivativeTensor[CurrentSystem].ax-=rab*Dab.x*fa.x+rbc*Dbc.x*fc.x+rbc*Dbc.x*fb.x+rcd*Dcd.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=rab*Dab.y*fa.x+rbc*Dbc.y*fc.x+rbc*Dbc.y*fb.x+rcd*Dcd.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=rab*Dab.z*fa.x+rbc*Dbc.z*fc.x+rbc*Dbc.z*fb.x+rcd*Dcd.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=rab*Dab.x*fa.y+rbc*Dbc.x*fc.y+rbc*Dbc.x*fb.y+rcd*Dcd.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=rab*Dab.y*fa.y+rbc*Dbc.y*fc.y+rbc*Dbc.y*fb.y+rcd*Dcd.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=rab*Dab.z*fa.y+rbc*Dbc.z*fc.y+rbc*Dbc.z*fb.y+rcd*Dcd.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=rab*Dab.x*fa.z+rbc*Dbc.x*fc.z+rbc*Dbc.x*fb.z+rcd*Dcd.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=rab*Dab.y*fa.z+rbc*Dbc.y*fc.z+rbc*Dbc.y*fb.z+rcd*Dcd.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=rab*Dab.z*fa.z+rbc*Dbc.z*fc.z+rbc*Dbc.z*fb.z+rcd*Dcd.z*fd.z;
  }
}

void CalculateAdsorbateBondTorsionForce(int m)
{
  int i,A,B,C,D,Type;
  POINT posA,posB,posC,posD;
  REAL d,e, energy,rab,rbc,rcd,temp;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2,DCos;
  VECTOR dtA,dtB,dtC,dtD,fa,fb,fc,fd;
  REAL gamsa,gamsb,gamsc;
  REAL *parms;

  Type=Adsorbates[CurrentSystem][m].Type;
  for(i=0;i<Components[Type].NumberOfBondTorsions;i++)
  {
    A=Components[Type].BondTorsions[i].A;
    B=Components[Type].BondTorsions[i].B;
    C=Components[Type].BondTorsions[i].C;
    D=Components[Type].BondTorsions[i].D;
    parms=(REAL*)&Components[Type].BondTorsionArguments[i];

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
    posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;
    rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

    Dcb.x=posC.x-posB.x;
    Dcb.y=posC.y-posB.y;
    Dcb.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
    Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

    Ddc.x=posD.x-posC.x;
    Ddc.y=posD.y-posC.y;
    Ddc.z=posD.z-posC.z;
    rcd=sqrt(SQR(Ddc.x)+SQR(Ddc.y)+SQR(Ddc.z));

    dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
    dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

    dr.x=Dab.x-dot_ab*Dcb.x;
    dr.y=Dab.y-dot_ab*Dcb.y;
    dr.z=Dab.z-dot_ab*Dcb.z;
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    dr.x/=r; dr.y/=r; dr.z/=r;

    ds.x=Ddc.x-dot_cd*Dcb.x;
    ds.y=Ddc.y-dot_cd*Dcb.y;
    ds.z=Ddc.z-dot_cd*Dcb.z;
    s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
    ds.x/=s; ds.y/=s; ds.z/=s;

    // compute Cos(Phi)
    // Phi is defined in protein convention Phi(trans)=Pi
    CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].BondTorsionType[i])
    {
      case MM3_BOND_TORSION_CROSS:
        // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
        // =====================================================================================
        // p_0     [kcal/A mole]
        // p_1     [kcal/A mole]
        // p_2     [kcal/A mole]
        // p_3     [A]
        temp=(rbc-parms[3]);
        energy=parms[0]*temp*CosPhi+parms[1]*temp*(2.0*CosPhi2-1.0)+parms[2]*temp*(4.0*CosPhi2*CosPhi-3.0*CosPhi);
        DCos=parms[0]*temp+4.0*parms[1]*temp*CosPhi+parms[2]*temp*(12.0*CosPhi2-3.0);
        gamsa=0.0;
        gamsb=-(parms[0]*CosPhi+parms[1]*(2.0*CosPhi2-1.0)+parms[2]*(4.0*CosPhi2*CosPhi-3.0*CosPhi));
        gamsc=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Bond-Torsion potential in routine 'CalculateAdsorbateBondTorsionForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UAdsorbateBondTorsion[CurrentSystem]+=energy;

    // Calculate the first derivative vectors.
    d=dot_ab/rbc;
    e=dot_cd/rbc;

    dtA.x=(ds.x-CosPhi*dr.x)/r;
    dtA.y=(ds.y-CosPhi*dr.y)/r;
    dtA.z=(ds.z-CosPhi*dr.z)/r;

    dtD.x=(dr.x-CosPhi*ds.x)/s;
    dtD.y=(dr.y-CosPhi*ds.y)/s;
    dtD.z=(dr.z-CosPhi*ds.z)/s;

    dtB.x=dtA.x*(d-1.0)+e*dtD.x;
    dtB.y=dtA.y*(d-1.0)+e*dtD.y;
    dtB.z=dtA.z*(d-1.0)+e*dtD.z;

    dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
    dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
    dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

    Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;
    Ddc.x/=rcd; Ddc.y/=rcd; Ddc.z/=rcd;

    fa.x=DCos*dtA.x-gamsa*Dab.x;
    fa.y=DCos*dtA.y-gamsa*Dab.y;
    fa.z=DCos*dtA.z-gamsa*Dab.z;

    fb.x=DCos*dtB.x+gamsb*Dcb.x+gamsa*Dab.x;
    fb.y=DCos*dtB.y+gamsb*Dcb.y+gamsa*Dab.y;
    fb.z=DCos*dtB.z+gamsb*Dcb.z+gamsa*Dab.z;

    fc.x=DCos*dtC.x-gamsb*Dcb.x+gamsc*Ddc.x;
    fc.y=DCos*dtC.y-gamsb*Dcb.y+gamsc*Ddc.y;
    fc.z=DCos*dtC.z-gamsb*Dcb.z+gamsc*Ddc.z;

    fd.x=DCos*dtD.x-gamsc*Ddc.x;
    fd.y=DCos*dtD.y-gamsc*Ddc.y;
    fd.z=DCos*dtD.z-gamsc*Ddc.z;

    // forces torsion
    Adsorbates[CurrentSystem][m].Atoms[A].Force.x-=fa.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y-=fa.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z-=fa.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x-=fb.x;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y-=fb.y;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z-=fb.z;

    Adsorbates[CurrentSystem][m].Atoms[C].Force.x-=fc.x;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.y-=fc.y;
    Adsorbates[CurrentSystem][m].Atoms[C].Force.z-=fc.z;

    Adsorbates[CurrentSystem][m].Atoms[D].Force.x-=fd.x;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.y-=fd.y;
    Adsorbates[CurrentSystem][m].Atoms[D].Force.z-=fd.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax+=rbc*Dcb.x*(fc.x+fd.x)+rab*Dab.x*fa.x+rcd*Ddc.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].ay+=rbc*Dcb.x*(fc.y+fd.y)+rab*Dab.x*fa.y+rcd*Ddc.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].az+=rbc*Dcb.x*(fc.z+fd.z)+rab*Dab.x*fa.z+rcd*Ddc.x*fd.z;

    StrainDerivativeTensor[CurrentSystem].bx+=rbc*Dcb.y*(fc.x+fd.x)+rab*Dab.y*fa.x+rcd*Ddc.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].by+=rbc*Dcb.y*(fc.y+fd.y)+rab*Dab.y*fa.y+rcd*Ddc.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].bz+=rbc*Dcb.y*(fc.z+fd.z)+rab*Dab.y*fa.z+rcd*Ddc.y*fd.z;

    StrainDerivativeTensor[CurrentSystem].cx+=rbc*Dcb.z*(fc.x+fd.x)+rab*Dab.z*fa.x+rcd*Ddc.z*fd.x;
    StrainDerivativeTensor[CurrentSystem].cy+=rbc*Dcb.z*(fc.y+fd.y)+rab*Dab.z*fa.y+rcd*Ddc.z*fd.y;
    StrainDerivativeTensor[CurrentSystem].cz+=rbc*Dcb.z*(fc.z+fd.z)+rab*Dab.z*fa.z+rcd*Ddc.z*fd.z;
  }
}

void CalculateCationBondTorsionForce(int m)
{
  int i,A,B,C,D,Type;
  POINT posA,posB,posC,posD;
  REAL d,e, energy,rab,rbc,rcd,temp;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2,DCos;
  VECTOR dtA,dtB,dtC,dtD,fa,fb,fc,fd;
  REAL gamsa,gamsb,gamsc;
  REAL *parms;

  Type=Cations[CurrentSystem][m].Type;
  for(i=0;i<Components[Type].NumberOfBondTorsions;i++)
  {
    A=Components[Type].BondTorsions[i].A;
    B=Components[Type].BondTorsions[i].B;
    C=Components[Type].BondTorsions[i].C;
    D=Components[Type].BondTorsions[i].D;
    parms=(REAL*)&Components[Type].BondTorsionArguments[i];

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    posC=Cations[CurrentSystem][m].Atoms[C].Position;
    posD=Cations[CurrentSystem][m].Atoms[D].Position;

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;
    Dab=ApplyBoundaryCondition(Dab);
    rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));

    Dcb.x=posC.x-posB.x;
    Dcb.y=posC.y-posB.y;
    Dcb.z=posC.z-posB.z;
    Dcb=ApplyBoundaryCondition(Dcb);
    rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
    Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

    Ddc.x=posD.x-posC.x;
    Ddc.y=posD.y-posC.y;
    Ddc.z=posD.z-posC.z;
    Ddc=ApplyBoundaryCondition(Ddc);
    rcd=sqrt(SQR(Ddc.x)+SQR(Ddc.y)+SQR(Ddc.z));

    dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
    dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

    dr.x=Dab.x-dot_ab*Dcb.x;
    dr.y=Dab.y-dot_ab*Dcb.y;
    dr.z=Dab.z-dot_ab*Dcb.z;
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    dr.x/=r; dr.y/=r; dr.z/=r;

    ds.x=Ddc.x-dot_cd*Dcb.x;
    ds.y=Ddc.y-dot_cd*Dcb.y;
    ds.z=Ddc.z-dot_cd*Dcb.z;
    s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
    ds.x/=s; ds.y/=s; ds.z/=s;

    // compute Cos(Phi)
    // Phi is defined in protein convention Phi(trans)=Pi
    CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    switch(Components[Type].BondTorsionType[i])
    {
      case MM3_BOND_TORSION_CROSS:
        // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
        // =====================================================================================
        // p_0     [kcal/A mole]
        // p_1     [kcal/A mole]
        // p_2     [kcal/A mole]
        // p_3     [A]
        temp=(rbc-parms[3]);
        energy=parms[0]*temp*CosPhi+parms[1]*temp*(2.0*CosPhi2-1.0)+parms[2]*temp*(4.0*CosPhi2*CosPhi-3.0*CosPhi);
        DCos=parms[0]*temp+4.0*parms[1]*temp*CosPhi+parms[2]*temp*(12.0*CosPhi2-3.0);
        gamsa=0.0;
        gamsb=-(parms[0]*CosPhi+parms[1]*(2.0*CosPhi2-1.0)+parms[2]*(4.0*CosPhi2*CosPhi-3.0*CosPhi));
        gamsc=0.0;
        break;
      default:
        fprintf(stderr, "Undefined Bond-Torsion potential in routine 'CalculateCationBondTorsionForce' ('internal_force.c')\n");
        exit(0);
        break;
    }

    // energy
    UCationBondTorsion[CurrentSystem]+=energy;

    // Calculate the first derivative vectors.
    d=dot_ab/rbc;
    e=dot_cd/rbc;

    dtA.x=(ds.x-CosPhi*dr.x)/r;
    dtA.y=(ds.y-CosPhi*dr.y)/r;
    dtA.z=(ds.z-CosPhi*dr.z)/r;

    dtD.x=(dr.x-CosPhi*ds.x)/s;
    dtD.y=(dr.y-CosPhi*ds.y)/s;
    dtD.z=(dr.z-CosPhi*ds.z)/s;

    dtB.x=dtA.x*(d-1.0)+e*dtD.x;
    dtB.y=dtA.y*(d-1.0)+e*dtD.y;
    dtB.z=dtA.z*(d-1.0)+e*dtD.z;

    dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
    dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
    dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

    Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;
    Ddc.x/=rcd; Ddc.y/=rcd; Ddc.z/=rcd;

    fa.x=DCos*dtA.x-gamsa*Dab.x;
    fa.y=DCos*dtA.y-gamsa*Dab.y;
    fa.z=DCos*dtA.z-gamsa*Dab.z;

    fb.x=DCos*dtB.x+gamsb*Dcb.x+gamsa*Dab.x;
    fb.y=DCos*dtB.y+gamsb*Dcb.y+gamsa*Dab.y;
    fb.z=DCos*dtB.z+gamsb*Dcb.z+gamsa*Dab.z;

    fc.x=DCos*dtC.x-gamsb*Dcb.x+gamsc*Ddc.x;
    fc.y=DCos*dtC.y-gamsb*Dcb.y+gamsc*Ddc.y;
    fc.z=DCos*dtC.z-gamsb*Dcb.z+gamsc*Ddc.z;

    fd.x=DCos*dtD.x-gamsc*Ddc.x;
    fd.y=DCos*dtD.y-gamsc*Ddc.y;
    fd.z=DCos*dtD.z-gamsc*Ddc.z;

    // forces torsion
    Cations[CurrentSystem][m].Atoms[A].Force.x-=fa.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y-=fa.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z-=fa.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x-=fb.x;
    Cations[CurrentSystem][m].Atoms[B].Force.y-=fb.y;
    Cations[CurrentSystem][m].Atoms[B].Force.z-=fb.z;

    Cations[CurrentSystem][m].Atoms[C].Force.x-=fc.x;
    Cations[CurrentSystem][m].Atoms[C].Force.y-=fc.y;
    Cations[CurrentSystem][m].Atoms[C].Force.z-=fc.z;

    Cations[CurrentSystem][m].Atoms[D].Force.x-=fd.x;
    Cations[CurrentSystem][m].Atoms[D].Force.y-=fd.y;
    Cations[CurrentSystem][m].Atoms[D].Force.z-=fd.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax+=rbc*Dcb.x*(fc.x+fd.x)+rab*Dab.x*fa.x+rcd*Ddc.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].ay+=rbc*Dcb.x*(fc.y+fd.y)+rab*Dab.x*fa.y+rcd*Ddc.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].az+=rbc*Dcb.x*(fc.z+fd.z)+rab*Dab.x*fa.z+rcd*Ddc.x*fd.z;

    StrainDerivativeTensor[CurrentSystem].bx+=rbc*Dcb.y*(fc.x+fd.x)+rab*Dab.y*fa.x+rcd*Ddc.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].by+=rbc*Dcb.y*(fc.y+fd.y)+rab*Dab.y*fa.y+rcd*Ddc.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].bz+=rbc*Dcb.y*(fc.z+fd.z)+rab*Dab.y*fa.z+rcd*Ddc.y*fd.z;

    StrainDerivativeTensor[CurrentSystem].cx+=rbc*Dcb.z*(fc.x+fd.x)+rab*Dab.z*fa.x+rcd*Ddc.z*fd.x;
    StrainDerivativeTensor[CurrentSystem].cy+=rbc*Dcb.z*(fc.y+fd.y)+rab*Dab.z*fa.y+rcd*Ddc.z*fd.y;
    StrainDerivativeTensor[CurrentSystem].cz+=rbc*Dcb.z*(fc.z+fd.z)+rab*Dab.z*fa.z+rcd*Ddc.z*fd.z;
  }
}




// calculate intra LJ-force (status: tested and working)
void CalculateAdsorbateIntraVDWForce(int m)
{
  int i,NumberOfIntraVDW,Type,typeA,typeB,A,B;
  REAL rr,energy,force_factor,Scaling;
  VECTOR dr,f;
  POINT posA,posB;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfIntraVDW=Components[Type].NumberOfIntraVDW;
  for(i=0;i<NumberOfIntraVDW;i++)
  {
    A=Components[Type].IntraVDW[i].A;
    B=Components[Type].IntraVDW[i].B;
    Scaling=Components[Type].IntraVDWScaling[i];
    typeA=Adsorbates[CurrentSystem][m].Atoms[A].Type;
    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

    // note: a cutoff is also used for intra-van der Waals forces
    if(rr<CutOffVDWSquared)
    {
      typeB=Adsorbates[CurrentSystem][m].Atoms[B].Type;
      PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);
      energy*=Scaling;
      force_factor*=Scaling;

      // energy
      UAdsorbateIntraVDW[CurrentSystem]+=energy;

      StrainDerivativeTensor[CurrentSystem].ax+=force_factor*dr.x*dr.x;
      StrainDerivativeTensor[CurrentSystem].bx+=force_factor*dr.y*dr.x;
      StrainDerivativeTensor[CurrentSystem].cx+=force_factor*dr.z*dr.x;

      StrainDerivativeTensor[CurrentSystem].ay+=force_factor*dr.x*dr.y;
      StrainDerivativeTensor[CurrentSystem].by+=force_factor*dr.y*dr.y;
      StrainDerivativeTensor[CurrentSystem].cy+=force_factor*dr.z*dr.y;

      StrainDerivativeTensor[CurrentSystem].az+=force_factor*dr.x*dr.z;
      StrainDerivativeTensor[CurrentSystem].bz+=force_factor*dr.y*dr.z;
      StrainDerivativeTensor[CurrentSystem].cz+=force_factor*dr.z*dr.z;

      // forces
      f.x=force_factor*dr.x;
      f.y=force_factor*dr.y;
      f.z=force_factor*dr.z;

      Adsorbates[CurrentSystem][m].Atoms[A].Force.x-=f.x;
      Adsorbates[CurrentSystem][m].Atoms[A].Force.y-=f.y;
      Adsorbates[CurrentSystem][m].Atoms[A].Force.z-=f.z;

      Adsorbates[CurrentSystem][m].Atoms[B].Force.x+=f.x;
      Adsorbates[CurrentSystem][m].Atoms[B].Force.y+=f.y;
      Adsorbates[CurrentSystem][m].Atoms[B].Force.z+=f.z;
    }
  }
}

// calculate intra LJ-force (status: tested and working)
void CalculateCationIntraVDWForce(int m)
{
  int i,NumberOfIntraVDW,Type,typeA,typeB,A,B;
  REAL rr,energy,force_factor,Scaling;
  VECTOR dr,f;
  POINT posA,posB;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfIntraVDW=Components[Type].NumberOfIntraVDW;
  for(i=0;i<NumberOfIntraVDW;i++)
  {
    A=Components[Type].IntraVDW[i].A;
    B=Components[Type].IntraVDW[i].B;
    Scaling=Components[Type].IntraVDWScaling[i];
    typeA=Cations[CurrentSystem][m].Atoms[A].Type;
    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

    // note: a cutoff is also used for intra-van der Waals forces
    if(rr<CutOffVDWSquared)
    {
      typeB=Cations[CurrentSystem][m].Atoms[B].Type;
      PotentialGradient(typeA,typeB,rr,&energy,&force_factor,1.0);
      energy*=Scaling;
      force_factor*=Scaling;

      // energy
      UCationIntraVDW[CurrentSystem]+=energy;

      StrainDerivativeTensor[CurrentSystem].ax+=force_factor*dr.x*dr.x;
      StrainDerivativeTensor[CurrentSystem].bx+=force_factor*dr.y*dr.x;
      StrainDerivativeTensor[CurrentSystem].cx+=force_factor*dr.z*dr.x;

      StrainDerivativeTensor[CurrentSystem].ay+=force_factor*dr.x*dr.y;
      StrainDerivativeTensor[CurrentSystem].by+=force_factor*dr.y*dr.y;
      StrainDerivativeTensor[CurrentSystem].cy+=force_factor*dr.z*dr.y;

      StrainDerivativeTensor[CurrentSystem].az+=force_factor*dr.x*dr.z;
      StrainDerivativeTensor[CurrentSystem].bz+=force_factor*dr.y*dr.z;
      StrainDerivativeTensor[CurrentSystem].cz+=force_factor*dr.z*dr.z;

      // forces
      f.x=force_factor*dr.x;
      f.y=force_factor*dr.y;
      f.z=force_factor*dr.z;

      Cations[CurrentSystem][m].Atoms[A].Force.x-=f.x;
      Cations[CurrentSystem][m].Atoms[A].Force.y-=f.y;
      Cations[CurrentSystem][m].Atoms[A].Force.z-=f.z;

      Cations[CurrentSystem][m].Atoms[B].Force.x+=f.x;
      Cations[CurrentSystem][m].Atoms[B].Force.y+=f.y;
      Cations[CurrentSystem][m].Atoms[B].Force.z+=f.z;
    }
  }
}

// calculate intra Coulomb-force
void CalculateAdsorbateIntraChargeChargeForce(int m)
{
  int i,NumberOfIntraChargeCharge,Type,typeA,typeB,A,B;
  REAL r,rr;
  VECTOR dr;
  POINT posA,posB;
  REAL chargeA,chargeB,temp,Scaling;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfIntraChargeCharge=Components[Type].NumberOfIntraChargeCharge;
  for(i=0;i<NumberOfIntraChargeCharge;i++)
  {
    A=Components[Type].IntraChargeCharge[i].A;
    B=Components[Type].IntraChargeCharge[i].B;
    Scaling=Components[Type].IntraChargeChargeScaling[i];

    typeA=Adsorbates[CurrentSystem][m].Atoms[A].Type;
    chargeA=Adsorbates[CurrentSystem][m].Atoms[A].Charge;
    typeB=Adsorbates[CurrentSystem][m].Atoms[B].Type;
    chargeB=Adsorbates[CurrentSystem][m].Atoms[B].Charge;

    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
    posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    // note: no cutoff used here
    switch(ChargeMethod)
    {
      case NONE:
        temp=0.0;
        break;
      case SHIFTED_COULOMB:
      case TRUNCATED_COULOMB:
      case EWALD:
      default:
        UAdsorbateIntraChargeCharge[CurrentSystem]+=Scaling*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
        temp=Scaling*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
        break;
    }

    StrainDerivativeTensor[CurrentSystem].ax-=temp*dr.x*dr.x;
    StrainDerivativeTensor[CurrentSystem].bx-=temp*dr.y*dr.x;
    StrainDerivativeTensor[CurrentSystem].cx-=temp*dr.z*dr.x;

    StrainDerivativeTensor[CurrentSystem].ay-=temp*dr.x*dr.y;
    StrainDerivativeTensor[CurrentSystem].by-=temp*dr.y*dr.y;
    StrainDerivativeTensor[CurrentSystem].cy-=temp*dr.z*dr.y;

    StrainDerivativeTensor[CurrentSystem].az-=temp*dr.x*dr.z;
    StrainDerivativeTensor[CurrentSystem].bz-=temp*dr.y*dr.z;
    StrainDerivativeTensor[CurrentSystem].cz-=temp*dr.z*dr.z;

    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=temp*dr.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=temp*dr.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=temp*dr.z;

    Adsorbates[CurrentSystem][m].Atoms[B].Force.x-=temp*dr.x;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.y-=temp*dr.y;
    Adsorbates[CurrentSystem][m].Atoms[B].Force.z-=temp*dr.z;
  }
}

// calculate intra Coulomb-force
void CalculateCationIntraChargeChargeForce(int m)
{
  int i,NumberOfIntraChargeCharge,Type,typeA,typeB,A,B;
  REAL r,rr;
  VECTOR dr;
  POINT posA,posB;
  REAL chargeA,chargeB,temp,Scaling;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfIntraChargeCharge=Components[Type].NumberOfIntraChargeCharge;
  for(i=0;i<NumberOfIntraChargeCharge;i++)
  {
    A=Components[Type].IntraChargeCharge[i].A;
    B=Components[Type].IntraChargeCharge[i].B;
    Scaling=Components[Type].IntraChargeChargeScaling[i];

    typeA=Cations[CurrentSystem][m].Atoms[A].Type;
    chargeA=Cations[CurrentSystem][m].Atoms[A].Charge;
    typeB=Cations[CurrentSystem][m].Atoms[B].Type;
    chargeB=Cations[CurrentSystem][m].Atoms[B].Charge;

    posA=Cations[CurrentSystem][m].Atoms[A].Position;
    posB=Cations[CurrentSystem][m].Atoms[B].Position;
    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    // note: no cutoff used here
    switch(ChargeMethod)
    {
      case NONE:
        temp=0.0;
        break;
      case SHIFTED_COULOMB:
      case TRUNCATED_COULOMB:
      case EWALD:
      default:
        UCationIntraChargeCharge[CurrentSystem]+=Scaling*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
        temp=Scaling*COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
        break;
    }

    StrainDerivativeTensor[CurrentSystem].ax-=temp*dr.x*dr.x;
    StrainDerivativeTensor[CurrentSystem].bx-=temp*dr.y*dr.x;
    StrainDerivativeTensor[CurrentSystem].cx-=temp*dr.z*dr.x;

    StrainDerivativeTensor[CurrentSystem].ay-=temp*dr.x*dr.y;
    StrainDerivativeTensor[CurrentSystem].by-=temp*dr.y*dr.y;
    StrainDerivativeTensor[CurrentSystem].cy-=temp*dr.z*dr.y;

    StrainDerivativeTensor[CurrentSystem].az-=temp*dr.x*dr.z;
    StrainDerivativeTensor[CurrentSystem].bz-=temp*dr.y*dr.z;
    StrainDerivativeTensor[CurrentSystem].cz-=temp*dr.z*dr.z;

    Cations[CurrentSystem][m].Atoms[A].Force.x+=temp*dr.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=temp*dr.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=temp*dr.z;

    Cations[CurrentSystem][m].Atoms[B].Force.x-=temp*dr.x;
    Cations[CurrentSystem][m].Atoms[B].Force.y-=temp*dr.y;
    Cations[CurrentSystem][m].Atoms[B].Force.z-=temp*dr.z;
  }
}

// calculate intra Coulomb-force
void CalculateAdsorbateIntraChargeBondDipoleForce(int m)
{
  int i,NumberOfIntraChargeBondDipole,Type,A,B,B1,B2;
  REAL r,rr,ri2,length;
  REAL Bt0,Bt1,Bt2,cosB;
  VECTOR dr,dipoleB,term,fa1,fb1,fb2;
  POINT posA,posB,posB1,posB2;
  REAL ChargeA,temp;
  REAL DipoleMagnitudeB,energy;
  REAL fac1,fac2;
  REAL_MATRIX3x3 v;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfIntraChargeBondDipole=Components[Type].NumberOfIntraChargeBondDipole;
  for(i=0;i<NumberOfIntraChargeBondDipole;i++)
  {
    A=Components[Type].IntraChargeBondDipole[i].A;
    B=Components[Type].IntraChargeBondDipole[i].B;

    ChargeA=Adsorbates[CurrentSystem][m].Atoms[A].Charge;
    posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;

    B1=Components[Type].BondDipoles[B].A;
    B2=Components[Type].BondDipoles[B].B;
    DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
    posB1=Adsorbates[CurrentSystem][m].Atoms[B1].Position;
    posB2=Adsorbates[CurrentSystem][m].Atoms[B2].Position;
    dipoleB.x=posB2.x-posB1.x;
    dipoleB.y=posB2.y-posB1.y;
    dipoleB.z=posB2.z-posB1.z;
    posB.x=posB1.x+0.5*dipoleB.x;
    posB.y=posB1.y+0.5*dipoleB.y;
    posB.z=posB1.z+0.5*dipoleB.z;
    ri2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
    length=sqrt(ri2);
    temp=DipoleMagnitudeB/length;
    dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

    dr.x=posB.x-posA.x;
    dr.y=posB.y-posA.y;
    dr.z=posB.z-posA.z;
    dr=ApplyBoundaryCondition(dr);
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    Bt0=1.0/(r);
    Bt1=1.0/(r*rr);
    Bt2=3.0/(r*rr*rr);
    cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
    energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
    UAdsorbateIntraChargeBondDipole[CurrentSystem]-=energy;

    term.x=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
    term.y=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
    term.z=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

    fa1.x=term.x;
    fa1.y=term.y;
    fa1.z=term.z;

    Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa1.x;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa1.y;
    Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa1.z;

    fac1=COULOMBIC_CONVERSION_FACTOR*ChargeA*Bt1*DipoleMagnitudeB/length;
    fac2=COULOMBIC_CONVERSION_FACTOR*ChargeA*Bt1*cosB/(DipoleMagnitudeB*length);

    fb1.x=0.5*term.x+fac1*dr.x-fac2*dipoleB.x;
    fb1.y=0.5*term.y+fac1*dr.y-fac2*dipoleB.y;
    fb1.z=0.5*term.z+fac1*dr.z-fac2*dipoleB.z;

    fb2.x=0.5*term.x-fac1*dr.x+fac2*dipoleB.x;
    fb2.y=0.5*term.y-fac1*dr.y+fac2*dipoleB.y;
    fb2.z=0.5*term.z-fac1*dr.z+fac2*dipoleB.z;

    Adsorbates[CurrentSystem][m].Atoms[B1].Force.x-=fb1.x;
    Adsorbates[CurrentSystem][m].Atoms[B1].Force.y-=fb1.y;
    Adsorbates[CurrentSystem][m].Atoms[B1].Force.z-=fb1.z;

    Adsorbates[CurrentSystem][m].Atoms[B2].Force.x-=fb2.x;
    Adsorbates[CurrentSystem][m].Atoms[B2].Force.y-=fb2.y;
    Adsorbates[CurrentSystem][m].Atoms[B2].Force.z-=fb2.z;

    // convert forces on atoms to molecular virial
    v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
    v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
    v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

    v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
    v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
    v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

    v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
    v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
    v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

    // the strain derivative
    StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
    StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
    StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

    StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
    StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
    StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

    StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
    StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
    StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
  }
}

// calculate intra Coulomb-force
void CalculateCationIntraChargeBondDipoleForce(int m)
{
  int i,NumberOfIntraChargeBondDipole,Type,A,B,B1,B2;
  REAL r,rr,ri2,length;
  REAL Bt0,Bt1,Bt2,cosB;
  VECTOR dr,dipoleB,term,fa1,fb1,fb2;
  POINT posA,posB,posB1,posB2;
  REAL ChargeA,temp;
  REAL DipoleMagnitudeB,energy;
  REAL fac1,fac2;
  REAL_MATRIX3x3 v;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfIntraChargeBondDipole=Components[Type].NumberOfIntraChargeBondDipole;
  for(i=0;i<NumberOfIntraChargeBondDipole;i++)
  {
    A=Components[Type].IntraChargeBondDipole[i].A;
    B=Components[Type].IntraChargeBondDipole[i].B;

    ChargeA=Cations[CurrentSystem][m].Atoms[A].Charge;
    posA=Cations[CurrentSystem][m].Atoms[A].Position;

    B1=Components[Type].BondDipoles[B].A;
    B2=Components[Type].BondDipoles[B].B;
    DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
    posB1=Cations[CurrentSystem][m].Atoms[B1].Position;
    posB2=Cations[CurrentSystem][m].Atoms[B2].Position;
    dipoleB.x=posB2.x-posB1.x;
    dipoleB.y=posB2.y-posB1.y;
    dipoleB.z=posB2.z-posB1.z;
    posB.x=posB1.x+0.5*dipoleB.x;
    posB.y=posB1.y+0.5*dipoleB.y;
    posB.z=posB1.z+0.5*dipoleB.z;
    ri2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
    length=sqrt(ri2);
    temp=DipoleMagnitudeB/length;
    dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

    dr.x=posB.x-posA.x;
    dr.y=posB.y-posA.y;
    dr.z=posB.z-posA.z;
    dr=ApplyBoundaryCondition(dr);
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    Bt0=1.0/(r);
    Bt1=1.0/(r*rr);
    Bt2=3.0/(r*rr*rr);
    cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
    energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
    UCationIntraChargeBondDipole[CurrentSystem]-=energy;

    term.x=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
    term.y=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
    term.z=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

    fa1.x=term.x;
    fa1.y=term.y;
    fa1.z=term.z;

    Cations[CurrentSystem][m].Atoms[A].Force.x+=fa1.x;
    Cations[CurrentSystem][m].Atoms[A].Force.y+=fa1.y;
    Cations[CurrentSystem][m].Atoms[A].Force.z+=fa1.z;

    fac1=COULOMBIC_CONVERSION_FACTOR*ChargeA*Bt1*DipoleMagnitudeB/length;
    fac2=COULOMBIC_CONVERSION_FACTOR*ChargeA*Bt1*cosB/(DipoleMagnitudeB*length);

    fb1.x=0.5*term.x+fac1*dr.x-fac2*dipoleB.x;
    fb1.y=0.5*term.y+fac1*dr.y-fac2*dipoleB.y;
    fb1.z=0.5*term.z+fac1*dr.z-fac2*dipoleB.z;

    fb2.x=0.5*term.x-fac1*dr.x+fac2*dipoleB.x;
    fb2.y=0.5*term.y-fac1*dr.y+fac2*dipoleB.y;
    fb2.z=0.5*term.z-fac1*dr.z+fac2*dipoleB.z;

    Cations[CurrentSystem][m].Atoms[B1].Force.x-=fb1.x;
    Cations[CurrentSystem][m].Atoms[B1].Force.y-=fb1.y;
    Cations[CurrentSystem][m].Atoms[B1].Force.z-=fb1.z;

    Cations[CurrentSystem][m].Atoms[B2].Force.x-=fb2.x;
    Cations[CurrentSystem][m].Atoms[B2].Force.y-=fb2.y;
    Cations[CurrentSystem][m].Atoms[B2].Force.z-=fb2.z;

    // convert forces on atoms to molecular virial
    v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
    v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
    v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

    v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
    v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
    v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

    v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
    v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
    v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

    // the strain derivative
    StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
    StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
    StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

    StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
    StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
    StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

    StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
    StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
    StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
  }
}

// calculate intra Coulomb-force
void CalculateAdsorbateIntraBondDipoleBondDipoleForce(int m)
{
  int i,NumberOfIntraBondDipoleBondDipole,Type,A,B,A1,A2,B1,B2;
  REAL r,rr,ri2,rk2,length;
  REAL Bt0,Bt1,Bt2,Bt3,cosA,cosB,cosAB;
  VECTOR dr,dipoleA,dipoleB,term,termA,termB,fa1,fa2,fb1,fb2;
  POINT posA,posB,posA1,posA2,posB1,posB2;
  REAL temp;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,energy;
  REAL_MATRIX3x3 v;

  Type=Adsorbates[CurrentSystem][m].Type;
  NumberOfIntraBondDipoleBondDipole=Components[Type].NumberOfIntraBondDipoleBondDipole;
  for(i=0;i<NumberOfIntraBondDipoleBondDipole;i++)
  {
    A=Components[Type].IntraBondDipoleBondDipole[i].A;
    B=Components[Type].IntraBondDipoleBondDipole[i].B;

    A1=Components[Type].BondDipoles[A].A;
    A2=Components[Type].BondDipoles[A].B;
    posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
    posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
    DipoleMagnitudeA=Components[Type].BondDipoleMagnitude[A];
    dipoleA.x=posA2.x-posA1.x;
    dipoleA.y=posA2.y-posA1.y;
    dipoleA.z=posA2.z-posA1.z;
    posA.x=posA1.x+0.5*dipoleA.x;
    posA.y=posA1.y+0.5*dipoleA.y;
    posA.z=posA1.z+0.5*dipoleA.z;
    ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
    length=sqrt(ri2);
    temp=DipoleMagnitudeA/length;
    dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

    B1=Components[Type].BondDipoles[B].A;
    B2=Components[Type].BondDipoles[B].B;
    posB1=Adsorbates[CurrentSystem][m].Atoms[B1].Position;
    posB2=Adsorbates[CurrentSystem][m].Atoms[B2].Position;
    DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
    dipoleB.x=posB2.x-posB1.x;
    dipoleB.y=posB2.y-posB1.y;
    dipoleB.z=posB2.z-posB1.z;
    posB.x=posB1.x+0.5*dipoleB.x;
    posB.y=posB1.y+0.5*dipoleB.y;
    posB.z=posB1.z+0.5*dipoleB.z;
    rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
    length=sqrt(rk2);
    temp=DipoleMagnitudeB/length;
    dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    dr=ApplyBoundaryCondition(dr);
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    Bt0=1.0/(r);
    Bt1=1.0/(r*rr);
    Bt2=3.0/(r*rr*rr);
    Bt3=15.0/(r*rr*rr*rr);

    cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
    cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
    cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
    energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
    UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=energy;

    term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
    term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
    term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

    termA.x=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt2*cosA*cosB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt1*DipoleMagnitudeA*dipoleB.x/sqrt(ri2)
            +Bt2*DipoleMagnitudeA*cosB*dr.x/sqrt(ri2));

    termA.y=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt2*cosA*cosB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt1*DipoleMagnitudeA*dipoleB.y/sqrt(ri2)
            +Bt2*DipoleMagnitudeA*cosB*dr.y/sqrt(ri2));

    termA.z=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt2*cosA*cosB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt1*DipoleMagnitudeA*dipoleB.z/sqrt(ri2)
            +Bt2*DipoleMagnitudeA*cosB*dr.z/sqrt(ri2));

    termB.x=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt2*cosA*cosB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt1*DipoleMagnitudeB*dipoleA.x/sqrt(rk2)
            +Bt2*DipoleMagnitudeB*cosA*dr.x/sqrt(rk2));

    termB.y=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt2*cosA*cosB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt1*DipoleMagnitudeB*dipoleA.y/sqrt(rk2)
            +Bt2*DipoleMagnitudeB*cosA*dr.y/sqrt(rk2));

    termB.z=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt2*cosA*cosB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt1*DipoleMagnitudeB*dipoleA.z/sqrt(rk2)
            +Bt2*DipoleMagnitudeB*cosA*dr.z/sqrt(rk2));

    fa1.x=0.5*term.x+termA.x;
    fa1.y=0.5*term.y+termA.y;
    fa1.z=0.5*term.z+termA.z;
    fa2.x=0.5*term.x-termA.x;
    fa2.y=0.5*term.y-termA.y;
    fa2.z=0.5*term.z-termA.z;

    fb1.x=-0.5*term.x+termB.x;
    fb1.y=-0.5*term.y+termB.y;
    fb1.z=-0.5*term.z+termB.z;
    fb2.x=-0.5*term.x-termB.x;
    fb2.y=-0.5*term.y-termB.y;
    fb2.z=-0.5*term.z-termB.z;

    Adsorbates[CurrentSystem][m].Atoms[A1].Force.x-=fa1.x;
    Adsorbates[CurrentSystem][m].Atoms[A1].Force.y-=fa1.y;
    Adsorbates[CurrentSystem][m].Atoms[A1].Force.z-=fa1.z;

    Adsorbates[CurrentSystem][m].Atoms[A2].Force.x-=fa2.x;
    Adsorbates[CurrentSystem][m].Atoms[A2].Force.y-=fa2.y;
    Adsorbates[CurrentSystem][m].Atoms[A2].Force.z-=fa2.z;

    Adsorbates[CurrentSystem][m].Atoms[B1].Force.x-=fb1.x;
    Adsorbates[CurrentSystem][m].Atoms[B1].Force.y-=fb1.y;
    Adsorbates[CurrentSystem][m].Atoms[B1].Force.z-=fb1.z;

    Adsorbates[CurrentSystem][m].Atoms[B2].Force.x-=fb2.x;
    Adsorbates[CurrentSystem][m].Atoms[B2].Force.y-=fb2.y;
    Adsorbates[CurrentSystem][m].Atoms[B2].Force.z-=fb2.z;

    // convert forces on atoms to molecular virial
    // usually this conversion produces a torque on the center of mass and an asymmetric stress
    // by making it symmetric we regain the 'atomic' strain derivative
    v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
    v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
    v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

    v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
    v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
    v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

    v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
    v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
    v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

    // the strain derivative
    StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
    StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
    StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

    StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
    StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
    StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

    StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
    StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
    StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
  }
}

// calculate intra Coulomb-force
void CalculateCationIntraBondDipoleBondDipoleForce(int m)
{
  int i,NumberOfIntraBondDipoleBondDipole,Type,A,B,A1,A2,B1,B2;
  REAL r,rr,ri2,rk2,length;
  REAL Bt0,Bt1,Bt2,Bt3,cosA,cosB,cosAB;
  VECTOR dr,dipoleA,dipoleB,term,termA,termB,fa1,fa2,fb1,fb2;
  POINT posA,posB,posA1,posA2,posB1,posB2;
  REAL temp;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,energy;
  REAL_MATRIX3x3 v;

  Type=Cations[CurrentSystem][m].Type;
  NumberOfIntraBondDipoleBondDipole=Components[Type].NumberOfIntraBondDipoleBondDipole;
  for(i=0;i<NumberOfIntraBondDipoleBondDipole;i++)
  {
    A=Components[Type].IntraBondDipoleBondDipole[i].A;
    B=Components[Type].IntraBondDipoleBondDipole[i].B;

    A1=Components[Type].BondDipoles[A].A;
    A2=Components[Type].BondDipoles[A].B;
    posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
    posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
    DipoleMagnitudeA=Components[Type].BondDipoleMagnitude[A];
    dipoleA.x=posA2.x-posA1.x;
    dipoleA.y=posA2.y-posA1.y;
    dipoleA.z=posA2.z-posA1.z;
    posA.x=posA1.x+0.5*dipoleA.x;
    posA.y=posA1.y+0.5*dipoleA.y;
    posA.z=posA1.z+0.5*dipoleA.z;
    ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
    length=sqrt(ri2);
    temp=DipoleMagnitudeA/length;
    dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

    B1=Components[Type].BondDipoles[B].A;
    B2=Components[Type].BondDipoles[B].B;
    posB1=Cations[CurrentSystem][m].Atoms[B1].Position;
    posB2=Cations[CurrentSystem][m].Atoms[B2].Position;
    DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
    dipoleB.x=posB2.x-posB1.x;
    dipoleB.y=posB2.y-posB1.y;
    dipoleB.z=posB2.z-posB1.z;
    posB.x=posB1.x+0.5*dipoleB.x;
    posB.y=posB1.y+0.5*dipoleB.y;
    posB.z=posB1.z+0.5*dipoleB.z;
    rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
    length=sqrt(rk2);
    temp=DipoleMagnitudeB/length;
    dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    dr=ApplyBoundaryCondition(dr);
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    Bt0=1.0/(r);
    Bt1=1.0/(r*rr);
    Bt2=3.0/(r*rr*rr);
    Bt3=15.0/(r*rr*rr*rr);

    cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
    cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
    cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
    energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
    UCationIntraBondDipoleBondDipole[CurrentSystem]+=energy;

    term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
    term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
    term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

    termA.x=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt2*cosA*cosB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt1*DipoleMagnitudeA*dipoleB.x/sqrt(ri2)
            +Bt2*DipoleMagnitudeA*cosB*dr.x/sqrt(ri2));

    termA.y=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt2*cosA*cosB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt1*DipoleMagnitudeA*dipoleB.y/sqrt(ri2)
            +Bt2*DipoleMagnitudeA*cosB*dr.y/sqrt(ri2));

    termA.z=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt2*cosA*cosB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
            -Bt1*DipoleMagnitudeA*dipoleB.z/sqrt(ri2)
            +Bt2*DipoleMagnitudeA*cosB*dr.z/sqrt(ri2));

    termB.x=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt2*cosA*cosB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt1*DipoleMagnitudeB*dipoleA.x/sqrt(rk2)
            +Bt2*DipoleMagnitudeB*cosA*dr.x/sqrt(rk2));

    termB.y=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt2*cosA*cosB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt1*DipoleMagnitudeB*dipoleA.y/sqrt(rk2)
            +Bt2*DipoleMagnitudeB*cosA*dr.y/sqrt(rk2));

    termB.z=COULOMBIC_CONVERSION_FACTOR*(
            Bt1*cosAB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt2*cosA*cosB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
            -Bt1*DipoleMagnitudeB*dipoleA.z/sqrt(rk2)
            +Bt2*DipoleMagnitudeB*cosA*dr.z/sqrt(rk2));

    fa1.x=0.5*term.x+termA.x;
    fa1.y=0.5*term.y+termA.y;
    fa1.z=0.5*term.z+termA.z;
    fa2.x=0.5*term.x-termA.x;
    fa2.y=0.5*term.y-termA.y;
    fa2.z=0.5*term.z-termA.z;

    fb1.x=-0.5*term.x+termB.x;
    fb1.y=-0.5*term.y+termB.y;
    fb1.z=-0.5*term.z+termB.z;
    fb2.x=-0.5*term.x-termB.x;
    fb2.y=-0.5*term.y-termB.y;
    fb2.z=-0.5*term.z-termB.z;

    Cations[CurrentSystem][m].Atoms[A1].Force.x-=fa1.x;
    Cations[CurrentSystem][m].Atoms[A1].Force.y-=fa1.y;
    Cations[CurrentSystem][m].Atoms[A1].Force.z-=fa1.z;

    Cations[CurrentSystem][m].Atoms[A2].Force.x-=fa2.x;
    Cations[CurrentSystem][m].Atoms[A2].Force.y-=fa2.y;
    Cations[CurrentSystem][m].Atoms[A2].Force.z-=fa2.z;

    Cations[CurrentSystem][m].Atoms[B1].Force.x-=fb1.x;
    Cations[CurrentSystem][m].Atoms[B1].Force.y-=fb1.y;
    Cations[CurrentSystem][m].Atoms[B1].Force.z-=fb1.z;

    Cations[CurrentSystem][m].Atoms[B2].Force.x-=fb2.x;
    Cations[CurrentSystem][m].Atoms[B2].Force.y-=fb2.y;
    Cations[CurrentSystem][m].Atoms[B2].Force.z-=fb2.z;

    // convert forces on atoms to molecular virial
    // usually this conversion produces a torque on the center of mass and an asymmetric stress
    // by making it symmetric we regain the 'atomic' strain derivative
    v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
    v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
    v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

    v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
    v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
    v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

    v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
    v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
    v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

    // the strain derivative
    StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
    StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
    StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

    StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
    StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
    StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

    StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
    StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
    StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
  }
}

void CalculateHarmonicBondConstraintForce(void)
{
  int m;
  REAL U,DF,r,rr;
  REAL parms0,parms1;
  POINT posA,posB;
  VECTOR dr,f;

  UDistanceConstraints[CurrentSystem]=0.0;
  for(m=0;m<NumberOfHarmonicDistanceConstraints[CurrentSystem];m++)
  {
    posA=HarmonicDistanceConstraints[CurrentSystem][m][0]->Position;
    posB=HarmonicDistanceConstraints[CurrentSystem][m][1]->Position;

    parms0=HarmonicDistanceConstraintParameters[CurrentSystem][m][0];
    parms1=HarmonicDistanceConstraintParameters[CurrentSystem][m][1];

    dr.x=posA.x-posB.x;
    dr.y=posA.y-posB.y;
    dr.z=posA.z-posB.z;
    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    r=sqrt(rr);

    U=0.5*parms0*SQR(r-parms1);
    DF=parms0*(r-parms1)/r;

    // add contribution to the energy
    UDistanceConstraints[CurrentSystem]+=U;

    // forces are oppositely directed to the gradient
    f.x=-DF*dr.x;
    f.y=-DF*dr.y;
    f.z=-DF*dr.z;

    // add contribution to the forces
    HarmonicDistanceConstraints[CurrentSystem][m][0]->Force.x+=f.x;
    HarmonicDistanceConstraints[CurrentSystem][m][0]->Force.y+=f.y;
    HarmonicDistanceConstraints[CurrentSystem][m][0]->Force.z+=f.z;

    HarmonicDistanceConstraints[CurrentSystem][m][1]->Force.x-=f.x;
    HarmonicDistanceConstraints[CurrentSystem][m][1]->Force.y-=f.y;
    HarmonicDistanceConstraints[CurrentSystem][m][1]->Force.z-=f.z;

    // add contribution to the stress tensor
    StrainDerivativeTensor[CurrentSystem].ax-=dr.x*f.x;
    StrainDerivativeTensor[CurrentSystem].bx-=dr.y*f.x;
    StrainDerivativeTensor[CurrentSystem].cx-=dr.z*f.x;

    StrainDerivativeTensor[CurrentSystem].ay-=dr.x*f.y;
    StrainDerivativeTensor[CurrentSystem].by-=dr.y*f.y;
    StrainDerivativeTensor[CurrentSystem].cy-=dr.z*f.y;

    StrainDerivativeTensor[CurrentSystem].az-=dr.x*f.z;
    StrainDerivativeTensor[CurrentSystem].bz-=dr.y*f.z;
    StrainDerivativeTensor[CurrentSystem].cz-=dr.z*f.z;
  }
}

void CalculateHarmonicAngleConstraintForce(void)
{
  int m;
  REAL DF,U;
  REAL CosTheta,Theta,SinTheta;
  REAL rab,rac,rbc,DTDX;
  POINT posA,posB,posC;
  VECTOR Rab,Rac,Rbc,fa,fb,fc,dtA,dtB,dtC;
  REAL parms0,parms1;

  UAngleConstraints[CurrentSystem]=0.0;
  for(m=0;m<NumberOfHarmonicAngleConstraints[CurrentSystem];m++)
  {
    posA=HarmonicAngleConstraints[CurrentSystem][m][0]->Position;
    posB=HarmonicAngleConstraints[CurrentSystem][m][1]->Position;
    posC=HarmonicAngleConstraints[CurrentSystem][m][2]->Position;

    parms0=HarmonicAngleConstraintParameters[CurrentSystem][m][0];
    parms1=HarmonicAngleConstraintParameters[CurrentSystem][m][1];

    Rab.x=posA.x-posB.x;
    Rab.y=posA.y-posB.y;
    Rab.z=posA.z-posB.z;
    rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
    Rab.x/=rab;
    Rab.y/=rab;
    Rab.z/=rab;

    Rbc.x=posC.x-posB.x;
    Rbc.y=posC.y-posB.y;
    Rbc.z=posC.z-posB.z;
    rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
    Rbc.x/=rbc;
    Rbc.y/=rbc;
    Rbc.z/=rbc;

    Rac.x=posC.x-posA.x;
    Rac.y=posC.y-posA.y;
    Rac.z=posC.z-posA.z;
    rac=sqrt(SQR(Rac.x)+SQR(Rac.y)+SQR(Rac.z));
    Rac.x/=rac;
    Rac.y/=rac;
    Rac.z/=rac;

    CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
    Theta=acos(CosTheta);
    SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));
    DTDX=-1.0/sqrt(1.0-SQR(CosTheta));

    U=0.5*parms0*SQR(Theta-parms1);
    DF=parms0*(Theta-parms1)*DTDX;

    // add contribution to the energy
    UAngleConstraints[CurrentSystem]+=U;

    // Calculate the components of the derivatives.
    dtA.x=(Rbc.x-CosTheta*Rab.x)/rab;
    dtA.y=(Rbc.y-CosTheta*Rab.y)/rab;
    dtA.z=(Rbc.z-CosTheta*Rab.z)/rab;

    dtC.x=(Rab.x-CosTheta*Rbc.x)/rbc;
    dtC.y=(Rab.y-CosTheta*Rbc.y)/rbc;
    dtC.z=(Rab.z-CosTheta*Rbc.z)/rbc;

    dtB.x=-(dtA.x+dtC.x);
    dtB.y=-(dtA.y+dtC.y);
    dtB.z=-(dtA.z+dtC.z);

    // forces are oppositely directed to the gradient
    fa.x=-DF*dtA.x;
    fa.y=-DF*dtA.y;
    fa.z=-DF*dtA.z;

    fb.x=-DF*dtB.x;
    fb.y=-DF*dtB.y;
    fb.z=-DF*dtB.z;

    fc.x=-DF*dtC.x;
    fc.y=-DF*dtC.y;
    fc.z=-DF*dtC.z;

    // add contribution to the forces
    HarmonicAngleConstraints[CurrentSystem][m][0]->Force.x+=fa.x;
    HarmonicAngleConstraints[CurrentSystem][m][0]->Force.y+=fa.y;
    HarmonicAngleConstraints[CurrentSystem][m][0]->Force.z+=fa.z;

    HarmonicAngleConstraints[CurrentSystem][m][1]->Force.x-=(fa.x+fc.x);
    HarmonicAngleConstraints[CurrentSystem][m][1]->Force.y-=(fa.y+fc.y);
    HarmonicAngleConstraints[CurrentSystem][m][1]->Force.z-=(fa.z+fc.z);

    HarmonicAngleConstraints[CurrentSystem][m][2]->Force.x+=fc.x;
    HarmonicAngleConstraints[CurrentSystem][m][2]->Force.y+=fc.y;
    HarmonicAngleConstraints[CurrentSystem][m][2]->Force.z+=fc.z;

    // add contribution to the stress tensor
    // Note: rab and rbc are here because the vectors were normalized before
    StrainDerivativeTensor[CurrentSystem].ax-=rab*Rab.x*fa.x+rbc*Rbc.x*fc.x;
    StrainDerivativeTensor[CurrentSystem].bx-=rab*Rab.y*fa.x+rbc*Rbc.y*fc.x;
    StrainDerivativeTensor[CurrentSystem].cx-=rab*Rab.z*fa.x+rbc*Rbc.z*fc.x;

    StrainDerivativeTensor[CurrentSystem].ay-=rab*Rab.x*fa.y+rbc*Rbc.x*fc.y;
    StrainDerivativeTensor[CurrentSystem].by-=rab*Rab.y*fa.y+rbc*Rbc.y*fc.y;
    StrainDerivativeTensor[CurrentSystem].cy-=rab*Rab.z*fa.y+rbc*Rbc.z*fc.y;

    StrainDerivativeTensor[CurrentSystem].az-=rab*Rab.x*fa.z+rbc*Rbc.x*fc.z;
    StrainDerivativeTensor[CurrentSystem].bz-=rab*Rab.y*fa.z+rbc*Rbc.y*fc.z;
    StrainDerivativeTensor[CurrentSystem].cz-=rab*Rab.z*fa.z+rbc*Rbc.z*fc.z;
  }
}

void CalculateHarmonicDihedralConstraintForce(void)
{
  int m;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi,U,DF;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  VECTOR fa,fb,fc,fd;
  REAL parms0,parms1;

  UDihedralConstraints[CurrentSystem]=0.0;
  for(m=0;m<NumberOfHarmonicDihedralConstraints[CurrentSystem];m++)
  {
    posA=HarmonicDihedralConstraints[CurrentSystem][m][0]->Position;
    posB=HarmonicDihedralConstraints[CurrentSystem][m][1]->Position;
    posC=HarmonicDihedralConstraints[CurrentSystem][m][2]->Position;
    posD=HarmonicDihedralConstraints[CurrentSystem][m][3]->Position;

    parms0=HarmonicDihedralConstraintParameters[CurrentSystem][m][0];
    parms1=HarmonicDihedralConstraintParameters[CurrentSystem][m][1];

    Dab.x=posA.x-posB.x;
    Dab.y=posA.y-posB.y;
    Dab.z=posA.z-posB.z;

    Dcb.x=posC.x-posB.x;
    Dcb.y=posC.y-posB.y;
    Dcb.z=posC.z-posB.z;
    rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
    Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

    Ddc.x=posD.x-posC.x;
    Ddc.y=posD.y-posC.y;
    Ddc.z=posD.z-posC.z;
    dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
    dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

    dr.x=Dab.x-dot_ab*Dcb.x;
    dr.y=Dab.y-dot_ab*Dcb.y;
    dr.z=Dab.z-dot_ab*Dcb.z;
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    dr.x/=r; dr.y/=r; dr.z/=r;

    ds.x=Ddc.x-dot_cd*Dcb.x;
    ds.y=Ddc.y-dot_cd*Dcb.y;
    ds.z=Ddc.z-dot_cd*Dcb.z;
    s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
    ds.x/=s; ds.y/=s; ds.z/=s;

    // compute Cos(Phi)
    // Phi is defined in protein convention Phi(trans)=Pi
    CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

    // Ensure CosPhi is between -1 and 1.
    CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
    CosPhi2=SQR(CosPhi);

    // potential defined in terms of 'phi' and therefore contains a singularity
    // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
    // same direction as Rbc, and negative otherwise
    Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
    Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
    Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
    Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
    Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
    Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
    sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
          +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
    Phi=SIGN(acos(CosPhi),sign);
    SinPhi=sin(Phi);
    Phi-=parms1;
    Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
    SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
    U=0.5*parms0*SQR(Phi);
    DF=-parms0*(Phi)/SinPhi;

    UDihedralConstraints[CurrentSystem]+=U;

    // Calculate the first derivative vectors.
    d=dot_ab/rbc;
    e=dot_cd/rbc;

    dtA.x=(ds.x-CosPhi*dr.x)/r;
    dtA.y=(ds.y-CosPhi*dr.y)/r;
    dtA.z=(ds.z-CosPhi*dr.z)/r;

    dtD.x=(dr.x-CosPhi*ds.x)/s;
    dtD.y=(dr.y-CosPhi*ds.y)/s;
    dtD.z=(dr.z-CosPhi*ds.z)/s;

    dtB.x=dtA.x*(d-1.0)+e*dtD.x;
    dtB.y=dtA.y*(d-1.0)+e*dtD.y;
    dtB.z=dtA.z*(d-1.0)+e*dtD.z;

    dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
    dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
    dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

    // forces are oppositely directed to the gradient
    fa.x=-DF*dtA.x;
    fa.y=-DF*dtA.y;
    fa.z=-DF*dtA.z;

    fb.x=-DF*dtB.x;
    fb.y=-DF*dtB.y;
    fb.z=-DF*dtB.z;

    fc.x=-DF*dtC.x;
    fc.y=-DF*dtC.y;
    fc.z=-DF*dtC.z;

    fd.x=-DF*dtD.x;
    fd.y=-DF*dtD.y;
    fd.z=-DF*dtD.z;

    // add contribution to the forces
    HarmonicDihedralConstraints[CurrentSystem][m][0]->Force.x+=fa.x;
    HarmonicDihedralConstraints[CurrentSystem][m][0]->Force.y+=fa.y;
    HarmonicDihedralConstraints[CurrentSystem][m][0]->Force.z+=fa.z;

    HarmonicDihedralConstraints[CurrentSystem][m][1]->Force.x+=fb.x;
    HarmonicDihedralConstraints[CurrentSystem][m][1]->Force.y+=fb.y;
    HarmonicDihedralConstraints[CurrentSystem][m][1]->Force.z+=fb.z;

    HarmonicDihedralConstraints[CurrentSystem][m][2]->Force.x+=fc.x;
    HarmonicDihedralConstraints[CurrentSystem][m][2]->Force.y+=fc.y;
    HarmonicDihedralConstraints[CurrentSystem][m][2]->Force.z+=fc.z;

    HarmonicDihedralConstraints[CurrentSystem][m][3]->Force.x+=fd.x;
    HarmonicDihedralConstraints[CurrentSystem][m][3]->Force.y+=fd.y;
    HarmonicDihedralConstraints[CurrentSystem][m][3]->Force.z+=fd.z;

    // add contribution to the stress tensor
    // Note: rbc is here because the vector was normalized before
    StrainDerivativeTensor[CurrentSystem].ax-=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
    StrainDerivativeTensor[CurrentSystem].bx-=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
    StrainDerivativeTensor[CurrentSystem].cx-=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

    StrainDerivativeTensor[CurrentSystem].ay-=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
    StrainDerivativeTensor[CurrentSystem].by-=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
    StrainDerivativeTensor[CurrentSystem].cy-=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

    StrainDerivativeTensor[CurrentSystem].az-=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
    StrainDerivativeTensor[CurrentSystem].bz-=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
    StrainDerivativeTensor[CurrentSystem].cz-=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;
  }
}
