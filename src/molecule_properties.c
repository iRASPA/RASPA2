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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "constants.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "simulation.h"
#include "ewald.h"
#include "potentials.h"
#include "utils.h"
#include "output.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "spacegroup.h"
#include "spectra.h"

REAL ComputeBondDistanceFramework(int index)
{
  int A,B;
  VECTOR dr;
  REAL r;

  if(index>=Framework[CurrentSystem].NumberOfBonds[CurrentFramework])
    fprintf(stderr, "Error: bond index too large\n");

  A=Framework[CurrentSystem].Bonds[CurrentFramework][index].A;
  B=Framework[CurrentSystem].Bonds[CurrentFramework][index].B;

  dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
       Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
  dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
       Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
  dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
       Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
  dr=ApplyBoundaryCondition(dr);
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  return r;
}

REAL ComputeBendAngleFramework(int index)
{
  int A,B,C;
  VECTOR posA,posB,posC,Rab,Rbc;
  REAL CosTheta,Theta,rab,rbc;

  if(index>=Framework[CurrentSystem].NumberOfBends[CurrentFramework])
    fprintf(stderr, "Error: framework bend index too large\n");

  A=Framework[CurrentSystem].Bends[CurrentFramework][index].A;
  B=Framework[CurrentSystem].Bends[CurrentFramework][index].B;
  C=Framework[CurrentSystem].Bends[CurrentFramework][index].C;

  posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
  posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;
  posC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position;

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;
  Rab=ApplyBoundaryCondition(Rab);
  rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
  Rab.x/=rab;
  Rab.y/=rab;
  Rab.z/=rab;

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  Rbc=ApplyBoundaryCondition(Rbc);
  rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x/=rbc;
  Rbc.y/=rbc;
  Rbc.z/=rbc;

  CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
  CosTheta=SIGN(MIN2(fabs(CosTheta),(REAL)1.0),CosTheta);
  Theta=acos(CosTheta);

  return (RAD2DEG*Theta);
}

REAL ComputeTorsionAngleFramework(int index)
{
  int A,B,C,D;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rcd,Pb,Pc;
  REAL pb2,rpb1,pc2,rpc1,pbpc,rrbc;
  REAL cost,sint,theta;

  if(index>=Framework[CurrentSystem].NumberOfTorsions[CurrentFramework])
    fprintf(stderr, "Error: framework torsion index too large\n");

  A=Framework[CurrentSystem].Torsions[CurrentFramework][index].A;
  B=Framework[CurrentSystem].Torsions[CurrentFramework][index].B;
  C=Framework[CurrentSystem].Torsions[CurrentFramework][index].C;
  D=Framework[CurrentSystem].Torsions[CurrentFramework][index].D;

  posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
  posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;
  posC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position;
  posD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Position;

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;

  Rbc.x=posB.x-posC.x;
  Rbc.y=posB.y-posC.y;
  Rbc.z=posB.z-posC.z;

  Rcd.x=posC.x-posD.x;
  Rcd.y=posC.y-posD.y;
  Rcd.z=posC.z-posD.z;

  rrbc=1.0/sqrt(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z);

  // construct first dihedral vector
  Pb.x=Rab.y*Rbc.z-Rab.z*Rbc.y;
  Pb.y=Rab.z*Rbc.x-Rab.x*Rbc.z;
  Pb.z=Rab.x*Rbc.y-Rab.y*Rbc.x;
  pb2=Pb.x*Pb.x+Pb.y*Pb.y+Pb.z*Pb.z;
  rpb1=1.0/sqrt(pb2);

  // construct second dihedral vector
  Pc.x=Rbc.y*Rcd.z-Rbc.z*Rcd.y;
  Pc.y=Rbc.z*Rcd.x-Rbc.x*Rcd.z;
  Pc.z=Rbc.x*Rcd.y-Rbc.y*Rcd.x;
  pc2=Pc.x*Pc.x+Pc.y*Pc.y+Pc.z*Pc.z;
  rpc1=1.0/sqrt(pc2);

  // determine dihedral angle
  pbpc=Pb.x*Pc.x+Pb.y*Pc.y+Pb.z*Pc.z;
  cost=pbpc*rpb1*rpc1;
  sint=(Rbc.x*(Pc.y*Pb.z-Pc.z*Pb.y)+Rbc.y*(Pb.x*Pc.z-Pb.z*Pc.x)
        +Rbc.z*(Pc.x*Pb.y-Pc.y*Pb.x))*(rpb1*rpc1*rrbc);
  theta=atan2(sint,cost);
  if(theta<0.0) theta+=2.0*M_PI;
  return theta*180.0/M_PI;
}

REAL ComputeBondDistanceAdsorbate(int m,int index)
{
  int Type,A,B;
  VECTOR posA,posB,dr;
  REAL r;


  Type=Adsorbates[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfBonds)
    fprintf(stderr, "Error: bond index too large\n");

  A=Components[Type].Bonds[index].A;
  B=Components[Type].Bonds[index].B;
  posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
  posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
  dr.x=posA.x-posB.x;
  dr.y=posA.y-posB.y;
  dr.z=posA.z-posB.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  return r;
}

REAL ComputeBondDistanceCation(int m,int index)
{
  int Type,A,B;
  VECTOR posA,posB,dr;
  REAL r;


  Type=Cations[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfBonds)
    fprintf(stderr, "Error: bond index too large\n");

  A=Components[Type].Bonds[index].A;
  B=Components[Type].Bonds[index].B;
  posA=Cations[CurrentSystem][m].Atoms[A].Position;
  posB=Cations[CurrentSystem][m].Atoms[B].Position;
  dr.x=posA.x-posB.x;
  dr.y=posA.y-posB.y;
  dr.z=posA.z-posB.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  return r;
}

REAL ComputeUreyBradleyDistanceAdsorbate(int m,int index)
{
  int Type,A,C;
  VECTOR posA,posC,dr;
  REAL r;

  Type=Adsorbates[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfUreyBradleys)
    fprintf(stderr, "Error: Urey-Brdley index too large\n");

  A=Components[Type].UreyBradleys[index].A;
  C=Components[Type].UreyBradleys[index].C;
  posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
  posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
  dr.x=posA.x-posC.x;
  dr.y=posA.y-posC.y;
  dr.z=posA.z-posC.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  return r;
}

REAL ComputeUreyBradleyDistanceCation(int m,int index)
{
  int Type,A,C;
  VECTOR posA,posC,dr;
  REAL r;

  Type=Cations[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfUreyBradleys)
    fprintf(stderr, "Error: Urey-Bradley index too large\n");

  A=Components[Type].UreyBradleys[index].A;
  C=Components[Type].UreyBradleys[index].C;
  posA=Cations[CurrentSystem][m].Atoms[A].Position;
  posC=Cations[CurrentSystem][m].Atoms[C].Position;
  dr.x=posA.x-posC.x;
  dr.y=posA.y-posC.y;
  dr.z=posA.z-posC.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  return r;
}

REAL ComputeBendAngleAdsorbate(int m,int index)
{
  int Type,A,B,C;
  VECTOR posA,posB,posC,Rab,Rbc;
  REAL rab,rbc,theta;

  Type=Adsorbates[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfBends)
    fprintf(stderr, "Error: bend index too large\n");

  A=Components[Type].Bends[index].A;
  B=Components[Type].Bends[index].B;
  C=Components[Type].Bends[index].C;

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

  theta=acos(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
  return (RAD2DEG*theta);
}

REAL ComputeBendAngleCation(int m,int index)
{
  int Type,A,B,C;
  VECTOR posA,posB,posC,Rab,Rbc;
  REAL rab,rbc,theta;

  Type=Cations[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfBends)
    fprintf(stderr, "Error: bend index too large\n");

  A=Components[Type].Bends[index].A;
  B=Components[Type].Bends[index].B;
  C=Components[Type].Bends[index].C;

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

  theta=acos(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
  return theta*180.0/M_PI;
}

REAL ComputeTorsionAngleAdsorbate(int m,int index)
{
  int Type,A,B,C,D;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rcd,Pb,Pc;
  REAL pb2,rpb1,pc2,rpc1,pbpc,rrbc;
  REAL cost,sint,theta;

  Type=Adsorbates[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfTorsions)
    fprintf(stderr, "Error: torsion index too large\n");

  A=Components[Type].Torsions[index].A;
  B=Components[Type].Torsions[index].B;
  C=Components[Type].Torsions[index].C;
  D=Components[Type].Torsions[index].D;

  posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
  posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
  posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
  posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;

  Rbc.x=posB.x-posC.x;
  Rbc.y=posB.y-posC.y;
  Rbc.z=posB.z-posC.z;

  Rcd.x=posC.x-posD.x;
  Rcd.y=posC.y-posD.y;
  Rcd.z=posC.z-posD.z;

  rrbc=1.0/sqrt(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z);

  // construct first dihedral vector
  Pb.x=Rab.y*Rbc.z-Rab.z*Rbc.y;
  Pb.y=Rab.z*Rbc.x-Rab.x*Rbc.z;
  Pb.z=Rab.x*Rbc.y-Rab.y*Rbc.x;
  pb2=Pb.x*Pb.x+Pb.y*Pb.y+Pb.z*Pb.z;
  rpb1=1.0/sqrt(pb2);

  // construct second dihedral vector
  Pc.x=Rbc.y*Rcd.z-Rbc.z*Rcd.y;
  Pc.y=Rbc.z*Rcd.x-Rbc.x*Rcd.z;
  Pc.z=Rbc.x*Rcd.y-Rbc.y*Rcd.x;
  pc2=Pc.x*Pc.x+Pc.y*Pc.y+Pc.z*Pc.z;
  rpc1=1.0/sqrt(pc2);

  // determine dihedral angle
  pbpc=Pb.x*Pc.x+Pb.y*Pc.y+Pb.z*Pc.z;
  cost=pbpc*rpb1*rpc1;
  sint=(Rbc.x*(Pc.y*Pb.z-Pc.z*Pb.y)+Rbc.y*(Pb.x*Pc.z-Pb.z*Pc.x)
        +Rbc.z*(Pc.x*Pb.y-Pc.y*Pb.x))*(rpb1*rpc1*rrbc);
  theta=atan2(sint,cost);
  if(theta<0.0) theta+=2.0*M_PI;
  return theta*180.0/M_PI;
}

REAL ComputeTorsionAngleCation(int m,int index)
{
  int Type,A,B,C,D;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rcd,Pb,Pc;
  REAL pb2,rpb1,pc2,rpc1,pbpc,rrbc;
  REAL cost,sint,theta;

  Type=Cations[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfTorsions)
    fprintf(stderr, "Error: torsion index too large\n");

  A=Components[Type].Torsions[index].A;
  B=Components[Type].Torsions[index].B;
  C=Components[Type].Torsions[index].C;
  D=Components[Type].Torsions[index].D;

  posA=Cations[CurrentSystem][m].Atoms[A].Position;
  posB=Cations[CurrentSystem][m].Atoms[B].Position;
  posC=Cations[CurrentSystem][m].Atoms[C].Position;
  posD=Cations[CurrentSystem][m].Atoms[D].Position;

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;

  Rbc.x=posB.x-posC.x;
  Rbc.y=posB.y-posC.y;
  Rbc.z=posB.z-posC.z;

  Rcd.x=posC.x-posD.x;
  Rcd.y=posC.y-posD.y;
  Rcd.z=posC.z-posD.z;

  rrbc=1.0/sqrt(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z);

  // construct first dihedral vector
  Pb.x=Rab.y*Rbc.z-Rab.z*Rbc.y;
  Pb.y=Rab.z*Rbc.x-Rab.x*Rbc.z;
  Pb.z=Rab.x*Rbc.y-Rab.y*Rbc.x;
  pb2=Pb.x*Pb.x+Pb.y*Pb.y+Pb.z*Pb.z;
  rpb1=1.0/sqrt(pb2);

  // construct second dihedral vector
  Pc.x=Rbc.y*Rcd.z-Rbc.z*Rcd.y;
  Pc.y=Rbc.z*Rcd.x-Rbc.x*Rcd.z;
  Pc.z=Rbc.x*Rcd.y-Rbc.y*Rcd.x;
  pc2=Pc.x*Pc.x+Pc.y*Pc.y+Pc.z*Pc.z;
  rpc1=1.0/sqrt(pc2);

  // determine dihedral angle
  pbpc=Pb.x*Pc.x+Pb.y*Pc.y+Pb.z*Pc.z;
  cost=pbpc*rpb1*rpc1;
  sint=(Rbc.x*(Pc.y*Pb.z-Pc.z*Pb.y)+Rbc.y*(Pb.x*Pc.z-Pb.z*Pc.x)
        +Rbc.z*(Pc.x*Pb.y-Pc.y*Pb.x))*(rpb1*rpc1*rrbc);
  theta=atan2(sint,cost);
  if(theta<0.0) theta+=2.0*M_PI;
  return theta*180.0/M_PI;
}

int ReturnTorsionConformation(REAL theta)
{
  int conf;

  conf=0;
  if(((theta>=0.0)&&(theta<30.0))||((theta>=330.0)&&(theta<=360.0)))
    conf=SYNPERIPLANAR;
  else if ((theta>=30.0)&&(theta<90.0))
    conf=SYNCLINAL_PLUS;
  else if ((theta>=90.0)&&(theta<150.0))
    conf=ANTICLINAL_PLUS;
  else if((theta>=150.0)&(theta<210.0))
    conf=ANTIPERIPLANAR_PLUS;
  else if((theta>=210.0)&&(theta<270.0))
    conf=ANTICLINAL_MIN;
  else if((theta>=270.0)&&(theta<330.0))
    conf=SYNCLINAL_MIN;

  return conf;
}
