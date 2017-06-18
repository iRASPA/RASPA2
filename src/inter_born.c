/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'inter_born.c' is part of RASPA-2.0

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
#include <string.h>
#include <sys/stat.h>
#include "constants.h"
#include "utils.h"
#include "simulation.h"
#include "potentials.h"
#include "output.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "grids.h"
#include "ewald.h"
#include "inter_energy.h"
#include "inter_force.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "spectra.h"


void ComputeInterVDWBornTerm(void)
{
  int I,J,i,j,ig,jg,ia,ja;
  int typeA,typeB;
  int TypeMolA,TypeMolB;
  REAL rr;
  REAL energy;
  VECTOR posA,posB,dr;
  VECTOR comA,comB;
  REAL DF,DDF;
  int RigidI,RigidJ;
  VECTOR drJI,f;
  REAL scalingA,scalingB;

  // first loop over adsorbate molecules
  for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
  {
    TypeMolA=Adsorbates[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      RigidI=Components[TypeMolA].Groups[ig].Rigid;
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(RigidI)
          comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;
        else
          comA=Adsorbates[CurrentSystem][I].Atoms[i].Position;

        typeA=Adsorbates[CurrentSystem][I].Atoms[i].Type;
        posA=Adsorbates[CurrentSystem][I].Atoms[i].Position;
        scalingA=Adsorbates[CurrentSystem][I].Atoms[i].CFVDWScalingParameter;

        if(!OmitAdsorbateAdsorbateVDWInteractions)
        {
          // second loop over adsorbates
          for(J=I+1;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
          {
            TypeMolB=Adsorbates[CurrentSystem][J].Type;
            for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
            {
              RigidJ=Components[TypeMolB].Groups[jg].Rigid;
              for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
              {
                j=Components[TypeMolB].Groups[jg].Atoms[ja];

                if(RigidJ)
                  comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                else
                  comB=Adsorbates[CurrentSystem][J].Atoms[j].Position;

                typeB=Adsorbates[CurrentSystem][J].Atoms[j].Type;
                posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;
                scalingB=Adsorbates[CurrentSystem][J].Atoms[j].CFVDWScalingParameter;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                {
                  PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,scalingA*scalingB);

                  // add contribution to the energy
                  UAdsorbateAdsorbateVDW[CurrentSystem]+=energy;

                  f.x=DF*dr.x;
                  f.y=DF*dr.y;
                  f.z=DF*dr.z;

                  // add contribution to the first derivatives
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.x-=f.x;
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.y-=f.y;
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.z-=f.z;

                  Adsorbates[CurrentSystem][J].Atoms[j].Force.x+=f.x;
                  Adsorbates[CurrentSystem][J].Atoms[j].Force.y+=f.y;
                  Adsorbates[CurrentSystem][J].Atoms[j].Force.z+=f.z;

                  // add contribution to the strain derivative tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                  VECTOR dJI;

                  dJI.x=(posB.x-comB.x)-(posA.x-comA.x);
                  dJI.y=(posB.y-comB.y)-(posA.y-comA.y);
                  dJI.z=(posB.z-comB.z)-(posA.z-comA.z);

                  // the corrected com-com distance vector
                  drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
                  drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
                  drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

                  // add the contributions to the Born-term

                  BornTerm[CurrentSystem].xxxx+=DDF*drJI.x*dr.x*drJI.x*dr.x+DF*drJI.x*dr.x+DF*drJI.x*drJI.x;
                  BornTerm[CurrentSystem].xxyy+=DDF*drJI.x*dr.x*drJI.y*dr.y;
                  BornTerm[CurrentSystem].xxzz+=DDF*drJI.x*dr.x*drJI.z*dr.z;
                  BornTerm[CurrentSystem].xxyz+=0.5*DDF*drJI.x*dr.x*(drJI.y*dr.z+drJI.z*dr.y);
                  BornTerm[CurrentSystem].xxzx+=0.5*DDF*drJI.x*dr.x*(drJI.x*dr.z+drJI.z*dr.x)+0.5*DF*(drJI.x*dJI.z)+DF*(drJI.x*dr.z);
                  BornTerm[CurrentSystem].xxxy+=0.5*DDF*drJI.x*dr.x*(drJI.x*dr.y+drJI.y*dr.x)+0.5*DF*(drJI.x*dJI.y)+DF*drJI.x*dr.y;

                  BornTerm[CurrentSystem].yyyy+=DDF*drJI.y*dr.y*dr.y*dr.y+2.0*DF*drJI.y*dr.y+DDF*drJI.y*dr.y*dJI.y*dr.y+DF*drJI.y*dJI.y;
                  BornTerm[CurrentSystem].yyzz+=0.5*DDF*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.z*dr.z+dJI.z*dr.z);
                  BornTerm[CurrentSystem].yyyz+=0.5*DDF*(drJI.y*dr.y*dr.y*dr.z+dr.y*drJI.y*dr.y*dr.z)+DF*(drJI.y*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.y*dr.z+dJI.z*dr.y)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].yyzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.y+dr.x*drJI.z*dr.y*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
                  BornTerm[CurrentSystem].yyxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.y+dr.x*drJI.y*dr.y*dr.y)+DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.y+dJI.y*dr.y)+0.5*DF*(drJI.x*dJI.y);

                  BornTerm[CurrentSystem].zzzz+=DDF*drJI.z*dr.z*dr.z*dr.z+2.0*DF*drJI.z*dr.z+DDF*drJI.z*dr.z*dJI.z*dr.z+DF*drJI.z*dJI.z;
                  BornTerm[CurrentSystem].zzyz+=0.5*DDF*(drJI.y*dr.z*dr.z*dr.z+dr.y*drJI.z*dr.z*dr.z)+DF*(drJI.y*dr.z)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].zzzx+=0.5*DDF*(drJI.x*dr.z*dr.z*dr.z+dr.x*drJI.z*dr.z*dr.z)+DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.x*dJI.z);
                  BornTerm[CurrentSystem].zzxy+=0.5*DDF*(drJI.x*dr.y*dr.z*dr.z+dr.x*drJI.y*dr.z*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);

                  BornTerm[CurrentSystem].yzyz+=0.5*DDF*(drJI.y*dr.z*dr.y*dr.z+dr.y*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.y*dr.y)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.z*dJI.z+drJI.y*dJI.y);
                  BornTerm[CurrentSystem].yzzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.z+dr.x*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.y);
                  BornTerm[CurrentSystem].yzxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.z+dr.x*drJI.y*dr.y*dr.z)+0.5*DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.z);

                  BornTerm[CurrentSystem].zxzx+=0.5*DDF*(drJI.x*dr.z*dr.x*dr.z+dr.x*drJI.z*dr.x*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.z*dJI.z+drJI.x*dJI.x);
                  BornTerm[CurrentSystem].zxxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.z+dr.x*drJI.y*dr.x*dr.z)+0.5*DF*(drJI.y*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.y*dJI.z);

                  BornTerm[CurrentSystem].xyxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.y+dr.x*drJI.y*dr.x*dr.y)+0.5*DF*(drJI.y*dr.y+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.y+dJI.y*dr.x)+0.25*DF*(drJI.y*dJI.y+drJI.x*dJI.x);

                }
              }
            }
          }
        }

        // second loop over cations
        if(!OmitAdsorbateCationVDWInteractions)
        {
          for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
          {
            TypeMolB=Cations[CurrentSystem][J].Type;
            for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
            {
              RigidJ=Components[TypeMolB].Groups[jg].Rigid;
              for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
              {
                j=Components[TypeMolB].Groups[jg].Atoms[ja];

                if(RigidJ)
                  comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                else
                  comB=Cations[CurrentSystem][J].Atoms[j].Position;

                typeB=Cations[CurrentSystem][J].Atoms[j].Type;
                posB=Cations[CurrentSystem][J].Atoms[j].Position;
                scalingB=Adsorbates[CurrentSystem][J].Atoms[j].CFVDWScalingParameter;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                {
                  PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,scalingA*scalingB);

                  // add contribution to the energy
                  UAdsorbateCationVDW[CurrentSystem]+=energy;

                  f.x=DF*dr.x;
                  f.y=DF*dr.y;
                  f.z=DF*dr.z;

                  // add contribution to the first derivatives
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.x-=f.x;
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.y-=f.y;
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.z-=f.z;

                  Cations[CurrentSystem][J].Atoms[j].Force.x+=f.x;
                  Cations[CurrentSystem][J].Atoms[j].Force.y+=f.y;
                  Cations[CurrentSystem][J].Atoms[j].Force.z+=f.z;

                  // add contribution to the strain derivative tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                  VECTOR dJI;

                  dJI.x=(posB.x-comB.x)-(posA.x-comA.x);
                  dJI.y=(posB.y-comB.y)-(posA.y-comA.y);
                  dJI.z=(posB.z-comB.z)-(posA.z-comA.z);

                  // the corrected com-com distance vector
                  drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
                  drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
                  drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

                  // add the contributions to the Born-term
                  BornTerm[CurrentSystem].xxxx+=DDF*drJI.x*dr.x*drJI.x*dr.x+DF*drJI.x*dr.x+DF*drJI.x*drJI.x;
                  BornTerm[CurrentSystem].xxyy+=DDF*drJI.x*dr.x*drJI.y*dr.y;
                  BornTerm[CurrentSystem].xxzz+=DDF*drJI.x*dr.x*drJI.z*dr.z;
                  BornTerm[CurrentSystem].xxyz+=0.5*DDF*drJI.x*dr.x*(drJI.y*dr.z+drJI.z*dr.y);
                  BornTerm[CurrentSystem].xxzx+=0.5*DDF*drJI.x*dr.x*(drJI.x*dr.z+drJI.z*dr.x)+0.5*DF*(drJI.x*dJI.z)+DF*(drJI.x*dr.z);
                  BornTerm[CurrentSystem].xxxy+=0.5*DDF*drJI.x*dr.x*(drJI.x*dr.y+drJI.y*dr.x)+0.5*DF*(drJI.x*dJI.y)+DF*drJI.x*dr.y;

                  BornTerm[CurrentSystem].yyyy+=DDF*drJI.y*dr.y*dr.y*dr.y+2.0*DF*drJI.y*dr.y+DDF*drJI.y*dr.y*dJI.y*dr.y+DF*drJI.y*dJI.y;
                  BornTerm[CurrentSystem].yyzz+=0.5*DDF*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.z*dr.z+dJI.z*dr.z);
                  BornTerm[CurrentSystem].yyyz+=0.5*DDF*(drJI.y*dr.y*dr.y*dr.z+dr.y*drJI.y*dr.y*dr.z)+DF*(drJI.y*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.y*dr.z+dJI.z*dr.y)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].yyzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.y+dr.x*drJI.z*dr.y*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
                  BornTerm[CurrentSystem].yyxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.y+dr.x*drJI.y*dr.y*dr.y)+DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.y+dJI.y*dr.y)+0.5*DF*(drJI.x*dJI.y);

                  BornTerm[CurrentSystem].zzzz+=DDF*drJI.z*dr.z*dr.z*dr.z+2.0*DF*drJI.z*dr.z+DDF*drJI.z*dr.z*dJI.z*dr.z+DF*drJI.z*dJI.z;
                  BornTerm[CurrentSystem].zzyz+=0.5*DDF*(drJI.y*dr.z*dr.z*dr.z+dr.y*drJI.z*dr.z*dr.z)+DF*(drJI.y*dr.z)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].zzzx+=0.5*DDF*(drJI.x*dr.z*dr.z*dr.z+dr.x*drJI.z*dr.z*dr.z)+DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.x*dJI.z);
                  BornTerm[CurrentSystem].zzxy+=0.5*DDF*(drJI.x*dr.y*dr.z*dr.z+dr.x*drJI.y*dr.z*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);

                  BornTerm[CurrentSystem].yzyz+=0.5*DDF*(drJI.y*dr.z*dr.y*dr.z+dr.y*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.y*dr.y)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.z*dJI.z+drJI.y*dJI.y);
                  BornTerm[CurrentSystem].yzzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.z+dr.x*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.y);
                  BornTerm[CurrentSystem].yzxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.z+dr.x*drJI.y*dr.y*dr.z)+0.5*DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.z);

                  BornTerm[CurrentSystem].zxzx+=0.5*DDF*(drJI.x*dr.z*dr.x*dr.z+dr.x*drJI.z*dr.x*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.z*dJI.z+drJI.x*dJI.x);
                  BornTerm[CurrentSystem].zxxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.z+dr.x*drJI.y*dr.x*dr.z)+0.5*DF*(drJI.y*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.y*dJI.z);

                  BornTerm[CurrentSystem].xyxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.y+dr.x*drJI.y*dr.x*dr.y)+0.5*DF*(drJI.y*dr.y+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.y+dJI.y*dr.x)+0.25*DF*(drJI.y*dJI.y+drJI.x*dJI.x);

                }
              }
            }
          }
        }
      }
    }
  }

  if(!OmitCationCationVDWInteractions)
  {
    // first loop over cation molecules
    for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
    {
      TypeMolA=Cations[CurrentSystem][I].Type;
      for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
      {
        RigidI=Components[TypeMolA].Groups[ig].Rigid;
        for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
        {
          i=Components[TypeMolA].Groups[ig].Atoms[ia];

          if(RigidI)
            comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;
          else
            comA=Cations[CurrentSystem][I].Atoms[i].Position;

          typeA=Cations[CurrentSystem][I].Atoms[i].Type;
          posA=Cations[CurrentSystem][I].Atoms[i].Position;
          scalingA=Cations[CurrentSystem][I].Atoms[i].CFVDWScalingParameter;

          // second loop over cation
          for(J=I+1;J<NumberOfCationMolecules[CurrentSystem];J++)
          {
            TypeMolB=Cations[CurrentSystem][J].Type;
            for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
            {
              RigidJ=Components[TypeMolB].Groups[jg].Rigid;
              for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
              {
                j=Components[TypeMolB].Groups[jg].Atoms[ja];

                if(RigidJ)
                  comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                else
                  comB=Cations[CurrentSystem][J].Atoms[j].Position;

                typeB=Cations[CurrentSystem][J].Atoms[j].Type;
                posB=Cations[CurrentSystem][J].Atoms[j].Position;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                {
                  scalingB=Adsorbates[CurrentSystem][J].Atoms[j].CFVDWScalingParameter;

                  PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,scalingA*scalingB);

                  // add contribution to the energy
                  UCationCationVDW[CurrentSystem]+=energy;

                  f.x=DF*dr.x;
                  f.y=DF*dr.y;
                  f.z=DF*dr.z;

                  // add contribution to the first derivatives
                  Cations[CurrentSystem][I].Atoms[i].Force.x-=f.x;
                  Cations[CurrentSystem][I].Atoms[i].Force.y-=f.y;
                  Cations[CurrentSystem][I].Atoms[i].Force.z-=f.z;

                  Cations[CurrentSystem][J].Atoms[j].Force.x+=f.x;
                  Cations[CurrentSystem][J].Atoms[j].Force.y+=f.y;
                  Cations[CurrentSystem][J].Atoms[j].Force.z+=f.z;

                  // add contribution to the strain derivative tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                  VECTOR dJI;

                  dJI.x=(posB.x-comB.x)-(posA.x-comA.x);
                  dJI.y=(posB.y-comB.y)-(posA.y-comA.y);
                  dJI.z=(posB.z-comB.z)-(posA.z-comA.z);

                  // the corrected com-com distance vector
                  drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
                  drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
                  drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

                  // add the contributions to the Born-term
                  BornTerm[CurrentSystem].xxxx+=DDF*drJI.x*dr.x*dr.x*dr.x+2.0*DF*drJI.x*dr.x+DDF*drJI.x*dr.x*dJI.x*dr.x+DF*drJI.x*dJI.x;
                  BornTerm[CurrentSystem].xxyy+=DDF*drJI.x*dr.x*dr.y*dr.y+0.5*DDF*drJI.x*dr.x*(dJI.y*dr.y+dJI.y*dr.y);
                  BornTerm[CurrentSystem].xxzz+=DDF*drJI.x*dr.x*dr.z*dr.z+0.5*DDF*drJI.x*dr.x*(dJI.z*dr.z+dJI.z*dr.z);
                  BornTerm[CurrentSystem].xxyz+=DDF*drJI.x*dr.x*dr.y*dr.z+0.5*DDF*drJI.x*dr.x*(dJI.y*dr.z+dJI.z*dr.y);
                  BornTerm[CurrentSystem].xxzx+=DDF*drJI.x*dr.x*dr.x*dr.z+DF*(drJI.x*dr.z)+0.5*DDF*drJI.x*dr.x*(dJI.x*dr.z+dJI.z*dr.x)+0.5*DF*(drJI.x*dJI.z);
                  BornTerm[CurrentSystem].xxxy+=DDF*drJI.x*dr.x*dr.x*dr.y+DF*drJI.x*dr.y+0.5*DDF*drJI.x*dr.x*(dJI.x*dr.y+dJI.y*dr.x)+0.5*DF*(drJI.x*dJI.y);

                  BornTerm[CurrentSystem].yyyy+=DDF*drJI.y*dr.y*dr.y*dr.y+2.0*DF*drJI.y*dr.y+DDF*drJI.y*dr.y*dJI.y*dr.y+DF*drJI.y*dJI.y;
                  BornTerm[CurrentSystem].yyzz+=0.5*DDF*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.z*dr.z+dJI.z*dr.z);
                  BornTerm[CurrentSystem].yyyz+=0.5*DDF*(drJI.y*dr.y*dr.y*dr.z+dr.y*drJI.y*dr.y*dr.z)+DF*(drJI.y*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.y*dr.z+dJI.z*dr.y)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].yyzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.y+dr.x*drJI.z*dr.y*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
                  BornTerm[CurrentSystem].yyxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.y+dr.x*drJI.y*dr.y*dr.y)+DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.y+dJI.y*dr.y)+0.5*DF*(drJI.x*dJI.y);

                  BornTerm[CurrentSystem].zzzz+=DDF*drJI.z*dr.z*dr.z*dr.z+2.0*DF*drJI.z*dr.z+DDF*drJI.z*dr.z*dJI.z*dr.z+DF*drJI.z*dJI.z;
                  BornTerm[CurrentSystem].zzyz+=0.5*DDF*(drJI.y*dr.z*dr.z*dr.z+dr.y*drJI.z*dr.z*dr.z)+DF*(drJI.y*dr.z)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].zzzx+=0.5*DDF*(drJI.x*dr.z*dr.z*dr.z+dr.x*drJI.z*dr.z*dr.z)+DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.x*dJI.z);
                  BornTerm[CurrentSystem].zzxy+=0.5*DDF*(drJI.x*dr.y*dr.z*dr.z+dr.x*drJI.y*dr.z*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);

                  BornTerm[CurrentSystem].yzyz+=0.5*DDF*(drJI.y*dr.z*dr.y*dr.z+dr.y*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.y*dr.y)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.z*dJI.z+drJI.y*dJI.y);
                  BornTerm[CurrentSystem].yzzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.z+dr.x*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.y);
                  BornTerm[CurrentSystem].yzxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.z+dr.x*drJI.y*dr.y*dr.z)+0.5*DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.z);

                  BornTerm[CurrentSystem].zxzx+=0.5*DDF*(drJI.x*dr.z*dr.x*dr.z+dr.x*drJI.z*dr.x*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.z*dJI.z+drJI.x*dJI.x);
                  BornTerm[CurrentSystem].zxxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.z+dr.x*drJI.y*dr.x*dr.z)+0.5*DF*(drJI.y*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.y*dJI.z);

                  BornTerm[CurrentSystem].xyxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.y+dr.x*drJI.y*dr.x*dr.y)+0.5*DF*(drJI.y*dr.y+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.y+dJI.y*dr.x)+0.25*DF*(drJI.y*dJI.y+drJI.x*dJI.x);
                }
              }
            }
          }
        }
      }
    }
  }
}


void ComputeInterChargeChargeBornTerm(void)
{
  int I,J,i,j,ig,jg,ia,ja;
  int typeA,typeB;
  int TypeMolA,TypeMolB;
  REAL rr,energy;
  REAL chargeA,chargeB;
  VECTOR posA,posB,dr;
  VECTOR comA,comB;
  REAL DF,DDF;
  int RigidI,RigidJ;
  VECTOR drJI,f;

  if(ChargeMethod==NONE) return;

  // first loop over adsorbate molecules
  for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
  {
    TypeMolA=Adsorbates[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      RigidI=Components[TypeMolA].Groups[ig].Rigid;
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(RigidI)
          comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;
        else
          comA=Adsorbates[CurrentSystem][I].Atoms[i].Position;

        typeA=Adsorbates[CurrentSystem][I].Atoms[i].Type;
        chargeA=Adsorbates[CurrentSystem][I].Atoms[i].Charge;
        posA=Adsorbates[CurrentSystem][I].Atoms[i].Position;

        if(!OmitAdsorbateAdsorbateCoulombInteractions)
        {
          // second loop over adsorbates
          for(J=I+1;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
          {
            TypeMolB=Adsorbates[CurrentSystem][J].Type;
            for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
            {
              RigidJ=Components[TypeMolB].Groups[jg].Rigid;
              for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
              {
                j=Components[TypeMolB].Groups[jg].Atoms[ja];

                if(RigidJ)
                  comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                else
                  comB=Adsorbates[CurrentSystem][J].Atoms[j].Position;

                typeB=Adsorbates[CurrentSystem][J].Atoms[j].Type;
                chargeB=Adsorbates[CurrentSystem][J].Atoms[j].Charge;
                posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeChargeSquared[CurrentSystem])
                {
                  PotentialSecondDerivativeCoulombic(chargeA,chargeB,rr,&energy,&DF,&DDF);

                  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=energy;

                  f.x=DF*dr.x;
                  f.y=DF*dr.y;
                  f.z=DF*dr.z;

                  // add contribution to the first derivatives
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.x-=f.x;
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.y-=f.y;
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.z-=f.z;

                  Adsorbates[CurrentSystem][J].Atoms[j].Force.x+=f.x;
                  Adsorbates[CurrentSystem][J].Atoms[j].Force.y+=f.y;
                  Adsorbates[CurrentSystem][J].Atoms[j].Force.z+=f.z;

                  // add contribution to the strain derivative tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                  VECTOR dJI;

                  dJI.x=(posB.x-comB.x)-(posA.x-comA.x);
                  dJI.y=(posB.y-comB.y)-(posA.y-comA.y);
                  dJI.z=(posB.z-comB.z)-(posA.z-comA.z);

                  // the corrected com-com distance vector
                  drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
                  drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
                  drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

                  // add the contributions to the Born-term
                  BornTerm[CurrentSystem].xxxx+=DDF*drJI.x*dr.x*drJI.x*dr.x+DF*drJI.x*dr.x+DF*drJI.x*drJI.x;
                  BornTerm[CurrentSystem].xxyy+=DDF*drJI.x*dr.x*drJI.y*dr.y;
                  BornTerm[CurrentSystem].xxzz+=DDF*drJI.x*dr.x*drJI.z*dr.z;
                  BornTerm[CurrentSystem].xxyz+=0.5*DDF*drJI.x*dr.x*(drJI.y*dr.z+drJI.z*dr.y);
                  BornTerm[CurrentSystem].xxzx+=0.5*DDF*drJI.x*dr.x*(drJI.x*dr.z+drJI.z*dr.x)+0.5*DF*(drJI.x*dJI.z)+DF*(drJI.x*dr.z);
                  BornTerm[CurrentSystem].xxxy+=0.5*DDF*drJI.x*dr.x*(drJI.x*dr.y+drJI.y*dr.x)+0.5*DF*(drJI.x*dJI.y)+DF*drJI.x*dr.y;

                  BornTerm[CurrentSystem].yyyy+=DDF*drJI.y*dr.y*dr.y*dr.y+2.0*DF*drJI.y*dr.y+DDF*drJI.y*dr.y*dJI.y*dr.y+DF*drJI.y*dJI.y;
                  BornTerm[CurrentSystem].yyzz+=0.5*DDF*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.z*dr.z+dJI.z*dr.z);
                  BornTerm[CurrentSystem].yyyz+=0.5*DDF*(drJI.y*dr.y*dr.y*dr.z+dr.y*drJI.y*dr.y*dr.z)+DF*(drJI.y*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.y*dr.z+dJI.z*dr.y)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].yyzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.y+dr.x*drJI.z*dr.y*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
                  BornTerm[CurrentSystem].yyxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.y+dr.x*drJI.y*dr.y*dr.y)+DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.y+dJI.y*dr.y)+0.5*DF*(drJI.x*dJI.y);

                  BornTerm[CurrentSystem].zzzz+=DDF*drJI.z*dr.z*dr.z*dr.z+2.0*DF*drJI.z*dr.z+DDF*drJI.z*dr.z*dJI.z*dr.z+DF*drJI.z*dJI.z;
                  BornTerm[CurrentSystem].zzyz+=0.5*DDF*(drJI.y*dr.z*dr.z*dr.z+dr.y*drJI.z*dr.z*dr.z)+DF*(drJI.y*dr.z)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].zzzx+=0.5*DDF*(drJI.x*dr.z*dr.z*dr.z+dr.x*drJI.z*dr.z*dr.z)+DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.x*dJI.z);
                  BornTerm[CurrentSystem].zzxy+=0.5*DDF*(drJI.x*dr.y*dr.z*dr.z+dr.x*drJI.y*dr.z*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);

                  BornTerm[CurrentSystem].yzyz+=0.5*DDF*(drJI.y*dr.z*dr.y*dr.z+dr.y*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.y*dr.y)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.z*dJI.z+drJI.y*dJI.y);
                  BornTerm[CurrentSystem].yzzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.z+dr.x*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.y);
                  BornTerm[CurrentSystem].yzxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.z+dr.x*drJI.y*dr.y*dr.z)+0.5*DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.z);

                  BornTerm[CurrentSystem].zxzx+=0.5*DDF*(drJI.x*dr.z*dr.x*dr.z+dr.x*drJI.z*dr.x*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.z*dJI.z+drJI.x*dJI.x);
                  BornTerm[CurrentSystem].zxxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.z+dr.x*drJI.y*dr.x*dr.z)+0.5*DF*(drJI.y*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.y*dJI.z);

                  BornTerm[CurrentSystem].xyxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.y+dr.x*drJI.y*dr.x*dr.y)+0.5*DF*(drJI.y*dr.y+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.y+dJI.y*dr.x)+0.25*DF*(drJI.y*dJI.y+drJI.x*dJI.x);

                }
              }
            }
          }
        }

        // second loop over cations
        if(!OmitAdsorbateCationCoulombInteractions)
        {
          for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
          {
            TypeMolB=Cations[CurrentSystem][J].Type;
            for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
            {
              RigidJ=Components[TypeMolB].Groups[jg].Rigid;
              for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
              {
                j=Components[TypeMolB].Groups[jg].Atoms[ja];

                if(RigidJ)
                  comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                else
                  comB=Cations[CurrentSystem][J].Atoms[j].Position;

                typeB=Cations[CurrentSystem][J].Atoms[j].Type;
                chargeB=Cations[CurrentSystem][J].Atoms[j].Charge;
                posB=Cations[CurrentSystem][J].Atoms[j].Position;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeChargeSquared[CurrentSystem])
                {
                  PotentialSecondDerivativeCoulombic(chargeA,chargeB,rr,&energy,&DF,&DDF);

                  UAdsorbateCationChargeChargeReal[CurrentSystem]+=energy;

                  f.x=DF*dr.x;
                  f.y=DF*dr.y;
                  f.z=DF*dr.z;

                  // add contribution to the first derivatives
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.x-=f.x;
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.y-=f.y;
                  Adsorbates[CurrentSystem][I].Atoms[i].Force.z-=f.z;

                  Cations[CurrentSystem][J].Atoms[j].Force.x+=f.x;
                  Cations[CurrentSystem][J].Atoms[j].Force.y+=f.y;
                  Cations[CurrentSystem][J].Atoms[j].Force.z+=f.z;

                  // add contribution to the strain derivative tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                  VECTOR dJI;

                  dJI.x=(posB.x-comB.x)-(posA.x-comA.x);
                  dJI.y=(posB.y-comB.y)-(posA.y-comA.y);
                  dJI.z=(posB.z-comB.z)-(posA.z-comA.z);

                  // the corrected com-com distance vector
                  drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
                  drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
                  drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

                  // add the contributions to the Born-term
                  BornTerm[CurrentSystem].xxxx+=DDF*drJI.x*dr.x*drJI.x*dr.x+DF*drJI.x*dr.x+DF*drJI.x*drJI.x;
                  BornTerm[CurrentSystem].xxyy+=DDF*drJI.x*dr.x*drJI.y*dr.y;
                  BornTerm[CurrentSystem].xxzz+=DDF*drJI.x*dr.x*drJI.z*dr.z;
                  BornTerm[CurrentSystem].xxyz+=0.5*DDF*drJI.x*dr.x*(drJI.y*dr.z+drJI.z*dr.y);
                  BornTerm[CurrentSystem].xxzx+=0.5*DDF*drJI.x*dr.x*(drJI.x*dr.z+drJI.z*dr.x)+0.5*DF*(drJI.x*dJI.z)+DF*(drJI.x*dr.z);
                  BornTerm[CurrentSystem].xxxy+=0.5*DDF*drJI.x*dr.x*(drJI.x*dr.y+drJI.y*dr.x)+0.5*DF*(drJI.x*dJI.y)+DF*drJI.x*dr.y;

                  BornTerm[CurrentSystem].yyyy+=DDF*drJI.y*dr.y*dr.y*dr.y+2.0*DF*drJI.y*dr.y+DDF*drJI.y*dr.y*dJI.y*dr.y+DF*drJI.y*dJI.y;
                  BornTerm[CurrentSystem].yyzz+=0.5*DDF*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.z*dr.z+dJI.z*dr.z);
                  BornTerm[CurrentSystem].yyyz+=0.5*DDF*(drJI.y*dr.y*dr.y*dr.z+dr.y*drJI.y*dr.y*dr.z)+DF*(drJI.y*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.y*dr.z+dJI.z*dr.y)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].yyzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.y+dr.x*drJI.z*dr.y*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
                  BornTerm[CurrentSystem].yyxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.y+dr.x*drJI.y*dr.y*dr.y)+DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.y+dJI.y*dr.y)+0.5*DF*(drJI.x*dJI.y);

                  BornTerm[CurrentSystem].zzzz+=DDF*drJI.z*dr.z*dr.z*dr.z+2.0*DF*drJI.z*dr.z+DDF*drJI.z*dr.z*dJI.z*dr.z+DF*drJI.z*dJI.z;
                  BornTerm[CurrentSystem].zzyz+=0.5*DDF*(drJI.y*dr.z*dr.z*dr.z+dr.y*drJI.z*dr.z*dr.z)+DF*(drJI.y*dr.z)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].zzzx+=0.5*DDF*(drJI.x*dr.z*dr.z*dr.z+dr.x*drJI.z*dr.z*dr.z)+DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.x*dJI.z);
                  BornTerm[CurrentSystem].zzxy+=0.5*DDF*(drJI.x*dr.y*dr.z*dr.z+dr.x*drJI.y*dr.z*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);

                  BornTerm[CurrentSystem].yzyz+=0.5*DDF*(drJI.y*dr.z*dr.y*dr.z+dr.y*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.y*dr.y)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.z*dJI.z+drJI.y*dJI.y);
                  BornTerm[CurrentSystem].yzzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.z+dr.x*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.y);
                  BornTerm[CurrentSystem].yzxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.z+dr.x*drJI.y*dr.y*dr.z)+0.5*DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.z);

                  BornTerm[CurrentSystem].zxzx+=0.5*DDF*(drJI.x*dr.z*dr.x*dr.z+dr.x*drJI.z*dr.x*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.z*dJI.z+drJI.x*dJI.x);
                  BornTerm[CurrentSystem].zxxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.z+dr.x*drJI.y*dr.x*dr.z)+0.5*DF*(drJI.y*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.y*dJI.z);

                  BornTerm[CurrentSystem].xyxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.y+dr.x*drJI.y*dr.x*dr.y)+0.5*DF*(drJI.y*dr.y+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.y+dJI.y*dr.x)+0.25*DF*(drJI.y*dJI.y+drJI.x*dJI.x);

                }
              }
            }
          }
        }
      }
    }
  }

  if(!OmitCationCationCoulombInteractions)
  {
    // first loop over cation molecules
    for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
    {
      TypeMolA=Cations[CurrentSystem][I].Type;
      for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
      {
        RigidI=Components[TypeMolA].Groups[ig].Rigid;
        for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
        {
          i=Components[TypeMolA].Groups[ig].Atoms[ia];

          if(RigidI)
            comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;
          else
            comA=Cations[CurrentSystem][I].Atoms[i].Position;

          typeA=Cations[CurrentSystem][I].Atoms[i].Type;
          chargeA=Cations[CurrentSystem][I].Atoms[i].Charge;
          posA=Cations[CurrentSystem][I].Atoms[i].Position;

          // second loop over cation
          for(J=I+1;J<NumberOfCationMolecules[CurrentSystem];J++)
          {
            TypeMolB=Cations[CurrentSystem][J].Type;
            for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
            {
              RigidJ=Components[TypeMolB].Groups[jg].Rigid;
              for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
              {
                j=Components[TypeMolB].Groups[jg].Atoms[ja];

                if(RigidJ)
                  comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                else
                  comB=Cations[CurrentSystem][J].Atoms[j].Position;

                typeB=Cations[CurrentSystem][J].Atoms[j].Type;
                chargeB=Cations[CurrentSystem][J].Atoms[j].Charge;
                posB=Cations[CurrentSystem][J].Atoms[j].Position;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeChargeSquared[CurrentSystem])
                {
                  PotentialSecondDerivativeCoulombic(chargeA,chargeB,rr,&energy,&DF,&DDF);

                  // add contribution to the energy
                  UCationCationChargeChargeReal[CurrentSystem]+=energy;

                  f.x=DF*dr.x;
                  f.y=DF*dr.y;
                  f.z=DF*dr.z;

                  // add contribution to the first derivatives
                  Cations[CurrentSystem][I].Atoms[i].Force.x-=f.x;
                  Cations[CurrentSystem][I].Atoms[i].Force.y-=f.y;
                  Cations[CurrentSystem][I].Atoms[i].Force.z-=f.z;

                  Cations[CurrentSystem][J].Atoms[j].Force.x+=f.x;
                  Cations[CurrentSystem][J].Atoms[j].Force.y+=f.y;
                  Cations[CurrentSystem][J].Atoms[j].Force.z+=f.z;

                  // add contribution to the strain derivative tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;

                  VECTOR dJI;

                  dJI.x=(posB.x-comB.x)-(posA.x-comA.x);
                  dJI.y=(posB.y-comB.y)-(posA.y-comA.y);
                  dJI.z=(posB.z-comB.z)-(posA.z-comA.z);

                  // the corrected com-com distance vector
                  drJI.x=dr.x+(posB.x-comB.x)-(posA.x-comA.x);
                  drJI.y=dr.y+(posB.y-comB.y)-(posA.y-comA.y);
                  drJI.z=dr.z+(posB.z-comB.z)-(posA.z-comA.z);

                  // add the contributions to the Born-term
                  BornTerm[CurrentSystem].xxxx+=DDF*drJI.x*dr.x*dr.x*dr.x+2.0*DF*drJI.x*dr.x+DDF*drJI.x*dr.x*dJI.x*dr.x+DF*drJI.x*dJI.x;
                  BornTerm[CurrentSystem].xxyy+=DDF*drJI.x*dr.x*dr.y*dr.y+0.5*DDF*drJI.x*dr.x*(dJI.y*dr.y+dJI.y*dr.y);
                  BornTerm[CurrentSystem].xxzz+=DDF*drJI.x*dr.x*dr.z*dr.z+0.5*DDF*drJI.x*dr.x*(dJI.z*dr.z+dJI.z*dr.z);
                  BornTerm[CurrentSystem].xxyz+=DDF*drJI.x*dr.x*dr.y*dr.z+0.5*DDF*drJI.x*dr.x*(dJI.y*dr.z+dJI.z*dr.y);
                  BornTerm[CurrentSystem].xxzx+=DDF*drJI.x*dr.x*dr.x*dr.z+DF*(drJI.x*dr.z)+0.5*DDF*drJI.x*dr.x*(dJI.x*dr.z+dJI.z*dr.x)+0.5*DF*(drJI.x*dJI.z);
                  BornTerm[CurrentSystem].xxxy+=DDF*drJI.x*dr.x*dr.x*dr.y+DF*drJI.x*dr.y+0.5*DDF*drJI.x*dr.x*(dJI.x*dr.y+dJI.y*dr.x)+0.5*DF*(drJI.x*dJI.y);

                  BornTerm[CurrentSystem].yyyy+=DDF*drJI.y*dr.y*dr.y*dr.y+2.0*DF*drJI.y*dr.y+DDF*drJI.y*dr.y*dJI.y*dr.y+DF*drJI.y*dJI.y;
                  BornTerm[CurrentSystem].yyzz+=0.5*DDF*(drJI.y*dr.y*dr.z*dr.z+dr.y*drJI.y*dr.z*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.z*dr.z+dJI.z*dr.z);
                  BornTerm[CurrentSystem].yyyz+=0.5*DDF*(drJI.y*dr.y*dr.y*dr.z+dr.y*drJI.y*dr.y*dr.z)+DF*(drJI.y*dr.z)+0.5*DDF*drJI.y*dr.y*(dJI.y*dr.z+dJI.z*dr.y)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].yyzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.y+dr.x*drJI.z*dr.y*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.y+dJI.y*dr.y);
                  BornTerm[CurrentSystem].yyxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.y+dr.x*drJI.y*dr.y*dr.y)+DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.y+dJI.y*dr.y)+0.5*DF*(drJI.x*dJI.y);

                  BornTerm[CurrentSystem].zzzz+=DDF*drJI.z*dr.z*dr.z*dr.z+2.0*DF*drJI.z*dr.z+DDF*drJI.z*dr.z*dJI.z*dr.z+DF*drJI.z*dJI.z;
                  BornTerm[CurrentSystem].zzyz+=0.5*DDF*(drJI.y*dr.z*dr.z*dr.z+dr.y*drJI.z*dr.z*dr.z)+DF*(drJI.y*dr.z)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.y*dJI.z);
                  BornTerm[CurrentSystem].zzzx+=0.5*DDF*(drJI.x*dr.z*dr.z*dr.z+dr.x*drJI.z*dr.z*dr.z)+DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.z*dr.z+dJI.z*dr.z)+0.5*DF*(drJI.x*dJI.z);
                  BornTerm[CurrentSystem].zzxy+=0.5*DDF*(drJI.x*dr.y*dr.z*dr.z+dr.x*drJI.y*dr.z*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.z*dr.z+dJI.z*dr.z);

                  BornTerm[CurrentSystem].yzyz+=0.5*DDF*(drJI.y*dr.z*dr.y*dr.z+dr.y*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.y*dr.y)+0.25*DDF*(drJI.y*dr.z+drJI.z*dr.y)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.z*dJI.z+drJI.y*dJI.y);
                  BornTerm[CurrentSystem].yzzx+=0.5*DDF*(drJI.x*dr.z*dr.y*dr.z+dr.x*drJI.z*dr.y*dr.z)+0.5*DF*(drJI.x*dr.y)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.y);
                  BornTerm[CurrentSystem].yzxy+=0.5*DDF*(drJI.x*dr.y*dr.y*dr.z+dr.x*drJI.y*dr.y*dr.z)+0.5*DF*(drJI.x*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.y*dr.z+dJI.z*dr.y)+0.25*DF*(drJI.x*dJI.z);

                  BornTerm[CurrentSystem].zxzx+=0.5*DDF*(drJI.x*dr.z*dr.x*dr.z+dr.x*drJI.z*dr.x*dr.z)+0.5*DF*(drJI.z*dr.z+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.z+drJI.z*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.z*dJI.z+drJI.x*dJI.x);
                  BornTerm[CurrentSystem].zxxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.z+dr.x*drJI.y*dr.x*dr.z)+0.5*DF*(drJI.y*dr.z)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.z+dJI.z*dr.x)+0.25*DF*(drJI.y*dJI.z);

                  BornTerm[CurrentSystem].xyxy+=0.5*DDF*(drJI.x*dr.y*dr.x*dr.y+dr.x*drJI.y*dr.x*dr.y)+0.5*DF*(drJI.y*dr.y+drJI.x*dr.x)+0.25*DDF*(drJI.x*dr.y+drJI.y*dr.x)*(dJI.x*dr.y+dJI.y*dr.x)+0.25*DF*(drJI.y*dJI.y+drJI.x*dJI.x);
                }
              }
            }
          }
        }
      }
    }
  }
}


