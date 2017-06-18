/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'recrossing.c' is part of RASPA-2.0

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
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <sys/stat.h>
#include "simulation.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "utils.h"
#include "molecule.h"
#include "input.h"
#include "output.h"
#include "mc_moves.h"
#include "statistics.h"
#include "potentials.h"
#include "cbmc.h"
#include "grids.h"
#include "thermo_baro_stats.h"
#include "integration.h"
#include "recrossing.h"
#include "sample.h"
#include "movies.h"
#include "rigid.h"

// BarrierRecrossing calculates the transmission coefficient. The transmission coefficient
// is the fraction of particles that start on top of the barrier and end up in the cage A
// when the initial velocity was directed to A. If the particle ends up in B it has recrossed
// the barrier and this will lower the transmission coefficient.

// This perticular implementation uses two particles (which do not notice eachother) with opposite
// velocities. Only if the particles end up in different cages do we get a contribution (the other
// two possibilities cancel exactly).

// As input the routine needs
// (1) a file "window_confs.dat" with center-of-mass coordinates of the particle on top of the barrier.
// (2) a windows-center point (the location of the barrier).
// (3) a vector describing the direction of reaction coordinate. This can often be chosen
//     as a vector normal to the cage window.

// As output the routine gives a file "kt.data" containing
//   column 1: the reaction coordinate
//   column 2: the transmission coeficient
//   column 3: the error of the transmission coefficient

extern bool STREAM;

VECTOR *BarrierPosition;
VECTOR *BarrierNormal;
VECTOR *BarrierAngle;
REAL *MaxBarrierDistance;
REAL *MaxBarrierTime;
int *NumberOfVelocities;
int *PutMoleculeOnBarrier;


int ReadMoleculePosition(FILE *fileptr)
{
  int i,j,l;
  VECTOR s;
  FILE *FilePtr;

  FilePtr=fileptr;

  fscanf(FilePtr,"#NumberOfAdsorbateMolecules:%d\n",&NumberOfAdsorbateMolecules[CurrentSystem]);
  fscanf(FilePtr,"#NumberOfCationMolecules:%d\n",&NumberOfCationMolecules[CurrentSystem]);
  fscanf(FilePtr,"#NumberOfFrameworkAtoms:%d\n",&Framework[CurrentSystem].TotalNumberOfAtoms);
  fscanf(FilePtr,"#Barrier:%lf %lf %lf\n",&s.x,&s.y,&s.z);

  // barrier position in Cartesian coordinates
  BarrierPosition[CurrentSystem].x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
  BarrierPosition[CurrentSystem].y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
  BarrierPosition[CurrentSystem].z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      if(fscanf(FilePtr,"%lf %lf %lf %lf %lf %lf\n",
        &Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition.x,&Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition.y,
        &Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition.z,&Adsorbates[CurrentSystem][i].Atoms[j].ReferenceVelocity.x,
        &Adsorbates[CurrentSystem][i].Atoms[j].ReferenceVelocity.y,&Adsorbates[CurrentSystem][i].Atoms[j].ReferenceVelocity.z)==EOF) return EOF;
      fprintf(stderr, "%g %g %g\n",
         Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition.x,
         Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition.y,
         Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition.z);
    }
  }
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      if(fscanf(FilePtr,"%lf %lf %lf %lf %lf %lf\n",
        &Cations[CurrentSystem][i].Atoms[j].ReferencePosition.x,&Cations[CurrentSystem][i].Atoms[j].ReferencePosition.y,
        &Cations[CurrentSystem][i].Atoms[j].ReferencePosition.z,&Cations[CurrentSystem][i].Atoms[j].ReferenceVelocity.x,
        &Cations[CurrentSystem][i].Atoms[j].ReferenceVelocity.y,&Cations[CurrentSystem][i].Atoms[j].ReferenceVelocity.z)==EOF) return EOF;
    }
  }
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(l=0;l<Framework[CurrentSystem].NumberOfFrameworks;l++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[l];i++)
      {
        if(fscanf(FilePtr,"%lf %lf %lf %lf %lf %lf\n",
          &Framework[CurrentSystem].Atoms[l][i].ReferencePosition.x,
          &Framework[CurrentSystem].Atoms[l][i].ReferencePosition.y,
          &Framework[CurrentSystem].Atoms[l][i].ReferencePosition.z,
          &Framework[CurrentSystem].Atoms[l][i].ReferenceVelocity.x,
          &Framework[CurrentSystem].Atoms[l][i].ReferenceVelocity.y,
          &Framework[CurrentSystem].Atoms[l][i].ReferenceVelocity.z)==EOF) return EOF;
      }
    }
  }
  return 0;
}

  // the unit-length barrier in fractional coordinates
void SetBarierNormal(void)
{
  VECTOR s;

  switch(FreeEnergyMappingType[CurrentSystem])
  {
    case A_MAPPING:
      s.x=1.0;
      s.y=0.0;
      s.z=0.0;
      break;
    case B_MAPPING:
      s.x=0.0;
      s.y=1.0;
      s.z=0.0;
      break;
    case C_MAPPING:
      s.x=0.0;
      s.y=0.0;
      s.z=1.0;
      break;
    case MAP_AB_DIAGONAL:
      s.x=M_SQRT1_2;
      s.y=-M_SQRT1_2;
      s.z=0.0;
      break;
    case MAP_AC_DIAGONAL:
      s.x=M_SQRT1_2;
      s.y=0.0;
      s.z=-M_SQRT1_2;
      break;
    case MAP_BC_DIAGONAL:
      s.x=0.0;
      s.y=M_SQRT1_2;
      s.z=-M_SQRT1_2;
    case MAP_O_AB_DIAGONAL:
      s.x=M_SQRT1_2;
      s.y=M_SQRT1_2;
      s.z=0.0;
      break;
    case MAP_O_AC_DIAGONAL:
      s.x=M_SQRT1_2;
      s.y=0.0;
      s.z=M_SQRT1_2;
      break;
    case MAP_O_BC_DIAGONAL:
      s.x=0.0;
      s.y=M_SQRT1_2;
      s.z=M_SQRT1_2;
      break;
    case MAP_A_BC_DIAGONAL:
      s.x=M_SQRT1_3;
      s.y=-M_SQRT1_3;
      s.z=-M_SQRT1_3;
      break;
    case MAP_B_AC_DIAGONAL:
      s.x=-M_SQRT1_3;
      s.y=M_SQRT1_3;
      s.z=-M_SQRT1_3;
      break;
    case MAP_C_AB_DIAGONAL:
      s.x=-M_SQRT1_3;
      s.y=-M_SQRT1_3;
      s.z=M_SQRT1_3;
      break;
    case MAP_O_ABC_DIAGONAL:
      s.x=-M_SQRT1_3;
      s.y=-M_SQRT1_3;
      s.z=-M_SQRT1_3;
      break;
    default:
      s.x=0.0;
      s.y=0.0;
      s.z=0.0;
      break;
  }

  BarrierNormal[CurrentSystem]=s;
}

void ReverseVelocities(void)
{
  int i,j,l,A,Type;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x=-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y=-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z=-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z;
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.x=-Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.x;
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.y=-Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.y;
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.z=-Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.z;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum=AngularVelocityToQuaternionMomentumAdsorbates(i,l);
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x=-Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y=-Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z=-Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z;
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x=-Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y=-Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z=-Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z;
        Cations[CurrentSystem][i].Groups[l].AngularVelocity.x=-Cations[CurrentSystem][i].Groups[l].AngularVelocity.x;
        Cations[CurrentSystem][i].Groups[l].AngularVelocity.y=-Cations[CurrentSystem][i].Groups[l].AngularVelocity.y;
        Cations[CurrentSystem][i].Groups[l].AngularVelocity.z=-Cations[CurrentSystem][i].Groups[l].AngularVelocity.z;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum=AngularVelocityToQuaternionMomentumCations(i,l);
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Cations[CurrentSystem][i].Atoms[A].Velocity.x=-Cations[CurrentSystem][i].Atoms[A].Velocity.x;
          Cations[CurrentSystem][i].Atoms[A].Velocity.y=-Cations[CurrentSystem][i].Atoms[A].Velocity.y;
          Cations[CurrentSystem][i].Atoms[A].Velocity.z=-Cations[CurrentSystem][i].Atoms[A].Velocity.z;
        }
      }
    }
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(l=0;l<Framework[CurrentSystem].NumberOfFrameworks;l++)
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[l];i++)
      {
        Framework[CurrentSystem].Atoms[l][i].Velocity.x=-Framework[CurrentSystem].Atoms[l][i].Velocity.x;
        Framework[CurrentSystem].Atoms[l][i].Velocity.y=-Framework[CurrentSystem].Atoms[l][i].Velocity.y;
        Framework[CurrentSystem].Atoms[l][i].Velocity.z=-Framework[CurrentSystem].Atoms[l][i].Velocity.z;
      }
  }

  CreateCartesianVelocities();
}

void StorePositionsAndVelocities(void)
{
  int i,j,l,A,Type;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassReferenceVelocity=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
        Adsorbates[CurrentSystem][i].Groups[l].AngularReferenceVelocity=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum=Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassReferencePosition=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition;
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Adsorbates[CurrentSystem][i].Atoms[A].ReferenceVelocity=Adsorbates[CurrentSystem][i].Atoms[A].Velocity;
          Adsorbates[CurrentSystem][i].Atoms[A].ReferencePosition=Adsorbates[CurrentSystem][i].Atoms[A].Position;
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Cations[CurrentSystem][i].Groups[l].CenterOfMassReferenceVelocity=Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
        Cations[CurrentSystem][i].Groups[l].AngularReferenceVelocity=Cations[CurrentSystem][i].Groups[l].AngularVelocity;
        Cations[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum=Cations[CurrentSystem][i].Groups[l].QuaternionMomentum;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassReferencePosition=Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition;
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Cations[CurrentSystem][i].Atoms[A].ReferenceVelocity=Cations[CurrentSystem][i].Atoms[A].Velocity;
          Cations[CurrentSystem][i].Atoms[A].ReferencePosition=Cations[CurrentSystem][i].Atoms[A].Position;
        }
      }
    }
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(l=0;l<Framework[CurrentSystem].NumberOfFrameworks;l++)
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[l];i++)
      {
        Framework[CurrentSystem].Atoms[l][i].ReferenceVelocity=Framework[CurrentSystem].Atoms[l][i].Velocity;
        Framework[CurrentSystem].Atoms[l][i].ReferencePosition=Framework[CurrentSystem].Atoms[l][i].Position;
      }
  }
}

// retrieve velocities and invert the directions
void RetrievePositionsAndVelocitiesAndInvertVelocities(void)
{
  int i,j,l,A,Type;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassReferenceVelocity;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum=Adsorbates[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassReferencePosition;

        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x=-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassReferenceVelocity.x;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y=-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassReferenceVelocity.y;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z=-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassReferenceVelocity.z;
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.x=-Adsorbates[CurrentSystem][i].Groups[l].AngularReferenceVelocity.x;
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.y=-Adsorbates[CurrentSystem][i].Groups[l].AngularReferenceVelocity.y;
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.z=-Adsorbates[CurrentSystem][i].Groups[l].AngularReferenceVelocity.z;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.r=-Adsorbates[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum.r;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.i=-Adsorbates[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum.i;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.j=-Adsorbates[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum.j;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.k=-Adsorbates[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum.k;
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Adsorbates[CurrentSystem][i].Atoms[A].Position=Adsorbates[CurrentSystem][i].Atoms[A].ReferencePosition;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x=-Adsorbates[CurrentSystem][i].Atoms[A].ReferenceVelocity.x;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y=-Adsorbates[CurrentSystem][i].Atoms[A].ReferenceVelocity.y;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z=-Adsorbates[CurrentSystem][i].Atoms[A].ReferenceVelocity.z;
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity=Cations[CurrentSystem][i].Groups[l].CenterOfMassReferenceVelocity;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum=Cations[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition=Cations[CurrentSystem][i].Groups[l].CenterOfMassReferencePosition;

        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x=-Cations[CurrentSystem][i].Groups[l].CenterOfMassReferenceVelocity.x;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y=-Cations[CurrentSystem][i].Groups[l].CenterOfMassReferenceVelocity.y;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z=-Cations[CurrentSystem][i].Groups[l].CenterOfMassReferenceVelocity.z;
        Cations[CurrentSystem][i].Groups[l].AngularVelocity.x=-Cations[CurrentSystem][i].Groups[l].AngularReferenceVelocity.x;
        Cations[CurrentSystem][i].Groups[l].AngularVelocity.y=-Cations[CurrentSystem][i].Groups[l].AngularReferenceVelocity.y;
        Cations[CurrentSystem][i].Groups[l].AngularVelocity.z=-Cations[CurrentSystem][i].Groups[l].AngularReferenceVelocity.z;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.r=-Cations[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum.r;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.i=-Cations[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum.i;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.j=-Cations[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum.j;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.k=-Cations[CurrentSystem][i].Groups[l].QuaternionReferenceMomentum.k;
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Cations[CurrentSystem][i].Atoms[A].Position=Cations[CurrentSystem][i].Atoms[A].ReferencePosition;
          Cations[CurrentSystem][i].Atoms[A].Velocity.x=-Cations[CurrentSystem][i].Atoms[A].ReferenceVelocity.x;
          Cations[CurrentSystem][i].Atoms[A].Velocity.y=-Cations[CurrentSystem][i].Atoms[A].ReferenceVelocity.y;
          Cations[CurrentSystem][i].Atoms[A].Velocity.z=-Cations[CurrentSystem][i].Atoms[A].ReferenceVelocity.z;
        }
      }
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    UpdateGroupCenterOfMassAdsorbate(i);
    ComputeQuaternionAdsorbate(i);
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    UpdateGroupCenterOfMassCation(i);
    ComputeQuaternionCation(i);
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(l=0;l<Framework[CurrentSystem].NumberOfFrameworks;l++)
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[l];i++)
      {
        Framework[CurrentSystem].Atoms[l][i].Position=Framework[CurrentSystem].Atoms[l][i].ReferencePosition;
        Framework[CurrentSystem].Atoms[l][i].Velocity.x=-Framework[CurrentSystem].Atoms[l][i].ReferenceVelocity.x;
        Framework[CurrentSystem].Atoms[l][i].Velocity.y=-Framework[CurrentSystem].Atoms[l][i].ReferenceVelocity.y;
        Framework[CurrentSystem].Atoms[l][i].Velocity.z=-Framework[CurrentSystem].Atoms[l][i].ReferenceVelocity.z;
      }
  }
}

int BarrierRecrossing(void)
{
  int i,j,k,count[CurrentSystem],index,*iterations;
  static FILE **FilePtr;
  FILE *FilePtrOut;
  VECTOR ra,rb,velA;
  REAL DistanceAToPlaneABC,DistanceBToPlaneABC;
  REAL DistanceAToPlaneXYZ,DistanceBToPlaneXYZ;
  REAL U0=0.0,drift;
  int ReactionBead;
  char buffer[1024];
  REAL **Kt;
  REAL Qdot,*QdotAv;
  REAL TempStartA,TempStartB,TempEndA,TempEndB;
  int stop,*Status;

  if (STREAM)
  {
    fprintf(stderr, "Streaming not implemented for BarrierRecrossing.");
    exit(0);
  }

  iterations=(int*)calloc(NumberOfSystems,sizeof(int));
  Status=(int*)calloc(NumberOfSystems,sizeof(int));
  QdotAv=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  FilePtr=(FILE**)calloc(NumberOfSystems,sizeof(FILE*));
  Kt=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  for(i=0;i<NumberOfSystems;i++)
  {
    iterations[i]=(int)(MaxBarrierTime[i]/DeltaT);
    Kt[i]=(REAL*)calloc(iterations[i],sizeof(REAL));
  }

  ReactionBead=Components[0].StartingBead;

  InitializeNoseHooverAllSystems();

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    SetBarierNormal();

    fprintf(OutputFilePtr[CurrentSystem],"barrier time is restricted to %lf [ps]\n",(double)((REAL)iterations[CurrentSystem]*DeltaT));
    fprintf(OutputFilePtr[CurrentSystem],"window: %lf %lf %lf\n",(double)BarrierPosition[CurrentSystem].x,
              (double)BarrierPosition[CurrentSystem].y,(double)BarrierPosition[CurrentSystem].z);
    fprintf(OutputFilePtr[CurrentSystem],"normal: %lf %lf %lf\n",(double)BarrierNormal[CurrentSystem].x,
              (double)BarrierNormal[CurrentSystem].y,(double)BarrierNormal[CurrentSystem].z);
    fprintf(OutputFilePtr[CurrentSystem],"angle: %lf %lf %lf\n",(double)BarrierAngle[CurrentSystem].x,
              (double)BarrierAngle[CurrentSystem].y,(double)BarrierAngle[CurrentSystem].z);

    sprintf(buffer,"dcTST_starting_configurations/System_%d/configurations.data",CurrentSystem);

    if(!(FilePtr[CurrentSystem]=fopen(buffer,"r")))
    {
      fprintf(stderr, "Error opening input-file 'configurations.data'\n");
      exit(1);
    }

    count[CurrentSystem]=0;
    QdotAv[CurrentSystem]=0.0;
    for(i=0;i<iterations[CurrentSystem];i++)
      Kt[CurrentSystem][i]=0.0;
  }

  for(;;)
  {
    // read the new configuration from files
    stop=TRUE;
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      Status[CurrentSystem]=ReadMoleculePosition(FilePtr[CurrentSystem]);
      stop&=(Status[CurrentSystem]==EOF);
    }

    // stop when all configurations from all files for all systems are read
    if(stop) break;

    // for all systems handle the read configuration
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      // continue as long as we can read configurations
      if(Status[CurrentSystem]!=EOF)
      {
        // print the reaction position in Cartesian coordinates
        fprintf(OutputFilePtr[CurrentSystem],"New configuration: %g %g %g\n",
            Adsorbates[CurrentSystem][0].Atoms[ReactionBead].ReferencePosition.x,
            Adsorbates[CurrentSystem][0].Atoms[ReactionBead].ReferencePosition.y,
            Adsorbates[CurrentSystem][0].Atoms[ReactionBead].ReferencePosition.z);

        for(k=0;k<NumberOfVelocities[CurrentSystem];k++)
        {
          // store the positions of the adsorbates to the reference positions
          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
            for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
              Adsorbates[CurrentSystem][i].Atoms[j].Position=Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition;

          // store the positions of the cations to the reference positions
          for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
            for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
              Cations[CurrentSystem][i].Atoms[j].Position=Cations[CurrentSystem][i].Atoms[j].ReferencePosition;

          // for a flexible framework: store the positions of the frameworks to the reference positions
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
              for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
                 Framework[CurrentSystem].Atoms[CurrentFramework][i].Position=
                   Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferencePosition;
          }

          // determine the center-of-mass and the quaternions for the adsorbates
          for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
          {
            UpdateGroupCenterOfMassAdsorbate(i);
            ComputeQuaternionAdsorbate(i);
          }

          // determine the center-of-mass and the quaternions for the cations
          for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
          {
            UpdateGroupCenterOfMassCation(i);
            ComputeQuaternionCation(i);
          }

          // initialize the velocity (linear and angular velocity drawn from a Maxwell-Boltzmann distribution)
          InitializeAdsorbateVelocities();
          InitializeCationVelocities();
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            InitializeFrameworkVelocities();

          // for rigid units: create atomic velocities from com-velocity and quaternion-velocity
          CreateCartesianVelocities();

          // compute initial forces
          InitializeForces();

          fprintf(OutputFilePtr[CurrentSystem],"Potential energy at start: %lf [K]\n",UTotal[CurrentSystem]*ENERGY_TO_KELVIN);

          // Force particle 1 to the "right" based on the atomic velocity of the reaction-bead
          velA=Adsorbates[CurrentSystem][0].Atoms[ReactionBead].Velocity;
          if(DotProduct(ConvertFromXYZtoABCUnitCell(velA),BarrierNormal[CurrentSystem])<=0.0)
            ReverseVelocities();

          StorePositionsAndVelocities();

          velA=Adsorbates[CurrentSystem][0].Atoms[ReactionBead].Velocity;
          Qdot=DotProduct(ConvertFromXYZtoABCUnitCell(velA),BarrierNormal[CurrentSystem]);
          QdotAv[CurrentSystem]+=fabs(Qdot);

          // Integrate A
          // ===========


          //InitializeForces();
          InitializeForces();
          TempStartA=2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem]);

          drift=0.0;
          DistanceAToPlaneABC=0.0;
          for(index=0;index<iterations[CurrentSystem];index++)
          {
            // compute the distance from the barrier
            ra.x=Adsorbates[CurrentSystem][0].Atoms[ReactionBead].Position.x-BarrierPosition[CurrentSystem].x;
            ra.y=Adsorbates[CurrentSystem][0].Atoms[ReactionBead].Position.y-BarrierPosition[CurrentSystem].y;
            ra.z=Adsorbates[CurrentSystem][0].Atoms[ReactionBead].Position.z-BarrierPosition[CurrentSystem].z;
            DistanceAToPlaneABC=DotProduct(ConvertFromXYZtoABCUnitCell(ra),BarrierNormal[CurrentSystem]);
            ra.x=BarrierNormal[CurrentSystem].x*DistanceAToPlaneABC;
            ra.y=BarrierNormal[CurrentSystem].y*DistanceAToPlaneABC;
            ra.z=BarrierNormal[CurrentSystem].z*DistanceAToPlaneABC;
            DistanceAToPlaneXYZ=DotProduct(ConvertFromABCtoXYZUnitCell(ra),BarrierNormal[CurrentSystem]);


            // only two contributions, the other two cancel by construction
            if((index==0)||(DistanceAToPlaneXYZ>0.0))
              Kt[CurrentSystem][index]+=0.5*Qdot;
            else
              Kt[CurrentSystem][index]-=0.5*Qdot;

            // keep track of the drift in the energy during integration
            if(index==0)
              U0=ConservedEnergy[CurrentSystem];
            else
              drift+=fabs(ConservedEnergy[CurrentSystem]-U0)/U0;

            // integrate as long as the distance is within the specified barrier-distance
            if(fabs(DistanceAToPlaneXYZ)<MaxBarrierDistance[CurrentSystem])
              Integration();
          }

          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            UKinetic[CurrentSystem]=GetFrameworkKineticEnergy()+GetAdsorbateKineticEnergy()+GetCationKineticEnergy();
          else
            UKinetic[CurrentSystem]=GetAdsorbateKineticEnergy()+GetCationKineticEnergy();

          TempEndA=2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem]);

          fprintf(OutputFilePtr[CurrentSystem],"Energy Conservation trajectory A: %lf %lf drift: %lf distance to plane: %lf\n",
              (double)U0,(double)ConservedEnergy[CurrentSystem],(double)(drift/(REAL)iterations[CurrentSystem]),(double)DistanceAToPlaneXYZ);
          fflush(OutputFilePtr[CurrentSystem]);

          // Integrate B
          // ===========

          RetrievePositionsAndVelocitiesAndInvertVelocities();

          InitializeForces();
          TempStartB=2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem]);

          drift=0.0;
          DistanceBToPlaneABC=0.0;
          for(index=0;index<iterations[CurrentSystem];index++)
          {
            // compute the distance from the barrier
            rb.x=Adsorbates[CurrentSystem][0].Atoms[ReactionBead].Position.x-BarrierPosition[CurrentSystem].x;
            rb.y=Adsorbates[CurrentSystem][0].Atoms[ReactionBead].Position.y-BarrierPosition[CurrentSystem].y;
            rb.z=Adsorbates[CurrentSystem][0].Atoms[ReactionBead].Position.z-BarrierPosition[CurrentSystem].z;
            DistanceBToPlaneABC=DotProduct(ConvertFromXYZtoABCUnitCell(rb),BarrierNormal[CurrentSystem]);
            rb.x=BarrierNormal[CurrentSystem].x*DistanceBToPlaneABC;
            rb.y=BarrierNormal[CurrentSystem].y*DistanceBToPlaneABC;
            rb.z=BarrierNormal[CurrentSystem].z*DistanceBToPlaneABC;
            DistanceBToPlaneXYZ=DotProduct(ConvertFromABCtoXYZUnitCell(rb),BarrierNormal[CurrentSystem]);


            // only two contributions, the other two cancel by construction
            if((index==0)||DistanceBToPlaneXYZ<=0.0)
              Kt[CurrentSystem][index]+=0.5*Qdot;
            else
              Kt[CurrentSystem][index]-=0.5*Qdot;

            // keep track of the drift in the energy during integration
            if(index==0)
              U0=ConservedEnergy[CurrentSystem];
            else
              drift+=fabs(ConservedEnergy[CurrentSystem]-U0)/U0;

            // integrate as long as the distance is within the specified barrier-distance
            if(fabs(DistanceBToPlaneXYZ)<MaxBarrierDistance[CurrentSystem])
              Integration();
          }

          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            UKinetic[CurrentSystem]=GetFrameworkKineticEnergy()+GetAdsorbateKineticEnergy()+GetCationKineticEnergy();
          else
            UKinetic[CurrentSystem]=GetAdsorbateKineticEnergy()+GetCationKineticEnergy();

          TempEndB=2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem]);

          fprintf(OutputFilePtr[CurrentSystem],"Energy Conservation trajectory B: %lf %lf drift: %lf distance to plane: %lf\n",
              (double)U0,(double)ConservedEnergy[CurrentSystem],(double)(drift/(REAL)iterations[CurrentSystem]),(double)DistanceBToPlaneXYZ);
          fflush(OutputFilePtr[CurrentSystem]);

          count[CurrentSystem]++;
        }
        fprintf(OutputFilePtr[CurrentSystem],"\n");

        // write output data to file
        mkdir("TransmissionCoefficient",S_IRWXU);
        sprintf(buffer,"TransmissionCoefficient/System_%d",CurrentSystem);
        mkdir(buffer,S_IRWXU);
        sprintf(buffer,"TransmissionCoefficient/System_%d/kt.data",CurrentSystem);
        FilePtrOut=fopen(buffer,"w");
        fprintf(FilePtrOut,"#iteration %d\n",count[CurrentSystem]);
        for(i=0;i<iterations[CurrentSystem];i++)
          fprintf(FilePtrOut,"%d %lf %lf\n",i,(double)((REAL)i*DeltaT),(double)(Kt[CurrentSystem][i]/(QdotAv[CurrentSystem])));
        fclose(FilePtrOut);
      }
    }
  }

  // close the reading files
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    fclose(FilePtr[CurrentSystem]);
  return 0;
}

void AllocateRecrossingMemory(void)
{
  BarrierPosition=(VECTOR*)calloc(NumberOfSystems,sizeof(VECTOR));
  BarrierNormal=(VECTOR*)calloc(NumberOfSystems,sizeof(VECTOR));
  BarrierAngle=(VECTOR*)calloc(NumberOfSystems,sizeof(VECTOR));

  MaxBarrierDistance=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  MaxBarrierTime=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  NumberOfVelocities=(int*)calloc(NumberOfSystems,sizeof(int));
  PutMoleculeOnBarrier=(int*)calloc(NumberOfSystems,sizeof(int));
}

