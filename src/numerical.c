/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'numerical.c' is part of RASPA-2.0

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
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "constants.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "framework_born.h"
#include "framework_hessian.h"
#include "inter_force.h"
#include "internal_force.h"
#include "internal_hessian.h"
#include "internal_born.h"
#include "inter_hessian.h"
#include "molecule.h"
#include "potentials.h"
#include "numerical.h"
#include "integration.h"
#include "matrix.h"
#include "thermo_baro_stats.h"
#include "ewald.h"
#include "simulation.h"
#include "output.h"
#include "spectra.h"
#include "cubic_spline_1d.h"
#include "minimization.h"
#include "rigid.h"
#include "mc_moves.h"

extern bool STREAM;

void NumericallyComputeGradient(int np,int nb,REAL *x,REAL *NumericalGradient)
{
  int i,j;
  REAL store,delta,Energy;
  REAL EnergyForward1,EnergyBackward1;
  REAL EnergyForward2,EnergyBackward2;
  REAL *Gradient=NULL;
  REAL_MATRIX Hessian;
  REAL_MATRIX3x3 StrainDerivative;
  REAL pressure,ExtPressure;

  CurrentSystem=0;
  BoundaryCondition[0]=TRICLINIC;

  delta=1e-7;
  for(i=0;i<np+nb;i++)
  {
    store=x[i];

    EnergyBackward2=0.0;
    x[i]=store-delta;
    CreatePositionsFromGeneralizedCoordinates(np,nb,x);
    EvaluateDerivatives(np,&EnergyBackward2,Gradient,Hessian,&StrainDerivative,FALSE,FALSE);

    EnergyBackward1=0.0;
    x[i]=store-0.5*delta;
    CreatePositionsFromGeneralizedCoordinates(np,nb,x);
    EvaluateDerivatives(np,&EnergyBackward1,Gradient,Hessian,&StrainDerivative,FALSE,FALSE);

    EnergyForward1=0.0;
    x[i]=store+0.5*delta;
    CreatePositionsFromGeneralizedCoordinates(np,nb,x);
    EvaluateDerivatives(np,&EnergyForward1,Gradient,Hessian,&StrainDerivative,FALSE,FALSE);

    EnergyForward2=0.0;
    x[i]=store+delta;
    CreatePositionsFromGeneralizedCoordinates(np,nb,x);
    EvaluateDerivatives(np,&EnergyForward2,Gradient,Hessian,&StrainDerivative,FALSE,FALSE);

    x[i]=store;
    CreatePositionsFromGeneralizedCoordinates(np,nb,x);
    EvaluateDerivatives(np,&Energy,Gradient,Hessian,&StrainDerivative,FALSE,FALSE);

    NumericalGradient[i]=(EnergyBackward2-8.0*EnergyBackward1+8.0*EnergyForward1-EnergyForward2)/(6.0*delta);
  }

  // correction for tail-correction
  pressure=0.0;
  for(i=0;i<NumberOfPseudoAtoms;i++)
    for(j=0;j<NumberOfPseudoAtoms;j++)
    {
      if(TailCorrection[i][j])
        pressure+=(2.0/3.0)*M_PI*CUBE(CutOffVDW)*(REAL)NumberOfPseudoAtomsType[CurrentSystem][i]*(REAL)NumberOfPseudoAtomsType[CurrentSystem][j]*
                  PotentialValue(i,j,CutOffVDWSquared,1.0);
    }

  // correction for fixed external pressure
  ExtPressure=therm_baro_stats.ExternalPressure[CurrentSystem][0];

  switch(Ensemble[CurrentSystem])
  {
    case NVE:
    case NVT:
      break;
    case NPT:
    case NPH:
      NumericalGradient[np]-=3.0*pressure/Volume[CurrentSystem];
      NumericalGradient[np]+=ExtPressure*3.0*Volume[CurrentSystem];
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          break;
        case ANISOTROPIC:
          switch(Dimension)
          {
            case 3:
              NumericalGradient[np+2]-=pressure/Volume[CurrentSystem];
              NumericalGradient[np+2]+=ExtPressure*Volume[CurrentSystem];
            case 2:
              NumericalGradient[np+1]-=pressure/Volume[CurrentSystem];
              NumericalGradient[np+1]+=ExtPressure*Volume[CurrentSystem];
            case 1:
              NumericalGradient[np]-=pressure/Volume[CurrentSystem];
              NumericalGradient[np]+=ExtPressure*Volume[CurrentSystem];
              break;
          }
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              NumericalGradient[np]-=pressure/Volume[CurrentSystem];
              NumericalGradient[np+1]-=pressure/Volume[CurrentSystem];
              NumericalGradient[np+3]-=pressure/Volume[CurrentSystem];

              NumericalGradient[np]+=ExtPressure*Volume[CurrentSystem];
              NumericalGradient[np+1]+=ExtPressure*Volume[CurrentSystem];
              NumericalGradient[np+3]+=ExtPressure*Volume[CurrentSystem];
              break;
            case MONOCLINIC_BETA_ANGLE:
              NumericalGradient[np]-=pressure/Volume[CurrentSystem];
              NumericalGradient[np+2]-=pressure/Volume[CurrentSystem];
              NumericalGradient[np+3]-=pressure/Volume[CurrentSystem];

              NumericalGradient[np]+=ExtPressure*Volume[CurrentSystem];
              NumericalGradient[np+2]+=ExtPressure*Volume[CurrentSystem];
              NumericalGradient[np+3]+=ExtPressure*Volume[CurrentSystem];
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              NumericalGradient[np]-=pressure/Volume[CurrentSystem];
              NumericalGradient[np+2]-=pressure/Volume[CurrentSystem];
              NumericalGradient[np+3]-=pressure/Volume[CurrentSystem];

              NumericalGradient[np]+=ExtPressure*Volume[CurrentSystem];
              NumericalGradient[np+2]+=ExtPressure*Volume[CurrentSystem];
              NumericalGradient[np+3]+=ExtPressure*Volume[CurrentSystem];
              break;
          }
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          switch(Dimension)
          {
            case 3:
             NumericalGradient[np+5]-=pressure/Volume[CurrentSystem];
            case 2:
              NumericalGradient[np+Dimension]-=pressure/Volume[CurrentSystem];
            case 1:
              NumericalGradient[np]-=pressure/Volume[CurrentSystem];
              break;
          }

            if(therm_baro_stats.UseExternalStress)
            {
              switch(Dimension)
              {
                case 3:
                  NumericalGradient[np]+=therm_baro_stats.ExternalStress[CurrentSystem].ax*Volume[CurrentSystem];
                  NumericalGradient[np+1]+=therm_baro_stats.ExternalStress[CurrentSystem].bx*Volume[CurrentSystem];
                  NumericalGradient[np+2]+=therm_baro_stats.ExternalStress[CurrentSystem].cx*Volume[CurrentSystem];
                  NumericalGradient[np+3]+=therm_baro_stats.ExternalStress[CurrentSystem].by*Volume[CurrentSystem];
                  NumericalGradient[np+4]+=therm_baro_stats.ExternalStress[CurrentSystem].cy*Volume[CurrentSystem];
                  NumericalGradient[np+5]+=therm_baro_stats.ExternalStress[CurrentSystem].cz*Volume[CurrentSystem];
                  break;
                case 2:
                  NumericalGradient[np]+=therm_baro_stats.ExternalStress[CurrentSystem].ax*Volume[CurrentSystem];
                  NumericalGradient[np+1]+=therm_baro_stats.ExternalStress[CurrentSystem].bx*Volume[CurrentSystem];
                  NumericalGradient[np+2]+=therm_baro_stats.ExternalStress[CurrentSystem].by*Volume[CurrentSystem];
                case 1:
                  NumericalGradient[np]+=therm_baro_stats.ExternalStress[CurrentSystem].ax*Volume[CurrentSystem];
                  break;
              }
            }
            else
            {
              switch(Dimension)
              {
                case 3:
                  NumericalGradient[np+5]+=ExtPressure*Volume[CurrentSystem];
                case 2:
                  NumericalGradient[np+Dimension]+=ExtPressure*Volume[CurrentSystem];
                case 1:
                  NumericalGradient[np]+=ExtPressure*Volume[CurrentSystem];
                  break;
              }
            }

          break;
      }
  }

}


void NumericallyComputeDerivatives(int np,int nb,REAL *x,REAL *GradientNumerical,REAL_MATRIX HessianNumerical)
{
  int i,j;
  REAL store,delta,Energy;
  REAL EnergyForward1,EnergyBackward1;
  REAL EnergyForward2,EnergyBackward2;
  REAL *GradientForward1,*GradientForward2;
  REAL *GradientBackward1,*GradientBackward2;
  REAL *Gradient;
  REAL_MATRIX Hessian;
  REAL_MATRIX3x3 StrainDerivative;

  delta=1e-4;
  Gradient=(REAL*)calloc(np+nb,sizeof(REAL));
  GradientForward1=(REAL*)calloc(np+nb,sizeof(REAL));
  GradientForward2=(REAL*)calloc(np+nb,sizeof(REAL));
  GradientBackward1=(REAL*)calloc(np+nb,sizeof(REAL));
  GradientBackward2=(REAL*)calloc(np+nb,sizeof(REAL));

  for(i=0;i<np+nb;i++)
  {
    store=x[i];

    EnergyForward1=0.0;
    for(j=0;j<np+nb;j++)
      GradientForward1[j]=0.0;
    x[i]=store+0.5*delta;

    CreatePositionsFromGeneralizedCoordinates(np,nb,x);
    EvaluateDerivatives(np,&EnergyForward1,GradientForward1,Hessian,&StrainDerivative,TRUE,FALSE);
    if(MinimizationVariables==FRACTIONAL)
      ConvertGradientFromCartesianToFractional(GradientForward1);

    EnergyForward2=0.0;
    for(j=0;j<np+nb;j++)
      GradientForward2[j]=0.0;
    x[i]=store+delta;
    CreatePositionsFromGeneralizedCoordinates(np,nb,x);
    EvaluateDerivatives(np,&EnergyForward2,GradientForward2,Hessian,&StrainDerivative,TRUE,FALSE);
    if(MinimizationVariables==FRACTIONAL)
      ConvertGradientFromCartesianToFractional(GradientForward2);

    EnergyBackward1=0.0;
    for(j=0;j<np+nb;j++)
      GradientBackward1[j]=0.0;
    x[i]=store-0.5*delta;
    CreatePositionsFromGeneralizedCoordinates(np,nb,x);
    EvaluateDerivatives(np,&EnergyBackward1,GradientBackward1,Hessian,&StrainDerivative,TRUE,FALSE);
    if(MinimizationVariables==FRACTIONAL)
      ConvertGradientFromCartesianToFractional(GradientBackward1);

    EnergyBackward2=0.0;
    for(j=0;j<np+nb;j++)
      GradientBackward2[j]=0.0;
    x[i]=store-delta;
    CreatePositionsFromGeneralizedCoordinates(np,nb,x);
    EvaluateDerivatives(np,&EnergyBackward2,GradientBackward2,Hessian,&StrainDerivative,TRUE,FALSE);
    if(MinimizationVariables==FRACTIONAL)
      ConvertGradientFromCartesianToFractional(GradientBackward2);

    x[i]=store;
    CreatePositionsFromGeneralizedCoordinates(np,nb,x);
    EvaluateDerivatives(np,&Energy,Gradient,Hessian,&StrainDerivative,TRUE,FALSE);

    for(j=0;j<np+nb;j++)
      HessianNumerical.element[i][j]=(GradientBackward2[j]-8.0*GradientBackward1[j]+8.0*GradientForward1[j]-GradientForward2[j])/(6.0*delta);
  }

  for(i=0;i<np+nb;i++)
    for(j=0;j<np+nb;j++)
      HessianNumerical.element[j][i]=HessianNumerical.element[i][j];

  free(GradientForward1);
  free(GradientForward2);
  free(GradientBackward1);
  free(GradientBackward2);
}


// Routines for numerical evaluation of the forces from the energy, hessian from the forces etc.
// Note that many problems arise at the cutoff boundary arise when the potentials are unshifted,

const REAL delta=1e-5;

void SaveFrameworkPositionsToReferenceValues(void)
{
  int i,f1;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      Framework[CurrentSystem].Atoms[f1][i].ReferencePosition=
        ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
}

void PlaceFrameworkInBoxFromReferenceValues(void)
{
  int i,f1;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Framework[CurrentSystem].Atoms[f1][i].Position=
          ConvertFromABCtoXYZ(Framework[CurrentSystem].Atoms[f1][i].ReferencePosition);
      Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition = Framework[CurrentSystem].Atoms[f1][i].Position;
    }
}

void SaveAdsorbateAtomPositionsToReferenceValues(void)
{
  int i,m;
  int Type;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfGroups;i++)
      Adsorbates[CurrentSystem][m].Groups[i].CenterOfMassReferencePosition=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][m].Groups[i].CenterOfMassPosition);
    for(i=0;i<Adsorbates[CurrentSystem][m].NumberOfAtoms;i++)
      Adsorbates[CurrentSystem][m].Atoms[i].ReferencePosition=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][m].Atoms[i].Position);
  }
}

void PlaceAdsorbateAtomsInBoxFromReferenceValues(void)
{
  int i,m;
  int Type;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfGroups;i++)
      Adsorbates[CurrentSystem][m].Groups[i].CenterOfMassPosition=ConvertFromABCtoXYZ(Adsorbates[CurrentSystem][m].Groups[i].CenterOfMassReferencePosition);
    for(i=0;i<Components[Type].NumberOfAtoms;i++)
    {
      Adsorbates[CurrentSystem][m].Atoms[i].Position=ConvertFromABCtoXYZ(Adsorbates[CurrentSystem][m].Atoms[i].ReferencePosition);
      Adsorbates[CurrentSystem][m].Atoms[i].AnisotropicPosition=Adsorbates[CurrentSystem][m].Atoms[i].Position;
    }
  }
}

void SaveCationAtomPositionsToReferenceValues(void)
{
  int i,m;
  int Type;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfGroups;i++)
      Cations[CurrentSystem][m].Groups[i].CenterOfMassReferencePosition=ConvertFromXYZtoABC(Cations[CurrentSystem][m].Groups[i].CenterOfMassPosition);
    for(i=0;i<Cations[CurrentSystem][m].NumberOfAtoms;i++)
      Cations[CurrentSystem][m].Atoms[i].ReferencePosition=ConvertFromXYZtoABC(Cations[CurrentSystem][m].Atoms[i].Position);
  }
}

void PlaceCationAtomsInBoxFromReferenceValues(void)
{
  int i,m;
  int Type;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfGroups;i++)
      Cations[CurrentSystem][m].Groups[i].CenterOfMassPosition=ConvertFromABCtoXYZ(Cations[CurrentSystem][m].Groups[i].CenterOfMassReferencePosition);
    for(i=0;i<Components[Type].NumberOfAtoms;i++)
    {
      Cations[CurrentSystem][m].Atoms[i].Position=ConvertFromABCtoXYZ(Cations[CurrentSystem][m].Atoms[i].ReferencePosition);
      Cations[CurrentSystem][m].Atoms[i].AnisotropicPosition=Cations[CurrentSystem][m].Atoms[i].Position;
    }
  }
}

void ComputeGradientsNumerically(REAL *Gradients)
{
  int i,k,m,A,l,Type,index;
  int type;
  VECTOR pos,force,derivative;
  REAL EnergyCentral,EnergyForward1,EnergyForward2,EnergyBackward1,EnergyBackward2;

  derivative.x=derivative.y=derivative.z=0.0;

  index=0;
  for(k=0;k<Framework[CurrentSystem].NumberOfFrameworks;k++)
  {
    if(Framework[CurrentSystem].FrameworkModels[k]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[k];i++)
      {
        pos=Framework[CurrentSystem].Atoms[k][i].Position;
        type=Framework[CurrentSystem].Atoms[k][i].Type;

        // x-direction
        Framework[CurrentSystem].Atoms[k][i].Position=pos;
        CalculateForce();
        force=Framework[CurrentSystem].Atoms[k][i].Force;
        EnergyCentral=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.x=pos.x+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.x=pos.x+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.x=pos.x-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.x=pos.x-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Framework[CurrentSystem].Atoms[k][i].Position=pos;

        derivative.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

        // y-direction
        CalculateForce();
        EnergyCentral=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.y=pos.y+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.y=pos.y+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.y=pos.y-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.y=pos.y-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Framework[CurrentSystem].Atoms[k][i].Position=pos;

        derivative.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

        // z-direction
        CalculateForce();
        EnergyCentral=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.z=pos.z+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.z=pos.z+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.z=pos.z-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.z=pos.z-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Framework[CurrentSystem].Atoms[k][i].Position=pos;

        derivative.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

        Gradients[index++]=derivative.x;
        Gradients[index++]=derivative.y;
        Gradients[index++]=derivative.z;
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;

    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        // store the original position 'pos'
        pos=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition;

        // reference energy
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        CalculateForce();
        force=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassForce;
        EnergyCentral=UTotal[CurrentSystem];

        // x-direction
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        derivative.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

        // y-direction
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        derivative.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

        // z-direction
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        derivative.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

        Gradients[index++]=derivative.x;
        Gradients[index++]=derivative.y;
        Gradients[index++]=derivative.z;
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          // store the original position 'pos'
          pos=Adsorbates[CurrentSystem][m].Atoms[A].Position;

          // reference energy
          Adsorbates[CurrentSystem][m].Atoms[A].Position=pos;
          CalculateForce();
          force=Adsorbates[CurrentSystem][m].Atoms[A].Force;
          EnergyCentral=UTotal[CurrentSystem];

          // x-direction
          Adsorbates[CurrentSystem][m].Atoms[A].Position.x=pos.x+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x=pos.x+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x=pos.x-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x=pos.x-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Adsorbates[CurrentSystem][m].Atoms[A].Position=pos;
          derivative.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

          // y-direction
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y=pos.y+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.y=pos.y+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.y=pos.y-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.y=pos.y-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Adsorbates[CurrentSystem][m].Atoms[A].Position=pos;
          derivative.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

          // z-direction
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z=pos.z+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.z=pos.z+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.z=pos.z-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.z=pos.z-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Adsorbates[CurrentSystem][m].Atoms[A].Position=pos;
          derivative.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

          Gradients[index++]=derivative.x;
          Gradients[index++]=derivative.y;
          Gradients[index++]=derivative.z;
        }
      }

    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;

    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        // store the original position 'pos'
        pos=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition;

        // reference energy
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        CalculateForce();
        force=Cations[CurrentSystem][m].Groups[l].CenterOfMassForce;
        EnergyCentral=UTotal[CurrentSystem];

        // x-direction
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        derivative.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

        // y-direction
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        derivative.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

        // z-direction
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        derivative.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

        Gradients[index++]=derivative.x;
        Gradients[index++]=derivative.y;
        Gradients[index++]=derivative.z;
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          // store the original position 'pos'
          pos=Cations[CurrentSystem][m].Atoms[A].Position;

          // reference energy
          Cations[CurrentSystem][m].Atoms[A].Position=pos;
          CalculateForce();
          force=Cations[CurrentSystem][m].Atoms[A].Force;
          EnergyCentral=UTotal[CurrentSystem];

          // x-direction
          Cations[CurrentSystem][m].Atoms[A].Position.x=pos.x+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.x=pos.x+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.x=pos.x-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.x=pos.x-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Cations[CurrentSystem][m].Atoms[A].Position=pos;
          derivative.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

          // y-direction
          Cations[CurrentSystem][m].Atoms[A].Position.y=pos.y+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.y=pos.y+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.y=pos.y-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.y=pos.y-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Cations[CurrentSystem][m].Atoms[A].Position=pos;
          derivative.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

          // z-direction
          Cations[CurrentSystem][m].Atoms[A].Position.z=pos.z+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.z=pos.z+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.z=pos.z-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.z=pos.z-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Cations[CurrentSystem][m].Atoms[A].Position=pos;
          derivative.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

          Gradients[index++]=derivative.x;
          Gradients[index++]=derivative.y;
          Gradients[index++]=derivative.z;
        }
      }
    }
  }
}


void TestGradientsNumerically(void)
{
  int i;
  REAL *x,Energy;
  REAL *Gradient,*NumericalGradient;
  REAL_MATRIX3x3 StrainFirstDerivative;
  REAL_MATRIX Hessian;
  char buffer[256];
  FILE *FilePtr;
  REAL largest_difference;

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet implemented for this function.");
    exit(0);
  }

  // loop over all the framwork-atoms and fix the order and 'index' into the Hessian
  NumberOfMinimizationVariables=0;
  NumberOfCellMinimizationVariables=0;
  OrderNumberOfMinimiationVariables();

  AllocateMinimizationLocalMemory();

  CurrentSystem=0;

  switch(Ensemble[CurrentSystem])
  {
    case NVE:
    case NVT:
      NumberOfCellMinimizationVariables=0;
      break;
    case NPT:
    case NPH:
      NumberOfCellMinimizationVariables=1;
      NumberOfMinimizationVariables+=1;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ANISOTROPIC:
        case ISOTROPIC:
          NumberOfCellMinimizationVariables=Dimension;
          NumberOfMinimizationVariables+=Dimension;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          NumberOfCellMinimizationVariables=4;
          NumberOfMinimizationVariables+=4;
          BoundaryCondition[CurrentSystem]=TRICLINIC;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          NumberOfCellMinimizationVariables=Dimension*(Dimension+1)/2;
          NumberOfMinimizationVariables+=Dimension*(Dimension+1)/2;
          BoundaryCondition[CurrentSystem]=TRICLINIC;
          break;
      }
  }


  // allocate memory for the generalized coordinates, and derivatives
  x=(REAL*)calloc(NumberOfMinimizationVariables,sizeof(REAL));
  Gradient=(REAL*)calloc(NumberOfMinimizationVariables,sizeof(REAL));
  NumericalGradient=(REAL*)calloc(NumberOfMinimizationVariables,sizeof(REAL));

  CreateGeneralizedCoordinatesFromPositions(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x);
  NumericallyComputeGradient(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x,NumericalGradient);

  CreateGeneralizedCoordinatesFromPositions(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x);
  EvaluateDerivatives(NumberOfCoordinatesMinimizationVariables,&Energy,Gradient,Hessian,&StrainFirstDerivative,TRUE,FALSE);

  mkdir("Numerical",S_IRWXU);

  sprintf(buffer,"Numerical/System_%d",CurrentSystem);
  mkdir(buffer,S_IRWXU);

  sprintf(buffer,"Numerical/System_%d/Gradient.dat",CurrentSystem);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# NumberOfCoordinatesMinimizationVariables: %d\n",NumberOfCoordinatesMinimizationVariables);
  fprintf(FilePtr,"# NumberOfCellMinimizationVariables: %d\n",NumberOfCellMinimizationVariables);
  fprintf(FilePtr,"# Energy: % 18.10f\n",Energy*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"# Strain derivative: % 18.10f % 18.10f % 18.10f\n",StrainFirstDerivative.ax,StrainFirstDerivative.bx,StrainFirstDerivative.cx);
  fprintf(FilePtr,"#                    % 18.10f % 18.10f % 18.10f\n",StrainFirstDerivative.ay,StrainFirstDerivative.by,StrainFirstDerivative.cy);
  fprintf(FilePtr,"#                    % 18.10f % 18.10f % 18.10f\n",StrainFirstDerivative.az,StrainFirstDerivative.bz,StrainFirstDerivative.cz);
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"# Gradients:\n");
  fprintf(FilePtr,"# Column 1: analytical\n");
  fprintf(FilePtr,"# Column 2: numerical\n");
  fprintf(FilePtr,"# Column 3: difference\n");
  fprintf(FilePtr,"# ==========\n");

  largest_difference=0.0;
  for(i=0;i<NumberOfCoordinatesMinimizationVariables;i++)
  {
    fprintf(FilePtr,"%5d %22.6f %22.6f %22.6f\n",i,Gradient[i],NumericalGradient[i],fabs(Gradient[i]-NumericalGradient[i]));
    if(fabs(Gradient[i]-NumericalGradient[i])>largest_difference) largest_difference=fabs(Gradient[i]-NumericalGradient[i]);
  }
  fprintf(FilePtr,"\n");
  if(NumberOfCellMinimizationVariables>0)
  {
    fprintf(FilePtr,"# strain derivatives\n");
    fprintf(FilePtr,"# ==================\n");
    for(;i<NumberOfMinimizationVariables;i++)
    {
      fprintf(FilePtr,"%5d %22.6f %22.6f %22.6f\n",i,Gradient[i],NumericalGradient[i],fabs(Gradient[i]-NumericalGradient[i]));
      if(fabs(Gradient[i]-NumericalGradient[i])>largest_difference) largest_difference=fabs(Gradient[i]-NumericalGradient[i]);
    }
  }


  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"# Largest difference: %22.10f\n\n",largest_difference);
  fflush(FilePtr);

  free(NumericalGradient);
  free(Gradient);
  free(x);
}

void TestHessianNumerically(void)
{
  int i,j;
  REAL *x,Energy;
  REAL *Gradient,*NumericalGradient;
  REAL_MATRIX3x3 StrainFirstDerivative;
  REAL_MATRIX Hessian,NumericalHessian;
  char buffer[256];
  FILE *FilePtr;
  REAL largest_difference,error;

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet implemented for this function.");
    exit(0);
  }

  CurrentSystem=0;
  NumberOfMinimizationVariables=0;
  NumberOfCellMinimizationVariables=0;
  OrderNumberOfMinimiationVariables();

  switch(Ensemble[CurrentSystem])
  {
    case NVE:
    case NVT:
      NumberOfCellMinimizationVariables=0;
      break;
    case NPT:
    case NPH:
      NumberOfCellMinimizationVariables=1;
      NumberOfMinimizationVariables+=1;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ANISOTROPIC:
        case ISOTROPIC:
          NumberOfCellMinimizationVariables=Dimension;
          NumberOfMinimizationVariables+=Dimension;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          NumberOfCellMinimizationVariables=4;
          NumberOfMinimizationVariables+=4;
          BoundaryCondition[CurrentSystem]=TRICLINIC;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          NumberOfCellMinimizationVariables=Dimension*(Dimension+1)/2;
          NumberOfMinimizationVariables+=Dimension*(Dimension+1)/2;
          BoundaryCondition[CurrentSystem]=TRICLINIC;
          break;
      }
  }

  // allocate memory for the generalized coordinates, and derivatives
  x=(REAL*)calloc(NumberOfMinimizationVariables,sizeof(REAL));
  Gradient=(REAL*)calloc(NumberOfMinimizationVariables,sizeof(REAL));
  NumericalGradient=(REAL*)calloc(NumberOfMinimizationVariables,sizeof(REAL));
  Hessian=CreateRealMatrix(NumberOfMinimizationVariables,NumberOfMinimizationVariables);
  NumericalHessian=CreateRealMatrix(NumberOfMinimizationVariables,NumberOfMinimizationVariables);

  CreateGeneralizedCoordinatesFromPositions(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x);
  NumericallyComputeDerivatives(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x,Gradient,NumericalHessian);

  CreateGeneralizedCoordinatesFromPositions(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,x);
  EvaluateDerivatives(NumberOfCoordinatesMinimizationVariables,&Energy,Gradient,Hessian,&StrainFirstDerivative,TRUE,TRUE);

  for(i=0;i<NumberOfMinimizationVariables;i++)
    for(j=0;j<NumberOfMinimizationVariables;j++)
      Hessian.element[j][i]=Hessian.element[i][j];


  mkdir("Numerical",S_IRWXU);

  sprintf(buffer,"Numerical/System_%d",CurrentSystem);
  mkdir(buffer,S_IRWXU);

  sprintf(buffer,"Numerical/System_%d/Hessian.dat",CurrentSystem);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# NumberOfCoordinatesMinimizationVariables: %d\n",NumberOfCoordinatesMinimizationVariables);
  fprintf(FilePtr,"# NumberOfCellMinimizationVariables: %d\n",NumberOfCellMinimizationVariables);
  fprintf(FilePtr,"# first value is the anaytical value, second the numerical value from finite difference\n");

  fprintf(FilePtr,"# Hessian matrix\n");
  fprintf(FilePtr,"# ==============\n");
  largest_difference=0.0;
  for(i=0;i<NumberOfCoordinatesMinimizationVariables;i++)
    for(j=0;j<NumberOfCoordinatesMinimizationVariables;j++)
    {
      fprintf(FilePtr,"%5d %5d %22.6f %22.6f %22.6f\n",i,j,Hessian.element[i][j],NumericalHessian.element[i][j],
              fabs(Hessian.element[i][j]-NumericalHessian.element[i][j]));
      error=fabs(Hessian.element[i][j]-NumericalHessian.element[i][j]);
      if(error>largest_difference) largest_difference=error;
    }
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"# Largest difference Hessian: % 18.10f\n\n",largest_difference);

  if(NumberOfCellMinimizationVariables>0)
  {
    fprintf(FilePtr,"# Cross-term\n");
    fprintf(FilePtr,"# ==========\n");
    largest_difference=0.0;
    for(i=NumberOfCoordinatesMinimizationVariables;i<NumberOfMinimizationVariables;i++)
      for(j=0;j<NumberOfCoordinatesMinimizationVariables;j++)
      {
        fprintf(FilePtr,"%5d %5d %22.6f %22.6f %22.6f\n",i,j,Hessian.element[i][j],NumericalHessian.element[i][j],
                fabs(Hessian.element[i][j]-NumericalHessian.element[i][j]));
        error=fabs(Hessian.element[i][j]-NumericalHessian.element[i][j]);
        if(error>largest_difference) largest_difference=error;
      }
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"# Largest difference cross-term: % 18.10f\n\n",largest_difference);

    fprintf(FilePtr,"# Born-term\n");
    fprintf(FilePtr,"# =========\n");
    largest_difference=0.0;
    for(i=NumberOfCoordinatesMinimizationVariables;i<NumberOfMinimizationVariables;i++)
      for(j=NumberOfCoordinatesMinimizationVariables;j<NumberOfMinimizationVariables;j++)
      {
        fprintf(FilePtr,"%5d %5d %22.6f %22.6f %22.6f\n",i,j,Hessian.element[i][j],NumericalHessian.element[i][j],
                fabs(Hessian.element[i][j]-NumericalHessian.element[i][j]));
        error=fabs(Hessian.element[i][j]-NumericalHessian.element[i][j]);
        if(error>largest_difference) largest_difference=error;
      }
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"# Largest difference Born-term: % 18.10f\n\n",largest_difference);
  }

  fclose(FilePtr);

  free(NumericalGradient);
  free(Gradient);
  free(x);
  DeleteRealMatrix(Hessian);
  DeleteRealMatrix(NumericalHessian);
}


void TestForcesNumerically(void)
{
  int i,k,m,A,l,Type;
  int type;
  VECTOR pos,force,derivative,second_derivative,largest_difference;
  REAL EnergyCentral,EnergyForward1,EnergyForward2,EnergyBackward1,EnergyBackward2;
  char buffer[256];
  FILE *FilePtr;

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet implemented for this function.");
    exit(0);
  }

  derivative.x=derivative.y=derivative.z=0.0;
  second_derivative.x=second_derivative.y=second_derivative.z=0.0;
  largest_difference.x=largest_difference.y=largest_difference.z=0.0;

  mkdir("Numerical",S_IRWXU);

  sprintf(buffer,"Numerical/System_%d",CurrentSystem);
  mkdir(buffer,S_IRWXU);

  sprintf(buffer,"Numerical/System_%d/Forces.dat",CurrentSystem);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# Force x,y,z, first line analytical, second line minus numerical value\n");
  fprintf(FilePtr,"# The force should be the exact opposite to the numerical derivative\n");
  fprintf(OutputFilePtr[CurrentSystem],"# Force x,y,z, first line analytical, second line numerical value\n");
  fprintf(OutputFilePtr[CurrentSystem],"# The force should be the exact opposite to the numerical derivative\n");

  fprintf(FilePtr,"# Framework atoms:\n");
  fprintf(FilePtr,"# ================\n");
  fprintf(OutputFilePtr[CurrentSystem],"# Framework atoms:\n");
  fprintf(OutputFilePtr[CurrentSystem],"# ================\n");

  ConstructBondDipolesFromBondsFramework();

  for(k=0;k<Framework[CurrentSystem].NumberOfFrameworks;k++)
  {
    if(Framework[CurrentSystem].FrameworkModels[k]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[k];i++)
      {
        pos=Framework[CurrentSystem].Atoms[k][i].Position;
        type=Framework[CurrentSystem].Atoms[k][i].Type;

        // x-direction
        Framework[CurrentSystem].Atoms[k][i].Position=pos;
        CalculateForce();
        force=Framework[CurrentSystem].Atoms[k][i].Force;
        EnergyCentral=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.x=pos.x+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.x=pos.x+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.x=pos.x-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.x=pos.x-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Framework[CurrentSystem].Atoms[k][i].Position=pos;

        derivative.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
        second_derivative.x=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

        // y-direction
        CalculateForce();
        EnergyCentral=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.y=pos.y+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.y=pos.y+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.y=pos.y-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.y=pos.y-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Framework[CurrentSystem].Atoms[k][i].Position=pos;

        derivative.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
        second_derivative.y=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

        // z-direction
        CalculateForce();
        EnergyCentral=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.z=pos.z+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.z=pos.z+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.z=pos.z-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Framework[CurrentSystem].Atoms[k][i].Position.z=pos.z-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Framework[CurrentSystem].Atoms[k][i].Position=pos;
        CalculateForce();

        derivative.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
        second_derivative.z=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

        if(fabs(force.x+derivative.x)>largest_difference.x) largest_difference.x=fabs(force.x+derivative.x);
        if(fabs(force.y+derivative.y)>largest_difference.y) largest_difference.y=fabs(force.y+derivative.y);
        if(fabs(force.z+derivative.z)>largest_difference.z) largest_difference.z=fabs(force.z+derivative.z);


         fprintf(FilePtr,"%4d %4d [%8s] % 18.10f % 18.10f % 18.10f\n",
           k,i,PseudoAtoms[type].Name,
           (double)force.x,(double)force.y,(double)force.z);
        fprintf(FilePtr,"                     % 18.10f % 18.10f % 18.10f\n",
           (double)(-derivative.x),(double)(-derivative.y),(double)(-derivative.z));
        fprintf(OutputFilePtr[CurrentSystem],"%4d %4d [%8s] % 18.10f % 18.10f % 18.10f\n",
           k,i,PseudoAtoms[type].Name,
            (double)force.x,(double)force.y,(double)force.z);
        fprintf(OutputFilePtr[CurrentSystem],"                     % 18.10f % 18.10f % 18.10f\n",
           (double)(-derivative.x),(double)(-derivative.y),(double)(-derivative.z));

        fflush(OutputFilePtr[CurrentSystem]);
        fflush(FilePtr);
      }
    }
  }
  fprintf(FilePtr,"Largest difference: % 18.10f % 18.10f % 18.10f\n\n",
       largest_difference.x,largest_difference.y,largest_difference.z);
  fprintf(OutputFilePtr[CurrentSystem],"Largest difference: % 18.10f % 18.10f % 18.10f\n\n",
       largest_difference.x,largest_difference.y,largest_difference.z);
  fflush(FilePtr);
  fflush(OutputFilePtr[CurrentSystem]);


  largest_difference.x=largest_difference.y=largest_difference.z=0.0;
  if(NumberOfAdsorbateMolecules[CurrentSystem]>0)
  {
    fprintf(FilePtr,"# Adsorbate atoms:\n");
    fprintf(FilePtr,"# ================\n");
    fprintf(OutputFilePtr[CurrentSystem],"# Adsorbate atoms:\n");
    fprintf(OutputFilePtr[CurrentSystem],"# ================\n");
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;

    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        // store the original position 'pos'
        pos=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition;

        // reference energy
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        CalculateForce();
        force=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassForce;
        EnergyCentral=UTotal[CurrentSystem];

        // x-direction
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        derivative.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
        second_derivative.x=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

        // y-direction
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        derivative.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
        second_derivative.y=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

        // z-direction
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        CalculateForce();

        derivative.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
        second_derivative.z=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

        if(fabs(force.x+derivative.x)>largest_difference.x) largest_difference.x=fabs(force.x+derivative.x);
        if(fabs(force.y+derivative.y)>largest_difference.y) largest_difference.y=fabs(force.y+derivative.y);
        if(fabs(force.z+derivative.z)>largest_difference.z) largest_difference.z=fabs(force.z+derivative.z);

          fprintf(FilePtr,"Rigid    mol=%4d group=%4d            % 18.10f % 18.10f % 18.10f\n",
             m,l,
             (double)force.x,(double)force.y,(double)force.z);
          fprintf(FilePtr,"                                        % 18.10f % 18.10f % 18.10f\n",
             (double)(-derivative.x),(double)(-derivative.y),(double)(-derivative.z));
          fprintf(OutputFilePtr[CurrentSystem],"Rigid    mol=%4d group=%4d            % 18.10f % 18.10f % 18.10f\n",
             m,l,
             (double)force.x,(double)force.y,(double)force.z);
          fprintf(OutputFilePtr[CurrentSystem],"                                        % 18.10f % 18.10f % 18.10f\n",
             (double)(-derivative.x),(double)(-derivative.y),(double)(-derivative.z));
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          // store the original position 'pos'
          pos=Adsorbates[CurrentSystem][m].Atoms[A].Position;

          // reference energy
          Adsorbates[CurrentSystem][m].Atoms[A].Position=pos;
          CalculateForce();
          force=Adsorbates[CurrentSystem][m].Atoms[A].Force;
          EnergyCentral=UTotal[CurrentSystem];

          // x-direction
          Adsorbates[CurrentSystem][m].Atoms[A].Position.x=pos.x+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x=pos.x+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x=pos.x-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.x=pos.x-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Adsorbates[CurrentSystem][m].Atoms[A].Position=pos;
          derivative.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
          second_derivative.x=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

          // y-direction
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y=pos.y+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.y=pos.y+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.y=pos.y-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.y=pos.y-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Adsorbates[CurrentSystem][m].Atoms[A].Position=pos;
          derivative.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
          second_derivative.y=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

          // z-direction
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z=pos.z+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.z=pos.z+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.z=pos.z-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Adsorbates[CurrentSystem][m].Atoms[A].Position.z=pos.z-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Adsorbates[CurrentSystem][m].Atoms[A].Position=pos;
          CalculateForce();

          derivative.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
          second_derivative.z=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);


          if(fabs(force.x+derivative.x)>largest_difference.x) largest_difference.x=fabs(force.x+derivative.x);
          if(fabs(force.y+derivative.y)>largest_difference.y) largest_difference.y=fabs(force.y+derivative.y);
          if(fabs(force.z+derivative.z)>largest_difference.z) largest_difference.z=fabs(force.z+derivative.z);

          type=Adsorbates[CurrentSystem][m].Atoms[A].Type;
          fprintf(FilePtr,"Flexible mol=%4d  atom=%4d [%8s] % 18.10f % 18.10f % 18.10f\n",
             m,A,PseudoAtoms[type].Name,
             (double)force.x,(double)force.y,(double)force.z);
          fprintf(FilePtr,"                                        % 18.10f % 18.10f % 18.10f\n",
             (double)(-derivative.x),(double)(-derivative.y),(double)(-derivative.z));
          fprintf(OutputFilePtr[CurrentSystem],"Flexible mol=%4d  atom=%4d [%8s] % 18.10f % 18.10f % 18.10f\n",
             m,A,PseudoAtoms[type].Name,
             (double)force.x,(double)force.y,(double)force.z);
          fprintf(OutputFilePtr[CurrentSystem],"                                        % 18.10f % 18.10f % 18.10f\n",
             (double)(-derivative.x),(double)(-derivative.y),(double)(-derivative.z));
        }
      }

    }
  }

  if(NumberOfAdsorbateMolecules[CurrentSystem]>0)
  {
    fprintf(FilePtr,"Largest difference: % 18.10f % 18.10f % 18.10f\n\n",
         largest_difference.x,largest_difference.y,largest_difference.z);
    fprintf(OutputFilePtr[0],"Largest difference: % 18.10f % 18.10f % 18.10f\n\n",
         largest_difference.x,largest_difference.y,largest_difference.z);
  }

  largest_difference.x=largest_difference.y=largest_difference.z=0.0;
  if(NumberOfCationMolecules[CurrentSystem]>0)
  {
    fprintf(FilePtr,"# Cation atoms:\n");
    fprintf(FilePtr,"# =============\n");
    fprintf(OutputFilePtr[CurrentSystem],"# Cation atoms:\n");
    fprintf(OutputFilePtr[CurrentSystem],"# =============\n");
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;

    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        // store the original position 'pos'
        pos=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition;

        // reference energy
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        CalculateForce();
        force=Cations[CurrentSystem][m].Groups[l].CenterOfMassForce;
        EnergyCentral=UTotal[CurrentSystem];

        // x-direction
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=pos.x-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        derivative.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
        second_derivative.x=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

        // y-direction
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=pos.y-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        derivative.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
        second_derivative.y=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

        // z-direction
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z+delta;
        CalculateForce();
        EnergyForward2=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z+0.5*delta;
        CalculateForce();
        EnergyForward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z-0.5*delta;
        CalculateForce();
        EnergyBackward1=UTotal[CurrentSystem];

        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=pos.z-delta;
        CalculateForce();
        EnergyBackward2=UTotal[CurrentSystem];

        // restore original position
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition=pos;
        CalculateForce();

        derivative.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
        second_derivative.z=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

        if(fabs(force.x+derivative.x)>largest_difference.x) largest_difference.x=fabs(force.x+derivative.x);
        if(fabs(force.y+derivative.y)>largest_difference.y) largest_difference.y=fabs(force.y+derivative.y);
        if(fabs(force.z+derivative.z)>largest_difference.z) largest_difference.z=fabs(force.z+derivative.z);

          fprintf(FilePtr,"Rigid    mol=%4d group=%4d            % 18.10f % 18.10f % 18.10f\n",
             m,l,
             (double)force.x,(double)force.y,(double)force.z);
          fprintf(FilePtr,"                                        % 18.10f % 18.10f % 18.10f\n",
             (double)(-derivative.x),(double)(-derivative.y),(double)(-derivative.z));
          fprintf(OutputFilePtr[CurrentSystem],"Rigid    mol=%4d group=%4d            % 18.10f % 18.10f % 18.10f\n",
             m,l,
             (double)force.x,(double)force.y,(double)force.z);
          fprintf(OutputFilePtr[CurrentSystem],"                                        % 18.10f % 18.10f % 18.10f\n",
             (double)(-derivative.x),(double)(-derivative.y),(double)(-derivative.z));
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          // store the original position 'pos'
          pos=Cations[CurrentSystem][m].Atoms[A].Position;

          // reference energy
          Cations[CurrentSystem][m].Atoms[A].Position=pos;
          CalculateForce();
          force=Cations[CurrentSystem][m].Atoms[A].Force;
          EnergyCentral=UTotal[CurrentSystem];

          // x-direction
          Cations[CurrentSystem][m].Atoms[A].Position.x=pos.x+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.x=pos.x+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.x=pos.x-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.x=pos.x-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Cations[CurrentSystem][m].Atoms[A].Position=pos;
          derivative.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
          second_derivative.x=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

          // y-direction
          Cations[CurrentSystem][m].Atoms[A].Position.y=pos.y+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.y=pos.y+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.y=pos.y-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.y=pos.y-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Cations[CurrentSystem][m].Atoms[A].Position=pos;
          derivative.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
          second_derivative.y=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

          // z-direction
          Cations[CurrentSystem][m].Atoms[A].Position.z=pos.z+delta;
          CalculateForce();
          EnergyForward2=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.z=pos.z+0.5*delta;
          CalculateForce();
          EnergyForward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.z=pos.z-0.5*delta;
          CalculateForce();
          EnergyBackward1=UTotal[CurrentSystem];

          Cations[CurrentSystem][m].Atoms[A].Position.z=pos.z-delta;
          CalculateForce();
          EnergyBackward2=UTotal[CurrentSystem];

          // restore original position
          Cations[CurrentSystem][m].Atoms[A].Position=pos;
          CalculateForce();

          derivative.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
          second_derivative.z=(EnergyForward2-2.0*EnergyCentral+EnergyBackward2)/SQR(delta);

          if(fabs(force.x+derivative.x)>largest_difference.x) largest_difference.x=fabs(force.x+derivative.x);
          if(fabs(force.y+derivative.y)>largest_difference.y) largest_difference.y=fabs(force.y+derivative.y);
          if(fabs(force.z+derivative.z)>largest_difference.z) largest_difference.z=fabs(force.z+derivative.z);

          type=Cations[CurrentSystem][m].Atoms[A].Type;
          fprintf(FilePtr,"Flexible mol=%4d  atom=%4d [%8s] % 18.10f % 18.10f % 18.10f\n",
             m,A,PseudoAtoms[type].Name,
             (double)force.x,(double)force.y,(double)force.z);
          fprintf(FilePtr,"                                        % 18.10f % 18.10f % 18.10f\n",
             (double)(-derivative.x),(double)(-derivative.y),(double)(-derivative.z));
          fprintf(OutputFilePtr[CurrentSystem],"Flexible mol=%4d  atom=%4d [%8s] % 18.10f % 18.10f % 18.10f\n",
             m,A,PseudoAtoms[type].Name,
             (double)force.x,(double)force.y,(double)force.z);
          fprintf(OutputFilePtr[CurrentSystem],"                                        % 18.10f % 18.10f % 18.10f\n",
             (double)(-derivative.x),(double)(-derivative.y),(double)(-derivative.z));
        }
      }
    }
  }

  if(NumberOfCationMolecules[CurrentSystem]>0)
  {
    fprintf(FilePtr,"Largest difference: % 18.10f % 18.10f % 18.10f\n",
         largest_difference.x,largest_difference.y,largest_difference.z);
    fprintf(OutputFilePtr[CurrentSystem],"Largest difference: % 18.10f % 18.10f % 18.10f\n",
         largest_difference.x,largest_difference.y,largest_difference.z);
  }

  fclose(FilePtr);
}

void TestElectricField(void)
{
  int i,j,l,f1;
  FILE *FilePtr;
  REAL Charge;
  char buffer[2024];

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet implemented for this function.");
    exit(0);
  }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Framework[CurrentSystem].Atoms[f1][i].Force.x=0.0;
      Framework[CurrentSystem].Atoms[f1][i].Force.y=0.0;
      Framework[CurrentSystem].Atoms[f1][i].Force.z=0.0;
      Framework[CurrentSystem].Atoms[f1][i].ElectricField.x=0.0;
      Framework[CurrentSystem].Atoms[f1][i].ElectricField.y=0.0;
      Framework[CurrentSystem].Atoms[f1][i].ElectricField.z=0.0;
      Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField.x=0.0;
      Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField.y=0.0;
      Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField.z=0.0;
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Force.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].Force.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].Force.z=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField.z=0.0;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Force.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].Force.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].Force.z=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField.z=0.0;
    }
  }

  // contribution of the intermolecular interactions (between adsorbates and/or cations)
  CalculateTotalInterChargeChargeCoulombForce();
  CalculateFrameworkAdsorbateChargeChargeForce();
  CalculateFrameworkCationChargeChargeForce();

  // the contribution of the charges present in the system (long-range interaction)
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    EwaldFourierForce();

/*
  CalculateFrameworkAdsorbateChargeInducedDipoleForce();
  CalculateFrameworkAdsorbateInducedDipoleInducedDipoleForce();

  CalculateInterChargeInducedDipoleForce();
  CalculateInterInducedDipoleInducedDipoleForce();

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    ComputeInducedDipolesForcesEwald();
*/

  CalculateElectricField();

  mkdir("Numerical",S_IRWXU);

  sprintf(buffer,"Numerical/System_%d",CurrentSystem);
  mkdir(buffer,S_IRWXU);

  sprintf(buffer,"Numerical/System_%d/ElectricField.dat",CurrentSystem);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# Framework:\n");
  fprintf(FilePtr,"#=================================================================================\n");
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      fprintf(FilePtr,"%d %18.10f %18.10f %18.10f <-> %18.10f %18.10f %18.10f\n",i,
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.x,
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.y,
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.z,
        Framework[CurrentSystem].Atoms[f1][i].Force.x/Charge,
        Framework[CurrentSystem].Atoms[f1][i].Force.y/Charge,
        Framework[CurrentSystem].Atoms[f1][i].Force.z/Charge);
    }
  }
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"# Adsorbates:\n");
  fprintf(FilePtr,"#=================================================================================\n");
  for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
    for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
    {
      Charge=Adsorbates[CurrentSystem][j].Atoms[l].Charge;
      fprintf(FilePtr,"%d %18.10f %18.10f %18.10f <-> %18.10f %18.10f %18.10f\n",i,
        Adsorbates[CurrentSystem][j].Atoms[l].ReferenceElectricField.x+Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.x,
        Adsorbates[CurrentSystem][j].Atoms[l].ReferenceElectricField.y+Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.y,
        Adsorbates[CurrentSystem][j].Atoms[l].ReferenceElectricField.z+Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.z,
        Adsorbates[CurrentSystem][j].Atoms[l].Force.x/Charge,
        Adsorbates[CurrentSystem][j].Atoms[l].Force.y/Charge,
        Adsorbates[CurrentSystem][j].Atoms[l].Force.z/Charge);
    }
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"# Cations:\n");
  fprintf(FilePtr,"#=================================================================================\n");
  for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
    for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
    {
      Charge=Cations[CurrentSystem][j].Atoms[l].Charge;
      fprintf(FilePtr,"%d %18.10f %18.10f %18.10f <-> %18.10f %18.10f %18.10f\n",i,
        Cations[CurrentSystem][j].Atoms[l].ElectricField.x,
        Cations[CurrentSystem][j].Atoms[l].ElectricField.y,
        Cations[CurrentSystem][j].Atoms[l].ElectricField.z,
        Cations[CurrentSystem][j].Atoms[l].Force.x/Charge,
        Cations[CurrentSystem][j].Atoms[l].Force.y/Charge,
        Cations[CurrentSystem][j].Atoms[l].Force.z/Charge);
    }
  fprintf(FilePtr,"\n");

  fclose(FilePtr);
}

REAL_MATRIX3x3 ComputeStrainDerivativeNumerically(void)
{
  int i,j;
  REAL pressure,temp;
  REAL EnergyCentral,EnergyForward1,EnergyForward2,EnergyBackward1,EnergyBackward2,det;
  REAL_MATRIX3x3 StoredBox,StoredReplicaBox,strain;
  REAL_MATRIX3x3 StrainDerivative;
  int StoredBoundaryCondition;
  const REAL delta=1e-7;
  int ncell,k1,k2,k3;

  // store the box, boundary-condition
  StoredBox=Box[CurrentSystem];
  StoredReplicaBox=ReplicaBox[CurrentSystem];
  StoredBoundaryCondition=BoundaryCondition[CurrentSystem];

  ComputeBornTerm=FALSE;

  // switch to boundary conditions
  BoundaryCondition[CurrentSystem]=TRICLINIC;

  SetupKVectors();
  CalculateForce();
  EnergyCentral=UTotal[CurrentSystem];

  SaveFrameworkPositionsToReferenceValues();
  SaveAdsorbateAtomPositionsToReferenceValues();
  SaveCationAtomPositionsToReferenceValues();

  // ax-element
  strain.ax=1.0+delta; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;       strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward2=UTotal[CurrentSystem];

  strain.ax=1.0+0.5*delta; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;           strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;           strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward1=UTotal[CurrentSystem];

  strain.ax=1.0-0.5*delta; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;           strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;           strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward1=UTotal[CurrentSystem];

  strain.ax=1.0-delta; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;       strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward2=UTotal[CurrentSystem];
  StrainDerivative.ax=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

  // bx-element
  strain.ax=1.0;       strain.bx=delta;     strain.cx=0.0;
  strain.ay=0.0;       strain.by=1.0;       strain.cy=0.0;
  strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward2=UTotal[CurrentSystem];

  strain.ax=1.0;        strain.bx=0.5*delta;  strain.cx=0.0;
  strain.ay=0.0;        strain.by=1.0;        strain.cy=0.0;
  strain.az=0.0;        strain.bz=0.0;        strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward1=UTotal[CurrentSystem];

  strain.ax=1.0;         strain.bx=-0.5*delta;  strain.cx=0.0;
  strain.ay=0.0;         strain.by=1.0;         strain.cy=0.0;
  strain.az=0.0;         strain.bz=0.0;         strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward1=UTotal[CurrentSystem];

  strain.ax=1.0;        strain.bx=-delta;     strain.cx=0.0;
  strain.ay=0.0;        strain.by=1.0;        strain.cy=0.0;
  strain.az=0.0;        strain.bz=0.0;        strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward2=UTotal[CurrentSystem];
  StrainDerivative.bx=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

  // ay-element
  strain.ax=1.0;       strain.bx=0.0;       strain.cx=0.0;
  strain.ay=delta;     strain.by=1.0;       strain.cy=0.0;
  strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward2=UTotal[CurrentSystem];

  strain.ax=1.0;        strain.bx=0.0;        strain.cx=0.0;
  strain.ay=0.5*delta;  strain.by=1.0;        strain.cy=0.0;
  strain.az=0.0;        strain.bz=0.0;        strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward1=UTotal[CurrentSystem];

  strain.ax=1.0;         strain.bx=0.0;         strain.cx=0.0;
  strain.ay=-0.5*delta;  strain.by=1.0;         strain.cy=0.0;
  strain.az=0.0;         strain.bz=0.0;         strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward1=UTotal[CurrentSystem];

  strain.ax=1.0;        strain.bx=0.0;        strain.cx=0.0;
  strain.ay=-delta;     strain.by=1.0;        strain.cy=0.0;
  strain.az=0.0;        strain.bz=0.0;        strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward2=UTotal[CurrentSystem];
  StrainDerivative.ay=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

  // cx-element
  strain.ax=1.0;       strain.bx=0.0; strain.cx=delta;
  strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;       strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward2=UTotal[CurrentSystem];

  strain.ax=1.0;        strain.bx=0.0; strain.cx=0.5*delta;
  strain.ay=0.0;        strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;        strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward1=UTotal[CurrentSystem];

  strain.ax=1.0;         strain.bx=0.0; strain.cx=-0.5*delta;
  strain.ay=0.0;         strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;         strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward1=UTotal[CurrentSystem];

  strain.ax=1.0;        strain.bx=0.0; strain.cx=-delta;
  strain.ay=0.0;        strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;        strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward2=UTotal[CurrentSystem];
  StrainDerivative.cx=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

  // az-element
  strain.ax=1.0;       strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
  strain.az=delta;     strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseBox[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward2=UTotal[CurrentSystem];

  strain.ax=1.0;        strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;        strain.by=1.0; strain.cy=0.0;
  strain.az=0.5*delta;  strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward1=UTotal[CurrentSystem];

  strain.ax=1.0;         strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;         strain.by=1.0; strain.cy=0.0;
  strain.az=-0.5*delta;  strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward1=UTotal[CurrentSystem];

  strain.ax=1.0;        strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;        strain.by=1.0; strain.cy=0.0;
  strain.az=-delta;     strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward2=UTotal[CurrentSystem];
  StrainDerivative.az=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

  // by-element
  strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0+delta; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0;       strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward2=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;           strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0+0.5*delta; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0;           strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward1=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;           strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0-0.5*delta; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0;           strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward1=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0-delta; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0;       strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward2=UTotal[CurrentSystem];
  StrainDerivative.by=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

  // cy-element
  strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;       strain.cy=delta;
  strain.az=0.0; strain.bz=0.0;       strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward2=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;        strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;        strain.cy=0.5*delta;
  strain.az=0.0; strain.bz=0.0;        strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward1=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;         strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;         strain.cy=-0.5*delta;
  strain.az=0.0; strain.bz=0.0;         strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward1=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;        strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;        strain.cy=-delta;
  strain.az=0.0; strain.bz=0.0;        strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward2=UTotal[CurrentSystem];
  StrainDerivative.cy=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);


  // bz-element
  strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;       strain.cy=0.0;
  strain.az=0.0; strain.bz=delta;     strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward2=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;        strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;        strain.cy=0.0;
  strain.az=0.0; strain.bz=0.5*delta;  strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward1=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;         strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;         strain.cy=0.0;
  strain.az=0.0; strain.bz=-0.5*delta;  strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward1=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;        strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;        strain.cy=0.0;
  strain.az=0.0; strain.bz=-delta;     strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward2=UTotal[CurrentSystem];
  StrainDerivative.bz=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

  // cz-element
  strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0; strain.cz=1.0+delta;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward2=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0; strain.cz=1.0+0.5*delta;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyForward1=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0; strain.cz=1.0-0.5*delta;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward1=UTotal[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0; strain.cz=1.0-delta;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  EnergyBackward2=UTotal[CurrentSystem];
  StrainDerivative.cz=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);

  // restore initial box, positions, and boundary condition
  Box[CurrentSystem]=StoredBox;
  ReplicaBox[CurrentSystem]=StoredReplicaBox;
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  BoundaryCondition[CurrentSystem]=StoredBoundaryCondition;
  SetupKVectors();


  // correction for tail-correction
  pressure=0.0;
  for(i=0;i<NumberOfPseudoAtoms;i++)
    for(j=0;j<NumberOfPseudoAtoms;j++)
    {
      if(TailCorrection[i][j])
        pressure+=(2.0/3.0)*M_PI*CUBE(CutOffVDW)*(REAL)NumberOfPseudoAtomsType[CurrentSystem][i]*(REAL)NumberOfPseudoAtomsType[CurrentSystem][j]*
                  PotentialValue(i,j,CutOffVDWSquared,1.0);
    }
  StrainDerivative.ax-=pressure/Volume[CurrentSystem];
  StrainDerivative.by-=pressure/Volume[CurrentSystem];
  StrainDerivative.cz-=pressure/Volume[CurrentSystem];


  // symmetrize strain-derivative (asymmetry originates from torque induced by rigid units)
  temp=0.5*(StrainDerivative.ay+StrainDerivative.bx);
  StrainDerivative.ay=StrainDerivative.bx=temp;
  temp=0.5*(StrainDerivative.az+StrainDerivative.cx);
  StrainDerivative.az=StrainDerivative.cx=temp;
  temp=0.5*(StrainDerivative.bz+StrainDerivative.cy);
  StrainDerivative.bz=StrainDerivative.cy=temp;

  return StrainDerivative;
}


REAL_MATRIX9x9 ComputeStrainSecondDerivativeNumerically(void)
{
  REAL det;
  REAL_MATRIX3x3 StrainDerivativeCentral,StrainDerivativeForward1,StrainDerivativeForward2;
  REAL_MATRIX3x3 StrainDerivativeBackward1,StrainDerivativeBackward2;
  REAL_MATRIX3x3 StoredBox,StoredReplicaBox,strain;
  REAL_MATRIX9x9 StrainSecondDerivative;
  int StoredBoundaryCondition;
  int ncell,k1,k2,k3;
  const REAL delta=1e-8;

  StoredBox=Box[CurrentSystem];
  StoredReplicaBox=ReplicaBox[CurrentSystem];
  StoredBoundaryCondition=BoundaryCondition[CurrentSystem];
  BoundaryCondition[CurrentSystem]=TRICLINIC;

  ComputeBornTerm=TRUE;

  InitializeMatrix9x9(&StrainSecondDerivative);

  SetupKVectors();
  CalculateForce();
  StrainDerivativeCentral=StrainDerivativeTensor[CurrentSystem];

  SaveFrameworkPositionsToReferenceValues();
  SaveAdsorbateAtomPositionsToReferenceValues();
  SaveCationAtomPositionsToReferenceValues();

  // ax-element
  strain.ax=1.0+delta; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;       strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0+0.5*delta; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;           strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;           strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0-0.5*delta; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;           strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;           strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0-delta; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
  strain.az=0.0;       strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

  StrainSecondDerivative.xxxx=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
  StrainSecondDerivative.xxyy=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
  StrainSecondDerivative.xxzz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
  StrainSecondDerivative.xxyz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
  StrainSecondDerivative.xxzx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
  StrainSecondDerivative.xxxy=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);

  // bx/ay-element
  strain.ax=1.0;       strain.bx=0.5*delta; strain.cx=0.0;
  strain.ay=0.5*delta; strain.by=1.0;       strain.cy=0.0;
  strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0;        strain.bx=0.25*delta; strain.cx=0.0;
  strain.ay=0.25*delta; strain.by=1.0;        strain.cy=0.0;
  strain.az=0.0;        strain.bz=0.0;        strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0;         strain.bx=-0.25*delta; strain.cx=0.0;
  strain.ay=-0.25*delta; strain.by=1.0;         strain.cy=0.0;
  strain.az=0.0;         strain.bz=0.0;         strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0;        strain.bx=-0.5*delta; strain.cx=0.0;
  strain.ay=-0.5*delta; strain.by=1.0;        strain.cy=0.0;
  strain.az=0.0;        strain.bz=0.0;        strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

  StrainSecondDerivative.xyxy=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
  StrainSecondDerivative.yyxy=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
  StrainSecondDerivative.yzxy=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
  StrainSecondDerivative.zxxy=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
  StrainSecondDerivative.zzxy=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

  // cx/az-element
  strain.ax=1.0;       strain.bx=0.0; strain.cx=0.5*delta;
  strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
  strain.az=0.5*delta; strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0;        strain.bx=0.0; strain.cx=0.25*delta;
  strain.ay=0.0;        strain.by=1.0; strain.cy=0.0;
  strain.az=0.25*delta; strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0;         strain.bx=0.0; strain.cx=-0.25*delta;
  strain.ay=0.0;         strain.by=1.0; strain.cy=0.0;
  strain.az=-0.25*delta; strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0;        strain.bx=0.0; strain.cx=-0.5*delta;
  strain.ay=0.0;        strain.by=1.0; strain.cy=0.0;
  strain.az=-0.5*delta; strain.bz=0.0; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

  StrainSecondDerivative.zxzx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
  StrainSecondDerivative.zzzx=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
  StrainSecondDerivative.yzzx=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);

  // by-element
  strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0+delta; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0;       strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;           strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0+0.5*delta; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0;           strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;           strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0-0.5*delta; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0;           strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0-delta; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0;       strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

  StrainSecondDerivative.yyyy=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
  StrainSecondDerivative.yyzz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
  StrainSecondDerivative.yyyz=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
  StrainSecondDerivative.yyzx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);

  // cy/bz-element
  strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;       strain.cy=0.5*delta;
  strain.az=0.0; strain.bz=0.5*delta; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;        strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;        strain.cy=0.25*delta;
  strain.az=0.0; strain.bz=0.25*delta; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;         strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;         strain.cy=-0.25*delta;
  strain.az=0.0; strain.bz=-0.25*delta; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0;        strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0;        strain.cy=-0.5*delta;
  strain.az=0.0; strain.bz=-0.5*delta; strain.cz=1.0;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

  StrainSecondDerivative.yzyz=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
  StrainSecondDerivative.zzyz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

  // cz-element
  strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0; strain.cz=1.0+delta;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0; strain.cz=1.0+0.5*delta;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0; strain.cz=1.0-0.5*delta;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

  strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
  strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
  strain.az=0.0; strain.bz=0.0; strain.cz=1.0-delta;
  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  ReplicaBox[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredReplicaBox);
  InverseBoxMatrix(&ReplicaBox[CurrentSystem],&InverseReplicaBox[CurrentSystem]);
  ncell=0;
  for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
      {
        ReplicaShift[ncell].x=Box[CurrentSystem].ax*k1+Box[CurrentSystem].bx*k2+Box[CurrentSystem].cx*k3;
        ReplicaShift[ncell].y=Box[CurrentSystem].ay*k1+Box[CurrentSystem].by*k2+Box[CurrentSystem].cy*k3;
        ReplicaShift[ncell].z=Box[CurrentSystem].az*k1+Box[CurrentSystem].bz*k2+Box[CurrentSystem].cz*k3;
        ncell++;
      }
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();
  SetupKVectors();
  CalculateForce();
  StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

  StrainSecondDerivative.zzzz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

  // use symmetry Cijkl = Cjikl = Cijlk = Cjilk = Cklij
  StrainSecondDerivative.yxxx=StrainSecondDerivative.xxxy;
  StrainSecondDerivative.xxyx=StrainSecondDerivative.xxxy;
  StrainSecondDerivative.yxyx=StrainSecondDerivative.xyxy;
  StrainSecondDerivative.zxyx=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.yxzx=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.zyxx=StrainSecondDerivative.xxyz;
  StrainSecondDerivative.xyyx=StrainSecondDerivative.xyxy;
  StrainSecondDerivative.yyyx=StrainSecondDerivative.yyxy;
  StrainSecondDerivative.zyyx=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.zyzx=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.xzxx=StrainSecondDerivative.xxzx;
  StrainSecondDerivative.xzyx=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.yzyx=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.zzyx=StrainSecondDerivative.zzxy;
  StrainSecondDerivative.xzzx=StrainSecondDerivative.zxzx;
  StrainSecondDerivative.yxxy=StrainSecondDerivative.xyxy;
  StrainSecondDerivative.yxyy=StrainSecondDerivative.yyxy;
  StrainSecondDerivative.xxzy=StrainSecondDerivative.xxyz;
  StrainSecondDerivative.yxzy=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.zxzy=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.zyxy=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.zyyy=StrainSecondDerivative.yyyz;
  StrainSecondDerivative.xyzy=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.yyzy=StrainSecondDerivative.yyyz;
  StrainSecondDerivative.zyzy=StrainSecondDerivative.yzyz;
  StrainSecondDerivative.xzxy=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.xzyy=StrainSecondDerivative.yyzx;
  StrainSecondDerivative.xzzy=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.yzzy=StrainSecondDerivative.yzyz;
  StrainSecondDerivative.zzzy=StrainSecondDerivative.zzyz;
  StrainSecondDerivative.xxxz=StrainSecondDerivative.xxzx;
  StrainSecondDerivative.yxxz=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.zxxz=StrainSecondDerivative.zxzx;
  StrainSecondDerivative.yxyz=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.yxzz=StrainSecondDerivative.zzxy;
  StrainSecondDerivative.xyxz=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.yyxz=StrainSecondDerivative.yyzx;
  StrainSecondDerivative.zyxz=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.zyyz=StrainSecondDerivative.yzyz;
  StrainSecondDerivative.zyzz=StrainSecondDerivative.zzyz;
  StrainSecondDerivative.xzxz=StrainSecondDerivative.zxzx;
  StrainSecondDerivative.yzxz=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.zzxz=StrainSecondDerivative.zzzx;
  StrainSecondDerivative.xzyz=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.xzzz=StrainSecondDerivative.zzzx;

  // use symmetry Cijkl = Cklij
  StrainSecondDerivative.zxxx=StrainSecondDerivative.xxzx;
  StrainSecondDerivative.xyxx=StrainSecondDerivative.xxxy;
  StrainSecondDerivative.yyxx=StrainSecondDerivative.xxyy;
  StrainSecondDerivative.xyzx=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.yzxx=StrainSecondDerivative.xxyz;
  StrainSecondDerivative.zzxx=StrainSecondDerivative.xxzz;
  StrainSecondDerivative.zxyy=StrainSecondDerivative.yyzx;
  StrainSecondDerivative.xyyy=StrainSecondDerivative.yyxy;
  StrainSecondDerivative.yzyy=StrainSecondDerivative.yyyz;
  StrainSecondDerivative.zzyy=StrainSecondDerivative.yyzz;
  StrainSecondDerivative.zxyz=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.zxzz=StrainSecondDerivative.zzzx;
  StrainSecondDerivative.xyyz=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.xyzz=StrainSecondDerivative.zzxy;
  StrainSecondDerivative.yzzz=StrainSecondDerivative.zzyz;


  // restore initial box, positions, and boundary condition
  Box[CurrentSystem]=StoredBox;
  ReplicaBox[CurrentSystem]=StoredReplicaBox;
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();

  BoundaryCondition[CurrentSystem]=StoredBoundaryCondition;
  return StrainSecondDerivative;
}


void TestStressTensorNumerically(void)
{
  int i,f1;
  REAL_MATRIX3x3 strain_derivative;
  char buffer[256];
  FILE *FilePtr;

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet implemented for this function.");
    exit(0);
  }

  mkdir("Numerical",S_IRWXU);

  sprintf(buffer,"Numerical/System_%d",CurrentSystem);
  mkdir(buffer,S_IRWXU);


  sprintf(buffer,"Numerical/System_%d/StrainDerivativeTensor.dat",CurrentSystem);
  FilePtr=fopen(buffer,"w");

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Framework[CurrentSystem].Atoms[f1][i].Velocity.x=0.0;
      Framework[CurrentSystem].Atoms[f1][i].Velocity.y=0.0;
      Framework[CurrentSystem].Atoms[f1][i].Velocity.z=0.0;
    }
  }


  // compute the first and second derivative of the energy with respect to strain analytically
  CalculateForce();
  fprintf(FilePtr,"# energy: %g [K] %g [kj/mol] %g [eV]\n\n",UTotal[CurrentSystem]*ENERGY_TO_KELVIN,
        UTotal[CurrentSystem]*ENERGY_TO_KJ_PER_MOL,UTotal[CurrentSystem]*ENERGY_TO_EV);

  fprintf(FilePtr,"# stress tensor (analytically)\n");
  fprintf(FilePtr,"% 18.10f % 18.10f % 18.10f\n",
    (double)StressTensor[CurrentSystem].ax,(double)StressTensor[CurrentSystem].bx,(double)StressTensor[CurrentSystem].cx);
  fprintf(FilePtr,"% 18.10f % 18.10f % 18.10f\n",
    (double)StressTensor[CurrentSystem].ay,(double)StressTensor[CurrentSystem].by,(double)StressTensor[CurrentSystem].cy);
  fprintf(FilePtr,"% 18.10f % 18.10f % 18.10f\n",
    (double)StressTensor[CurrentSystem].az,(double)StressTensor[CurrentSystem].bz,(double)StressTensor[CurrentSystem].cz);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"# strain derivative tensor (analytically)\n");
  fprintf(FilePtr,"% 18.10f % 18.10f % 18.10f\n",
    (double)StrainDerivativeTensor[CurrentSystem].ax,(double)StrainDerivativeTensor[CurrentSystem].bx,(double)StrainDerivativeTensor[CurrentSystem].cx);
  fprintf(FilePtr,"% 18.10f % 18.10f % 18.10f\n",
    (double)StrainDerivativeTensor[CurrentSystem].ay,(double)StrainDerivativeTensor[CurrentSystem].by,(double)StrainDerivativeTensor[CurrentSystem].cy);
  fprintf(FilePtr,"% 18.10f % 18.10f % 18.10f\n",
    (double)StrainDerivativeTensor[CurrentSystem].az,(double)StrainDerivativeTensor[CurrentSystem].bz,(double)StrainDerivativeTensor[CurrentSystem].cz);
  fprintf(FilePtr,"\n");

  // compute the first derivative of the energy with respect to strain numerically
  strain_derivative=ComputeStrainDerivativeNumerically();

  fprintf(FilePtr,"# strain derivative tensor (numerically)\n");
  fprintf(FilePtr,"% 18.10f % 18.10f % 18.10f\n",
    (double)strain_derivative.ax,(double)strain_derivative.bx,(double)strain_derivative.cx);
  fprintf(FilePtr,"% 18.10f % 18.10f % 18.10f\n",
    (double)strain_derivative.ay,(double)strain_derivative.by,(double)strain_derivative.cy);
  fprintf(FilePtr,"% 18.10f % 18.10f % 18.10f\n",
    (double)strain_derivative.az,(double)strain_derivative.bz,(double)strain_derivative.cz);
  fprintf(FilePtr,"\n");

  fclose(FilePtr);
}

void TestStrainSecondDerivativeNumerically(void)
{
  REAL_MATRIX9x9 strain_second_derivative;
  int StoredComputeBornTerm;
  char buffer[256];
  FILE *FilePtr;

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet implemented for this function.");
    exit(0);
  }

  StoredComputeBornTerm=ComputeBornTerm;

  InitializeMatrix9x9(&strain_second_derivative);

  mkdir("Numerical",S_IRWXU);

  sprintf(buffer,"Numerical/System_%d",CurrentSystem);
  mkdir(buffer,S_IRWXU);

  sprintf(buffer,"Numerical/System_%d/StrainSecondDerivativeTensor.dat",CurrentSystem);
  FilePtr=fopen(buffer,"w");

  ComputeBornTerm=TRUE;
  CalculateForce();
  fprintf(FilePtr,"# energy: %g [K] %g [kj/mol] %g [eV]\n\n",UTotal[CurrentSystem]*ENERGY_TO_KELVIN,
        UTotal[CurrentSystem]*ENERGY_TO_KJ_PER_MOL,UTotal[CurrentSystem]*ENERGY_TO_EV);
  ComputeBornTerm=FALSE;

  // compute the second derivative of the energy with respect to strain numerically
  strain_second_derivative=ComputeStrainSecondDerivativeNumerically();

  fprintf(FilePtr,"# strain second derivative tensor (numerically from analytical first derivative)\n");
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xxxx,(double)strain_second_derivative.xxyy,(double)strain_second_derivative.xxzz,
    (double)strain_second_derivative.xxyz,(double)strain_second_derivative.xxzx,(double)strain_second_derivative.xxxy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.yyxx,(double)strain_second_derivative.yyyy,(double)strain_second_derivative.yyzz,
    (double)strain_second_derivative.yyyz,(double)strain_second_derivative.yyzx,(double)strain_second_derivative.yyxy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.zzxx,(double)strain_second_derivative.zzyy,(double)strain_second_derivative.zzzz,
    (double)strain_second_derivative.zzyz,(double)strain_second_derivative.zzzx,(double)strain_second_derivative.zzxy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.yzxx,(double)strain_second_derivative.yzyy,(double)strain_second_derivative.yzzz,
    (double)strain_second_derivative.yzyz,(double)strain_second_derivative.yzzx,(double)strain_second_derivative.yzxy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.zxxx,(double)strain_second_derivative.zxyy,(double)strain_second_derivative.zxzz,
    (double)strain_second_derivative.zxyz,(double)strain_second_derivative.zxzx,(double)strain_second_derivative.zxxy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xyxx,(double)strain_second_derivative.xyyy,(double)strain_second_derivative.xyzz,
    (double)strain_second_derivative.xyyz,(double)strain_second_derivative.xyzx,(double)strain_second_derivative.xyxy);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"# strain second derivative tensor (analytically)\n");
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xxxx,(double)BornTerm[CurrentSystem].xxyy,(double)BornTerm[CurrentSystem].xxzz,
    (double)BornTerm[CurrentSystem].xxyz,(double)BornTerm[CurrentSystem].xxzx,(double)BornTerm[CurrentSystem].xxxy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].yyxx,(double)BornTerm[CurrentSystem].yyyy,(double)BornTerm[CurrentSystem].yyzz,
    (double)BornTerm[CurrentSystem].yyyz,(double)BornTerm[CurrentSystem].yyzx,(double)BornTerm[CurrentSystem].yyxy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].zzxx,(double)BornTerm[CurrentSystem].zzyy,(double)BornTerm[CurrentSystem].zzzz,
    (double)BornTerm[CurrentSystem].zzyz,(double)BornTerm[CurrentSystem].zzzx,(double)BornTerm[CurrentSystem].zzxy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].yzxx,(double)BornTerm[CurrentSystem].yzyy,(double)BornTerm[CurrentSystem].yzzz,
    (double)BornTerm[CurrentSystem].yzyz,(double)BornTerm[CurrentSystem].yzzx,(double)BornTerm[CurrentSystem].yzxy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].zxxx,(double)BornTerm[CurrentSystem].zxyy,(double)BornTerm[CurrentSystem].zxzz,
    (double)BornTerm[CurrentSystem].zxyz,(double)BornTerm[CurrentSystem].zxzx,(double)BornTerm[CurrentSystem].zxxy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xyxx,(double)BornTerm[CurrentSystem].xyyy,(double)BornTerm[CurrentSystem].xyzz,
    (double)BornTerm[CurrentSystem].xyyz,(double)BornTerm[CurrentSystem].xyzx,(double)BornTerm[CurrentSystem].xyxy);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"# difference strain second derivative tensor\n");
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xxxx-strain_second_derivative.xxxx),(double)(BornTerm[CurrentSystem].xxyy-strain_second_derivative.xxyy),(double)(BornTerm[CurrentSystem].xxzz-strain_second_derivative.xxzz),
    (double)(BornTerm[CurrentSystem].xxyz-strain_second_derivative.xxyz),(double)(BornTerm[CurrentSystem].xxzx-strain_second_derivative.xxzx),(double)(BornTerm[CurrentSystem].xxxy-strain_second_derivative.xxxy));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].yyxx-strain_second_derivative.yyxx),(double)(BornTerm[CurrentSystem].yyyy-strain_second_derivative.yyyy),(double)(BornTerm[CurrentSystem].yyzz-strain_second_derivative.yyzz),
    (double)(BornTerm[CurrentSystem].yyyz-strain_second_derivative.yyyz),(double)(BornTerm[CurrentSystem].yyzx-strain_second_derivative.yyzx),(double)(BornTerm[CurrentSystem].yyxy-strain_second_derivative.yyxy));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].zzxx-strain_second_derivative.zzxx),(double)(BornTerm[CurrentSystem].zzyy-strain_second_derivative.zzyy),(double)(BornTerm[CurrentSystem].zzzz-strain_second_derivative.zzzz),
    (double)(BornTerm[CurrentSystem].zzyz-strain_second_derivative.zzyz),(double)(BornTerm[CurrentSystem].zzzx-strain_second_derivative.zzzx),(double)(BornTerm[CurrentSystem].zzxy-strain_second_derivative.zzxy));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].yzxx-strain_second_derivative.yzxx),(double)(BornTerm[CurrentSystem].yzyy-strain_second_derivative.yzyy),(double)(BornTerm[CurrentSystem].yzzz-strain_second_derivative.yzzz),
    (double)(BornTerm[CurrentSystem].yzyz-strain_second_derivative.yzyz),(double)(BornTerm[CurrentSystem].yzzx-strain_second_derivative.yzzx),(double)(BornTerm[CurrentSystem].yzxy-strain_second_derivative.yzxy));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].zxxx-strain_second_derivative.zxxx),(double)(BornTerm[CurrentSystem].zxyy-strain_second_derivative.zxyy),(double)(BornTerm[CurrentSystem].zxzz-strain_second_derivative.zxzz),
    (double)(BornTerm[CurrentSystem].zxyz-strain_second_derivative.zxyz),(double)(BornTerm[CurrentSystem].zxzx-strain_second_derivative.zxzx),(double)(BornTerm[CurrentSystem].zxxy-strain_second_derivative.zxxy));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xyxx-strain_second_derivative.xyxx),(double)(BornTerm[CurrentSystem].xyyy-strain_second_derivative.xyyy),(double)(BornTerm[CurrentSystem].xyzz-strain_second_derivative.xyzz),
    (double)(BornTerm[CurrentSystem].xyyz-strain_second_derivative.xyyz),(double)(BornTerm[CurrentSystem].xyzx-strain_second_derivative.xyzx),(double)(BornTerm[CurrentSystem].xyxy-strain_second_derivative.xyxy));
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"# strain second derivative tensor (analytically)\n");
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xxxx,(double)BornTerm[CurrentSystem].yxxx,(double)BornTerm[CurrentSystem].zxxx,
    (double)BornTerm[CurrentSystem].xxyx,(double)BornTerm[CurrentSystem].yxyx,(double)BornTerm[CurrentSystem].zxyx,
    (double)BornTerm[CurrentSystem].xxzx,(double)BornTerm[CurrentSystem].yxzx,(double)BornTerm[CurrentSystem].zxzx);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xyxx,(double)BornTerm[CurrentSystem].yyxx,(double)BornTerm[CurrentSystem].zyxx,
    (double)BornTerm[CurrentSystem].xyyx,(double)BornTerm[CurrentSystem].yyyx,(double)BornTerm[CurrentSystem].zyyx,
    (double)BornTerm[CurrentSystem].xyzx,(double)BornTerm[CurrentSystem].yyzx,(double)BornTerm[CurrentSystem].zyzx);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xzxx,(double)BornTerm[CurrentSystem].yzxx,(double)BornTerm[CurrentSystem].zzxx,
    (double)BornTerm[CurrentSystem].xzyx,(double)BornTerm[CurrentSystem].yzyx,(double)BornTerm[CurrentSystem].zzyx,
    (double)BornTerm[CurrentSystem].xzzx,(double)BornTerm[CurrentSystem].yzzx,(double)BornTerm[CurrentSystem].zzzx);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xxxy,(double)BornTerm[CurrentSystem].yxxy,(double)BornTerm[CurrentSystem].zxxy,
    (double)BornTerm[CurrentSystem].xxyy,(double)BornTerm[CurrentSystem].yxyy,(double)BornTerm[CurrentSystem].zxyy,
    (double)BornTerm[CurrentSystem].xxzy,(double)BornTerm[CurrentSystem].yxzy,(double)BornTerm[CurrentSystem].zxzy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xyxy,(double)BornTerm[CurrentSystem].yyxy,(double)BornTerm[CurrentSystem].zyxy,
    (double)BornTerm[CurrentSystem].xyyy,(double)BornTerm[CurrentSystem].yyyy,(double)BornTerm[CurrentSystem].zyyy,
    (double)BornTerm[CurrentSystem].xyzy,(double)BornTerm[CurrentSystem].yyzy,(double)BornTerm[CurrentSystem].zyzy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xzxy,(double)BornTerm[CurrentSystem].yzxy,(double)BornTerm[CurrentSystem].zzxy,
    (double)BornTerm[CurrentSystem].xzyy,(double)BornTerm[CurrentSystem].yzyy,(double)BornTerm[CurrentSystem].zzyy,
    (double)BornTerm[CurrentSystem].xzzy,(double)BornTerm[CurrentSystem].yzzy,(double)BornTerm[CurrentSystem].zzzy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xxxz,(double)BornTerm[CurrentSystem].yxxz,(double)BornTerm[CurrentSystem].zxxz,
    (double)BornTerm[CurrentSystem].xxyz,(double)BornTerm[CurrentSystem].yxyz,(double)BornTerm[CurrentSystem].zxyz,
    (double)BornTerm[CurrentSystem].xxzz,(double)BornTerm[CurrentSystem].yxzz,(double)BornTerm[CurrentSystem].zxzz);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xyxz,(double)BornTerm[CurrentSystem].yyxz,(double)BornTerm[CurrentSystem].zyxz,
    (double)BornTerm[CurrentSystem].xyyz,(double)BornTerm[CurrentSystem].yyyz,(double)BornTerm[CurrentSystem].zyyz,
    (double)BornTerm[CurrentSystem].xyzz,(double)BornTerm[CurrentSystem].yyzz,(double)BornTerm[CurrentSystem].zyzz);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xzxz,(double)BornTerm[CurrentSystem].yzxz,(double)BornTerm[CurrentSystem].zzxz,
    (double)BornTerm[CurrentSystem].xzyz,(double)BornTerm[CurrentSystem].yzyz,(double)BornTerm[CurrentSystem].zzyz,
    (double)BornTerm[CurrentSystem].xzzz,(double)BornTerm[CurrentSystem].yzzz,(double)BornTerm[CurrentSystem].zzzz);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"# strain second derivative tensor (numerically from analytical first derivative)\n");
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xxxx,(double)strain_second_derivative.yxxx,(double)strain_second_derivative.zxxx,
    (double)strain_second_derivative.xxyx,(double)strain_second_derivative.yxyx,(double)strain_second_derivative.zxyx,
    (double)strain_second_derivative.xxzx,(double)strain_second_derivative.yxzx,(double)strain_second_derivative.zxzx);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xyxx,(double)strain_second_derivative.yyxx,(double)strain_second_derivative.zyxx,
    (double)strain_second_derivative.xyyx,(double)strain_second_derivative.yyyx,(double)strain_second_derivative.zyyx,
    (double)strain_second_derivative.xyzx,(double)strain_second_derivative.yyzx,(double)strain_second_derivative.zyzx);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xzxx,(double)strain_second_derivative.yzxx,(double)strain_second_derivative.zzxx,
    (double)strain_second_derivative.xzyx,(double)strain_second_derivative.yzyx,(double)strain_second_derivative.zzyx,
    (double)strain_second_derivative.xzzx,(double)strain_second_derivative.yzzx,(double)strain_second_derivative.zzzx);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xxxy,(double)strain_second_derivative.yxxy,(double)strain_second_derivative.zxxy,
    (double)strain_second_derivative.xxyy,(double)strain_second_derivative.yxyy,(double)strain_second_derivative.zxyy,
    (double)strain_second_derivative.xxzy,(double)strain_second_derivative.yxzy,(double)strain_second_derivative.zxzy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xyxy,(double)strain_second_derivative.yyxy,(double)strain_second_derivative.zyxy,
    (double)strain_second_derivative.xyyy,(double)strain_second_derivative.yyyy,(double)strain_second_derivative.zyyy,
    (double)strain_second_derivative.xyzy,(double)strain_second_derivative.yyzy,(double)strain_second_derivative.zyzy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xzxy,(double)strain_second_derivative.yzxy,(double)strain_second_derivative.zzxy,
    (double)strain_second_derivative.xzyy,(double)strain_second_derivative.yzyy,(double)strain_second_derivative.zzyy,
    (double)strain_second_derivative.xzzy,(double)strain_second_derivative.yzzy,(double)strain_second_derivative.zzzy);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xxxz,(double)strain_second_derivative.yxxz,(double)strain_second_derivative.zxxz,
    (double)strain_second_derivative.xxyz,(double)strain_second_derivative.yxyz,(double)strain_second_derivative.zxyz,
    (double)strain_second_derivative.xxzz,(double)strain_second_derivative.yxzz,(double)strain_second_derivative.zxzz);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xyxz,(double)strain_second_derivative.yyxz,(double)strain_second_derivative.zyxz,
    (double)strain_second_derivative.xyyz,(double)strain_second_derivative.yyyz,(double)strain_second_derivative.zyyz,
    (double)strain_second_derivative.xyzz,(double)strain_second_derivative.yyzz,(double)strain_second_derivative.zyzz);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)strain_second_derivative.xzxz,(double)strain_second_derivative.yzxz,(double)strain_second_derivative.zzxz,
    (double)strain_second_derivative.xzyz,(double)strain_second_derivative.yzyz,(double)strain_second_derivative.zzyz,
    (double)strain_second_derivative.xzzz,(double)strain_second_derivative.yzzz,(double)strain_second_derivative.zzzz);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"# difference strain second derivative tensor\n");
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xxxx-strain_second_derivative.xxxx),(double)(BornTerm[CurrentSystem].yxxx-strain_second_derivative.yxxx),(double)(BornTerm[CurrentSystem].zxxx-strain_second_derivative.zxxx),
    (double)(BornTerm[CurrentSystem].xxyx-strain_second_derivative.xxyx),(double)(BornTerm[CurrentSystem].yxyx-strain_second_derivative.yxyx),(double)(BornTerm[CurrentSystem].zxyx-strain_second_derivative.zxyx),
    (double)(BornTerm[CurrentSystem].xxzx-strain_second_derivative.xxzx),(double)(BornTerm[CurrentSystem].yxzx-strain_second_derivative.yxzx),(double)(BornTerm[CurrentSystem].zxzx-strain_second_derivative.zxzx));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xyxx-strain_second_derivative.xyxx),(double)(BornTerm[CurrentSystem].yyxx-strain_second_derivative.yyxx),(double)(BornTerm[CurrentSystem].zyxx-strain_second_derivative.zyxx),
    (double)(BornTerm[CurrentSystem].xyyx-strain_second_derivative.xyyx),(double)(BornTerm[CurrentSystem].yyyx-strain_second_derivative.yyyx),(double)(BornTerm[CurrentSystem].zyyx-strain_second_derivative.zyyx),
    (double)(BornTerm[CurrentSystem].xyzx-strain_second_derivative.xyzx),(double)(BornTerm[CurrentSystem].yyzx-strain_second_derivative.yyzx),(double)(BornTerm[CurrentSystem].zyzx-strain_second_derivative.zyzx));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xzxx-strain_second_derivative.xzxx),(double)(BornTerm[CurrentSystem].yzxx-strain_second_derivative.yzxx),(double)(BornTerm[CurrentSystem].zzxx-strain_second_derivative.zzxx),
    (double)(BornTerm[CurrentSystem].xzyx-strain_second_derivative.xzyx),(double)(BornTerm[CurrentSystem].yzyx-strain_second_derivative.yzyx),(double)(BornTerm[CurrentSystem].zzyx-strain_second_derivative.zzyx),
    (double)(BornTerm[CurrentSystem].xzzx-strain_second_derivative.xzzx),(double)(BornTerm[CurrentSystem].yzzx-strain_second_derivative.yzzx),(double)(BornTerm[CurrentSystem].zzzx-strain_second_derivative.zzzx));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xxxy-strain_second_derivative.xxxy),(double)(BornTerm[CurrentSystem].yxxy-strain_second_derivative.yxxy),(double)(BornTerm[CurrentSystem].zxxy-strain_second_derivative.zxxy),
    (double)(BornTerm[CurrentSystem].xxyy-strain_second_derivative.xxyy),(double)(BornTerm[CurrentSystem].yxyy-strain_second_derivative.yxyy),(double)(BornTerm[CurrentSystem].zxyy-strain_second_derivative.zxyy),
    (double)(BornTerm[CurrentSystem].xxzy-strain_second_derivative.xxzy),(double)(BornTerm[CurrentSystem].yxzy-strain_second_derivative.yxzy),(double)(BornTerm[CurrentSystem].zxzy-strain_second_derivative.zxzy));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xyxy-strain_second_derivative.xyxy),(double)(BornTerm[CurrentSystem].yyxy-strain_second_derivative.yyxy),(double)(BornTerm[CurrentSystem].zyxy-strain_second_derivative.zyxy),
    (double)(BornTerm[CurrentSystem].xyyy-strain_second_derivative.xyyy),(double)(BornTerm[CurrentSystem].yyyy-strain_second_derivative.yyyy),(double)(BornTerm[CurrentSystem].zyyy-strain_second_derivative.zyyy),
    (double)(BornTerm[CurrentSystem].xyzy-strain_second_derivative.xyzy),(double)(BornTerm[CurrentSystem].yyzy-strain_second_derivative.yyzy),(double)(BornTerm[CurrentSystem].zyzy-strain_second_derivative.zyzy));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xzxy-strain_second_derivative.xzxy),(double)(BornTerm[CurrentSystem].yzxy-strain_second_derivative.yzxy),(double)(BornTerm[CurrentSystem].zzxy-strain_second_derivative.zzxy),
    (double)(BornTerm[CurrentSystem].xzyy-strain_second_derivative.xzyy),(double)(BornTerm[CurrentSystem].yzyy-strain_second_derivative.yzyy),(double)(BornTerm[CurrentSystem].zzyy-strain_second_derivative.zzyy),
    (double)(BornTerm[CurrentSystem].xzzy-strain_second_derivative.xzzy),(double)(BornTerm[CurrentSystem].yzzy-strain_second_derivative.yzzy),(double)(BornTerm[CurrentSystem].zzzy-strain_second_derivative.zzzy));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xxxz-strain_second_derivative.xxxz),(double)(BornTerm[CurrentSystem].yxxz-strain_second_derivative.yxxz),(double)(BornTerm[CurrentSystem].zxxz-strain_second_derivative.zxxz),
    (double)(BornTerm[CurrentSystem].xxyz-strain_second_derivative.xxyz),(double)(BornTerm[CurrentSystem].yxyz-strain_second_derivative.yxyz),(double)(BornTerm[CurrentSystem].zxyz-strain_second_derivative.zxyz),
    (double)(BornTerm[CurrentSystem].xxzz-strain_second_derivative.xxzz),(double)(BornTerm[CurrentSystem].yxzz-strain_second_derivative.yxzz),(double)(BornTerm[CurrentSystem].zxzz-strain_second_derivative.zxzz));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xyxz-strain_second_derivative.xyxz),(double)(BornTerm[CurrentSystem].yyxz-strain_second_derivative.yyxz),(double)(BornTerm[CurrentSystem].zyxz-strain_second_derivative.zyxz),
    (double)(BornTerm[CurrentSystem].xyyz-strain_second_derivative.xyyz),(double)(BornTerm[CurrentSystem].yyyz-strain_second_derivative.yyyz),(double)(BornTerm[CurrentSystem].zyyz-strain_second_derivative.zyyz),
    (double)(BornTerm[CurrentSystem].xyzz-strain_second_derivative.xyzz),(double)(BornTerm[CurrentSystem].yyzz-strain_second_derivative.yyzz),(double)(BornTerm[CurrentSystem].zyzz-strain_second_derivative.zyzz));
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)(BornTerm[CurrentSystem].xzxz-strain_second_derivative.xzxz),(double)(BornTerm[CurrentSystem].yzxz-strain_second_derivative.yzxz),(double)(BornTerm[CurrentSystem].zzxz-strain_second_derivative.zzxz),
    (double)(BornTerm[CurrentSystem].xzyz-strain_second_derivative.xzyz),(double)(BornTerm[CurrentSystem].yzyz-strain_second_derivative.yzyz),(double)(BornTerm[CurrentSystem].zzyz-strain_second_derivative.zzyz),
    (double)(BornTerm[CurrentSystem].xzzz-strain_second_derivative.xzzz),(double)(BornTerm[CurrentSystem].yzzz-strain_second_derivative.yzzz),(double)(BornTerm[CurrentSystem].zzzz-strain_second_derivative.zzzz));
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"# strain second derivative tensor (analytically, same order as Hessian-matrix)\n");
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xxxx,(double)BornTerm[CurrentSystem].xxxy,(double)BornTerm[CurrentSystem].xxxz,
    (double)BornTerm[CurrentSystem].xxyy,(double)BornTerm[CurrentSystem].xxyz,(double)BornTerm[CurrentSystem].xxzz);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xyxx,(double)BornTerm[CurrentSystem].xyxy,(double)BornTerm[CurrentSystem].xyxz,
    (double)BornTerm[CurrentSystem].xyyy,(double)BornTerm[CurrentSystem].xyyz,(double)BornTerm[CurrentSystem].xyzz);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].xzxx,(double)BornTerm[CurrentSystem].xzxy,(double)BornTerm[CurrentSystem].xzxz,
    (double)BornTerm[CurrentSystem].xzyy,(double)BornTerm[CurrentSystem].xzyz,(double)BornTerm[CurrentSystem].xzzz);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].yyxx,(double)BornTerm[CurrentSystem].yyxy,(double)BornTerm[CurrentSystem].yyxz,
    (double)BornTerm[CurrentSystem].yyyy,(double)BornTerm[CurrentSystem].yyyz,(double)BornTerm[CurrentSystem].yyzz);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].yzxx,(double)BornTerm[CurrentSystem].yzxy,(double)BornTerm[CurrentSystem].yzxz,
    (double)BornTerm[CurrentSystem].yzyy,(double)BornTerm[CurrentSystem].yzyz,(double)BornTerm[CurrentSystem].yzzz);
  fprintf(FilePtr,"% 18.4f % 18.4f % 18.4f % 18.4f % 18.4f % 18.4f\n",
    (double)BornTerm[CurrentSystem].zzxx,(double)BornTerm[CurrentSystem].zzxy,(double)BornTerm[CurrentSystem].zzxz,
    (double)BornTerm[CurrentSystem].zzyy,(double)BornTerm[CurrentSystem].zzyz,(double)BornTerm[CurrentSystem].zzzz);
  fprintf(FilePtr,"\n");


  fclose(FilePtr);

  ComputeBornTerm=StoredComputeBornTerm;
}



void TestEnergyForcesHessian(void)
{
  int i,j,m,l,A,f1;
  int MolType,AtomType;
  VECTOR pos,vel,force;
  char buffer[256];
  FILE *FilePtr;

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet implemented for this function.");
    exit(0);
  }

  mkdir("Numerical",S_IRWXU);

  sprintf(buffer,"Numerical/System_%d",CurrentSystem);
  mkdir(buffer,S_IRWXU);

  sprintf(buffer,"Numerical/System_%d/EnergyForces.dat",CurrentSystem);
  FilePtr=fopen(buffer,"w");

  CalculateForce();
  fprintf(FilePtr,"# energy: %18.10f\n",UTotal[CurrentSystem]*ENERGY_TO_KELVIN);

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        fprintf(FilePtr,"Framework Atom: %-5d %-2d Position: %24.20f %24.20f %24.20f   Force: %18.12f %18.12f %18.12f Fixed: %d %d %d\n",
                  j,
                  CurrentFramework,
                  (double)Framework[CurrentSystem].Atoms[f1][j].Position.x,
                  (double)Framework[CurrentSystem].Atoms[f1][j].Position.y,
                  (double)Framework[CurrentSystem].Atoms[f1][j].Position.z,
                  (double)Framework[CurrentSystem].Atoms[f1][j].Force.x,
                  (double)Framework[CurrentSystem].Atoms[f1][j].Force.y,
                  (double)Framework[CurrentSystem].Atoms[f1][j].Force.z,
                  (int)Framework[CurrentSystem].Atoms[f1][j].Fixed.x,
                  (int)Framework[CurrentSystem].Atoms[f1][j].Fixed.y,
                  (int)Framework[CurrentSystem].Atoms[f1][j].Fixed.z);
    }
    fprintf(FilePtr,"\n");
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    MolType=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[MolType].NumberOfGroups;l++)
    {
      if(Components[MolType].Groups[l].Rigid) // rigid unit
      {
        pos=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition;
        vel=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity;
        force=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassForce;

        fprintf(FilePtr,"Molecule: %-5d  Group: %-5d  Position: %18.12f %18.12f %18.12f  Force: %18.12f %18.12f %18.12f\n",
                 m,l,pos.x,pos.y,pos.z,force.x,force.y,force.z);

        for(i=0;i<Components[MolType].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[MolType].Groups[l].Atoms[i];
          AtomType=Adsorbates[CurrentSystem][m].Atoms[A].Type;
          pos=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          vel=Adsorbates[CurrentSystem][m].Atoms[A].Velocity;
          force=Adsorbates[CurrentSystem][m].Atoms[A].Force;

          fprintf(FilePtr,"\tAtom: %-5d  Position: %18.12f %18.12f %18.12f  Force: %18.12f %18.12f %18.12f\n",
                  A,pos.x,pos.y,pos.z,force.x,force.y,force.z);
        }
      }
      else
      {
        for(i=0;i<Components[MolType].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[MolType].Groups[l].Atoms[i];
          AtomType=Adsorbates[CurrentSystem][m].Atoms[A].Type;
          pos=Adsorbates[CurrentSystem][m].Atoms[A].Position;
          vel=Adsorbates[CurrentSystem][m].Atoms[A].Velocity;
          force=Adsorbates[CurrentSystem][m].Atoms[A].Force;

          fprintf(FilePtr,"Molecule: %-5d  Atom: %-5d  Position: %18.12f %18.12f %18.12f  Force: %18.12f %18.12f %18.12f\n",
                  i,A,pos.x,pos.y,pos.z,force.x,force.y,force.z);

        }
      }
    }
  }
  fclose(FilePtr);
}

// HERE
void TestElectricField2(void)
{
  int i,j,m,l;
  VECTOR pos;
  VECTOR der;
  REAL EnergyCentral,EnergyForward2,EnergyForward1,EnergyBackward1,EnergyBackward2;
  REAL delta,charge;

  delta=1e-5;

  CurrentSystem=0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Force.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].Force.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].Force.z=0.0;

      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Force.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].Force.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].Force.z=0.0;

      Cations[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }
  EwaldFourierForce();

  EwaldFourierStaticElectricField();
  CurrentSystem=0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    for(l=0;l<Adsorbates[CurrentSystem][m].NumberOfAtoms;l++)
    {
      charge=Adsorbates[CurrentSystem][m].Atoms[l].Charge;

      // x-direction
      pos=Adsorbates[CurrentSystem][m].Atoms[l].Position;
      EnergyCentral=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      Adsorbates[CurrentSystem][m].Atoms[l].Position.x=pos.x+delta;
      EnergyForward2=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      Adsorbates[CurrentSystem][m].Atoms[l].Position.x=pos.x+0.5*delta;
      EnergyForward1=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      Adsorbates[CurrentSystem][m].Atoms[l].Position.x=pos.x-0.5*delta;
      EnergyBackward1=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      Adsorbates[CurrentSystem][m].Atoms[l].Position.x=pos.x-delta;
      EnergyBackward2=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      der.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
      Adsorbates[CurrentSystem][m].Atoms[l].Position=pos;

      // y-direction
      pos=Adsorbates[CurrentSystem][m].Atoms[l].Position;
      Adsorbates[CurrentSystem][m].Atoms[l].Position.y=pos.y+delta;
      EnergyForward2=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      Adsorbates[CurrentSystem][m].Atoms[l].Position.y=pos.y+0.5*delta;
      EnergyForward1=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      Adsorbates[CurrentSystem][m].Atoms[l].Position.y=pos.y-0.5*delta;
      EnergyBackward1=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      Adsorbates[CurrentSystem][m].Atoms[l].Position.y=pos.y-delta;
      EnergyBackward2=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      der.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
      Adsorbates[CurrentSystem][m].Atoms[l].Position=pos;

      // z-direction
      pos=Adsorbates[CurrentSystem][m].Atoms[l].Position;
      Adsorbates[CurrentSystem][m].Atoms[l].Position.z=pos.z+delta;
      EnergyForward2=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      Adsorbates[CurrentSystem][m].Atoms[l].Position.z=pos.z+0.5*delta;
      EnergyForward1=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      Adsorbates[CurrentSystem][m].Atoms[l].Position.z=pos.z-0.5*delta;
      EnergyBackward1=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      Adsorbates[CurrentSystem][m].Atoms[l].Position.z=pos.z-delta;
      EnergyBackward2=EwaldFourierElectrostaticPotentialAdsorbate(m,l);
      der.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
      Adsorbates[CurrentSystem][m].Atoms[l].Position=pos;

      //charge=1.0;
      fprintf(stderr, "%d %d, electric field numerically: %g %g %g Force/q: %g %g %g Elec: %g %g %g diff: %g %g %g\n",m,l,-der.x,-der.y,-der.z,
         Adsorbates[CurrentSystem][m].Atoms[l].Force.x/charge,
         Adsorbates[CurrentSystem][m].Atoms[l].Force.y/charge,
         Adsorbates[CurrentSystem][m].Atoms[l].Force.z/charge,
         Adsorbates[CurrentSystem][m].Atoms[l].ElectricField.x,Adsorbates[CurrentSystem][m].Atoms[l].ElectricField.y,Adsorbates[CurrentSystem][m].Atoms[l].ElectricField.z,
         fabs(-der.x-Adsorbates[CurrentSystem][m].Atoms[l].ElectricField.x),
         fabs(-der.y-Adsorbates[CurrentSystem][m].Atoms[l].ElectricField.y),
         fabs(-der.z-Adsorbates[CurrentSystem][m].Atoms[l].ElectricField.z));
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    for(l=0;l<Cations[CurrentSystem][m].NumberOfAtoms;l++)
    {
      charge=Cations[CurrentSystem][m].Atoms[l].Charge;

      // x-direction
      pos=Cations[CurrentSystem][m].Atoms[l].Position;
      EnergyCentral=EwaldFourierElectrostaticPotentialCation(m,l);
      Cations[CurrentSystem][m].Atoms[l].Position.x=pos.x+delta;
      EnergyForward2=EwaldFourierElectrostaticPotentialCation(m,l);
      Cations[CurrentSystem][m].Atoms[l].Position.x=pos.x+0.5*delta;
      EnergyForward1=EwaldFourierElectrostaticPotentialCation(m,l);
      Cations[CurrentSystem][m].Atoms[l].Position.x=pos.x-0.5*delta;
      EnergyBackward1=EwaldFourierElectrostaticPotentialCation(m,l);
      Cations[CurrentSystem][m].Atoms[l].Position.x=pos.x-delta;
      EnergyBackward2=EwaldFourierElectrostaticPotentialCation(m,l);
      der.x=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
      Cations[CurrentSystem][m].Atoms[l].Position=pos;

      // y-direction
      pos=Cations[CurrentSystem][m].Atoms[l].Position;
      Cations[CurrentSystem][m].Atoms[l].Position.y=pos.y+delta;
      EnergyForward2=EwaldFourierElectrostaticPotentialCation(m,l);
      Cations[CurrentSystem][m].Atoms[l].Position.y=pos.y+0.5*delta;
      EnergyForward1=EwaldFourierElectrostaticPotentialCation(m,l);
      Cations[CurrentSystem][m].Atoms[l].Position.y=pos.y-0.5*delta;
      EnergyBackward1=EwaldFourierElectrostaticPotentialCation(m,l);
      Cations[CurrentSystem][m].Atoms[l].Position.y=pos.y-delta;
      EnergyBackward2=EwaldFourierElectrostaticPotentialCation(m,l);
      der.y=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
      Cations[CurrentSystem][m].Atoms[l].Position=pos;

      // z-direction
      pos=Cations[CurrentSystem][m].Atoms[l].Position;
      Cations[CurrentSystem][m].Atoms[l].Position.z=pos.z+delta;
      EnergyForward2=EwaldFourierElectrostaticPotentialCation(m,l);
      Cations[CurrentSystem][m].Atoms[l].Position.z=pos.z+0.5*delta;
      EnergyForward1=EwaldFourierElectrostaticPotentialCation(m,l);
      Cations[CurrentSystem][m].Atoms[l].Position.z=pos.z-0.5*delta;
      EnergyBackward1=EwaldFourierElectrostaticPotentialCation(m,l);
      Cations[CurrentSystem][m].Atoms[l].Position.z=pos.z-delta;
      EnergyBackward2=EwaldFourierElectrostaticPotentialCation(m,l);
      der.z=(-EnergyForward2+8.0*EnergyForward1-8.0*EnergyBackward1+EnergyBackward2)/(6.0*delta);
      Cations[CurrentSystem][m].Atoms[l].Position=pos;

      charge=1.0;
      fprintf(stderr, "%d %d, electric field numerically: %g %g %g Force/q: %g %g %g Elec: %g %g %g diff: %g %g %g\n",m,l,-der.x,-der.y,-der.z,
         Cations[CurrentSystem][m].Atoms[l].Force.x/charge,
         Cations[CurrentSystem][m].Atoms[l].Force.y/charge,
         Cations[CurrentSystem][m].Atoms[l].Force.z/charge,
         Cations[CurrentSystem][m].Atoms[l].ElectricField.x,Cations[CurrentSystem][m].Atoms[l].ElectricField.y,Cations[CurrentSystem][m].Atoms[l].ElectricField.z,
         fabs(-der.x-Cations[CurrentSystem][m].Atoms[l].ElectricField.x),
         fabs(-der.y-Cations[CurrentSystem][m].Atoms[l].ElectricField.y),
         fabs(-der.z-Cations[CurrentSystem][m].Atoms[l].ElectricField.z));
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Force.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].Force.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].Force.z=0.0;

      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Force.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].Force.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].Force.z=0.0;

      Cations[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }

  for(l=0;l<Cations[CurrentSystem][NumberOfCationMolecules[CurrentSystem]-1].NumberOfAtoms;l++)
  {
    TrialPosition[CurrentSystem][l]=Cations[CurrentSystem][NumberOfCationMolecules[CurrentSystem]-1].Atoms[l].Position;
    ElectricFieldAtTrialPosition[CurrentSystem][l].x=0.0;
    ElectricFieldAtTrialPosition[CurrentSystem][l].y=0.0;
    ElectricFieldAtTrialPosition[CurrentSystem][l].z=0.0;
  }
  //NumberOfCationMolecules[CurrentSystem]--;

  CurrentComponent=4;
  ComputeStaticElectricFieldEwald(FALSE,-1,10);

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    for(l=0;l<Adsorbates[CurrentSystem][m].NumberOfAtoms;l++)
    {
       fprintf(stderr, "%d %d, electric field: %g %g %g\n",m,l,
          Adsorbates[CurrentSystem][m].Atoms[l].ElectricField.x,
          Adsorbates[CurrentSystem][m].Atoms[l].ElectricField.y,
          Adsorbates[CurrentSystem][m].Atoms[l].ElectricField.z);
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    for(l=0;l<Cations[CurrentSystem][m].NumberOfAtoms;l++)
    {
       fprintf(stderr, "%d %d, electric field: %g %g %g\n",m,l,
          Cations[CurrentSystem][m].Atoms[l].ElectricField.x,
          Cations[CurrentSystem][m].Atoms[l].ElectricField.y,
          Cations[CurrentSystem][m].Atoms[l].ElectricField.z);
    }
  }
  for(l=0;l<Components[CurrentComponent].NumberOfAtoms;l++)
  {
    fprintf(stderr, "trial: %g %g %g\n",ElectricFieldAtTrialPosition[CurrentSystem][l].x,
                               ElectricFieldAtTrialPosition[CurrentSystem][l].y,
                               ElectricFieldAtTrialPosition[CurrentSystem][l].z);
  }

}

void TestElectricField3(void)
{
  CurrentSystem=0;
  NumberOfAdsorbateMolecules[CurrentSystem]=0;
  NumberOfCationMolecules[CurrentSystem]=0;

  Framework[CurrentSystem].NumberOfAtoms[0]=2;

  Framework[CurrentSystem].Atoms[0][0].Position.x=-0.5;
  Framework[CurrentSystem].Atoms[0][0].Position.y=0.0;
  Framework[CurrentSystem].Atoms[0][0].Position.z=0.0;
  Framework[CurrentSystem].Atoms[0][0].ElectricField.x=0.0;
  Framework[CurrentSystem].Atoms[0][0].ElectricField.y=0.0;
  Framework[CurrentSystem].Atoms[0][0].ElectricField.z=0.0;
  Framework[CurrentSystem].Atoms[0][0].InducedDipole.x=1.0;
  Framework[CurrentSystem].Atoms[0][0].InducedDipole.y=0.0;
  Framework[CurrentSystem].Atoms[0][0].InducedDipole.z=0.0;

  Framework[CurrentSystem].Atoms[0][1].Position.x=0.5;
  Framework[CurrentSystem].Atoms[0][1].Position.y=0.0;
  Framework[CurrentSystem].Atoms[0][1].Position.z=0.0;
  Framework[CurrentSystem].Atoms[0][1].ElectricField.x=0.0;
  Framework[CurrentSystem].Atoms[0][1].ElectricField.y=0.0;
  Framework[CurrentSystem].Atoms[0][1].ElectricField.z=0.0;
  Framework[CurrentSystem].Atoms[0][1].InducedDipole.x=1.0;
  Framework[CurrentSystem].Atoms[0][1].InducedDipole.y=0.0;
  Framework[CurrentSystem].Atoms[0][1].InducedDipole.z=0.0;


  NumberOfAdsorbateMolecules[CurrentSystem]=0;
  Adsorbates[CurrentSystem][0].Atoms[0].Position.x=-0.5;
  Adsorbates[CurrentSystem][0].Atoms[0].Position.y=0.0;
  Adsorbates[CurrentSystem][0].Atoms[0].Position.z=0.0;
  Adsorbates[CurrentSystem][0].Atoms[0].ElectricField.x=0.0;
  Adsorbates[CurrentSystem][0].Atoms[0].ElectricField.y=0.0;
  Adsorbates[CurrentSystem][0].Atoms[0].ElectricField.z=0.0;
  Adsorbates[CurrentSystem][0].Atoms[0].InducedDipole.x=1.0;
  Adsorbates[CurrentSystem][0].Atoms[0].InducedDipole.y=0.0;
  Adsorbates[CurrentSystem][0].Atoms[0].InducedDipole.z=0.0;

  Adsorbates[CurrentSystem][1].Atoms[0].Position.x=0.5;
  Adsorbates[CurrentSystem][1].Atoms[0].Position.y=0.0;
  Adsorbates[CurrentSystem][1].Atoms[0].Position.z=0.0;
  Adsorbates[CurrentSystem][1].Atoms[0].ElectricField.x=0.0;
  Adsorbates[CurrentSystem][1].Atoms[0].ElectricField.y=0.0;
  Adsorbates[CurrentSystem][1].Atoms[0].ElectricField.z=0.0;
  Adsorbates[CurrentSystem][1].Atoms[0].InducedDipole.x=1.0;
  Adsorbates[CurrentSystem][1].Atoms[0].InducedDipole.y=0.0;
  Adsorbates[CurrentSystem][1].Atoms[0].InducedDipole.z=0.0;

  NumberOfCationMolecules[CurrentSystem]=2;
  Cations[CurrentSystem][0].Atoms[0].Position.x=-0.5;
  Cations[CurrentSystem][0].Atoms[0].Position.y=1.0;
  Cations[CurrentSystem][0].Atoms[0].Position.z=0.0;
  Cations[CurrentSystem][0].Atoms[0].ElectricField.x=0.0;
  Cations[CurrentSystem][0].Atoms[0].ElectricField.y=0.0;
  Cations[CurrentSystem][0].Atoms[0].ElectricField.z=0.0;
  Cations[CurrentSystem][0].Atoms[0].InducedDipole.x=1.0;
  Cations[CurrentSystem][0].Atoms[0].InducedDipole.y=0.0;
  Cations[CurrentSystem][0].Atoms[0].InducedDipole.z=0.0;

  Cations[CurrentSystem][1].Atoms[0].Position.x=0.5;
  Cations[CurrentSystem][1].Atoms[0].Position.y=1.0;
  Cations[CurrentSystem][1].Atoms[0].Position.z=0.0;
  Cations[CurrentSystem][1].Atoms[0].ElectricField.x=0.0;
  Cations[CurrentSystem][1].Atoms[0].ElectricField.y=0.0;
  Cations[CurrentSystem][1].Atoms[0].ElectricField.z=0.0;
  Cations[CurrentSystem][1].Atoms[0].InducedDipole.x=1.0;
  Cations[CurrentSystem][1].Atoms[0].InducedDipole.y=0.0;
  Cations[CurrentSystem][1].Atoms[0].InducedDipole.z=0.0;


  ComputeElectricFieldFromInducedDipolesEwald();

  fprintf(stderr, "Framework Atom 0: %g %g %g\n",
     Framework[CurrentSystem].Atoms[0][0].ElectricField.x,
     Framework[CurrentSystem].Atoms[0][0].ElectricField.y,
     Framework[CurrentSystem].Atoms[0][0].ElectricField.z);
  fprintf(stderr, "Framework Atom 0: %g %g %g\n",
     Framework[CurrentSystem].Atoms[0][1].ElectricField.x,
     Framework[CurrentSystem].Atoms[0][1].ElectricField.y,
     Framework[CurrentSystem].Atoms[0][1].ElectricField.z);

  fprintf(stderr, "Atom 0: %g %g %g\n",
     Adsorbates[CurrentSystem][0].Atoms[0].ElectricField.x,
     Adsorbates[CurrentSystem][0].Atoms[0].ElectricField.y,
     Adsorbates[CurrentSystem][0].Atoms[0].ElectricField.z);
  fprintf(stderr, "Atom 0: %g %g %g\n\n",
     Adsorbates[CurrentSystem][1].Atoms[0].ElectricField.x,
     Adsorbates[CurrentSystem][1].Atoms[0].ElectricField.y,
     Adsorbates[CurrentSystem][1].Atoms[0].ElectricField.z);

  fprintf(stderr, "Atom 0: %g %g %g\n",
     Cations[CurrentSystem][0].Atoms[0].ElectricField.x,
     Cations[CurrentSystem][0].Atoms[0].ElectricField.y,
     Cations[CurrentSystem][0].Atoms[0].ElectricField.z);
  fprintf(stderr, "Atom 0: %g %g %g\n",
     Cations[CurrentSystem][1].Atoms[0].ElectricField.x,
     Cations[CurrentSystem][1].Atoms[0].ElectricField.y,
     Cations[CurrentSystem][1].Atoms[0].ElectricField.z);



  exit(0);
}

void CheckStatusNumerically(void)
{

/*
  // randomize the framework positions
  int i,f1;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Framework[CurrentSystem].Atoms[f1][i].Position.x+=0.1*(0.5-RandomNumber());
      Framework[CurrentSystem].Atoms[f1][i].Position.y+=0.1*(0.5-RandomNumber());
      Framework[CurrentSystem].Atoms[f1][i].Position.z+=0.1*(0.5-RandomNumber());
    }
*/

  ComputeBornTerm=FALSE;
  BoundaryCondition[CurrentSystem]=TRICLINIC;

  CurrentSystem=0;
  StoredBox=Box[CurrentSystem];
  StoredReplicaBox=ReplicaBox[CurrentSystem];
  StoredInverseBox=InverseBox[CurrentSystem];

  AllocateMinimizationLocalMemory();

  // open output-file for systems
  OpenOutputFile();

  // print simulation settings to the output-file
  PrintPreSimulationStatus();

  // compute initial energy/forces
  InitializeForcesAllSystems();

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    TestEnergyForcesHessian();

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    TestForcesNumerically();

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    TestGradientsNumerically();

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    TestStressTensorNumerically();

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    TestStrainSecondDerivativeNumerically();

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    TestHessianNumerically();

  // recompute the energy
  InitializeForcesAllSystems();

  PrintPostSimulationStatus();

  CloseOutputFile();
}

void ComputeCrossTermNumerically(REAL_MATRIX CrossTerm)
{
  int i,j,n,index;
  REAL_MATRIX3x3 StrainDerivativeForward1,StrainDerivativeForward2;
  REAL_MATRIX3x3 StrainDerivativeBackward1,StrainDerivativeBackward2;
  REAL_MATRIX3x3 StrainSecondDerivativeX,StrainSecondDerivativeY,StrainSecondDerivativeZ;
  const REAL delta=1e-5;
  VECTOR pos;
  int type;

  n=CrossTerm.m;
  for(i=0;i<n;i++)
  {
    CrossTerm.element[i][0]=CrossTerm.element[i][1]=CrossTerm.element[i][2]=0.0;
    CrossTerm.element[i][3]=CrossTerm.element[i][4]=CrossTerm.element[i][5]=0.0;
    CrossTerm.element[i][6]=CrossTerm.element[i][7]=CrossTerm.element[i][8]=0.0;
  }

  index=0;
  for(i=0;i<Framework[0].NumberOfAtoms[0];i++)
  {
    pos=Framework[CurrentSystem].Atoms[0][i].Position;
    type=Framework[CurrentSystem].Atoms[0][i].Type;

    // x-direction
    Framework[CurrentSystem].Atoms[0][i].Position.x=pos.x+delta;
    CalculateForce();
    StrainDerivativeForward2=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.x=pos.x+0.5*delta;
    CalculateForce();
    StrainDerivativeForward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.x=pos.x-0.5*delta;
    CalculateForce();
    StrainDerivativeBackward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.x=pos.x-delta;
    CalculateForce();
    StrainDerivativeBackward2=StrainDerivativeTensor[0];

    StrainSecondDerivativeX.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
    StrainSecondDerivativeX.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
    StrainSecondDerivativeX.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
    StrainSecondDerivativeX.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
    StrainSecondDerivativeX.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivativeX.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
    StrainSecondDerivativeX.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivativeX.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
    StrainSecondDerivativeX.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

    // restore original position
    Framework[CurrentSystem].Atoms[0][i].Position.x=pos.x;

    // y-direction
    Framework[CurrentSystem].Atoms[0][i].Position.y=pos.y+delta;
    CalculateForce();
    StrainDerivativeForward2=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.y=pos.y+0.5*delta;
    CalculateForce();
    StrainDerivativeForward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.y=pos.y-0.5*delta;
    CalculateForce();
    StrainDerivativeBackward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.y=pos.y-delta;
    CalculateForce();
    StrainDerivativeBackward2=StrainDerivativeTensor[0];

    StrainSecondDerivativeY.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
    StrainSecondDerivativeY.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
    StrainSecondDerivativeY.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
    StrainSecondDerivativeY.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
    StrainSecondDerivativeY.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivativeY.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
    StrainSecondDerivativeY.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivativeY.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
    StrainSecondDerivativeY.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

    // restore original position
    Framework[CurrentSystem].Atoms[0][i].Position.y=pos.y;

    // z-direction
    Framework[CurrentSystem].Atoms[0][i].Position.z=pos.z+delta;
    CalculateForce();
    StrainDerivativeForward2=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.z=pos.z+0.5*delta;
    CalculateForce();
    StrainDerivativeForward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.z=pos.z-0.5*delta;
    CalculateForce();
    StrainDerivativeBackward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.z=pos.z-delta;
    CalculateForce();
    StrainDerivativeBackward2=StrainDerivativeTensor[0];

    StrainSecondDerivativeZ.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
    StrainSecondDerivativeZ.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
    StrainSecondDerivativeZ.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
    StrainSecondDerivativeZ.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
    StrainSecondDerivativeZ.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivativeZ.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
    StrainSecondDerivativeZ.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivativeZ.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
    StrainSecondDerivativeZ.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

    // restore original position
    Framework[CurrentSystem].Atoms[0][i].Position.z=pos.z;

    CrossTerm.element[3*index][0]=StrainSecondDerivativeX.ax; CrossTerm.element[3*index][3]=StrainSecondDerivativeY.ax; CrossTerm.element[3*index][6]=StrainSecondDerivativeZ.ax;
    CrossTerm.element[3*index][1]=StrainSecondDerivativeX.bx; CrossTerm.element[3*index][4]=StrainSecondDerivativeY.bx; CrossTerm.element[3*index][7]=StrainSecondDerivativeZ.bx;
    CrossTerm.element[3*index][2]=StrainSecondDerivativeX.cx; CrossTerm.element[3*index][5]=StrainSecondDerivativeY.cx; CrossTerm.element[3*index][8]=StrainSecondDerivativeZ.cx;

    CrossTerm.element[3*index+1][0]=StrainSecondDerivativeX.ay; CrossTerm.element[3*index+1][3]=StrainSecondDerivativeY.ay; CrossTerm.element[3*index+1][6]=StrainSecondDerivativeZ.ay;
    CrossTerm.element[3*index+1][1]=StrainSecondDerivativeX.by; CrossTerm.element[3*index+1][4]=StrainSecondDerivativeY.by; CrossTerm.element[3*index+1][7]=StrainSecondDerivativeZ.by;
    CrossTerm.element[3*index+1][2]=StrainSecondDerivativeX.cy; CrossTerm.element[3*index+1][5]=StrainSecondDerivativeY.cy; CrossTerm.element[3*index+1][8]=StrainSecondDerivativeZ.cy;

    CrossTerm.element[3*index+2][0]=StrainSecondDerivativeX.az; CrossTerm.element[3*index+2][3]=StrainSecondDerivativeY.az; CrossTerm.element[3*index+2][6]=StrainSecondDerivativeZ.az;
    CrossTerm.element[3*index+2][1]=StrainSecondDerivativeX.bz; CrossTerm.element[3*index+2][4]=StrainSecondDerivativeY.bz; CrossTerm.element[3*index+2][7]=StrainSecondDerivativeZ.bz;
    CrossTerm.element[3*index+2][2]=StrainSecondDerivativeX.cz; CrossTerm.element[3*index+2][5]=StrainSecondDerivativeY.cz; CrossTerm.element[3*index+2][8]=StrainSecondDerivativeZ.cz;

    index++;
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      pos=Adsorbates[CurrentSystem][i].Atoms[j].Position;
      type=Adsorbates[CurrentSystem][i].Atoms[j].Type;

      // x-direction
      Adsorbates[CurrentSystem][i].Atoms[j].Position.x=pos.x+delta;
      CalculateForce();
      StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

      Adsorbates[CurrentSystem][i].Atoms[j].Position.x=pos.x+0.5*delta;
      CalculateForce();
      StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

      Adsorbates[CurrentSystem][i].Atoms[j].Position.x=pos.x-0.5*delta;
      CalculateForce();
      StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

      Adsorbates[CurrentSystem][i].Atoms[j].Position.x=pos.x-delta;
      CalculateForce();
      StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

      StrainSecondDerivativeX.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
      StrainSecondDerivativeX.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
      StrainSecondDerivativeX.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
      StrainSecondDerivativeX.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
      StrainSecondDerivativeX.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
      StrainSecondDerivativeX.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
      StrainSecondDerivativeX.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
      StrainSecondDerivativeX.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
      StrainSecondDerivativeX.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

      // restore original position
      Adsorbates[CurrentSystem][i].Atoms[j].Position.x=pos.x;

      // y-direction
      Adsorbates[CurrentSystem][i].Atoms[j].Position.y=pos.y+delta;
      CalculateForce();
      StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

      Adsorbates[CurrentSystem][i].Atoms[j].Position.y=pos.y+0.5*delta;
      CalculateForce();
      StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

      Adsorbates[CurrentSystem][i].Atoms[j].Position.y=pos.y-0.5*delta;
      CalculateForce();
      StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

      Adsorbates[CurrentSystem][i].Atoms[j].Position.y=pos.y-delta;
      CalculateForce();
      StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

      StrainSecondDerivativeY.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
      StrainSecondDerivativeY.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
      StrainSecondDerivativeY.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
      StrainSecondDerivativeY.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
      StrainSecondDerivativeY.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
      StrainSecondDerivativeY.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
      StrainSecondDerivativeY.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
      StrainSecondDerivativeY.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
      StrainSecondDerivativeY.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

      // restore original position
      Adsorbates[CurrentSystem][i].Atoms[j].Position.y=pos.y;

      // z-direction
      Adsorbates[CurrentSystem][i].Atoms[j].Position.z=pos.z+delta;
      CalculateForce();
      StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

      Adsorbates[CurrentSystem][i].Atoms[j].Position.z=pos.z+0.5*delta;
      CalculateForce();
      StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

      Adsorbates[CurrentSystem][i].Atoms[j].Position.z=pos.z-0.5*delta;
      CalculateForce();
      StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

      Adsorbates[CurrentSystem][i].Atoms[j].Position.z=pos.z-delta;
      CalculateForce();
      StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

      StrainSecondDerivativeZ.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
      StrainSecondDerivativeZ.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
      StrainSecondDerivativeZ.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
      StrainSecondDerivativeZ.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
      StrainSecondDerivativeZ.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
      StrainSecondDerivativeZ.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
      StrainSecondDerivativeZ.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
      StrainSecondDerivativeZ.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
      StrainSecondDerivativeZ.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

      // restore original position
      Adsorbates[CurrentSystem][i].Atoms[j].Position.z=pos.z;

      CrossTerm.element[3*index][0]=StrainSecondDerivativeX.ax; CrossTerm.element[3*index][3]=StrainSecondDerivativeY.ax; CrossTerm.element[3*index][6]=StrainSecondDerivativeZ.ax;
      CrossTerm.element[3*index][1]=StrainSecondDerivativeX.bx; CrossTerm.element[3*index][4]=StrainSecondDerivativeY.bx; CrossTerm.element[3*index][7]=StrainSecondDerivativeZ.bx;
      CrossTerm.element[3*index][2]=StrainSecondDerivativeX.cx; CrossTerm.element[3*index][5]=StrainSecondDerivativeY.cx; CrossTerm.element[3*index][8]=StrainSecondDerivativeZ.cx;

      CrossTerm.element[3*index+1][0]=StrainSecondDerivativeX.ay; CrossTerm.element[3*index+1][3]=StrainSecondDerivativeY.ay; CrossTerm.element[3*index+1][6]=StrainSecondDerivativeZ.ay;
      CrossTerm.element[3*index+1][1]=StrainSecondDerivativeX.by; CrossTerm.element[3*index+1][4]=StrainSecondDerivativeY.by; CrossTerm.element[3*index+1][7]=StrainSecondDerivativeZ.by;
      CrossTerm.element[3*index+1][2]=StrainSecondDerivativeX.cy; CrossTerm.element[3*index+1][5]=StrainSecondDerivativeY.cy; CrossTerm.element[3*index+1][8]=StrainSecondDerivativeZ.cy;

      CrossTerm.element[3*index+2][0]=StrainSecondDerivativeX.az; CrossTerm.element[3*index+2][3]=StrainSecondDerivativeY.az; CrossTerm.element[3*index+2][6]=StrainSecondDerivativeZ.az;
      CrossTerm.element[3*index+2][1]=StrainSecondDerivativeX.bz; CrossTerm.element[3*index+2][4]=StrainSecondDerivativeY.bz; CrossTerm.element[3*index+2][7]=StrainSecondDerivativeZ.bz;
      CrossTerm.element[3*index+2][2]=StrainSecondDerivativeX.cz; CrossTerm.element[3*index+2][5]=StrainSecondDerivativeY.cz; CrossTerm.element[3*index+2][8]=StrainSecondDerivativeZ.cz;

      index++;
    }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      pos=Cations[CurrentSystem][i].Atoms[j].Position;
      type=Cations[CurrentSystem][i].Atoms[j].Type;

      // x-direction
      Cations[CurrentSystem][i].Atoms[j].Position.x=pos.x+delta;
      CalculateForce();
      StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

      Cations[CurrentSystem][i].Atoms[j].Position.x=pos.x+0.5*delta;
      CalculateForce();
      StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

      Cations[CurrentSystem][i].Atoms[j].Position.x=pos.x-0.5*delta;
      CalculateForce();
      StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

      Cations[CurrentSystem][i].Atoms[j].Position.x=pos.x-delta;
      CalculateForce();
      StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

      StrainSecondDerivativeX.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
      StrainSecondDerivativeX.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
      StrainSecondDerivativeX.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
      StrainSecondDerivativeX.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
      StrainSecondDerivativeX.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
      StrainSecondDerivativeX.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
      StrainSecondDerivativeX.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
      StrainSecondDerivativeX.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
      StrainSecondDerivativeX.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

      // restore original position
      Cations[CurrentSystem][i].Atoms[j].Position.x=pos.x;

      // y-direction
      Cations[CurrentSystem][i].Atoms[j].Position.y=pos.y+delta;
      CalculateForce();
      StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

      Cations[CurrentSystem][i].Atoms[j].Position.y=pos.y+0.5*delta;
      CalculateForce();
      StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

      Cations[CurrentSystem][i].Atoms[j].Position.y=pos.y-0.5*delta;
      CalculateForce();
      StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

      Cations[CurrentSystem][i].Atoms[j].Position.y=pos.y-delta;
      CalculateForce();
      StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

      StrainSecondDerivativeY.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
      StrainSecondDerivativeY.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
      StrainSecondDerivativeY.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
      StrainSecondDerivativeY.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
      StrainSecondDerivativeY.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
      StrainSecondDerivativeY.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
      StrainSecondDerivativeY.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
      StrainSecondDerivativeY.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
      StrainSecondDerivativeY.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

      // restore original position
      Cations[CurrentSystem][i].Atoms[j].Position.y=pos.y;

      // z-direction
      Cations[CurrentSystem][i].Atoms[j].Position.z=pos.z+delta;
      CalculateForce();
      StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

      Cations[CurrentSystem][i].Atoms[j].Position.z=pos.z+0.5*delta;
      CalculateForce();
      StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

      Cations[CurrentSystem][i].Atoms[j].Position.z=pos.z-0.5*delta;
      CalculateForce();
      StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

      Cations[CurrentSystem][i].Atoms[j].Position.z=pos.z-delta;
      CalculateForce();
      StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

      StrainSecondDerivativeZ.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
      StrainSecondDerivativeZ.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
      StrainSecondDerivativeZ.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
      StrainSecondDerivativeZ.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
      StrainSecondDerivativeZ.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
      StrainSecondDerivativeZ.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
      StrainSecondDerivativeZ.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
      StrainSecondDerivativeZ.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
      StrainSecondDerivativeZ.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

      // restore original position
      Cations[CurrentSystem][i].Atoms[j].Position.z=pos.z;

      CrossTerm.element[3*index][0]=StrainSecondDerivativeX.ax; CrossTerm.element[3*index][3]=StrainSecondDerivativeY.ax; CrossTerm.element[3*index][6]=StrainSecondDerivativeZ.ax;
      CrossTerm.element[3*index][1]=StrainSecondDerivativeX.bx; CrossTerm.element[3*index][4]=StrainSecondDerivativeY.bx; CrossTerm.element[3*index][7]=StrainSecondDerivativeZ.bx;
      CrossTerm.element[3*index][2]=StrainSecondDerivativeX.cx; CrossTerm.element[3*index][5]=StrainSecondDerivativeY.cx; CrossTerm.element[3*index][8]=StrainSecondDerivativeZ.cx;

      CrossTerm.element[3*index+1][0]=StrainSecondDerivativeX.ay; CrossTerm.element[3*index+1][3]=StrainSecondDerivativeY.ay; CrossTerm.element[3*index+1][6]=StrainSecondDerivativeZ.ay;
      CrossTerm.element[3*index+1][1]=StrainSecondDerivativeX.by; CrossTerm.element[3*index+1][4]=StrainSecondDerivativeY.by; CrossTerm.element[3*index+1][7]=StrainSecondDerivativeZ.by;
      CrossTerm.element[3*index+1][2]=StrainSecondDerivativeX.cy; CrossTerm.element[3*index+1][5]=StrainSecondDerivativeY.cy; CrossTerm.element[3*index+1][8]=StrainSecondDerivativeZ.cy;

      CrossTerm.element[3*index+2][0]=StrainSecondDerivativeX.az; CrossTerm.element[3*index+2][3]=StrainSecondDerivativeY.az; CrossTerm.element[3*index+2][6]=StrainSecondDerivativeZ.az;
      CrossTerm.element[3*index+2][1]=StrainSecondDerivativeX.bz; CrossTerm.element[3*index+2][4]=StrainSecondDerivativeY.bz; CrossTerm.element[3*index+2][7]=StrainSecondDerivativeZ.bz;
      CrossTerm.element[3*index+2][2]=StrainSecondDerivativeX.cz; CrossTerm.element[3*index+2][5]=StrainSecondDerivativeY.cz; CrossTerm.element[3*index+2][8]=StrainSecondDerivativeZ.cz;

      index++;
    }
}


void ConstrucGeneralizedHessianMatrix(REAL_MATRIX GeneralizedHessianMatrix)
{
/*
  int i,j,n;
  REAL Energy,*Gradients;
  REAL_MATRIX3x3 StrainDerivativeTensor;
  REAL_MATRIX HessianMatrix,CrossTerm,BornMatrix;

  n=GeneralizedHessianMatrix.m;

  CrossTerm=CreateRealMatrix(n,9);
  Gradients=(REAL*)calloc(n,sizeof(REAL));

  HessianMatrix=CreateRealMatrix(n,n);
  BornMatrix=CreateRealMatrix(9,9);

  ComputeCrossTerm=TRUE;
  ComputeEnergyGradientHessian(&Energy,Gradients,&StrainDerivativeTensor,HessianMatrix,CrossTerm);

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      GeneralizedHessianMatrix.element[i][j]=HessianMatrix.element[i][j];

  ComputeCrossTermNumerically(CrossTerm);

  for(i=0;i<n;i++)
    for(j=0;j<9;j++)
      GeneralizedHessianMatrix.element[i][j+n]=GeneralizedHessianMatrix.element[j+n][i]=CrossTerm.element[i][j];

  ComputeBornTerm=TRUE;
  CalculateForce();
  Convert9x9ToRealMatrix(&BornTerm[CurrentSystem],BornMatrix);

  for(i=0;i<9;i++)
    for(j=0;j<9;j++)
      GeneralizedHessianMatrix.element[i+n][j+n]=BornMatrix.element[i][j];

  free(Gradients);
  DeleteRealMatrix(BornMatrix);
  DeleteRealMatrix(CrossTerm);
  DeleteRealMatrix(HessianMatrix);
*/
}

void ComputeNormalModeDerivativeNumerically(REAL_MATRIX GeneralizedHessianMatrix,REAL *Eigenvalues,REAL *Positions,
               CUBIC_SPLINE *Splines,REAL_MATRIX3x3 StoredBox,REAL_MATRIX3x3 StoredInverseBox)
{
  int i,n,mode,size,frames;
  REAL factor,amplitude;
  REAL *displacement;
  REAL frequency,temperature,ReferenceEnergy;
  REAL *x,*y;

  size=200;

  n=GeneralizedHessianMatrix.n;
  displacement=(REAL*)calloc(n,sizeof(REAL));
  x=(REAL*)calloc(size,sizeof(REAL));
  y=(REAL*)calloc(size,sizeof(REAL));

  for(i=0;i<n;i++)
    displacement[i]=0.0;

  //StoredPositionsToRealPositions(Positions,displacement,GeneralizedHessianMatrix,StoredBox,StoredInverseBox);
  CalculateForce();
  ReferenceEnergy=UTotal[CurrentSystem];

  fprintf(stderr, "Creating normal mode derivatives\n");
  for(mode=6;mode<27;mode++)
  {
    frequency=SIGN(sqrt(fabs(Eigenvalues[mode]))*TO_WAVENUMBERS,Eigenvalues[mode]);
    temperature=1000.0;
    amplitude=sqrt(2.0*MOLAR_GAS_CONSTANT*temperature)*TO_WAVENUMBERS/fabs(frequency);

    for(frames=0;frames<size;frames++)
    {
      factor=-amplitude+2.0*amplitude*frames/size;

      for(i=0;i<n;i++)
        displacement[i]=factor*GeneralizedHessianMatrix.element[mode][i];

      //StoredPositionsToRealPositions(Positions,displacement,GeneralizedHessianMatrix,StoredBox,StoredInverseBox);
      CalculateForce();

      x[frames]=factor;
      y[frames]=UTotal[CurrentSystem]-ReferenceEnergy;
    }
    fprintf(stderr, "Creating spline for mode: %d\n",mode);
    Splines[mode]=CreateCubicSpline(size,x,y,SECOND_DERIVATIVE,0,0);
  }
  fprintf(stderr, "Done (Creating normal mode derivatives)\n");
  free(displacement);
}

// eventually need to be implemented anaytically
void CalculateStrainDerivativeOfNumericalParts(void)
{
  int m;

  StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
  StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
  StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    CalculateFrameworkInversionBendForce();
    CalculateFrameworkBondBondForce();
    CalculateFrameworkBondBendForce();
    CalculateFrameworkBendBendForce();
    CalculateFrameworkBondTorsionForce();
    CalculateFrameworkBendTorsionForce();
  }

  // contribution of the intra molecular interactions of the adsorbates
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    CalculateAdsorbateInversionBendForce(m);
    CalculateAdsorbateBondBondForce(m);
    CalculateAdsorbateBondBendForce(m);
    CalculateAdsorbateBendBendForce(m);
    CalculateAdsorbateBondTorsionForce(m);
    CalculateAdsorbateBendTorsionForce(m);
  }

  // contribution of the intra molecular interactions of the cations
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    CalculateCationInversionBendForce(m);
    CalculateCationBondBondForce(m);
    CalculateCationBondBendForce(m);
    CalculateCationBendBendForce(m);
    CalculateCationBondTorsionForce(m);
    CalculateCationBendTorsionForce(m);
  }
}

void AddRemainderOfCrossTermNumerically(REAL_MATRIX HessianMatrix)
{
  int i,n;
  INT_VECTOR3 index;
  REAL_MATRIX3x3 StrainDerivativeForward1,StrainDerivativeForward2;
  REAL_MATRIX3x3 StrainDerivativeBackward1,StrainDerivativeBackward2;
  REAL_MATRIX3x3 StrainSecondDerivativeX,StrainSecondDerivativeY,StrainSecondDerivativeZ;
  const REAL delta=1e-6;
  VECTOR pos;
  int type;

  for(i=0;i<Framework[0].NumberOfAtoms[0];i++)
  {
    pos=Framework[CurrentSystem].Atoms[0][i].Position;
    type=Framework[CurrentSystem].Atoms[0][i].Type;
    index=Framework[CurrentSystem].Atoms[0][i].HessianIndex;

    // x-direction
    Framework[CurrentSystem].Atoms[0][i].Position.x=pos.x+delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeForward2=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.x=pos.x+0.5*delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeForward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.x=pos.x-0.5*delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeBackward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.x=pos.x-delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeBackward2=StrainDerivativeTensor[0];

    StrainSecondDerivativeX.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
    StrainSecondDerivativeX.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
    StrainSecondDerivativeX.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
    StrainSecondDerivativeX.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
    StrainSecondDerivativeX.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivativeX.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
    StrainSecondDerivativeX.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivativeX.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
    StrainSecondDerivativeX.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

    // restore original position
    Framework[CurrentSystem].Atoms[0][i].Position.x=pos.x;

    // y-direction
    Framework[CurrentSystem].Atoms[0][i].Position.y=pos.y+delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeForward2=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.y=pos.y+0.5*delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeForward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.y=pos.y-0.5*delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeBackward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.y=pos.y-delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeBackward2=StrainDerivativeTensor[0];

    StrainSecondDerivativeY.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
    StrainSecondDerivativeY.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
    StrainSecondDerivativeY.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
    StrainSecondDerivativeY.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
    StrainSecondDerivativeY.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivativeY.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
    StrainSecondDerivativeY.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivativeY.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
    StrainSecondDerivativeY.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

    // restore original position
    Framework[CurrentSystem].Atoms[0][i].Position.y=pos.y;

    // z-direction
    Framework[CurrentSystem].Atoms[0][i].Position.z=pos.z+delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeForward2=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.z=pos.z+0.5*delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeForward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.z=pos.z-0.5*delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeBackward1=StrainDerivativeTensor[0];

    Framework[CurrentSystem].Atoms[0][i].Position.z=pos.z-delta;
    CalculateStrainDerivativeOfNumericalParts();
    StrainDerivativeBackward2=StrainDerivativeTensor[0];

    StrainSecondDerivativeZ.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
    StrainSecondDerivativeZ.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
    StrainSecondDerivativeZ.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
    StrainSecondDerivativeZ.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
    StrainSecondDerivativeZ.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivativeZ.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
    StrainSecondDerivativeZ.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivativeZ.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
    StrainSecondDerivativeZ.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

    // restore original position
    Framework[CurrentSystem].Atoms[0][i].Position.z=pos.z;

    n=NumberOfCoordinatesMinimizationVariables;
    switch(Ensemble[CurrentSystem])
    {
      case NPT:
        HessianMatrix.element[index.x][n]+=StrainSecondDerivativeX.ax;
        HessianMatrix.element[index.y][n]+=StrainSecondDerivativeY.ax;
        HessianMatrix.element[index.z][n]+=StrainSecondDerivativeZ.ax;
        break;
      case NPTPR:
      case NPHPR:
        switch(NPTPRCellType[CurrentSystem])
        {
          case ISOTROPIC:
          case ANISOTROPIC:
            HessianMatrix.element[index.x][n]+=StrainSecondDerivativeX.ax;
            HessianMatrix.element[index.x][n+1]+=StrainSecondDerivativeX.by;
            HessianMatrix.element[index.x][n+2]+=StrainSecondDerivativeX.cz;
            HessianMatrix.element[index.y][n]+=StrainSecondDerivativeY.ax;
            HessianMatrix.element[index.y][n+1]+=StrainSecondDerivativeY.by;
            HessianMatrix.element[index.y][n+2]+=StrainSecondDerivativeY.cz;
            HessianMatrix.element[index.z][n]+=StrainSecondDerivativeZ.ax;
            HessianMatrix.element[index.z][n+1]+=StrainSecondDerivativeZ.by;
            HessianMatrix.element[index.z][n+2]+=StrainSecondDerivativeZ.cz;
            break;
            break;
          case REGULAR:
          case REGULAR_UPPER_TRIANGLE:
            HessianMatrix.element[index.x][n]+=StrainSecondDerivativeX.ax;
            HessianMatrix.element[index.x][n+1]+=StrainSecondDerivativeX.ay;
            HessianMatrix.element[index.x][n+2]+=StrainSecondDerivativeX.az;
            HessianMatrix.element[index.x][n+3]+=StrainSecondDerivativeX.by;
            HessianMatrix.element[index.x][n+4]+=StrainSecondDerivativeX.bz;
            HessianMatrix.element[index.x][n+5]+=StrainSecondDerivativeX.cz;
            HessianMatrix.element[index.y][n]+=StrainSecondDerivativeY.ax;
            HessianMatrix.element[index.y][n+1]+=StrainSecondDerivativeY.ay;
            HessianMatrix.element[index.y][n+2]+=StrainSecondDerivativeY.az;
            HessianMatrix.element[index.y][n+3]+=StrainSecondDerivativeY.by;
            HessianMatrix.element[index.y][n+4]+=StrainSecondDerivativeY.bz;
            HessianMatrix.element[index.y][n+5]+=StrainSecondDerivativeY.cz;
            HessianMatrix.element[index.z][n]+=StrainSecondDerivativeZ.ax;
            HessianMatrix.element[index.z][n+1]+=StrainSecondDerivativeZ.ay;
            HessianMatrix.element[index.z][n+2]+=StrainSecondDerivativeZ.az;
            HessianMatrix.element[index.z][n+3]+=StrainSecondDerivativeZ.by;
            HessianMatrix.element[index.z][n+4]+=StrainSecondDerivativeZ.bz;
            HessianMatrix.element[index.z][n+5]+=StrainSecondDerivativeZ.cz;
            break;
          case MONOCLINIC:
          case MONOCLINIC_UPPER_TRIANGLE:
            switch(MonoclinicAngleType[CurrentSystem])
            {
              case MONOCLINIC_ALPHA_ANGLE:
                HessianMatrix.element[index.x][n]+=StrainSecondDerivativeX.ax;
                HessianMatrix.element[index.x][n+1]+=StrainSecondDerivativeX.by;
                HessianMatrix.element[index.x][n+2]+=StrainSecondDerivativeX.bz;
                HessianMatrix.element[index.x][n+3]+=StrainSecondDerivativeX.cz;
                HessianMatrix.element[index.y][n]+=StrainSecondDerivativeY.ax;
                HessianMatrix.element[index.y][n+1]+=StrainSecondDerivativeY.by;
                HessianMatrix.element[index.y][n+2]+=StrainSecondDerivativeY.bz;
                HessianMatrix.element[index.y][n+3]+=StrainSecondDerivativeY.cz;
                HessianMatrix.element[index.z][n]+=StrainSecondDerivativeZ.ax;
                HessianMatrix.element[index.z][n+1]+=StrainSecondDerivativeZ.by;
                HessianMatrix.element[index.z][n+2]+=StrainSecondDerivativeZ.bz;
                HessianMatrix.element[index.z][n+3]+=StrainSecondDerivativeZ.cz;
                break;
              case MONOCLINIC_BETA_ANGLE:
                HessianMatrix.element[index.x][n]+=StrainSecondDerivativeX.ax;
                HessianMatrix.element[index.x][n+1]+=StrainSecondDerivativeX.az;
                HessianMatrix.element[index.x][n+2]+=StrainSecondDerivativeX.by;
                HessianMatrix.element[index.x][n+3]+=StrainSecondDerivativeX.cz;
                HessianMatrix.element[index.y][n]+=StrainSecondDerivativeY.ax;
                HessianMatrix.element[index.y][n+1]+=StrainSecondDerivativeY.az;
                HessianMatrix.element[index.y][n+2]+=StrainSecondDerivativeY.by;
                HessianMatrix.element[index.y][n+3]+=StrainSecondDerivativeY.cz;
                HessianMatrix.element[index.z][n]+=StrainSecondDerivativeZ.ax;
                HessianMatrix.element[index.z][n+1]+=StrainSecondDerivativeZ.az;
                HessianMatrix.element[index.z][n+2]+=StrainSecondDerivativeZ.by;
                HessianMatrix.element[index.z][n+3]+=StrainSecondDerivativeZ.cz;
                break;
              case MONOCLINIC_GAMMA_ANGLE:
                HessianMatrix.element[index.x][n]+=StrainSecondDerivativeX.ax;
                HessianMatrix.element[index.x][n+1]+=StrainSecondDerivativeX.ay;
                HessianMatrix.element[index.x][n+2]+=StrainSecondDerivativeX.by;
                HessianMatrix.element[index.x][n+3]+=StrainSecondDerivativeX.cz;
                HessianMatrix.element[index.y][n]+=StrainSecondDerivativeY.ax;
                HessianMatrix.element[index.y][n+1]+=StrainSecondDerivativeY.ay;
                HessianMatrix.element[index.y][n+2]+=StrainSecondDerivativeY.by;
                HessianMatrix.element[index.y][n+3]+=StrainSecondDerivativeY.cz;
                HessianMatrix.element[index.z][n]+=StrainSecondDerivativeZ.ax;
                HessianMatrix.element[index.z][n+1]+=StrainSecondDerivativeZ.ay;
                HessianMatrix.element[index.z][n+2]+=StrainSecondDerivativeZ.by;
                HessianMatrix.element[index.z][n+3]+=StrainSecondDerivativeZ.cz;
                break;
            }
            break;
          default:
            fprintf(stderr, "Unknown NPTPRCellType\n");
            exit(0);
            break;
        }
        break;
      case NVT:
      case NVE:
        break;
    }
  }

}

void AddRemainderOfBornTermNumerically(REAL_MATRIX HessianMatrix)
{
  REAL det;
  REAL_MATRIX3x3 StrainDerivativeCentral,StrainDerivativeForward1,StrainDerivativeForward2;
  REAL_MATRIX3x3 StrainDerivativeBackward1,StrainDerivativeBackward2;
  REAL_MATRIX3x3 StoredBox,StoredReplicaBox,strain;
  REAL_MATRIX9x9 StrainSecondDerivative;
  int StoredBoundaryCondition,StoreComputeBornTerm;
  int n,ncell,k1,k2,k3;
  const REAL delta=1e-8;
  //int StoredBornTerm;

  StoredBox=Box[CurrentSystem];
  StoredBoundaryCondition=BoundaryCondition[CurrentSystem];
  BoundaryCondition[CurrentSystem]=TRICLINIC;
  StoreComputeBornTerm=ComputeBornTerm;

  ComputeBornTerm=FALSE;

  InitializeMatrix9x9(&StrainSecondDerivative);
  InitializeMatrix9x9(&BornTerm[CurrentSystem]);

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    CalculateFrameworkBendTorsionForce();
  }
  StrainDerivativeCentral=StrainDerivativeTensor[CurrentSystem];

  SaveFrameworkPositionsToReferenceValues();
  SaveAdsorbateAtomPositionsToReferenceValues();
  SaveCationAtomPositionsToReferenceValues();


  if((Ensemble[CurrentSystem]==NPTPR||Ensemble[CurrentSystem]==NPHPR)&&(NPTPRCellType[CurrentSystem]==REGULAR))
  {
    // ax-element
    strain.ax=1.0+delta; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;       strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0+0.5*delta; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0;           strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;           strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0-0.5*delta; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0;           strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;           strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0-delta; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;       strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.xxxx=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
    StrainSecondDerivative.xxyy=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivative.xxzz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
    StrainSecondDerivative.xxyz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
    StrainSecondDerivative.xxzx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivative.xxxy=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);

    // bx/ay-element
    strain.ax=1.0;       strain.bx=0.5*delta; strain.cx=0.0;
    strain.ay=0.5*delta; strain.by=1.0;       strain.cy=0.0;
    strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;        strain.bx=0.25*delta; strain.cx=0.0;
    strain.ay=0.25*delta; strain.by=1.0;        strain.cy=0.0;
    strain.az=0.0;        strain.bz=0.0;        strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;         strain.bx=-0.25*delta; strain.cx=0.0;
    strain.ay=-0.25*delta; strain.by=1.0;         strain.cy=0.0;
    strain.az=0.0;         strain.bz=0.0;         strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;        strain.bx=-0.5*delta; strain.cx=0.0;
    strain.ay=-0.5*delta; strain.by=1.0;        strain.cy=0.0;
    strain.az=0.0;        strain.bz=0.0;        strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.xyxy=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
    StrainSecondDerivative.yyxy=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivative.yzxy=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
    StrainSecondDerivative.zxxy=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivative.zzxy=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

    // cx/az-element
    strain.ax=1.0;       strain.bx=0.0; strain.cx=0.5*delta;
    strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
    strain.az=0.5*delta; strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;        strain.bx=0.0; strain.cx=0.25*delta;
    strain.ay=0.0;        strain.by=1.0; strain.cy=0.0;
    strain.az=0.25*delta; strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;         strain.bx=0.0; strain.cx=-0.25*delta;
    strain.ay=0.0;         strain.by=1.0; strain.cy=0.0;
    strain.az=-0.25*delta; strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;        strain.bx=0.0; strain.cx=-0.5*delta;
    strain.ay=0.0;        strain.by=1.0; strain.cy=0.0;
    strain.az=-0.5*delta; strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.zxzx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivative.zzzx=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
    StrainSecondDerivative.yzzx=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);

    // by-element
    strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0+delta; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0;       strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;           strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0+0.5*delta; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0;           strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;           strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0-0.5*delta; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0;           strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0-delta; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0;       strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.yyyy=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivative.yyzz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
    StrainSecondDerivative.yyyz=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
    StrainSecondDerivative.yyzx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);

    // cy/bz-element
    strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0;       strain.cy=0.5*delta;
    strain.az=0.0; strain.bz=0.5*delta; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;        strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0;        strain.cy=0.25*delta;
    strain.az=0.0; strain.bz=0.25*delta; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;         strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0;         strain.cy=-0.25*delta;
    strain.az=0.0; strain.bz=-0.25*delta; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;        strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0;        strain.cy=-0.5*delta;
    strain.az=0.0; strain.bz=-0.5*delta; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.yzyz=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
    StrainSecondDerivative.zzyz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

    // cz-element
    strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0; strain.cz=1.0+delta;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0; strain.cz=1.0+0.5*delta;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0; strain.cz=1.0-0.5*delta;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0; strain.cz=1.0-delta;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.zzzz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
  }
  else
  {
    // ax-element
    strain.ax=1.0+delta; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;       strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
      StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
      StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0+0.5*delta; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0;           strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;           strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0-0.5*delta; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0;           strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;           strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0-delta; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;       strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.xxxx=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
    StrainSecondDerivative.xxyy=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivative.xxzz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
    StrainSecondDerivative.xxyz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
    StrainSecondDerivative.xxzx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivative.xxxy=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);

    // bx/ay-element
    strain.ax=1.0;       strain.bx=delta;     strain.cx=0.0;
    strain.ay=0.0;       strain.by=1.0;       strain.cy=0.0;
    strain.az=0.0;       strain.bz=0.0;       strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;        strain.bx=0.5*delta;  strain.cx=0.0;
    strain.ay=0.0;        strain.by=1.0;        strain.cy=0.0;
    strain.az=0.0;        strain.bz=0.0;        strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;         strain.bx=-0.5*delta;  strain.cx=0.0;
    strain.ay=0.0;         strain.by=1.0;         strain.cy=0.0;
    strain.az=0.0;         strain.bz=0.0;         strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;        strain.bx=-delta;     strain.cx=0.0;
    strain.ay=0.0;        strain.by=1.0;        strain.cy=0.0;
    strain.az=0.0;        strain.bz=0.0;        strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.xyxy=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
    StrainSecondDerivative.yyxy=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivative.yzxy=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
    StrainSecondDerivative.zxxy=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivative.zzxy=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

    // cx/az-element
    strain.ax=1.0;       strain.bx=0.0; strain.cx=delta;
    strain.ay=0.0;       strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;       strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;        strain.bx=0.0; strain.cx=0.5*delta;
    strain.ay=0.0;        strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;        strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;         strain.bx=0.0; strain.cx=-0.5*delta;
    strain.ay=0.0;         strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;         strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0;        strain.bx=0.0; strain.cx=-delta;
    strain.ay=0.0;        strain.by=1.0; strain.cy=0.0;
    strain.az=0.0;        strain.bz=0.0; strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.zxzx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
    StrainSecondDerivative.zzzx=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
    StrainSecondDerivative.yzzx=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);

    // by-element
    strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0+delta; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0;       strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;           strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0+0.5*delta; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0;           strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;           strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0-0.5*delta; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0;           strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0-delta; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0;       strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.yyyy=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
    StrainSecondDerivative.yyzz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
    StrainSecondDerivative.yyyz=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
    StrainSecondDerivative.yyzx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);

    // cy/bz-element
    strain.ax=1.0; strain.bx=0.0;       strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0;       strain.cy=delta;
    strain.az=0.0; strain.bz=0.0;       strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;        strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0;        strain.cy=0.5*delta;
    strain.az=0.0; strain.bz=0.0;        strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;         strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0;         strain.cy=-0.5*delta;
    strain.az=0.0; strain.bz=0.0;         strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0;        strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0;        strain.cy=-delta;
    strain.az=0.0; strain.bz=0.0;        strain.cz=1.0;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.yzyz=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
    StrainSecondDerivative.zzyz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

    // cz-element
    strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0; strain.cz=1.0+delta;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0; strain.cz=1.0+0.5*delta;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0; strain.cz=1.0-0.5*delta;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

    strain.ax=1.0; strain.bx=0.0; strain.cx=0.0;
    strain.ay=0.0; strain.by=1.0; strain.cy=0.0;
    strain.az=0.0; strain.bz=0.0; strain.cz=1.0-delta;
    Box[CurrentSystem]=MatrixMatrixMultiplication3x3(strain,StoredBox);
    InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();
    StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
    StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
    StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBendTorsionForce();
    }
    StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

    StrainSecondDerivative.zzzz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);
  }


  // use symmetry Cijkl = Cjikl = Cijlk = Cjilk = Cklij
  StrainSecondDerivative.yxxx=StrainSecondDerivative.xxxy;
  StrainSecondDerivative.xxyx=StrainSecondDerivative.xxxy;
  StrainSecondDerivative.yxyx=StrainSecondDerivative.xyxy;
  StrainSecondDerivative.zxyx=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.yxzx=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.zyxx=StrainSecondDerivative.xxyz;
  StrainSecondDerivative.xyyx=StrainSecondDerivative.xyxy;
  StrainSecondDerivative.yyyx=StrainSecondDerivative.yyxy;
  StrainSecondDerivative.zyyx=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.zyzx=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.xzxx=StrainSecondDerivative.xxzx;
  StrainSecondDerivative.xzyx=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.yzyx=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.zzyx=StrainSecondDerivative.zzxy;
  StrainSecondDerivative.xzzx=StrainSecondDerivative.zxzx;
  StrainSecondDerivative.yxxy=StrainSecondDerivative.xyxy;
  StrainSecondDerivative.yxyy=StrainSecondDerivative.yyxy;
  StrainSecondDerivative.xxzy=StrainSecondDerivative.xxyz;
  StrainSecondDerivative.yxzy=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.zxzy=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.zyxy=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.zyyy=StrainSecondDerivative.yyyz;
  StrainSecondDerivative.xyzy=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.yyzy=StrainSecondDerivative.yyyz;
  StrainSecondDerivative.zyzy=StrainSecondDerivative.yzyz;
  StrainSecondDerivative.xzxy=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.xzyy=StrainSecondDerivative.yyzx;
  StrainSecondDerivative.xzzy=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.yzzy=StrainSecondDerivative.yzyz;
  StrainSecondDerivative.zzzy=StrainSecondDerivative.zzyz;
  StrainSecondDerivative.xxxz=StrainSecondDerivative.xxzx;
  StrainSecondDerivative.yxxz=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.zxxz=StrainSecondDerivative.zxzx;
  StrainSecondDerivative.yxyz=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.yxzz=StrainSecondDerivative.zzxy;
  StrainSecondDerivative.xyxz=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.yyxz=StrainSecondDerivative.yyzx;
  StrainSecondDerivative.zyxz=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.zyyz=StrainSecondDerivative.yzyz;
  StrainSecondDerivative.zyzz=StrainSecondDerivative.zzyz;
  StrainSecondDerivative.xzxz=StrainSecondDerivative.zxzx;
  StrainSecondDerivative.yzxz=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.zzxz=StrainSecondDerivative.zzzx;
  StrainSecondDerivative.xzyz=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.xzzz=StrainSecondDerivative.zzzx;

  // use symmetry Cijkl = Cklij
  StrainSecondDerivative.zxxx=StrainSecondDerivative.xxzx;
  StrainSecondDerivative.xyxx=StrainSecondDerivative.xxxy;
  StrainSecondDerivative.yyxx=StrainSecondDerivative.xxyy;
  StrainSecondDerivative.xyzx=StrainSecondDerivative.zxxy;
  StrainSecondDerivative.yzxx=StrainSecondDerivative.xxyz;
  StrainSecondDerivative.zzxx=StrainSecondDerivative.xxzz;
  StrainSecondDerivative.zxyy=StrainSecondDerivative.yyzx;
  StrainSecondDerivative.xyyy=StrainSecondDerivative.yyxy;
  StrainSecondDerivative.yzyy=StrainSecondDerivative.yyyz;
  StrainSecondDerivative.zzyy=StrainSecondDerivative.yyzz;
  StrainSecondDerivative.zxyz=StrainSecondDerivative.yzzx;
  StrainSecondDerivative.zxzz=StrainSecondDerivative.zzzx;
  StrainSecondDerivative.xyyz=StrainSecondDerivative.yzxy;
  StrainSecondDerivative.xyzz=StrainSecondDerivative.zzxy;
  StrainSecondDerivative.yzzz=StrainSecondDerivative.zzyz;


  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          break;
        case ANISOTROPIC:
          HessianMatrix.element[n][n]+=StrainSecondDerivative.xxxx;
          HessianMatrix.element[n][n+1]+=StrainSecondDerivative.xxyy;
          HessianMatrix.element[n][n+2]+=StrainSecondDerivative.xxzz;
          HessianMatrix.element[n+1][n+1]+=StrainSecondDerivative.yyyy;
          HessianMatrix.element[n+1][n+2]+=StrainSecondDerivative.yyzz;
          HessianMatrix.element[n+2][n+2]+=StrainSecondDerivative.zzzz;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          HessianMatrix.element[n][n]+=StrainSecondDerivative.xxxx;
          HessianMatrix.element[n][n+1]+=StrainSecondDerivative.xxxy;
          HessianMatrix.element[n][n+2]+=StrainSecondDerivative.xxxz;
          HessianMatrix.element[n][n+3]+=StrainSecondDerivative.xxyy;
          HessianMatrix.element[n][n+4]+=StrainSecondDerivative.xxyz;
          HessianMatrix.element[n][n+5]+=StrainSecondDerivative.xxzz;

          HessianMatrix.element[n+1][n+1]+=StrainSecondDerivative.xyxy;
          HessianMatrix.element[n+1][n+2]+=StrainSecondDerivative.xyxz;
          HessianMatrix.element[n+1][n+3]+=StrainSecondDerivative.xyyy;
          HessianMatrix.element[n+1][n+4]+=StrainSecondDerivative.xyyz;
          HessianMatrix.element[n+1][n+5]+=StrainSecondDerivative.xyzz;

          HessianMatrix.element[n+2][n+2]+=StrainSecondDerivative.xzxz;
          HessianMatrix.element[n+2][n+3]+=StrainSecondDerivative.xzyy;
          HessianMatrix.element[n+2][n+4]+=StrainSecondDerivative.xzyz;
          HessianMatrix.element[n+2][n+5]+=StrainSecondDerivative.xzzz;

          HessianMatrix.element[n+3][n+3]+=StrainSecondDerivative.yyyy;
          HessianMatrix.element[n+3][n+4]+=StrainSecondDerivative.yyyz;
          HessianMatrix.element[n+3][n+5]+=StrainSecondDerivative.yyzz;

          HessianMatrix.element[n+4][n+4]+=StrainSecondDerivative.yzyz;
          HessianMatrix.element[n+4][n+5]+=StrainSecondDerivative.yzzz;

          HessianMatrix.element[n+5][n+5]+=StrainSecondDerivative.zzzz;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
               HessianMatrix.element[n][n]+=StrainSecondDerivative.xxxx;
               HessianMatrix.element[n][n+1]+=StrainSecondDerivative.xxyy;
               HessianMatrix.element[n][n+2]+=StrainSecondDerivative.xxyz;
               HessianMatrix.element[n][n+3]+=StrainSecondDerivative.xxzz;

               HessianMatrix.element[n+1][n+1]+=StrainSecondDerivative.yyyy;
               HessianMatrix.element[n+1][n+2]+=StrainSecondDerivative.yyyz;
               HessianMatrix.element[n+1][n+3]+=StrainSecondDerivative.yyzz;

               HessianMatrix.element[n+2][n+2]+=StrainSecondDerivative.yzyz;
               HessianMatrix.element[n+2][n+3]+=StrainSecondDerivative.yzzz;

               HessianMatrix.element[n+3][n+3]+=StrainSecondDerivative.zzzz;
              break;
            case MONOCLINIC_BETA_ANGLE:
               HessianMatrix.element[n][n]+=StrainSecondDerivative.xxxx;
               HessianMatrix.element[n][n+1]+=StrainSecondDerivative.xxxz;
               HessianMatrix.element[n][n+2]+=StrainSecondDerivative.xxyy;
               HessianMatrix.element[n][n+3]+=StrainSecondDerivative.xxzz;

               HessianMatrix.element[n+1][n+1]+=StrainSecondDerivative.xzxz;
               HessianMatrix.element[n+1][n+2]+=StrainSecondDerivative.xzyy;
               HessianMatrix.element[n+1][n+3]+=StrainSecondDerivative.xzzz;

               HessianMatrix.element[n+2][n+2]+=StrainSecondDerivative.yyyy;
               HessianMatrix.element[n+2][n+3]+=StrainSecondDerivative.yyzz;

               HessianMatrix.element[n+3][n+3]+=StrainSecondDerivative.zzzz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
               HessianMatrix.element[n][n]+=StrainSecondDerivative.xxxx;
               HessianMatrix.element[n][n+1]+=StrainSecondDerivative.xxxy;
               HessianMatrix.element[n][n+2]+=StrainSecondDerivative.xxyy;
               HessianMatrix.element[n][n+3]+=StrainSecondDerivative.xxzz;

               HessianMatrix.element[n+1][n+1]+=StrainSecondDerivative.xyxy;
               HessianMatrix.element[n+1][n+2]+=StrainSecondDerivative.xyyy;
               HessianMatrix.element[n+1][n+3]+=StrainSecondDerivative.xyzz;

               HessianMatrix.element[n+2][n+2]+=StrainSecondDerivative.yyyy;
               HessianMatrix.element[n+2][n+3]+=StrainSecondDerivative.yyzz;

               HessianMatrix.element[n+3][n+3]+=StrainSecondDerivative.zzzz;
            default:
              break;
          }
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }


  // restore initial box, positions, and boundary condition
  Box[CurrentSystem]=StoredBox;
  InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);

  PlaceFrameworkInBoxFromReferenceValues();
  PlaceAdsorbateAtomsInBoxFromReferenceValues();
  PlaceCationAtomsInBoxFromReferenceValues();

  BoundaryCondition[CurrentSystem]=StoredBoundaryCondition;
  ComputeBornTerm=StoreComputeBornTerm;
}

void ComputeForcesForMinimalSet(void)
{

  CalculateFrameworkBendForce();
  CalculateFrameworkInversionBendForce();
  CalculateFrameworkTorsionForce();
  CalculateFrameworkImproperTorsionForce();
  CalculateFrameworkBondBondForce();
  CalculateFrameworkBondBendForce();
  CalculateFrameworkBendBendForce();
  CalculateFrameworkBondTorsionForce();
  CalculateFrameworkBendTorsionForce();

  CalculateAdsorbateBendBornTerm();
  CalculateAdsorbateTorsionBornTerm();
  CalculateAdsorbateImproperTorsionBornTerm();

  CalculateCationBendBornTerm();
  CalculateCationTorsionBornTerm();
  CalculateCationImproperTorsionBornTerm();
}

void ComputeCrossTermNumericallyMinimalSet(REAL_MATRIX CrossTerm)
{
  int i,n,index,f1;
  REAL_MATRIX3x3 StrainDerivativeForward1,StrainDerivativeForward2;
  REAL_MATRIX3x3 StrainDerivativeBackward1,StrainDerivativeBackward2;
  REAL_MATRIX3x3 StrainSecondDerivativeX,StrainSecondDerivativeY,StrainSecondDerivativeZ;
  const REAL delta=1e-5;
  VECTOR pos;
  int type;

  n=CrossTerm.m;
  for(i=0;i<n;i++)
  {
    CrossTerm.element[i][0]=CrossTerm.element[i][1]=CrossTerm.element[i][2]=0.0;
    CrossTerm.element[i][3]=CrossTerm.element[i][4]=CrossTerm.element[i][5]=0.0;
    CrossTerm.element[i][6]=CrossTerm.element[i][7]=CrossTerm.element[i][8]=0.0;
  }

  index=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        pos=Framework[CurrentSystem].Atoms[f1][i].Position;
        type=Framework[CurrentSystem].Atoms[f1][i].Type;

        // x-direction
        Framework[CurrentSystem].Atoms[f1][i].Position.x=pos.x+delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

        Framework[CurrentSystem].Atoms[f1][i].Position.x=pos.x+0.5*delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

        Framework[CurrentSystem].Atoms[f1][i].Position.x=pos.x-0.5*delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

        Framework[CurrentSystem].Atoms[f1][i].Position.x=pos.x-delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

        StrainSecondDerivativeX.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
        StrainSecondDerivativeX.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
        StrainSecondDerivativeX.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
        StrainSecondDerivativeX.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
        StrainSecondDerivativeX.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
        StrainSecondDerivativeX.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
        StrainSecondDerivativeX.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
        StrainSecondDerivativeX.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
        StrainSecondDerivativeX.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

        // restore original position
        Framework[CurrentSystem].Atoms[f1][i].Position.x=pos.x;

        // y-direction
        Framework[CurrentSystem].Atoms[f1][i].Position.y=pos.y+delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

        Framework[CurrentSystem].Atoms[f1][i].Position.y=pos.y+0.5*delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

        Framework[CurrentSystem].Atoms[f1][i].Position.y=pos.y-0.5*delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

        Framework[CurrentSystem].Atoms[f1][i].Position.y=pos.y-delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

        StrainSecondDerivativeY.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
        StrainSecondDerivativeY.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
        StrainSecondDerivativeY.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
        StrainSecondDerivativeY.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
        StrainSecondDerivativeY.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
        StrainSecondDerivativeY.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
        StrainSecondDerivativeY.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
        StrainSecondDerivativeY.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
        StrainSecondDerivativeY.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

        // restore original position
        Framework[CurrentSystem].Atoms[f1][i].Position.y=pos.y;

        // z-direction
        Framework[CurrentSystem].Atoms[f1][i].Position.z=pos.z+delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeForward2=StrainDerivativeTensor[CurrentSystem];

        Framework[CurrentSystem].Atoms[f1][i].Position.z=pos.z+0.5*delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeForward1=StrainDerivativeTensor[CurrentSystem];

        Framework[CurrentSystem].Atoms[f1][i].Position.z=pos.z-0.5*delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeBackward1=StrainDerivativeTensor[CurrentSystem];

        Framework[CurrentSystem].Atoms[f1][i].Position.z=pos.z-delta;

        StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
        StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
        StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

        ComputeForcesForMinimalSet();

        StrainDerivativeBackward2=StrainDerivativeTensor[CurrentSystem];

        StrainSecondDerivativeZ.ax=(-StrainDerivativeForward2.ax+8.0*StrainDerivativeForward1.ax-8.0*StrainDerivativeBackward1.ax+StrainDerivativeBackward2.ax)/(6.0*delta);
        StrainSecondDerivativeZ.ay=(-StrainDerivativeForward2.ay+8.0*StrainDerivativeForward1.ay-8.0*StrainDerivativeBackward1.ay+StrainDerivativeBackward2.ay)/(6.0*delta);
        StrainSecondDerivativeZ.az=(-StrainDerivativeForward2.az+8.0*StrainDerivativeForward1.az-8.0*StrainDerivativeBackward1.az+StrainDerivativeBackward2.az)/(6.0*delta);
        StrainSecondDerivativeZ.bx=(-StrainDerivativeForward2.bx+8.0*StrainDerivativeForward1.bx-8.0*StrainDerivativeBackward1.bx+StrainDerivativeBackward2.bx)/(6.0*delta);
        StrainSecondDerivativeZ.by=(-StrainDerivativeForward2.by+8.0*StrainDerivativeForward1.by-8.0*StrainDerivativeBackward1.by+StrainDerivativeBackward2.by)/(6.0*delta);
        StrainSecondDerivativeZ.bz=(-StrainDerivativeForward2.bz+8.0*StrainDerivativeForward1.bz-8.0*StrainDerivativeBackward1.bz+StrainDerivativeBackward2.bz)/(6.0*delta);
        StrainSecondDerivativeZ.cx=(-StrainDerivativeForward2.cx+8.0*StrainDerivativeForward1.cx-8.0*StrainDerivativeBackward1.cx+StrainDerivativeBackward2.cx)/(6.0*delta);
        StrainSecondDerivativeZ.cy=(-StrainDerivativeForward2.cy+8.0*StrainDerivativeForward1.cy-8.0*StrainDerivativeBackward1.cy+StrainDerivativeBackward2.cy)/(6.0*delta);
        StrainSecondDerivativeZ.cz=(-StrainDerivativeForward2.cz+8.0*StrainDerivativeForward1.cz-8.0*StrainDerivativeBackward1.cz+StrainDerivativeBackward2.cz)/(6.0*delta);

        // restore original position
        Framework[CurrentSystem].Atoms[f1][i].Position.z=pos.z;

        CrossTerm.element[3*index][0]=StrainSecondDerivativeX.ax; CrossTerm.element[3*index][3]=StrainSecondDerivativeY.ax; CrossTerm.element[3*index][6]=StrainSecondDerivativeZ.ax;
        CrossTerm.element[3*index][1]=StrainSecondDerivativeX.bx; CrossTerm.element[3*index][4]=StrainSecondDerivativeY.bx; CrossTerm.element[3*index][7]=StrainSecondDerivativeZ.bx;
        CrossTerm.element[3*index][2]=StrainSecondDerivativeX.cx; CrossTerm.element[3*index][5]=StrainSecondDerivativeY.cx; CrossTerm.element[3*index][8]=StrainSecondDerivativeZ.cx;

        CrossTerm.element[3*index+1][0]=StrainSecondDerivativeX.ay; CrossTerm.element[3*index+1][3]=StrainSecondDerivativeY.ay; CrossTerm.element[3*index+1][6]=StrainSecondDerivativeZ.ay;
        CrossTerm.element[3*index+1][1]=StrainSecondDerivativeX.by; CrossTerm.element[3*index+1][4]=StrainSecondDerivativeY.by; CrossTerm.element[3*index+1][7]=StrainSecondDerivativeZ.by;
        CrossTerm.element[3*index+1][2]=StrainSecondDerivativeX.cy; CrossTerm.element[3*index+1][5]=StrainSecondDerivativeY.cy; CrossTerm.element[3*index+1][8]=StrainSecondDerivativeZ.cy;

        CrossTerm.element[3*index+2][0]=StrainSecondDerivativeX.az; CrossTerm.element[3*index+2][3]=StrainSecondDerivativeY.az; CrossTerm.element[3*index+2][6]=StrainSecondDerivativeZ.az;
        CrossTerm.element[3*index+2][1]=StrainSecondDerivativeX.bz; CrossTerm.element[3*index+2][4]=StrainSecondDerivativeY.bz; CrossTerm.element[3*index+2][7]=StrainSecondDerivativeZ.bz;
        CrossTerm.element[3*index+2][2]=StrainSecondDerivativeX.cz; CrossTerm.element[3*index+2][5]=StrainSecondDerivativeY.cz; CrossTerm.element[3*index+2][8]=StrainSecondDerivativeZ.cz;

        index++;
      }
    }
  }
}

