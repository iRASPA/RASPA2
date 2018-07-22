/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'equations_of_state.c' is part of RASPA-2.0

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
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include "simulation.h"
#include "ewald.h"
#include "utils.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "internal_energy.h"
#include "inter_energy.h"
#include "integration.h"
#include "potentials.h"
#include "mc_moves.h"
#include "input.h"
#include "output.h"
#include "cbmc.h"
#include "ewald.h"
#include "grids.h"
#include "statistics.h"
#include "thermo_baro_stats.h"
#include "equations_of_state.h"

// The fugacity coefficity phi is defined through the relationship RT ln phi = G - G^{ideal gas}.
// The deviation of phi from unity is a measure of the deviation of the fluid from ideal gas behaviour.
// Note that in the special case that phi^liq=phi^vap, the liquid and vapour phases will be at equilibrium
// with each other. This enables saturated vapour pressures to be predicted: specify the T of interest and
// then vary P until phi^liq=phi^vap. Alternatively boiling points may be predicted by specifying the P of
// interest and varying T until the fugacity coefficient condition is obeyed.

int EquationOfState;
int MultiComponentMixingRules;
REAL **BinaryInteractionParameter;
int **ComputeFugacityCoefficient;

REAL *HeliumVoidFraction;
REAL *ExcessVolume;

int *FluidState;
REAL Compressibility[3];

static REAL *a;
static REAL *A;
static REAL *b;
static REAL *B;
static REAL **aij;
static REAL **Aij;
static REAL (*FugacityCoefficients)[3];

void RescaleMolarFractions(void)
{
  int i,j;
  REAL total;

  for(j=0;j<NumberOfSystems;j++)
  {
    total=0.0;
    for(i=0;i<NumberOfComponents;i++)
      if(Components[i].Swapable)
        total+=Components[i].MolFraction[j];

    if(total>0.0)
    {
      for(i=0;i<NumberOfComponents;i++)
        if(Components[i].Swapable)
          Components[i].MolFraction[j]/=total;
    }
    else
    {
      for(i=0;i<NumberOfComponents;i++)
        if(Components[i].Swapable)
          Components[i].MolFraction[j]=1.0/NumberOfComponents;
    }
  }
}

void ComputeGasPropertiesForAllSystems(void)
{
  int i,j,k;
  REAL alpha,w,kappa,Tc,Pc,Tr;
  REAL Amix,Bmix;
  REAL excess_volume,Pressure,Temperature,FrameworkMass;
  int NumberOfSolutions;
  REAL Coefficients[4];
  REAL number_of_unit_cells,temp;

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    FrameworkMass=0.0;
    for(i=0;i<Framework[CurrentSystem].NumberOfFrameworks;i++)
      FrameworkMass+=Framework[CurrentSystem].FrameworkMassPerComponent[i];

    // add the mass of cations to the framework
    for(i=0;i<NumberOfComponents;i++)
      if((!Components[i].Swapable)&&(Components[i].ExtraFrameworkMolecule))
        FrameworkMass+=Components[i].NumberOfMolecules[CurrentSystem]*
                       Components[i].Mass;

    Framework[CurrentSystem].FrameworkMass=FrameworkMass;
    Framework[CurrentSystem].FrameworkDensity=FrameworkMass/(Volume[CurrentSystem]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT);

    number_of_unit_cells=NumberOfUnitCells[CurrentSystem].x*NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z;
    for(i=0;i<NumberOfComponents;i++)
    {
      // molec/uc -> mol/kg
      Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]=1000.0*number_of_unit_cells/FrameworkMass;
      // molec/uc -> g/g
      Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]=Components[i].Mass*1000.0*number_of_unit_cells/FrameworkMass;
    }

    // use global pressure to compute partial pressures
    if(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]>=0.0)
    {
      for(i=0;i<NumberOfComponents;i++)
        Components[i].PartialPressure[CurrentSystem]=therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*
                   Components[i].MolFraction[CurrentSystem];
    }
    else // check that all partial pressure are properly defined
    {
      for(i=0;i<NumberOfComponents;i++)
      {
        if(Components[i].Swapable)
        {
          if(Components[i].PartialPressure[CurrentSystem]<0.0)
          {
            fprintf(stderr, "Error: Partial pressure of component %d [%s] is NOT set !!\n",
              i,Components[i].Name);
            exit(-1);
          }
        }
      }
    }

    // check: critical T and/or P not set, and also fugacity coefficient is not initialized
    for(i=0;i<NumberOfComponents;i++)
    {
      if(Components[i].Swapable)
      {
        Tc=Components[i].CriticalTemperature;
        Pc=Components[i].CriticalPressure;
        if((Tc<1e-10)||(Pc<1e-10))
        {
          if(Components[i].FugacityCoefficient[CurrentSystem]<0.0)
          {
            fprintf(stderr, "Error: Fugacity coefficient of component %d [%s] is NOT set for system %d,\n"
                   "       Critial pressure/temperature NOT set either (set to zero in molecule-definition)\n",
              i,Components[i].Name,CurrentSystem);
            exit(-1);
          }
        }
      }
    }

    // here, use pressure in Pascal and temperature in Kelvin
    Pressure=therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR;
    Temperature=therm_baro_stats.ExternalTemperature[CurrentSystem];

    for(i=0;i<NumberOfComponents;i++)
    {
      if(Components[i].Swapable)
      {
        Tc=Components[i].CriticalTemperature;
        Pc=Components[i].CriticalPressure;
        w=Components[i].AcentricFactor;
        Tr=Temperature/Tc;
        switch(EquationOfState)
        {
          case SOAVE_REDLICH_KWONG:
            kappa=0.480+1.574*w-0.176*SQR(w);
            alpha=SQR(1.0+kappa*(1.0-sqrt(Tr)));
            a[i]=0.42748*alpha*SQR(MOLAR_GAS_CONSTANT*Tc)/Pc;
            b[i]=0.08664*MOLAR_GAS_CONSTANT*Tc/Pc;
            A[i]=a[i]*Pressure/SQR(MOLAR_GAS_CONSTANT*Temperature);
            B[i]=b[i]*Pressure/(MOLAR_GAS_CONSTANT*Temperature);
            break;
          case PENG_ROBINSON_GASEM:
           kappa=0.134+0.508*w-0.0467*SQR(w);
           alpha=exp((2.0+0.836*Tr)*(1.0-pow(Tr,kappa)));
           a[i]=0.45724*alpha*SQR(MOLAR_GAS_CONSTANT*Tc)/Pc;
           b[i]=0.07780*MOLAR_GAS_CONSTANT*Tc/Pc;
           A[i]=a[i]*Pressure/SQR(MOLAR_GAS_CONSTANT*Temperature);
           B[i]=b[i]*Pressure/(MOLAR_GAS_CONSTANT*Temperature);
           break;
          case PENG_ROBINSON:
          default:
            kappa=0.37464+1.54226*w-0.26992*SQR(w);
            alpha=SQR(1.0+kappa*(1.0-sqrt(Tr)));
            a[i]=0.45724*alpha*SQR(MOLAR_GAS_CONSTANT*Tc)/Pc;
            b[i]=0.07780*MOLAR_GAS_CONSTANT*Tc/Pc;
            A[i]=a[i]*Pressure/SQR(MOLAR_GAS_CONSTANT*Temperature);
            B[i]=b[i]*Pressure/(MOLAR_GAS_CONSTANT*Temperature);
            break;
        }
      }
    }

    switch(MultiComponentMixingRules)
    {
      default:
      case VAN_DER_WAALS:
        for(i=0;i<NumberOfComponents;i++)
          for(j=0;j<NumberOfComponents;j++)
          {
            if(Components[i].Swapable&&Components[j].Swapable)
            {
              aij[i][j]=(1.0-BinaryInteractionParameter[i][j])*sqrt(a[i]*a[j]);
              Aij[i][j]=(1.0-BinaryInteractionParameter[i][j])*sqrt(A[i]*A[j]);
            }
          }
        Amix=0.0;
        Bmix=0.0;
        for(i=0;i<NumberOfComponents;i++)
        {
          if(Components[i].Swapable)
          {
            Bmix+=Components[i].MolFraction[CurrentSystem]*b[i];
            for(j=0;j<NumberOfComponents;j++)
              if(Components[j].Swapable)
                Amix+=Components[i].MolFraction[CurrentSystem]*Components[j].MolFraction[CurrentSystem]*aij[i][j];
          }
        }
        Amix*=Pressure/SQR(MOLAR_GAS_CONSTANT*Temperature);
        Bmix*=Pressure/(MOLAR_GAS_CONSTANT*Temperature);
        break;
    }

    switch(EquationOfState)
    {
      case SOAVE_REDLICH_KWONG:
        Coefficients[3]=1.0;
        Coefficients[2]=-1.0;
        Coefficients[1]=Amix-Bmix-SQR(Bmix);
        Coefficients[0]=-Amix*Bmix;
        break;
      case PENG_ROBINSON:
      case PENG_ROBINSON_GASEM:
      default:
        Coefficients[3]=1.0;
        Coefficients[2]=Bmix-1.0;
        Coefficients[1]=Amix-3.0*SQR(Bmix)-2.0*Bmix;
        Coefficients[0]=-(Amix*Bmix-SQR(Bmix)-CUBE(Bmix));
        break;
    }

    // solve equation-of-state
    cubic(Coefficients,Compressibility,&NumberOfSolutions);

    // sort the compressibilities, Compressibility[0] the highest, Compressibility[2] the lowest
    if (Compressibility[0]<Compressibility[1])
    {
      temp=Compressibility[0];
      Compressibility[0]=Compressibility[1];
      Compressibility[1]=temp;
    }
    if (Compressibility[1]<Compressibility[2])
    {
      temp=Compressibility[1];
      Compressibility[1]=Compressibility[2];
      Compressibility[2]=temp;
    }
    if(Compressibility[0]<Compressibility[1])
    {
      temp=Compressibility[0];
      Compressibility[0]=Compressibility[1];
      Compressibility[1]=temp;
    }

    // In cases where three roots are obtained, the highest value corresponds to the solution for vapour phase,
    // and the lowest value corresponds to the solution for liquid phase. (The intermediate solution has no
    // physical meaning). Note that the solutions in this case don't imply that the liquid and vapour phases
    // are in equilibrium with each other. Rather, one of the phases will be stable and the other will be a
    // metastable state. Only at the special condition that T = Tsat and P = Psat are vapour and liquid in
    // equilibrium with each other.

    // In other cases, only a single root is obtained. This may correspond to liquid, vapour or supercritical fluid,
    // depending on the relative magnitudes of P and Pc, and T and Tc.

    switch(EquationOfState)
    {
      case SOAVE_REDLICH_KWONG:
        for(i=0;i<NumberOfComponents;i++)
        {
          if(Components[i].Swapable)
          {
            for(j=0;j<NumberOfSolutions;j++)
            {
              temp=0.0;
              for(k=0;k<NumberOfComponents;k++)
                temp+=2.0*Components[k].MolFraction[CurrentSystem]*Aij[i][k];

              FugacityCoefficients[i][j]=exp((B[i]/Bmix)*(Compressibility[j]-1.0)-log(Compressibility[j]-Bmix)-
                  Amix*(temp/Amix-B[i]/Bmix)*log(1.0+Bmix/Compressibility[j])/Bmix);
            }
          }
        }
        break;
      case PENG_ROBINSON_GASEM:
      case PENG_ROBINSON:
      default:
        for(i=0;i<NumberOfComponents;i++)
        {
          if(Components[i].Swapable)
          {
            for(j=0;j<NumberOfSolutions;j++)
            {
              temp=0.0;
              for(k=0;k<NumberOfComponents;k++)
                temp+=2.0*Components[k].MolFraction[CurrentSystem]*Aij[i][k];

              FugacityCoefficients[i][j]=exp((B[i]/Bmix)*(Compressibility[j]-1.0)-log(Compressibility[j]-Bmix)
                        -(Amix/(2.0*sqrt(2.0)*Bmix))*(temp/Amix-B[i]/Bmix)*
                  log((Compressibility[j]+(1.0+sqrt(2.0))*Bmix)/(Compressibility[j]+(1.0-sqrt(2))*Bmix)));
            }
          }
        }
        break;
    }

    if(ExcessVolume[CurrentSystem]>0.0)
      excess_volume=ExcessVolume[CurrentSystem]*HeliumVoidFraction[CurrentSystem];
    else
      excess_volume=Volume[CurrentSystem]*HeliumVoidFraction[CurrentSystem];

    // get the gas-phase fugacity coefficient
    for(i=0;i<NumberOfComponents;i++)
    {

      if(Components[i].Swapable)
      {
        if(NumberOfSolutions==1)
        {
          if(ComputeFugacityCoefficient[CurrentSystem][i])
            Components[i].FugacityCoefficient[CurrentSystem]=FugacityCoefficients[i][0];
          Components[i].AmountOfExcessMolecules[CurrentSystem]=
             Components[i].MolFraction[CurrentSystem]*AVOGADRO_CONSTANT*
             excess_volume*CUBE(ANGSTROM)*Pressure/(Compressibility[0]*MOLAR_GAS_CONSTANT*Temperature);
          Components[i].BulkFluidDensity[CurrentSystem]=(Pressure/(Compressibility[0]*MOLAR_GAS_CONSTANT*Temperature))*
                                                         Components[i].Mass/1000.0;
          Components[i].Compressibility[CurrentSystem]=Compressibility[0];
          if((Temperature>Components[i].CriticalTemperature)&&(Pressure>Components[i].CriticalPressure))
            FluidState[CurrentSystem]=SUPER_CRITICAL_FLUID;
          else if((Temperature<Components[i].CriticalTemperature)&&(Pressure<Components[i].CriticalPressure))
            FluidState[CurrentSystem]=VAPOUR;
          else if((Temperature<Components[i].CriticalTemperature)&&(Pressure>Components[i].CriticalPressure))
            FluidState[CurrentSystem]=LIQUID;
        }
        else
        {
          if(Compressibility[2]>0.0)
          {
            if(FugacityCoefficients[i][0]<FugacityCoefficients[i][2])
            {
              // vapour (stable) and liquid (metastable)
              FluidState[CurrentSystem]=VAPOUR_STABLE;
              if(ComputeFugacityCoefficient[CurrentSystem][i])
                Components[i].FugacityCoefficient[CurrentSystem]=FugacityCoefficients[i][0];
              Components[i].AmountOfExcessMolecules[CurrentSystem]=Components[i].MolFraction[CurrentSystem]*AVOGADRO_CONSTANT*
                excess_volume*CUBE(ANGSTROM)*Pressure/(Compressibility[0]*MOLAR_GAS_CONSTANT*Temperature);
              Components[i].BulkFluidDensity[CurrentSystem]=(Pressure/(Compressibility[0]*MOLAR_GAS_CONSTANT*Temperature))*
                                                            Components[i].Mass/1000.0;
              Components[i].Compressibility[CurrentSystem]=Compressibility[0];
            }
            else if(FugacityCoefficients[i][0]>FugacityCoefficients[i][2])
            {
              // vapour (metastable) and liquid (stable)
              FluidState[CurrentSystem]=LIQUID_STABLE;
              if(ComputeFugacityCoefficient[CurrentSystem][i])
                Components[i].FugacityCoefficient[CurrentSystem]=FugacityCoefficients[i][2];
              Components[i].AmountOfExcessMolecules[CurrentSystem]=Components[i].MolFraction[CurrentSystem]*AVOGADRO_CONSTANT*
                excess_volume*CUBE(ANGSTROM)*Pressure/(Compressibility[2]*MOLAR_GAS_CONSTANT*Temperature);
              Components[i].BulkFluidDensity[CurrentSystem]=(Pressure/(Compressibility[2]*MOLAR_GAS_CONSTANT*Temperature))*
                                                            Components[i].Mass/1000.0;
              Components[i].Compressibility[CurrentSystem]=Compressibility[2];
            }
            else
            {
              // vapour (stable) and liquid (stable)
              FluidState[CurrentSystem]=VAPOUR_LIQUID_STABLE;
              if(ComputeFugacityCoefficient[CurrentSystem][i])
                Components[i].FugacityCoefficient[CurrentSystem]=FugacityCoefficients[i][0];
              Components[i].AmountOfExcessMolecules[CurrentSystem]=Components[i].MolFraction[CurrentSystem]*AVOGADRO_CONSTANT*
                excess_volume*CUBE(ANGSTROM)*Pressure/(Compressibility[0]*MOLAR_GAS_CONSTANT*Temperature);
              Components[i].BulkFluidDensity[CurrentSystem]=(Pressure/(Compressibility[0]*MOLAR_GAS_CONSTANT*Temperature))*
                                                            Components[i].Mass/1000.0;
              Components[i].Compressibility[CurrentSystem]=Compressibility[0];
            }
          }
          else
          {
            if(ComputeFugacityCoefficient[CurrentSystem][i])
              Components[i].FugacityCoefficient[CurrentSystem]=FugacityCoefficients[i][0];
            Components[i].AmountOfExcessMolecules[CurrentSystem]=Components[i].MolFraction[CurrentSystem]*AVOGADRO_CONSTANT*
              excess_volume*CUBE(ANGSTROM)*Pressure/(Compressibility[0]*MOLAR_GAS_CONSTANT*Temperature);
            Components[i].BulkFluidDensity[CurrentSystem]=(Pressure/(Compressibility[0]*MOLAR_GAS_CONSTANT*Temperature))*
                                                          Components[i].Mass/1000.0;
            Components[i].Compressibility[CurrentSystem]=Compressibility[0];
            if((Temperature>Components[i].CriticalTemperature)&&(Pressure>Components[i].CriticalPressure))
              FluidState[CurrentSystem]=SUPER_CRITICAL_FLUID;
            else if((Temperature<Components[i].CriticalTemperature)&&(Pressure<Components[i].CriticalPressure))
              FluidState[CurrentSystem]=VAPOUR;
            else if((Temperature<Components[i].CriticalTemperature)&&(Pressure>Components[i].CriticalPressure))
              FluidState[CurrentSystem]=LIQUID;
          }
        }
      }
      else
      {
        Components[i].FugacityCoefficient[CurrentSystem]=-1.0;
        Components[i].PartialPressure[CurrentSystem]=-1.0;
      }
    }

    // The term cm^3 (STP) is not a unit of volume but a unit to express the number of gas molecules
    // It is essentially a term in units of pV, e.g. cm^3 at 1 atm and 273 K.
    // The volume of 1 mole of a gas is 22.4 liters

    for(i=0;i<NumberOfComponents;i++)
    {
      // molec/uc -> cm^3 (STP)/g
      Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]=1e6*(MOLAR_GAS_CONSTANT*273.15/ATM_TO_PA)*number_of_unit_cells/FrameworkMass;

      // molec/uc -> cm^3 (STP)/cm^3
      Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]=number_of_unit_cells*(MOLAR_GAS_CONSTANT*273.15/ATM_TO_PA)
                                                             /(Volume[CurrentSystem]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT);

      // mol/kg -> cm^3 (STP)/g
      Components[i].MOL_PER_KG_TO_CC_STP_G[CurrentSystem]=1e3*(MOLAR_GAS_CONSTANT*273.15/ATM_TO_PA);

      // mol/kg -> cm^3 (STP)/cm^3
      Components[i].MOL_PER_KG_TO_CC_STP_CC[CurrentSystem]=(MOLAR_GAS_CONSTANT*273.15/ATM_TO_PA)*FrameworkMass
                                                            /(1e3*Volume[CurrentSystem]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT);
    }

    for(i=0;i<NumberOfComponents;i++)
    {
      if(Components[i].Swapable)
      {
        if(Components[i].PartialPressure[CurrentSystem]<0.0)
        {
          fprintf(stderr, "Error: Partial pressure of component %d [%s] is NOT set !!\n",
            i,Components[i].Name);
          exit(-1);
        }
      }
    }
  }
}

static int versionNumber=1;

void WriteRestartEquationOfState(FILE *FilePtr)
{
  int i;
  REAL Check;

  fwrite(&versionNumber,sizeof(int),1,FilePtr);

  for(i=0;i<NumberOfComponents;i++)
    fwrite(BinaryInteractionParameter[i],NumberOfComponents,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
    fwrite(ComputeFugacityCoefficient[i],NumberOfComponents,sizeof(int),FilePtr);

  fwrite(HeliumVoidFraction,NumberOfSystems,sizeof(REAL),FilePtr);
  fwrite(ExcessVolume,NumberOfSystems,sizeof(REAL),FilePtr);
  fwrite(FluidState,NumberOfSystems,sizeof(int),FilePtr);

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}


void AllocateEquationOfStateMemory(void)
{
  int i;

  BinaryInteractionParameter=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
  a=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
  A=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
  b=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
  B=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
  aij=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
  Aij=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
  for(i=0;i<NumberOfComponents;i++)
  {
    BinaryInteractionParameter[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    aij[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    Aij[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
  }
  FugacityCoefficients=(REAL(*)[3])calloc(NumberOfComponents,sizeof(REAL[3]));

  ComputeFugacityCoefficient=(int**)calloc(NumberOfSystems,sizeof(int*));
  for(i=0;i<NumberOfSystems;i++)
    ComputeFugacityCoefficient[i]=(int*)calloc(NumberOfComponents,sizeof(int));

  HeliumVoidFraction=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  ExcessVolume=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  FluidState=(int*)calloc(NumberOfSystems,sizeof(int));
}

void ReadRestartEquationOfState(FILE *FilePtr)
{
  int i;
  REAL Check;
  int readversionNumber=0;

  AllocateEquationOfStateMemory();

  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumber)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  for(i=0;i<NumberOfComponents;i++)
    fread(BinaryInteractionParameter[i],NumberOfComponents,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
    fread(ComputeFugacityCoefficient[i],NumberOfComponents,sizeof(int),FilePtr);

  fread(HeliumVoidFraction,NumberOfSystems,sizeof(REAL),FilePtr);
  fread(ExcessVolume,NumberOfSystems,sizeof(REAL),FilePtr);
  fread(FluidState,NumberOfSystems,sizeof(int),FilePtr);

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartEquationOfState)\n");
    ContinueAfterCrash=FALSE;
  }
}

