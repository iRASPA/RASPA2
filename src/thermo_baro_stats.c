/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'thermo_baro_stats.c' is part of RASPA-2.0

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
#include "simulation.h"
#include "output.h"
#include "molecule.h"
#include "framework.h"
#include "inter_force.h"
#include "framework_energy.h"
#include "grids.h"
#include "thermo_baro_stats.h"
#include "complex.h"
#include "matrix.h"
#include "ewald.h"

/*************************************************************************************
  Thermo-and barostats

  algorithms from:

  M. E. Tuckerman, J. Alejandre, R. Lopez-Rendon, A. L. Jochem, and G. J. Martyna
  A Liouville-operator derived measure-preserving integrator for molecular dynamics
  simulations in the isothermal-isobaric ensemble
  Journal of Physics A: Mathematical and General, volume 39, pages 5629-5651, 2006

  The NVE, NVT, and NPT integrators are
  1) preserving the correct phase-space volume element, i.e. measure-preserving
  2) reversible

  The algorithms are modified for rigid molecules.

  With the Andersen method (1980) of pressure control, the volume of the cell can change,
  but its shape is preserved by allowing the cell to change isotropically.
  The Anderson method is useful for liquid simulations since the box could become quite
  elongated in the absence of restoring forces if the shape of the cell were allowed to
  change. A constant shape also makes the dynamics analysis easier.
  However, this method is not very useful for studying materials under nonisotropic
  stress or phase transitions, which involve changes in both cell lengths and cell angles
  (for these conditions, the Parrinello-Rahman method should be used).

  The Parrinello-Rahman method of pressure and stress control can allow simulation of a
  model under externally applied stress. This is useful for studying the stress-strain
  relationship of materials. Both the shape and the volume of the cell can change, so
  that the internal stress of the system can match the externally applied stress.

  Raspa performs atomic scaling. The major advantage of using atomic scaling of the
  coordinates is that atom overlapping can be avoided. Such overlap can occur if centers
  of mass are moved instead of individual atoms. In addition, for large models having
  internal flexibility, atomic scaling yields a smoother response to pressure changes.
  For rigid molecules the molecular scaling and molecular virial is computed.
 *************************************************************************************/

int NumberOfIsothermPressures;
int CurrentIsothermPressure;

// data structure containing the parameters for the thermo-and barostats
THERMO_BAROSTATS therm_baro_stats;

static REAL w[10];

static REAL_MATRIX3x3 CellForcePressurePart;
static REAL_MATRIX3x3 CellForceKineticPressurePart;
static REAL_MATRIX3x3 CellForce;   // NEW
static REAL_MATRIX3x3 *CellVelocity;
static REAL *CellMass;

static REAL *LnVolumePosition;
REAL *LnVolumeVelocity;
REAL *LnVolumeMass;

static REAL **ThermostatDegreesOfFreedom;
static REAL **ThermostatForce;
static REAL **ThermostatVelocity;
static REAL **ThermostatPosition;
static REAL **ThermostatMass;

static REAL **BarostatDegreesOfFreedom;
static REAL **BarostatForce;
static REAL **BarostatVelocity;
static REAL **BarostatPosition;
static REAL **BarostatMass;

// the translational degrees of freedom are splitted into
// a) framework
// b) adsorbates
// c) cations
// this results in better (faster) equipartition for flexible frameworks, the high frequency oscillations
// of the framework and the low frequency oscillations of the adsorbates leads to adiabatic-type behaviour

static REAL **ThermostatDegreesOfFreedomTranslation;
static REAL **ThermostatForceTranslation;
static REAL **ThermostatVelocityTranslation;
static REAL **ThermostatPositionTranslation;
static REAL **ThermostatMassTranslation;

static REAL **ThermostatDegreesOfFreedomTranslationFramework;
static REAL **ThermostatForceTranslationFramework;
static REAL **ThermostatVelocityTranslationFramework;
static REAL **ThermostatPositionTranslationFramework;
static REAL **ThermostatMassTranslationFramework;

static REAL **ThermostatDegreesOfFreedomTranslationAdsorbates;
static REAL **ThermostatForceTranslationAdsorbates;
static REAL **ThermostatVelocityTranslationAdsorbates;
static REAL **ThermostatPositionTranslationAdsorbates;
static REAL **ThermostatMassTranslationAdsorbates;

static REAL **ThermostatDegreesOfFreedomTranslationCations;
static REAL **ThermostatForceTranslationCations;
static REAL **ThermostatVelocityTranslationCations;
static REAL **ThermostatPositionTranslationCations;
static REAL **ThermostatMassTranslationCations;

static REAL **ThermostatDegreesOfFreedomRotation;
static REAL **ThermostatForceRotation;
static REAL **ThermostatVelocityRotation;
static REAL **ThermostatPositionRotation;
static REAL **ThermostatMassRotation;

REAL_MATRIX3x3 GetKineticStressTensor(void)
{
  int i,j,k,l,f1;
  int Type,A;
  REAL Mass;
  REAL_MATRIX3x3 KineticPressureTensor;
  VECTOR Velocity;
  INT_VECTOR3 Fixed;

  KineticPressureTensor.ax=KineticPressureTensor.bx=KineticPressureTensor.cx=0.0;
  KineticPressureTensor.ay=KineticPressureTensor.by=KineticPressureTensor.cy=0.0;
  KineticPressureTensor.az=KineticPressureTensor.bz=KineticPressureTensor.cz=0.0;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][k].Type].Mass;
          Fixed=Framework[CurrentSystem].Atoms[f1][k].Fixed;

          if(!Fixed.x)               KineticPressureTensor.ax+=Mass*Framework[CurrentSystem].Atoms[f1][k].Velocity.x*Framework[CurrentSystem].Atoms[f1][k].Velocity.x;
          if((!Fixed.x)&&(!Fixed.y)) KineticPressureTensor.ay+=Mass*Framework[CurrentSystem].Atoms[f1][k].Velocity.x*Framework[CurrentSystem].Atoms[f1][k].Velocity.y;
          if((!Fixed.x)&&(!Fixed.z)) KineticPressureTensor.az+=Mass*Framework[CurrentSystem].Atoms[f1][k].Velocity.x*Framework[CurrentSystem].Atoms[f1][k].Velocity.z;
          if((!Fixed.y)&&(!Fixed.x)) KineticPressureTensor.bx+=Mass*Framework[CurrentSystem].Atoms[f1][k].Velocity.y*Framework[CurrentSystem].Atoms[f1][k].Velocity.x;
          if(!Fixed.y)               KineticPressureTensor.by+=Mass*Framework[CurrentSystem].Atoms[f1][k].Velocity.y*Framework[CurrentSystem].Atoms[f1][k].Velocity.y;
          if((!Fixed.y)&&(!Fixed.z)) KineticPressureTensor.bz+=Mass*Framework[CurrentSystem].Atoms[f1][k].Velocity.y*Framework[CurrentSystem].Atoms[f1][k].Velocity.z;
          if((!Fixed.z)&&(!Fixed.x)) KineticPressureTensor.cx+=Mass*Framework[CurrentSystem].Atoms[f1][k].Velocity.z*Framework[CurrentSystem].Atoms[f1][k].Velocity.x;
          if((!Fixed.z)&&(!Fixed.y)) KineticPressureTensor.cy+=Mass*Framework[CurrentSystem].Atoms[f1][k].Velocity.z*Framework[CurrentSystem].Atoms[f1][k].Velocity.y;
          if(!Fixed.z)               KineticPressureTensor.cz+=Mass*Framework[CurrentSystem].Atoms[f1][k].Velocity.z*Framework[CurrentSystem].Atoms[f1][k].Velocity.z;
        }
      }
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Velocity=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
        Mass=Components[Type].Groups[l].Mass;

        Fixed=Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x)               KineticPressureTensor.ax+=Mass*Velocity.x*Velocity.x;
        if((!Fixed.x)&&(!Fixed.y)) KineticPressureTensor.ay+=Mass*Velocity.x*Velocity.y;
        if((!Fixed.x)&&(!Fixed.z)) KineticPressureTensor.az+=Mass*Velocity.x*Velocity.z;
        if((!Fixed.y)&&(!Fixed.x)) KineticPressureTensor.bx+=Mass*Velocity.y*Velocity.x;
        if(!Fixed.y)               KineticPressureTensor.by+=Mass*Velocity.y*Velocity.y;
        if((!Fixed.y)&&(!Fixed.z)) KineticPressureTensor.bz+=Mass*Velocity.y*Velocity.z;
        if((!Fixed.z)&&(!Fixed.x)) KineticPressureTensor.cx+=Mass*Velocity.z*Velocity.x;
        if((!Fixed.z)&&(!Fixed.y)) KineticPressureTensor.cy+=Mass*Velocity.z*Velocity.y;
        if(!Fixed.z)               KineticPressureTensor.cz+=Mass*Velocity.z*Velocity.z;
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];

          Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
          Fixed=Adsorbates[CurrentSystem][i].Atoms[A].Fixed;

          if(!Fixed.x)               KineticPressureTensor.ax+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x;
          if((!Fixed.x)&&(!Fixed.y)) KineticPressureTensor.ay+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y;
          if((!Fixed.x)&&(!Fixed.z)) KineticPressureTensor.az+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z;
          if((!Fixed.y)&&(!Fixed.x)) KineticPressureTensor.bx+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x;
          if(!Fixed.y)               KineticPressureTensor.by+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y;
          if((!Fixed.y)&&(!Fixed.z)) KineticPressureTensor.bz+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z;
          if((!Fixed.z)&&(!Fixed.x)) KineticPressureTensor.cx+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x;
          if((!Fixed.z)&&(!Fixed.y)) KineticPressureTensor.cy+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y;
          if(!Fixed.z)               KineticPressureTensor.cz+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z;
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
        Velocity=Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
        Mass=Components[Type].Groups[l].Mass;

        Fixed=Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x)               KineticPressureTensor.ax+=Mass*Velocity.x*Velocity.x;
        if((!Fixed.x)&&(!Fixed.y)) KineticPressureTensor.ay+=Mass*Velocity.x*Velocity.y;
        if((!Fixed.x)&&(!Fixed.z)) KineticPressureTensor.az+=Mass*Velocity.x*Velocity.z;
        if((!Fixed.y)&&(!Fixed.x)) KineticPressureTensor.bx+=Mass*Velocity.y*Velocity.x;
        if(!Fixed.y)               KineticPressureTensor.by+=Mass*Velocity.y*Velocity.y;
        if((!Fixed.y)&&(!Fixed.z)) KineticPressureTensor.bz+=Mass*Velocity.y*Velocity.z;
        if((!Fixed.z)&&(!Fixed.x)) KineticPressureTensor.cx+=Mass*Velocity.z*Velocity.x;
        if((!Fixed.z)&&(!Fixed.y)) KineticPressureTensor.cy+=Mass*Velocity.z*Velocity.y;
        if(!Fixed.z)               KineticPressureTensor.cz+=Mass*Velocity.z*Velocity.z;
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];

          Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[A].Type].Mass;
          Fixed=Cations[CurrentSystem][i].Atoms[A].Fixed;

          if(!Fixed.x)               KineticPressureTensor.ax+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.x*Cations[CurrentSystem][i].Atoms[A].Velocity.x;
          if((!Fixed.x)&&(!Fixed.y)) KineticPressureTensor.ay+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.x*Cations[CurrentSystem][i].Atoms[A].Velocity.y;
          if((!Fixed.x)&&(!Fixed.z)) KineticPressureTensor.az+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.x*Cations[CurrentSystem][i].Atoms[A].Velocity.z;
          if((!Fixed.y)&&(!Fixed.x)) KineticPressureTensor.bx+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.y*Cations[CurrentSystem][i].Atoms[A].Velocity.x;
          if(!Fixed.y)               KineticPressureTensor.by+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.y*Cations[CurrentSystem][i].Atoms[A].Velocity.y;
          if((!Fixed.y)&&(!Fixed.z)) KineticPressureTensor.bz+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.y*Cations[CurrentSystem][i].Atoms[A].Velocity.z;
          if((!Fixed.z)&&(!Fixed.x)) KineticPressureTensor.cx+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.z*Cations[CurrentSystem][i].Atoms[A].Velocity.x;
          if((!Fixed.z)&&(!Fixed.y)) KineticPressureTensor.cy+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.z*Cations[CurrentSystem][i].Atoms[A].Velocity.y;
          if(!Fixed.z)               KineticPressureTensor.cz+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.z*Cations[CurrentSystem][i].Atoms[A].Velocity.z;
        }
      }
    }
  }

  return KineticPressureTensor;
}


REAL GetCoreShellTemperature(void)
{
  int f1;
  int A,B,total_nr_core_shells;
  REAL MassA,MassB;
  REAL CoreShellKineticEnergy,ccc,sss,ppp;
  VECTOR VelA,VelB;
  int nr_core_shells,nr_atoms;

  CoreShellKineticEnergy=0;
  total_nr_core_shells=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    nr_core_shells=Framework[CurrentSystem].NumberOfCoreShells[f1];
    total_nr_core_shells+=nr_core_shells;
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    if(nr_core_shells>0)
    {
      // loop over shells
      for(A=0;A<nr_atoms-nr_core_shells;A++)
      {
        // get the core for this shell
        B=Framework[CurrentSystem].CoreShellConnectivity[f1][A];
        if(B>0)
        {
          MassA=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Mass;
          MassB=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Mass;
          VelA=Framework[CurrentSystem].Atoms[f1][A].Velocity;
          VelB=Framework[CurrentSystem].Atoms[f1][B].Velocity;
          ppp=(SQR(MassA*VelA.x+MassB*VelB.x)+SQR(MassA*VelA.y+MassB*VelB.y)+SQR(MassA*VelA.z+MassB*VelB.z))/(MassA+MassB);
          ccc=MassA*(SQR(VelA.x)+SQR(VelA.y)+SQR(VelA.z));
          sss=MassB*(SQR(VelB.x)+SQR(VelB.y)+SQR(VelB.z));
          CoreShellKineticEnergy+=0.5*(ccc+sss-ppp);

        }
      }
    }
  }
  return CoreShellKineticEnergy/total_nr_core_shells;
}

void AdjustCoreShellVelocities(void)
{
  int f1;
  int A,B;
  REAL MassA,MassB,pke,scale,rmu;
  VECTOR VelA,VelB,DeltaV,comv;
  int nr_core_shells,nr_atoms;

  pke=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*1.0e-5;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    nr_core_shells=Framework[CurrentSystem].NumberOfCoreShells[f1];
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    if(nr_core_shells>0)
    {
      // loop over shells
      for(A=0;A<nr_atoms-nr_core_shells;A++)
      {
        // get the core for this shell
        B=Framework[CurrentSystem].CoreShellConnectivity[f1][A];
        if(B>0)
        {
          MassA=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Mass;
          MassB=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Mass;
          rmu=(MassA*MassB)/(MassA+MassB);
          VelA=Framework[CurrentSystem].Atoms[f1][A].Velocity;
          VelB=Framework[CurrentSystem].Atoms[f1][B].Velocity;
          DeltaV.x=VelB.x-VelA.x;
          DeltaV.y=VelB.y-VelA.y;
          DeltaV.z=VelB.z-VelA.z;
          scale=sqrt(pke/(rmu*(SQR(DeltaV.x)+SQR(DeltaV.y)+SQR(DeltaV.z))));

          comv.x=(MassA*VelA.x+MassB*VelB.x)/(MassA+MassB);
          comv.y=(MassA*VelA.y+MassB*VelB.y)/(MassA+MassB);
          comv.z=(MassA*VelA.z+MassB*VelB.z)/(MassA+MassB);

          VelA.x=comv.x-scale*rmu*DeltaV.x/MassA;
          VelA.y=comv.y-scale*rmu*DeltaV.y/MassA;
          VelA.z=comv.z-scale*rmu*DeltaV.z/MassA;
          Framework[CurrentSystem].Atoms[f1][A].Velocity=VelA;

          VelB.x=comv.x+scale*rmu*DeltaV.x/MassB;
          VelB.y=comv.y+scale*rmu*DeltaV.y/MassB;
          VelB.z=comv.z+scale*rmu*DeltaV.z/MassB;
          Framework[CurrentSystem].Atoms[f1][B].Velocity=VelB;
        }
      }
    }
  }
}


void InitializeBoxVelocities(void)
{
  switch(Ensemble[CurrentSystem])
  {
    case NPH:
      LnVolumeVelocity[CurrentSystem]=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/LnVolumeMass[CurrentSystem]);
      break;
    case NPHPR:
      CellVelocity[CurrentSystem].ax=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].ay=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].az=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].bx=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].by=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].bz=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].cx=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].cy=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].cz=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      break;
    default:
      break;
  }
}



void InitializeNoseHooverCurrentSystem(void)
{
  int i;
  int M,N;

  M=therm_baro_stats.ThermostatChainLength;
  N=therm_baro_stats.BarostatChainLength;

  CellMass[CurrentSystem]=((DegreesOfFreedomTranslation[CurrentSystem]+3.0)/3.0)*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
                          SQR(therm_baro_stats.time_scale_parameter_barostat);

  LnVolumeMass[CurrentSystem]=(DegreesOfFreedom[CurrentSystem]+3.0)*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
                    SQR(therm_baro_stats.time_scale_parameter_barostat);

  LnVolumePosition[CurrentSystem]=log(Volume[CurrentSystem]);

  switch(Ensemble[CurrentSystem])
  {
    case NPTPR:
    case MuPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          // 6 degrees of freedom for a fully flexible box
          BarostatDegreesOfFreedom[CurrentSystem][0]=6.0*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          // 4 degrees of freedom for a fully flexible box
          BarostatDegreesOfFreedom[CurrentSystem][0]=4.0*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
          break;
        case ANISOTROPIC:
          // 3 degrees of freedom for an anisotropic box
          BarostatDegreesOfFreedom[CurrentSystem][0]=3.0*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
          break;
        case ISOTROPIC:
          // 1 degree of freedom for an isotropic box
          BarostatDegreesOfFreedom[CurrentSystem][0]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
          break;
      }
      break;
    case NPT:
    case MuPT:
      // 1 degree of freedom for an isotropic box
      BarostatDegreesOfFreedom[CurrentSystem][0]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
      break;
    default:
      BarostatDegreesOfFreedom[CurrentSystem][0]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
      break;
  }

  for(i=1;i<N;i++)
    BarostatDegreesOfFreedom[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

  for(i=0;i<N;i++)
    BarostatMass[CurrentSystem][i]=BarostatDegreesOfFreedom[CurrentSystem][i]*SQR(therm_baro_stats.time_scale_parameter_barostat);

  // general
  ThermostatDegreesOfFreedom[CurrentSystem][0]=DegreesOfFreedom[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
  for(i=1;i<M;i++)
    ThermostatDegreesOfFreedom[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

  for(i=0;i<M;i++)
    ThermostatMass[CurrentSystem][i]=ThermostatDegreesOfFreedom[CurrentSystem][i]*SQR(therm_baro_stats.time_scale_parameter_thermostat);

  // translation framework
  ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][0]=DegreesOfFreedomFramework[CurrentSystem]*
                              K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
  for(i=1;i<M;i++)
    ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

  for(i=0;i<M;i++)
    ThermostatMassTranslationFramework[CurrentSystem][i]=ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][i]*
                                                           SQR(therm_baro_stats.time_scale_parameter_thermostat);

  // translation adsorbates
  ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][0]=DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]*
                                                                      K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

  for(i=1;i<M;i++)
    ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

  for(i=0;i<M;i++)
    ThermostatMassTranslationAdsorbates[CurrentSystem][i]=ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][i]*
                                                            SQR(therm_baro_stats.time_scale_parameter_thermostat);

  // translation cations
  ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][0]=DegreesOfFreedomTranslationalCations[CurrentSystem]*
                                                                   K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
  for(i=1;i<M;i++)
    ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

  for(i=0;i<M;i++)
    ThermostatMassTranslationCations[CurrentSystem][i]=ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][i]*
                                                         SQR(therm_baro_stats.time_scale_parameter_thermostat);

  // translation
  ThermostatDegreesOfFreedomTranslation[CurrentSystem][0]=DegreesOfFreedomTranslation[CurrentSystem]*
                                                K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
  for(i=1;i<M;i++)
    ThermostatDegreesOfFreedomTranslation[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

  for(i=0;i<M;i++)
    ThermostatMassTranslation[CurrentSystem][i]=ThermostatDegreesOfFreedomTranslation[CurrentSystem][i]*
                                                    SQR(therm_baro_stats.time_scale_parameter_thermostat);

  // rotation
  ThermostatDegreesOfFreedomRotation[CurrentSystem][0]=DegreesOfFreedomRotation[CurrentSystem]*
                                                K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
  for(i=1;i<M;i++)
    ThermostatDegreesOfFreedomRotation[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

  for(i=0;i<M;i++)
  {
    if(ThermostatDegreesOfFreedomRotation[CurrentSystem][i]<1e-8)
       ThermostatMassRotation[CurrentSystem][i]=1.0;
    else
      ThermostatMassRotation[CurrentSystem][i]=ThermostatDegreesOfFreedomRotation[CurrentSystem][i]*
                                                  SQR(therm_baro_stats.time_scale_parameter_thermostat);
  }

  // initialize the cell velocities
  switch(Ensemble[CurrentSystem])
  {
    case NPH:
      LnVolumeVelocity[CurrentSystem]=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/LnVolumeMass[CurrentSystem]);
      break;
    case NPTPR:
    case MuPTPR:
    case NPHPR:
      CellVelocity[CurrentSystem].ax=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].ay=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].az=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].bx=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].by=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].bz=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].cx=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].cy=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
      CellVelocity[CurrentSystem].cz=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
    default:
      break;
  }

}

void InitializeNoseHooverAllSystems(void)
{
  int i,k;
  int M,N;
  REAL UKineticCell;
  REAL avg;

  M=therm_baro_stats.ThermostatChainLength;
  N=therm_baro_stats.BarostatChainLength;

  if(M<1)
  {
    fprintf(stderr, "Error: Number of thermostats in the NHC chain %d has to be larger than 0\n",M);
    exit(0);
  }

  if(M>=MAXIMUM_LENGTH_THERMOSTATS)
  {
    fprintf(stderr, "Number of thermostats in the NHC chain %d exceed the limit of %d\n",
      M,MAXIMUM_LENGTH_THERMOSTATS);
    exit(0);
  }

  if(N<1)
  {
    fprintf(stderr, "Error: Number of barostats in the NHC chain %d has to be larger than 0\n",N);
    exit(0);
  }

  if(N>=MAXIMUM_LENGTH_BAROSTATS)
  {
    fprintf(stderr, "Number of barostats in the NHC chain %d exceed the limit of %d\n",
      M,MAXIMUM_LENGTH_BAROSTATS);
    exit(0);
  }

  switch(therm_baro_stats.NumberOfYoshidaSuzukiSteps)
  {
    case 1:
      w[0]=1.0;
      break;
    case 3:
      w[0]=1.0/(2.0-pow(2.0,1.0/3.0));
      //w[1]=-pow(2.0,1.0/3.0)/(2.0-pow(2.0,1.0/3.0));
      w[1]=1.0-2.0*w[0];
      w[2]=w[0];
      break;
    case 5:
      w[0]=1.0/(4.0-pow(4.0,1.0/3.0));
      w[1]=w[0];
      w[2]=1.0-4.0*w[0];
      w[3]=w[0];
      w[4]=w[0];
      break;
    case 7:
      w[0]=0.784513610477560;
      w[1]=0.235573213359357;
      w[2]=-1.17767998417887;
      w[3]=1.0-2.0*(w[0]+w[1]+w[2]);
      w[4]=-1.17767998417887;
      w[5]=0.235573213359357;
      w[6]=0.784513610477560;
      break;
    case 9:
      w[0]=0.192;
      w[1]=0.554910818409783619692725006662999;
      w[2]=0.124659619941888644216504240951585;
      w[3]=-0.843182063596933505315033808282941;
      w[4]=1.0-2.0*(w[0]+w[1]+w[2]+w[3]);
      w[5]=-0.843182063596933505315033808282941;
      w[6]=0.124659619941888644216504240951585;
      w[7]=0.554910818409783619692725006662999;
      w[8]=0.192;
      break;
    default:
      fprintf(stderr, "Error: Yoshida-Suzuki-steps should be: 1,3,5,7 or 9\n");
      exit(0);
      break;
  }

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    UNoseHoover[CurrentSystem]=0.0;
    UNoseHooverAdsorbates[CurrentSystem]=0.0;
    UNoseHooverCations[CurrentSystem]=0.0;
    UNoseHooverFramework[CurrentSystem]=0.0;

    CellMass[CurrentSystem]=((DegreesOfFreedomTranslation[CurrentSystem]+3.0)/3.0)*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
                 SQR(therm_baro_stats.time_scale_parameter_barostat);

    LnVolumeMass[CurrentSystem]=(DegreesOfFreedomTranslation[CurrentSystem]+3.0)*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
                    SQR(therm_baro_stats.time_scale_parameter_barostat);

    LnVolumePosition[CurrentSystem]=log(Volume[CurrentSystem]);

    switch(Ensemble[CurrentSystem])
    {
      case NPTPR:
      case MuPTPR:
      case NPHPR:
        switch(NPTPRCellType[CurrentSystem])
        {
          case REGULAR:
          case REGULAR_UPPER_TRIANGLE:
            // 6 degrees of freedom for a fully flexible box
            BarostatDegreesOfFreedom[CurrentSystem][0]=6.0*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
            break;
          case MONOCLINIC:
          case MONOCLINIC_UPPER_TRIANGLE:
            // 4 degrees of freedom for a fully flexible box
            BarostatDegreesOfFreedom[CurrentSystem][0]=4.0*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
            break;
          case ANISOTROPIC:
            // 3 degrees of freedom for an anisotropic box
            BarostatDegreesOfFreedom[CurrentSystem][0]=3.0*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
            break;
          case ISOTROPIC:
            // 1 degree of freedom for an isotropic box
            BarostatDegreesOfFreedom[CurrentSystem][0]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
            break;
        }
        break;
      case NPT:
      case MuPT:
        // 1 degree of freedom for an isotropic box
        BarostatDegreesOfFreedom[CurrentSystem][0]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
        break;
      default:
        BarostatDegreesOfFreedom[CurrentSystem][0]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
        break;
    }

    for(i=1;i<N;i++)
      BarostatDegreesOfFreedom[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

    for(i=0;i<N;i++)
      BarostatMass[CurrentSystem][i]=BarostatDegreesOfFreedom[CurrentSystem][i]*SQR(therm_baro_stats.time_scale_parameter_barostat);


    // general
    ThermostatDegreesOfFreedom[CurrentSystem][0]=DegreesOfFreedom[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
    for(i=1;i<M;i++)
      ThermostatDegreesOfFreedom[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

    for(i=0;i<M;i++)
      ThermostatMass[CurrentSystem][i]=ThermostatDegreesOfFreedom[CurrentSystem][i]*SQR(therm_baro_stats.time_scale_parameter_thermostat);

    // translation framework
    ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][0]=DegreesOfFreedomFramework[CurrentSystem]*
                              K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
    for(i=1;i<M;i++)
      ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

    for(i=0;i<M;i++)
      ThermostatMassTranslationFramework[CurrentSystem][i]=ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][i]*
                                                           SQR(therm_baro_stats.time_scale_parameter_thermostat);

    // translation adsorbates
    ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][0]=DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]*
                                                                      K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
    for(i=1;i<M;i++)
      ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

    for(i=0;i<M;i++)
      ThermostatMassTranslationAdsorbates[CurrentSystem][i]=ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][i]*
                                                            SQR(therm_baro_stats.time_scale_parameter_thermostat);

    // translation cations
    ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][0]=DegreesOfFreedomTranslationalCations[CurrentSystem]*
                                                                   K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
    for(i=1;i<M;i++)
      ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

    for(i=0;i<M;i++)
      ThermostatMassTranslationCations[CurrentSystem][i]=ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][i]*
                                                         SQR(therm_baro_stats.time_scale_parameter_thermostat);


    // translation
    ThermostatDegreesOfFreedomTranslation[CurrentSystem][0]=DegreesOfFreedomTranslation[CurrentSystem]*
                                                K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
    for(i=1;i<M;i++)
      ThermostatDegreesOfFreedomTranslation[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

    for(i=0;i<M;i++)
      ThermostatMassTranslation[CurrentSystem][i]=ThermostatDegreesOfFreedomTranslation[CurrentSystem][i]*
                                                    SQR(therm_baro_stats.time_scale_parameter_thermostat);

    // rotation
    ThermostatDegreesOfFreedomRotation[CurrentSystem][0]=DegreesOfFreedomRotation[CurrentSystem]*
                                                K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];
    for(i=1;i<M;i++)
      ThermostatDegreesOfFreedomRotation[CurrentSystem][i]=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem];

    for(i=0;i<M;i++)
    {
      if(ThermostatDegreesOfFreedomRotation[CurrentSystem][i]<1e-8)
         ThermostatMassRotation[CurrentSystem][i]=1.0;
      else
        ThermostatMassRotation[CurrentSystem][i]=ThermostatDegreesOfFreedomRotation[CurrentSystem][i]*
                                                    SQR(therm_baro_stats.time_scale_parameter_thermostat);
    }

    for(k=0;k<M;k++)
    {
      ThermostatForce[CurrentSystem][0]=0.0;
      ThermostatVelocity[CurrentSystem][0]=0.0;
      ThermostatPosition[CurrentSystem][0]=0.0;
      ThermostatForceTranslation[CurrentSystem][0]=0.0;
      ThermostatForceRotation[CurrentSystem][0]=0.0;
      ThermostatVelocityTranslation[CurrentSystem][0]=0.0;
      ThermostatVelocityRotation[CurrentSystem][0]=0.0;
      ThermostatPositionTranslation[CurrentSystem][0]=0.0;
      ThermostatPositionRotation[CurrentSystem][0]=0.0;
    }
    for(k=0;k<N;k++)
    {
      BarostatForce[CurrentSystem][0]=0.0;
      BarostatVelocity[CurrentSystem][0]=0.0;
      BarostatPosition[CurrentSystem][0]=0.0;
    }

    for(k=0;k<M;k++)
    {
      ThermostatVelocity[CurrentSystem][k]=RandomGaussianNumber()*sqrt(ThermostatDegreesOfFreedom[CurrentSystem][k]/ThermostatMass[CurrentSystem][k]);
      ThermostatVelocityTranslation[CurrentSystem][k]=RandomGaussianNumber()*sqrt(ThermostatDegreesOfFreedomTranslation[CurrentSystem][k]/ThermostatMassTranslation[CurrentSystem][k]);
      ThermostatVelocityRotation[CurrentSystem][k]=RandomGaussianNumber()*sqrt(ThermostatDegreesOfFreedomRotation[CurrentSystem][k]/ThermostatMassTranslation[CurrentSystem][k]);

      ThermostatVelocityTranslationFramework[CurrentSystem][k]=RandomGaussianNumber()*sqrt(ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][k]/ThermostatMassTranslationFramework[CurrentSystem][k]);
      ThermostatVelocityTranslationAdsorbates[CurrentSystem][k]=RandomGaussianNumber()*sqrt(ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][k]/ThermostatMassTranslationAdsorbates[CurrentSystem][k]);
      ThermostatVelocityTranslationCations[CurrentSystem][k]=RandomGaussianNumber()*sqrt(ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][k]/ThermostatMassTranslationCations[CurrentSystem][k]);
    }

    for(k=0;k<N;k++)
      BarostatVelocity[CurrentSystem][k]=RandomGaussianNumber()*sqrt(BarostatDegreesOfFreedom[CurrentSystem][k]/BarostatMass[CurrentSystem][k]);


    // initialize the cell velocities
    switch(Ensemble[CurrentSystem])
    {
      case NPT:
      case MuPT:
      case NPH:
        LnVolumeVelocity[CurrentSystem]=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/LnVolumeMass[CurrentSystem]);
        break;
      case NPTPR:
      case MuPTPR:
      case NPHPR:
        switch(NPTPRCellType[CurrentSystem])
        {
          case ISOTROPIC:
            CellVelocity[CurrentSystem].ay=CellVelocity[CurrentSystem].az=CellVelocity[CurrentSystem].bx=0.0;
            CellVelocity[CurrentSystem].bz=CellVelocity[CurrentSystem].cx=CellVelocity[CurrentSystem].cz=0.0;
            CellVelocity[CurrentSystem].ax=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].by=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].cz=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            avg=(CellVelocity[CurrentSystem].ax+CellVelocity[CurrentSystem].by+CellVelocity[CurrentSystem].cz)/3.0;
            CellVelocity[CurrentSystem].ax=avg;
            CellVelocity[CurrentSystem].by=avg;
            CellVelocity[CurrentSystem].cz=avg;
            break;
          case ANISOTROPIC:
            CellVelocity[CurrentSystem].ax=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].by=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].cz=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            break;
          case REGULAR:
          default:
            CellVelocity[CurrentSystem].ax=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].ay=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].az=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].bx=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].by=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].bz=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].cx=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].cy=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            CellVelocity[CurrentSystem].cz=RandomGaussianNumber()*sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/CellMass[CurrentSystem]);
            break;
        }

        ThermostatForceTranslationFramework[CurrentSystem][0]=
             (2.0*GetTranslationKineticEnergyFramework()-ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][0])/
                    ThermostatMassTranslationFramework[CurrentSystem][0];
        for(k=0;k<M-1;k++)
        {
          ThermostatForceTranslationFramework[CurrentSystem][k+1]=
              (ThermostatMassTranslationFramework[CurrentSystem][k]*SQR(ThermostatVelocityTranslationFramework[CurrentSystem][k])-
               ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][k+1])/ThermostatMassTranslationFramework[CurrentSystem][k+1];
        }
        ThermostatForceTranslationAdsorbates[CurrentSystem][0]=
             (2.0*GetTranslationKineticEnergyAdsorbates()-ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][0])/
                    ThermostatMassTranslationAdsorbates[CurrentSystem][0];
        for(k=0;k<M-1;k++)
        {
          ThermostatForceTranslationAdsorbates[CurrentSystem][k+1]=
              (ThermostatMassTranslationAdsorbates[CurrentSystem][k]*SQR(ThermostatVelocityTranslationAdsorbates[CurrentSystem][k])-
               ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][k+1])/ThermostatMassTranslationAdsorbates[CurrentSystem][k+1];
        }

        UKineticCell=CellMass[CurrentSystem]*
           (SQR(CellVelocity[CurrentSystem].ax)+SQR(CellVelocity[CurrentSystem].ay)+SQR(CellVelocity[CurrentSystem].az)+
            SQR(CellVelocity[CurrentSystem].bx)+SQR(CellVelocity[CurrentSystem].by)+SQR(CellVelocity[CurrentSystem].bz)+
            SQR(CellVelocity[CurrentSystem].cx)+SQR(CellVelocity[CurrentSystem].cy)+SQR(CellVelocity[CurrentSystem].cz));

        BarostatForce[CurrentSystem][0]=(UKineticCell-BarostatDegreesOfFreedom[CurrentSystem][0])/BarostatMass[CurrentSystem][0];

        for(k=0;k<N-1;k++)
          BarostatForce[CurrentSystem][k+1]=(BarostatMass[CurrentSystem][k]*SQR(BarostatVelocity[CurrentSystem][k])-
                                    BarostatDegreesOfFreedom[CurrentSystem][k+1])/BarostatMass[CurrentSystem][k+1];
      default:
        break;
  }

  }
}

void ComputeNoseHooverEnergySystem(void)
{
  int i,M,N;
  REAL ExtPressure;

  M=therm_baro_stats.ThermostatChainLength;
  N=therm_baro_stats.BarostatChainLength;
  ExtPressure=therm_baro_stats.ExternalPressure[CurrentSystem][0];

  UNoseHoover[CurrentSystem]=0.0;
  switch(Ensemble[CurrentSystem])
  {
    case MuPT:
    case NPT:
      UNoseHoover[CurrentSystem]=ExtPressure*Volume[CurrentSystem]+
               0.5*LnVolumeMass[CurrentSystem]*SQR(LnVolumeVelocity[CurrentSystem]);

      for(i=0;i<M;i++)
          UNoseHoover[CurrentSystem]+=0.5*SQR(BarostatVelocity[CurrentSystem][i])*BarostatMass[CurrentSystem][i]+
                 K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*BarostatPosition[CurrentSystem][i];
      // fall through
    case NVT:
    case MuVT:
      if(DegreesOfFreedomFramework[CurrentSystem]>0)
      {
        UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationFramework[CurrentSystem][0])*
            ThermostatMassTranslationFramework[CurrentSystem][0]+
            DegreesOfFreedomFramework[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
            ThermostatPositionTranslationFramework[CurrentSystem][0];
        for(i=1;i<M;i++)
          UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationFramework[CurrentSystem][i])*
                 ThermostatMassTranslationFramework[CurrentSystem][i]+
                 K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionTranslationFramework[CurrentSystem][i];
      }

      if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
      {
        UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationAdsorbates[CurrentSystem][0])*
            ThermostatMassTranslationAdsorbates[CurrentSystem][0]+
            DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
            ThermostatPositionTranslationAdsorbates[CurrentSystem][0];

        for(i=1;i<M;i++)
          UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationAdsorbates[CurrentSystem][i])*
               ThermostatMassTranslationAdsorbates[CurrentSystem][i]+
                K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionTranslationAdsorbates[CurrentSystem][i];
      }

      if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
      {
        UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationCations[CurrentSystem][0])*
                ThermostatMassTranslationCations[CurrentSystem][0]+
                DegreesOfFreedomTranslationalCations[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
                ThermostatPositionTranslationCations[CurrentSystem][0];

        for(i=1;i<M;i++)
          UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationCations[CurrentSystem][i])*
                ThermostatMassTranslationCations[CurrentSystem][i]+
                K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionTranslationCations[CurrentSystem][i];
      }

      if(DegreesOfFreedomRotation[CurrentSystem]>0)
      {
         UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityRotation[CurrentSystem][0])*ThermostatMassRotation[CurrentSystem][0]+
                     DegreesOfFreedomRotation[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
                     ThermostatPositionRotation[CurrentSystem][0];
         for(i=1;i<M;i++)
           UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityRotation[CurrentSystem][i])*ThermostatMassRotation[CurrentSystem][i]+
                 K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionRotation[CurrentSystem][i];
      }
      break;
    case NPTPR:
    case MuPTPR:
      UNoseHoover[CurrentSystem]=therm_baro_stats.ExternalPressure[CurrentSystem][0]*Volume[CurrentSystem];

      UNoseHoover[CurrentSystem]+=0.5*CellMass[CurrentSystem]*
          (SQR(CellVelocity[CurrentSystem].ax)+SQR(CellVelocity[CurrentSystem].ay)+SQR(CellVelocity[CurrentSystem].az)+
           SQR(CellVelocity[CurrentSystem].bx)+SQR(CellVelocity[CurrentSystem].by)+SQR(CellVelocity[CurrentSystem].bz)+
           SQR(CellVelocity[CurrentSystem].cx)+SQR(CellVelocity[CurrentSystem].cy)+SQR(CellVelocity[CurrentSystem].cz));

      for(i=0;i<N;i++)
        UNoseHoover[CurrentSystem]+=0.5*SQR(BarostatVelocity[CurrentSystem][i])*BarostatMass[CurrentSystem][i]+
            BarostatDegreesOfFreedom[CurrentSystem][i]*BarostatPosition[CurrentSystem][i];

      if(DegreesOfFreedomFramework[CurrentSystem]>0)
      {
        UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationFramework[CurrentSystem][0])*
            ThermostatMassTranslationFramework[CurrentSystem][0]+
            DegreesOfFreedomFramework[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
            ThermostatPositionTranslationFramework[CurrentSystem][0];
        for(i=1;i<M;i++)
          UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationFramework[CurrentSystem][i])*
                 ThermostatMassTranslationFramework[CurrentSystem][i]+
                 K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionTranslationFramework[CurrentSystem][i];
      }

      if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
      {
        UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationAdsorbates[CurrentSystem][0])*
            ThermostatMassTranslationAdsorbates[CurrentSystem][0]+
            DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
            ThermostatPositionTranslationAdsorbates[CurrentSystem][0];

        for(i=1;i<M;i++)
          UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationAdsorbates[CurrentSystem][i])*
               ThermostatMassTranslationAdsorbates[CurrentSystem][i]+
                K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionTranslationAdsorbates[CurrentSystem][i];
      }

      if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
      {
        UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationCations[CurrentSystem][0])*
                ThermostatMassTranslationCations[CurrentSystem][0]+
                DegreesOfFreedomTranslationalCations[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
                ThermostatPositionTranslationCations[CurrentSystem][0];

        for(i=1;i<M;i++)
          UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationCations[CurrentSystem][i])*
                ThermostatMassTranslationCations[CurrentSystem][i]+
                K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionTranslationCations[CurrentSystem][i];
      }

      if(DegreesOfFreedomRotation[CurrentSystem]>0)
      {
         UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityRotation[CurrentSystem][0])*ThermostatMassRotation[CurrentSystem][0]+
                     DegreesOfFreedomRotation[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
                     ThermostatPositionRotation[CurrentSystem][0];
         for(i=1;i<M;i++)
           UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityRotation[CurrentSystem][i])*ThermostatMassRotation[CurrentSystem][i]+
                 K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionRotation[CurrentSystem][i];
      }
      break;
    case NPH:
      UNoseHoover[CurrentSystem]=ExtPressure*Volume[CurrentSystem]+
               0.5*LnVolumeMass[CurrentSystem]*SQR(LnVolumeVelocity[CurrentSystem]);
      break;
    case NPHPR:
      UNoseHoover[CurrentSystem]=ExtPressure*Volume[CurrentSystem];

      UNoseHoover[CurrentSystem]+=0.5*CellMass[CurrentSystem]*
          (SQR(CellVelocity[CurrentSystem].ax)+SQR(CellVelocity[CurrentSystem].ay)+SQR(CellVelocity[CurrentSystem].az)+
           SQR(CellVelocity[CurrentSystem].bx)+SQR(CellVelocity[CurrentSystem].by)+SQR(CellVelocity[CurrentSystem].bz)+
           SQR(CellVelocity[CurrentSystem].cx)+SQR(CellVelocity[CurrentSystem].cy)+SQR(CellVelocity[CurrentSystem].cz));
      break;
    default:
      break;
  }
}

void UpdateVelocities(void)
{
  int i,k,l,A,Type,f1;
  REAL E2,E4,E6,E8,aa,bb,aa2,arg2,mass,alpha;
  VECTOR vel,force;
  INT_VECTOR3 Fixed;

  // Taylor-expansion coefficients of sinh(x)/x
  E2=1.0/6.0;
  E4=E2/20.0;
  E6=E4/42.0;
  E8=E6/72.0;

  switch(Ensemble[CurrentSystem])
  {
    case NPH:
    case NPT:
    case MuPT:
      alpha=(1.0+3.0/DegreesOfFreedomTranslation[CurrentSystem]);
      aa=exp(-0.25*DeltaT*alpha*LnVolumeVelocity[CurrentSystem]);
      aa2=SQR(aa);
      arg2=SQR(0.25*DeltaT*alpha*LnVolumeVelocity[CurrentSystem]);
      bb=aa*((((E8*arg2+E6)*arg2+E4)*arg2+E2)*arg2+1.0)*DeltaT;
      break;
    case NVE:
    case NVT:
    case MuVT:
    case NPTPR:
    case MuPTPR:
    case NPHPR:
    default:
      aa=1.0;
      aa2=1.0;
      bb=DeltaT;
      break;
  }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      vel=Framework[CurrentSystem].Atoms[f1][i].Velocity;
      force=Framework[CurrentSystem].Atoms[f1][i].Force;
      mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
      Fixed=Framework[CurrentSystem].Atoms[f1][i].Fixed;
      if(!Fixed.x) vel.x=vel.x*aa2+0.5*bb*force.x/mass;
      if(!Fixed.y) vel.y=vel.y*aa2+0.5*bb*force.y/mass;
      if(!Fixed.z) vel.z=vel.z*aa2+0.5*bb*force.z/mass;
      Framework[CurrentSystem].Atoms[f1][i].Velocity=vel;
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        mass=Components[Type].Groups[l].Mass;
        vel=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
        force=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassForce;
        Fixed=Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x) vel.x=vel.x*aa2+0.5*bb*force.x/mass;
        if(!Fixed.y) vel.y=vel.y*aa2+0.5*bb*force.y/mass;
        if(!Fixed.z) vel.z=vel.z*aa2+0.5*bb*force.z/mass;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity=vel;

        Fixed=Adsorbates[CurrentSystem][i].Groups[l].FixedOrientation;
        if((!Fixed.x)||(!Fixed.y)||(!Fixed.z))
        {
          Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.r+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.r;
          Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.i+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.i;
          Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.j+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.j;
          Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.k+=0.5*DeltaT*Adsorbates[CurrentSystem][i].Groups[l].QuaternionForce.k;
        }
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];

          vel=Adsorbates[CurrentSystem][i].Atoms[A].Velocity;
          force=Adsorbates[CurrentSystem][i].Atoms[A].Force;
          mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
          Fixed=Adsorbates[CurrentSystem][i].Atoms[A].Fixed;
          if(!Fixed.x) vel.x=vel.x*aa2+0.5*bb*force.x/mass;
          if(!Fixed.y) vel.y=vel.y*aa2+0.5*bb*force.y/mass;
          if(!Fixed.z) vel.z=vel.z*aa2+0.5*bb*force.z/mass;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity=vel;
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
        mass=Components[Type].Groups[l].Mass;
        vel=Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
        force=Cations[CurrentSystem][i].Groups[l].CenterOfMassForce;
        Fixed=Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x) vel.x=vel.x*aa2+0.5*bb*force.x/mass;
        if(!Fixed.y) vel.y=vel.y*aa2+0.5*bb*force.y/mass;
        if(!Fixed.z) vel.z=vel.z*aa2+0.5*bb*force.z/mass;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity=vel;

        // FIX
        Fixed=Cations[CurrentSystem][i].Groups[l].FixedOrientation;
        if((!Fixed.x)||(!Fixed.y)||(!Fixed.z))
        {
          Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.r+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.r;
          Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.i+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.i;
          Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.j+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.j;
          Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.k+=0.5*DeltaT*Cations[CurrentSystem][i].Groups[l].QuaternionForce.k;
        }
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];

          vel=Cations[CurrentSystem][i].Atoms[A].Velocity;
          force=Cations[CurrentSystem][i].Atoms[A].Force;
          mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[A].Type].Mass;
          Fixed=Cations[CurrentSystem][i].Atoms[A].Fixed;
          if(!Fixed.x) vel.x=vel.x*aa2+0.5*bb*force.x/mass;
          if(!Fixed.y) vel.y=vel.y*aa2+0.5*bb*force.y/mass;
          if(!Fixed.z) vel.z=vel.z*aa2+0.5*bb*force.z/mass;
          Cations[CurrentSystem][i].Atoms[A].Velocity=vel;
        }
      }
    }
  }
}


void UpdatePositions(void)
{
  int i,k,l,A,Type,f1;
  REAL E2,E4,E6,E8,aa,bb,aa2,arg2,det;
  VECTOR pos,vel;
  INT_VECTOR3 Fixed;

  E2=1.0/6.0;
  E4=E2/20.0;
  E6=E4/42.0;
  E8=E6/72.0;

  switch(Ensemble[CurrentSystem])
  {
    case NPH:
    case NPT:
    case MuPT:
      aa=exp(0.5*DeltaT*LnVolumeVelocity[CurrentSystem]);
      aa2=SQR(aa);
      arg2=SQR(0.5*DeltaT*LnVolumeVelocity[CurrentSystem]);
      bb=aa*((((E8*arg2+E6)*arg2+E4)*arg2+E2)*arg2+1.0)*DeltaT;
      break;
    case NVE:
    case NVT:
    case MuVT:
    default:
      aa=1.0;
      aa2=1.0;
      bb=DeltaT;
      break;
  }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      pos=Framework[CurrentSystem].Atoms[f1][i].Position;
      Framework[CurrentSystem].Atoms[f1][i].RattleReferencePosition=pos;
      vel=Framework[CurrentSystem].Atoms[f1][i].Velocity;
      Fixed=Framework[CurrentSystem].Atoms[f1][i].Fixed;
      if(!Fixed.x) pos.x=pos.x*aa2+vel.x*bb;
      if(!Fixed.y) pos.y=pos.y*aa2+vel.y*bb;
      if(!Fixed.z) pos.z=pos.z*aa2+vel.z*bb;
      Framework[CurrentSystem].Atoms[f1][i].Position=pos;
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        pos=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition;
        vel=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
        Fixed=Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x) pos.x=pos.x*aa2+vel.x*bb;
        if(!Fixed.y) pos.y=pos.y*aa2+vel.y*bb;
        if(!Fixed.z) pos.z=pos.z*aa2+vel.z*bb;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition=pos;
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];

          pos=Adsorbates[CurrentSystem][i].Atoms[A].Position;
          Adsorbates[CurrentSystem][i].Atoms[A].RattleReferencePosition=pos;
          vel=Adsorbates[CurrentSystem][i].Atoms[A].Velocity;
          Fixed=Adsorbates[CurrentSystem][i].Atoms[A].Fixed;
          if(!Fixed.x) pos.x=pos.x*aa2+vel.x*bb;
          if(!Fixed.y) pos.y=pos.y*aa2+vel.y*bb;
          if(!Fixed.z) pos.z=pos.z*aa2+vel.z*bb;
          Adsorbates[CurrentSystem][i].Atoms[A].Position=pos;
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
        pos=Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition;
        vel=Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
        Fixed=Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x) pos.x=pos.x*aa2+vel.x*bb;
        if(!Fixed.y) pos.y=pos.y*aa2+vel.y*bb;
        if(!Fixed.z) pos.z=pos.z*aa2+vel.z*bb;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition=pos;
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];

          pos=Cations[CurrentSystem][i].Atoms[A].Position;
          Cations[CurrentSystem][i].Atoms[A].RattleReferencePosition=pos;
          vel=Cations[CurrentSystem][i].Atoms[A].Velocity;
          Fixed=Cations[CurrentSystem][i].Atoms[A].Fixed;
          if(!Fixed.x) pos.x=pos.x*aa2+vel.x*bb;
          if(!Fixed.y) pos.y=pos.y*aa2+vel.y*bb;
          if(!Fixed.z) pos.z=pos.z*aa2+vel.z*bb;
          Cations[CurrentSystem][i].Atoms[A].Position=pos;
        }
      }
    }
  }

  if((Ensemble[CurrentSystem]==NPT)||(Ensemble[CurrentSystem]==MuPT)||(Ensemble[CurrentSystem]==NPH))
  {
    LnVolumePosition[CurrentSystem]+=LnVolumeVelocity[CurrentSystem]*DeltaT;

    Box[CurrentSystem].ax*=aa2;  Box[CurrentSystem].ay*=aa2;  Box[CurrentSystem].az*=aa2;
    Box[CurrentSystem].bx*=aa2;  Box[CurrentSystem].by*=aa2;  Box[CurrentSystem].bz*=aa2;
    Box[CurrentSystem].cx*=aa2;  Box[CurrentSystem].cy*=aa2;  Box[CurrentSystem].cz*=aa2;
    CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
    Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
    CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);

    if(MIN3(BoxProperties[CurrentSystem].cx,BoxProperties[CurrentSystem].cy,
            BoxProperties[CurrentSystem].cz)<2.0*CutOffVDW)
    {
       fprintf(stderr, "ERROR: (System (%d) Cutoff smaller than half of one of the perpendicular boxlengths !!!\n",CurrentSystem);
       fprintf(stderr, "       Cutoff: %lf perpendicular boxlengths: %lf %lf %lf\n",(double)CutOffVDW,
               (double)BoxProperties[CurrentSystem].cx,(double)BoxProperties[CurrentSystem].cy,(double)BoxProperties[CurrentSystem].cz);
       exit(0);
    }

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      SetupKVectors();
  }
}

void UpdateCellVelocity(void)
{
  REAL IntPressure,ExtPressure;
  REAL LnVolumeKineticPressurePart,LnVolumePressurePart;
  REAL UKineticTranslation,UKineticTranslationFramework,UKineticTranslationAdsorbates;
  REAL UKineticTranslationCations,UKineticRotation;

  UKineticTranslationFramework=2.0*GetTranslationKineticEnergyFramework();
  UKineticTranslationAdsorbates=2.0*GetTranslationKineticEnergyAdsorbates();
  UKineticTranslationCations=2.0*GetTranslationKineticEnergyCations();
  UKineticRotation=2.0*GetRotationalKineticEnergy();

  UKineticTranslation=UKineticTranslationFramework+UKineticTranslationAdsorbates+
         UKineticTranslationCations;

  IntPressure=Trace3x3Matrix(&ConfigurationalStressTensor[CurrentSystem])/3.0;
  ExtPressure=therm_baro_stats.ExternalPressure[CurrentSystem][0];

  LnVolumePressurePart=3.0*(IntPressure-ExtPressure)*Volume[CurrentSystem];
  LnVolumeKineticPressurePart=(UKineticTranslation)*(1.0+3.0/DegreesOfFreedomTranslation[CurrentSystem]);

  LnVolumeVelocity[CurrentSystem]+=0.5*(LnVolumeKineticPressurePart+LnVolumePressurePart)*DeltaT/LnVolumeMass[CurrentSystem];
}


void NoseHooverNPT(void)
{
  int i,j,k;
  int M,nc,nyosh;
  REAL AA,scale,UKineticBaro;

  M=therm_baro_stats.ThermostatChainLength;
  nc=therm_baro_stats.NumberOfRespaSteps;
  nyosh=therm_baro_stats.NumberOfYoshidaSuzukiSteps;

  scale=1.0;

  UKineticBaro=LnVolumeMass[CurrentSystem]*SQR(LnVolumeVelocity[CurrentSystem]);
  BarostatForce[CurrentSystem][0]=(UKineticBaro-
                                   BarostatDegreesOfFreedom[CurrentSystem][0])/BarostatMass[CurrentSystem][0];

  for(i=0;i<nc;i++)
  {
    for(j=0;j<nyosh;j++)
    {
      BarostatVelocity[CurrentSystem][M-1]+=BarostatForce[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);

      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*BarostatVelocity[CurrentSystem][M-k-1]);
        BarostatVelocity[CurrentSystem][M-k-2]=BarostatVelocity[CurrentSystem][M-k-2]*SQR(AA)+
                                               BarostatForce[CurrentSystem][M-k-2]*AA*w[j]*DeltaT/(4.0*nc);
      }

      for(k=0;k<M;k++)
        BarostatPosition[CurrentSystem][k]+=BarostatVelocity[CurrentSystem][k]*w[j]*DeltaT/(2.0*nc);

      AA=exp(-(w[j]*DeltaT/(2.0*nc))*BarostatVelocity[CurrentSystem][0]);
      scale*=AA;
      BarostatForce[CurrentSystem][0]=(SQR(scale)*UKineticBaro-
                                       BarostatDegreesOfFreedom[CurrentSystem][0])/BarostatMass[CurrentSystem][0];
      for(k=0;k<M-1;k++)
      {

        AA=exp(-(w[j]*DeltaT/(8.0*nc))*BarostatVelocity[CurrentSystem][k+1]);
        BarostatVelocity[CurrentSystem][k]=BarostatVelocity[CurrentSystem][k]*SQR(AA)+
                               BarostatForce[CurrentSystem][k]*AA*(w[j]*DeltaT/(4.0*nc));

        BarostatForce[CurrentSystem][k+1]=(BarostatMass[CurrentSystem][k]*SQR(BarostatVelocity[CurrentSystem][k])-
                                           BarostatDegreesOfFreedom[CurrentSystem][k+1])/BarostatMass[CurrentSystem][k+1];

      }

      BarostatVelocity[CurrentSystem][M-1]+=BarostatForce[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);
    }
  }

  LnVolumeVelocity[CurrentSystem]*=scale;
}


REAL GetTranslationKineticEnergyFramework(void)
{
  int i,f1;
  REAL Mass,TranslationEnergy;
  INT_VECTOR3 Fixed;

  TranslationEnergy=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
      Fixed=Framework[CurrentSystem].Atoms[f1][i].Fixed;
      if(!Fixed.x) TranslationEnergy+=0.5*Mass*SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.x);
      if(!Fixed.y) TranslationEnergy+=0.5*Mass*SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.y);
      if(!Fixed.z) TranslationEnergy+=0.5*Mass*SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.z);
    }
  }
  return TranslationEnergy;
}

REAL GetTranslationKineticEnergyAdsorbates(void)
{
  int i,j,l;
  int Type,A;
  REAL Mass,TranslationEnergy;
  INT_VECTOR3 Fixed;

  TranslationEnergy=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Groups[l].Mass;
        Fixed=Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x);
        if(!Fixed.y) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y);
        if(!Fixed.z) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z);
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];

          Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
          Fixed=Adsorbates[CurrentSystem][i].Atoms[A].Fixed;
          if(!Fixed.x) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x);
          if(!Fixed.y) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y);
          if(!Fixed.z) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z);
        }
      }
    }
  }
  return TranslationEnergy;

}

REAL GetTranslationKineticEnergyCations(void)
{
  int i,j,l;
  int Type,A;
  REAL Mass,TranslationEnergy;
  INT_VECTOR3 Fixed;

  TranslationEnergy=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Groups[l].Mass;
        Fixed=Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x);
        if(!Fixed.y) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y);
        if(!Fixed.z) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z);
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];

          Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[A].Type].Mass;
          Fixed=Cations[CurrentSystem][i].Atoms[A].Fixed;
          if(!Fixed.x) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Atoms[A].Velocity.x);
          if(!Fixed.y) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Atoms[A].Velocity.y);
          if(!Fixed.z) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Atoms[A].Velocity.z);
        }
      }
    }
  }
  return TranslationEnergy;
}


REAL GetTranslationKineticEnergy(void)
{
  int i,j,l;
  int Type,A,f1;
  REAL Mass,TranslationEnergy;
  INT_VECTOR3 Fixed;

  TranslationEnergy=0.0;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
        Fixed=Framework[CurrentSystem].Atoms[f1][i].Fixed;

        if(!Fixed.x) TranslationEnergy+=0.5*Mass*SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.x);
        if(!Fixed.y) TranslationEnergy+=0.5*Mass*SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.y);
        if(!Fixed.z) TranslationEnergy+=0.5*Mass*SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.z);
      }
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Groups[l].Mass;
        Fixed=Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x);
        if(!Fixed.y) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x);
        if(!Fixed.z) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x);
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];

          Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
          Fixed=Adsorbates[CurrentSystem][i].Atoms[A].Fixed;
          if(!Fixed.x) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x);
          if(!Fixed.y) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y);
          if(!Fixed.z) TranslationEnergy+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z);
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
        Mass=Components[Type].Groups[l].Mass;
        Fixed=Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x);
        if(!Fixed.y) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y);
        if(!Fixed.z) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z);
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];

          Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[A].Type].Mass;
          Fixed=Cations[CurrentSystem][i].Atoms[A].Fixed;
          if(!Fixed.x) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Atoms[A].Velocity.x);
          if(!Fixed.y) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Atoms[A].Velocity.y);
          if(!Fixed.z) TranslationEnergy+=0.5*Mass*SQR(Cations[CurrentSystem][i].Atoms[A].Velocity.z);
        }
      }
    }
  }

  return TranslationEnergy;
}

REAL GetRotationalKineticEnergy(void)
{
  int i,l;
  int Type;
  QUATERNION p,q;
  VECTOR I;
  REAL RotationalEnergy;

  RotationalEnergy=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        I=Components[Type].Groups[l].InverseInertiaVector;
        p=Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum;
        q=Adsorbates[CurrentSystem][i].Groups[l].Quaternion;
        RotationalEnergy+=SQR(-p.r*q.i+p.i*q.r+p.j*q.k-p.k*q.j)*I.x/8.0;
        RotationalEnergy+=SQR(-p.r*q.j-p.i*q.k+p.j*q.r+p.k*q.i)*I.y/8.0;
        RotationalEnergy+=SQR(-p.r*q.k+p.i*q.j-p.j*q.i+p.k*q.r)*I.z/8.0;
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
        I=Components[Type].Groups[l].InverseInertiaVector;
        p=Cations[CurrentSystem][i].Groups[l].QuaternionMomentum;
        q=Cations[CurrentSystem][i].Groups[l].Quaternion;
        RotationalEnergy+=SQR(-p.r*q.i+p.i*q.r+p.j*q.k-p.k*q.j)*I.x/8.0;
        RotationalEnergy+=SQR(-p.r*q.j-p.i*q.k+p.j*q.r+p.k*q.i)*I.y/8.0;
        RotationalEnergy+=SQR(-p.r*q.k+p.i*q.j-p.j*q.i+p.k*q.r)*I.z/8.0;
      }
    }
  }
  return RotationalEnergy;
}

REAL GetCellTemperature(void)
{
  REAL UKineticCell,T;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
    case MuPT:
    case NPH:
      UKineticCell=0.5*LnVolumeMass[CurrentSystem]*SQR(LnVolumeVelocity[CurrentSystem]);
      T=2.0*UKineticCell/K_B;
      break;
    case NPTPR:
    case MuPTPR:
    case NPHPR:
      UKineticCell=0.5*CellMass[CurrentSystem]*
         (SQR(CellVelocity[CurrentSystem].ax)+SQR(CellVelocity[CurrentSystem].ay)+SQR(CellVelocity[CurrentSystem].az)+
          SQR(CellVelocity[CurrentSystem].bx)+SQR(CellVelocity[CurrentSystem].by)+SQR(CellVelocity[CurrentSystem].bz)+
          SQR(CellVelocity[CurrentSystem].cx)+SQR(CellVelocity[CurrentSystem].cy)+SQR(CellVelocity[CurrentSystem].cz));
      switch(NPTPRCellType[CurrentSystem])
      {
        default:
        case REGULAR:
          // 6 degrees of freedom for a fully flexible box
          T=2.0*UKineticCell/(9.0*K_B);
          break;
        case REGULAR_UPPER_TRIANGLE:
          // 6 degrees of freedom for a fully flexible box
          T=2.0*UKineticCell/(6.0*K_B);
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          // 4 degrees of freedom for a fully flexible box
          T=2.0*UKineticCell/(4.0*K_B);
          break;
        case ANISOTROPIC:
          // 3 degrees of freedom for an anisotropic box
          T=2.0*UKineticCell/(3.0*K_B);
          break;
        case ISOTROPIC:
          // 1 degree of freedom for an isotropic box
          T=2.0*UKineticCell/(3.0*K_B);
          break;
      }
      break;
    default:
      T=0.0;
      break;
  }

  return T;
}

REAL GetCellKineticEnergy(void)
{
  REAL U;

  U=0.5*LnVolumeMass[CurrentSystem]*SQR(LnVolumeVelocity[CurrentSystem]);
  return U;
}

void NoseHooverNVTFramework(void)
{
  int i,j,k,f1;
  int M,nc,nyosh;
  REAL AA;
  REAL ScaleTranslationFramework;
  REAL UKineticTranslationFramework;
  INT_VECTOR3 Fixed;

  M=therm_baro_stats.ThermostatChainLength;
  nc=therm_baro_stats.NumberOfRespaSteps;
  nyosh=therm_baro_stats.NumberOfYoshidaSuzukiSteps;

  ScaleTranslationFramework=1.0;

  UKineticTranslationFramework=2.0*GetTranslationKineticEnergyFramework();

  ThermostatForceTranslationFramework[CurrentSystem][0]=
        (UKineticTranslationFramework-ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][0])/
                ThermostatMassTranslationFramework[CurrentSystem][0];

  for(i=0;i<nc;i++)
    for(j=0;j<nyosh;j++)
    {
      ThermostatVelocityTranslationFramework[CurrentSystem][M-1]+=ThermostatForceTranslationFramework[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);
      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationFramework[CurrentSystem][M-k-1]);
        ThermostatVelocityTranslationFramework[CurrentSystem][M-k-2]=ThermostatVelocityTranslationFramework[CurrentSystem][M-k-2]*SQR(AA)+
                                     ThermostatForceTranslationFramework[CurrentSystem][M-k-2]*AA*w[j]*DeltaT/(4.0*nc);
      }

      AA=exp(-(w[j]*DeltaT/(2.0*nc))*ThermostatVelocityTranslationFramework[CurrentSystem][0]);
      ScaleTranslationFramework*=AA;
      ThermostatForceTranslationFramework[CurrentSystem][0]=(SQR(ScaleTranslationFramework)*UKineticTranslationFramework-
            ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][0])/
                                        ThermostatMassTranslationFramework[CurrentSystem][0];

      for(k=0;k<M;k++)
        ThermostatPositionTranslationFramework[CurrentSystem][k]+=ThermostatVelocityTranslationFramework[CurrentSystem][k]*w[j]*DeltaT/(2.0*nc);

      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationFramework[CurrentSystem][k+1]);
        ThermostatVelocityTranslationFramework[CurrentSystem][k]=ThermostatVelocityTranslationFramework[CurrentSystem][k]*SQR(AA)+
                               ThermostatForceTranslationFramework[CurrentSystem][k]*AA*(w[j]*DeltaT/(4.0*nc));

        ThermostatForceTranslationFramework[CurrentSystem][k+1]=
              (ThermostatMassTranslationFramework[CurrentSystem][k]*SQR(ThermostatVelocityTranslationFramework[CurrentSystem][k])-
               ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][k+1])/ThermostatMassTranslationFramework[CurrentSystem][k+1];
      }

      ThermostatVelocityTranslationFramework[CurrentSystem][M-1]+=ThermostatForceTranslationFramework[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);
    }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Fixed=Framework[CurrentSystem].Atoms[f1][i].Fixed;

      if(!Fixed.x) Framework[CurrentSystem].Atoms[f1][i].Velocity.x*=ScaleTranslationFramework;
      if(!Fixed.y) Framework[CurrentSystem].Atoms[f1][i].Velocity.y*=ScaleTranslationFramework;
      if(!Fixed.z) Framework[CurrentSystem].Atoms[f1][i].Velocity.z*=ScaleTranslationFramework;
    }
  }

  UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationFramework[CurrentSystem][0])*
      ThermostatMassTranslationFramework[CurrentSystem][0]+
      DegreesOfFreedomFramework[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
      ThermostatPositionTranslationFramework[CurrentSystem][0];
  for(i=1;i<M;i++)
    UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationFramework[CurrentSystem][i])*
           ThermostatMassTranslationFramework[CurrentSystem][i]+
           K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionTranslationFramework[CurrentSystem][i];
}

void NoseHooverNVTAdsorbates(void)
{
  int i,j,k,l;
  int A,M,nc,nyosh,Type;
  REAL AA;
  REAL ScaleTranslationAdsorbates;
  REAL UKineticTranslationAdsorbates;

  M=therm_baro_stats.ThermostatChainLength;
  nc=therm_baro_stats.NumberOfRespaSteps;
  nyosh=therm_baro_stats.NumberOfYoshidaSuzukiSteps;

  ScaleTranslationAdsorbates=1.0;

  UKineticTranslationAdsorbates=2.0*GetTranslationKineticEnergyAdsorbates();

  ThermostatForceTranslationAdsorbates[CurrentSystem][0]=
        (UKineticTranslationAdsorbates-ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][0])/
                ThermostatMassTranslationAdsorbates[CurrentSystem][0];

  for(i=0;i<nc;i++)
  {
    for(j=0;j<nyosh;j++)
    {
      ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-1]+=ThermostatForceTranslationAdsorbates[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);

      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-k-1]);
        ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-k-2]=ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-k-2]*SQR(AA)+
                                       ThermostatForceTranslationAdsorbates[CurrentSystem][M-k-2]*AA*w[j]*DeltaT/(4.0*nc);
      }

      AA=exp(-(w[j]*DeltaT/(2.0*nc))*ThermostatVelocityTranslationAdsorbates[CurrentSystem][0]);
      ScaleTranslationAdsorbates*=AA;
      ThermostatForceTranslationAdsorbates[CurrentSystem][0]=(SQR(ScaleTranslationAdsorbates)*UKineticTranslationAdsorbates-
              ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][0])/
                                          ThermostatMassTranslationAdsorbates[CurrentSystem][0];

      for(k=0;k<M;k++)
        ThermostatPositionTranslationAdsorbates[CurrentSystem][k]+=ThermostatVelocityTranslationAdsorbates[CurrentSystem][k]*w[j]*DeltaT/(2.0*nc);

      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationAdsorbates[CurrentSystem][k+1]);
        ThermostatVelocityTranslationAdsorbates[CurrentSystem][k]=ThermostatVelocityTranslationAdsorbates[CurrentSystem][k]*SQR(AA)+
                                 ThermostatForceTranslationAdsorbates[CurrentSystem][k]*AA*(w[j]*DeltaT/(4.0*nc));

        ThermostatForceTranslationAdsorbates[CurrentSystem][k+1]=
                (ThermostatMassTranslationAdsorbates[CurrentSystem][k]*SQR(ThermostatVelocityTranslationAdsorbates[CurrentSystem][k])-
                ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][k+1])/ThermostatMassTranslationAdsorbates[CurrentSystem][k+1];
      }
      ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-1]+=ThermostatForceTranslationAdsorbates[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x*=ScaleTranslationAdsorbates;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y*=ScaleTranslationAdsorbates;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z*=ScaleTranslationAdsorbates;
      }
      else
      {
        if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x*=ScaleTranslationAdsorbates;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y*=ScaleTranslationAdsorbates;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z*=ScaleTranslationAdsorbates;
        }
      }
    }
  }

  UNoseHoover[CurrentSystem]=0.5*SQR(ThermostatVelocityTranslationAdsorbates[CurrentSystem][0])*
        ThermostatMassTranslationAdsorbates[CurrentSystem][0]+
        DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
        ThermostatPositionTranslationAdsorbates[CurrentSystem][0];

  for(i=1;i<M;i++)
    UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationAdsorbates[CurrentSystem][i])*
            ThermostatMassTranslationAdsorbates[CurrentSystem][i]+
            K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionTranslationAdsorbates[CurrentSystem][i];
}

void NoseHooverNVTCations(void)
{
  int i,j,k,l;
  int A,M,nc,nyosh,Type;
  REAL AA;
  REAL ScaleTranslationCations;
  REAL UKineticTranslationCations;

  M=therm_baro_stats.ThermostatChainLength;
  nc=therm_baro_stats.NumberOfRespaSteps;
  nyosh=therm_baro_stats.NumberOfYoshidaSuzukiSteps;

  ScaleTranslationCations=1.0;

  UKineticTranslationCations=2.0*GetTranslationKineticEnergyCations();

  ThermostatForceTranslationCations[CurrentSystem][0]=
        (UKineticTranslationCations-ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][0])/
                ThermostatMassTranslationCations[CurrentSystem][0];

  for(i=0;i<nc;i++)
    for(j=0;j<nyosh;j++)
    {
      ThermostatVelocityTranslationCations[CurrentSystem][M-1]+=ThermostatForceTranslationCations[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);

      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationCations[CurrentSystem][M-k-1]);
        ThermostatVelocityTranslationCations[CurrentSystem][M-k-2]=ThermostatVelocityTranslationCations[CurrentSystem][M-k-2]*SQR(AA)+
                                     ThermostatForceTranslationCations[CurrentSystem][M-k-2]*AA*w[j]*DeltaT/(4.0*nc);
      }

      AA=exp(-(w[j]*DeltaT/(2.0*nc))*ThermostatVelocityTranslationCations[CurrentSystem][0]);
      ScaleTranslationCations*=AA;
      ThermostatForceTranslationCations[CurrentSystem][0]=(SQR(ScaleTranslationCations)*UKineticTranslationCations-
            ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][0])/
                                        ThermostatMassTranslationCations[CurrentSystem][0];

      for(k=0;k<M;k++)
        ThermostatPositionTranslationCations[CurrentSystem][k]+=ThermostatVelocityTranslationCations[CurrentSystem][k]*w[j]*DeltaT/(2.0*nc);

      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationCations[CurrentSystem][k+1]);
        ThermostatVelocityTranslationCations[CurrentSystem][k]=ThermostatVelocityTranslationCations[CurrentSystem][k]*SQR(AA)+
                               ThermostatForceTranslationCations[CurrentSystem][k]*AA*(w[j]*DeltaT/(4.0*nc));

        ThermostatForceTranslationCations[CurrentSystem][k+1]=
              (ThermostatMassTranslationCations[CurrentSystem][k]*SQR(ThermostatVelocityTranslationCations[CurrentSystem][k])-
               ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][k+1])/ThermostatMassTranslationCations[CurrentSystem][k+1];
      }
      ThermostatVelocityTranslationCations[CurrentSystem][M-1]+=ThermostatForceTranslationCations[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);
    }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x*=ScaleTranslationCations;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y*=ScaleTranslationCations;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z*=ScaleTranslationCations;
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Cations[CurrentSystem][i].Atoms[A].Velocity.x*=ScaleTranslationCations;
          Cations[CurrentSystem][i].Atoms[A].Velocity.y*=ScaleTranslationCations;
          Cations[CurrentSystem][i].Atoms[A].Velocity.z*=ScaleTranslationCations;
        }
     }
    }
  }

  UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationCations[CurrentSystem][0])*
            ThermostatMassTranslationCations[CurrentSystem][0]+
            DegreesOfFreedomTranslationalCations[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
            ThermostatPositionTranslationCations[CurrentSystem][0];

  for(i=1;i<M;i++)
    UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityTranslationCations[CurrentSystem][i])*
            ThermostatMassTranslationCations[CurrentSystem][i]+
            K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionTranslationCations[CurrentSystem][i];
}

void NoseHooverNVTRotation(void)
{
  int i,j,k,l;
  int M,nc,nyosh,Type;
  REAL AA;
  REAL ScaleRotation;
  REAL UKineticRotation;

  M=therm_baro_stats.ThermostatChainLength;
  nc=therm_baro_stats.NumberOfRespaSteps;
  nyosh=therm_baro_stats.NumberOfYoshidaSuzukiSteps;

  ScaleRotation=1.0;

  UKineticRotation=2.0*GetRotationalKineticEnergy();

  ThermostatForceRotation[CurrentSystem][0]=(UKineticRotation-ThermostatDegreesOfFreedomRotation[CurrentSystem][0])/
                ThermostatMassRotation[CurrentSystem][0];

  for(i=0;i<nc;i++)
    for(j=0;j<nyosh;j++)
    {
      ThermostatVelocityRotation[CurrentSystem][M-1]+=ThermostatForceRotation[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);

      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityRotation[CurrentSystem][M-k-1]);
        ThermostatVelocityRotation[CurrentSystem][M-k-2]=ThermostatVelocityRotation[CurrentSystem][M-k-2]*SQR(AA)+
                                     ThermostatForceRotation[CurrentSystem][M-k-2]*AA*w[j]*DeltaT/(4.0*nc);
      }

      AA=exp(-(w[j]*DeltaT/(2.0*nc))*ThermostatVelocityRotation[CurrentSystem][0]);
      ScaleRotation*=AA;
      ThermostatForceRotation[CurrentSystem][0]=(SQR(ScaleRotation)*UKineticRotation-ThermostatDegreesOfFreedomRotation[CurrentSystem][0])/
                                        ThermostatMassRotation[CurrentSystem][0];

      for(k=0;k<M;k++)
        ThermostatPositionRotation[CurrentSystem][k]+=ThermostatVelocityRotation[CurrentSystem][k]*w[j]*DeltaT/(2.0*nc);

      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityRotation[CurrentSystem][k+1]);
        ThermostatVelocityRotation[CurrentSystem][k]=ThermostatVelocityRotation[CurrentSystem][k]*SQR(AA)+
                               ThermostatForceRotation[CurrentSystem][k]*AA*(w[j]*DeltaT/(4.0*nc));

        ThermostatForceRotation[CurrentSystem][k+1]=(ThermostatMassRotation[CurrentSystem][k]*SQR(ThermostatVelocityRotation[CurrentSystem][k])-
                                      ThermostatDegreesOfFreedomRotation[CurrentSystem][k+1])/ThermostatMassRotation[CurrentSystem][k+1];
      }
      ThermostatVelocityRotation[CurrentSystem][M-1]+=ThermostatForceRotation[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);
    }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.r*=ScaleRotation;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.i*=ScaleRotation;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.j*=ScaleRotation;
        Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum.k*=ScaleRotation;
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
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.r*=ScaleRotation;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.i*=ScaleRotation;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.j*=ScaleRotation;
        Cations[CurrentSystem][i].Groups[l].QuaternionMomentum.k*=ScaleRotation;
      }
    }
  }

  UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityRotation[CurrentSystem][0])*ThermostatMassRotation[CurrentSystem][0]+
                              DegreesOfFreedomRotation[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*
                              ThermostatPositionRotation[CurrentSystem][0];
  for(i=1;i<M;i++)
    UNoseHoover[CurrentSystem]+=0.5*SQR(ThermostatVelocityRotation[CurrentSystem][i])*ThermostatMassRotation[CurrentSystem][i]+
            K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*ThermostatPositionRotation[CurrentSystem][i];
}

// NoseHoover baro/thermostat (NPT-Ensemble, fully flexible cell) by Glenn J. Martyna et al.
// Ref: Molecular Physics 1996, Vol. 87, No. 5, 1117-1157

void NoseHooverNPTPR(void)
{
  int i,j,k,l,m,f1,A;
  int M,N,nc,nyosh,Type;
  REAL AA,trace,avg;
  REAL_MATRIX3x3 vtemps;
  REAL_MATRIX3x3 eigenvectors;
  VECTOR eigenvalues,tempv;
  VECTOR vel,vexpdt;
  REAL UKineticTranslation,UKineticCell;
  REAL UKineticTranslationFramework;
  REAL UKineticTranslationAdsorbates;
  REAL UKineticTranslationCations;
  INT_VECTOR3 Fixed;


  M=therm_baro_stats.ThermostatChainLength;
  N=therm_baro_stats.BarostatChainLength;
  nc=therm_baro_stats.NumberOfRespaSteps;
  nyosh=therm_baro_stats.NumberOfYoshidaSuzukiSteps;

  // compute force on the first thermostat
  UKineticTranslation=2.0*GetTranslationKineticEnergy();
  UKineticCell=CellMass[CurrentSystem]*
      (SQR(CellVelocity[CurrentSystem].ax)+SQR(CellVelocity[CurrentSystem].ay)+SQR(CellVelocity[CurrentSystem].az)+
       SQR(CellVelocity[CurrentSystem].bx)+SQR(CellVelocity[CurrentSystem].by)+SQR(CellVelocity[CurrentSystem].bz)+
       SQR(CellVelocity[CurrentSystem].cx)+SQR(CellVelocity[CurrentSystem].cy)+SQR(CellVelocity[CurrentSystem].cz));


  // thermostat force on the framework
  UKineticTranslationFramework=2.0*GetTranslationKineticEnergyFramework();
  ThermostatForceTranslationFramework[CurrentSystem][0]=
        (UKineticTranslationFramework-ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][0])/
                ThermostatMassTranslationFramework[CurrentSystem][0];

  // thermostat force on the adsorbates
  UKineticTranslationAdsorbates=2.0*GetTranslationKineticEnergyAdsorbates();
  ThermostatForceTranslationAdsorbates[CurrentSystem][0]=
        (UKineticTranslationAdsorbates-ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][0])/
                ThermostatMassTranslationAdsorbates[CurrentSystem][0];

  // thermostat force on the cations
  UKineticTranslationCations=2.0*GetTranslationKineticEnergyCations();
  ThermostatForceTranslationCations[CurrentSystem][0]=
        (UKineticTranslationCations-ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][0])/
                ThermostatMassTranslationCations[CurrentSystem][0];


  BarostatForce[CurrentSystem][0]=(UKineticCell-BarostatDegreesOfFreedom[CurrentSystem][0])/BarostatMass[CurrentSystem][0];

  // compute force on the box
  CellForce=GetKineticStressTensor();
  trace=UKineticTranslation/(REAL)DegreesOfFreedomTranslation[CurrentSystem];
  CellForce.ax+=trace;
  CellForce.by+=trace;
  CellForce.cz+=trace;

  CellForce.ax+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].ax-therm_baro_stats.ExternalPressure[CurrentSystem][0]);
  CellForce.ay+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].ay);
  CellForce.az+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].az);
  CellForce.bx+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].bx);
  CellForce.by+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].by-therm_baro_stats.ExternalPressure[CurrentSystem][0]);
  CellForce.bz+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].bz);
  CellForce.cx+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cx);
  CellForce.cy+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cy);
  CellForce.cz+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cz-therm_baro_stats.ExternalPressure[CurrentSystem][0]);

  CellForce.ax/=CellMass[CurrentSystem]; CellForce.bx/=CellMass[CurrentSystem]; CellForce.cx/=CellMass[CurrentSystem];
  CellForce.ay/=CellMass[CurrentSystem]; CellForce.by/=CellMass[CurrentSystem]; CellForce.cy/=CellMass[CurrentSystem];
  CellForce.az/=CellMass[CurrentSystem]; CellForce.bz/=CellMass[CurrentSystem]; CellForce.cz/=CellMass[CurrentSystem];

  switch(NPTPRCellType[CurrentSystem])
  {
    case REGULAR_UPPER_TRIANGLE:
      CellForce.ay=CellForce.bz=CellForce.az=0.0;
      break;
    case ISOTROPIC:
      CellForce.ay=CellForce.bz=CellForce.az=0.0;
      CellForce.bx=CellForce.cy=CellForce.cx=0.0;
      avg=(CellForce.ax+CellForce.by+CellForce.cz)/3.0;
      CellForce.ax=CellForce.by=CellForce.cz=avg;
      break;
    case ANISOTROPIC:
      CellForce.ay=CellForce.bz=CellForce.az=0.0;
      CellForce.bx=CellForce.cy=CellForce.cx=0.0;
      break;
    case MONOCLINIC:
      CellForce.bx=CellForce.cy=CellForce.ay=CellForce.bz=0.0;
      avg=0.5*(CellForce.az+CellForce.cx);
      CellForce.az=CellForce.cx=avg;
      break;
    case MONOCLINIC_UPPER_TRIANGLE:
      CellForce.bx=CellForce.cy=CellForce.ay=CellForce.az=CellForce.bz=0.0;
      break;
    case REGULAR:
    default:
      avg=0.5*(CellForce.ay+CellForce.bx);
      CellForce.ay=CellForce.bx=avg;
      avg=0.5*(CellForce.az+CellForce.cx);
      CellForce.az=CellForce.cx=avg;
      avg=0.5*(CellForce.bz+CellForce.cy);
      CellForce.bz=CellForce.cy=avg;
      break;
  }

  // start the multiple time-step procedure
  for(i=0;i<nc;i++)
    for(j=0;j<nyosh;j++)
    {
      // update the framework thermostat velocities
      ThermostatVelocityTranslationFramework[CurrentSystem][M-1]+=ThermostatForceTranslationFramework[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);
      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationFramework[CurrentSystem][M-k-1]);
        ThermostatVelocityTranslationFramework[CurrentSystem][M-k-2]=ThermostatVelocityTranslationFramework[CurrentSystem][M-k-2]*SQR(AA)+
                                     ThermostatForceTranslationFramework[CurrentSystem][M-k-2]*AA*w[j]*DeltaT/(4.0*nc);
      }
      // update the adsorbate thermostat velocities
      ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-1]+=ThermostatForceTranslationAdsorbates[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);
      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-k-1]);
        ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-k-2]=ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-k-2]*SQR(AA)+
                                     ThermostatForceTranslationAdsorbates[CurrentSystem][M-k-2]*AA*w[j]*DeltaT/(4.0*nc);
      }
      // update the cation thermostat velocities
      ThermostatVelocityTranslationCations[CurrentSystem][M-1]+=ThermostatForceTranslationCations[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);
      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationCations[CurrentSystem][M-k-1]);
        ThermostatVelocityTranslationCations[CurrentSystem][M-k-2]=ThermostatVelocityTranslationCations[CurrentSystem][M-k-2]*SQR(AA)+
                                     ThermostatForceTranslationCations[CurrentSystem][M-k-2]*AA*w[j]*DeltaT/(4.0*nc);
      }



      // update the barostat velocities
      BarostatVelocity[CurrentSystem][N-1]+=BarostatForce[CurrentSystem][N-1]*w[j]*DeltaT/(4.0*nc);
      for(k=0;k<N-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*BarostatVelocity[CurrentSystem][N-k-1]);
        BarostatVelocity[CurrentSystem][N-k-2]=BarostatVelocity[CurrentSystem][N-k-2]*SQR(AA)+
                                   BarostatForce[CurrentSystem][N-k-2]*AA*w[j]*DeltaT/(4.0*nc);
      }

      // update the box velocities
      AA=exp(-w[j]*BarostatVelocity[CurrentSystem][0]*DeltaT/(8.0*(REAL)nc));
      CellVelocity[CurrentSystem].ax=CellVelocity[CurrentSystem].ax*SQR(AA)+CellForce.ax*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].ay=CellVelocity[CurrentSystem].ay*SQR(AA)+CellForce.ay*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].az=CellVelocity[CurrentSystem].az*SQR(AA)+CellForce.az*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].bx=CellVelocity[CurrentSystem].bx*SQR(AA)+CellForce.bx*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].by=CellVelocity[CurrentSystem].by*SQR(AA)+CellForce.by*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].bz=CellVelocity[CurrentSystem].bz*SQR(AA)+CellForce.bz*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].cx=CellVelocity[CurrentSystem].cx*SQR(AA)+CellForce.cx*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].cy=CellVelocity[CurrentSystem].cy*SQR(AA)+CellForce.cy*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].cz=CellVelocity[CurrentSystem].cz*SQR(AA)+CellForce.cz*AA*w[j]*DeltaT/(4.0*(REAL)nc);

      // update the particle velocities with NHC
      AA=exp(-ThermostatVelocityTranslationFramework[CurrentSystem][0]*w[j]*DeltaT/(4.0*nc));
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          Fixed=Framework[CurrentSystem].Atoms[f1][k].Fixed;
          if(!Fixed.x) Framework[CurrentSystem].Atoms[f1][k].Velocity.x*=AA;
          if(!Fixed.y) Framework[CurrentSystem].Atoms[f1][k].Velocity.y*=AA;
          if(!Fixed.z) Framework[CurrentSystem].Atoms[f1][k].Velocity.z*=AA;
        }
      }

      AA=exp(-ThermostatVelocityTranslationAdsorbates[CurrentSystem][0]*w[j]*DeltaT/(4.0*nc));
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        Type=Adsorbates[CurrentSystem][k].Type;
        for(l=0;l<Components[Type].NumberOfGroups;l++)
        {
          if(Components[Type].Groups[l].Rigid)
          {
            Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity.x*=AA;
            Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity.y*=AA;
            Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity.z*=AA;
          }
          else
          {
            if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
              for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
              {
                A=Components[Type].Groups[l].Atoms[m];
                Adsorbates[CurrentSystem][k].Atoms[A].Velocity.x*=AA;
                Adsorbates[CurrentSystem][k].Atoms[A].Velocity.y*=AA;
                Adsorbates[CurrentSystem][k].Atoms[A].Velocity.z*=AA;
              }
          }
        }
      }

      AA=exp(-ThermostatVelocityTranslationCations[CurrentSystem][0]*w[j]*DeltaT/(4.0*nc));
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        Type=Cations[CurrentSystem][k].Type;
        for(l=0;l<Components[Type].NumberOfGroups;l++)
        {
          if(Components[Type].Groups[l].Rigid)
          {
            Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity.x*=AA;
            Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity.y*=AA;
            Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity.z*=AA;
          }
          else
          {
            if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
              for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
              {
                A=Components[Type].Groups[l].Atoms[m];
                Cations[CurrentSystem][k].Atoms[A].Velocity.x*=AA;
                Cations[CurrentSystem][k].Atoms[A].Velocity.y*=AA;
                Cations[CurrentSystem][k].Atoms[A].Velocity.z*=AA;
              }
          }
        }
      }

      // update the thermostat positions
      for(k=0;k<M;k++)
      {
        ThermostatPositionTranslationFramework[CurrentSystem][k]+=ThermostatVelocityTranslationFramework[CurrentSystem][k]*w[j]*DeltaT/(2.0*nc);
        ThermostatPositionTranslationAdsorbates[CurrentSystem][k]+=ThermostatVelocityTranslationAdsorbates[CurrentSystem][k]*w[j]*DeltaT/(2.0*nc);
        ThermostatPositionTranslationCations[CurrentSystem][k]+=ThermostatVelocityTranslationCations[CurrentSystem][k]*w[j]*DeltaT/(2.0*nc);
      }

      // update the barostat positions
      for(k=0;k<N;k++)
        BarostatPosition[CurrentSystem][k]+=BarostatVelocity[CurrentSystem][k]*w[j]*DeltaT/(2.0*nc);

      // update the particle velocities
      vtemps=CellVelocity[CurrentSystem];
      trace=(CellVelocity[CurrentSystem].ax+CellVelocity[CurrentSystem].by+CellVelocity[CurrentSystem].cz)/(REAL)DegreesOfFreedomTranslation[CurrentSystem];
      vtemps.ax+=trace; vtemps.by+=trace; vtemps.bz+=trace;

      switch(NPTPRCellType[CurrentSystem])
      {
        case REGULAR_UPPER_TRIANGLE:
        case MONOCLINIC_UPPER_TRIANGLE:
          vexpdt.x=exp(-vtemps.ax*w[j]*DeltaT/(2.0*nc));
          vexpdt.y=exp(-vtemps.by*w[j]*DeltaT/(2.0*nc));
          vexpdt.z=exp(-vtemps.cz*w[j]*DeltaT/(2.0*nc));

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              vel=Framework[CurrentSystem].Atoms[f1][k].Velocity;
              vel.x=vel.x*vexpdt.x
                    -CellVelocity[CurrentSystem].cx*vel.z*(vexpdt.x+vexpdt.z)*w[j]*DeltaT/(4.0*nc)
                    -CellVelocity[CurrentSystem].bx*vel.y*(vexpdt.x+vexpdt.y)*w[j]*DeltaT/(4.0*nc)
                    +CellVelocity[CurrentSystem].cy*vel.z*CellVelocity[CurrentSystem].bx*(0.5*(vexpdt.x+vexpdt.z)+vexpdt.y)*SQR(w[j]*DeltaT/(4.0*nc));
              vel.y=vel.y*vexpdt.y-CellVelocity[CurrentSystem].cy*vel.z*(vexpdt.z+vexpdt.y)*w[j]*DeltaT/(4.0*nc);
              vel.z=vel.z*vexpdt.z;
              Framework[CurrentSystem].Atoms[f1][k].Velocity=vel;
            }
          }

          for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
          {
            Type=Adsorbates[CurrentSystem][k].Type;
            for(l=0;l<Components[Type].NumberOfGroups;l++)
            {
              if(Components[Type].Groups[l].Rigid)
              {
                vel=Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity;
                vel.x=vel.x*vexpdt.x-CellVelocity[CurrentSystem].cx*vel.z*(vexpdt.x+vexpdt.z)*w[j]*DeltaT/(4.0*nc)
                      -CellVelocity[CurrentSystem].bx*vel.y*(vexpdt.x+vexpdt.y)*w[j]*DeltaT/(4.0*nc)
                      +CellVelocity[CurrentSystem].cy*vel.z*CellVelocity[CurrentSystem].bx*(0.5*(vexpdt.x+vexpdt.z)+vexpdt.y)*SQR(w[j]*DeltaT/(4.0*nc));
                vel.y=vel.y*vexpdt.y-CellVelocity[CurrentSystem].cy*vel.z*(vexpdt.z+vexpdt.y)*w[j]*DeltaT/(4.0*nc);
                vel.z=vel.z*vexpdt.z;
                Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity=vel;
              }
              else
              {
                if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
                  for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
                  {
                    A=Components[Type].Groups[l].Atoms[m];
                    vel=Adsorbates[CurrentSystem][k].Atoms[A].Velocity;
                    vel.x=vel.x*vexpdt.x-CellVelocity[CurrentSystem].cx*vel.z*(vexpdt.x+vexpdt.z)*w[j]*DeltaT/(4.0*nc)
                          -CellVelocity[CurrentSystem].bx*vel.y*(vexpdt.x+vexpdt.y)*w[j]*DeltaT/(4.0*nc)
                          +CellVelocity[CurrentSystem].cy*vel.z*CellVelocity[CurrentSystem].bx*(0.5*(vexpdt.x+vexpdt.z)+vexpdt.y)*SQR(w[j]*DeltaT/(4.0*nc));
                    vel.y=vel.y*vexpdt.y-CellVelocity[CurrentSystem].cy*vel.z*(vexpdt.z+vexpdt.y)*w[j]*DeltaT/(4.0*nc);
                    vel.z=vel.z*vexpdt.z;
                    Adsorbates[CurrentSystem][k].Atoms[A].Velocity=vel;
                  }
              }
            }
          }

         for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
          {
            Type=Cations[CurrentSystem][k].Type;
            for(l=0;l<Components[Type].NumberOfGroups;l++)
            {
              if(Components[Type].Groups[l].Rigid)
              {
                vel=Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity;
                vel.x=vel.x*vexpdt.x-CellVelocity[CurrentSystem].cx*vel.z*(vexpdt.x+vexpdt.z)*w[j]*DeltaT/(4.0*nc)
                      -CellVelocity[CurrentSystem].bx*vel.y*(vexpdt.x+vexpdt.y)*w[j]*DeltaT/(4.0*nc)
                      +CellVelocity[CurrentSystem].cy*vel.z*CellVelocity[CurrentSystem].bx*(0.5*(vexpdt.x+vexpdt.z)+vexpdt.y)*SQR(w[j]*DeltaT/(4.0*nc));
                vel.y=vel.y*vexpdt.y-CellVelocity[CurrentSystem].cy*vel.z*(vexpdt.z+vexpdt.y)*w[j]*DeltaT/(4.0*nc);
                vel.z=vel.z*vexpdt.z;
                Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity=vel;
              }
              else
              {
                if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
                  for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
                  {
                    A=Components[Type].Groups[l].Atoms[m];
                    vel=Cations[CurrentSystem][k].Atoms[A].Velocity;
                    vel=Adsorbates[CurrentSystem][k].Atoms[A].Velocity;
                    vel.x=vel.x*vexpdt.x-CellVelocity[CurrentSystem].cx*vel.z*(vexpdt.x+vexpdt.z)*w[j]*DeltaT/(4.0*nc)
                          -CellVelocity[CurrentSystem].bx*vel.y*(vexpdt.x+vexpdt.y)*w[j]*DeltaT/(4.0*nc)
                          +CellVelocity[CurrentSystem].cy*vel.z*CellVelocity[CurrentSystem].bx*(0.5*(vexpdt.x+vexpdt.z)+vexpdt.y)*SQR(w[j]*DeltaT/(4.0*nc));
                    vel.y=vel.y*vexpdt.y-CellVelocity[CurrentSystem].cy*vel.z*(vexpdt.z+vexpdt.y)*w[j]*DeltaT/(4.0*nc);
                    vel.z=vel.z*vexpdt.z;
                    Cations[CurrentSystem][k].Atoms[A].Velocity=vel;
                  }
              }
            }
          }
          break;
        default:
        case REGULAR:
        case ISOTROPIC:
        case ANISOTROPIC:
        case MONOCLINIC:
          EigenSystem3x3(vtemps,&eigenvectors,&eigenvalues);

          vexpdt.x=exp(-eigenvalues.x*w[j]*DeltaT/(2.0*nc));
          vexpdt.y=exp(-eigenvalues.y*w[j]*DeltaT/(2.0*nc));
          vexpdt.z=exp(-eigenvalues.z*w[j]*DeltaT/(2.0*nc));

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              vel=Framework[CurrentSystem].Atoms[f1][k].Velocity;
              tempv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
              tempv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
              tempv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
              tempv.x*=vexpdt.x;
              tempv.y*=vexpdt.y;
              tempv.z*=vexpdt.z;
              vel.x=tempv.x*eigenvectors.ax+tempv.y*eigenvectors.ay+tempv.z*eigenvectors.az;
              vel.y=tempv.x*eigenvectors.bx+tempv.y*eigenvectors.by+tempv.z*eigenvectors.bz;
              vel.z=tempv.x*eigenvectors.cx+tempv.y*eigenvectors.cy+tempv.z*eigenvectors.cz;
              Framework[CurrentSystem].Atoms[f1][k].Velocity=vel;
            }
          }

          for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
          {
            Type=Adsorbates[CurrentSystem][k].Type;
            for(l=0;l<Components[Type].NumberOfGroups;l++)
            {
              if(Components[Type].Groups[l].Rigid)
              {
                vel=Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity;
                tempv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
                tempv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
                tempv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
                tempv.x*=vexpdt.x;
                tempv.y*=vexpdt.y;
                tempv.z*=vexpdt.z;
                vel.x=tempv.x*eigenvectors.ax+tempv.y*eigenvectors.ay+tempv.z*eigenvectors.az;
                vel.y=tempv.x*eigenvectors.bx+tempv.y*eigenvectors.by+tempv.z*eigenvectors.bz;
                vel.z=tempv.x*eigenvectors.cx+tempv.y*eigenvectors.cy+tempv.z*eigenvectors.cz;
                Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity=vel;
              }
              else
              {
                if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
                  for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
                  {
                    A=Components[Type].Groups[l].Atoms[m];
                    vel=Adsorbates[CurrentSystem][k].Atoms[A].Velocity;
                    tempv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
                    tempv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
                    tempv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
                    tempv.x*=vexpdt.x;
                    tempv.y*=vexpdt.y;
                    tempv.z*=vexpdt.z;
                    vel.x=tempv.x*eigenvectors.ax+tempv.y*eigenvectors.ay+tempv.z*eigenvectors.az;
                    vel.y=tempv.x*eigenvectors.bx+tempv.y*eigenvectors.by+tempv.z*eigenvectors.bz;
                    vel.z=tempv.x*eigenvectors.cx+tempv.y*eigenvectors.cy+tempv.z*eigenvectors.cz;
                    Adsorbates[CurrentSystem][k].Atoms[A].Velocity=vel;
                  }
              }
            }
          }

          for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
          {
            Type=Cations[CurrentSystem][k].Type;
            for(l=0;l<Components[Type].NumberOfGroups;l++)
            {
              if(Components[Type].Groups[l].Rigid)
              {
                vel=Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity;
                tempv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
                tempv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
                tempv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
                tempv.x*=vexpdt.x;
                tempv.y*=vexpdt.y;
                tempv.z*=vexpdt.z;
                vel.x=tempv.x*eigenvectors.ax+tempv.y*eigenvectors.ay+tempv.z*eigenvectors.az;
                vel.y=tempv.x*eigenvectors.bx+tempv.y*eigenvectors.by+tempv.z*eigenvectors.bz;
                vel.z=tempv.x*eigenvectors.cx+tempv.y*eigenvectors.cy+tempv.z*eigenvectors.cz;
                Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity=vel;
              }
              else
              {
                if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
                  for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
                  {
                    A=Components[Type].Groups[l].Atoms[m];
                    vel=Cations[CurrentSystem][k].Atoms[A].Velocity;
                    tempv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
                    tempv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
                    tempv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
                    tempv.x*=vexpdt.x;
                    tempv.y*=vexpdt.y;
                    tempv.z*=vexpdt.z;
                    vel.x=tempv.x*eigenvectors.ax+tempv.y*eigenvectors.ay+tempv.z*eigenvectors.az;
                    vel.y=tempv.x*eigenvectors.bx+tempv.y*eigenvectors.by+tempv.z*eigenvectors.bz;
                    vel.z=tempv.x*eigenvectors.cx+tempv.y*eigenvectors.cy+tempv.z*eigenvectors.cz;
                    Cations[CurrentSystem][k].Atoms[A].Velocity=vel;
                  }
              }
            }
          }
          break;
      }

      // update the particle velocities with NHC
      AA=exp(-ThermostatVelocityTranslationFramework[CurrentSystem][0]*w[j]*DeltaT/(4.0*nc));
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          Framework[CurrentSystem].Atoms[f1][k].Velocity.x*=AA;
          Framework[CurrentSystem].Atoms[f1][k].Velocity.y*=AA;
          Framework[CurrentSystem].Atoms[f1][k].Velocity.z*=AA;
        }
      }

      AA=exp(-ThermostatVelocityTranslationAdsorbates[CurrentSystem][0]*w[j]*DeltaT/(4.0*nc));
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        Type=Adsorbates[CurrentSystem][k].Type;
        for(l=0;l<Components[Type].NumberOfGroups;l++)
        {
          if(Components[Type].Groups[l].Rigid)
          {
            Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity.x*=AA;
            Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity.y*=AA;
            Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity.z*=AA;
          }
          else
          {
            if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
              for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
              {
                A=Components[Type].Groups[l].Atoms[m];
                Adsorbates[CurrentSystem][k].Atoms[A].Velocity.x*=AA;
                Adsorbates[CurrentSystem][k].Atoms[A].Velocity.y*=AA;
                Adsorbates[CurrentSystem][k].Atoms[A].Velocity.z*=AA;
              }
          }
        }
      }

      AA=exp(-ThermostatVelocityTranslationCations[CurrentSystem][0]*w[j]*DeltaT/(4.0*nc));
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        Type=Cations[CurrentSystem][k].Type;
        for(l=0;l<Components[Type].NumberOfGroups;l++)
        {
          if(Components[Type].Groups[l].Rigid)
          {
            Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity.x*=AA;
            Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity.y*=AA;
            Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity.z*=AA;
          }
          else
          {
            if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
              for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
              {
                A=Components[Type].Groups[l].Atoms[m];
                Cations[CurrentSystem][k].Atoms[A].Velocity.x*=AA;
                Cations[CurrentSystem][k].Atoms[A].Velocity.y*=AA;
                Cations[CurrentSystem][k].Atoms[A].Velocity.z*=AA;
              }
          }
        }
      }

      // update the force on the box
      UKineticTranslation=2.0*GetTranslationKineticEnergy();
      CellForce=GetKineticStressTensor();
      trace=UKineticTranslation/(REAL)DegreesOfFreedomTranslation[CurrentSystem];
      CellForce.ax+=trace;
      CellForce.by+=trace;
      CellForce.cz+=trace;

      CellForce.ax+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].ax-therm_baro_stats.ExternalPressure[CurrentSystem][0]);
      CellForce.ay+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].ay);
      CellForce.az+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].az);
      CellForce.bx+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].bx);
      CellForce.by+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].by-therm_baro_stats.ExternalPressure[CurrentSystem][0]);
      CellForce.bz+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].bz);
      CellForce.cx+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cx);
      CellForce.cy+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cy);
      CellForce.cz+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cz-therm_baro_stats.ExternalPressure[CurrentSystem][0]);

      CellForce.ax/=CellMass[CurrentSystem]; CellForce.bx/=CellMass[CurrentSystem]; CellForce.cx/=CellMass[CurrentSystem];
      CellForce.ay/=CellMass[CurrentSystem]; CellForce.by/=CellMass[CurrentSystem]; CellForce.cy/=CellMass[CurrentSystem];
      CellForce.az/=CellMass[CurrentSystem]; CellForce.bz/=CellMass[CurrentSystem]; CellForce.cz/=CellMass[CurrentSystem];

      switch(NPTPRCellType[CurrentSystem])
      {
        case REGULAR_UPPER_TRIANGLE:
          CellForce.ay=CellForce.bz=CellForce.az=0.0;
          break;
        case ISOTROPIC:
          CellForce.ay=CellForce.bz=CellForce.az=0.0;
          CellForce.bx=CellForce.cy=CellForce.cx=0.0;
          avg=(CellForce.ax+CellForce.by+CellForce.cz)/3.0;
          CellForce.ax=CellForce.by=CellForce.cz=avg;
          break;
        case ANISOTROPIC:
          CellForce.ay=CellForce.bz=CellForce.az=0.0;
          CellForce.bx=CellForce.cy=CellForce.cx=0.0;
          break;
        case MONOCLINIC:
          CellForce.bx=CellForce.cy=CellForce.ay=CellForce.bz=0.0;
          avg=0.5*(CellForce.az+CellForce.cx);
          CellForce.az=CellForce.cx=avg;
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          CellForce.bx=CellForce.cy=CellForce.ay=CellForce.az=CellForce.bz=0.0;
          break;
        case REGULAR:
        default:
          avg=0.5*(CellForce.ay+CellForce.bx);
          CellForce.ay=CellForce.bx=avg;
          avg=0.5*(CellForce.az+CellForce.cx);
          CellForce.az=CellForce.cx=avg;
          avg=0.5*(CellForce.bz+CellForce.cy);
          CellForce.bz=CellForce.cy=avg;
          break;
      }

      // update the box velocities
      AA=exp(-w[j]*BarostatVelocity[CurrentSystem][0]*DeltaT/(8.0*(REAL)nc));
      CellVelocity[CurrentSystem].ax=CellVelocity[CurrentSystem].ax*SQR(AA)+CellForce.ax*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].ay=CellVelocity[CurrentSystem].ay*SQR(AA)+CellForce.ay*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].az=CellVelocity[CurrentSystem].az*SQR(AA)+CellForce.az*AA*w[j]*DeltaT/(4.0*(REAL)nc);

      CellVelocity[CurrentSystem].bx=CellVelocity[CurrentSystem].bx*SQR(AA)+CellForce.bx*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].by=CellVelocity[CurrentSystem].by*SQR(AA)+CellForce.by*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].bz=CellVelocity[CurrentSystem].bz*SQR(AA)+CellForce.bz*AA*w[j]*DeltaT/(4.0*(REAL)nc);

      CellVelocity[CurrentSystem].cx=CellVelocity[CurrentSystem].cx*SQR(AA)+CellForce.cx*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].cy=CellVelocity[CurrentSystem].cy*SQR(AA)+CellForce.cy*AA*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].cz=CellVelocity[CurrentSystem].cz*SQR(AA)+CellForce.cz*AA*w[j]*DeltaT/(4.0*(REAL)nc);

      UKineticCell=CellMass[CurrentSystem]*
          (SQR(CellVelocity[CurrentSystem].ax)+SQR(CellVelocity[CurrentSystem].ay)+SQR(CellVelocity[CurrentSystem].az)+
           SQR(CellVelocity[CurrentSystem].bx)+SQR(CellVelocity[CurrentSystem].by)+SQR(CellVelocity[CurrentSystem].bz)+
           SQR(CellVelocity[CurrentSystem].cx)+SQR(CellVelocity[CurrentSystem].cy)+SQR(CellVelocity[CurrentSystem].cz));


      // update the force on the first framework thermostat
      UKineticTranslationFramework=2.0*GetTranslationKineticEnergyFramework();
      ThermostatForceTranslationFramework[CurrentSystem][0]=
            (UKineticTranslationFramework-ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][0])/
                   ThermostatMassTranslationFramework[CurrentSystem][0];

      // update the force on the first adsorbate thermostat
      UKineticTranslationAdsorbates=2.0*GetTranslationKineticEnergyAdsorbates();
      ThermostatForceTranslationAdsorbates[CurrentSystem][0]=
             (UKineticTranslationAdsorbates-ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][0])/
                     ThermostatMassTranslationAdsorbates[CurrentSystem][0];

      // update the force on the first cation thermostat
      UKineticTranslationCations=2.0*GetTranslationKineticEnergyCations();
      ThermostatForceTranslationCations[CurrentSystem][0]=
            (UKineticTranslationCations-ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][0])/
                    ThermostatMassTranslationCations[CurrentSystem][0];

      BarostatForce[CurrentSystem][0]=(UKineticCell-BarostatDegreesOfFreedom[CurrentSystem][0])/BarostatMass[CurrentSystem][0];


      // update the framework thermostat velocities
      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationFramework[CurrentSystem][k+1]);
        ThermostatVelocityTranslationFramework[CurrentSystem][k]=ThermostatVelocityTranslationFramework[CurrentSystem][k]*SQR(AA)+
                               ThermostatForceTranslationFramework[CurrentSystem][k]*AA*(w[j]*DeltaT/(4.0*nc));
        ThermostatForceTranslationFramework[CurrentSystem][k+1]=(ThermostatMassTranslationFramework[CurrentSystem][k]*SQR(ThermostatVelocityTranslationFramework[CurrentSystem][k])-
                                      ThermostatDegreesOfFreedomTranslationFramework[CurrentSystem][k+1])/ThermostatMassTranslationFramework[CurrentSystem][k+1];
      }
      ThermostatVelocityTranslationFramework[CurrentSystem][M-1]+=ThermostatForceTranslationFramework[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);

      // update the adsorbate thermostat velocities
      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationAdsorbates[CurrentSystem][k+1]);
        ThermostatVelocityTranslationAdsorbates[CurrentSystem][k]=ThermostatVelocityTranslationAdsorbates[CurrentSystem][k]*SQR(AA)+
                               ThermostatForceTranslationAdsorbates[CurrentSystem][k]*AA*(w[j]*DeltaT/(4.0*nc));
        ThermostatForceTranslationAdsorbates[CurrentSystem][k+1]=(ThermostatMassTranslationAdsorbates[CurrentSystem][k]*SQR(ThermostatVelocityTranslationAdsorbates[CurrentSystem][k])-
                                      ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][k+1])/ThermostatMassTranslationAdsorbates[CurrentSystem][k+1];
      }
      ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-1]+=ThermostatForceTranslationAdsorbates[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);

      // update the cations thermostat velocities
      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationCations[CurrentSystem][k+1]);
        ThermostatVelocityTranslationCations[CurrentSystem][k]=ThermostatVelocityTranslationCations[CurrentSystem][k]*SQR(AA)+
                               ThermostatForceTranslationCations[CurrentSystem][k]*AA*(w[j]*DeltaT/(4.0*nc));
        ThermostatForceTranslationCations[CurrentSystem][k+1]=(ThermostatMassTranslationCations[CurrentSystem][k]*SQR(ThermostatVelocityTranslationCations[CurrentSystem][k])-
                                      ThermostatDegreesOfFreedomTranslationCations[CurrentSystem][k+1])/ThermostatMassTranslationCations[CurrentSystem][k+1];
      }
      ThermostatVelocityTranslationCations[CurrentSystem][M-1]+=ThermostatForceTranslationCations[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);

      // update the barostat velocities
      for(k=0;k<N-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*BarostatVelocity[CurrentSystem][k+1]);
        BarostatVelocity[CurrentSystem][k]=BarostatVelocity[CurrentSystem][k]*SQR(AA)+
                             BarostatForce[CurrentSystem][k]*AA*(w[j]*DeltaT/(4.0*nc));
        BarostatForce[CurrentSystem][k+1]=(BarostatMass[CurrentSystem][k]*SQR(BarostatVelocity[CurrentSystem][k])-
                                    BarostatDegreesOfFreedom[CurrentSystem][k+1])/BarostatMass[CurrentSystem][k+1];
      }
      BarostatVelocity[CurrentSystem][N-1]+=BarostatForce[CurrentSystem][N-1]*w[j]*DeltaT/(4.0*nc);
    }
}

void UpdatePositionsVelocitiesNPTPR(void)
{
  int i,l,A,Type,f1,j;
  REAL E2,E4,E6,E8,aa,arg2,det;
  VECTOR aa2,bb;
  REAL_MATRIX3x3 eigenvectors,vtemps;
  VECTOR eigenvalues;
  VECTOR pos,vel,u,uv;
  int nres,ires;
  REAL e2,e4,e6,e8;
  REAL dt2,dt24,vg48,dt_res,dt2_res,dt24_res,dt;
  REAL gamma11,gamma22,gamma33;
  REAL sinh11,sinh22,sinh33;
  REAL ax,ay,az,bx,by,bz;
  REAL vg15,vg59,vg19,sindt;
  REAL_MATRIX3x3 OldBox;


  E2=1.0/6.0;
  E4=E2/20.0;
  E6=E4/42.0;
  E8=E6/72.0;

  e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
  e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);


  switch(NPTPRCellType[CurrentSystem])
  {
    case REGULAR_UPPER_TRIANGLE:
    case MONOCLINIC_UPPER_TRIANGLE:
      dt=DeltaT;
      dt2     = dt/2.0;
      dt24    = dt2*dt2;
      vg48    = CellVelocity[CurrentSystem].bx*CellVelocity[CurrentSystem].cy;

      nres     = 5;
      dt_res   = dt  /((double)nres);
      dt2_res  = dt2 /((double)nres);
      dt24_res = dt24/((double)(nres*nres));

      gamma11  = (CellVelocity[CurrentSystem].ax*dt_res)/2.0;
      arg2     = gamma11*gamma11;
      sinh11   = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

      gamma22  = (CellVelocity[CurrentSystem].by*dt_res)/2.0;
      arg2     = gamma22*gamma22;
      sinh22   = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

      gamma33  = (CellVelocity[CurrentSystem].cz*dt_res)/2.0;
      arg2     = gamma33*gamma33;
      sinh33   = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

      ax       = dt_res*sinh11*exp(CellVelocity[CurrentSystem].ax*dt2_res);
      ay       = dt_res*sinh22*exp(CellVelocity[CurrentSystem].by*dt2_res);
      az       = dt_res*sinh33*exp(CellVelocity[CurrentSystem].cz*dt2_res);

      bx       = exp(CellVelocity[CurrentSystem].ax*dt_res);
      by       = exp(CellVelocity[CurrentSystem].by*dt_res);
      bz       = exp(CellVelocity[CurrentSystem].cz*dt_res);

      // Evolve the particles from 0 to dt
      for(ires=1;ires<=nres;ires++)
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
          {
            pos=Framework[CurrentSystem].Atoms[f1][i].Position;
            vel=Framework[CurrentSystem].Atoms[f1][i].Velocity;

            //  Evolve the postitions from 0 to dt with exp(iL_corr*dt/2)
            pos.x+=(CellVelocity[CurrentSystem].bx*pos.y*dt2_res+0.5*vg48*dt24_res*pos.z
                    +CellVelocity[CurrentSystem].cx*pos.z*dt2_res);
            pos.y+=(CellVelocity[CurrentSystem].cy*pos.z*dt2_res);

            // Evolve the postitions from 0 to dt with exp(iL_reference*dt)
            pos.x=ax*vel.x+bx*pos.x;
            pos.y=ay*vel.y+by*pos.y;
            pos.z=az*vel.z+bz*pos.z;

            // Evolve the postitions from dt2 to dt with exp(iL_correction*dt2)
            pos.x+=(CellVelocity[CurrentSystem].bx*pos.y*dt2_res+0.5*vg48*dt24_res*pos.z
                   +CellVelocity[CurrentSystem].cx*pos.z*dt2_res);
            pos.y+=(CellVelocity[CurrentSystem].cy*pos.z*dt2_res);
            Framework[CurrentSystem].Atoms[f1][i].Position=pos;
          }
        }

        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          Type=Adsorbates[CurrentSystem][i].Type;
          for(l=0;l<Components[Type].NumberOfGroups;l++)
          {
            if(Components[Type].Groups[l].Rigid)
            {
              vel=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
              pos=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition;

              //  Evolve the postitions from 0 to dt with exp(iL_corr*dt/2)
              pos.x+=(CellVelocity[CurrentSystem].bx*pos.y*dt2_res+0.5*vg48*dt24_res*pos.z
                      +CellVelocity[CurrentSystem].cx*pos.z*dt2_res);
              pos.y+=(CellVelocity[CurrentSystem].cy*pos.z*dt2_res);

              // Evolve the postitions from 0 to dt with exp(iL_reference*dt)
              pos.x=ax*vel.x+bx*pos.x;
              pos.y=ay*vel.y+by*pos.y;
              pos.z=az*vel.z+bz*pos.z;

              // Evolve the postitions from dt2 to dt with exp(iL_correction*dt2)
              pos.x+=(CellVelocity[CurrentSystem].bx*pos.y*dt2_res+0.5*vg48*dt24_res*pos.z
                     +CellVelocity[CurrentSystem].cx*pos.z*dt2_res);
              pos.y+=(CellVelocity[CurrentSystem].cy*pos.z*dt2_res);
              Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition=pos;
            }
            else
            {
              if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
              {
                for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
                {
                  A=Components[Type].Groups[l].Atoms[j];
                  vel=Adsorbates[CurrentSystem][i].Atoms[A].Velocity;
                  pos=Adsorbates[CurrentSystem][i].Atoms[A].Position;

                  //  Evolve the postitions from 0 to dt with exp(iL_corr*dt/2)
                  pos.x+=(CellVelocity[CurrentSystem].bx*pos.y*dt2_res+0.5*vg48*dt24_res*pos.z
                          +CellVelocity[CurrentSystem].cx*pos.z*dt2_res);
                  pos.y+=(CellVelocity[CurrentSystem].cy*pos.z*dt2_res);

                  // Evolve the postitions from 0 to dt with exp(iL_reference*dt)
                  pos.x=ax*vel.x+bx*pos.x;
                  pos.y=ay*vel.y+by*pos.y;
                  pos.z=az*vel.z+bz*pos.z;

                  // Evolve the postitions from dt2 to dt with exp(iL_correction*dt2)
                  pos.x+=(CellVelocity[CurrentSystem].bx*pos.y*dt2_res+0.5*vg48*dt24_res*pos.z
                          +CellVelocity[CurrentSystem].cx*pos.z*dt2_res);
                  pos.y+=(CellVelocity[CurrentSystem].cy*pos.z*dt2_res);
                  Adsorbates[CurrentSystem][i].Atoms[A].Position=pos;
                }
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
              vel=Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
              pos=Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition;

              //  Evolve the postitions from 0 to dt with exp(iL_corr*dt/2)
              pos.x+=(CellVelocity[CurrentSystem].bx*pos.y*dt2_res+0.5*vg48*dt24_res*pos.z
                      +CellVelocity[CurrentSystem].cx*pos.z*dt2_res);
              pos.y+=(CellVelocity[CurrentSystem].cy*pos.z*dt2_res);

              // Evolve the postitions from 0 to dt with exp(iL_reference*dt)
              pos.x=ax*vel.x+bx*pos.x;
              pos.y=ay*vel.y+by*pos.y;
              pos.z=az*vel.z+bz*pos.z;

              // Evolve the postitions from dt2 to dt with exp(iL_correction*dt2)
              pos.x+=(CellVelocity[CurrentSystem].bx*pos.y*dt2_res+0.5*vg48*dt24_res*pos.z
                     +CellVelocity[CurrentSystem].cx*pos.z*dt2_res);
              pos.y+=(CellVelocity[CurrentSystem].cy*pos.z*dt2_res);
              Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition=pos;
            }
            else
            {
              if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
              {
                for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
                {
                  A=Components[Type].Groups[l].Atoms[j];
                  vel=Cations[CurrentSystem][i].Atoms[A].Velocity;
                  pos=Cations[CurrentSystem][i].Atoms[A].Position;

                  //  Evolve the postitions from 0 to dt with exp(iL_corr*dt/2)
                  pos.x+=(CellVelocity[CurrentSystem].bx*pos.y*dt2_res+0.5*vg48*dt24_res*pos.z
                          +CellVelocity[CurrentSystem].cx*pos.z*dt2_res);
                  pos.y+=(CellVelocity[CurrentSystem].cy*pos.z*dt2_res);

                  // Evolve the postitions from 0 to dt with exp(iL_reference*dt)
                  pos.x=ax*vel.x+bx*pos.x;
                  pos.y=ay*vel.y+by*pos.y;
                  pos.z=az*vel.z+bz*pos.z;

                  // Evolve the postitions from dt2 to dt with exp(iL_correction*dt2)
                  pos.x+=(CellVelocity[CurrentSystem].bx*pos.y*dt2_res+0.5*vg48*dt24_res*pos.z
                          +CellVelocity[CurrentSystem].cx*pos.z*dt2_res);
                  pos.y+=(CellVelocity[CurrentSystem].cy*pos.z*dt2_res);
                  Cations[CurrentSystem][i].Atoms[A].Position=pos;
                }
              }
            }
          }
        }
      }

      // Evolve the cell matrix
      nres     = 5;
      dt_res   = dt/((double)nres);
      dt24_res = dt24/((double)(nres*nres));

      for(ires=1;ires<=nres;ires++)
      {
        // Evolve the matrix from 0 to dt2 with exp(iL_correction*dt2)

        Box[CurrentSystem].cx+=(CellVelocity[CurrentSystem].bx*Box[CurrentSystem].cy*dt_res);

        // Save the old h-matrix
        OldBox=Box[CurrentSystem];

        // Evolve the matrix from 0 to dt with exp(iL_reference*dt)
        Box[CurrentSystem].ax=OldBox.ax*exp(CellVelocity[CurrentSystem].ax*dt_res);
        Box[CurrentSystem].by=OldBox.by*exp(CellVelocity[CurrentSystem].by*dt_res);
        Box[CurrentSystem].cz=OldBox.cz*exp(CellVelocity[CurrentSystem].cz*dt_res);

        vg15=(CellVelocity[CurrentSystem].ax-CellVelocity[CurrentSystem].by)*dt_res/2.0;
        arg2=vg15*vg15;
        sindt=(((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

        Box[CurrentSystem].bx=OldBox.bx*exp(CellVelocity[CurrentSystem].ax*dt_res)
                +OldBox.by*CellVelocity[CurrentSystem].bx*dt_res
                *exp((CellVelocity[CurrentSystem].ax+CellVelocity[CurrentSystem].by)*dt_res/2.0)*sindt;

        vg59=(CellVelocity[CurrentSystem].by-CellVelocity[CurrentSystem].cz)*dt_res/2.0;
        arg2=vg59*vg59;
        sindt=(((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

        Box[CurrentSystem].cy=OldBox.cy*exp(CellVelocity[CurrentSystem].by*dt_res)
                +OldBox.cz*CellVelocity[CurrentSystem].cy*dt_res
                *exp((CellVelocity[CurrentSystem].by+CellVelocity[CurrentSystem].cz)*dt_res/2.0)*sindt;

        vg19=(CellVelocity[CurrentSystem].ax-CellVelocity[CurrentSystem].cz)*dt_res;
        arg2=vg19*vg19;
        sindt=(((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

        Box[CurrentSystem].cx=OldBox.cx*exp(CellVelocity[CurrentSystem].ax*dt_res)
                +OldBox.cz*CellVelocity[CurrentSystem].cx*dt_res
                *exp((CellVelocity[CurrentSystem].ax+CellVelocity[CurrentSystem].cz)*dt_res)*sindt;

        // Evolve the matrix from dt2 to dt with exp(iL_correction*dt)
        Box[CurrentSystem].cx += CellVelocity[CurrentSystem].bx*Box[CurrentSystem].cy*dt_res;
      }
      break;
    case REGULAR:
    case MONOCLINIC:
    case ISOTROPIC:
    case ANISOTROPIC:
    default:
      EigenSystem3x3(CellVelocity[CurrentSystem],&eigenvectors,&eigenvalues);

      aa=exp(0.5*DeltaT*eigenvalues.x);   aa2.x=SQR(aa);
      arg2=SQR(0.5*DeltaT*eigenvalues.x); bb.x=aa*((((E8*arg2+E6)*arg2+E4)*arg2+E2)*arg2+1.0)*DeltaT;

      aa=exp(0.5*DeltaT*eigenvalues.y);   aa2.y=SQR(aa);
      arg2=SQR(0.5*DeltaT*eigenvalues.y); bb.y=aa*((((E8*arg2+E6)*arg2+E4)*arg2+E2)*arg2+1.0)*DeltaT;

      aa=exp(0.5*DeltaT*eigenvalues.z);   aa2.z=SQR(aa);
      arg2=SQR(0.5*DeltaT*eigenvalues.z); bb.z=aa*((((E8*arg2+E6)*arg2+E4)*arg2+E2)*arg2+1.0)*DeltaT;


      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          vel=Framework[CurrentSystem].Atoms[f1][i].Velocity;
          pos=Framework[CurrentSystem].Atoms[f1][i].Position;
          u.x=eigenvectors.ax*pos.x+eigenvectors.bx*pos.y+eigenvectors.cx*pos.z;
          u.y=eigenvectors.ay*pos.x+eigenvectors.by*pos.y+eigenvectors.cy*pos.z;
          u.z=eigenvectors.az*pos.x+eigenvectors.bz*pos.y+eigenvectors.cz*pos.z;
          uv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
          uv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
          uv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
          u.x=u.x*aa2.x+uv.x*bb.x;
          u.y=u.y*aa2.y+uv.y*bb.y;
          u.z=u.z*aa2.z+uv.z*bb.z;
          pos.x=u.x*eigenvectors.ax+u.y*eigenvectors.ay+u.z*eigenvectors.az;
          pos.y=u.x*eigenvectors.bx+u.y*eigenvectors.by+u.z*eigenvectors.bz;
          pos.z=u.x*eigenvectors.cx+u.y*eigenvectors.cy+u.z*eigenvectors.cz;

          Framework[CurrentSystem].Atoms[f1][i].Position=pos;
        }
      }

      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        Type=Adsorbates[CurrentSystem][i].Type;
        for(l=0;l<Components[Type].NumberOfGroups;l++)
        {
          if(Components[Type].Groups[l].Rigid)
          {
            vel=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
            pos=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition;
            u.x=eigenvectors.ax*pos.x+eigenvectors.bx*pos.y+eigenvectors.cx*pos.z;
            u.y=eigenvectors.ay*pos.x+eigenvectors.by*pos.y+eigenvectors.cy*pos.z;
            u.z=eigenvectors.az*pos.x+eigenvectors.bz*pos.y+eigenvectors.cz*pos.z;
            uv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
            uv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
            uv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
            u.x=u.x*aa2.x+uv.x*bb.x;
            u.y=u.y*aa2.y+uv.y*bb.y;
            u.z=u.z*aa2.z+uv.z*bb.z;
            pos.x=u.x*eigenvectors.ax+u.y*eigenvectors.ay+u.z*eigenvectors.az;
            pos.y=u.x*eigenvectors.bx+u.y*eigenvectors.by+u.z*eigenvectors.bz;
            pos.z=u.x*eigenvectors.cx+u.y*eigenvectors.cy+u.z*eigenvectors.cz;
            Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition=pos;
          }
          else
          {
            if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
              for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
              {
                A=Components[Type].Groups[l].Atoms[j];
                vel=Adsorbates[CurrentSystem][i].Atoms[A].Velocity;
                pos=Adsorbates[CurrentSystem][i].Atoms[A].Position;
                u.x=eigenvectors.ax*pos.x+eigenvectors.bx*pos.y+eigenvectors.cx*pos.z;
                u.y=eigenvectors.ay*pos.x+eigenvectors.by*pos.y+eigenvectors.cy*pos.z;
                u.z=eigenvectors.az*pos.x+eigenvectors.bz*pos.y+eigenvectors.cz*pos.z;
                uv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
                uv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
                uv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
                u.x=u.x*aa2.x+uv.x*bb.x;
                u.y=u.y*aa2.y+uv.y*bb.y;
                u.z=u.z*aa2.z+uv.z*bb.z;
                pos.x=u.x*eigenvectors.ax+u.y*eigenvectors.ay+u.z*eigenvectors.az;
                pos.y=u.x*eigenvectors.bx+u.y*eigenvectors.by+u.z*eigenvectors.bz;
                pos.z=u.x*eigenvectors.cx+u.y*eigenvectors.cy+u.z*eigenvectors.cz;
                Adsorbates[CurrentSystem][i].Atoms[A].Position=pos;
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
            vel=Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
            pos=Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition;
            u.x=eigenvectors.ax*pos.x+eigenvectors.bx*pos.y+eigenvectors.cx*pos.z;
            u.y=eigenvectors.ay*pos.x+eigenvectors.by*pos.y+eigenvectors.cy*pos.z;
            u.z=eigenvectors.az*pos.x+eigenvectors.bz*pos.y+eigenvectors.cz*pos.z;
            uv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
            uv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
            uv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
            u.x=u.x*aa2.x+uv.x*bb.x;
            u.y=u.y*aa2.y+uv.y*bb.y;
            u.z=u.z*aa2.z+uv.z*bb.z;
            pos.x=u.x*eigenvectors.ax+u.y*eigenvectors.ay+u.z*eigenvectors.az;
            pos.y=u.x*eigenvectors.bx+u.y*eigenvectors.by+u.z*eigenvectors.bz;
            pos.z=u.x*eigenvectors.cx+u.y*eigenvectors.cy+u.z*eigenvectors.cz;
            Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition=pos;
          }
          else
          {
            if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
              for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
              {
                A=Components[Type].Groups[l].Atoms[j];
                vel=Cations[CurrentSystem][i].Atoms[A].Velocity;
                pos=Cations[CurrentSystem][i].Atoms[A].Position;
                u.x=eigenvectors.ax*pos.x+eigenvectors.bx*pos.y+eigenvectors.cx*pos.z;
                u.y=eigenvectors.ay*pos.x+eigenvectors.by*pos.y+eigenvectors.cy*pos.z;
                u.z=eigenvectors.az*pos.x+eigenvectors.bz*pos.y+eigenvectors.cz*pos.z;
                uv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
                uv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
                uv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
                u.x=u.x*aa2.x+uv.x*bb.x;
                u.y=u.y*aa2.y+uv.y*bb.y;
                u.z=u.z*aa2.z+uv.z*bb.z;
                pos.x=u.x*eigenvectors.ax+u.y*eigenvectors.ay+u.z*eigenvectors.az;
                pos.y=u.x*eigenvectors.bx+u.y*eigenvectors.by+u.z*eigenvectors.bz;
                pos.z=u.x*eigenvectors.cx+u.y*eigenvectors.cy+u.z*eigenvectors.cz;
                Cations[CurrentSystem][i].Atoms[A].Position=pos;
              }
          }
        }
      }

      // update box
      vtemps.ax=eigenvectors.ax*Box[CurrentSystem].ax+eigenvectors.bx*Box[CurrentSystem].ay+eigenvectors.cx*Box[CurrentSystem].az;
      vtemps.ay=eigenvectors.ax*Box[CurrentSystem].bx+eigenvectors.bx*Box[CurrentSystem].by+eigenvectors.cx*Box[CurrentSystem].bz;
      vtemps.az=eigenvectors.ax*Box[CurrentSystem].cx+eigenvectors.bx*Box[CurrentSystem].cy+eigenvectors.cx*Box[CurrentSystem].cz;

      vtemps.bx=eigenvectors.ay*Box[CurrentSystem].ax+eigenvectors.by*Box[CurrentSystem].ay+eigenvectors.cy*Box[CurrentSystem].az;
      vtemps.by=eigenvectors.ay*Box[CurrentSystem].bx+eigenvectors.by*Box[CurrentSystem].by+eigenvectors.cy*Box[CurrentSystem].bz;
      vtemps.bz=eigenvectors.ay*Box[CurrentSystem].cx+eigenvectors.by*Box[CurrentSystem].cy+eigenvectors.cy*Box[CurrentSystem].cz;

      vtemps.cx=eigenvectors.az*Box[CurrentSystem].ax+eigenvectors.bz*Box[CurrentSystem].ay+eigenvectors.cz*Box[CurrentSystem].az;
      vtemps.cy=eigenvectors.az*Box[CurrentSystem].bx+eigenvectors.bz*Box[CurrentSystem].by+eigenvectors.cz*Box[CurrentSystem].bz;
      vtemps.cz=eigenvectors.az*Box[CurrentSystem].cx+eigenvectors.bz*Box[CurrentSystem].cy+eigenvectors.cz*Box[CurrentSystem].cz;

      vtemps.ax*=aa2.x;
      vtemps.ay*=aa2.x;
      vtemps.az*=aa2.x;

      vtemps.bx*=aa2.y;
      vtemps.by*=aa2.y;
      vtemps.bz*=aa2.y;

      vtemps.cx*=aa2.z;
      vtemps.cy*=aa2.z;
      vtemps.cz*=aa2.z;

      Box[CurrentSystem].ax=eigenvectors.ax*vtemps.ax+eigenvectors.ay*vtemps.bx+eigenvectors.az*vtemps.cx;
      Box[CurrentSystem].bx=eigenvectors.ax*vtemps.ay+eigenvectors.ay*vtemps.by+eigenvectors.az*vtemps.cy;
      Box[CurrentSystem].cx=eigenvectors.ax*vtemps.az+eigenvectors.ay*vtemps.bz+eigenvectors.az*vtemps.cz;

      Box[CurrentSystem].ay=eigenvectors.bx*vtemps.ax+eigenvectors.by*vtemps.bx+eigenvectors.bz*vtemps.cx;
      Box[CurrentSystem].by=eigenvectors.bx*vtemps.ay+eigenvectors.by*vtemps.by+eigenvectors.bz*vtemps.cy;
      Box[CurrentSystem].cy=eigenvectors.bx*vtemps.az+eigenvectors.by*vtemps.bz+eigenvectors.bz*vtemps.cz;

      Box[CurrentSystem].az=eigenvectors.cx*vtemps.ax+eigenvectors.cy*vtemps.bx+eigenvectors.cz*vtemps.cx;
      Box[CurrentSystem].bz=eigenvectors.cx*vtemps.ay+eigenvectors.cy*vtemps.by+eigenvectors.cz*vtemps.cy;
      Box[CurrentSystem].cz=eigenvectors.cx*vtemps.az+eigenvectors.cy*vtemps.bz+eigenvectors.cz*vtemps.cz;
      break;
  }

  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);

  AlphaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bx);
  BetaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].by);
  GammaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bz);

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    SetupKVectors();
}

// isoenthalpic-isobaric ensemble
void NoseHooverNPHPR(void)
{
  int i,j,k,l,m,f1,A;
  int nc,nyosh,Type;
  REAL trace,avg;
  REAL_MATRIX3x3 vtemps;
  REAL_MATRIX3x3 eigenvectors;
  VECTOR eigenvalues,tempv;
  VECTOR vel,vexpdt;
  REAL UKineticTranslation;

  nc=therm_baro_stats.NumberOfRespaSteps;
  nyosh=therm_baro_stats.NumberOfYoshidaSuzukiSteps;

  // compute force on the box
  UKineticTranslation=2.0*GetTranslationKineticEnergy();
  CellForce=GetKineticStressTensor();
  trace=UKineticTranslation/(REAL)DegreesOfFreedomTranslation[CurrentSystem];
  CellForce.ax+=trace;
  CellForce.by+=trace;
  CellForce.cz+=trace;

  CellForce.ax+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].ax-therm_baro_stats.ExternalPressure[CurrentSystem][0]);
  CellForce.ay+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].ay);
  CellForce.az+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].az);
  CellForce.bx+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].bx);
  CellForce.by+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].by-therm_baro_stats.ExternalPressure[CurrentSystem][0]);
  CellForce.bz+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].bz);
  CellForce.cx+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cx);
  CellForce.cy+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cy);
  CellForce.cz+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cz-therm_baro_stats.ExternalPressure[CurrentSystem][0]);

  CellForce.ax/=CellMass[CurrentSystem]; CellForce.bx/=CellMass[CurrentSystem]; CellForce.cx/=CellMass[CurrentSystem];
  CellForce.ay/=CellMass[CurrentSystem]; CellForce.by/=CellMass[CurrentSystem]; CellForce.cy/=CellMass[CurrentSystem];
  CellForce.az/=CellMass[CurrentSystem]; CellForce.bz/=CellMass[CurrentSystem]; CellForce.cz/=CellMass[CurrentSystem];

  switch(NPTPRCellType[CurrentSystem])
  {
    case REGULAR_UPPER_TRIANGLE:
      CellForce.ay=CellForce.bz=CellForce.az=0.0;
      break;
    case ISOTROPIC:
      CellForce.ay=CellForce.bz=CellForce.az=0.0;
      CellForce.bx=CellForce.cy=CellForce.cx=0.0;
      avg=(CellForce.ax+CellForce.by+CellForce.cz)/3.0;
      CellForce.ax=CellForce.by=CellForce.cz=avg;
      break;
    case ANISOTROPIC:
      CellForce.ay=CellForce.bz=CellForce.az=0.0;
      CellForce.bx=CellForce.cy=CellForce.cx=0.0;
      break;
    case MONOCLINIC:
      CellForce.bx=CellForce.cy=CellForce.ay=CellForce.bz=0.0;
      avg=0.5*(CellForce.az+CellForce.cx);
      CellForce.az=CellForce.cx=avg;
      break;
    case MONOCLINIC_UPPER_TRIANGLE:
      CellForce.bx=CellForce.cy=CellForce.ay=CellForce.az=CellForce.bz=0.0;
      break;
    case REGULAR:
    default:
      avg=0.5*(CellForce.ay+CellForce.bx);
      CellForce.ay=CellForce.bx=avg;
      avg=0.5*(CellForce.az+CellForce.cx);
      CellForce.az=CellForce.cx=avg;
      avg=0.5*(CellForce.bz+CellForce.cy);
      CellForce.bz=CellForce.cy=avg;
      break;
  }

  // start the multiple time-step procedure
  for(i=0;i<nc;i++)
    for(j=0;j<nyosh;j++)
    {
      // update the box velocities
      CellVelocity[CurrentSystem].ax+=CellForce.ax*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].ay+=CellForce.ay*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].az+=CellForce.az*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].bx+=CellForce.bx*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].by+=CellForce.by*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].bz+=CellForce.bz*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].cx+=CellForce.cx*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].cy+=CellForce.cy*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].cz+=CellForce.cz*w[j]*DeltaT/(4.0*(REAL)nc);

      // update the particle velocities
      vtemps=CellVelocity[CurrentSystem];
      trace=(CellVelocity[CurrentSystem].ax+CellVelocity[CurrentSystem].by+CellVelocity[CurrentSystem].cz)/(REAL)DegreesOfFreedomTranslation[CurrentSystem];
      vtemps.ax+=trace; vtemps.by+=trace; vtemps.bz+=trace;

      switch(NPTPRCellType[CurrentSystem])
      {
        case REGULAR_UPPER_TRIANGLE:
        case MONOCLINIC_UPPER_TRIANGLE:
          vexpdt.x=exp(-vtemps.ax*w[j]*DeltaT/(2.0*nc));
          vexpdt.y=exp(-vtemps.by*w[j]*DeltaT/(2.0*nc));
          vexpdt.z=exp(-vtemps.cz*w[j]*DeltaT/(2.0*nc));

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              vel=Framework[CurrentSystem].Atoms[f1][k].Velocity;
              vel.x=vel.x*vexpdt.x-CellVelocity[CurrentSystem].cx*vel.z*(vexpdt.x+vexpdt.z)*w[j]*DeltaT/(4.0*nc)
                    -CellVelocity[CurrentSystem].bx*vel.y*(vexpdt.x+vexpdt.y)*w[j]*DeltaT/(4.0*nc)
                    +CellVelocity[CurrentSystem].cy*vel.z*CellVelocity[CurrentSystem].bx*(0.5*(vexpdt.x+vexpdt.z)+vexpdt.y)*SQR(w[j]*DeltaT/(4.0*nc));
              vel.y=vel.y*vexpdt.y-CellVelocity[CurrentSystem].cy*vel.z*(vexpdt.z+vexpdt.y)*w[j]*DeltaT/(4.0*nc);
              vel.z=vel.z*vexpdt.z;
              Framework[CurrentSystem].Atoms[f1][k].Velocity=vel;
            }
          }

          for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
          {
            Type=Adsorbates[CurrentSystem][k].Type;
            for(l=0;l<Components[Type].NumberOfGroups;l++)
            {
              if(Components[Type].Groups[l].Rigid)
              {
                vel=Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity;
                vel.x=vel.x*vexpdt.x-CellVelocity[CurrentSystem].cx*vel.z*(vexpdt.x+vexpdt.z)*w[j]*DeltaT/(4.0*nc)
                      -CellVelocity[CurrentSystem].bx*vel.y*(vexpdt.x+vexpdt.y)*w[j]*DeltaT/(4.0*nc)
                      +CellVelocity[CurrentSystem].cy*vel.z*CellVelocity[CurrentSystem].bx*(0.5*(vexpdt.x+vexpdt.z)+vexpdt.y)*SQR(w[j]*DeltaT/(4.0*nc));
                vel.y=vel.y*vexpdt.y-CellVelocity[CurrentSystem].cy*vel.z*(vexpdt.z+vexpdt.y)*w[j]*DeltaT/(4.0*nc);
                vel.z=vel.z*vexpdt.z;
                Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity=vel;
              }
              else
              {
                if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
                  for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
                  {
                    A=Components[Type].Groups[l].Atoms[m];
                    vel=Adsorbates[CurrentSystem][k].Atoms[A].Velocity;
                    vel.x=vel.x*vexpdt.x-CellVelocity[CurrentSystem].cx*vel.z*(vexpdt.x+vexpdt.z)*w[j]*DeltaT/(4.0*nc)
                          -CellVelocity[CurrentSystem].bx*vel.y*(vexpdt.x+vexpdt.y)*w[j]*DeltaT/(4.0*nc)
                          +CellVelocity[CurrentSystem].cy*vel.z*CellVelocity[CurrentSystem].bx*(0.5*(vexpdt.x+vexpdt.z)+vexpdt.y)*SQR(w[j]*DeltaT/(4.0*nc));
                    vel.y=vel.y*vexpdt.y-CellVelocity[CurrentSystem].cy*vel.z*(vexpdt.z+vexpdt.y)*w[j]*DeltaT/(4.0*nc);
                    vel.z=vel.z*vexpdt.z;
                    Adsorbates[CurrentSystem][k].Atoms[A].Velocity=vel;
                  }
              }
            }
          }

         for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
          {
            Type=Cations[CurrentSystem][k].Type;
            for(l=0;l<Components[Type].NumberOfGroups;l++)
            {
              if(Components[Type].Groups[l].Rigid)
              {
                vel=Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity;
                vel.x=vel.x*vexpdt.x-CellVelocity[CurrentSystem].cx*vel.z*(vexpdt.x+vexpdt.z)*w[j]*DeltaT/(4.0*nc)
                      -CellVelocity[CurrentSystem].bx*vel.y*(vexpdt.x+vexpdt.y)*w[j]*DeltaT/(4.0*nc)
                      +CellVelocity[CurrentSystem].cy*vel.z*CellVelocity[CurrentSystem].bx*(0.5*(vexpdt.x+vexpdt.z)+vexpdt.y)*SQR(w[j]*DeltaT/(4.0*nc));
                vel.y=vel.y*vexpdt.y-CellVelocity[CurrentSystem].cy*vel.z*(vexpdt.z+vexpdt.y)*w[j]*DeltaT/(4.0*nc);
                vel.z=vel.z*vexpdt.z;
                Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity=vel;
              }
              else
              {
                if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
                  for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
                  {
                    A=Components[Type].Groups[l].Atoms[m];
                    vel=Cations[CurrentSystem][k].Atoms[A].Velocity;
                    vel=Adsorbates[CurrentSystem][k].Atoms[A].Velocity;
                    vel.x=vel.x*vexpdt.x-CellVelocity[CurrentSystem].cx*vel.z*(vexpdt.x+vexpdt.z)*w[j]*DeltaT/(4.0*nc)
                          -CellVelocity[CurrentSystem].bx*vel.y*(vexpdt.x+vexpdt.y)*w[j]*DeltaT/(4.0*nc)
                          +CellVelocity[CurrentSystem].cy*vel.z*CellVelocity[CurrentSystem].bx*(0.5*(vexpdt.x+vexpdt.z)+vexpdt.y)*SQR(w[j]*DeltaT/(4.0*nc));
                    vel.y=vel.y*vexpdt.y-CellVelocity[CurrentSystem].cy*vel.z*(vexpdt.z+vexpdt.y)*w[j]*DeltaT/(4.0*nc);
                    vel.z=vel.z*vexpdt.z;
                    Cations[CurrentSystem][k].Atoms[A].Velocity=vel;
                  }
              }
            }
          }
          break;
        default:
        case REGULAR:
        case ISOTROPIC:
        case ANISOTROPIC:
        case MONOCLINIC:
          EigenSystem3x3(vtemps,&eigenvectors,&eigenvalues);

          vexpdt.x=exp(-eigenvalues.x*w[j]*DeltaT/(2.0*nc));
          vexpdt.y=exp(-eigenvalues.y*w[j]*DeltaT/(2.0*nc));
          vexpdt.z=exp(-eigenvalues.z*w[j]*DeltaT/(2.0*nc));

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              vel=Framework[CurrentSystem].Atoms[f1][k].Velocity;
              tempv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
              tempv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
              tempv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
              tempv.x*=vexpdt.x;
              tempv.y*=vexpdt.y;
              tempv.z*=vexpdt.z;
              vel.x=tempv.x*eigenvectors.ax+tempv.y*eigenvectors.ay+tempv.z*eigenvectors.az;
              vel.y=tempv.x*eigenvectors.bx+tempv.y*eigenvectors.by+tempv.z*eigenvectors.bz;
              vel.z=tempv.x*eigenvectors.cx+tempv.y*eigenvectors.cy+tempv.z*eigenvectors.cz;
              Framework[CurrentSystem].Atoms[f1][k].Velocity=vel;
            }
          }

          for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
          {
            Type=Adsorbates[CurrentSystem][k].Type;
            for(l=0;l<Components[Type].NumberOfGroups;l++)
            {
              if(Components[Type].Groups[l].Rigid)
              {
                vel=Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity;
                tempv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
                tempv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
                tempv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
                tempv.x*=vexpdt.x;
                tempv.y*=vexpdt.y;
                tempv.z*=vexpdt.z;
                vel.x=tempv.x*eigenvectors.ax+tempv.y*eigenvectors.ay+tempv.z*eigenvectors.az;
                vel.y=tempv.x*eigenvectors.bx+tempv.y*eigenvectors.by+tempv.z*eigenvectors.bz;
                vel.z=tempv.x*eigenvectors.cx+tempv.y*eigenvectors.cy+tempv.z*eigenvectors.cz;
                Adsorbates[CurrentSystem][k].Groups[l].CenterOfMassVelocity=vel;
              }
              else
              {
                if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
                  for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
                  {
                    A=Components[Type].Groups[l].Atoms[m];
                    vel=Adsorbates[CurrentSystem][k].Atoms[A].Velocity;
                    tempv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
                    tempv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
                    tempv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
                    tempv.x*=vexpdt.x;
                    tempv.y*=vexpdt.y;
                    tempv.z*=vexpdt.z;
                    vel.x=tempv.x*eigenvectors.ax+tempv.y*eigenvectors.ay+tempv.z*eigenvectors.az;
                    vel.y=tempv.x*eigenvectors.bx+tempv.y*eigenvectors.by+tempv.z*eigenvectors.bz;
                    vel.z=tempv.x*eigenvectors.cx+tempv.y*eigenvectors.cy+tempv.z*eigenvectors.cz;
                    Adsorbates[CurrentSystem][k].Atoms[A].Velocity=vel;
                  }
              }
            }
          }

          for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
          {
            Type=Cations[CurrentSystem][k].Type;
            for(l=0;l<Components[Type].NumberOfGroups;l++)
            {
              if(Components[Type].Groups[l].Rigid)
              {
                vel=Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity;
                tempv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
                tempv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
                tempv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
                tempv.x*=vexpdt.x;
                tempv.y*=vexpdt.y;
                tempv.z*=vexpdt.z;
                vel.x=tempv.x*eigenvectors.ax+tempv.y*eigenvectors.ay+tempv.z*eigenvectors.az;
                vel.y=tempv.x*eigenvectors.bx+tempv.y*eigenvectors.by+tempv.z*eigenvectors.bz;
                vel.z=tempv.x*eigenvectors.cx+tempv.y*eigenvectors.cy+tempv.z*eigenvectors.cz;
                Cations[CurrentSystem][k].Groups[l].CenterOfMassVelocity=vel;
              }
              else
              {
                if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
                  for(m=0;m<Components[Type].Groups[l].NumberOfGroupAtoms;m++)
                  {
                    A=Components[Type].Groups[l].Atoms[m];
                    vel=Cations[CurrentSystem][k].Atoms[A].Velocity;
                    tempv.x=eigenvectors.ax*vel.x+eigenvectors.bx*vel.y+eigenvectors.cx*vel.z;
                    tempv.y=eigenvectors.ay*vel.x+eigenvectors.by*vel.y+eigenvectors.cy*vel.z;
                    tempv.z=eigenvectors.az*vel.x+eigenvectors.bz*vel.y+eigenvectors.cz*vel.z;
                    tempv.x*=vexpdt.x;
                    tempv.y*=vexpdt.y;
                    tempv.z*=vexpdt.z;
                    vel.x=tempv.x*eigenvectors.ax+tempv.y*eigenvectors.ay+tempv.z*eigenvectors.az;
                    vel.y=tempv.x*eigenvectors.bx+tempv.y*eigenvectors.by+tempv.z*eigenvectors.bz;
                    vel.z=tempv.x*eigenvectors.cx+tempv.y*eigenvectors.cy+tempv.z*eigenvectors.cz;
                    Cations[CurrentSystem][k].Atoms[A].Velocity=vel;
                  }
              }
            }
          }
          break;
      }

      // update the force on the box
      UKineticTranslation=2.0*GetTranslationKineticEnergy();
      CellForce=GetKineticStressTensor();
      trace=UKineticTranslation/(REAL)DegreesOfFreedomTranslation[CurrentSystem];
      CellForce.ax+=trace;
      CellForce.by+=trace;
      CellForce.cz+=trace;

      CellForce.ax+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].ax-therm_baro_stats.ExternalPressure[CurrentSystem][0]);
      CellForce.ay+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].ay);
      CellForce.az+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].az);
      CellForce.bx+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].bx);
      CellForce.by+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].by-therm_baro_stats.ExternalPressure[CurrentSystem][0]);
      CellForce.bz+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].bz);
      CellForce.cx+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cx);
      CellForce.cy+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cy);
      CellForce.cz+=Volume[CurrentSystem]*(ConfigurationalStressTensor[CurrentSystem].cz-therm_baro_stats.ExternalPressure[CurrentSystem][0]);

      CellForce.ax/=CellMass[CurrentSystem]; CellForce.bx/=CellMass[CurrentSystem]; CellForce.cx/=CellMass[CurrentSystem];
      CellForce.ay/=CellMass[CurrentSystem]; CellForce.by/=CellMass[CurrentSystem]; CellForce.cy/=CellMass[CurrentSystem];
      CellForce.az/=CellMass[CurrentSystem]; CellForce.bz/=CellMass[CurrentSystem]; CellForce.cz/=CellMass[CurrentSystem];

      switch(NPTPRCellType[CurrentSystem])
      {
        case REGULAR_UPPER_TRIANGLE:
          CellForce.ay=CellForce.bz=CellForce.az=0.0;
          break;
        case ISOTROPIC:
          CellForce.ay=CellForce.bz=CellForce.az=0.0;
          CellForce.bx=CellForce.cy=CellForce.cx=0.0;
          avg=(CellForce.ax+CellForce.by+CellForce.cz)/3.0;
          CellForce.ax=CellForce.by=CellForce.cz=avg;
          break;
        case ANISOTROPIC:
          CellForce.ay=CellForce.bz=CellForce.az=0.0;
          CellForce.bx=CellForce.cy=CellForce.cx=0.0;
          break;
        case MONOCLINIC:
          CellForce.bx=CellForce.cy=CellForce.ay=CellForce.bz=0.0;
          avg=0.5*(CellForce.az+CellForce.cx);
          CellForce.az=CellForce.cx=avg;
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          CellForce.bx=CellForce.cy=CellForce.ay=CellForce.az=CellForce.bz=0.0;
          break;
        case REGULAR:
        default:
          avg=0.5*(CellForce.ay+CellForce.bx);
          CellForce.ay=CellForce.bx=avg;
          avg=0.5*(CellForce.az+CellForce.cx);
          CellForce.az=CellForce.cx=avg;
          avg=0.5*(CellForce.bz+CellForce.cy);
          CellForce.bz=CellForce.cy=avg;
          break;
      }

      // update the box velocities
      CellVelocity[CurrentSystem].ax+=CellForce.ax*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].ay+=CellForce.ay*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].az+=CellForce.az*w[j]*DeltaT/(4.0*(REAL)nc);

      CellVelocity[CurrentSystem].bx+=CellForce.bx*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].by+=CellForce.by*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].bz+=CellForce.bz*w[j]*DeltaT/(4.0*(REAL)nc);

      CellVelocity[CurrentSystem].cx+=CellForce.cx*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].cy+=CellForce.cy*w[j]*DeltaT/(4.0*(REAL)nc);
      CellVelocity[CurrentSystem].cz+=CellForce.cz*w[j]*DeltaT/(4.0*(REAL)nc);

    }
}

int versionNumber=1;

void WriteRestartThermoBarostats(FILE *FilePtr)
{
  int i;
  REAL Check;

  fwrite(&versionNumber,sizeof(int),1,FilePtr);

  fwrite(&therm_baro_stats,1,sizeof(therm_baro_stats),FilePtr);
  fwrite(&NumberOfIsothermPressures,1,sizeof(int),FilePtr);
  fwrite(&CurrentIsothermPressure,1,sizeof(int),FilePtr);

  fwrite(therm_baro_stats.ExternalTemperature,NumberOfSystems,sizeof(REAL),FilePtr);
  fwrite(therm_baro_stats.ExternalPressureTensor,NumberOfSystems,sizeof(REAL_MATRIX3x3),FilePtr);
  fwrite(therm_baro_stats.ExternalSurfaceTension,NumberOfSystems,sizeof(VECTOR),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
    fwrite(therm_baro_stats.ExternalPressure[i],NumberOfIsothermPressures,sizeof(REAL),FilePtr);
  fwrite(therm_baro_stats.ExternalStress,NumberOfSystems,sizeof(REAL_MATRIX3x3),FilePtr);

  fwrite(&w,1,sizeof(w),FilePtr);

  fwrite(&CellForcePressurePart,1,sizeof(CellForcePressurePart),FilePtr);
  fwrite(&CellForceKineticPressurePart,1,sizeof(CellForceKineticPressurePart),FilePtr);

  fwrite(CellVelocity,NumberOfSystems,sizeof(REAL_MATRIX3x3),FilePtr);
  fwrite(CellMass,NumberOfSystems,sizeof(REAL),FilePtr);

  fwrite(LnVolumePosition,NumberOfSystems,sizeof(REAL),FilePtr);
  fwrite(LnVolumeVelocity,NumberOfSystems,sizeof(REAL),FilePtr);
  fwrite(LnVolumeMass,NumberOfSystems,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    fwrite(ThermostatDegreesOfFreedom[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatForce[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatVelocity[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatPosition[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatMass[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fwrite(ThermostatDegreesOfFreedomTranslation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatForceTranslation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatVelocityTranslation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatPositionTranslation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatMassTranslation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fwrite(ThermostatDegreesOfFreedomTranslationFramework[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatForceTranslationFramework[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatVelocityTranslationFramework[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatPositionTranslationFramework[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatMassTranslationFramework[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fwrite(ThermostatDegreesOfFreedomTranslationAdsorbates[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatForceTranslationAdsorbates[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatVelocityTranslationAdsorbates[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatPositionTranslationAdsorbates[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatMassTranslationAdsorbates[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fwrite(ThermostatDegreesOfFreedomTranslationCations[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatForceTranslationCations[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatVelocityTranslationCations[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatPositionTranslationCations[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatMassTranslationCations[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fwrite(ThermostatDegreesOfFreedomRotation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatForceRotation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatVelocityRotation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatPositionRotation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fwrite(ThermostatMassRotation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fwrite(BarostatDegreesOfFreedom[i],therm_baro_stats.BarostatChainLength,sizeof(REAL),FilePtr);
    fwrite(BarostatForce[i],therm_baro_stats.BarostatChainLength,sizeof(REAL),FilePtr);
    fwrite(BarostatVelocity[i],therm_baro_stats.BarostatChainLength,sizeof(REAL),FilePtr);
    fwrite(BarostatPosition[i],therm_baro_stats.BarostatChainLength,sizeof(REAL),FilePtr);
    fwrite(BarostatMass[i],therm_baro_stats.BarostatChainLength,sizeof(REAL),FilePtr);
  }

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void AllocateThermoBaroStatMemory(void)
{
  int i;

  therm_baro_stats.ExternalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  therm_baro_stats.ExternalPressureTensor=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  therm_baro_stats.ExternalSurfaceTension=(VECTOR*)calloc(NumberOfSystems,sizeof(VECTOR));
  therm_baro_stats.ExternalPressure=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  for(i=0;i<NumberOfSystems;i++)
    therm_baro_stats.ExternalPressure[i]=(REAL*)calloc(NumberOfIsothermPressures,sizeof(REAL));
  therm_baro_stats.ExternalStress=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));

  CellVelocity=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  CellMass=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  LnVolumePosition=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  LnVolumeVelocity=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  LnVolumeMass=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  ThermostatDegreesOfFreedom=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatForce=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatVelocity=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatPosition=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatMass=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  ThermostatDegreesOfFreedomTranslation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatForceTranslation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatVelocityTranslation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatPositionTranslation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatMassTranslation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  ThermostatDegreesOfFreedomTranslationFramework=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatForceTranslationFramework=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatVelocityTranslationFramework=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatPositionTranslationFramework=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatMassTranslationFramework=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  ThermostatDegreesOfFreedomTranslationAdsorbates=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatForceTranslationAdsorbates=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatVelocityTranslationAdsorbates=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatPositionTranslationAdsorbates=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatMassTranslationAdsorbates=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  ThermostatDegreesOfFreedomTranslationCations=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatForceTranslationCations=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatVelocityTranslationCations=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatPositionTranslationCations=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatMassTranslationCations=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  ThermostatDegreesOfFreedomRotation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatForceRotation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatVelocityRotation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatPositionRotation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ThermostatMassRotation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  for(i=0;i<NumberOfSystems;i++)
  {
    ThermostatDegreesOfFreedom[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatForce[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatVelocity[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatPosition[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatMass[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));

    ThermostatDegreesOfFreedomTranslation[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatForceTranslation[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatVelocityTranslation[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatPositionTranslation[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatMassTranslation[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));

    ThermostatDegreesOfFreedomTranslationFramework[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatForceTranslationFramework[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatVelocityTranslationFramework[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatPositionTranslationFramework[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatMassTranslationFramework[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));

    ThermostatDegreesOfFreedomTranslationAdsorbates[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatForceTranslationAdsorbates[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatVelocityTranslationAdsorbates[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatPositionTranslationAdsorbates[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatMassTranslationAdsorbates[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));

    ThermostatDegreesOfFreedomTranslationCations[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatForceTranslationCations[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatVelocityTranslationCations[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatPositionTranslationCations[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatMassTranslationCations[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));

    ThermostatDegreesOfFreedomRotation[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatForceRotation[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatVelocityRotation[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatPositionRotation[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
    ThermostatMassRotation[i]=(REAL*)calloc(therm_baro_stats.ThermostatChainLength,sizeof(REAL));
  }

  BarostatDegreesOfFreedom=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BarostatForce=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BarostatVelocity=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BarostatPosition=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BarostatMass=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  for(i=0;i<NumberOfSystems;i++)
  {
    BarostatDegreesOfFreedom[i]=(REAL*)calloc(therm_baro_stats.BarostatChainLength,sizeof(REAL));
    BarostatForce[i]=(REAL*)calloc(therm_baro_stats.BarostatChainLength,sizeof(REAL));
    BarostatVelocity[i]=(REAL*)calloc(therm_baro_stats.BarostatChainLength,sizeof(REAL));
    BarostatPosition[i]=(REAL*)calloc(therm_baro_stats.BarostatChainLength,sizeof(REAL));
    BarostatMass[i]=(REAL*)calloc(therm_baro_stats.BarostatChainLength,sizeof(REAL));
  }
}

void ReadRestartThermoBarostats(FILE *FilePtr)
{
  int i;
  REAL Check;
  int readversionNumber=0;

  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumber)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  fread(&therm_baro_stats,1,sizeof(therm_baro_stats),FilePtr);
  fread(&NumberOfIsothermPressures,1,sizeof(int),FilePtr);
  fread(&CurrentIsothermPressure,1,sizeof(int),FilePtr);

  AllocateThermoBaroStatMemory();

  fread(therm_baro_stats.ExternalTemperature,NumberOfSystems,sizeof(REAL),FilePtr);
  fread(therm_baro_stats.ExternalPressureTensor,NumberOfSystems,sizeof(REAL_MATRIX3x3),FilePtr);
  fread(therm_baro_stats.ExternalSurfaceTension,NumberOfSystems,sizeof(VECTOR),FilePtr);
  for(i=0;i<NumberOfSystems;i++)
    fread(therm_baro_stats.ExternalPressure[i],NumberOfIsothermPressures,sizeof(REAL),FilePtr);
  fread(therm_baro_stats.ExternalStress,NumberOfSystems,sizeof(REAL_MATRIX3x3),FilePtr);

  fread(&w,1,sizeof(w),FilePtr);

  fread(&CellForcePressurePart,1,sizeof(CellForcePressurePart),FilePtr);
  fread(&CellForceKineticPressurePart,1,sizeof(CellForceKineticPressurePart),FilePtr);

  fread(CellVelocity,NumberOfSystems,sizeof(REAL_MATRIX3x3),FilePtr);
  fread(CellMass,NumberOfSystems,sizeof(REAL),FilePtr);

  fread(LnVolumePosition,NumberOfSystems,sizeof(REAL),FilePtr);
  fread(LnVolumeVelocity,NumberOfSystems,sizeof(REAL),FilePtr);
  fread(LnVolumeMass,NumberOfSystems,sizeof(REAL),FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    fread(ThermostatDegreesOfFreedom[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatForce[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatVelocity[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatPosition[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatMass[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fread(ThermostatDegreesOfFreedomTranslation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatForceTranslation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatVelocityTranslation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatPositionTranslation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatMassTranslation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fread(ThermostatDegreesOfFreedomTranslationFramework[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatForceTranslationFramework[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatVelocityTranslationFramework[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatPositionTranslationFramework[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatMassTranslationFramework[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fread(ThermostatDegreesOfFreedomTranslationAdsorbates[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatForceTranslationAdsorbates[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatVelocityTranslationAdsorbates[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatPositionTranslationAdsorbates[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatMassTranslationAdsorbates[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fread(ThermostatDegreesOfFreedomTranslationCations[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatForceTranslationCations[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatVelocityTranslationCations[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatPositionTranslationCations[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatMassTranslationCations[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fread(ThermostatDegreesOfFreedomRotation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatForceRotation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatVelocityRotation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatPositionRotation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);
    fread(ThermostatMassRotation[i],therm_baro_stats.ThermostatChainLength,sizeof(REAL),FilePtr);

    fread(BarostatDegreesOfFreedom[i],therm_baro_stats.BarostatChainLength,sizeof(REAL),FilePtr);
    fread(BarostatForce[i],therm_baro_stats.BarostatChainLength,sizeof(REAL),FilePtr);
    fread(BarostatVelocity[i],therm_baro_stats.BarostatChainLength,sizeof(REAL),FilePtr);
    fread(BarostatPosition[i],therm_baro_stats.BarostatChainLength,sizeof(REAL),FilePtr);
    fread(BarostatMass[i],therm_baro_stats.BarostatChainLength,sizeof(REAL),FilePtr);
  }

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartThermoBarostats)\n");
    ContinueAfterCrash=FALSE;
  }
}
