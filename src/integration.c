/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'integration.c' is part of RASPA-2.0

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
#include <math.h>
#include "constants.h"
#include "utils.h"
#include "simulation.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "framework_born.h"
#include "integration.h"
#include "ewald.h"
#include "internal_energy.h"
#include "internal_force.h"
#include "internal_born.h"
#include "inter_energy.h"
#include "inter_force.h"
#include "inter_born.h"
#include "thermo_baro_stats.h"
#include "rigid.h"
#include "mc_moves.h"
#include "potentials.h"
#include "sample.h"

void ComputeKineticEnergySystem(void)
{
  int i,j,l,f1;
  int Type,A;
  REAL Mass;
  VECTOR I,Velocity;
  QUATERNION p,q;
  INT_VECTOR3 Fixed;

  UHostKinetic[CurrentSystem]=0.0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
          if(!Framework[CurrentSystem].Atoms[f1][i].Fixed.x)
            UHostKinetic[CurrentSystem]+=0.5*Mass*SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.x);
          if(!Framework[CurrentSystem].Atoms[f1][i].Fixed.y)
            UHostKinetic[CurrentSystem]+=0.5*Mass*SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.y);
          if(!Framework[CurrentSystem].Atoms[f1][i].Fixed.z)
            UHostKinetic[CurrentSystem]+=0.5*Mass*SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.z);
        }
      }
    }
  }

  UAdsorbateTranslationalKinetic[CurrentSystem]=0.0;
  UAdsorbateRotationalKinetic[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        I=Components[Type].Groups[l].InverseInertiaVector;
        Velocity=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
        Mass=Components[Type].Groups[l].Mass;
        Fixed=Adsorbates[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x) UAdsorbateTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Velocity.x);
        if(!Fixed.y) UAdsorbateTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Velocity.y);
        if(!Fixed.z) UAdsorbateTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Velocity.z);

        Fixed=Adsorbates[CurrentSystem][i].Groups[l].FixedOrientation;
        if((!Fixed.x)||(!Fixed.y)||(!Fixed.z))
        {
          p=Adsorbates[CurrentSystem][i].Groups[l].QuaternionMomentum;
          q=Adsorbates[CurrentSystem][i].Groups[l].Quaternion;

          UAdsorbateRotationalKinetic[CurrentSystem]+=SQR(-p.r*q.i+p.i*q.r+p.j*q.k-p.k*q.j)*I.x/8.0;
          UAdsorbateRotationalKinetic[CurrentSystem]+=SQR(-p.r*q.j-p.i*q.k+p.j*q.r+p.k*q.i)*I.y/8.0;
          UAdsorbateRotationalKinetic[CurrentSystem]+=SQR(-p.r*q.k+p.i*q.j-p.j*q.i+p.k*q.r)*I.z/8.0;
        }
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;

          Fixed=Adsorbates[CurrentSystem][i].Atoms[A].Fixed;
          if(!Fixed.x) UAdsorbateTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x);
          if(!Fixed.y) UAdsorbateTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y);
          if(!Fixed.z) UAdsorbateTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z);
        }
      }
    }
  }

  UCationTranslationalKinetic[CurrentSystem]=0.0;
  UCationRotationalKinetic[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        I=Components[Cations[CurrentSystem][i].Type].Groups[l].InverseInertiaVector;
        Velocity=Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity;
        Mass=Components[Type].Groups[l].Mass;
        Fixed=Cations[CurrentSystem][i].Groups[l].FixedCenterOfMass;
        if(!Fixed.x) UCationTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Velocity.x);
        if(!Fixed.y) UCationTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Velocity.y);
        if(!Fixed.z) UCationTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Velocity.z);

        Fixed=Cations[CurrentSystem][i].Groups[l].FixedOrientation;
        if((!Fixed.x)||(!Fixed.y)||(!Fixed.z))
        {
          p=Cations[CurrentSystem][i].Groups[l].QuaternionMomentum;
          q=Cations[CurrentSystem][i].Groups[l].Quaternion;

          UCationRotationalKinetic[CurrentSystem]+=SQR(-p.r*q.i+p.i*q.r+p.j*q.k-p.k*q.j)*I.x/8.0;
          UCationRotationalKinetic[CurrentSystem]+=SQR(-p.r*q.j-p.i*q.k+p.j*q.r+p.k*q.i)*I.y/8.0;
          UCationRotationalKinetic[CurrentSystem]+=SQR(-p.r*q.k+p.i*q.j-p.j*q.i+p.k*q.r)*I.z/8.0;
        }
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[A].Type].Mass;

          Fixed=Cations[CurrentSystem][i].Atoms[A].Fixed;
          if(!Fixed.x) UCationTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Cations[CurrentSystem][i].Atoms[A].Velocity.x);
          if(!Fixed.y) UCationTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Cations[CurrentSystem][i].Atoms[A].Velocity.y);
          if(!Fixed.z) UCationTranslationalKinetic[CurrentSystem]+=0.5*Mass*SQR(Cations[CurrentSystem][i].Atoms[A].Velocity.z);
        }
      }
    }
  }

  UAdsorbateKinetic[CurrentSystem]=UAdsorbateTranslationalKinetic[CurrentSystem]+
                                   UAdsorbateRotationalKinetic[CurrentSystem];
  UCationKinetic[CurrentSystem]=UCationTranslationalKinetic[CurrentSystem]+
                                UCationRotationalKinetic[CurrentSystem];
  UKinetic[CurrentSystem]=UHostKinetic[CurrentSystem]+
       UAdsorbateKinetic[CurrentSystem]+UCationKinetic[CurrentSystem];
}

void InitializeForces(void)
{
  int i;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    ComputeQuaternionAdsorbate(i);

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    ComputeQuaternionCation(i);

  CalculateForce();

  CreateCartesianVelocities();

  ComputeKineticEnergySystem();

  ComputeNoseHooverEnergySystem();

  ConservedEnergy[CurrentSystem]=
           UTotal[CurrentSystem]+
           UKinetic[CurrentSystem]+
           UNoseHoover[CurrentSystem];
}


void InitializeForcesAllSystems(void)
{
  int i;

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      ComputeQuaternionAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      ComputeQuaternionCation(i);

    CreateCartesianVelocities();

    CalculateForce();

    ComputeKineticEnergySystem();

    ComputeNoseHooverEnergySystem();

    ConservedEnergy[CurrentSystem]=
             UTotal[CurrentSystem]+
             UKinetic[CurrentSystem]+
             UNoseHoover[CurrentSystem];

  }
  CurrentSystem=0;
}

void CalculateTotalEnergyAllSystems(void)
{
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    CalculateEnergy();
  CurrentSystem=0;
}

void CalculateEnergy(void)
{

  // initialize energies
  UAdsorbateBond[CurrentSystem]=0.0;
  UAdsorbateUreyBradley[CurrentSystem]=0.0;
  UAdsorbateBend[CurrentSystem]=0.0;
  UAdsorbateInversionBend[CurrentSystem]=0.0;
  UAdsorbateTorsion[CurrentSystem]=0.0;
  UAdsorbateImproperTorsion[CurrentSystem]=0.0;
  UAdsorbateBondBond[CurrentSystem]=0.0;
  UAdsorbateBondBend[CurrentSystem]=0.0;
  UAdsorbateBendBend[CurrentSystem]=0.0;
  UAdsorbateBondTorsion[CurrentSystem]=0.0;
  UAdsorbateBendTorsion[CurrentSystem]=0.0;
  UAdsorbateIntraVDW[CurrentSystem]=0.0;
  UAdsorbateIntraChargeCharge[CurrentSystem]=0.0;
  UAdsorbateIntraChargeBondDipole[CurrentSystem]=0.0;
  UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]=0.0;

  UCationBond[CurrentSystem]=0.0;
  UCationUreyBradley[CurrentSystem]=0.0;
  UCationBend[CurrentSystem]=0.0;
  UCationInversionBend[CurrentSystem]=0.0;
  UCationTorsion[CurrentSystem]=0.0;
  UCationImproperTorsion[CurrentSystem]=0.0;
  UCationBondBond[CurrentSystem]=0.0;
  UCationBondBend[CurrentSystem]=0.0;
  UCationBendBend[CurrentSystem]=0.0;
  UCationBondTorsion[CurrentSystem]=0.0;
  UCationBendTorsion[CurrentSystem]=0.0;
  UCationIntraVDW[CurrentSystem]=0.0;
  UCationIntraChargeCharge[CurrentSystem]=0.0;
  UCationIntraChargeBondDipole[CurrentSystem]=0.0;
  UCationIntraBondDipoleBondDipole[CurrentSystem]=0.0;

  UHostBond[CurrentSystem]=0.0;
  UHostUreyBradley[CurrentSystem]=0.0;
  UHostBend[CurrentSystem]=0.0;
  UHostInversionBend[CurrentSystem]=0.0;
  UHostTorsion[CurrentSystem]=0.0;
  UHostImproperTorsion[CurrentSystem]=0.0;
  UHostBondBond[CurrentSystem]=0.0;
  UHostBondBend[CurrentSystem]=0.0;
  UHostBendBend[CurrentSystem]=0.0;
  UHostBondTorsion[CurrentSystem]=0.0;
  UHostBendTorsion[CurrentSystem]=0.0;

  UHostHost[CurrentSystem]=0.0;
  UAdsorbateAdsorbate[CurrentSystem]=0.0;
  UCationCation[CurrentSystem]=0.0;
  UHostAdsorbate[CurrentSystem]=0.0;
  UHostCation[CurrentSystem]=0.0;
  UAdsorbateCation[CurrentSystem]=0.0;

  UHostHostVDW[CurrentSystem]=0.0;
  UAdsorbateAdsorbateVDW[CurrentSystem]=0.0;
  UCationCationVDW[CurrentSystem]=0.0;
  UHostAdsorbateVDW[CurrentSystem]=0.0;
  UHostCationVDW[CurrentSystem]=0.0;
  UAdsorbateCationVDW[CurrentSystem]=0.0;

  UHostHostCoulomb[CurrentSystem]=0.0;
  UAdsorbateAdsorbateCoulomb[CurrentSystem]=0.0;
  UCationCationCoulomb[CurrentSystem]=0.0;
  UHostAdsorbateCoulomb[CurrentSystem]=0.0;
  UHostCationCoulomb[CurrentSystem]=0.0;
  UAdsorbateCationCoulomb[CurrentSystem]=0.0;

  UHostHostChargeChargeReal[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  UCationCationChargeChargeReal[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  UHostCationChargeChargeReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeReal[CurrentSystem]=0.0;

  UHostHostChargeBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleReal[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleReal[CurrentSystem]=0.0;

  UHostHostBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;

  UHostHostChargeChargeFourier[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UCationCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UHostCationChargeChargeFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourier[CurrentSystem]=0.0;

  UHostHostChargeBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=0.0;

  UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;

  UHostPolarization[CurrentSystem]=0.0;
  UAdsorbatePolarization[CurrentSystem]=0.0;
  UCationPolarization[CurrentSystem]=0.0;

  UHostBackPolarization[CurrentSystem]=0.0;
  UAdsorbateBackPolarization[CurrentSystem]=0.0;
  UCationBackPolarization[CurrentSystem]=0.0;

  UTailCorrection[CurrentSystem]=0.0;
  UDistanceConstraints[CurrentSystem]=0.0;
  UAngleConstraints[CurrentSystem]=0.0;
  UDihedralConstraints[CurrentSystem]=0.0;
  UInversionBendConstraints[CurrentSystem]=0.0;
  UOutOfPlaneDistanceConstraints[CurrentSystem]=0.0;
  UExclusionConstraints[CurrentSystem]=0.0;
  UNoseHoover[CurrentSystem]=0.0;
  UTotal[CurrentSystem]=0.0;

  // compute the modified VDW sites for anisotropic models
  CalculateAnisotropicSites();

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    UHostBond[CurrentSystem]=CalculateFrameworkBondEnergy(TRUE,0,0);
    UHostUreyBradley[CurrentSystem]=CalculateFrameworkUreyBradleyEnergy(TRUE,0,0);
    UHostBend[CurrentSystem]=CalculateFrameworkBendEnergy(TRUE,0,0);
    UHostInversionBend[CurrentSystem]=CalculateFrameworkInversionBendEnergy(TRUE,0,0);
    UHostTorsion[CurrentSystem]=CalculateFrameworkTorsionEnergy(TRUE,0,0);
    UHostImproperTorsion[CurrentSystem]=CalculateFrameworkImproperTorsionEnergy(TRUE,0,0);
    UHostBondBond[CurrentSystem]=CalculateFrameworkBondBondEnergy(TRUE,0,0);
    UHostBondBend[CurrentSystem]=CalculateFrameworkBondBendEnergy(TRUE,0,0);
    UHostBendBend[CurrentSystem]=CalculateFrameworkBendBendEnergy(TRUE,0,0);
    UHostBondTorsion[CurrentSystem]=CalculateFrameworkBondTorsionEnergy(TRUE,0,0);
    UHostBendTorsion[CurrentSystem]=CalculateFrameworkBendTorsionEnergy(TRUE,0,0);

    if(UseReplicas[CurrentSystem])
    {
      CalculateFrameworkIntraReplicaVDWEnergy();
      CalculateFrameworkIntraReplicaChargeChargeEnergy();
    }
    else
    {
      CalculateFrameworkIntraVDWEnergy();
      CalculateFrameworkIntraChargeChargeEnergy();
    }
    CalculateFrameworkIntraChargeBondDipoleEnergy();
    CalculateFrameworkIntraBondDipoleBondDipoleEnergy();
  }

  CalculateBondEnergyAdsorbates();
  CalculateUreyBradleyEnergyAdsorbates();
  CalculateBendEnergyAdsorbates();
  CalculateInversionBendEnergyAdsorbates();
  CalculateTorsionEnergyAdsorbates();
  CalculateImproperTorsionEnergyAdsorbates();
  CalculateBondBondEnergyAdsorbates();
  CalculateBondBendEnergyAdsorbates();
  CalculateBendBendEnergyAdsorbates();
  CalculateBendTorsionEnergyAdsorbates();
  CalculateBondTorsionEnergyAdsorbates();

  CalculateIntraVDWEnergyAdsorbates();
  CalculateIntraChargeChargeEnergyAdsorbates();
  CalculateIntraChargeBondDipoleEnergyAdsorbates();
  CalculateIntraBondDipoleBondDipoleEnergyAdsorbates();

  CalculateBondEnergyCations();
  CalculateUreyBradleyEnergyCations();
  CalculateBendEnergyCations();
  CalculateInversionBendEnergyCations();
  CalculateTorsionEnergyCations();
  CalculateImproperTorsionEnergyCations();
  CalculateBondBondEnergyCations();
  CalculateBondBendEnergyCations();
  CalculateBendBendEnergyCations();
  CalculateBondTorsionEnergyCations();
  CalculateBendTorsionEnergyCations();

  CalculateIntraVDWEnergyCations();

  CalculateIntraChargeChargeEnergyCations();
  CalculateIntraChargeBondDipoleEnergyCations();
  CalculateIntraBondDipoleBondDipoleEnergyCations();

  // contribution of the intermolecular interactions (between adsorbates and/or cations)
  CalculateTotalInterVDWEnergy();
  CalculateTotalInterChargeChargeCoulombEnergy();
  CalculateTotalInterChargeBondDipoleCoulombEnergy();
  CalculateTotalInterBondDipoleBondDipoleCoulombEnergy();

  // contribution of the interactions of the adsorbates with the framework
  if(UseReplicas[CurrentSystem])
  {
    CalculateFrameworkAdsorbateReplicaVDWEnergy();
    CalculateFrameworkAdsorbateReplicaChargeChargeEnergy();
  }
  else
  {
    CalculateFrameworkAdsorbateVDWEnergy();
    CalculateFrameworkAdsorbateChargeChargeEnergy();
  }
  CalculateFrameworkAdsorbateChargeBondDipoleEnergy();
  CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergy();

  // contribution of the interactions of the cations with the framework
  if(UseReplicas[CurrentSystem])
  {
    CalculateFrameworkCationReplicaVDWEnergy();
    CalculateFrameworkCationReplicaChargeChargeEnergy();
  }
  else
  {
    CalculateFrameworkCationVDWEnergy();
    CalculateFrameworkCationChargeChargeEnergy();
  }
  CalculateFrameworkCationChargeBondDipoleEnergy();
  CalculateFrameworkCationBondDipoleBondDipoleEnergy();

  // the contribution of the charges present in the system (long-range interaction)
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    EwaldFourierEnergy();

  CalculateTailCorrection();

  CalculateHarmonicBondConstraintEnergy();
  CalculateHarmonicAngleConstraintEnergy();
  CalculateHarmonicDihedralConstraintEnergy();
  CalculateConstraintsExclusionEnergy();

  if(ComputePolarization)
    CalculateElectricField();

  UAdsorbateAdsorbateCoulomb[CurrentSystem]=
      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
  UAdsorbateAdsorbate[CurrentSystem]=UAdsorbateAdsorbateVDW[CurrentSystem]+UAdsorbateAdsorbateCoulomb[CurrentSystem];

  UCationCationCoulomb[CurrentSystem]=
      UCationCationChargeChargeReal[CurrentSystem]+UCationCationChargeChargeFourier[CurrentSystem]+
      UCationCationChargeBondDipoleReal[CurrentSystem]+UCationCationChargeBondDipoleFourier[CurrentSystem]+
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+UCationCationBondDipoleBondDipoleFourier[CurrentSystem];
  UCationCation[CurrentSystem]=UCationCationVDW[CurrentSystem]+UCationCationCoulomb[CurrentSystem];

  UAdsorbateCationCoulomb[CurrentSystem]=
      UAdsorbateCationChargeChargeReal[CurrentSystem]+UAdsorbateCationChargeChargeFourier[CurrentSystem]+
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem];
  UAdsorbateCation[CurrentSystem]=UAdsorbateCationVDW[CurrentSystem]+UAdsorbateCationCoulomb[CurrentSystem];

  UHostAdsorbateCoulomb[CurrentSystem]=
      UHostAdsorbateChargeChargeReal[CurrentSystem]+UHostAdsorbateChargeChargeFourier[CurrentSystem]+
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
  UHostAdsorbate[CurrentSystem]=UHostAdsorbateVDW[CurrentSystem]+UHostAdsorbateCoulomb[CurrentSystem];

  UHostCationCoulomb[CurrentSystem]=
      UHostCationChargeChargeReal[CurrentSystem]+UHostCationChargeChargeFourier[CurrentSystem]+
      UHostCationChargeBondDipoleReal[CurrentSystem]+UHostCationChargeBondDipoleFourier[CurrentSystem]+
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+UHostCationBondDipoleBondDipoleFourier[CurrentSystem];
  UHostCation[CurrentSystem]=UHostCationVDW[CurrentSystem]+UHostCationCoulomb[CurrentSystem];

  UHostHostCoulomb[CurrentSystem]=
      UHostHostChargeChargeReal[CurrentSystem]+UHostHostChargeChargeFourier[CurrentSystem]+
      UHostHostChargeBondDipoleReal[CurrentSystem]+UHostHostChargeBondDipoleFourier[CurrentSystem]+
      UHostHostBondDipoleBondDipoleReal[CurrentSystem]+UHostHostBondDipoleBondDipoleFourier[CurrentSystem];
  UHostHost[CurrentSystem]=UHostHostVDW[CurrentSystem]+UHostHostCoulomb[CurrentSystem];

  UTotal[CurrentSystem]=
         // inter
         UHostHost[CurrentSystem]+UCationCation[CurrentSystem]+UAdsorbateAdsorbate[CurrentSystem]+
         UHostCation[CurrentSystem]+UHostAdsorbate[CurrentSystem]+UAdsorbateCation[CurrentSystem]+
         // intra framework
         UHostBond[CurrentSystem]+UHostUreyBradley[CurrentSystem]+UHostBend[CurrentSystem]+
         UHostInversionBend[CurrentSystem]+UHostTorsion[CurrentSystem]+UHostImproperTorsion[CurrentSystem]+
         UHostBondBond[CurrentSystem]+UHostBondBend[CurrentSystem]+UHostBendBend[CurrentSystem]+
         UHostBondTorsion[CurrentSystem]+UHostBendTorsion[CurrentSystem]+
         // intra adsorbates
         UAdsorbateBond[CurrentSystem]+UAdsorbateUreyBradley[CurrentSystem]+UAdsorbateBend[CurrentSystem]+
         UAdsorbateInversionBend[CurrentSystem]+UAdsorbateTorsion[CurrentSystem]+UAdsorbateImproperTorsion[CurrentSystem]+
         UAdsorbateBondBond[CurrentSystem]+UAdsorbateBondBend[CurrentSystem]+UAdsorbateBendBend[CurrentSystem]+
         UAdsorbateBondTorsion[CurrentSystem]+UAdsorbateBendTorsion[CurrentSystem]+
         UAdsorbateIntraVDW[CurrentSystem]+UAdsorbateIntraChargeCharge[CurrentSystem]+
         UAdsorbateIntraChargeBondDipole[CurrentSystem]+UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+
         // intra cations
         UCationBond[CurrentSystem]+UCationUreyBradley[CurrentSystem]+UCationBend[CurrentSystem]+
         UCationInversionBend[CurrentSystem]+UCationTorsion[CurrentSystem]+UCationImproperTorsion[CurrentSystem]+
         UCationBondBond[CurrentSystem]+UCationBondBend[CurrentSystem]+UCationBendBend[CurrentSystem]+
         UCationBondTorsion[CurrentSystem]+UCationBendTorsion[CurrentSystem]+
         UCationIntraVDW[CurrentSystem]+UCationIntraChargeCharge[CurrentSystem]+
         UCationIntraChargeBondDipole[CurrentSystem]+UCationIntraBondDipoleBondDipole[CurrentSystem]+
         // polarization
         UHostPolarization[CurrentSystem]+UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem]+
         UHostBackPolarization[CurrentSystem]+UAdsorbateBackPolarization[CurrentSystem]+UCationBackPolarization[CurrentSystem]+
         // tailcorrection
         UTailCorrection[CurrentSystem]+
         // harmonic constraint energy
         UDistanceConstraints[CurrentSystem]+UAngleConstraints[CurrentSystem]+UDihedralConstraints[CurrentSystem]+
         UInversionBendConstraints[CurrentSystem]+UOutOfPlaneDistanceConstraints[CurrentSystem]+UExclusionConstraints[CurrentSystem];
}

void CalculateForce(void)
{
  int i,j,f1;
  VECTOR com,pos,force;
  REAL temp;
  int Type,A;
  int k,l;

  // initialize stress
  StressTensor[CurrentSystem].ax=StressTensor[CurrentSystem].bx=StressTensor[CurrentSystem].cx=0.0;
  StressTensor[CurrentSystem].ay=StressTensor[CurrentSystem].by=StressTensor[CurrentSystem].cy=0.0;
  StressTensor[CurrentSystem].az=StressTensor[CurrentSystem].bz=StressTensor[CurrentSystem].cz=0.0;

  StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
  StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
  StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

  StrainDerivativeRigidCorrection[CurrentSystem].ax=StrainDerivativeRigidCorrection[CurrentSystem].bx=StrainDerivativeRigidCorrection[CurrentSystem].cx=0.0;
  StrainDerivativeRigidCorrection[CurrentSystem].ay=StrainDerivativeRigidCorrection[CurrentSystem].by=StrainDerivativeRigidCorrection[CurrentSystem].cy=0.0;
  StrainDerivativeRigidCorrection[CurrentSystem].az=StrainDerivativeRigidCorrection[CurrentSystem].bz=StrainDerivativeRigidCorrection[CurrentSystem].cz=0.0;

  // initialize energies
  UAdsorbateBond[CurrentSystem]=0.0;
  UAdsorbateUreyBradley[CurrentSystem]=0.0;
  UAdsorbateBend[CurrentSystem]=0.0;
  UAdsorbateInversionBend[CurrentSystem]=0.0;
  UAdsorbateTorsion[CurrentSystem]=0.0;
  UAdsorbateImproperTorsion[CurrentSystem]=0.0;
  UAdsorbateBondBond[CurrentSystem]=0.0;
  UAdsorbateBondBend[CurrentSystem]=0.0;
  UAdsorbateBendBend[CurrentSystem]=0.0;
  UAdsorbateBondTorsion[CurrentSystem]=0.0;
  UAdsorbateBendTorsion[CurrentSystem]=0.0;
  UAdsorbateIntraVDW[CurrentSystem]=0.0;
  UAdsorbateIntraChargeCharge[CurrentSystem]=0.0;
  UAdsorbateIntraChargeBondDipole[CurrentSystem]=0.0;
  UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]=0.0;

  UCationBond[CurrentSystem]=0.0;
  UCationUreyBradley[CurrentSystem]=0.0;
  UCationBend[CurrentSystem]=0.0;
  UCationInversionBend[CurrentSystem]=0.0;
  UCationTorsion[CurrentSystem]=0.0;
  UCationImproperTorsion[CurrentSystem]=0.0;
  UCationBondBond[CurrentSystem]=0.0;
  UCationBondBend[CurrentSystem]=0.0;
  UCationBendBend[CurrentSystem]=0.0;
  UCationBondTorsion[CurrentSystem]=0.0;
  UCationBendTorsion[CurrentSystem]=0.0;
  UCationIntraVDW[CurrentSystem]=0.0;
  UCationIntraChargeCharge[CurrentSystem]=0.0;
  UCationIntraChargeBondDipole[CurrentSystem]=0.0;
  UCationIntraBondDipoleBondDipole[CurrentSystem]=0.0;

  UHostBond[CurrentSystem]=0.0;
  UHostUreyBradley[CurrentSystem]=0.0;
  UHostBend[CurrentSystem]=0.0;
  UHostInversionBend[CurrentSystem]=0.0;
  UHostTorsion[CurrentSystem]=0.0;
  UHostImproperTorsion[CurrentSystem]=0.0;
  UHostBondBond[CurrentSystem]=0.0;
  UHostBondBend[CurrentSystem]=0.0;
  UHostBendBend[CurrentSystem]=0.0;
  UHostBondTorsion[CurrentSystem]=0.0;
  UHostBendTorsion[CurrentSystem]=0.0;

  UHostHost[CurrentSystem]=0.0;
  UAdsorbateAdsorbate[CurrentSystem]=0.0;
  UCationCation[CurrentSystem]=0.0;
  UHostAdsorbate[CurrentSystem]=0.0;
  UHostCation[CurrentSystem]=0.0;
  UAdsorbateCation[CurrentSystem]=0.0;

  UHostHostVDW[CurrentSystem]=0.0;
  UAdsorbateAdsorbateVDW[CurrentSystem]=0.0;
  UCationCationVDW[CurrentSystem]=0.0;
  UHostAdsorbateVDW[CurrentSystem]=0.0;
  UHostCationVDW[CurrentSystem]=0.0;
  UAdsorbateCationVDW[CurrentSystem]=0.0;

  UHostHostCoulomb[CurrentSystem]=0.0;
  UAdsorbateAdsorbateCoulomb[CurrentSystem]=0.0;
  UCationCationCoulomb[CurrentSystem]=0.0;
  UHostAdsorbateCoulomb[CurrentSystem]=0.0;
  UHostCationCoulomb[CurrentSystem]=0.0;
  UAdsorbateCationCoulomb[CurrentSystem]=0.0;

  UHostHostChargeChargeReal[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  UCationCationChargeChargeReal[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  UHostCationChargeChargeReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeReal[CurrentSystem]=0.0;

  UHostHostChargeBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleReal[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleReal[CurrentSystem]=0.0;

  UHostHostBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;

  UHostHostChargeChargeFourier[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UCationCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UHostCationChargeChargeFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourier[CurrentSystem]=0.0;

  UHostHostChargeBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=0.0;

  UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;

  UHostPolarization[CurrentSystem]=0.0;
  UAdsorbatePolarization[CurrentSystem]=0.0;
  UCationPolarization[CurrentSystem]=0.0;

  UHostBackPolarization[CurrentSystem]=0.0;
  UAdsorbateBackPolarization[CurrentSystem]=0.0;
  UCationBackPolarization[CurrentSystem]=0.0;

  UTailCorrection[CurrentSystem]=0.0;
  UDistanceConstraints[CurrentSystem]=0.0;
  UAngleConstraints[CurrentSystem]=0.0;
  UDihedralConstraints[CurrentSystem]=0.0;
  UInversionBendConstraints[CurrentSystem]=0.0;
  UOutOfPlaneDistanceConstraints[CurrentSystem]=0.0;
  UExclusionConstraints[CurrentSystem]=0.0;
  UNoseHoover[CurrentSystem]=0.0;
  UTotal[CurrentSystem]=0.0;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Framework[CurrentSystem].Atoms[f1][i].Force.x=0.0;
        Framework[CurrentSystem].Atoms[f1][i].Force.y=0.0;
        Framework[CurrentSystem].Atoms[f1][i].Force.z=0.0;
      }
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

  // create atomic positions from com-positions (needed to compute the atomicly defined forces)
  CreateCartesianPositions();

  CalculateAnisotropicSites();

  // compute the current bond-dipole vectors from the positions of the individual bond atoms
  ConstructBondDipolesFromBondsFramework();

  if(ComputeBornTerm)
  {
    InitializeMatrix9x9(&BornTerm[CurrentSystem]);

    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBondBornTerm();
      CalculateFrameworkUreyBradleyBornTerm();
      CalculateFrameworkBendBornTerm();
      CalculateFrameworkTorsionBornTerm();
      CalculateFrameworkImproperTorsionBornTerm();

      if(UseReplicas[CurrentSystem])
      {
        CalculateFrameworkIntraReplicaVDWBornTerm();
        CalculateFrameworkIntraReplicaChargeChargeBornTerm();
      }
      else
      {
        CalculateFrameworkIntraVDWBornTerm();
        CalculateFrameworkIntraChargeChargeBornTerm();
      }
    }

    CalculateAdsorbateBondBornTerm();
    CalculateAdsorbateUreyBradleyBornTerm();
    CalculateAdsorbateBendBornTerm();
    CalculateAdsorbateTorsionBornTerm();
    CalculateAdsorbateImproperTorsionBornTerm();
    CalculateAdsorbateIntraVDWBornTerm();
    CalculateAdsorbateIntraCoulombBornTerm();

    CalculateCationBondBornTerm();
    CalculateCationUreyBradleyBornTerm();
    CalculateCationBendBornTerm();
    CalculateCationTorsionBornTerm();
    CalculateCationImproperTorsionBornTerm();
    CalculateCationIntraVDWBornTerm();
    CalculateCationIntraCoulombBornTerm();

    ComputeInterVDWBornTerm();
    ComputeInterChargeChargeBornTerm();

    if(UseReplicas[CurrentSystem])
    {
      FrameworkAdsorbateReplicaVDWBornTerm();
      FrameworkAdsorbateReplicaChargeChargeBornTerm();
      FrameworkCationReplicaVDWBornTerm();
      FrameworkCationReplicaChargeChargeBornTerm();
    }
    else
    {
      FrameworkAdsorbateVDWBornTerm();
      FrameworkAdsorbateChargeChargeBornTerm();
      FrameworkCationVDWBornTerm();
      FrameworkCationChargeChargeBornTerm();
    }

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      EwaldFourierBornTerm();

    FillInRemainingTermsInBornTerm();
  }
  else
  {
    // contribution of the intra molecular interactions of the framework
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      CalculateFrameworkBondForce();
      CalculateFrameworkUreyBradleyForce();
      CalculateFrameworkBendForce();
      CalculateFrameworkInversionBendForce();
      CalculateFrameworkTorsionForce();
      CalculateFrameworkImproperTorsionForce();
      CalculateFrameworkBondBondForce();
      CalculateFrameworkBondBendForce();
      CalculateFrameworkBendBendForce();
      CalculateFrameworkBondTorsionForce();
      CalculateFrameworkBendTorsionForce();

      if(UseReplicas[CurrentSystem])
      {
        CalculateFrameworkIntraReplicaVDWForce();
        CalculateFrameworkIntraReplicaChargeChargeForce();
      }
      else
      {
        CalculateFrameworkIntraVDWForce();
        CalculateFrameworkIntraChargeChargeForce();
      }
      CalculateFrameworkIntraChargeBondDipoleForce();
      CalculateFrameworkIntraBondDipoleBondDipoleForce();
    }

    // contribution of the intra molecular interactions of the adsorbates
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      CalculateAdsorbateBondForce(i);
      CalculateAdsorbateUreyBradleyForce(i);
      CalculateAdsorbateBendForce(i);
      CalculateAdsorbateInversionBendForce(i);
      CalculateAdsorbateTorsionForce(i);
      CalculateAdsorbateImproperTorsionForce(i);
      CalculateAdsorbateBondBondForce(i);
      CalculateAdsorbateBondBendForce(i);
      CalculateAdsorbateBendBendForce(i);
      CalculateAdsorbateBondTorsionForce(i);
      CalculateAdsorbateBendTorsionForce(i);

      CalculateAdsorbateIntraVDWForce(i);
      CalculateAdsorbateIntraChargeChargeForce(i);
      CalculateAdsorbateIntraChargeBondDipoleForce(i);
      CalculateAdsorbateIntraBondDipoleBondDipoleForce(i);
    }
    // contribution of the intra molecular interactions of the cations
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      CalculateCationBondForce(i);
      CalculateCationUreyBradleyForce(i);
      CalculateCationBendForce(i);
      CalculateCationInversionBendForce(i);
      CalculateCationTorsionForce(i);
      CalculateCationImproperTorsionForce(i);
      CalculateCationBondBondForce(i);
      CalculateCationBondBendForce(i);
      CalculateCationBendBendForce(i);
      CalculateCationBondTorsionForce(i);
      CalculateCationBendTorsionForce(i);

      CalculateCationIntraVDWForce(i);

      CalculateCationIntraChargeChargeForce(i);
      CalculateCationIntraChargeBondDipoleForce(i);
      CalculateCationIntraBondDipoleBondDipoleForce(i);
    }

    // contribution of the intermolecular interactions (between adsorbates and/or cations)
    if(UseReplicas[CurrentSystem])
    {
      CalculateTotalInterReplicaVDWForce();
      CalculateTotalInterReplicaChargeChargeCoulombForce();
    }
    else
    {
      CalculateTotalInterVDWForce();
      CalculateTotalInterChargeChargeCoulombForce();

    }
    CalculateTotalInterChargeBondDipoleCoulombForce();
    CalculateTotalInterBondDipoleBondDipoleCoulombForce();

    // contribution of the interactions of the adsorbates with the framework
    if(UseReplicas[CurrentSystem])
    {
      CalculateFrameworkAdsorbateReplicaVDWForce();
      CalculateFrameworkAdsorbateReplicaChargeChargeForce();
    }
    else
    {
      CalculateFrameworkAdsorbateVDWForce();
      CalculateFrameworkAdsorbateChargeChargeForce();
    }
    CalculateFrameworkAdsorbateChargeBondDipoleForce();
    CalculateFrameworkAdsorbateBondDipoleBondDipoleForce();

    // contribution of the interactions of the cations with the framework
    if(UseReplicas[CurrentSystem])
    {
      CalculateFrameworkCationReplicaVDWForce();
      CalculateFrameworkCationReplicaChargeChargeForce();
    }
    else
    {
      CalculateFrameworkCationVDWForce();
      CalculateFrameworkCationChargeChargeForce();
    }
    CalculateFrameworkCationChargeBondDipoleForce();
    CalculateFrameworkCationBondDipoleBondDipoleForce();

    // the contribution of the charges present in the system (long-range interaction)
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      EwaldFourierForce();
  }

  if(ComputePolarization)
  {
    CalculateElectricField();

    CalculateFrameworkAdsorbateChargeInducedDipoleForce();
    CalculateFrameworkAdsorbateInducedDipoleInducedDipoleForce();

    CalculateInterChargeInducedDipoleForce();
    CalculateInterInducedDipoleInducedDipoleForce();

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      ComputeInducedDipolesForcesEwald();
  }

  CreateQuaternionAndComForces();

  CalculateTailCorrection();

  CalculateHarmonicBondConstraintForce();
  CalculateHarmonicAngleConstraintForce();
  CalculateHarmonicDihedralConstraintForce();
  CalculateConstraintsExclusionForce();


  UAdsorbateAdsorbateCoulomb[CurrentSystem]=
      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
  UAdsorbateAdsorbate[CurrentSystem]=UAdsorbateAdsorbateVDW[CurrentSystem]+UAdsorbateAdsorbateCoulomb[CurrentSystem];

  UCationCationCoulomb[CurrentSystem]=
      UCationCationChargeChargeReal[CurrentSystem]+UCationCationChargeChargeFourier[CurrentSystem]+
      UCationCationChargeBondDipoleReal[CurrentSystem]+UCationCationChargeBondDipoleFourier[CurrentSystem]+
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+UCationCationBondDipoleBondDipoleFourier[CurrentSystem];
  UCationCation[CurrentSystem]=UCationCationVDW[CurrentSystem]+UCationCationCoulomb[CurrentSystem];

  UAdsorbateCationCoulomb[CurrentSystem]=
      UAdsorbateCationChargeChargeReal[CurrentSystem]+UAdsorbateCationChargeChargeFourier[CurrentSystem]+
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem];
  UAdsorbateCation[CurrentSystem]=UAdsorbateCationVDW[CurrentSystem]+UAdsorbateCationCoulomb[CurrentSystem];

  UHostAdsorbateCoulomb[CurrentSystem]=
      UHostAdsorbateChargeChargeReal[CurrentSystem]+UHostAdsorbateChargeChargeFourier[CurrentSystem]+
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
  UHostAdsorbate[CurrentSystem]=UHostAdsorbateVDW[CurrentSystem]+UHostAdsorbateCoulomb[CurrentSystem];

  UHostCationCoulomb[CurrentSystem]=
      UHostCationChargeChargeReal[CurrentSystem]+UHostCationChargeChargeFourier[CurrentSystem]+
      UHostCationChargeBondDipoleReal[CurrentSystem]+UHostCationChargeBondDipoleFourier[CurrentSystem]+
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+UHostCationBondDipoleBondDipoleFourier[CurrentSystem];
  UHostCation[CurrentSystem]=UHostCationVDW[CurrentSystem]+UHostCationCoulomb[CurrentSystem];

  UHostHostCoulomb[CurrentSystem]=
      UHostHostChargeChargeReal[CurrentSystem]+UHostHostChargeChargeFourier[CurrentSystem]+
      UHostHostChargeBondDipoleReal[CurrentSystem]+UHostHostChargeBondDipoleFourier[CurrentSystem]+
      UHostHostBondDipoleBondDipoleReal[CurrentSystem]+UHostHostBondDipoleBondDipoleFourier[CurrentSystem];
  UHostHost[CurrentSystem]=UHostHostVDW[CurrentSystem]+UHostHostCoulomb[CurrentSystem];

  // correct the strain derivative tensor for rigid units
  // first the tensor was computed atomicly, here the rigid interactions are subtracted
  // note that instantaneous molecular strain derivative tensor is NOT symmetric
  // the total force is not directed along the line of the center of the molecules and the
  // contribution subtracted here introduces a torque on the molecules
  // the average pressure tensor however, is symmetric

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        com=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition;

        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          pos=Adsorbates[CurrentSystem][i].Atoms[A].Position;
          force=Adsorbates[CurrentSystem][i].Atoms[A].Force;

          StrainDerivativeRigidCorrection[CurrentSystem].ax+=(pos.x-com.x)*force.x;
          StrainDerivativeRigidCorrection[CurrentSystem].ay+=(pos.x-com.x)*force.y;
          StrainDerivativeRigidCorrection[CurrentSystem].az+=(pos.x-com.x)*force.z;

          StrainDerivativeRigidCorrection[CurrentSystem].bx+=(pos.y-com.y)*force.x;
          StrainDerivativeRigidCorrection[CurrentSystem].by+=(pos.y-com.y)*force.y;
          StrainDerivativeRigidCorrection[CurrentSystem].bz+=(pos.y-com.y)*force.z;

          StrainDerivativeRigidCorrection[CurrentSystem].cx+=(pos.z-com.z)*force.x;
          StrainDerivativeRigidCorrection[CurrentSystem].cy+=(pos.z-com.z)*force.y;
          StrainDerivativeRigidCorrection[CurrentSystem].cz+=(pos.z-com.z)*force.z;
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
        com=Cations[CurrentSystem][i].Groups[l].CenterOfMassPosition;

        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          pos=Cations[CurrentSystem][i].Atoms[A].Position;
          force=Cations[CurrentSystem][i].Atoms[A].Force;

          StrainDerivativeRigidCorrection[CurrentSystem].ax+=(pos.x-com.x)*force.x;
          StrainDerivativeRigidCorrection[CurrentSystem].ay+=(pos.x-com.x)*force.y;
          StrainDerivativeRigidCorrection[CurrentSystem].az+=(pos.x-com.x)*force.z;

          StrainDerivativeRigidCorrection[CurrentSystem].bx+=(pos.y-com.y)*force.x;
          StrainDerivativeRigidCorrection[CurrentSystem].by+=(pos.y-com.y)*force.y;
          StrainDerivativeRigidCorrection[CurrentSystem].bz+=(pos.y-com.y)*force.z;

          StrainDerivativeRigidCorrection[CurrentSystem].cx+=(pos.z-com.z)*force.x;
          StrainDerivativeRigidCorrection[CurrentSystem].cy+=(pos.z-com.z)*force.y;
          StrainDerivativeRigidCorrection[CurrentSystem].cz+=(pos.z-com.z)*force.z;
        }
      }
    }
  }

  // correct the first strain-derivative for constraints (rigid units)
  StrainDerivativeTensor[CurrentSystem].ax+=StrainDerivativeRigidCorrection[CurrentSystem].ax;
  StrainDerivativeTensor[CurrentSystem].ay+=StrainDerivativeRigidCorrection[CurrentSystem].ay;
  StrainDerivativeTensor[CurrentSystem].az+=StrainDerivativeRigidCorrection[CurrentSystem].az;
  StrainDerivativeTensor[CurrentSystem].bx+=StrainDerivativeRigidCorrection[CurrentSystem].bx;
  StrainDerivativeTensor[CurrentSystem].by+=StrainDerivativeRigidCorrection[CurrentSystem].by;
  StrainDerivativeTensor[CurrentSystem].bz+=StrainDerivativeRigidCorrection[CurrentSystem].bz;
  StrainDerivativeTensor[CurrentSystem].cx+=StrainDerivativeRigidCorrection[CurrentSystem].cx;
  StrainDerivativeTensor[CurrentSystem].cy+=StrainDerivativeRigidCorrection[CurrentSystem].cy;
  StrainDerivativeTensor[CurrentSystem].cz+=StrainDerivativeRigidCorrection[CurrentSystem].cz;

  // correct the first strain-derivative for tail-corrections
  StrainDerivativeTensor[CurrentSystem].ax+=StrainDerivativeTailCorrection[CurrentSystem];
  StrainDerivativeTensor[CurrentSystem].by+=StrainDerivativeTailCorrection[CurrentSystem];
  StrainDerivativeTensor[CurrentSystem].cz+=StrainDerivativeTailCorrection[CurrentSystem];

  // symmetrize strain-derivative (asymmetry originates from torque induced by rigid units and bond-dipoles)
  temp=0.5*(StrainDerivativeTensor[CurrentSystem].ay+StrainDerivativeTensor[CurrentSystem].bx);
  StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].bx=temp;
  temp=0.5*(StrainDerivativeTensor[CurrentSystem].az+StrainDerivativeTensor[CurrentSystem].cx);
  StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].cx=temp;
  temp=0.5*(StrainDerivativeTensor[CurrentSystem].bz+StrainDerivativeTensor[CurrentSystem].cy);
  StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cy=temp;


  // correct the second strain-derivative for tail-corrections
  if(ComputeBornTerm)
  {
    BornTerm[CurrentSystem].xxxx-=StrainDerivativeTailCorrection[CurrentSystem];
    BornTerm[CurrentSystem].yyxx-=StrainDerivativeTailCorrection[CurrentSystem];
    BornTerm[CurrentSystem].zzxx-=StrainDerivativeTailCorrection[CurrentSystem];
    BornTerm[CurrentSystem].xxyy-=StrainDerivativeTailCorrection[CurrentSystem];
    BornTerm[CurrentSystem].yyyy-=StrainDerivativeTailCorrection[CurrentSystem];
    BornTerm[CurrentSystem].zzyy-=StrainDerivativeTailCorrection[CurrentSystem];
    BornTerm[CurrentSystem].xxzz-=StrainDerivativeTailCorrection[CurrentSystem];
    BornTerm[CurrentSystem].yyzz-=StrainDerivativeTailCorrection[CurrentSystem];
    BornTerm[CurrentSystem].zzzz-=StrainDerivativeTailCorrection[CurrentSystem];
  }

  // compute the stress-tensor from the strain-derivative tensor
  ConfigurationalStressTensor[CurrentSystem].ax=-StrainDerivativeTensor[CurrentSystem].ax/Volume[CurrentSystem];
  ConfigurationalStressTensor[CurrentSystem].bx=-StrainDerivativeTensor[CurrentSystem].bx/Volume[CurrentSystem];
  ConfigurationalStressTensor[CurrentSystem].cx=-StrainDerivativeTensor[CurrentSystem].cx/Volume[CurrentSystem];

  ConfigurationalStressTensor[CurrentSystem].ay=-StrainDerivativeTensor[CurrentSystem].ay/Volume[CurrentSystem];
  ConfigurationalStressTensor[CurrentSystem].by=-StrainDerivativeTensor[CurrentSystem].by/Volume[CurrentSystem];
  ConfigurationalStressTensor[CurrentSystem].cy=-StrainDerivativeTensor[CurrentSystem].cy/Volume[CurrentSystem];

  ConfigurationalStressTensor[CurrentSystem].az=-StrainDerivativeTensor[CurrentSystem].az/Volume[CurrentSystem];
  ConfigurationalStressTensor[CurrentSystem].bz=-StrainDerivativeTensor[CurrentSystem].bz/Volume[CurrentSystem];
  ConfigurationalStressTensor[CurrentSystem].cz=-StrainDerivativeTensor[CurrentSystem].cz/Volume[CurrentSystem];

  // compute the pressure arising from interactions from the trace of the stress tensor
  Pressure[CurrentSystem]=(StressTensor[CurrentSystem].ax+StressTensor[CurrentSystem].by+StressTensor[CurrentSystem].cz)/3.0;

  KineticStressTensor[CurrentSystem]=GetKineticStressTensor();

  KineticStressTensor[CurrentSystem].ax/=Volume[CurrentSystem];
  KineticStressTensor[CurrentSystem].ay/=Volume[CurrentSystem];
  KineticStressTensor[CurrentSystem].az/=Volume[CurrentSystem];

  KineticStressTensor[CurrentSystem].bx/=Volume[CurrentSystem];
  KineticStressTensor[CurrentSystem].by/=Volume[CurrentSystem];
  KineticStressTensor[CurrentSystem].bz/=Volume[CurrentSystem];

  KineticStressTensor[CurrentSystem].cx/=Volume[CurrentSystem];
  KineticStressTensor[CurrentSystem].cy/=Volume[CurrentSystem];
  KineticStressTensor[CurrentSystem].cz/=Volume[CurrentSystem];

  StressTensor[CurrentSystem].ax=ConfigurationalStressTensor[CurrentSystem].ax+KineticStressTensor[CurrentSystem].ax;
  StressTensor[CurrentSystem].ay=ConfigurationalStressTensor[CurrentSystem].ay+KineticStressTensor[CurrentSystem].ay;
  StressTensor[CurrentSystem].az=ConfigurationalStressTensor[CurrentSystem].az+KineticStressTensor[CurrentSystem].az;
  StressTensor[CurrentSystem].bx=ConfigurationalStressTensor[CurrentSystem].bx+KineticStressTensor[CurrentSystem].bx;
  StressTensor[CurrentSystem].by=ConfigurationalStressTensor[CurrentSystem].by+KineticStressTensor[CurrentSystem].by;
  StressTensor[CurrentSystem].bz=ConfigurationalStressTensor[CurrentSystem].bz+KineticStressTensor[CurrentSystem].bz;
  StressTensor[CurrentSystem].cx=ConfigurationalStressTensor[CurrentSystem].cx+KineticStressTensor[CurrentSystem].cx;
  StressTensor[CurrentSystem].cy=ConfigurationalStressTensor[CurrentSystem].cy+KineticStressTensor[CurrentSystem].cy;
  StressTensor[CurrentSystem].cz=ConfigurationalStressTensor[CurrentSystem].cz+KineticStressTensor[CurrentSystem].cz;

  UTotal[CurrentSystem]=
         // inter
         UHostHost[CurrentSystem]+UCationCation[CurrentSystem]+UAdsorbateAdsorbate[CurrentSystem]+
         UHostCation[CurrentSystem]+UHostAdsorbate[CurrentSystem]+UAdsorbateCation[CurrentSystem]+
         // intra framework
         UHostBond[CurrentSystem]+UHostUreyBradley[CurrentSystem]+UHostBend[CurrentSystem]+
         UHostInversionBend[CurrentSystem]+UHostTorsion[CurrentSystem]+UHostImproperTorsion[CurrentSystem]+
         UHostBondBond[CurrentSystem]+UHostBondBend[CurrentSystem]+UHostBendBend[CurrentSystem]+
         UHostBondTorsion[CurrentSystem]+UHostBendTorsion[CurrentSystem]+
         // intra adsorbates
         UAdsorbateBond[CurrentSystem]+UAdsorbateUreyBradley[CurrentSystem]+UAdsorbateBend[CurrentSystem]+
         UAdsorbateInversionBend[CurrentSystem]+UAdsorbateTorsion[CurrentSystem]+UAdsorbateImproperTorsion[CurrentSystem]+
         UAdsorbateBondBond[CurrentSystem]+UAdsorbateBondBend[CurrentSystem]+UAdsorbateBendBend[CurrentSystem]+
         UAdsorbateBondTorsion[CurrentSystem]+UAdsorbateBendTorsion[CurrentSystem]+
         UAdsorbateIntraVDW[CurrentSystem]+UAdsorbateIntraChargeCharge[CurrentSystem]+
         UAdsorbateIntraChargeBondDipole[CurrentSystem]+UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+
         // intra cations
         UCationBond[CurrentSystem]+UCationUreyBradley[CurrentSystem]+UCationBend[CurrentSystem]+
         UCationInversionBend[CurrentSystem]+UCationTorsion[CurrentSystem]+UCationImproperTorsion[CurrentSystem]+
         UCationBondBond[CurrentSystem]+UCationBondBend[CurrentSystem]+UCationBendBend[CurrentSystem]+
         UCationBondTorsion[CurrentSystem]+UCationBendTorsion[CurrentSystem]+
         UCationIntraVDW[CurrentSystem]+UCationIntraChargeCharge[CurrentSystem]+
         UCationIntraChargeBondDipole[CurrentSystem]+UCationIntraBondDipoleBondDipole[CurrentSystem]+
         // polarization
         UHostPolarization[CurrentSystem]+UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem]+
         UHostBackPolarization[CurrentSystem]+UAdsorbateBackPolarization[CurrentSystem]+UCationBackPolarization[CurrentSystem]+
         // tailcorrection
         UTailCorrection[CurrentSystem]+
         // harmonic constraint energy
         UDistanceConstraints[CurrentSystem]+UAngleConstraints[CurrentSystem]+UDihedralConstraints[CurrentSystem]+
         UInversionBendConstraints[CurrentSystem]+UOutOfPlaneDistanceConstraints[CurrentSystem]+UExclusionConstraints[CurrentSystem];

}

void CalculateMolecularExcessPressure(void)
{
  int i,j,l;
  VECTOR pos,force,com;
  REAL temp;
  int Type;

  StrainDerivativeTensor[CurrentSystem].ax=StrainDerivativeTensor[CurrentSystem].bx=StrainDerivativeTensor[CurrentSystem].cx=0.0;
  StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].by=StrainDerivativeTensor[CurrentSystem].cy=0.0;
  StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cz=0.0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Force.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].Force.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].Force.z=0.0;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Force.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].Force.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].Force.z=0.0;
    }
  }

  CalculateTotalInterVDWForce();
  CalculateTotalInterChargeChargeCoulombForce();
  CalculateTotalInterChargeBondDipoleCoulombForce();
  CalculateTotalInterBondDipoleBondDipoleCoulombForce();

  // the contribution of the charges present in the system (long-range interaction)
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    EwaldFourierForce();

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    com=GetAdsorbateCenterOfMass(i);
    for(l=0;l<Components[Type].NumberOfAtoms;l++)
    {
      pos=Adsorbates[CurrentSystem][i].Atoms[l].Position;
      force=Adsorbates[CurrentSystem][i].Atoms[l].Force;

      StrainDerivativeTensor[CurrentSystem].ax+=(pos.x-com.x)*force.x;
      StrainDerivativeTensor[CurrentSystem].ay+=(pos.x-com.x)*force.y;
      StrainDerivativeTensor[CurrentSystem].az+=(pos.x-com.x)*force.z;

      StrainDerivativeTensor[CurrentSystem].bx+=(pos.y-com.y)*force.x;
      StrainDerivativeTensor[CurrentSystem].by+=(pos.y-com.y)*force.y;
      StrainDerivativeTensor[CurrentSystem].bz+=(pos.y-com.y)*force.z;

      StrainDerivativeTensor[CurrentSystem].cx+=(pos.z-com.z)*force.x;
      StrainDerivativeTensor[CurrentSystem].cy+=(pos.z-com.z)*force.y;
      StrainDerivativeTensor[CurrentSystem].cz+=(pos.z-com.z)*force.z;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    com=GetCationCenterOfMass(i);
    for(l=0;l<Components[Type].NumberOfAtoms;l++)
    {
      pos=Cations[CurrentSystem][i].Atoms[l].Position;
      force=Cations[CurrentSystem][i].Atoms[l].Force;

      StrainDerivativeTensor[CurrentSystem].ax+=(pos.x-com.x)*force.x;
      StrainDerivativeTensor[CurrentSystem].ay+=(pos.x-com.x)*force.y;
      StrainDerivativeTensor[CurrentSystem].az+=(pos.x-com.x)*force.z;

      StrainDerivativeTensor[CurrentSystem].bx+=(pos.y-com.y)*force.x;
      StrainDerivativeTensor[CurrentSystem].by+=(pos.y-com.y)*force.y;
      StrainDerivativeTensor[CurrentSystem].bz+=(pos.y-com.y)*force.z;

      StrainDerivativeTensor[CurrentSystem].cx+=(pos.z-com.z)*force.x;
      StrainDerivativeTensor[CurrentSystem].cy+=(pos.z-com.z)*force.y;
      StrainDerivativeTensor[CurrentSystem].cz+=(pos.z-com.z)*force.z;
    }
  }

  temp=0.5*(StrainDerivativeTensor[CurrentSystem].ay+StrainDerivativeTensor[CurrentSystem].bx);
  StrainDerivativeTensor[CurrentSystem].ay=StrainDerivativeTensor[CurrentSystem].bx=temp;
  temp=0.5*(StrainDerivativeTensor[CurrentSystem].az+StrainDerivativeTensor[CurrentSystem].cx);
  StrainDerivativeTensor[CurrentSystem].az=StrainDerivativeTensor[CurrentSystem].cx=temp;
  temp=0.5*(StrainDerivativeTensor[CurrentSystem].bz+StrainDerivativeTensor[CurrentSystem].cy);
  StrainDerivativeTensor[CurrentSystem].bz=StrainDerivativeTensor[CurrentSystem].cy=temp;
}

// the main routine for MD: the integrator
void Integration(void)
{
  switch(Ensemble[CurrentSystem])
  {
    case NPTPR:
    case MuPTPR:
      NoseHooverNPTPR();

      if(DegreesOfFreedomRotation[CurrentSystem]>0)
        NoseHooverNVTRotation();

      // evolve the positions a half timestep
      UpdateVelocities();

      UpdatePositionsVelocitiesNPTPR();

      // evolve the part of rigid bodies involving free rotation
      NoSquishFreeRotorOrderTwo();

      // compute the forces on all the atoms
      CalculateForce();

      // evolve the positions a half timestep
      UpdateVelocities();

      if(DegreesOfFreedomRotation[CurrentSystem]>0)
        NoseHooverNVTRotation();

      NoseHooverNPTPR();
      break;
    case NPHPR:
      NoseHooverNPHPR();

      // evolve the positions a half timestep
      UpdateVelocities();

      UpdatePositionsVelocitiesNPTPR();

      // evolve the part of rigid bodies involving free rotation
      NoSquishFreeRotorOrderTwo();

      // store the structure-factors for the Ewald-summations
      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      {
        SetupKVectors();
        EwaldEnergyIon();
      }

      // compute the forces on all the atoms
      CalculateForce();

      // evolve the positions a half timestep
      UpdateVelocities();

      NoseHooverNPHPR();
      break;
    case NPH:
      // update the velocity of the cell
      UpdateCellVelocity();

      // evolve the positions a half timestep
      UpdateVelocities();

      // evolve the positions a full timestep
      UpdatePositions();

      // evolve the part of rigid bodies involving free rotation
      NoSquishFreeRotorOrderTwo();

      // store the structure-factors for the Ewald-summations
      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      {
        SetupKVectors();
        EwaldEnergyIon();
      }

      // compute the forces on all the atoms
      CalculateForce();

      // evolve the positions a half timestep
      UpdateVelocities();

      // update the velocity of the cell
      UpdateCellVelocity();
      break;
    case NPT:
    case MuPT:
      // apply thermo/barostats for temperature, volume, and/or box-shape control
      // integrate the fictive velocities another half a timestep
      NoseHooverNPT();

      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
        NoseHooverNVTFramework();

      if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
        NoseHooverNVTAdsorbates();

      if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
        NoseHooverNVTCations();

      if(DegreesOfFreedomRotation[CurrentSystem]>0)
        NoseHooverNVTRotation();

      // update the velocity of the cell
      UpdateCellVelocity();

      // evolve the positions a half timestep
      UpdateVelocities();

      // evolve the positions a full timestep
      UpdatePositions();

      // evolve the part of rigid bodies involving free rotation
      NoSquishFreeRotorOrderTwo();

      // compute the forces on all the atoms
      CalculateForce();

      // evolve the positions a half timestep
      UpdateVelocities();

      // update the velocity of the cell
      UpdateCellVelocity();

      // apply thermo/barostats for temperature, volume, and/or box-shape control
      // integrate the fictive velocities another half a timestep
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
        NoseHooverNVTFramework();

      if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
        NoseHooverNVTAdsorbates();

      if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
        NoseHooverNVTCations();

      if(DegreesOfFreedomRotation[CurrentSystem]>0)
        NoseHooverNVTRotation();

      NoseHooverNPT();
      break;
    case NVT:
    case MuVT:
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
        NoseHooverNVTFramework();

      if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
        NoseHooverNVTAdsorbates();

      if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
        NoseHooverNVTCations();

      if(DegreesOfFreedomRotation[CurrentSystem]>0)
        NoseHooverNVTRotation();

      // evolve the positions a half timestep
      UpdateVelocities();

      // evolve the positions a full timestep
      UpdatePositions();

      // first part of bond-constraints
      RattleStageOne();

      // evolve the part of rigid bodies involving free rotation
      NoSquishFreeRotorOrderTwo();

      // compute the forces on all the atoms
      CalculateForce();

      // evolve the positions a half timestep
      UpdateVelocities();

      // second part of bond-constraints
      RattleStageTwo();

      // apply thermo/barostats for temperature, volume, and/or box-shape control
      // integrate the fictive velocities another half a timestep
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
        NoseHooverNVTFramework();

      if(DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]>0)
        NoseHooverNVTAdsorbates();

      if(DegreesOfFreedomTranslationalCations[CurrentSystem]>0)
        NoseHooverNVTCations();

      if(DegreesOfFreedomRotation[CurrentSystem]>0)
        NoseHooverNVTRotation();
      break;
    case NVE:
    default:
      // evolve the positions a half timestep
      UpdateVelocities();

      // evolve the positions a full timestep
      UpdatePositions();

      // first part of bond-constraints
      RattleStageOne();

      // evolve the part of rigid bodies involving free rotation
      NoSquishFreeRotorOrderTwo();

      // compute the forces on all the atoms
      CalculateForce();

      // evolve the positions a half timestep
      UpdateVelocities();

      // second part of bond-constraints
      RattleStageTwo();
      break;
  }

  // for rigid units: create atomic velocities from com-velocity and quaternion-velocity
  CreateCartesianVelocities();

  ComputeNoseHooverEnergySystem();

  ComputeKineticEnergySystem();

  ConservedEnergy[CurrentSystem]=UTotal[CurrentSystem]+UKinetic[CurrentSystem]+UNoseHoover[CurrentSystem];
}

int CalculateElectricField(void)
{
  int i,j,k,f1;
  REAL_MATRIX3x3 pol_factor;
  VECTOR electric_field,electric_field1,electric_field2;
  VECTOR induced_dipole;
  int Type;

  // set "static" polarization to zero
  UHostPolarization[CurrentSystem]=0.0;
  UAdsorbatePolarization[CurrentSystem]=0.0;
  UCationPolarization[CurrentSystem]=0.0;

  // set "back" polarization to zero
  UHostBackPolarization[CurrentSystem]=0.0;
  UAdsorbateBackPolarization[CurrentSystem]=0.0;
  UCationBackPolarization[CurrentSystem]=0.0;

  // set framework electric fields to zero
  if((BackPolarization)||(Framework[CurrentSystem].FrameworkModel==FLEXIBLE))
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField.x=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField.y=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField.z=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.x=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.y=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.z=0.0;
      }
    }
  }

  // set adsorbate electric fields to zero
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField.z=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }

  // set cation electric fields to zero
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField.z=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }

  // calculate static electric field
  CalculateTotalInterChargeChargeCoulombElectricFieldMC(FALSE,-1,-1);
  CalculateFrameworkChargeChargeElectricFieldMC(FALSE,-1,-1);
  ComputeStaticElectricFieldEwald(FALSE,-1,-1);

  // set the framework induces dipoles when backpolarization is used or the framework is flexible
  if((BackPolarization)||(Framework[CurrentSystem].FrameworkModel==FLEXIBLE))
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Type=Framework[CurrentSystem].Atoms[f1][i].Type;
        pol_factor=PseudoAtoms[Type].Polarization;
        electric_field=Framework[CurrentSystem].Atoms[f1][i].ElectricField;

        // Eq. 10 from T.M. Nymand and P. Linse, J. Chem. Phys., 112(14), 2000
        Framework[CurrentSystem].Atoms[f1][i].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
        Framework[CurrentSystem].Atoms[f1][i].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
        Framework[CurrentSystem].Atoms[f1][i].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;

        Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField=electric_field;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.x=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.y=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.z=0.0;
      }
    }
  }

  // set the adsorbate induces dipoles
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      pol_factor=PseudoAtoms[Type].Polarization;
      electric_field=Adsorbates[CurrentSystem][i].Atoms[j].ElectricField;

      // Eq. 10 from T.M. Nymand and P. Linse, J. Chem. Phys., 112(14), 2000
      Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
      Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
      Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;

      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField=electric_field;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }

  // set the cation induces dipoles
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][i].Atoms[j].Type;
      pol_factor=PseudoAtoms[Type].Polarization;
      electric_field=Cations[CurrentSystem][i].Atoms[j].ElectricField;

      // Eq. 10 from T.M. Nymand and P. Linse, J. Chem. Phys., 112(14), 2000
      Cations[CurrentSystem][i].Atoms[j].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
      Cations[CurrentSystem][i].Atoms[j].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
      Cations[CurrentSystem][i].Atoms[j].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;

      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField=electric_field;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }

  // initialize the polarization energies
  UHostPolarization[CurrentSystem]=0.0;
  UAdsorbatePolarization[CurrentSystem]=0.0;
  UCationPolarization[CurrentSystem]=0.0;

  // calculate the host-polarization from the induces dipoles and electric field, Eq. 13 from T.M. Nymand and P. Linse, J. Chem. Phys., 112(14), 2000
  if((BackPolarization)||(Framework[CurrentSystem].FrameworkModel==FLEXIBLE))
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          Type=Framework[CurrentSystem].Atoms[f1][i].Type;
          electric_field=Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField;
          induced_dipole=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;
          UHostPolarization[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
        }
      }
    }
  }

  // calculate the adsorbate-polarization from the induces dipoles and electric field, Eq. 13 from T.M. Nymand and P. Linse, J. Chem. Phys., 112(14), 2000
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      induced_dipole=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;
      electric_field=Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField;
      UAdsorbatePolarization[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }

  // calculate the cation-polarization from the induces dipoles and electric field, Eq. 13 from T.M. Nymand and P. Linse, J. Chem. Phys., 112(14), 2000
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][i].Atoms[j].Type;
      induced_dipole=Cations[CurrentSystem][i].Atoms[j].InducedDipole;
      electric_field=Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField;
      UCationPolarization[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }

  if(!BackPolarization) return 0;

  // iterate to convergence

  for(k=0;k<NumberOfBackPolarizationSteps;k++)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.x=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.y=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.z=0.0;
      }
    }

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
        Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
        Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
      }
    }

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Cations[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
        Cations[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
        Cations[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
      }
    }

    // compute electric field generated by the induced dipoles
    //CalculateInterElectricFieldFromInducedDipoles();
    CalculateInterElectricFieldFromInducedDipoleMC(FALSE,-1,-1);
    CalculateFrameworkMoleculeElectricFieldFromInducedDipoleMC(FALSE,-1,-1);
    ComputeElectricFieldFromInducedDipolesEwaldMC(FALSE,-1,-1);

    // computer new induced dipoles from electric field
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Type=Framework[CurrentSystem].Atoms[f1][i].Type;
        pol_factor=PseudoAtoms[Type].Polarization;
        electric_field1=Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField;
        electric_field2=Framework[CurrentSystem].Atoms[f1][i].ElectricField;
        electric_field.x=electric_field1.x+electric_field2.x;
        electric_field.y=electric_field1.y+electric_field2.y;
        electric_field.z=electric_field1.z+electric_field2.z;

        // update the new induced dipoles from the statc electric field plus the induced field
        Framework[CurrentSystem].Atoms[f1][i].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
        Framework[CurrentSystem].Atoms[f1][i].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
        Framework[CurrentSystem].Atoms[f1][i].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;
      }
    }

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        pol_factor=PseudoAtoms[Type].Polarization;
        electric_field1=Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField;
        electric_field2=Adsorbates[CurrentSystem][i].Atoms[j].ElectricField;
        electric_field.x=electric_field1.x+electric_field2.x;
        electric_field.y=electric_field1.y+electric_field2.y;
        electric_field.z=electric_field1.z+electric_field2.z;

        // update the new induced dipoles from the statc electric field plus the induced field
        Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
        Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
        Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;
      }
    }

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Type=Cations[CurrentSystem][i].Atoms[j].Type;
        pol_factor=PseudoAtoms[Type].Polarization;
        electric_field1=Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField;
        electric_field2=Cations[CurrentSystem][i].Atoms[j].ElectricField;
        electric_field.x=electric_field1.x+electric_field2.x;
        electric_field.y=electric_field1.y+electric_field2.y;
        electric_field.z=electric_field1.z+electric_field2.z;

        // update the new induced dipoles from the statc electric field plus the induced field
        Cations[CurrentSystem][i].Atoms[j].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
        Cations[CurrentSystem][i].Atoms[j].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
        Cations[CurrentSystem][i].Atoms[j].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;
      }
    }
  }

  // sum the energy
  UHostBackPolarization[CurrentSystem]=0.0;
  UAdsorbateBackPolarization[CurrentSystem]=0.0;
  UCationBackPolarization[CurrentSystem]=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Type=Framework[CurrentSystem].Atoms[f1][i].Type;
      pol_factor=PseudoAtoms[Type].Polarization;
      electric_field=Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField;
      induced_dipole=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;
      UHostBackPolarization[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      pol_factor=PseudoAtoms[Type].Polarization;
      induced_dipole=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;
      electric_field=Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField;
      UAdsorbateBackPolarization[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][i].Atoms[j].Type;
      pol_factor=PseudoAtoms[Type].Polarization;
      induced_dipole=Cations[CurrentSystem][i].Atoms[j].InducedDipole;
      electric_field=Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField;
      UCationBackPolarization[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }

  UHostBackPolarization[CurrentSystem]-=UHostPolarization[CurrentSystem];
  UAdsorbateBackPolarization[CurrentSystem]-=UAdsorbatePolarization[CurrentSystem];
  UCationBackPolarization[CurrentSystem]-=UCationPolarization[CurrentSystem];

  return 0;
}


// integrate using a leap-frog algorithm
// never use (unless for comparison or testing)
void IntegrationLeapFrog(void)
{
  int i,j,l;
  int Type,A;
  VECTOR Force,Velocity;
  VECTOR Torque,dr,op,tr,inertia,inv,oq;
  REAL Mass,tmp,eps,rnorm;
  REAL_MATRIX3x3 M;
  QUATERNION q,qn,qna,qnb;

  CurrentSystem=0;
  UAdsorbateKinetic[CurrentSystem]=0.0;

  UAdsorbateRotationalKinetic[CurrentSystem]=0.0;
  UAdsorbateTranslationalKinetic[CurrentSystem]=0.0;

  // rotational motion
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        // compute the torque
        Torque.x=0.0;
        Torque.y=0.0;
        Torque.z=0.0;
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Force=Adsorbates[CurrentSystem][i].Atoms[A].Force;
          dr.x=Adsorbates[CurrentSystem][i].Atoms[A].Position.x-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.x;
          dr.y=Adsorbates[CurrentSystem][i].Atoms[A].Position.y-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.y;
          dr.z=Adsorbates[CurrentSystem][i].Atoms[A].Position.z-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.z;
          Torque.x+=dr.y*Force.z-dr.z*Force.y;
          Torque.y+=dr.z*Force.x-dr.x*Force.z;
          Torque.z+=dr.x*Force.y-dr.y*Force.x;
        }
        BuildRotationMatrix(&M,Adsorbates[CurrentSystem][i].Groups[l].Quaternion);

        // angular velocity at time step n (first guess)
        op=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity;
        inv=Components[Type].Groups[l].InverseInertiaVector;
        inertia=Components[Type].Groups[l].InertiaVector;
        for(j=0;j<5;j++)
        {
          tr.x=(Torque.x*M.ax+Torque.y*M.bx+Torque.z*M.cx)*inv.x
               +(inertia.y-inertia.z)*op.y*op.z*inv.x;
          tr.y=(Torque.x*M.ay+Torque.y*M.by+Torque.z*M.cy)*inv.y
               +(inertia.z-inertia.x)*op.z*op.x*inv.y;
          tr.z=(Torque.x*M.az+Torque.y*M.bz+Torque.z*M.cz)*inv.z
               +(inertia.x-inertia.y)*op.x*op.y*inv.z;

          // improved angular velocity at time step n
          op.x=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.x+0.5*DeltaT*tr.x;
          op.y=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.y+0.5*DeltaT*tr.y;
          op.z=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.z+0.5*DeltaT*tr.z;
        }

        // angular velocity at time step n+1  (needed for quat algorithm)
        oq.x=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.x+1.5*DeltaT*tr.x;
        oq.y=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.y+1.5*DeltaT*tr.y;
        oq.z=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.z+1.5*DeltaT*tr.z;

        // angular velocity at time step n+1/2
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.x+=DeltaT*tr.x;
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.y+=DeltaT*tr.y;
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.z+=DeltaT*tr.z;

        // rotational kinetic energy
        UAdsorbateRotationalKinetic[CurrentSystem]+=
              0.5*(inertia.x*SQR(op.x)+inertia.y*SQR(op.y)+inertia.z*SQR(op.z));

        q=Adsorbates[CurrentSystem][i].Groups[l].Quaternion;
        qn.r=q.r+(-q.i*op.x-q.j*op.y-q.k*op.z)*0.5*DeltaT;
        qn.i=q.i+(q.r*op.x-q.k*op.y+q.j*op.z)*0.5*DeltaT;
        qn.j=q.j+(q.k*op.x+q.r*op.y-q.i*op.z)*0.5*DeltaT;
        qn.k=q.k+(-q.j*op.x+q.i*op.y+q.r*op.z)*0.5*DeltaT;

        qnb.r=qnb.i=qnb.j=qnb.k=0.0;
        do
        {
          qna.r=0.5*(-q.i*op.x-q.j*op.y-q.k*op.z)
               +0.5*(-qn.i*oq.x-qn.j*oq.y-qn.k*oq.z);
          qna.i=0.5*(q.r*op.x-q.k*op.y+q.j*op.z)
               +0.5*(qn.r*oq.x-qn.k*oq.y+qn.j*oq.z);
          qna.j=0.5*(q.k*op.x+q.r*op.y-q.i*op.z)
               +0.5*(qn.k*oq.x+qn.r*oq.y-qn.i*oq.z);
          qna.k=0.5*(-q.j*op.x+q.i*op.y+q.r*op.z)
               +0.5*(-qn.j*oq.x+qn.i*oq.y+qn.r*oq.z);

          qn.r=q.r+0.5*qna.r*DeltaT;
          qn.i=q.i+0.5*qna.i*DeltaT;
          qn.j=q.j+0.5*qna.j*DeltaT;
          qn.k=q.k+0.5*qna.k*DeltaT;

          rnorm=1.0/sqrt(SQR(qn.r)+SQR(qn.i)+SQR(qn.j)+SQR(qn.k));
          qn.r*=rnorm;
          qn.i*=rnorm;
          qn.j*=rnorm;
          qn.k*=rnorm;

          //convergence test
          eps=(SQR(qna.r-qnb.r)+SQR(qna.i-qnb.i)+SQR(qna.j-qnb.j)+SQR(qna.k-qnb.k))*SQR(DeltaT);
          qnb=qna;
        }while(eps>1e-5);

        // store new quaternions
        Adsorbates[CurrentSystem][i].Groups[l].Quaternion=qn;
      }
    }
  }

  // translational motion
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Mass=Components[Adsorbates[CurrentSystem][i].Type].Mass;
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Force.x=Force.y=Force.z=0.0;
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Force.x+=Adsorbates[CurrentSystem][i].Atoms[A].Force.x;
          Force.y+=Adsorbates[CurrentSystem][i].Atoms[A].Force.y;
          Force.z+=Adsorbates[CurrentSystem][i].Atoms[A].Force.z;
        }

        // centre of mass velocities at half-step to compute the kinetic energy
        tmp=0.5*DeltaT/Mass;
        Velocity.x=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x+tmp*Force.x;
        Velocity.y=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y+tmp*Force.y;
        Velocity.z=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z+tmp*Force.z;

        UAdsorbateTranslationalKinetic[CurrentSystem]+=0.5*Mass*(SQR(Velocity.x)+SQR(Velocity.y)+SQR(Velocity.z));

        // advance velocity by leapfrog
        tmp=DeltaT/Mass;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x+=tmp*Force.x;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y+=tmp*Force.y;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z+=tmp*Force.z;

        // advance position by leapfrog
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.x+=DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.y+=DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.z+=DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z;
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Mass=PseudoAtoms[Components[Type].Type[A]].Mass;

          Force=Adsorbates[CurrentSystem][i].Atoms[A].Force;

          // centre of mass velocities at half-step to compute the kinetic energy
          tmp=0.5*DeltaT/Mass;
          Velocity.x=Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x+tmp*Force.x;
          Velocity.y=Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y+tmp*Force.y;
          Velocity.z=Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z+tmp*Force.z;

          UAdsorbateTranslationalKinetic[CurrentSystem]+=0.5*Mass*(SQR(Velocity.x)+SQR(Velocity.y)+SQR(Velocity.z));


          // advance velocity by leapfrog
          tmp=DeltaT/Mass;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x+=tmp*Force.x;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y+=tmp*Force.y;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z+=tmp*Force.z;

          // advance position by leapfrog
          Adsorbates[CurrentSystem][i].Atoms[A].Position.x+=DeltaT*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x;
          Adsorbates[CurrentSystem][i].Atoms[A].Position.y+=DeltaT*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y;
          Adsorbates[CurrentSystem][i].Atoms[A].Position.z+=DeltaT*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z;
        }
      }
    }
  }

  UKinetic[CurrentSystem]=UAdsorbateTranslationalKinetic[CurrentSystem]+UAdsorbateRotationalKinetic[CurrentSystem];

  ConservedEnergy[CurrentSystem]=
             UTotal[CurrentSystem]+
             UKinetic[CurrentSystem]+
             UNoseHoover[CurrentSystem]+
             UNoseHooverAdsorbates[CurrentSystem]+
             UNoseHooverCations[CurrentSystem]+
             UNoseHooverFramework[CurrentSystem];

    // new atomic positions for atoms in rigid bodies - relative to c.o.m
  CreateCartesianPositions();

  CalculateForce();
}

void IntegrationLeapFrogAtomicPositions(void)
{
  int i,j;

  CurrentSystem=0;

  // advance position by leapfrog
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Position.x+=DeltaT*Adsorbates[CurrentSystem][i].Atoms[j].Velocity.x;
      Adsorbates[CurrentSystem][i].Atoms[j].Position.y+=DeltaT*Adsorbates[CurrentSystem][i].Atoms[j].Velocity.y;
      Adsorbates[CurrentSystem][i].Atoms[j].Position.z+=DeltaT*Adsorbates[CurrentSystem][i].Atoms[j].Velocity.z;
    }
  }
}

void IntegrationLeapFrogAtomicVelocities(void)
{
  int i,j;
  int Type;
  REAL Mass,tmp;
  VECTOR Force,Velocity;

  CurrentSystem=0;
  UKinetic[CurrentSystem]=0.0;

  // translational motion
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Mass=Components[Adsorbates[CurrentSystem][i].Type].Mass;
    Type=Adsorbates[CurrentSystem][i].Type;

    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[j].Type].Mass;
      Force=Adsorbates[CurrentSystem][i].Atoms[j].Force;

      // centre of mass velocities at half-step to compute the kinetic energy
      tmp=0.5*DeltaT/Mass;
      Velocity.x=Adsorbates[CurrentSystem][i].Atoms[j].Velocity.x+tmp*Force.x;
      Velocity.y=Adsorbates[CurrentSystem][i].Atoms[j].Velocity.y+tmp*Force.y;
      Velocity.z=Adsorbates[CurrentSystem][i].Atoms[j].Velocity.z+tmp*Force.z;

      UKinetic[CurrentSystem]+=0.5*Mass*(SQR(Velocity.x)+SQR(Velocity.y)+SQR(Velocity.z));

      // advance velocity by leapfrog
      tmp=DeltaT/Mass;
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.x+=tmp*Force.x;
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.y+=tmp*Force.y;
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.z+=tmp*Force.z;
    }
  }

  ConservedEnergy[CurrentSystem]=
             UTotal[CurrentSystem]+
             UKinetic[CurrentSystem]+
             UNoseHoover[CurrentSystem]+
             UNoseHooverAdsorbates[CurrentSystem]+
             UNoseHooverCations[CurrentSystem]+
             UNoseHooverFramework[CurrentSystem];
}


void IntegrationLeapFrogAtomicPositionsAndQuaternions(void)
{
  int i,j,l;
  int Type,A;
  VECTOR Force;
  VECTOR Torque,dr,op,tr,inertia,inv,oq;
  REAL Mass,eps,rnorm;
  REAL_MATRIX3x3 M;
  QUATERNION q,qn,qna,qnb;

  CurrentSystem=0;

  // rotational motion
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        // compute the torque
        Torque.x=0.0;
        Torque.y=0.0;
        Torque.z=0.0;
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Force=Adsorbates[CurrentSystem][i].Atoms[A].Force;
          dr.x=Adsorbates[CurrentSystem][i].Atoms[A].Position.x-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.x;
          dr.y=Adsorbates[CurrentSystem][i].Atoms[A].Position.y-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.y;
          dr.z=Adsorbates[CurrentSystem][i].Atoms[A].Position.z-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.z;
          Torque.x+=dr.y*Force.z-dr.z*Force.y;
          Torque.y+=dr.z*Force.x-dr.x*Force.z;
          Torque.z+=dr.x*Force.y-dr.y*Force.x;
        }
        BuildRotationMatrix(&M,Adsorbates[CurrentSystem][i].Groups[l].Quaternion);

        // angular velocity at time step n (first guess)
        op=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity;
        inv=Components[Type].Groups[l].InverseInertiaVector;
        inertia=Components[Type].Groups[l].InertiaVector;
        for(j=0;j<5;j++)
        {
          tr.x=(Torque.x*M.ax+Torque.y*M.bx+Torque.z*M.cx)*inv.x
               +(inertia.y-inertia.z)*op.y*op.z*inv.x;
          tr.y=(Torque.x*M.ay+Torque.y*M.by+Torque.z*M.cy)*inv.y
               +(inertia.z-inertia.x)*op.z*op.x*inv.y;
          tr.z=(Torque.x*M.az+Torque.y*M.bz+Torque.z*M.cz)*inv.z
               +(inertia.x-inertia.y)*op.x*op.y*inv.z;

          // improved angular velocity at time step n
          op.x=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.x-0.5*DeltaT*tr.x;
          op.y=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.y-0.5*DeltaT*tr.y;
          op.z=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.z-0.5*DeltaT*tr.z;
        }

        // angular velocity at time step n+1  (needed for quat algorithm)
        oq.x=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.x+0.5*DeltaT*tr.x;
        oq.y=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.y+0.5*DeltaT*tr.y;
        oq.z=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.z+0.5*DeltaT*tr.z;

        // rotational kinetic energy
        UAdsorbateRotationalKinetic[CurrentSystem]+=
              0.5*(inertia.x*SQR(op.x)+inertia.y*SQR(op.y)+inertia.z*SQR(op.z));

        q=Adsorbates[CurrentSystem][i].Groups[l].Quaternion;
        qn.r=q.r+(-q.i*op.x-q.j*op.y-q.k*op.z)*0.5*DeltaT;
        qn.i=q.i+(q.r*op.x-q.k*op.y+q.j*op.z)*0.5*DeltaT;
        qn.j=q.j+(q.k*op.x+q.r*op.y-q.i*op.z)*0.5*DeltaT;
        qn.k=q.k+(-q.j*op.x+q.i*op.y+q.r*op.z)*0.5*DeltaT;

        qnb.r=qnb.i=qnb.j=qnb.k=0.0;
        do
        {
          qna.r=0.5*(-q.i*op.x-q.j*op.y-q.k*op.z)
               +0.5*(-qn.i*oq.x-qn.j*oq.y-qn.k*oq.z);
          qna.i=0.5*(q.r*op.x-q.k*op.y+q.j*op.z)
               +0.5*(qn.r*oq.x-qn.k*oq.y+qn.j*oq.z);
          qna.j=0.5*(q.k*op.x+q.r*op.y-q.i*op.z)
               +0.5*(qn.k*oq.x+qn.r*oq.y-qn.i*oq.z);
          qna.k=0.5*(-q.j*op.x+q.i*op.y+q.r*op.z)
               +0.5*(-qn.j*oq.x+qn.i*oq.y+qn.r*oq.z);

          qn.r=q.r+0.5*qna.r*DeltaT;
          qn.i=q.i+0.5*qna.i*DeltaT;
          qn.j=q.j+0.5*qna.j*DeltaT;
          qn.k=q.k+0.5*qna.k*DeltaT;

          rnorm=1.0/sqrt(SQR(qn.r)+SQR(qn.i)+SQR(qn.j)+SQR(qn.k));
          qn.r*=rnorm;
          qn.i*=rnorm;
          qn.j*=rnorm;
          qn.k*=rnorm;

          //convergence test
          eps=(SQR(qna.r-qnb.r)+SQR(qna.i-qnb.i)+SQR(qna.j-qnb.j)+SQR(qna.k-qnb.k))*SQR(DeltaT);
          qnb=qna;
        }while(eps>1e-5);

        // store new quaternions
        Adsorbates[CurrentSystem][i].Groups[l].Quaternion=qn;
      }
    }
  }

  // translational motion
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
    {
      Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[0][i].Type].Mass;

      // advance velocity by leapfrog
      Framework[CurrentSystem].Atoms[0][i].Position.x+=DeltaT*Framework[CurrentSystem].Atoms[0][i].Velocity.x;
      Framework[CurrentSystem].Atoms[0][i].Position.y+=DeltaT*Framework[CurrentSystem].Atoms[0][i].Velocity.y;
      Framework[CurrentSystem].Atoms[0][i].Position.z+=DeltaT*Framework[CurrentSystem].Atoms[0][i].Velocity.z;
    }
  }
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Mass=Components[Adsorbates[CurrentSystem][i].Type].Mass;
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        // advance position by leapfrog
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.x+=DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.y+=DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.z+=DeltaT*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z;
      }
      else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];

          // advance position by leapfrog
          Adsorbates[CurrentSystem][i].Atoms[A].Position.x+=DeltaT*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x;
          Adsorbates[CurrentSystem][i].Atoms[A].Position.y+=DeltaT*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y;
          Adsorbates[CurrentSystem][i].Atoms[A].Position.z+=DeltaT*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z;
        }
      }
    }
  }
}

void IntegrationLeapFrogAtomicVelocitiesAndQuaternions(void)
{
  int i,j,l;
  int Type,A;
  VECTOR Force,Velocity;
  VECTOR Torque,dr,op,tr,inertia,inv;
  REAL Mass,tmp;
  REAL_MATRIX3x3 M;

  CurrentSystem=0;
  UKinetic[CurrentSystem]=0.0;

  // rotational motion
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        // compute the torque
        Torque.x=0.0;
        Torque.y=0.0;
        Torque.z=0.0;
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Force=Adsorbates[CurrentSystem][i].Atoms[A].Force;
          dr.x=Adsorbates[CurrentSystem][i].Atoms[A].Position.x-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.x;
          dr.y=Adsorbates[CurrentSystem][i].Atoms[A].Position.y-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.y;
          dr.z=Adsorbates[CurrentSystem][i].Atoms[A].Position.z-Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassPosition.z;
          Torque.x+=dr.y*Force.z-dr.z*Force.y;
          Torque.y+=dr.z*Force.x-dr.x*Force.z;
          Torque.z+=dr.x*Force.y-dr.y*Force.x;
        }
        BuildRotationMatrix(&M,Adsorbates[CurrentSystem][i].Groups[l].Quaternion);

        // angular velocity at time step n (first guess)
        op=Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity;
        inv=Components[Type].Groups[l].InverseInertiaVector;
        inertia=Components[Type].Groups[l].InertiaVector;
        for(j=0;j<5;j++)
        {
          tr.x=(Torque.x*M.ax+Torque.y*M.bx+Torque.z*M.cx)*inv.x
               +(inertia.y-inertia.z)*op.y*op.z*inv.x;
          tr.y=(Torque.x*M.ay+Torque.y*M.by+Torque.z*M.cy)*inv.y
               +(inertia.z-inertia.x)*op.z*op.x*inv.y;
          tr.z=(Torque.x*M.az+Torque.y*M.bz+Torque.z*M.cz)*inv.z
               +(inertia.x-inertia.y)*op.x*op.y*inv.z;
        }

        // angular velocity at time step n+1/2
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.x+=DeltaT*tr.x;
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.y+=DeltaT*tr.y;
        Adsorbates[CurrentSystem][i].Groups[l].AngularVelocity.z+=DeltaT*tr.z;

        // rotational kinetic energy
        UAdsorbateRotationalKinetic[CurrentSystem]+=
              0.5*(inertia.x*SQR(op.x)+inertia.y*SQR(op.y)+inertia.z*SQR(op.z));
      }
    }
  }

  // translational motion

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
    {
      Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[0][i].Type].Mass;
      Force=Framework[CurrentSystem].Atoms[0][i].Force;

      // centre of mass velocities at half-step to compute the kinetic energy
      tmp=0.5*DeltaT/Mass;
      Velocity.x=Framework[CurrentSystem].Atoms[0][i].Velocity.x+tmp*Force.x;
      Velocity.y=Framework[CurrentSystem].Atoms[0][i].Velocity.y+tmp*Force.y;
      Velocity.z=Framework[CurrentSystem].Atoms[0][i].Velocity.z+tmp*Force.z;

      UKinetic[CurrentSystem]+=0.5*Mass*(SQR(Velocity.x)+SQR(Velocity.y)+SQR(Velocity.z));

      // advance velocity by leapfrog
      tmp=DeltaT/Mass;
      Framework[CurrentSystem].Atoms[0][i].Velocity.x+=tmp*Force.x;
      Framework[CurrentSystem].Atoms[0][i].Velocity.y+=tmp*Force.y;
      Framework[CurrentSystem].Atoms[0][i].Velocity.z+=tmp*Force.z;
    }
  }
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Mass=Components[Adsorbates[CurrentSystem][i].Type].Mass;
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Force.x=Force.y=Force.z=0.0;
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Force.x+=Adsorbates[CurrentSystem][i].Atoms[A].Force.x;
          Force.y+=Adsorbates[CurrentSystem][i].Atoms[A].Force.y;
          Force.z+=Adsorbates[CurrentSystem][i].Atoms[A].Force.z;
        }

        // centre of mass velocities at half-step to compute the kinetic energy
        tmp=0.5*DeltaT/Mass;
        Velocity.x=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x+tmp*Force.x;
        Velocity.y=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y+tmp*Force.y;
        Velocity.z=Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z+tmp*Force.z;

        UKinetic[CurrentSystem]+=0.5*Mass*(SQR(Velocity.x)+SQR(Velocity.y)+SQR(Velocity.z));

        // advance velocity by leapfrog
        tmp=DeltaT/Mass;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x+=tmp*Force.x;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y+=tmp*Force.y;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z+=tmp*Force.z;
      }
     else
      {
        for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
        {
          A=Components[Type].Groups[l].Atoms[j];
          Mass=PseudoAtoms[Components[Type].Type[A]].Mass;

          Force=Adsorbates[CurrentSystem][i].Atoms[A].Force;

          // centre of mass velocities at half-step to compute the kinetic energy
          tmp=0.5*DeltaT/Mass;
          Velocity.x=Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x+tmp*Force.x;
          Velocity.y=Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y+tmp*Force.y;
          Velocity.z=Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z+tmp*Force.z;

          UKinetic[CurrentSystem]+=0.5*Mass*(SQR(Velocity.x)+SQR(Velocity.y)+SQR(Velocity.z));

          // advance velocity by leapfrog
          tmp=DeltaT/Mass;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x+=tmp*Force.x;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y+=tmp*Force.y;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z+=tmp*Force.z;
        }
      }
    }
  }

  ConservedEnergy[CurrentSystem]=
             UTotal[CurrentSystem]+
             UKinetic[CurrentSystem]+
             UNoseHoover[CurrentSystem]+
             UNoseHooverAdsorbates[CurrentSystem]+
             UNoseHooverCations[CurrentSystem]+
             UNoseHooverFramework[CurrentSystem];
}

REAL NumericalElectrostaticPotential(VECTOR pos)
{
  REAL U1,U2;
  REAL charge1,charge2;
  int type1,type2;
  int nr;

  nr=32;

  CurrentSystem=0;
  Adsorbates[CurrentSystem][nr].Atoms[0].Position=pos;
  type1=ReturnPseudoAtomNumber("TEST1");
  Adsorbates[CurrentSystem][nr].Atoms[0].Type=type1;
  charge1=PseudoAtoms[type1].Charge1;
  CalculateEnergy();
  U1=UTotal[0];

  type2=ReturnPseudoAtomNumber("TEST2");
  Adsorbates[CurrentSystem][nr].Atoms[0].Type=type2;
  charge2=PseudoAtoms[type2].Charge1;
  CalculateEnergy();
  U2=UTotal[0];

  return (U2-U1)/(charge2-charge1);
}


void TestElectrostaticPotential(void)
{
  int i;
  REAL U;
  VECTOR pos;

  for(i=0;i<200;i++)
  {
    pos.x=RandomNumber()*Box[CurrentSystem].ax;
    pos.y=RandomNumber()*Box[CurrentSystem].by;
    pos.z=RandomNumber()*Box[CurrentSystem].cz;

    U=0.0;
    //U+=EwaldFourierElectrostaticPotential(pos);
    U+=CalculateFrameworkElectrostaticPotential(pos);
    U+=CalculateInterChargeElectrostaticPotential(pos);

    fprintf(stderr, "(%g,%g,%g) Final U= %18.10f <-> %18.10f\n",pos.x,pos.y,pos.z,U*ENERGY_TO_KELVIN,NumericalElectrostaticPotential(pos)*ENERGY_TO_KELVIN);
  }
}
