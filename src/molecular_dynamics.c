/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'molecular_dynamics.c' is part of RASPA-2.0

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
#include <string.h>
#include <math.h>
#include <time.h>
#include "integration.h"
#include "simulation.h"
#include "molecule.h"
#include "framework_energy.h"
#include "framework.h"
#include "molecule.h"
#include "input.h"
#include "output.h"
#include "mc_moves.h"
#include "statistics.h"
#include "potentials.h"
#include "cbmc.h"
#include "grids.h"
#include "movies.h"
#include "spectra.h"
#include "sample.h"
#include "thermo_baro_stats.h"
#include "recrossing.h"
#include "molecular_dynamics.h"
#include "rigid.h"
#include "equations_of_state.h"

void MolecularDynamicsSimulation(void)
{
  int i,j,k;
  REAL ran;

  // for a crash-recovery we skip the initialization and jump to the
  // position right after the crash-file was written in the production run
  if(ContinueAfterCrash)
  {
    if(SimulationStage==POSITION_INITIALIZATION) goto ContinueAfterCrashLabel1;
    else if(SimulationStage==VELOCITY_EQUILIBRATION) goto ContinueAfterCrashLabel2;
    else goto ContinueAfterCrashLabel3;
  }

  // compute gas properties
  // here the pressure is converted to fugacities
  ComputeGasPropertiesForAllSystems();

  // open output-file
  OpenOutputFile();

  // print all the pre-simulations data
  PrintPreSimulationStatus();

  // allocate memory
  SampleRadialDistributionFunction(ALLOCATE);
  SampleProjectedLengthsDistributionFunction(ALLOCATE);
  SampleProjectedAnglesDistributionFunction(ALLOCATE);
  SampleNumberOfMoleculesHistogram(ALLOCATE);
  SamplePositionHistogram(ALLOCATE);
  SampleFreeEnergyProfile(ALLOCATE);
  SampleEndToEndDistanceHistogram(ALLOCATE);
  SampleEnergyHistogram(ALLOCATE);
  SampleFrameworkSpacingHistogram(ALLOCATE);
  SampleResidenceTimes(ALLOCATE);
  SampleDistanceHistogram(ALLOCATE);
  SampleBendAngleHistogram(ALLOCATE);
  SampleDihedralAngleHistogram(ALLOCATE);
  SampleAngleBetweenPlanesHistogram(ALLOCATE);
  SampleMoleculePropertyHistogram(ALLOCATE);
  SampleInfraRedSpectra(ALLOCATE);
  SampleMeanSquaredDisplacementOrderN(ALLOCATE);
  SampleVelocityAutoCorrelationFunctionOrderN(ALLOCATE);
  SampleRotationalVelocityAutoCorrelationFunctionOrderN(ALLOCATE);
  SampleMolecularOrientationAutoCorrelationFunctionOrderN(ALLOCATE);
  SampleBondOrientationAutoCorrelationFunctionOrderN(ALLOCATE);
  SampleMeanSquaredDisplacement(ALLOCATE);
  SampleVelocityAutoCorrelationFunction(ALLOCATE);
  SampleDensityProfile3DVTKGrid(ALLOCATE);
  SampleCationAndAdsorptionSites(ALLOCATE);
  SampleDcTSTConfigurationFiles(ALLOCATE);
  SamplePDBMovies(ALLOCATE,-1);


  // initialize
  InitializesEnergiesAllSystems();
  InitializeSmallMCStatisticsAllSystems();
  InitializeMCMovesStatisticsAllSystems();

  // compute initial energy
  CalculateTotalEnergyAllSystems();


  // Monte-Carlo initializing period to achieve a rapid equilibration of the positions
  SimulationStage=POSITION_INITIALIZATION;
  for(CurrentCycle=0;CurrentCycle<NumberOfInitializationCycles;CurrentCycle++)
  {
    if((CurrentCycle%PrintEvery)==0)
      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        PrintIntervalStatusInit(CurrentCycle,NumberOfInitializationCycles,OutputFilePtr[CurrentSystem]);

    for(j=0;j<NumberOfSystems*NumberOfComponents;j++)
    {
      // choose a random system
      CurrentSystem=(int)(RandomNumber()*(REAL)NumberOfSystems);

      for(k=0;k<MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[CurrentSystem]);k++)
      {
        // choose a random component
        CurrentComponent=(int)(RandomNumber()*(REAL)NumberOfComponents);

        // choose any of the MC moves randomly with the selected probability
        ran=RandomNumber();

        if(ran<Components[CurrentComponent].ProbabilityTranslationMove)
          TranslationMove();
        else if(ran<Components[CurrentComponent].ProbabilityRandomTranslationMove)
          RandomTranslationMove();
        else if(ran<Components[CurrentComponent].ProbabilityRotationMove)
          RotationMove();
        else if(ran<Components[CurrentComponent].ProbabilityPartialReinsertionMove)
          PartialReinsertionMove();
        else if(ran<Components[CurrentComponent].ProbabilityReinsertionMove)
          ReinsertionMove();
        else if(ran<Components[CurrentComponent].ProbabilityReinsertionInPlaceMove)
          ReinsertionInPlaceMove();
        else if(ran<Components[CurrentComponent].ProbabilityReinsertionInPlaneMove)
          ReinsertionInPlaneMove();
        else if(ran<Components[CurrentComponent].ProbabilityIdentityChangeMove)
          IdentityChangeMove();
        else if(ran<Components[CurrentComponent].ProbabilitySwapMove)
        {
          if(RandomNumber()<0.5) SwapAddMove();
          else SwapRemoveMove();
        }
        else if(ran<Components[CurrentComponent].ProbabilityCFSwapLambdaMove)
          CFSwapLambaMove();
        else if(ran<Components[CurrentComponent].ProbabilityCBCFSwapLambdaMove)
          CBCFSwapLambaMove();
        else if(ran<Components[CurrentComponent].ProbabilityWidomMove)
          ;
        else if(ran<Components[CurrentComponent].ProbabilitySurfaceAreaMove)
          ;
        else if(ran<Components[CurrentComponent].ProbabilityGibbsChangeMove)
          GibbsParticleTransferMove();
        else if(ran<Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove)
          GibbsIdentityChangeMove();
        else if(ran<Components[CurrentComponent].ProbabilityCFGibbsChangeMove)
          CFGibbsParticleTransferMove();
        else if(ran<Components[CurrentComponent].ProbabilityCBCFGibbsChangeMove)
          CBCFGibbsParticleTransferMove();
        else if(ran<Components[CurrentComponent].ProbabilityParallelTemperingMove)
          ParallelTemperingMove();
        else if(ran<Components[CurrentComponent].ProbabilityHyperParallelTemperingMove)
          HyperParallelTemperingMove();
        else if(ran<Components[CurrentComponent].ProbabilityParallelMolFractionMove)
          ParallelMolFractionMove();
        else if(ran<Components[CurrentComponent].ProbabilityChiralInversionMove)
          ChiralInversionMove();
        else if(ran<Components[CurrentComponent].ProbabilityHybridNVEMove)
          HybridNVEMove();
        else if(ran<Components[CurrentComponent].ProbabilityHybridNPHMove)
          HybridNPHMove();
        else if(ran<Components[CurrentComponent].ProbabilityHybridNPHPRMove)
          HybridNPHPRMove();
        else if(ran<Components[CurrentComponent].ProbabilityVolumeChangeMove)
          VolumeMove();
        else if(ran<Components[CurrentComponent].ProbabilityBoxShapeChangeMove)
          BoxShapeChangeMove();
        else if(ran<Components[CurrentComponent].ProbabilityGibbsVolumeChangeMove)
          GibbsVolumeMove();
        else if(ran<Components[CurrentComponent].ProbabilityFrameworkChangeMove)
          FrameworkChangeMove();
        else if(ran<Components[CurrentComponent].ProbabilityFrameworkShiftMove)
          FrameworkShiftMove();
      }
    }

    if(CurrentCycle%PrintEvery==0)
    {
      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      {
        OptimizeVolumeChangeAcceptence();
        OptimizeTranslationAcceptence();
        if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
        {
          OptimizeFrameworkChangeAcceptence();
          OptimizeFrameworkShiftAcceptence();
        }
        RescaleMaximumRotationAnglesSmallMC();
      }
    }

    if((CurrentCycle>0)&&(WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
      WriteBinaryRestartFiles();

    ContinueAfterCrashLabel1:;
  }

  // initialize
  InitializesEnergiesAllSystems();
  InitializeSmallMCStatisticsAllSystems();
  InitializeMCMovesStatisticsAllSystems();

  SimulationStage=VELOCITY_SCALING;
  if(ReinitializeVelocities)
  {
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      InitializeAdsorbateVelocities();
      InitializeCationVelocities();
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
        InitializeFrameworkVelocities();

      RemoveVelocityDrift();
      //AdjustSystemAngularRotationToZero();
    }
  }

  //InitializeNoseHooverCurrentSystem();
  InitializeNoseHooverAllSystems();

  // compute initial energy
  InitializeForcesAllSystems();

  // Molecular-Dynamics initializing period to remove excessive kinetic energy due to
  // poor equilibration of the positions (use with caution)
  // Velocity-scaling scales the velocities at each time step to the desired temperature
/*
  for(CurrentCycle=0;CurrentCycle<NumberOfVelocityScalingCycles;CurrentCycle++)
  {
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      // scake the velocity to the exact temperature
      AdjustVelocitiesToTemperature();

      // regularly output system status and restart files
      if(CurrentCycle%PrintEvery==0)
      {
        PrintIntervalStatusEquilibration(CurrentCycle,NumberOfEquilibrationCycles,OutputFilePtr[CurrentSystem]);
        PrintRestartFile();
      }
      // evolve the system a full time-step
      Integration();
    }
  }
*/

  // set the current ensemble to the initialization ensemble
  for(i=0;i<NumberOfSystems;i++)
    Ensemble[i]=InitEnsemble[i];

  InitializesEnergyAveragesAllSystems();

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    ReferenceEnergy[CurrentSystem]=ConservedEnergy[CurrentSystem];
    Drift[CurrentSystem]=0.0;
  }

  // Molecular-Dynamics initializing period to achieve a rapid equilibration of the velocities
  SimulationStage=VELOCITY_EQUILIBRATION;
  for(CurrentCycle=0;CurrentCycle<NumberOfEquilibrationCycles;CurrentCycle++)
  {
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      // regularly output system status and restart files
      if(CurrentCycle%PrintEvery==0)
      {
        PrintIntervalStatusEquilibration(CurrentCycle,NumberOfEquilibrationCycles,OutputFilePtr[CurrentSystem]);
        PrintRestartFile();
      }

      // particle insertion/deletion move for the osmotic-ensemble(MuPT), anisotropic osmotic-ensemble (MuPTPR) and MuVT ensemble
      if((Ensemble[CurrentSystem]==MuPT)||(Ensemble[CurrentSystem]==MuPTPR)||(Ensemble[CurrentSystem]==MuVT))
      {
        CurrentComponent=(int)(RandomNumber()*NumberOfComponents);
        if((CurrentCycle%(Components[CurrentComponent].SwapEvery)==0)&&
           (Components[CurrentComponent].ProbabilitySwapMove>0.0))
        {
          if(RandomNumber()<0.5) SwapAddMove();
          else SwapRemoveMove();
        }
      }

      // evolve the system a full time-step
      Integration();

      // prevent heating of the core-shells
      AdjustCoreShellVelocities();

      // update the current energy-drift
      Drift[CurrentSystem]+=fabs((ConservedEnergy[CurrentSystem]-ReferenceEnergy[CurrentSystem])/ReferenceEnergy[CurrentSystem]);
    }

    if((CurrentCycle>0)&&(WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
      WriteBinaryRestartFiles();

    ContinueAfterCrashLabel2:;
  }


  // initialize sampling-routines at the start of the production run
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    Ensemble[CurrentSystem]=RunEnsemble[CurrentSystem];

    ReferenceEnergy[CurrentSystem]=ConservedEnergy[CurrentSystem];
    Drift[CurrentSystem]=0.0;

    SampleRadialDistributionFunction(INITIALIZE);
    SampleProjectedLengthsDistributionFunction(INITIALIZE);
    SampleProjectedAnglesDistributionFunction(INITIALIZE);
    SampleNumberOfMoleculesHistogram(INITIALIZE);
    SamplePositionHistogram(INITIALIZE);
    SampleFreeEnergyProfile(INITIALIZE);
    SampleEndToEndDistanceHistogram(INITIALIZE);
    SampleEnergyHistogram(INITIALIZE);
    SampleFrameworkSpacingHistogram(INITIALIZE);
    SampleResidenceTimes(INITIALIZE);
    SampleDistanceHistogram(INITIALIZE);
    SampleBendAngleHistogram(INITIALIZE);
    SampleDihedralAngleHistogram(INITIALIZE);
    SampleAngleBetweenPlanesHistogram(INITIALIZE);
    SampleMoleculePropertyHistogram(INITIALIZE);
    SampleInfraRedSpectra(INITIALIZE);
    SampleMeanSquaredDisplacementOrderN(INITIALIZE);
    SampleVelocityAutoCorrelationFunctionOrderN(INITIALIZE);
    SampleRotationalVelocityAutoCorrelationFunctionOrderN(INITIALIZE);
    SampleMolecularOrientationAutoCorrelationFunctionOrderN(INITIALIZE);
    SampleBondOrientationAutoCorrelationFunctionOrderN(INITIALIZE);
    SampleMeanSquaredDisplacement(INITIALIZE);
    SampleVelocityAutoCorrelationFunction(INITIALIZE);
    SampleDensityProfile3DVTKGrid(INITIALIZE);
    SampleCationAndAdsorptionSites(INITIALIZE);
    SampleDcTSTConfigurationFiles(INITIALIZE);
    SamplePDBMovies(INITIALIZE,-1);
  }

  // Molecular-Dynamics production run
  // loop over the amount of production cycles (MD integration steps)
  SimulationStage=PRODUCTION;
  for(CurrentCycle=0;CurrentCycle<NumberOfCycles;CurrentCycle++)
  {
    // detect erroneous chirality changes
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      CheckChiralityMolecules();

    // loop over all the systems and handle one by one
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      // update all the average energies
      UpdateEnergyAveragesCurrentSystem();

      if(CurrentCycle%PrintPropertiesEvery==0)
        PrintPropertyStatus(CurrentCycle,NumberOfCycles,OutputFilePtr[CurrentSystem]);

      if(CurrentCycle%PrintEvery==0)
      {
        PrintIntervalStatusProduction(CurrentCycle,NumberOfCycles,OutputFilePtr[CurrentSystem]);
        PrintRestartFile();
      }

      // particle insertion/deletion move for the osmotic-ensemble(MuPT), anisotropic osmotic-ensemble (MuPTPR) and MuVT ensemble
      if((Ensemble[CurrentSystem]==MuPT)||(Ensemble[CurrentSystem]==MuPTPR)||(Ensemble[CurrentSystem]==MuVT))
      {
        CurrentComponent=(int)(RandomNumber()*NumberOfComponents);
        if((CurrentCycle%(Components[CurrentComponent].SwapEvery)==0)&&
           (Components[CurrentComponent].FractionOfSwapMove>0.0))
        {
          if(RandomNumber()<0.5) SwapAddMove();
          else SwapRemoveMove();
        }
      }

      for(CurrentComponent=0;CurrentComponent<NumberOfComponents;CurrentComponent++)
      {
        if(Components[CurrentComponent].FractionOfWidomMove>0.0)
          WidomMove();
      }

      // evolve the system a full time step
      Integration();

      // update the current energy-drift
      Drift[CurrentSystem]+=fabs((ConservedEnergy[CurrentSystem]-ReferenceEnergy[CurrentSystem])/
            ReferenceEnergy[CurrentSystem]);

      if(Drift[CurrentSystem]/(CurrentCycle+1.0)>1e2)
      {
        fprintf(OutputFilePtr[CurrentSystem],"\n\nERROR: unstable integration (energy drift %g > 1e2, simulation stopped)\n\n",Drift[CurrentSystem]/(CurrentCycle+1.0));
        fprintf(OutputFilePtr[CurrentSystem],"Check that: 1) all interaction are defined properly\n");
        fprintf(OutputFilePtr[CurrentSystem],"            2) the time step is not too large\n");
        fprintf(OutputFilePtr[CurrentSystem],"            3) the system is equilibrated\n\n");
        fflush(OutputFilePtr[CurrentSystem]);
        exit(0);
      }
    }

    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      SampleRadialDistributionFunction(SAMPLE);
      SampleProjectedLengthsDistributionFunction(SAMPLE);
      SampleProjectedAnglesDistributionFunction(SAMPLE);
      SampleNumberOfMoleculesHistogram(SAMPLE);
      SamplePositionHistogram(SAMPLE);
      SampleFreeEnergyProfile(SAMPLE);
      SampleEndToEndDistanceHistogram(SAMPLE);
      SampleEnergyHistogram(SAMPLE);
      SampleFrameworkSpacingHistogram(SAMPLE);
      SampleResidenceTimes(SAMPLE);
      SampleDistanceHistogram(SAMPLE);
      SampleBendAngleHistogram(SAMPLE);
      SampleDihedralAngleHistogram(SAMPLE);
      SampleAngleBetweenPlanesHistogram(SAMPLE);
      SampleMoleculePropertyHistogram(SAMPLE);
      SampleInfraRedSpectra(SAMPLE);
      SampleMeanSquaredDisplacementOrderN(SAMPLE);
      SampleVelocityAutoCorrelationFunctionOrderN(SAMPLE);
      SampleRotationalVelocityAutoCorrelationFunctionOrderN(SAMPLE);
      SampleMolecularOrientationAutoCorrelationFunctionOrderN(SAMPLE);
      SampleBondOrientationAutoCorrelationFunctionOrderN(SAMPLE);
      SampleMeanSquaredDisplacement(SAMPLE);
      SampleVelocityAutoCorrelationFunction(SAMPLE);
      SampleDensityProfile3DVTKGrid(SAMPLE);
      SampleCationAndAdsorptionSites(SAMPLE);
      SampleDcTSTConfigurationFiles(SAMPLE);
      SamplePDBMovies(SAMPLE,-1);
    }

    if((CurrentCycle>0)&&(WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
      WriteBinaryRestartFiles();

    ContinueAfterCrashLabel3:;

    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      // regulary output sampling function
      SampleRadialDistributionFunction(PRINT);
      SampleProjectedLengthsDistributionFunction(PRINT);
      SampleProjectedAnglesDistributionFunction(PRINT);
      SampleNumberOfMoleculesHistogram(PRINT);
      SamplePositionHistogram(PRINT);
      SampleFreeEnergyProfile(PRINT);
      SampleEndToEndDistanceHistogram(PRINT);
      SampleEnergyHistogram(PRINT);
      SampleFrameworkSpacingHistogram(PRINT);
      SampleResidenceTimes(PRINT);
      SampleDistanceHistogram(PRINT);
      SampleBendAngleHistogram(PRINT);
      SampleDihedralAngleHistogram(PRINT);
      SampleAngleBetweenPlanesHistogram(PRINT);
      SampleMoleculePropertyHistogram(PRINT);
      SampleInfraRedSpectra(PRINT);
      SampleMeanSquaredDisplacementOrderN(PRINT);
      SampleVelocityAutoCorrelationFunctionOrderN(PRINT);
      SampleRotationalVelocityAutoCorrelationFunctionOrderN(PRINT);
      SampleMolecularOrientationAutoCorrelationFunctionOrderN(PRINT);
      SampleBondOrientationAutoCorrelationFunctionOrderN(PRINT);
      SampleMeanSquaredDisplacement(PRINT);
      SampleVelocityAutoCorrelationFunction(PRINT);
      SampleDensityProfile3DVTKGrid(PRINT);
      SampleCationAndAdsorptionSites(PRINT);
      SampleDcTSTConfigurationFiles(PRINT);
      if(Movies[CurrentSystem]&&(CurrentCycle%WriteMoviesEvery[CurrentSystem]==0))
        SamplePDBMovies(PRINT,-1);
    }
  }

  if(WriteBinaryRestartFileEvery>0)
    WriteBinaryRestartFiles();

  // finalize and clean up
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    SampleRadialDistributionFunction(FINALIZE);
    SampleProjectedLengthsDistributionFunction(FINALIZE);
    SampleProjectedAnglesDistributionFunction(FINALIZE);
    SampleNumberOfMoleculesHistogram(FINALIZE);
    SamplePositionHistogram(FINALIZE);
    SampleFreeEnergyProfile(FINALIZE);
    SampleEndToEndDistanceHistogram(FINALIZE);
    SampleEnergyHistogram(FINALIZE);
    SampleFrameworkSpacingHistogram(FINALIZE);
    SampleResidenceTimes(FINALIZE);
    SampleDistanceHistogram(FINALIZE);
    SampleBendAngleHistogram(FINALIZE);
    SampleDihedralAngleHistogram(FINALIZE);
    SampleAngleBetweenPlanesHistogram(FINALIZE);
    SampleMoleculePropertyHistogram(FINALIZE);
    SampleInfraRedSpectra(FINALIZE);
    SampleMeanSquaredDisplacementOrderN(FINALIZE);
    SampleVelocityAutoCorrelationFunctionOrderN(FINALIZE);
    SampleRotationalVelocityAutoCorrelationFunctionOrderN(FINALIZE);
    SampleMolecularOrientationAutoCorrelationFunctionOrderN(FINALIZE);
    SampleBondOrientationAutoCorrelationFunctionOrderN(FINALIZE);
    SampleMeanSquaredDisplacement(FINALIZE);
    SampleVelocityAutoCorrelationFunction(FINALIZE);
    SampleDensityProfile3DVTKGrid(FINALIZE);
    SampleCationAndAdsorptionSites(FINALIZE);
    SampleDcTSTConfigurationFiles(FINALIZE);
    SamplePDBMovies(FINALIZE,-1);
  }

  // print post-status
  PrintPostSimulationStatus();
  CloseOutputFile();
}
