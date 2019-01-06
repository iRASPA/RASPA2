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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "simulation.h"
#include "molecule.h"
#include "framework_energy.h"
#include "framework.h"
#include "utils.h"
#include "molecule.h"
#include "input.h"
#include "output.h"
#include "mc_moves.h"
#include "statistics.h"
#include "potentials.h"
#include "spacegroup.h"
#include "cbmc.h"
#include "grids.h"
#include "movies.h"
#include "sample.h"
#include "recrossing.h"
#include "monte_carlo.h"
#include "integration.h"
#include "thermo_baro_stats.h"
#include "equations_of_state.h"
#include "ewald.h"

void DebugEnergyStatus(void);

/*********************************************************************************************************
 * Name       | MonteCarloSimulation                                                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Main routine to do a Monte-Carlo (MC) simulation.                                        *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void MonteCarloSimulation(void)
{
  int i,j;
  int NumberOfParticleMoves;
  int NumberOfSteps;
  int SelectedSystem;
  double cpu_start,cpu_end;
  double cpu_before,cpu_after;
  REAL ran;


  // for a crash-recovery we skip the initialization/equilibration and jump to the
  // position right after the crash-file was written in the production run
  if(ContinueAfterCrash)
  {
    if(SimulationStage==POSITION_INITIALIZATION) goto ContinueAfterCrashLabel1;
    else if(SimulationStage==CF_WANG_LANDAU_EQUILIBRATION) goto ContinueAfterCrashLabel2;
    else if(SimulationStage==PRODUCTION) goto ContinueAfterCrashLabel3;
    else goto ContinueAfterCrashLabel4;
  }

  // allocate memory for sampling routines
  SampleRadialDistributionFunction(ALLOCATE);
  SampleProjectedLengthsDistributionFunction(ALLOCATE);
  SampleProjectedAnglesDistributionFunction(ALLOCATE);
  SampleNumberOfMoleculesHistogram(ALLOCATE);
  SamplePositionHistogram(ALLOCATE);
  SampleFreeEnergyProfile(ALLOCATE);
  SamplePoreSizeDistribution(ALLOCATE);
  SampleEndToEndDistanceHistogram(ALLOCATE);
  SampleEnergyHistogram(ALLOCATE);
  SampleThermoDynamicsFactor(ALLOCATE);
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
  SampleCOMDensityProfile3DVTKGrid(ALLOCATE); // testing
  SampleCationAndAdsorptionSites(ALLOCATE);
  SampleDcTSTConfigurationFiles(ALLOCATE);
  SamplePDBMovies(ALLOCATE,-1);


  // loop over all the pressures of the isotherm
  for(CurrentIsothermPressure=0;CurrentIsothermPressure<NumberOfIsothermPressures;CurrentIsothermPressure++)
  {
    // compute gas properties
    // here the pressure is converted to fugacities
    ComputeGasPropertiesForAllSystems();

    // open output-file for systems
    OpenOutputFile();

    // print simulation settings to the output-file
    PrintPreSimulationStatus();

    // initialize the energies and compute the total energies for all systems
    InitializesEnergiesAllSystems();
    InitializesEnergyAveragesAllSystems();

    InitializeSmallMCStatisticsAllSystems();
    InitializeMCMovesStatisticsAllSystems();

    // compute total energy for all systems
    CalculateTotalEnergyAllSystems();

    CFWangLandauIteration(INITIALIZE);
    CFRXMCWangLandauIteration(INITIALIZE);

    // initialization to reach equilibration of positions (no averages are computed yet)
    SimulationStage=POSITION_INITIALIZATION;
    for(CurrentCycle=0;CurrentCycle<NumberOfInitializationCycles;CurrentCycle++)
    {
      if((WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
        WriteBinaryRestartFiles();

      // a label to jump to for a restart, everything before here is skipped
      // this is the point where the previous binary restart file was written
      ContinueAfterCrashLabel1: ;

      cpu_start=get_cpu_time();

      // detect erroneous chirality changes
      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        CheckChiralityMolecules();

      // Print at 'PrintEvery' intervals the status and a restart-file
      if((CurrentCycle%PrintEvery)==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          PrintIntervalStatusInit(CurrentCycle,NumberOfInitializationCycles,OutputFilePtr[CurrentSystem]);
          PrintRestartFile();
        }
      }

      for(i=0;i<NumberOfSystems;i++)
      {
        // choose a random system
        SelectedSystem=(int)(RandomNumber()*(REAL)NumberOfSystems);

        //NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[SelectedSystem]+NumberOfCationMolecules[SelectedSystem]);
        NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[0]+NumberOfCationMolecules[0]);
        for(j=0;j<NumberOfSystems;j++)
        {
          NumberOfParticleMoves=MAX2(NumberOfParticleMoves,NumberOfAdsorbateMolecules[j]+NumberOfCationMolecules[j]);
        }
        NumberOfSteps=NumberOfParticleMoves*(NumberOfComponents==0?1:NumberOfComponents);

        for(j=0;j<NumberOfSteps;j++)
        {
          // set the selected system
          CurrentSystem=SelectedSystem;

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
          else if(ran<Components[CurrentComponent].ProbabilityRandomRotationMove)
            RandomRotationMove();
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
          else if(ran<Components[CurrentComponent].ProbabilityCFWidomLambdaMove)
            CFWidomLambaMove();
          else if(ran<Components[CurrentComponent].ProbabilityGibbsWidomMove)
            ;
          else if(ran<Components[CurrentComponent].ProbabilityCFWidomLambdaMove)
            CFWidomLambaMove();
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
          else if(ran<Components[CurrentComponent].ProbabilityExchangeFractionalParticleMove)
            ExchangeFractionalParticleMove();
          else if(ran<Components[CurrentComponent].ProbabilityCFGibbsSwapFractionalMoleculeToOtherBoxMove)
            CFGibbsSwapFractionalMoleculeToOtherBoxMove();
          else if(ran<Components[CurrentComponent].ProbabilityCFGibbsLambdaChangeMove)
            CFGibbsLambdaChangeMove();
          else if(ran<Components[CurrentComponent].ProbabilityCFGibbsFractionalToIntegerMove)
            CFGibbsFractionalToIntegerMove();
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
          else if(ran<Components[CurrentComponent].ProbabilityCFCRXMCLambdaChangeMove)
            CFCRXMCLambdaChangeMove();

          #ifdef DEBUG
            DebugEnergyStatus();
          #endif
        }
      }

      if(CurrentCycle%OptimizeAcceptenceEvery==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          OptimizeVolumeChangeAcceptence();
          OptimizeGibbsVolumeChangeAcceptence();
          OptimizeTranslationAcceptence();
          OptimizeRotationAcceptence();
          OptimizeFrameworkChangeAcceptence();
          OptimizeFrameworkShiftAcceptence();
          OptimizeCFLambdaAcceptence();
          OptimizeCFGibbsLambdaAcceptence();
          OptimizeCBCFLambdaChangeAcceptence();
          OptimizeCBCFGibbsLambdaChangeAcceptence();
          OptimizeRXMCLambdaChangeAcceptence();
          RescaleMaximumRotationAnglesSmallMC();
          OptimizeCFGibbsLambdaChangeAcceptence();
        }
      }
      cpu_end=get_cpu_time();
      CpuTotal+=(cpu_end-cpu_start);
      CpuTimeInitialization+=(cpu_end-cpu_start);

    }


    // initialize the energies and compute the total energies for all systems
    InitializesEnergiesAllSystems();
    InitializesEnergyAveragesAllSystems();

    InitializeSmallMCStatisticsAllSystems();
    InitializeMCMovesStatisticsAllSystems();

    // compute total energy for all systems
    CalculateTotalEnergyAllSystems();

    if(NumberOfEquilibrationCycles>0)
    {
      CFWangLandauIteration(INITIALIZE);
      CFRXMCWangLandauIteration(INITIALIZE);

      // initialization to reach equilibration of positions (no averages are computed yet)
      SimulationStage=CF_WANG_LANDAU_EQUILIBRATION;
      for(CurrentCycle=0;CurrentCycle<NumberOfEquilibrationCycles;CurrentCycle++)
      {
        if((WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
          WriteBinaryRestartFiles();

        // a label to jump to for a restart, everything before here is skipped
        // this is the point where the previous binary restart file was written
        ContinueAfterCrashLabel2: ;

        cpu_start=get_cpu_time();

        // detect erroneous chirality changes
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
          CheckChiralityMolecules();

        // Print at 'PrintEvery' intervals the status and a restart-file
        if((CurrentCycle%PrintEvery)==0)
        {
          for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
          {
            PrintIntervalStatusEquilibration(CurrentCycle,NumberOfEquilibrationCycles,OutputFilePtr[CurrentSystem]);
            PrintRestartFile();
          }
        }

        // select MC moves
        for(i=0;i<NumberOfSystems;i++)
        {
          // choose a random system
          SelectedSystem=(int)(RandomNumber()*(REAL)NumberOfSystems);

          NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[0]+NumberOfCationMolecules[0]);
          for(j=0;j<NumberOfSystems;j++)
          {
            NumberOfParticleMoves=MAX2(NumberOfParticleMoves,NumberOfAdsorbateMolecules[j]+NumberOfCationMolecules[j]);
          }
          NumberOfSteps=NumberOfParticleMoves*(NumberOfComponents==0?1:NumberOfComponents);

          for(j=0;j<NumberOfSteps;j++)
          {
            // set the selected system
            CurrentSystem=SelectedSystem;

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
            else if(ran<Components[CurrentComponent].ProbabilityRandomRotationMove)
              RandomRotationMove();
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
            {
              CFSwapLambaMove();
              CFWangLandauIteration(SAMPLE);
            }
            else if(ran<Components[CurrentComponent].ProbabilityCBCFSwapLambdaMove)
            {
              CFWangLandauIteration(SAMPLE);
              CBCFSwapLambaMove();
            }
            else if(ran<Components[CurrentComponent].ProbabilityWidomMove)
              ;
            else if(ran<Components[CurrentComponent].ProbabilityCFWidomLambdaMove)
            {
              CFWidomLambaMove();
              CFWangLandauIteration(SAMPLE);
            }
            else if(ran<Components[CurrentComponent].ProbabilityGibbsWidomMove)
              ;
            else if(ran<Components[CurrentComponent].ProbabilitySurfaceAreaMove)
              ;
            else if(ran<Components[CurrentComponent].ProbabilityGibbsChangeMove)
              GibbsParticleTransferMove();
            else if(ran<Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove)
              GibbsIdentityChangeMove();
            else if(ran<Components[CurrentComponent].ProbabilityCFGibbsChangeMove)
            {
              CFWangLandauIteration(SAMPLE);
              CFGibbsParticleTransferMove();
            }
            else if(ran<Components[CurrentComponent].ProbabilityCBCFGibbsChangeMove)
            {
              CFWangLandauIteration(SAMPLE);
              CBCFGibbsParticleTransferMove();
            }
            else if(ran<Components[CurrentComponent].ProbabilityExchangeFractionalParticleMove)
              ExchangeFractionalParticleMove();
            else if(ran<Components[CurrentComponent].ProbabilityCFGibbsSwapFractionalMoleculeToOtherBoxMove)
            {
              CFGibbsSwapFractionalMoleculeToOtherBoxMove();
            }
            else if(ran<Components[CurrentComponent].ProbabilityCFGibbsLambdaChangeMove)
            {
              CFGibbsLambdaChangeMove();
              CFWangLandauIteration(SAMPLE);
            }
            else if(ran<Components[CurrentComponent].ProbabilityCFGibbsFractionalToIntegerMove)
            {
              CFGibbsFractionalToIntegerMove();
            }
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
            else if(ran<Components[CurrentComponent].ProbabilityCFCRXMCLambdaChangeMove)
            {
              CFCRXMCLambdaChangeMove();
              CFRXMCWangLandauIteration(SAMPLE);
            }


            #ifdef DEBUG
              DebugEnergyStatus();
            #endif
          }
        }

        if(CurrentCycle%OptimizeAcceptenceEvery==0)
        {
          for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
          {
            OptimizeVolumeChangeAcceptence();
            OptimizeGibbsVolumeChangeAcceptence();
            OptimizeTranslationAcceptence();
            OptimizeRotationAcceptence();
            OptimizeFrameworkChangeAcceptence();
            OptimizeFrameworkShiftAcceptence();
            OptimizeCFLambdaAcceptence();
            OptimizeCFGibbsLambdaAcceptence();
            OptimizeCBCFLambdaChangeAcceptence();
            OptimizeCBCFGibbsLambdaChangeAcceptence();
            OptimizeRXMCLambdaChangeAcceptence();
            RescaleMaximumRotationAnglesSmallMC();
            OptimizeCFWidomAcceptence();
            OptimizeCFGibbsLambdaChangeAcceptence();
          }
        }

        if(CurrentCycle%CFWangLandauEvery==0)
        {
          CFWangLandauIteration(PRINT);
          CFRXMCWangLandauIteration(PRINT);
        }

        cpu_end=get_cpu_time();
        CpuTotal+=(cpu_end-cpu_start);
        CpuTimeEquilibration+=(cpu_end-cpu_start);
      }

      CFWangLandauIteration(FINALIZE);
      CFRXMCWangLandauIteration(FINALIZE);

      // after equilibration, recompute all the energies
      InitializesEnergiesAllSystems();
      InitializesEnergyAveragesAllSystems();

      InitializeSmallMCStatisticsAllSystems();
      InitializeMCMovesStatisticsAllSystems();

      CalculateTotalEnergyAllSystems();
    }

    // initialize sampling-routines at the start of the production run
    SampleRadialDistributionFunction(INITIALIZE);
    SampleProjectedLengthsDistributionFunction(INITIALIZE);
    SampleProjectedAnglesDistributionFunction(INITIALIZE);
    SampleNumberOfMoleculesHistogram(INITIALIZE);
    SamplePositionHistogram(INITIALIZE);
    SampleFreeEnergyProfile(INITIALIZE);
    SamplePoreSizeDistribution(INITIALIZE);
    SampleEndToEndDistanceHistogram(INITIALIZE);
    SampleEnergyHistogram(INITIALIZE);
    SampleThermoDynamicsFactor(INITIALIZE);
    SampleFrameworkSpacingHistogram(INITIALIZE);
    SampleResidenceTimes(INITIALIZE);
    SampleDistanceHistogram(INITIALIZE);
    SampleBendAngleHistogram(INITIALIZE);
    SampleDihedralAngleHistogram(INITIALIZE);
    SampleAngleBetweenPlanesHistogram(INITIALIZE);
    SampleMoleculePropertyHistogram(INITIALIZE);
    SampleDensityProfile3DVTKGrid(INITIALIZE);
    SampleCOMDensityProfile3DVTKGrid(INITIALIZE);
    SampleCationAndAdsorptionSites(INITIALIZE);
    SampleDcTSTConfigurationFiles(INITIALIZE);
    SamplePDBMovies(INITIALIZE,-1);

    ClearLambdaHistogram();


    SimulationStage=PRODUCTION;
    for(CurrentCycle=0;CurrentCycle<NumberOfCycles;CurrentCycle++)
    {
      cpu_start=get_cpu_time();

      // detect erroneous chirality changes
      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        CheckChiralityMolecules();

      // sample energy average and system/particle properties
      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      {
        UpdateEnergyAveragesCurrentSystem();

        SampleRadialDistributionFunction(SAMPLE);
        SampleProjectedLengthsDistributionFunction(SAMPLE);
        SampleProjectedAnglesDistributionFunction(SAMPLE);
        SampleNumberOfMoleculesHistogram(SAMPLE);
        SamplePositionHistogram(SAMPLE);
        SampleFreeEnergyProfile(SAMPLE);
        SamplePoreSizeDistribution(SAMPLE);
        SampleEndToEndDistanceHistogram(SAMPLE);
        SampleEnergyHistogram(SAMPLE);
        SampleThermoDynamicsFactor(SAMPLE);
        SampleFrameworkSpacingHistogram(SAMPLE);
        SampleResidenceTimes(SAMPLE);
        SampleDistanceHistogram(SAMPLE);
        SampleBendAngleHistogram(SAMPLE);
        SampleDihedralAngleHistogram(SAMPLE);
        SampleAngleBetweenPlanesHistogram(SAMPLE);
        SampleMoleculePropertyHistogram(SAMPLE);
        SampleDensityProfile3DVTKGrid(SAMPLE);
        SampleCOMDensityProfile3DVTKGrid(SAMPLE);
        SampleCationAndAdsorptionSites(SAMPLE);
        SampleDcTSTConfigurationFiles(SAMPLE);
        SamplePDBMovies(SAMPLE,-1);
      }

      if(CurrentCycle%PrintPropertiesEvery==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
          PrintPropertyStatus(CurrentCycle,NumberOfCycles,OutputFilePtr[CurrentSystem]);
      }

      // Print at 'PrintEvery' intervals the status and a restart-file
      if((CurrentCycle%PrintEvery)==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          PrintIntervalStatusProduction(CurrentCycle,NumberOfCycles,OutputFilePtr[CurrentSystem]);
          PrintRestartFile();
        }
      }

      // select MC moves
      for(i=0;i<NumberOfSystems;i++)
      {
        // choose a random system
        SelectedSystem=(int)(RandomNumber()*(REAL)NumberOfSystems);

        //NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[SelectedSystem]+NumberOfCationMolecules[SelectedSystem]);
        NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[0]+NumberOfCationMolecules[0]);
        for(j=0;j<NumberOfSystems;j++)
        {
          NumberOfParticleMoves=MAX2(NumberOfParticleMoves,NumberOfAdsorbateMolecules[j]+NumberOfCationMolecules[j]);
        }
        NumberOfSteps=NumberOfParticleMoves*(NumberOfComponents==0?1:NumberOfComponents);

        // loop over the MC 'steps' per MC 'cycle'
        for(j=0;j<NumberOfSteps;j++)
        {
          // set the selected system
          CurrentSystem=SelectedSystem;

          // choose a random component
          CurrentComponent=(int)(RandomNumber()*(REAL)NumberOfComponents);

          // sample Lambda-histogram
          SampleLambdaHistogram();

          // choose any of the MC moves randomly with the selected probability
          ran=RandomNumber();

          if(ran<Components[CurrentComponent].ProbabilityTranslationMove)
          {
            cpu_before=get_cpu_time();
            TranslationMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeTranslationMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityRandomTranslationMove)
          {
            cpu_before=get_cpu_time();
            RandomTranslationMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeRandomTranslationMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityRotationMove)
          {
            cpu_before=get_cpu_time();
            RotationMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeRotationMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityRandomRotationMove)
          {
            cpu_before=get_cpu_time();
            RandomRotationMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeRandomRotationMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityPartialReinsertionMove)
          {
            cpu_before=get_cpu_time();
            PartialReinsertionMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimePartialReinsertionMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityReinsertionMove)
          {
            cpu_before=get_cpu_time();
            ReinsertionMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeReinsertionMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityReinsertionInPlaceMove)
          {
            cpu_before=get_cpu_time();
            ReinsertionInPlaceMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeReinsertionInPlaceMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityReinsertionInPlaneMove)
          {
            cpu_before=get_cpu_time();
            ReinsertionInPlaneMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeReinsertionInPlaneMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityIdentityChangeMove)
          {
            cpu_before=get_cpu_time();
            IdentityChangeMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeIdentityChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilitySwapMove)
          {
            if(RandomNumber()<0.5) 
            {
              cpu_before=get_cpu_time();
              SwapAddMove();
              cpu_after=get_cpu_time();
              Components[CurrentComponent].CpuTimeSwapMoveInsertion[CurrentSystem]+=(cpu_after-cpu_before);
            }
            else 
            {
              cpu_before=get_cpu_time();
              SwapRemoveMove();
              cpu_after=get_cpu_time();
              Components[CurrentComponent].CpuTimeSwapMoveDeletion[CurrentSystem]+=(cpu_after-cpu_before);
            }
          }
          else if(ran<Components[CurrentComponent].ProbabilityCFSwapLambdaMove)
          {
            cpu_before=get_cpu_time();
            CFSwapLambaMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCFSwapLambdaMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCBCFSwapLambdaMove)
          {
            cpu_before=get_cpu_time();
            CBCFSwapLambaMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCBCFSwapLambdaMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityWidomMove)
          {
            cpu_before=get_cpu_time();
            WidomMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeWidomMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCFWidomLambdaMove)
          {
            cpu_before=get_cpu_time();
            CFWidomLambaMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCFWidomLambdaMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityGibbsWidomMove)
          {
            cpu_before=get_cpu_time();
            GibbsWidomMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeGibbsWidomMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilitySurfaceAreaMove)
          {
            cpu_before=get_cpu_time();
            SurfaceAreaMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeSurfaceAreaMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityGibbsChangeMove)
          {
            cpu_before=get_cpu_time();
            GibbsParticleTransferMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeGibbsChangeMove[0]+=0.5*(cpu_after-cpu_before);
            Components[CurrentComponent].CpuTimeGibbsChangeMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove)
          {
            cpu_before=get_cpu_time();
            GibbsIdentityChangeMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeGibbsIdentityChangeMove[0]+=0.5*(cpu_after-cpu_before);
            Components[CurrentComponent].CpuTimeGibbsIdentityChangeMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCFGibbsChangeMove)
          {
            cpu_before=get_cpu_time();
            CFGibbsParticleTransferMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCFGibbsChangeMove[0]+=0.5*(cpu_after-cpu_before);
            Components[CurrentComponent].CpuTimeCFGibbsChangeMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCBCFGibbsChangeMove)
          {
            cpu_before=get_cpu_time();
            CBCFGibbsParticleTransferMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCBCFGibbsChangeMove[0]+=0.5*(cpu_after-cpu_before);
            Components[CurrentComponent].CpuTimeCBCFGibbsChangeMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityExchangeFractionalParticleMove)
          {
            cpu_before=get_cpu_time();
            ExchangeFractionalParticleMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeExchangeFractionalParticleMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCFGibbsSwapFractionalMoleculeToOtherBoxMove)
          {
            cpu_before=get_cpu_time();
            CFGibbsSwapFractionalMoleculeToOtherBoxMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove[0]+=0.5*(cpu_after-cpu_before);
            Components[CurrentComponent].CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove[0]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCFGibbsLambdaChangeMove)
          {
            cpu_before=get_cpu_time();
            CFGibbsLambdaChangeMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCFGibbsLambdaChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCFGibbsFractionalToIntegerMove)
          {
            cpu_before=get_cpu_time();
            CFGibbsFractionalToIntegerMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCFGibbsFractionalToIntegerMove[0]+=0.5*(cpu_after-cpu_before);
            Components[CurrentComponent].CpuTimeCFGibbsFractionalToIntegerMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityParallelTemperingMove) 
          {
            // note: cpu-timings done in the move itself
            ParallelTemperingMove();
          }
          else if(ran<Components[CurrentComponent].ProbabilityHyperParallelTemperingMove) 
          {
            // note: cpu-timings done in the move itself
            HyperParallelTemperingMove();
          }
          else if(ran<Components[CurrentComponent].ProbabilityParallelMolFractionMove) 
          {
            // note: cpu-timings done in the move itself
            ParallelMolFractionMove();
          }
          else if(ran<Components[CurrentComponent].ProbabilityChiralInversionMove) 
          {
            cpu_before=get_cpu_time();
            ChiralInversionMove();
            cpu_after=get_cpu_time();
            CpuTimeChiralInversionMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityHybridNVEMove) 
          {
            cpu_before=get_cpu_time();
            HybridNVEMove();
            cpu_after=get_cpu_time();
            CpuTimeHybridNVEMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityHybridNPHMove) 
          {
            cpu_before=get_cpu_time();
            HybridNPHMove();
            cpu_after=get_cpu_time();
            CpuTimeHybridNPHMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityHybridNPHPRMove) 
          {
            cpu_before=get_cpu_time();
            HybridNPHPRMove();
            cpu_after=get_cpu_time();
            CpuTimeHybridNPHPRMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityVolumeChangeMove) 
          {
            cpu_before=get_cpu_time();
            VolumeMove();
            cpu_after=get_cpu_time();
            CpuTimeVolumeChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityBoxShapeChangeMove) 
          {
            cpu_before=get_cpu_time();
            BoxShapeChangeMove();
            cpu_after=get_cpu_time();
            CpuTimeBoxShapeChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityGibbsVolumeChangeMove) 
          {
            cpu_before=get_cpu_time();
            GibbsVolumeMove();
            cpu_after=get_cpu_time();
            CpuTimeGibbsVolumeChangeMove[0]+=0.5*(cpu_after-cpu_before);
            CpuTimeGibbsVolumeChangeMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityFrameworkChangeMove) 
          {
            cpu_before=get_cpu_time();
            FrameworkChangeMove();
            cpu_after=get_cpu_time();
            CpuTimeFrameworkChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityFrameworkShiftMove) 
          {
            cpu_before=get_cpu_time();
            FrameworkShiftMove();
            cpu_after=get_cpu_time();
            CpuTimeFrameworkShiftMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCFCRXMCLambdaChangeMove)
          {
            cpu_before=get_cpu_time();
            CFCRXMCLambdaChangeMove();
            cpu_after=get_cpu_time();
            CpuTimeCFCRXMCLambdaChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          
          #ifdef DEBUG
            DebugEnergyStatus();
          #endif
        }
      }


      if(CurrentCycle%OptimizeAcceptenceEvery==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          OptimizeVolumeChangeAcceptence();
          OptimizeGibbsVolumeChangeAcceptence();
          OptimizeTranslationAcceptence();
          OptimizeRotationAcceptence();
          OptimizeFrameworkChangeAcceptence();
          OptimizeFrameworkShiftAcceptence();
          OptimizeCFLambdaAcceptence();
          OptimizeCFGibbsLambdaAcceptence();
          OptimizeCBCFLambdaChangeAcceptence();
          OptimizeCBCFGibbsLambdaChangeAcceptence();
          OptimizeRXMCLambdaChangeAcceptence();
          RescaleMaximumRotationAnglesSmallMC();
          OptimizeCFGibbsLambdaChangeAcceptence();
        }
      }

      cpu_end=get_cpu_time();
      CpuTotal+=(cpu_end-cpu_start);
      CpuTimeProductionRun+=(cpu_end-cpu_start);

      if((WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
        WriteBinaryRestartFiles();

      // a label to jump to for a restart, everything before here is skipped
      // this is the point where the previous binary restart file was written
      ContinueAfterCrashLabel3: ;

      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      {
        SampleRadialDistributionFunction(PRINT);
        SampleProjectedLengthsDistributionFunction(PRINT);
        SampleProjectedAnglesDistributionFunction(PRINT);
        SampleNumberOfMoleculesHistogram(PRINT);
        SamplePositionHistogram(PRINT);
        SampleFreeEnergyProfile(PRINT);
        SamplePoreSizeDistribution(PRINT);
        SampleEndToEndDistanceHistogram(PRINT);
        SampleEnergyHistogram(PRINT);
        SampleThermoDynamicsFactor(PRINT);
        SampleFrameworkSpacingHistogram(PRINT);
        SampleResidenceTimes(PRINT);
        SampleDistanceHistogram(PRINT);
        SampleBendAngleHistogram(PRINT);
        SampleDihedralAngleHistogram(PRINT);
        SampleAngleBetweenPlanesHistogram(PRINT);
        SampleMoleculePropertyHistogram(PRINT);
        SampleDensityProfile3DVTKGrid(PRINT);
	SampleCOMDensityProfile3DVTKGrid(PRINT);
        SampleCationAndAdsorptionSites(PRINT);
        SampleDcTSTConfigurationFiles(PRINT);
        if(Movies[CurrentSystem]&&(CurrentCycle%WriteMoviesEvery[CurrentSystem]==0))
          SamplePDBMovies(PRINT,-1);
      }
    }

    SimulationStage=FINISHED;

    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      OptimizeVolumeChangeAcceptence();
      OptimizeGibbsVolumeChangeAcceptence();
      OptimizeTranslationAcceptence();
      OptimizeRotationAcceptence();
      OptimizeFrameworkChangeAcceptence();
      OptimizeFrameworkShiftAcceptence();
      OptimizeCFLambdaAcceptence();
      OptimizeCFGibbsLambdaAcceptence();
      OptimizeCBCFLambdaChangeAcceptence();
      OptimizeCBCFGibbsLambdaChangeAcceptence();
      OptimizeRXMCLambdaChangeAcceptence();
      RescaleMaximumRotationAnglesSmallMC();
      OptimizeCFWidomAcceptence();
      OptimizeCFGibbsLambdaChangeAcceptence();
    }


    // write last binary restart-file
    // sampled properties can be remade from the binary restart-file
    if(WriteBinaryRestartFileEvery>0)
       WriteBinaryRestartFiles();

    ContinueAfterCrashLabel4: ;

    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      SampleRadialDistributionFunction(PRINT);
      SampleProjectedLengthsDistributionFunction(PRINT);
      SampleProjectedAnglesDistributionFunction(PRINT);
      SampleNumberOfMoleculesHistogram(PRINT);
      SamplePositionHistogram(PRINT);
      SampleFreeEnergyProfile(PRINT);
      SamplePoreSizeDistribution(PRINT);
      SampleEndToEndDistanceHistogram(PRINT);
      SampleEnergyHistogram(PRINT);
      SampleThermoDynamicsFactor(PRINT);
      SampleFrameworkSpacingHistogram(PRINT);
      SampleResidenceTimes(PRINT);
      SampleDistanceHistogram(PRINT);
      SampleBendAngleHistogram(PRINT);
      SampleDihedralAngleHistogram(PRINT);
      SampleAngleBetweenPlanesHistogram(PRINT);
      SampleMoleculePropertyHistogram(PRINT);
      SampleDensityProfile3DVTKGrid(PRINT);
      SampleCOMDensityProfile3DVTKGrid(PRINT);
      SampleCationAndAdsorptionSites(PRINT);
      SampleDcTSTConfigurationFiles(PRINT);
    }

    PrintPostSimulationStatus();

    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      PrintRestartFile();

    CloseOutputFile();
  }

  // set current prssure to the last one
  CurrentIsothermPressure=NumberOfIsothermPressures-1;


  // finalize output
  SampleRadialDistributionFunction(FINALIZE);
  SampleProjectedLengthsDistributionFunction(FINALIZE);
  SampleProjectedAnglesDistributionFunction(FINALIZE);
  SampleNumberOfMoleculesHistogram(FINALIZE);
  SamplePositionHistogram(FINALIZE);
  SampleFreeEnergyProfile(FINALIZE);
  SamplePoreSizeDistribution(FINALIZE);
  SampleEndToEndDistanceHistogram(FINALIZE);
  SampleEnergyHistogram(FINALIZE);
  SampleThermoDynamicsFactor(FINALIZE);
  SampleFrameworkSpacingHistogram(FINALIZE);
  SampleResidenceTimes(FINALIZE);
  SampleDistanceHistogram(FINALIZE);
  SampleBendAngleHistogram(FINALIZE);
  SampleDihedralAngleHistogram(FINALIZE);
  SampleAngleBetweenPlanesHistogram(FINALIZE);
  SampleMoleculePropertyHistogram(FINALIZE);
  SampleDensityProfile3DVTKGrid(FINALIZE);
  SampleCOMDensityProfile3DVTKGrid(FINALIZE);
  SampleCationAndAdsorptionSites(FINALIZE);
  SampleDcTSTConfigurationFiles(FINALIZE);
  SamplePDBMovies(FINALIZE,-1);
}

/*********************************************************************************************************
 * Name       | CheckEnergyOfSystem                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Routine to check the energy state of the system                                          *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void DebugEnergyStatus(void)
{
  int i;

  REAL StoredUHostBond,StoredUHostUreyBradley,StoredUHostBend,StoredUHostInversionBend;
  REAL StoredUHostTorsion,StoredUHostImproperTorsion,StoredUHostBondBond;
  REAL StoredUHostBendBend,StoredUHostBondBend,StoredUHostBondTorsion,StoredUHostBendTorsion;

  REAL StoredUAdsorbateBond,StoredUAdsorbateUreyBradley,StoredUAdsorbateBend,StoredUAdsorbateInversionBend;
  REAL StoredUAdsorbateTorsion,StoredUAdsorbateImproperTorsion,StoredUAdsorbateBondBond;
  REAL StoredUAdsorbateBendBend,StoredUAdsorbateBondBend,StoredUAdsorbateBondTorsion,StoredUAdsorbateBendTorsion;
  REAL StoredUAdsorbateIntraVDW,StoredUAdsorbateIntraChargeCharge;
  REAL StoredUAdsorbateIntraChargeBondDipole,StoredUAdsorbateIntraBondDipoleBondDipole;

  REAL StoredUCationBond,StoredUCationUreyBradley,StoredUCationBend,StoredUCationInversionBend;
  REAL StoredUCationTorsion,StoredUCationImproperTorsion,StoredUCationBondBond;
  REAL StoredUCationBendBend,StoredUCationBondBend,StoredUCationBondTorsion,StoredUCationBendTorsion;
  REAL StoredUCationIntraVDW,StoredUCationIntraChargeCharge;
  REAL StoredUCationIntraChargeBondDipole,StoredUCationIntraBondDipoleBondDipole;

  REAL StoredUHostHost,StoredUHostHostVDW,StoredUHostHostChargeChargeReal;
  REAL StoredUHostHostChargeBondDipoleReal,StoredUHostHostBondDipoleBondDipoleReal;
  REAL StoredUHostHostChargeChargeFourier,StoredUHostHostCoulomb;
  REAL StoredUHostHostChargeBondDipoleFourier,StoredUHostHostBondDipoleBondDipoleFourier;
  REAL StoredUHostAdsorbate,StoredUHostAdsorbateVDW,StoredUHostAdsorbateChargeChargeReal;
  REAL StoredUHostAdsorbateChargeBondDipoleReal,StoredUHostAdsorbateBondDipoleBondDipoleReal;
  REAL StoredUHostAdsorbateChargeChargeFourier,StoredUHostAdsorbateCoulomb;
  REAL StoredUHostAdsorbateChargeBondDipoleFourier,StoredUHostAdsorbateBondDipoleBondDipoleFourier;
  REAL StoredUHostCation,StoredUHostCationVDW,StoredUHostCationChargeChargeReal;
  REAL StoredUHostCationChargeBondDipoleReal,StoredUHostCationBondDipoleBondDipoleReal;
  REAL StoredUHostCationChargeChargeFourier,StoredUHostCationCoulomb;
  REAL StoredUHostCationChargeBondDipoleFourier,StoredUHostCationBondDipoleBondDipoleFourier;

  REAL StoredUAdsorbateAdsorbate,StoredUAdsorbateAdsorbateVDW,StoredUAdsorbateAdsorbateChargeChargeReal;
  REAL StoredUAdsorbateAdsorbateChargeBondDipoleReal,StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal;
  REAL StoredUAdsorbateAdsorbateChargeChargeFourier,StoredUAdsorbateAdsorbateCoulomb;
  REAL StoredUAdsorbateAdsorbateChargeBondDipoleFourier,StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier;
  REAL StoredUAdsorbateCation,StoredUAdsorbateCationVDW,StoredUAdsorbateCationChargeChargeReal;
  REAL StoredUAdsorbateCationChargeBondDipoleReal,StoredUAdsorbateCationBondDipoleBondDipoleReal;
  REAL StoredUAdsorbateCationChargeChargeFourier,StoredUAdsorbateCationCoulomb;
  REAL StoredUAdsorbateCationChargeBondDipoleFourier,StoredUAdsorbateCationBondDipoleBondDipoleFourier;
  REAL StoredUCationCation,StoredUCationCationVDW,StoredUCationCationChargeChargeReal;
  REAL StoredUCationCationChargeBondDipoleReal,StoredUCationCationBondDipoleBondDipoleReal;
  REAL StoredUCationCationChargeChargeFourier,StoredUCationCationCoulomb;
  REAL StoredUCationCationChargeBondDipoleFourier,StoredUCationCationBondDipoleBondDipoleFourier;
  REAL StoredUTotal,StoredUTailCorrection;

  REAL UHostPolarizationStored,UAdsorbatePolarizationStored,UCationPolarizationStored;
  REAL UHostBackPolarizationStored,UAdsorbateBackPolarizationStored,UCationBackPolarizationStored;


  // store all energies
  StoredUTotal=UTotal[CurrentSystem];
  StoredUTailCorrection=UTailCorrection[CurrentSystem];

  StoredUHostBond=UHostBond[CurrentSystem];
  StoredUHostUreyBradley=UHostUreyBradley[CurrentSystem];
  StoredUHostBend=UHostBend[CurrentSystem];
  StoredUHostInversionBend=UHostInversionBend[CurrentSystem];
  StoredUHostTorsion=UHostTorsion[CurrentSystem];
  StoredUHostImproperTorsion=UHostImproperTorsion[CurrentSystem];
  StoredUHostBondBond=UHostBondBond[CurrentSystem];
  StoredUHostBendBend=UHostBendBend[CurrentSystem];
  StoredUHostBondBend=UHostBondBend[CurrentSystem];
  StoredUHostBondTorsion=UHostBondTorsion[CurrentSystem];
  StoredUHostBendTorsion=UHostBendTorsion[CurrentSystem];

  StoredUAdsorbateBond=UAdsorbateBond[CurrentSystem];
  StoredUAdsorbateUreyBradley=UAdsorbateUreyBradley[CurrentSystem];
  StoredUAdsorbateBend=UAdsorbateBend[CurrentSystem];
  StoredUAdsorbateInversionBend=UAdsorbateInversionBend[CurrentSystem];
  StoredUAdsorbateTorsion=UAdsorbateTorsion[CurrentSystem];
  StoredUAdsorbateImproperTorsion=UAdsorbateImproperTorsion[CurrentSystem];
  StoredUAdsorbateBondBond=UAdsorbateBondBond[CurrentSystem];
  StoredUAdsorbateBendBend=UAdsorbateBendBend[CurrentSystem];
  StoredUAdsorbateBondBend=UAdsorbateBondBend[CurrentSystem];
  StoredUAdsorbateBondTorsion=UAdsorbateBondTorsion[CurrentSystem];
  StoredUAdsorbateBendTorsion=UAdsorbateBendTorsion[CurrentSystem];
  StoredUAdsorbateIntraVDW=UAdsorbateIntraVDW[CurrentSystem];
  StoredUAdsorbateIntraChargeCharge=UAdsorbateIntraChargeCharge[CurrentSystem];
  StoredUAdsorbateIntraChargeBondDipole=UAdsorbateIntraChargeBondDipole[CurrentSystem];
  StoredUAdsorbateIntraBondDipoleBondDipole=UAdsorbateIntraBondDipoleBondDipole[CurrentSystem];

  StoredUCationBond=UCationBond[CurrentSystem];
  StoredUCationUreyBradley=UCationUreyBradley[CurrentSystem];
  StoredUCationBend=UCationBend[CurrentSystem];
  StoredUCationInversionBend=UCationInversionBend[CurrentSystem];
  StoredUCationTorsion=UCationTorsion[CurrentSystem];
  StoredUCationImproperTorsion=UCationImproperTorsion[CurrentSystem];
  StoredUCationBondBond=UCationBondBond[CurrentSystem];
  StoredUCationBendBend=UCationBendBend[CurrentSystem];
  StoredUCationBondBend=UCationBondBend[CurrentSystem];
  StoredUCationBondTorsion=UCationBondTorsion[CurrentSystem];
  StoredUCationBendTorsion=UCationBendTorsion[CurrentSystem];
  StoredUCationIntraVDW=UCationIntraVDW[CurrentSystem];
  StoredUCationIntraChargeCharge=UCationIntraChargeCharge[CurrentSystem];
  StoredUCationIntraChargeBondDipole=UCationIntraChargeBondDipole[CurrentSystem];
  StoredUCationIntraBondDipoleBondDipole=UCationIntraBondDipoleBondDipole[CurrentSystem];

  StoredUHostHost=UHostHost[CurrentSystem];
  StoredUHostHostVDW=UHostHostVDW[CurrentSystem];
  StoredUHostHostChargeChargeReal=UHostHostChargeChargeReal[CurrentSystem];
  StoredUHostHostChargeChargeFourier=UHostHostChargeChargeFourier[CurrentSystem];
  StoredUHostHostChargeBondDipoleReal=UHostHostChargeBondDipoleReal[CurrentSystem];
  StoredUHostHostChargeBondDipoleFourier=UHostHostChargeBondDipoleFourier[CurrentSystem];
  StoredUHostHostBondDipoleBondDipoleReal=UHostHostBondDipoleBondDipoleReal[CurrentSystem];
  StoredUHostHostBondDipoleBondDipoleFourier=UHostHostBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUHostHostCoulomb=UHostHostCoulomb[CurrentSystem];

  StoredUHostAdsorbate=UHostAdsorbate[CurrentSystem];
  StoredUHostAdsorbateVDW=UHostAdsorbateVDW[CurrentSystem];
  StoredUHostAdsorbateChargeChargeReal=UHostAdsorbateChargeChargeReal[CurrentSystem];
  StoredUHostAdsorbateChargeBondDipoleReal=UHostAdsorbateChargeBondDipoleReal[CurrentSystem];
  StoredUHostAdsorbateBondDipoleBondDipoleReal=UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
  StoredUHostAdsorbateChargeChargeFourier=UHostAdsorbateChargeChargeFourier[CurrentSystem];
  StoredUHostAdsorbateChargeBondDipoleFourier=UHostAdsorbateChargeBondDipoleFourier[CurrentSystem];
  StoredUHostAdsorbateBondDipoleBondDipoleFourier=UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUHostAdsorbateCoulomb=UHostAdsorbateCoulomb[CurrentSystem];

  StoredUHostCation=UHostCation[CurrentSystem];
  StoredUHostCationVDW=UHostCationVDW[CurrentSystem];
  StoredUHostCationChargeChargeReal=UHostCationChargeChargeReal[CurrentSystem];
  StoredUHostCationChargeBondDipoleReal=UHostCationChargeBondDipoleReal[CurrentSystem];
  StoredUHostCationBondDipoleBondDipoleReal=UHostCationBondDipoleBondDipoleReal[CurrentSystem];
  StoredUHostCationChargeChargeFourier=UHostCationChargeChargeFourier[CurrentSystem];
  StoredUHostCationChargeBondDipoleFourier=UHostCationChargeBondDipoleFourier[CurrentSystem];
  StoredUHostCationBondDipoleBondDipoleFourier=UHostCationBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUHostCationCoulomb=UHostCationCoulomb[CurrentSystem];

  StoredUAdsorbateAdsorbate=UAdsorbateAdsorbate[CurrentSystem];
  StoredUAdsorbateAdsorbateVDW=UAdsorbateAdsorbateVDW[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeChargeReal=UAdsorbateAdsorbateChargeChargeReal[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeBondDipoleReal=UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem];
  StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal=UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeChargeFourier=UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeBondDipoleFourier=UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier=UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateCoulomb=UAdsorbateAdsorbateCoulomb[CurrentSystem];

  StoredUAdsorbateCation=UAdsorbateCation[CurrentSystem];
  StoredUAdsorbateCationVDW=UAdsorbateCationVDW[CurrentSystem];
  StoredUAdsorbateCationChargeChargeReal=UAdsorbateCationChargeChargeReal[CurrentSystem];
  StoredUAdsorbateCationChargeBondDipoleReal=UAdsorbateCationChargeBondDipoleReal[CurrentSystem];
  StoredUAdsorbateCationBondDipoleBondDipoleReal=UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem];
  StoredUAdsorbateCationChargeChargeFourier=UAdsorbateCationChargeChargeFourier[CurrentSystem];
  StoredUAdsorbateCationChargeBondDipoleFourier=UAdsorbateCationChargeBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateCationBondDipoleBondDipoleFourier=UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateCationCoulomb=UAdsorbateCationCoulomb[CurrentSystem];

  StoredUCationCation=UCationCation[CurrentSystem];
  StoredUCationCationVDW=UCationCationVDW[CurrentSystem];
  StoredUCationCationChargeChargeReal=UCationCationChargeChargeReal[CurrentSystem];
  StoredUCationCationChargeBondDipoleReal=UCationCationChargeBondDipoleReal[CurrentSystem];
  StoredUCationCationBondDipoleBondDipoleReal=UCationCationBondDipoleBondDipoleReal[CurrentSystem];
  StoredUCationCationChargeChargeFourier=UCationCationChargeChargeFourier[CurrentSystem];
  StoredUCationCationChargeBondDipoleFourier=UCationCationChargeBondDipoleFourier[CurrentSystem];
  StoredUCationCationBondDipoleBondDipoleFourier=UCationCationBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUCationCationCoulomb=UCationCationCoulomb[CurrentSystem];

  UHostPolarizationStored=UHostPolarization[CurrentSystem];
  UAdsorbatePolarizationStored=UAdsorbatePolarization[CurrentSystem];
  UCationPolarizationStored=UCationPolarization[CurrentSystem];

  UHostBackPolarizationStored=UHostBackPolarization[CurrentSystem];
  UAdsorbateBackPolarizationStored=UAdsorbateBackPolarization[CurrentSystem];
  UCationBackPolarizationStored=UCationBackPolarization[CurrentSystem];

  // store the structure-factors for the Ewald-summations
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    SaveCurrentEwaldStructureFactors(0,CurrentSystem);

 // NOTE: CalculateEnergy does not recompute the positions from the center-of-mass and orientation (otherwise the type needs to be switched too).
  CalculateEnergy();

  if(fabs(StoredUTotal-UTotal[CurrentSystem])<1e-4)
  {
    fprintf(stderr, "Energy status system [%d]: okay, %18.10f vs. %18.10f\n",CurrentSystem,StoredUTotal,UTotal[CurrentSystem]);
  }
  else
  {
    for(i=0;i<NumberOfSystems;i++)
    {
      fprintf(stderr, "\n\n");
      fprintf(stderr, "ERROR [energy status]: current energy (%18.10f) and true energy (%18.10f) are different in system %d\n",StoredUTotal*ENERGY_TO_KELVIN,UTotal[CurrentSystem]*ENERGY_TO_KELVIN,CurrentSystem);
      if(fabs(StoredUTotal)>1e-8) fprintf(stderr, "UTotal: %18.10f %18.10f\n",StoredUTotal*ENERGY_TO_KELVIN,UTotal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUTailCorrection)>1e-8) fprintf(stderr, "UTailCorrection: %18.10f %18.10f\n",StoredUTailCorrection*ENERGY_TO_KELVIN,UTailCorrection[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(StoredUHostBond)>1e-8) fprintf(stderr, "UHostBond: %18.10f %18.10f\n",StoredUHostBond*ENERGY_TO_KELVIN,UHostBond[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostUreyBradley)>1e-8) fprintf(stderr, "UHostUreyBradley: %18.10f %18.10f\n",StoredUHostUreyBradley*ENERGY_TO_KELVIN,UHostUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostBend)>1e-8) fprintf(stderr, "UHostBend: %18.10f %18.10f\n",StoredUHostBend*ENERGY_TO_KELVIN,UHostBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostInversionBend)>1e-8) fprintf(stderr, "UHostInversionBend: %18.10f %18.10f\n",StoredUHostInversionBend*ENERGY_TO_KELVIN,UHostInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostTorsion)>1e-8) fprintf(stderr, "UHostTorsion: %18.10f %18.10f\n",StoredUHostTorsion*ENERGY_TO_KELVIN,UHostTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostImproperTorsion)>1e-8) fprintf(stderr, "UHostImproperTorsion: %18.10f %18.10f\n",StoredUHostImproperTorsion*ENERGY_TO_KELVIN,UHostImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostBondBond)>1e-8) fprintf(stderr, "UHostBondBond: %18.10f %18.10f\n",StoredUHostBondBond*ENERGY_TO_KELVIN,UHostBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostBendBend)>1e-8) fprintf(stderr, "UHostBendBend: %18.10f %18.10f\n",StoredUHostBendBend*ENERGY_TO_KELVIN,UHostBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostBondBend)>1e-8) fprintf(stderr, "UHostBondBend: %18.10f %18.10f\n",StoredUHostBondBend*ENERGY_TO_KELVIN,UHostBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostBondTorsion)>1e-8) fprintf(stderr, "UHostBondTorsion: %18.10f %18.10f\n",StoredUHostBondTorsion*ENERGY_TO_KELVIN,UHostBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostBendTorsion)>1e-8) fprintf(stderr, "UHostBendTorsion: %18.10f %18.10f\n",StoredUHostBendTorsion*ENERGY_TO_KELVIN,UHostBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(StoredUAdsorbateBond)>1e-8) fprintf(stderr, "UAdsorbateBond: %18.10f %18.10f\n",StoredUAdsorbateBond*ENERGY_TO_KELVIN,UAdsorbateBond[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateUreyBradley)>1e-8) fprintf(stderr, "UAdsorbateUreyBradley: %18.10f %18.10f\n",StoredUAdsorbateUreyBradley*ENERGY_TO_KELVIN,UAdsorbateUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateBend)>1e-8) fprintf(stderr, "UAdsorbateBend: %18.10f %18.10f\n",StoredUAdsorbateBend*ENERGY_TO_KELVIN,UAdsorbateBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateInversionBend)>1e-8) fprintf(stderr, "UAdsorbateInversionBend: %18.10f %18.10f\n",StoredUAdsorbateInversionBend*ENERGY_TO_KELVIN,UAdsorbateInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateTorsion)>1e-8) fprintf(stderr, "UAdsorbateTorsion: %18.10f %18.10f\n",StoredUAdsorbateTorsion*ENERGY_TO_KELVIN,UAdsorbateTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateImproperTorsion)>1e-8) fprintf(stderr, "UAdsorbateImproperTorsion: %18.10f %18.10f\n",StoredUAdsorbateImproperTorsion*ENERGY_TO_KELVIN,UAdsorbateImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateBondBond)>1e-8) fprintf(stderr, "UAdsorbateBondBond: %18.10f %18.10f\n",StoredUAdsorbateBondBond*ENERGY_TO_KELVIN,UAdsorbateBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateBendBend)>1e-8) fprintf(stderr, "UAdsorbateBendBend: %18.10f %18.10f\n",StoredUAdsorbateBendBend*ENERGY_TO_KELVIN,UAdsorbateBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateBondBend)>1e-8) fprintf(stderr, "UAdsorbateBondBend: %18.10f %18.10f\n",StoredUAdsorbateBondBend*ENERGY_TO_KELVIN,UAdsorbateBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateBondTorsion)>1e-8) fprintf(stderr, "UAdsorbateBondTorsion: %18.10f %18.10f\n",StoredUAdsorbateBondTorsion*ENERGY_TO_KELVIN,UAdsorbateBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateBendTorsion)>1e-8) fprintf(stderr, "UAdsorbateBendTorsion: %18.10f %18.10f\n",StoredUAdsorbateBendTorsion*ENERGY_TO_KELVIN,UAdsorbateBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateIntraVDW)>1e-8) fprintf(stderr, "UAdsorbateIntraVDW: %18.10f %18.10f\n",StoredUAdsorbateIntraVDW*ENERGY_TO_KELVIN,UAdsorbateIntraVDW[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateIntraChargeCharge)>1e-8) fprintf(stderr, "UAdsorbateIntraChargeCharge: %18.10f %18.10f\n",StoredUAdsorbateIntraChargeCharge*ENERGY_TO_KELVIN,UAdsorbateIntraChargeCharge[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateIntraChargeBondDipole)>1e-8) fprintf(stderr, "UAdsorbateIntraChargeBondDipole: %18.10f %18.10f\n",StoredUAdsorbateIntraChargeBondDipole*ENERGY_TO_KELVIN,UAdsorbateIntraChargeBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateIntraBondDipoleBondDipole)>1e-8) fprintf(stderr, "UAdsorbateIntraBondDipoleBondDipole: %18.10f %18.10f\n",StoredUAdsorbateIntraBondDipoleBondDipole*ENERGY_TO_KELVIN,UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(StoredUCationBond)>1e-8) fprintf(stderr, "UCationBond: %18.10f %18.10f\n",StoredUCationBond*ENERGY_TO_KELVIN,UCationBond[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationUreyBradley)>1e-8) fprintf(stderr, "UCationUreyBradley: %18.10f %18.10f\n",StoredUCationUreyBradley*ENERGY_TO_KELVIN,UCationUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationBend)>1e-8) fprintf(stderr, "UCationBend: %18.10f %18.10f\n",StoredUCationBend*ENERGY_TO_KELVIN,UCationBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationInversionBend)>1e-8) fprintf(stderr, "UCationInversionBend: %18.10f %18.10f\n",StoredUCationInversionBend*ENERGY_TO_KELVIN,UCationInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationTorsion)>1e-8) fprintf(stderr, "UCationTorsion: %18.10f %18.10f\n",StoredUCationTorsion*ENERGY_TO_KELVIN,UCationTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationImproperTorsion)>1e-8) fprintf(stderr, "UCationImproperTorsion: %18.10f %18.10f\n",StoredUCationImproperTorsion*ENERGY_TO_KELVIN,UCationImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationBondBond)>1e-8) fprintf(stderr, "UCationBondBond: %18.10f %18.10f\n",StoredUCationBondBond*ENERGY_TO_KELVIN,UCationBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationBendBend)>1e-8) fprintf(stderr, "UCationBendBend: %18.10f %18.10f\n",StoredUCationBendBend*ENERGY_TO_KELVIN,UCationBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationBondBend)>1e-8) fprintf(stderr, "UCationBondBend: %18.10f %18.10f\n",StoredUCationBondBend*ENERGY_TO_KELVIN,UCationBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationBondTorsion)>1e-8) fprintf(stderr, "UCationBondTorsion: %18.10f %18.10f\n",StoredUCationBondTorsion*ENERGY_TO_KELVIN,UCationBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationBendTorsion)>1e-8) fprintf(stderr, "UCationBendTorsion: %18.10f %18.10f\n",StoredUCationBendTorsion*ENERGY_TO_KELVIN,UCationBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationIntraVDW)>1e-8) fprintf(stderr, "UCationIntraVDW: %18.10f %18.10f\n",StoredUCationIntraVDW*ENERGY_TO_KELVIN,UCationIntraVDW[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationIntraChargeCharge)>1e-8) fprintf(stderr, "UCationIntraChargeCharge: %18.10f %18.10f\n",StoredUCationIntraChargeCharge*ENERGY_TO_KELVIN,UCationIntraChargeCharge[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationIntraChargeBondDipole)>1e-8) fprintf(stderr, "UCationIntraChargeBondDipole: %18.10f %18.10f\n",StoredUCationIntraChargeBondDipole*ENERGY_TO_KELVIN,UCationIntraChargeBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationIntraBondDipoleBondDipole)>1e-8) fprintf(stderr, "UCationIntraBondDipoleBondDipole: %18.10f %18.10f\n",StoredUCationIntraBondDipoleBondDipole*ENERGY_TO_KELVIN,UCationIntraBondDipoleBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(StoredUHostHost)>1e-8) fprintf(stderr, "UHostHost: %18.10f %18.10f\n",StoredUHostHost*ENERGY_TO_KELVIN,UHostHost[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostHostVDW)>1e-8) fprintf(stderr, "UHostHostVDW: %18.10f %18.10f\n",StoredUHostHostVDW*ENERGY_TO_KELVIN,UHostHostVDW[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostHostChargeChargeReal)>1e-8) fprintf(stderr, "UHostHostChargeChargeReal: %18.10f %18.10f\n",StoredUHostHostChargeChargeReal*ENERGY_TO_KELVIN,UHostHostChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostHostChargeChargeFourier)>1e-8) fprintf(stderr, "UHostHostChargeChargeFourier: %18.10f %18.10f\n",StoredUHostHostChargeChargeFourier*ENERGY_TO_KELVIN,UHostHostChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostHostChargeBondDipoleReal)>1e-8) fprintf(stderr, "UHostHostChargeBondDipoleReal: %18.10f %18.10f\n",StoredUHostHostChargeBondDipoleReal*ENERGY_TO_KELVIN,UHostHostChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostHostChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UHostHostChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUHostHostChargeBondDipoleFourier*ENERGY_TO_KELVIN,UHostHostChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostHostBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UHostHostBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUHostHostBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,UHostHostBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostHostBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UHostHostBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUHostHostBondDipoleBondDipoleFourier*ENERGY_TO_KELVIN,UHostHostBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostHostCoulomb)>1e-8) fprintf(stderr, "UHostHostCoulomb: %18.10f %18.10f\n",StoredUHostHostCoulomb*ENERGY_TO_KELVIN,UHostHostCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(StoredUHostAdsorbate)>1e-8) fprintf(stderr, "UHostAdsorbate: %18.10f %18.10f\n",StoredUHostAdsorbate*ENERGY_TO_KELVIN,UHostAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostAdsorbateVDW)>1e-8) fprintf(stderr, "UHostAdsorbateVDW: %18.10f %18.10f\n",StoredUHostAdsorbateVDW*ENERGY_TO_KELVIN,UHostAdsorbateVDW[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostAdsorbateChargeChargeReal)>1e-8) fprintf(stderr, "UHostAdsorbateChargeChargeReal: %18.10f %18.10f\n",StoredUHostAdsorbateChargeChargeReal*ENERGY_TO_KELVIN,UHostAdsorbateChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostAdsorbateChargeBondDipoleReal)>1e-8) fprintf(stderr, "UHostAdsorbateChargeBondDipoleReal: %18.10f %18.10f\n",StoredUHostAdsorbateChargeBondDipoleReal*ENERGY_TO_KELVIN,UHostAdsorbateChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostAdsorbateBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UHostAdsorbateBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUHostAdsorbateBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostAdsorbateChargeChargeFourier)>1e-8) fprintf(stderr, "UHostAdsorbateChargeChargeFourier: %18.10f %18.10f\n",StoredUHostAdsorbateChargeChargeFourier*ENERGY_TO_KELVIN,UHostAdsorbateChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostAdsorbateChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UHostAdsorbateChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUHostAdsorbateChargeBondDipoleFourier*ENERGY_TO_KELVIN,UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostAdsorbateBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UHostAdsorbateBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUHostAdsorbateBondDipoleBondDipoleFourier*ENERGY_TO_KELVIN,UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostAdsorbateCoulomb)>1e-8) fprintf(stderr, "UHostAdsorbateCoulomb: %18.10f %18.10f\n",StoredUHostAdsorbateCoulomb*ENERGY_TO_KELVIN,UHostAdsorbateCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(StoredUHostCation)>1e-8) fprintf(stderr, "UHostCation: %18.10f %18.10f\n",StoredUHostCation*ENERGY_TO_KELVIN,UHostCation[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostCationVDW)>1e-8) fprintf(stderr, "UHostCationVDW: %18.10f %18.10f\n",StoredUHostCationVDW*ENERGY_TO_KELVIN,UHostCationVDW[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostCationChargeChargeReal)>1e-8) fprintf(stderr, "UHostCationChargeChargeReal: %18.10f %18.10f\n",StoredUHostCationChargeChargeReal*ENERGY_TO_KELVIN,UHostCationChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostCationChargeBondDipoleReal)>1e-8) fprintf(stderr, "UHostCationChargeBondDipoleReal: %18.10f %18.10f\n",StoredUHostCationChargeBondDipoleReal*ENERGY_TO_KELVIN,UHostCationChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostCationBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UHostCationBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUHostCationBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,UHostCationBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostCationChargeChargeFourier)>1e-8) fprintf(stderr, "UHostCationChargeChargeFourier: %18.10f %18.10f\n",StoredUHostCationChargeChargeFourier*ENERGY_TO_KELVIN,UHostCationChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostCationChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UHostCationChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUHostCationChargeBondDipoleFourier*ENERGY_TO_KELVIN,UHostCationChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostCationBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UHostCationBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUHostCationBondDipoleBondDipoleFourier*ENERGY_TO_KELVIN,UHostCationBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUHostCationCoulomb)>1e-8) fprintf(stderr, "UHostCationCoulomb: %18.10f %18.10f\n",StoredUHostCationCoulomb*ENERGY_TO_KELVIN,UHostCationCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(StoredUAdsorbateAdsorbate)>1e-8) fprintf(stderr, "UAdsorbateAdsorbate: %18.10f %18.10f\n",StoredUAdsorbateAdsorbate*ENERGY_TO_KELVIN,UAdsorbateAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateAdsorbateVDW)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateVDW: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateVDW*ENERGY_TO_KELVIN,UAdsorbateAdsorbateVDW[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateAdsorbateChargeChargeReal)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateChargeChargeReal: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateChargeChargeReal*ENERGY_TO_KELVIN,UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateAdsorbateChargeBondDipoleReal)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateChargeBondDipoleReal: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateChargeBondDipoleReal*ENERGY_TO_KELVIN,UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateAdsorbateChargeChargeFourier)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateChargeChargeFourier: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateChargeChargeFourier*ENERGY_TO_KELVIN,UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateAdsorbateChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateChargeBondDipoleFourier*ENERGY_TO_KELVIN,UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier*ENERGY_TO_KELVIN,UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateAdsorbateCoulomb)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateCoulomb: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateCoulomb*ENERGY_TO_KELVIN,UAdsorbateAdsorbateCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(StoredUAdsorbateCation)>1e-8) fprintf(stderr, "UAdsorbateCation: %18.10f %18.10f\n",StoredUAdsorbateCation*ENERGY_TO_KELVIN,UAdsorbateCation[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateCationVDW)>1e-8) fprintf(stderr, "UAdsorbateCationVDW: %18.10f %18.10f\n",StoredUAdsorbateCationVDW*ENERGY_TO_KELVIN,UAdsorbateCationVDW[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateCationChargeChargeReal)>1e-8) fprintf(stderr, "UAdsorbateCationChargeChargeReal: %18.10f %18.10f\n",StoredUAdsorbateCationChargeChargeReal*ENERGY_TO_KELVIN,UAdsorbateCationChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateCationChargeBondDipoleReal)>1e-8) fprintf(stderr, "UAdsorbateCationChargeBondDipoleReal: %18.10f %18.10f\n",StoredUAdsorbateCationChargeBondDipoleReal*ENERGY_TO_KELVIN,UAdsorbateCationChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateCationBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UAdsorbateCationBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUAdsorbateCationBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateCationChargeChargeFourier)>1e-8) fprintf(stderr, "UAdsorbateCationChargeChargeFourier: %18.10f %18.10f\n",StoredUAdsorbateCationChargeChargeFourier*ENERGY_TO_KELVIN,UAdsorbateCationChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateCationChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UAdsorbateCationChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUAdsorbateCationChargeBondDipoleFourier*ENERGY_TO_KELVIN,UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateCationBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UAdsorbateCationBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUAdsorbateCationBondDipoleBondDipoleFourier*ENERGY_TO_KELVIN,UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUAdsorbateCationCoulomb)>1e-8) fprintf(stderr, "UAdsorbateCationCoulomb: %18.10f %18.10f\n",StoredUAdsorbateCationCoulomb*ENERGY_TO_KELVIN,UAdsorbateCationCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(StoredUCationCation)>1e-8) fprintf(stderr, "UCationCation: %18.10f %18.10f\n",StoredUCationCation*ENERGY_TO_KELVIN,UCationCation[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationCationVDW)>1e-8) fprintf(stderr, "UCationCationVDW: %18.10f %18.10f\n",StoredUCationCationVDW*ENERGY_TO_KELVIN,UCationCationVDW[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationCationChargeChargeReal)>1e-8) fprintf(stderr, "UCationCationChargeChargeReal: %18.10f %18.10f\n",StoredUCationCationChargeChargeReal*ENERGY_TO_KELVIN,UCationCationChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationCationChargeBondDipoleReal)>1e-8) fprintf(stderr, "UCationCationChargeBondDipoleReal: %18.10f %18.10f\n",StoredUCationCationChargeBondDipoleReal*ENERGY_TO_KELVIN,UCationCationChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationCationBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UCationCationBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUCationCationBondDipoleBondDipoleReal*ENERGY_TO_KELVIN,UCationCationBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationCationChargeChargeFourier)>1e-8) fprintf(stderr, "UCationCationChargeChargeFourier: %18.10f %18.10f\n",StoredUCationCationChargeChargeFourier*ENERGY_TO_KELVIN,UCationCationChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationCationChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UCationCationChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUCationCationChargeBondDipoleFourier*ENERGY_TO_KELVIN,UCationCationChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationCationBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UCationCationBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUCationCationBondDipoleBondDipoleFourier*ENERGY_TO_KELVIN,UCationCationBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(StoredUCationCationCoulomb)>1e-8) fprintf(stderr, "UCationCationCoulomb: %18.10f %18.10f\n",StoredUCationCationCoulomb*ENERGY_TO_KELVIN,UCationCationCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(UHostPolarizationStored)>1e-8) fprintf(stderr, "UHostPolarization: %18.10f %18.10f\n",UHostPolarizationStored*ENERGY_TO_KELVIN,UHostPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(UAdsorbatePolarizationStored)>1e-8) fprintf(stderr, "UAdsorbatePolarization: %18.10f %18.10f\n",UAdsorbatePolarizationStored*ENERGY_TO_KELVIN,UAdsorbatePolarization[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(UCationPolarizationStored)>1e-8) fprintf(stderr, "UCationPolarization: %18.10f %18.10f\n",UCationPolarizationStored*ENERGY_TO_KELVIN,UCationPolarization[CurrentSystem]*ENERGY_TO_KELVIN);

      if(fabs(UHostBackPolarizationStored)>1e-8) fprintf(stderr, "UHostBackPolarization: %18.10f %18.10f\n",UHostBackPolarizationStored*ENERGY_TO_KELVIN,UHostBackPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(UAdsorbateBackPolarizationStored)>1e-8) fprintf(stderr, "UAdsorbateBackPolarization: %18.10f %18.10f\n",UAdsorbateBackPolarizationStored*ENERGY_TO_KELVIN,UAdsorbateBackPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
      if(fabs(UCationBackPolarizationStored)>1e-8) fprintf(stderr, "UCationBackPolarization: %18.10f %18.10f\n",UCationBackPolarizationStored*ENERGY_TO_KELVIN,UCationBackPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
      fclose(OutputFilePtr[i]);
    }
    exit(0);
  }

  UHostBond[CurrentSystem]=StoredUHostBond;
  UHostUreyBradley[CurrentSystem]=StoredUHostUreyBradley;
  UHostBend[CurrentSystem]=StoredUHostBend;
  UHostInversionBend[CurrentSystem]=StoredUHostInversionBend;
  UHostTorsion[CurrentSystem]=StoredUHostTorsion;
  UHostImproperTorsion[CurrentSystem]=StoredUHostImproperTorsion;
  UHostBondBond[CurrentSystem]=StoredUHostBondBond;
  UHostBendBend[CurrentSystem]=StoredUHostBendBend;
  UHostBondBend[CurrentSystem]=StoredUHostBondBend;
  UHostBondTorsion[CurrentSystem]=StoredUHostBondTorsion;
  UHostBendTorsion[CurrentSystem]=StoredUHostBendTorsion;

  UAdsorbateBond[CurrentSystem]=StoredUAdsorbateBond;
  UAdsorbateUreyBradley[CurrentSystem]=StoredUAdsorbateUreyBradley;
  UAdsorbateBend[CurrentSystem]=StoredUAdsorbateBend;
  UAdsorbateInversionBend[CurrentSystem]=StoredUAdsorbateInversionBend;
  UAdsorbateTorsion[CurrentSystem]=StoredUAdsorbateTorsion;
  UAdsorbateImproperTorsion[CurrentSystem]=StoredUAdsorbateImproperTorsion;
  UAdsorbateBondBond[CurrentSystem]=StoredUAdsorbateBondBond;
  UAdsorbateBendBend[CurrentSystem]=StoredUAdsorbateBendBend;
  UAdsorbateBondTorsion[CurrentSystem]=StoredUAdsorbateBondTorsion;
  UAdsorbateBondBend[CurrentSystem]=StoredUAdsorbateBondBend;
  UAdsorbateBendTorsion[CurrentSystem]=StoredUAdsorbateBendTorsion;
  UAdsorbateIntraVDW[CurrentSystem]=StoredUAdsorbateIntraVDW;
  UAdsorbateIntraChargeCharge[CurrentSystem]=StoredUAdsorbateIntraChargeCharge;
  UAdsorbateIntraChargeBondDipole[CurrentSystem]=StoredUAdsorbateIntraChargeBondDipole;
  UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]=StoredUAdsorbateIntraBondDipoleBondDipole;

  UCationBond[CurrentSystem]=StoredUCationBond;
  UCationUreyBradley[CurrentSystem]=StoredUCationUreyBradley;
  UCationBend[CurrentSystem]=StoredUCationBend;
  UCationInversionBend[CurrentSystem]=StoredUCationInversionBend;
  UCationTorsion[CurrentSystem]=StoredUCationTorsion;
  UCationImproperTorsion[CurrentSystem]=StoredUCationImproperTorsion;
  UCationBondBond[CurrentSystem]=StoredUCationBondBond;
  UCationBendBend[CurrentSystem]=StoredUCationBendBend;
  UCationBondBend[CurrentSystem]=StoredUCationBondBend;
  UCationBondTorsion[CurrentSystem]=StoredUCationBondTorsion;
  UCationBendTorsion[CurrentSystem]=StoredUCationBendTorsion;
  UCationIntraVDW[CurrentSystem]=StoredUCationIntraVDW;
  UCationIntraChargeCharge[CurrentSystem]=StoredUCationIntraChargeCharge;
  UCationIntraChargeBondDipole[CurrentSystem]=StoredUCationIntraChargeBondDipole;
  UCationIntraBondDipoleBondDipole[CurrentSystem]=StoredUCationIntraBondDipoleBondDipole;

  UHostHost[CurrentSystem]=StoredUHostHost;
  UHostHostVDW[CurrentSystem]=StoredUHostHostVDW;
  UHostHostChargeChargeReal[CurrentSystem]=StoredUHostHostChargeChargeReal;
  UHostHostChargeBondDipoleReal[CurrentSystem]=StoredUHostHostChargeBondDipoleReal;
  UHostHostBondDipoleBondDipoleReal[CurrentSystem]=StoredUHostHostBondDipoleBondDipoleReal;
  UHostHostChargeChargeFourier[CurrentSystem]=StoredUHostHostChargeChargeFourier;
  UHostHostChargeBondDipoleFourier[CurrentSystem]=StoredUHostHostChargeBondDipoleFourier;
  UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=StoredUHostHostBondDipoleBondDipoleFourier;
  UHostHostCoulomb[CurrentSystem]=StoredUHostHostCoulomb;

  UHostAdsorbate[CurrentSystem]=StoredUHostAdsorbate;
  UHostAdsorbateVDW[CurrentSystem]=StoredUHostAdsorbateVDW;
  UHostAdsorbateChargeChargeReal[CurrentSystem]=StoredUHostAdsorbateChargeChargeReal;
  UHostAdsorbateChargeBondDipoleReal[CurrentSystem]=StoredUHostAdsorbateChargeBondDipoleReal;
  UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=StoredUHostAdsorbateBondDipoleBondDipoleReal;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=StoredUHostAdsorbateChargeChargeFourier;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=StoredUHostAdsorbateChargeBondDipoleFourier;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=StoredUHostAdsorbateBondDipoleBondDipoleFourier;
  UHostAdsorbateCoulomb[CurrentSystem]=StoredUHostAdsorbateCoulomb;

  UHostCation[CurrentSystem]=StoredUHostCation;
  UHostCationVDW[CurrentSystem]=StoredUHostCationVDW;
  UHostCationChargeChargeReal[CurrentSystem]=StoredUHostCationChargeChargeReal;
  UHostCationChargeBondDipoleReal[CurrentSystem]=StoredUHostCationChargeBondDipoleReal;
  UHostCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUHostCationBondDipoleBondDipoleReal;
  UHostCationChargeChargeFourier[CurrentSystem]=StoredUHostCationChargeChargeFourier;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=StoredUHostCationChargeBondDipoleFourier;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUHostCationBondDipoleBondDipoleFourier;
  UHostCationCoulomb[CurrentSystem]=StoredUHostCationCoulomb;

  UAdsorbateAdsorbate[CurrentSystem]=StoredUAdsorbateAdsorbate;
  UAdsorbateAdsorbateVDW[CurrentSystem]=StoredUAdsorbateAdsorbateVDW;
  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=StoredUAdsorbateAdsorbateChargeChargeReal;
  UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=StoredUAdsorbateAdsorbateChargeBondDipoleReal;
  UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal;
  UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=StoredUAdsorbateAdsorbateChargeChargeFourier;
  UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=StoredUAdsorbateAdsorbateChargeBondDipoleFourier;
  UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier;
  UAdsorbateAdsorbateCoulomb[CurrentSystem]=StoredUAdsorbateAdsorbateCoulomb;

  UAdsorbateCation[CurrentSystem]=StoredUAdsorbateCation;
  UAdsorbateCationVDW[CurrentSystem]=StoredUAdsorbateCationVDW;
  UAdsorbateCationChargeChargeReal[CurrentSystem]=StoredUAdsorbateCationChargeChargeReal;
  UAdsorbateCationChargeBondDipoleReal[CurrentSystem]=StoredUAdsorbateCationChargeBondDipoleReal;
  UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUAdsorbateCationBondDipoleBondDipoleReal;
  UAdsorbateCationChargeChargeFourier[CurrentSystem]=StoredUAdsorbateCationChargeChargeFourier;
  UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=StoredUAdsorbateCationChargeBondDipoleFourier;
  UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUAdsorbateCationBondDipoleBondDipoleFourier;
  UAdsorbateCationCoulomb[CurrentSystem]=StoredUAdsorbateCationCoulomb;

  UCationCation[CurrentSystem]=StoredUCationCation;
  UCationCationVDW[CurrentSystem]=StoredUCationCationVDW;
  UCationCationChargeChargeReal[CurrentSystem]=StoredUCationCationChargeChargeReal;
  UCationCationChargeBondDipoleReal[CurrentSystem]=StoredUCationCationChargeBondDipoleReal;
  UCationCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUCationCationBondDipoleBondDipoleReal;
  UCationCationChargeChargeFourier[CurrentSystem]=StoredUCationCationChargeChargeFourier;
  UCationCationChargeBondDipoleFourier[CurrentSystem]=StoredUCationCationChargeBondDipoleFourier;
  UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUCationCationBondDipoleBondDipoleFourier;
  UCationCationCoulomb[CurrentSystem]=StoredUCationCationCoulomb;

  UHostPolarization[CurrentSystem]=UHostPolarizationStored;
  UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationStored;
  UCationPolarization[CurrentSystem]=UCationPolarizationStored;

  UHostBackPolarization[CurrentSystem]=UHostBackPolarizationStored;
  UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationStored;
  UCationBackPolarization[CurrentSystem]=UCationBackPolarizationStored;

  UTotal[CurrentSystem]=StoredUTotal;
  UTailCorrection[CurrentSystem]=StoredUTailCorrection;

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    RetrieveStoredEwaldStructureFactors(0,CurrentSystem);
}


/*********************************************************************************************************
 * Name       | ComputePoreSizeDistribution                                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Main routine to compute the Pore-Size Distribution function (PSD).                       *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void ComputePoreSizeDistribution(void)
{
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    ComputePSDHistogram[CurrentSystem]=TRUE;

  SamplePoreSizeDistribution(ALLOCATE);
  SamplePoreSizeDistribution(INITIALIZE);

  for(CurrentCycle=0;CurrentCycle<NumberOfCycles;CurrentCycle++)
  {
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      if((CurrentCycle%PrintEvery)==0)
      {
        fprintf(OutputFilePtr[CurrentSystem],"Iteration: %lld\n",CurrentCycle);
        fflush(OutputFilePtr[CurrentSystem]);
      }
      SamplePoreSizeDistribution(SAMPLE);
      SamplePoreSizeDistribution(PRINT);
    }
  }

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
     SamplePoreSizeDistribution(PRINT);

  SamplePoreSizeDistribution(FINALIZE);
}
