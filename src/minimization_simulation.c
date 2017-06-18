/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'minimization_simulation.c' is part of RASPA-2.0

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
#include "minimization_simulation.h"
#include "integration.h"
#include "simulation.h"
#include "molecule.h"
#include "framework_energy.h"
#include "framework_force.h"
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
#include "minimization.h"
#include "numerical.h"

void MinimalizationSimulation(void)
{
  int j,k;
  REAL Drift,ReferenceEnergy,ran;
  int NumberOfSystemMoves,NumberOfParticleMoves,NumberOfSteps;
  int ran_int;

  Drift=0.0;
  ReferenceEnergy=0.0;
  CurrentSystem=0;

  fprintf(stderr, "Starting minimization\n");

  OpenOutputFile();

  //InitializeNoseHooverAllSystems();

  InitializesEnergiesAllSystems();
  InitializesEnergyAveragesAllSystems();

  InitializeSmallMCStatisticsAllSystems();
  InitializeMCMovesStatisticsAllSystems();

  CalculateTotalEnergyAllSystems();

  PrintPreSimulationStatus();

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    for(k=0;k<NumberOfCycles;k++)
    {
      InitializesEnergiesAllSystems();
      InitializesEnergyAveragesAllSystems();
      CalculateTotalEnergyAllSystems();

      for(CurrentCycle=0;CurrentCycle<NumberOfInitializationCycles;CurrentCycle++)
      {
        if((CurrentCycle%PrintEvery)==0)
          PrintIntervalStatusInit(CurrentCycle,NumberOfInitializationCycles,OutputFilePtr[CurrentSystem]);

        NumberOfSystemMoves=12;
        NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem]);
        NumberOfSteps=(NumberOfSystemMoves+NumberOfParticleMoves)*NumberOfComponents;

        for(j=0;j<NumberOfSteps;j++)
        {
          ran_int=(int)(RandomNumber()*NumberOfSteps);
          switch(ran_int)
          {
            case 0:
              if(RandomNumber()<ProbabilityParallelTemperingMove) ParallelTemperingMove();
              break;
            case 1:
              if(RandomNumber()<ProbabilityHyperParallelTemperingMove) HyperParallelTemperingMove();
              break;
            case 2:
              if(RandomNumber()<ProbabilityParallelMolFractionMove) ParallelMolFractionMove();
              break;
            case 3:
              if(RandomNumber()<ProbabilityChiralInversionMove) ChiralInversionMove();
              break;
            case 4:
              if(RandomNumber()<ProbabilityHybridNVEMove) HybridNVEMove();
              break;
            case 5:
              if(RandomNumber()<ProbabilityHybridNPHMove) HybridNPHMove();
              break;
            case 6:
              if(RandomNumber()<ProbabilityHybridNPHPRMove) HybridNPHPRMove();
              break;
            case 7:
              if(RandomNumber()<ProbabilityVolumeChangeMove) VolumeMove();
              break;
            case 8:
              if(RandomNumber()<ProbabilityBoxShapeChangeMove) BoxShapeChangeMove();
              break;
            case 9:
              if(RandomNumber()<ProbabilityGibbsVolumeChangeMove) GibbsVolumeMove();
              break;
            case 10:
              if(RandomNumber()<ProbabilityFrameworkChangeMove) FrameworkChangeMove();
              break;
            case 11:
              if(RandomNumber()<ProbabilityFrameworkShiftMove) FrameworkShiftMove();
              break;
            default:
              // choose component at random
              CurrentComponent=(int)(RandomNumber()*(REAL)NumberOfComponents);

              // choose the Monte Carlo move at random
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
              else if(ran<Components[CurrentComponent].ProbabilityWidomMove)
               ;
              else if(ran<Components[CurrentComponent].ProbabilitySurfaceAreaMove)
               ;
              else if(ran<Components[CurrentComponent].ProbabilityGibbsChangeMove)
                GibbsParticleTransferMove();
              else if(ran<Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove)
               GibbsIdentityChangeMove();
          }
        }

        if(CurrentCycle%1000==0)
          OptimizeTranslationAcceptence();
      }

      // do the minimization
      Minimization(k);

      InitializesEnergiesCurrentSystem();
      CalculateEnergy();
      PrintEnergyStatus(OutputFilePtr[CurrentSystem],"minimization full energy");
    }

    PrintRestartFile();
  }

  PrintPostSimulationStatus();
  CloseOutputFile();
}

void GlobalMinimumSimulation(void)
{
  int i,j,k,l;
  int success;
  REAL ran,GlobalMinimumEnergy,UInitial,UAfterMC;

  success=0;

  GlobalMinimumEnergy=EnergyOverlapCriteria;

  OpenOutputFile();

  PrintPreSimulationStatus();

  for(k=0;k<NumberOfCycles;k++)
  {
    do
    {
      CurrentSystem=0;
      for(i=0;i<NumberOfAdsorbateMolecules[0];i++)
        RemoveAdsorbateMolecule();

      for(i=0;i<NumberOfCationMolecules[0];i++)
        RemoveCationMolecule();

      for(j=0;j<NumberOfComponents;j++)
      {
        if(Components[j].CreateNumberOfMolecules[0]>0)
        {
          CurrentSystem=0;
          if(Components[j].ExtraFrameworkMolecule)
            MakeInitialCations(Components[j].CreateNumberOfMolecules[0],j);
          else
            MakeInitialAdsorbates(Components[j].CreateNumberOfMolecules[0],j);
        }
      }
      CalculateForce();
      UInitial=UTotal[0];

      for(l=0;l<NumberOfInitializationCycles;l++)
      {
        for(j=0;j<MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem]);j++)
        {
          // choose component at random
          CurrentComponent=(int)(RandomNumber()*(REAL)NumberOfComponents);

          // choose the Monte Carlo move at random
          ran=RandomNumber();

          if(Components[CurrentComponent].ExtraFrameworkMolecule)
          {
            if(ran<Components[CurrentComponent].ProbabilityTranslationMove)
              TranslationMoveCation();
            if(ran<Components[CurrentComponent].ProbabilityRotationMove)
              RotationMoveCation();
            else if(ran<Components[CurrentComponent].ProbabilityReinsertionMove)
              ReinsertionCationMove();
          }
          else
          {
            if(ran<Components[CurrentComponent].ProbabilityTranslationMove)
              TranslationMoveAdsorbate();
            if(ran<Components[CurrentComponent].ProbabilityRotationMove)
              RotationMoveAdsorbate();
            else if(ran<Components[CurrentComponent].ProbabilityReinsertionMove)
              ReinsertionAdsorbateMove();
          }
        }
      }


      CalculateForce();
      UAfterMC=UTotal[0];

      // do the minimization as
      //success=BakerMinimizationNoOutput();
      if(success==0)
        fprintf(OutputFilePtr[0],"Minimization failed to convergence within %d steps\n",MaximumNumberOfMinimizationSteps);
    } while(success==0);

    // recompute energy
    CalculateForce();
    CurrentSystem=0;
    fprintf(OutputFilePtr[0],"iteration: %d Energy before Monte-Carlo: %18.10f [K] Energy after Monte-Carlo: %18.10f [K]\n",
      k,UInitial*ENERGY_TO_KELVIN,UAfterMC*ENERGY_TO_KELVIN);
    fflush(OutputFilePtr[0]);

    if(UTotal[0]<GlobalMinimumEnergy)
    {
      fprintf(OutputFilePtr[0],"found lower new minimum, positions saved to restart-file\n");
      GlobalMinimumEnergy=UTotal[0];
      PrintRestartFile();
    }

    fprintf(OutputFilePtr[0],"\t\tMinimization steps: %d Minimized energy: %18.10f [K] (%18.10f [kJ/mol]) Global minimum: %18.10f (%18.10f [kJ/mol])\n",
        success,UTotal[0]*ENERGY_TO_KELVIN,UTotal[0]*ENERGY_TO_KJ_PER_MOL,GlobalMinimumEnergy*ENERGY_TO_KELVIN,GlobalMinimumEnergy*ENERGY_TO_KJ_PER_MOL);
    fflush(OutputFilePtr[0]);
  }


  PrintPostSimulationStatus();
  CloseOutputFile();
}
