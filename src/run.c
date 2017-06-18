/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'run.c' is part of RASPA-2.0

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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include "simulation.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "utils.h"
#include "molecule.h"
#include "output.h"
#include "mc_moves.h"
#include "statistics.h"
#include "potentials.h"
#include "cbmc.h"
#include "monte_carlo.h"
#include "grids.h"
#include "ewald.h"
#include "input.h"
#include "recrossing.h"
#include "minimization_simulation.h"
#include "molecular_dynamics.h"
#include "spectra.h"
#include "spacegroup.h"
#include "inter_energy.h"
#include "numerical.h"
#include "matrix.h"
#include "vector.h"
#include "movies.h"
#include "status.h"
#include "framework.h"

extern bool STREAM = false;
extern char *INPUT_CRYSTAL = NULL;
extern char *PORE_SIZE_DISTRIBUTION_OUTPUT = NULL;
extern size_t PORE_SIZE_DISTRIBUTION_OUTPUT_SIZE = NULL;
extern char **FILE_CONTENTS = NULL;
extern size_t *FILE_SIZES = NULL;

/**
 * The core logic is separated from main to simplify wrapper functionality
 */
char* run(char *inputData, char *inputCrystal, char *raspaDir, bool stream)
{
  int i = 0, j = 0, k = 0, maxAtomsBonds = 0;
  size_t chars = 0;
  REAL energy, force_factor;
  char *output = NULL, *temp = NULL, *delimiter = NULL;

  // There are a lot of globals kicking around.
  INPUT_CRYSTAL = strdup(inputCrystal);
  RASPA_DIRECTORY = strdup(raspaDir);
  STREAM = stream;

  if (STREAM)
    ReadInput(inputData);
  else
    ReadInputFile(inputData);

  // write the initial positions to files
  if (!STREAM)
  {
    WriteFrameworkDefinitionShell("initial");
    WriteFrameworkDefinitionCSSR("initial");
    WriteFrameworkDefinitionGulp("initial");
    WriteFrameworkDefinitionVASP("initial");
    WriteFrameworkDefinitionPDB("initial");
    WriteFrameworkDefinitionTinker("initial");
    WriteFrameworkDefinitionMOL("initial");
    WriteFrameworkDefinitionCIF("initial");

    if(CreateTinkerInput)
    {
      WriteSnapshotTinker("initial");
      WriteFrameworkDefinitionTinker("tinker");
      WriteTinkerParameterFile();
      WriteTinkerKeyFile();
    }
  }

  // compute powder diffraction pattern
  if(ComputePowderDiffractionPattern)
    PowderDiffraction();

  // run a simulation-type specified in the input-file
  switch(SimulationType)
  {
    case NUMERICAL:
      CheckStatusNumerically();
      break;
    case MONTE_CARLO:
      MonteCarloSimulation();
      break;
    case MOLECULAR_DYNAMICS:
      MolecularDynamicsSimulation();
      break;
    case SPECTRA:
      VibrationalAnalysis();
      break;
    case MINIMIZATION:
      MinimalizationSimulation();
      break;
    case GLOBAL_MINIMIZATION:
      GlobalMinimumSimulation();
      break;
    case MAKE_GRID:
      OpenOutputFile();
      PrintPreSimulationStatus();
      Framework[0].FrameworkModel=FULL;
      CurrentSystem=0;
      MakeGrid();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case MAKE_ASCI_GRID:
      OpenOutputFile();
      PrintPreSimulationStatus();
      Framework[0].FrameworkModel=FULL;
      CurrentSystem=0;
      MakeASCIGrid();
      PrintPostSimulationStatus();
      CloseOutputFile();
    case VISUALIZATION:
      OpenOutputFile();
      PrintPreSimulationStatus();
      FreeEnergyProfile3D();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case PORE_SIZE_DISTRIBUTION:
      OpenOutputFile();
      PrintPreSimulationStatus();
      ComputePoreSizeDistribution();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case BARRIER_RECROSSING:
      OpenOutputFile();
      PrintPreSimulationStatus();
      BarrierRecrossing();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case TEST_SPACEGROUP:
      TestSpacegroup();
      break;
    case PLOT_POTENTIAL:
      OpenOutputFile();
      PrintPreSimulationStatus();
      CurrentSystem=0;
      for(i=0;i<=1000;i++)
      {
        PotentialGradient(ReturnPseudoAtomNumber("CH4_sp3"),ReturnPseudoAtomNumber("CH4_sp3"),SQR(i*CutOffVDW/1000),&energy,&force_factor,1.0);
        //TestFunction(i*CutOffChargeCharge/1000,&energy,&force_factor);
        fprintf(stderr, "%g %g %g\n",i*CutOffVDW/1000,energy,force_factor);
      }
      CloseOutputFile();
      break;
    case STATUS:
      OpenOutputFile();
      PrintPreSimulationStatus();
      CurrentSystem=0;
      Status();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case GENERATE_FRAMEWORK:
      OpenOutputFile();
      PrintPreSimulationStatus();
      CurrentSystem=0;
      GenerateFramework();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case USER:
      OpenOutputFile();
      PrintPreSimulationStatus();
      CurrentSystem=0;
      VibrationalAnalysis();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
  }

  // If streaming, merge the simulation contents into one string
  if (STREAM)
  {
    // Just returning relevant data on PSD sim. Note: Should we return full
    // simulation data, or is the overhead not worth it?
    if (SimulationType == PORE_SIZE_DISTRIBUTION)
    {
      output = calloc(PORE_SIZE_DISTRIBUTION_OUTPUT_SIZE + 1, sizeof(char));
      strcat(output, PORE_SIZE_DISTRIBUTION_OUTPUT);
      free(PORE_SIZE_DISTRIBUTION_OUTPUT);
    }
    // Returning full output on standard case.
    else
    {
      delimiter = "\n\nEND OF SIMULATION\n\n";

      for (i = 0; i < NumberOfSystems; i++)
      {
        chars += FILE_SIZES[i];
      }
      chars += strlen(delimiter) * (NumberOfSystems - 1) + 1;
      output = calloc(chars, sizeof(char));

      for (i = 0; i < NumberOfSystems-1; i++)
      {
        strcat(output, FILE_CONTENTS[i]);
        strcat(output, delimiter);
      }
      strcat(output, FILE_CONTENTS[NumberOfSystems - 1]);

      for (i = 0; i < NumberOfSystems; i++)
      {
        free(FILE_CONTENTS[i]);
      }
    }

  // Write the final positions to files
  } else {
    WriteFrameworkDefinitionCSSR("final");
    WriteFrameworkDefinitionGulp("final");
    WriteFrameworkDefinitionVASP("final");
    WriteFrameworkDefinitionPDB("final");
    WriteFrameworkDefinitionMOL("final");
    WriteFrameworkDefinitionCIF("final");
    if(CreateTinkerInput)
      WriteSnapshotTinker("final");
  }

  free(FILE_CONTENTS);
  free(FILE_SIZES);
  free(INPUT_CRYSTAL);

  // Fighting memory leaks caused by unfreed globals.
  // This all came from valgrind.
  free(Eikx0);
  free(Eiky0);
  free(Eikz0);
  free(Eikx_bd0);
  free(Eiky_bd0);
  free(Eikz_bd0);
  free(Eikr_bd);
  free(Eikr_xy_bd);

  for (i = 0; i < NumberOfPseudoAtoms; i++)
  {
    free(PotentialParms[i]);
  }
  free(PotentialParms);

  for (i = 0; i < 80; i++)
  {
    free(CFVDWScalingRXMC[i]);
    free(CFChargeScalingRXMC[i]);
  }
  free(CFVDWScalingRXMC);
  free(CFChargeScalingRXMC);

  for (i = 0; i < NumberOfSystems; i++)
  {
    for (j = 0; j < Framework[i].NumberOfFrameworks; j++)
    {
      maxAtomsBonds = MAX2(Framework[i].NumberOfAtoms[j],
                           Framework[i].NumberOfBonds[j]);
      for (k = 0; k < maxAtomsBonds; k++)
      {
        free(Framework[i].ExclusionMatrix[j][k]);
      }
      free(Framework[i].ExclusionMatrix[j]);

      for (k = 0; k < TotalNumberOfReplicaCells[i] * Framework[i].NumberOfAtoms[j]; k++)
      {
        free(Framework[i].Neighbours[j][k]);
      }
      free(Framework[i].Neighbours[j]);

      free(Framework[i].Atoms[j]);
    }
    free(Framework[i].Atoms);
    free(Framework[i].ExclusionMatrix);
    free(Framework[i].Neighbours);
  }

  return output;
}
