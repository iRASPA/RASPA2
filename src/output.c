/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2015 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'output.c' is part of RASPA-2.0

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
#include <float.h>
#include <sys/stat.h>
#include <sys/sysctl.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "simulation.h"
#include "ewald.h"
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
#include "movies.h"
#include "sample.h"
#include "spectra.h"
#include "grids.h"
#include "statistics.h"
#include "thermo_baro_stats.h"
#include "equations_of_state.h"
#include "minimization.h"
#include "recrossing.h"
#include "warnings.h"
#include "rigid.h"
#include "spacegroup.h"

extern bool STREAM;
extern char **FILE_CONTENTS;
extern size_t *FILE_SIZES;

static REAL UHostHostRunning,UHostHostVDWRunning,UHostHostCoulombRunning;
static REAL UHostHostChargeChargeRealRunning,UHostHostChargeChargeFourierRunning;
static REAL UHostHostChargeBondDipoleRealRunning,UHostHostChargeBondDipoleFourierRunning;
static REAL UHostHostBondDipoleBondDipoleRealRunning,UHostHostBondDipoleBondDipoleFourierRunning;

static REAL UHostAdsorbateRunning,UHostAdsorbateVDWRunning,UHostAdsorbateCoulombRunning;
static REAL UHostAdsorbateChargeChargeRealRunning,UHostAdsorbateChargeChargeFourierRunning;
static REAL UHostAdsorbateChargeBondDipoleRealRunning,UHostAdsorbateChargeBondDipoleFourierRunning;
static REAL UHostAdsorbateBondDipoleBondDipoleRealRunning,UHostAdsorbateBondDipoleBondDipoleFourierRunning;

static REAL UHostCationRunning,UHostCationVDWRunning,UHostCationCoulombRunning;
static REAL UHostCationChargeChargeRealRunning,UHostCationChargeChargeFourierRunning;
static REAL UHostCationChargeBondDipoleRealRunning,UHostCationChargeBondDipoleFourierRunning;
static REAL UHostCationBondDipoleBondDipoleRealRunning,UHostCationBondDipoleBondDipoleFourierRunning;

static REAL UAdsorbateAdsorbateRunning,UAdsorbateAdsorbateVDWRunning,UAdsorbateAdsorbateCoulombRunning;
static REAL UAdsorbateAdsorbateChargeChargeRealRunning,UAdsorbateAdsorbateChargeChargeFourierRunning;
static REAL UAdsorbateAdsorbateChargeBondDipoleRealRunning,UAdsorbateAdsorbateChargeBondDipoleFourierRunning;
static REAL UAdsorbateAdsorbateBondDipoleBondDipoleRealRunning,UAdsorbateAdsorbateBondDipoleBondDipoleFourierRunning;

static REAL UCationCationRunning,UCationCationVDWRunning,UCationCationCoulombRunning;
static REAL UCationCationChargeChargeRealRunning,UCationCationChargeChargeFourierRunning;
static REAL UCationCationChargeBondDipoleRealRunning,UCationCationChargeBondDipoleFourierRunning;
static REAL UCationCationBondDipoleBondDipoleRealRunning,UCationCationBondDipoleBondDipoleFourierRunning;

static REAL UAdsorbateCationRunning,UAdsorbateCationVDWRunning,UAdsorbateCationCoulombRunning;
static REAL UAdsorbateCationChargeChargeRealRunning,UAdsorbateCationChargeChargeFourierRunning;
static REAL UAdsorbateCationChargeBondDipoleRealRunning,UAdsorbateCationChargeBondDipoleFourierRunning;
static REAL UAdsorbateCationBondDipoleBondDipoleRealRunning,UAdsorbateCationBondDipoleBondDipoleFourierRunning;

static REAL UHostBondRunning,UHostUreyBradleyRunning,UHostBendRunning;
static REAL UHostInversionBendRunning,UHostTorsionRunning,UHostImproperTorsionRunning,UHostOutOfPlaneRunning;
static REAL UHostBondBondRunning,UHostBondBendRunning,UHostBendBendRunning;
static REAL UHostBondTorsionRunning,UHostBendTorsionRunning;

static REAL UAdsorbateBondRunning,UAdsorbateUreyBradleyRunning,UAdsorbateBendRunning;
static REAL UAdsorbateInversionBendRunning,UAdsorbateTorsionRunning,UAdsorbateImproperTorsionRunning,UAdsorbateOutOfPlaneRunning;
static REAL UAdsorbateBondBondRunning,UAdsorbateBondBendRunning,UAdsorbateBendBendRunning;
static REAL UAdsorbateBondTorsionRunning,UAdsorbateBendTorsionRunning;
static REAL UAdsorbateIntraVDWRunning,UAdsorbateIntraChargeChargeCoulombRunning;
static REAL UAdsorbateIntraChargeBondDipoleCoulombRunning,UAdsorbateIntraBondDipoleBondDipoleCoulombRunning;

static REAL UCationBondRunning,UCationUreyBradleyRunning,UCationBendRunning;
static REAL UCationInversionBendRunning,UCationTorsionRunning,UCationImproperTorsionRunning,UCationOutOfPlaneRunning;
static REAL UCationBondBondRunning,UCationBondBendRunning,UCationBendBendRunning;
static REAL UCationBondTorsionRunning,UCationBendTorsionRunning;
static REAL UCationIntraVDWRunning,UCationIntraChargeChargeCoulombRunning;
static REAL UCationIntraChargeBondDipoleCoulombRunning,UCationIntraBondDipoleBondDipoleCoulombRunning;

static REAL UHostPolarizationRunning,UAdsorbatePolarizationRunning,UCationPolarizationRunning;
static REAL UHostBackPolarizationRunning,UAdsorbateBackPolarizationRunning,UCationBackPolarizationRunning;
static REAL UTailCorrectionRunning,UTotalRunning;
static REAL UDistanceConstraintsRunning,UAngleConstraintsRunning,UDihedralConstraintsRunning;
static REAL UInversionBendConstraintsRunning,UOutOfPlaneDistanceConstraintsRunning,UExclusionConstraintsRunning;

FILE **OutputFilePtr;

void OpenOutputFile(void)
{
  int i;
  char buffer[1024],buffer2[256];

  OutputFilePtr=(FILE**)calloc(NumberOfSystems,sizeof(FILE*));
  FILE_CONTENTS = (char**)malloc(NumberOfSystems * sizeof(char*));
  FILE_SIZES=(size_t*)malloc(NumberOfSystems * sizeof(size_t));

  if (STREAM)
  {
#ifdef __unix__
    // Loads output contents into a global
    for(i=0;i<NumberOfSystems;i++)
      OutputFilePtr[i]=open_memstream(&FILE_CONTENTS[i], &FILE_SIZES[i]);
#else
    fprintf(stderr, "Streaming only allowed on POSIX systems (for now)\n.");
    exit(1);
#endif
  }
  else
  {
    mkdir("Output",S_IRWXU);
    for(i=0;i<NumberOfSystems;i++)
    {
      sprintf(buffer,"Output/System_%d",i);
      mkdir(buffer,S_IRWXU);
    }
    for(i=0;i<NumberOfSystems;i++)
    {
      sprintf(buffer,"Output/System_%d/output_%s_%d.%d.%d_%lf_%lg%s",
              i,
              Framework[i].Name[0],
              NumberOfUnitCells[i].x,
              NumberOfUnitCells[i].y,
              NumberOfUnitCells[i].z,
              (double)therm_baro_stats.ExternalTemperature[i],
              (double)(therm_baro_stats.ExternalPressure[i][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
              FileNameAppend);

      // limit length of file-name
      strncpy(buffer2,buffer,250);
      sprintf(buffer2,"%s.data",buffer2);
      OutputFilePtr[i]=fopen(buffer2,"w");
    }
  }
}

void CloseOutputFile(void)
{
  int i;

  for(i=0;i<NumberOfSystems;i++)
    fclose(OutputFilePtr[i]);
}

void PrintPreSimulationStatus(void)
{
  int i;
  for(i=0;i<NumberOfSystems;i++)
    PrintPreSimulationStatusCurrentSystem(i);

  CurrentSystem=0;
}

void PrintPreSimulationStatusCurrentSystem(int system)
{
  int i,j,k,m,l;
  int ncell,k1,k2,k3;
  int A,B,Type,nr_args;
  int nr_free,nr_fixed;
  REAL charge,NetCharge;
  REAL smallest_charge,largest_charge;
  VECTOR Dipole,eigenvalues;
  REAL_MATRIX3x3 Quadrupole,Eigenvectors;
  char charge_string[256];
  char framework_charge_string[256];
  char polarization_string[256];
  FILE *FilePtr;
  char my_date[] = "Compile Date = " __DATE__;
  char my_time[] = "Compile Time = " __TIME__;
  time_t curtime;
  struct tm *loctime;
  char buffer[256],buffer1[256],buffer2[256],buffer3[256],buffer4[256],buffer5[256],buffer6[256];
  #if defined (__linux__)|| defined(__linux)
    int mib[2];
    size_t len;
  #endif
  #if defined (__APPLE__)
    size_t len;
    char cpudata[128],cpumodel[128],hostname[128];
    char osrelease[128],ostype[128],osversion[128];
  #endif

  // see what compiler is used
  #if defined(__GNUC__)
    #define COMPILER_STRING "gcc " __VERSION__
  #elif defined(__INTEL_COMPILER)
    #define COMPILER_STRING "Intel icc " __INTEL_COMPILER_BUILD_DATE
  #elif defined(__PGI)
    #define COMPILER_STRING "PGI compiler "
  #elif defined(__SUNPRO_C)
    #define COMPILER_STRING "Sun compiler "
  #elif defined(__sgi)
    #define COMPILER_STRING "SGI compiler " _COMPILER_VERSION
  #elif defined(__HP_aCC)
    #define COMPILER_STRING "HP compiler "
  #elif defined(__DECC)
    #define COMPILER_STRING "HP-Compaq-Digital compiler " __DECC_VER
  #elif defined(__xlC__)
    #define COMPILER_STRING "IBM compiler " __IBMC__
  #elif defined(_MSC_VER)
    #define COMPILER_STRING "Microsoft compiler "
  #elif defined(__BORLANDC__)
    #define COMPILER_STRING "Borland compiler "
  #elif defined(__MWERKS__)
    #define COMPILER_STRING "Metrowerks CodeWarrior "
  #else
    #define COMPILER_STRING "unknown compiler/version"
  #endif

  FilePtr=OutputFilePtr[system];

  fprintf(FilePtr,"Compiler and run-time data\n");
  fprintf(FilePtr,"===========================================================================\n");

  fprintf(FilePtr,"%s\n","RASPA 2.0");

  #if defined (__LP64__) || defined (__64BIT__) || defined (_LP64) || (__WORDSIZE == 64)
    fprintf(FilePtr,"Compiled as a 64-bits application\n");
    fprintf(FilePtr,"Compiler: %s\n",COMPILER_STRING);
    fprintf(FilePtr,"%s, %s\n\n",my_date,my_time);
  #else
    fprintf(FilePtr,"Compiled as a 32-bits application\n");
    fprintf(FilePtr,"Compiler: %s\n",COMPILER_STRING);
    fprintf(FilePtr,"%s, %s\n\n",my_date,my_time);
  #endif

  /* Get the current time.  */
  curtime = time (NULL);

  /* Convert it to local time representation.  */
  loctime = localtime (&curtime);

  /* Print out the date and time in the standard format.  */
  fprintf(FilePtr,"%s",asctime(loctime));

  /* Print it out in a nice format.  */
  strftime (buffer, 256, "Simulation started on %A, %B %d.", loctime);
  fprintf(FilePtr,"%s\n",buffer);
  strftime (buffer, 256, "The start time was %I:%M %p.", loctime);
  fprintf(FilePtr,"%s\n\n",buffer);

  // get hostname and cpu-info for linux
  #if defined (__linux__)|| defined(__linux)
    len = sizeof(buffer);

    mib[0] = CTL_KERN;
    mib[1] = KERN_NODENAME;
    sysctl(mib, 2, &buffer, &len, NULL, 0);
    fprintf(FilePtr,"Hostname:    %s\n",buffer);

    mib[0] = CTL_KERN;
    mib[1] = KERN_OSTYPE;
    sysctl(mib, 2, &buffer, &len, NULL, 0);
    fprintf(FilePtr,"OS type:     %s\n",buffer);

    mib[0] = CTL_KERN;
    mib[1] = KERN_OSRELEASE;
    sysctl(mib, 2, &buffer, &len, NULL, 0);
    fprintf(FilePtr,"OS release:  %s\n",buffer);

    mib[0] = CTL_KERN;
    mib[1] = KERN_VERSION;
    sysctl(mib, 2, &buffer, &len, NULL, 0);
    fprintf(FilePtr,"OS version:  %s\n",buffer);

    fprintf(FilePtr,"\n");
  #endif

  // get hostname and cpu-info for mac osx
  #if defined (__APPLE__)
    len = sizeof(cpudata);
    sysctlbyname("hw.machine", &cpudata, &len, NULL, 0);
    fprintf(FilePtr,"Cpu data:    %s\n",cpudata);

    len = sizeof(cpumodel);
    sysctlbyname("hw.model", &cpumodel, &len, NULL, 0);
    fprintf(FilePtr,"Cpu Model:   %s\n",cpumodel);

    len = sizeof(hostname);
    sysctlbyname("kern.hostname", &hostname, &len, NULL, 0);
    fprintf(FilePtr,"Host name:   %s\n",hostname);

    len = sizeof(osrelease);
    sysctlbyname("kern.osrelease", &osrelease, &len, NULL, 0);
    fprintf(FilePtr,"OS release:  %s\n",osrelease);

    len = sizeof(ostype);
    sysctlbyname("kern.ostype", &ostype, &len, NULL, 0);
    fprintf(FilePtr,"OS type:     %s\n",ostype);

    len = sizeof(osversion);
    sysctlbyname("kern.osversion", &osversion, &len, NULL, 0);
    fprintf(FilePtr,"OS version:  %s\n\n\n",osversion);
  #endif

  fprintf(FilePtr,"Simulation\n");
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"Dimensions: %d\n",Dimension);
  fprintf(FilePtr,"Random number seed: %lu\n",seed);
  fprintf(FilePtr,"RASPA directory set to: %s\n",RASPA_DIRECTORY);
  fprintf(FilePtr,"String appended to output-files: %s\n",FileNameAppend);
  fprintf(FilePtr,"Number of cycles: %lld\n",NumberOfCycles);
  fprintf(FilePtr,"Number of initializing cycles: %lld\n",NumberOfInitializationCycles);
  fprintf(FilePtr,"Number of equilibration cycles: %lld\n",NumberOfEquilibrationCycles);
  fprintf(FilePtr,"Print every: %d\n",PrintEvery);
  switch(BoundaryCondition[system])
  {
    case UNINITIALIZED_BOUNDARY_CONDITION:
      fprintf(FilePtr,"ERROR: Boundary Condition not set\n");
      exit(0);
      break;
    case FINITE:
      fprintf(FilePtr,"Finite system, no boundary condition applied\n");
      break;
    case CUBIC:
      fprintf(FilePtr,"Cubic boundary condition applied\n");
      break;
    case RECTANGULAR:
      fprintf(FilePtr,"Rectangular boundary condition applied\n");
      break;
    case TRICLINIC:
      fprintf(FilePtr,"Triclinic boundary condition applied\n");
      break;
  }
  fprintf(FilePtr,"Timestep: %lf\n",DeltaT);
  switch(InitEnsemble[system])
  {
    case NVE:
      fprintf(FilePtr,"Initialization Ensemble: NVE (constant number of particles N, constant volume V, constant energy E)\n");
      break;
    case NVT:
      fprintf(FilePtr,"Initialization Ensemble: NVT (constant number of particles N, constant volume V, constant average temperature T)\n");
      break;
    case NPT:
      fprintf(FilePtr,"Initialization Ensemble: NPT (constant number of particles N, constant average pressure P, constant average temperature T)\n");
      break;
    case NPH:
      fprintf(FilePtr,"Initialization Ensemble: NPH (constant number of particles N, constant average pressure P, constant enthalpy H)\n");
      break;
    case MuPT:
      fprintf(FilePtr,"Initialization Ensemble: MuPT (constant chemical potential mu, constant average pressure P, constant average temperature T)\n");
      break;
    case NPTPR:
      switch(NPTPRCellType[system])
      {
        case REGULAR:
          fprintf(FilePtr,"Initialization Ensemble: NPT (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
          break;
        case ANISOTROPIC:
          fprintf(FilePtr,"Initialization Ensemble: NPT angles fixed, box-lengths are allowed to vary anisotropically (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
          break;
        case ISOTROPIC:
          fprintf(FilePtr,"Initialization Ensemble: NPT angles fixed, box-lengths are allowed to vary isotropically (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPT beta/gamma angles fixed, alpha angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_BETA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPT alpha/gamma angles fixed, beta angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPT alpha/beta angles fixed, gamma angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          fprintf(FilePtr,"Initialization Ensemble: NPT full cell fluctuations (Parinello-Rahman) using the upper triangle of the box-matrix\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPT beta/gamma angles fixed, alpha angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_BETA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPT alpha/gamma angles fixed, beta angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPT alpha/beta angles fixed, gamma angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
          }
          break;
        default:
          fprintf(stderr, "Undefined cell type for NPTPR\n");
          exit(0);
          break;
      }
      break;
    case NPHPR:
      switch(NPTPRCellType[system])
      {
        case REGULAR:
          fprintf(FilePtr,"Initialization Ensemble: NPH (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant enthalpy H)\n");
          break;
        case ANISOTROPIC:
          fprintf(FilePtr,"Initialization Ensemble: NPH angles fixed, box-lengths are allowed to vary anisotropically (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant enthalpy H)\n");
          break;
        case ISOTROPIC:
          fprintf(FilePtr,"Initialization Ensemble: NPH angles fixed, box-lengths are allowed to vary isotropically (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant enthalpy H)\n");
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPH beta/gamma angles fixed, alpha angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_BETA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPH alpha/gamma angles fixed, beta angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPH alpha/beta angles fixed, gamma angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          fprintf(FilePtr,"Initialization Ensemble: NPH full cell fluctuations (Parinello-Rahman) using the upper triangle of the box-matrix\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant enthalpy)\n");
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPH beta/gamma angles fixed, alpha angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_BETA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPH alpha/gamma angles fixed, beta angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              fprintf(FilePtr,"Initialization Ensemble: NPH alpha/beta angles fixed, gamma angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
          }
          break;
        default:
          fprintf(stderr, "Undefined cell type for NPHPR\n");
          exit(0);
          break;
      }
      break;
  }

  switch(RunEnsemble[system])
  {
    case NVE:
      fprintf(FilePtr,"Production run ensemble: NVE (constant number of particles N, constant volume V, constant energy E)\n");
      break;
    case NVT:
      fprintf(FilePtr,"Production run ensemble: NVT (constant number of particles N, constant volume V, constant average temperature T)\n");
      break;
    case NPT:
      fprintf(FilePtr,"Production run ensemble: NPT (constant number of particles N, constant average pressure P, constant average temperature T)\n");
      break;
    case NPH:
      fprintf(FilePtr,"Production run ensemble: NPH (constant number of particles N, constant average pressure P, constant enthalpy H)\n");
      break;
    case MuPT:
      fprintf(FilePtr,"Production run ensemble: MuPT (constant chemical potential mu, constant average pressure P, constant average temperature T)\n");
      break;
    case NPTPR:
      switch(NPTPRCellType[system])
      {
        case REGULAR:
          fprintf(FilePtr,"Production run ensemble: NPT (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
          break;
        case ANISOTROPIC:
          fprintf(FilePtr,"Production run ensemble: NPT angles fixed, box-lengths are allowed to vary anisotropically (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
          break;          case ISOTROPIC:
          fprintf(FilePtr,"Production run ensemble: NPT angles fixed, box-lengths are allowed to vary isotropically (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPT beta/gamma angles fixed, alpha angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_BETA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPT alpha/gamma angles fixed, beta angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPT alpha/beta angles fixed, gamma angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          fprintf(FilePtr,"Production run ensemble: NPT full cell fluctuations (Parinello-Rahman) using the upper triangle of the box-matrix\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPT beta/gamma angles fixed, alpha angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_BETA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPT alpha/gamma angles fixed, beta angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPT alpha/beta angles fixed, gamma angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
          }
          break;
        default:
          fprintf(stderr, "Undefined cell type for NPTPR\n");
          exit(0);
          break;
      }
      break;
    case NPHPR:
      switch(NPTPRCellType[system])
      {
        case REGULAR:
          fprintf(FilePtr,"Production run ensemble: NPH (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant enthalpy H)\n");
          break;
        case ANISOTROPIC:
          fprintf(FilePtr,"Production run ensemble: NPH angles fixed, box-lengths are allowed to vary anisotropically (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant enthalpy H)\n");
          break;
        case ISOTROPIC:
          fprintf(FilePtr,"Production run ensemble: NPH angles fixed, box-lengths are allowed to vary isotropically (Parinello-Rahman)\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant enthalpy H)\n");
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPH beta/gamma angles fixed, alpha angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_BETA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPH alpha/gamma angles fixed, beta angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPH alpha/beta angles fixed, gamma angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          fprintf(FilePtr,"Production run ensemble: NPH full cell fluctuations (Parinello-Rahman) using the upper triangle of the box-matrix\n");
          fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant enthalpy)\n");
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPH beta/gamma angles fixed, alpha angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_BETA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPH alpha/gamma angles fixed, beta angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              fprintf(FilePtr,"Production run Ensemble: NPH alpha/beta angles fixed, gamma angle and box-lengths are allowed to vary (monoclinic Parinello-Rahman)  using the upper triangle of the box-matrix\n");
              fprintf(FilePtr,"          (constant number of particles N, constant average pressure P, constant average temperature T)\n");
              break;
          }
          break;
        default:
          fprintf(stderr, "Undefined cell type for NPHPR\n");
          exit(0);
          break;
      }
      break;
  }

  fprintf(FilePtr,"\tDegrees of freedom:                        %d\n",DegreesOfFreedom[system]);
  fprintf(FilePtr,"\tTranslational Degrees of freedom:          %d\n",DegreesOfFreedomTranslation[system]);
  fprintf(FilePtr,"\tRotational Degrees of freedom:             %d\n",DegreesOfFreedomRotation[system]);
  fprintf(FilePtr,"\tDegrees of freedom Framework:              %d\n",DegreesOfFreedomFramework[system]);
  fprintf(FilePtr,"\n\n");

  fprintf(FilePtr,"Mutual consistent basic set of units:\n");
  fprintf(FilePtr,"======================================\n");
  if(UseReducedUnits)
    fprintf(FilePtr,"Reduced units\n");
  else
  {
    fprintf(FilePtr,"Unit of temperature: Kelvin\n");
    fprintf(FilePtr,"Unit of length:      %g [m]\n",LENGTH_UNIT);
    fprintf(FilePtr,"Unit of time:        %g [s]\n",TIME_UNIT);
    fprintf(FilePtr,"Unit of mass:        %g [kg]\n",MASS_UNIT);
    fprintf(FilePtr,"Unit of charge:      %g [C/particle]\n",CHARGE_UNIT);
    fprintf(FilePtr,"\n");

    fprintf(FilePtr,"Derived units and their conversion factors:\n");
    fprintf(FilePtr,"===========================================\n");
    fprintf(FilePtr,"Unit of energy:              %g [J]\n",ENERGY_CONVERSION_FACTOR);
    fprintf(FilePtr,"Unit of force:               %g [N]\n",FORCE_CONVERSION_FACTOR);
    fprintf(FilePtr,"Unit of pressure:            %g [Pa]\n",PRESSURE_CONVERSION_FACTOR);
    fprintf(FilePtr,"Unit of velocity:            %g [m/s]\n",VELOCITY_CONVERSION_FACTOR);
    fprintf(FilePtr,"Unit of acceleration:        %g [m^2/s]\n",ACCELERATION_CONVERSION_FACTOR);
    fprintf(FilePtr,"Unit of diffusion:           %g [m^2/s]\n",DIFFUSION_CONVERSION_FACTOR);
    fprintf(FilePtr,"Unit of dipole moment:       %g [C.m]\n",DIPOLE_MOMENT_CONVERSION_FACTOR);
    fprintf(FilePtr,"Unit of electric potential:  %g [V]\n",ELECTRIC_POTENTIAL_CONVERSION_FACTOR);
    fprintf(FilePtr,"Unit of electric field:      %g [V]\n",ELECTRIC_FIELD_CONVERSION_FACTOR);
    fprintf(FilePtr,"Unit of polarizability:      %g [-]\n",POLARIZABILITY);
    fprintf(FilePtr,"Unit of Coulomb potential:   %-18.10f [K]\n",COULOMBIC_CONVERSION_FACTOR*ENERGY_TO_KELVIN);
    fprintf(FilePtr,"Unit of dielectric constant: %-18.10f [s^2 C^2/(kg m^3)]\n",DIELECTRIC_CONSTANT_CONVERSION_FACTOR);
    fprintf(FilePtr,"Unit of wave vectors:        %-18.10f [cm^1]\n",TO_WAVENUMBERS);
    fprintf(FilePtr,"Boltzmann constant:          %-18.10f [-]\n",K_B);
  }
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Internal conversion factors:\n");
  fprintf(FilePtr,"===========================================\n");
  fprintf(FilePtr,"Energy to Kelvin:                              %18.10f\n",ENERGY_TO_KELVIN);
  fprintf(FilePtr,"FH correction factor                           %18.10f\n",FH_CONVERSION_FACTOR);
  fprintf(FilePtr,"Heat capacity conversion factor:               %18.10f\n",HEAT_CAPACITY_CONVERSION_FACTOR);
  fprintf(FilePtr,"From Debye to internal units:                  %18.10f\n",DEBYE_CONVERSION_FACTOR);
  fprintf(FilePtr,"Isothermal compressibility conversion factor:  %18.10f\n",ISOTHERMAL_COMPRESSIBILITY_CONVERSION_FACTOR);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Energy conversion factors:\n");
  fprintf(FilePtr,"===========================================\n");
  fprintf(FilePtr,"From mdyne/A to kcal/mol/A^2:           %g\n",MDYNE_PER_ANGSTROM_TO_KCAL_PER_MOL_PER_ANGSTROM_SQUARED);
  fprintf(FilePtr,"From mdyne/A to kj/mol/A^2:             %g\n",MDYNE_PER_ANGSTROM_TO_KJ_PER_MOL_PER_ANGSTROM_SQUARED);
  fprintf(FilePtr,"From mdyne/A to K/A^2:                  %g\n",MDYNE_PER_ANGSTROM_TO_KELVIN_PER_ANGSTROM_SQUARED);
  fprintf(FilePtr,"From mdyne A/rad^2 to kcal/mol/deg^2:   %g\n",MDYNE_ANGSTROM_PER_RAD_TO_KCAL_PER_MOL_PER_DEGREE_SQUARED);
  fprintf(FilePtr,"\n\n");

  fprintf(FilePtr,"Properties computed\n");
  fprintf(FilePtr,"===========================================================================\n");

  fprintf(FilePtr,"Movies: %s\n",Movies[system]?"yes":"no");
  if(Movies[system])
    fprintf(FilePtr,"\tMovie snapshots are appended to file every %d cycles\n",WriteMoviesEvery[system]);

  // sampling the radial distribution function (RDF)
  fprintf(FilePtr,"Radial Distribution Function: %s\n",ComputeRDF[system]?"yes":"no");
  if(ComputeRDF[system])
  {
    fprintf(FilePtr,"\tRDF is written to file every %d cycles\n",WriteRDFEvery[system]);
    fprintf(FilePtr,"\tUpper limit of the rdf histogram: %f [A]\n",(double)RDFRange[system]);
    fprintf(FilePtr,"\tSize of the RDF histogram: %d\n",RDFHistogramSize[system]);
  }

  // sampling the number-of-molecules histogram
  fprintf(FilePtr,"Number of molecules GCMC histogram: %s\n",ComputeNumberOfMoleculesHistogram[system]?"yes":"no");
  if(ComputeNumberOfMoleculesHistogram[system])
  {
    fprintf(FilePtr,"\tNumber of molecules histograms are written to file every %d cycles\n",WriteNumberOfMoleculesHistogramEvery[system]);
    fprintf(FilePtr,"\tSize of the number-of-molecule histogram: %d\n",NumberOfMoleculesHistogramSize[system]);
    fprintf(FilePtr,"\tRange of the number-of-molecule histogram: %lf\n",NumberOfMoleculesRange[system]);
  }

  // sampling position histograms/free energies
  fprintf(FilePtr,"Histogram of the molecule positions: %s\n",ComputePositionHistogram[system]?"yes":"no");
  if(ComputePositionHistogram[system])
  {
    fprintf(FilePtr,"\tPosition histograms are written to file every %d cycles\n",WritePositionHistogramEvery[system]);
    fprintf(FilePtr,"\tSize of the position histograms: %d\n",PositionHistogramSize[system]);
  }

  // samples the free energy profiles in a,b,c directions
  fprintf(FilePtr,"Free energy profiles: %s\n",ComputeFreeEnergyProfile[system]?"yes":"no");
  if(ComputeFreeEnergyProfile[system])
  {
    fprintf(FilePtr,"\tFree energy profiles are written to file every %d cycles\n",WriteFreeEnergyProfileEvery[system]);
    fprintf(FilePtr,"\tSize of the free energy profile:  %d\n",FreeEnergyHistogramSize[system]);
  }

  // sampling the pore-size distribution (PSD)
  fprintf(FilePtr,"Pore Size Distribution Function: %s\n",((ComputePSDHistogram[system])||(SimulationType==PORE_SIZE_DISTRIBUTION))?"yes":"no");
  if((ComputePSDHistogram[system])||(SimulationType==PORE_SIZE_DISTRIBUTION))
  {
    fprintf(FilePtr,"\tPSD is written to file every %d cycles\n",WritePSDHistogramEvery[system]);
    fprintf(FilePtr,"\tPSD probe distance: %18.10f\n", Framework[system].PoreSizeDistributionProbeDistance);
    fprintf(FilePtr,"\tPSD maximum range: %18.10f\n",PSDRange[system]);
    fprintf(FilePtr,"\tPSD number of elements: %d\n",PSDHistogramSize[system]);
  }

  // sampling the end-to-end histograms
  fprintf(FilePtr,"End-to-end distance: %s\n",ComputeEndToEndDistanceHistogram[system]?"yes":"no");
  if(ComputeEndToEndDistanceHistogram[system])
  {
    fprintf(FilePtr,"\tEnd-to-end distance histograms are written to file every %d cycles\n",WriteEndToEndDistanceHistogramEvery[system]);
    fprintf(FilePtr,"\tEnd-to-end maximum range: %18.10f\n",EndToEndRange[system]);
    fprintf(FilePtr,"\tEnd-to-end number of elements: %d\n",EndToEndHistogramSize[system]);
  }

  // sampling the energy histogram
  fprintf(FilePtr,"Histogram of the energy of the system: %s\n",ComputeEnergyHistogram[system]?"yes":"no");
  if(ComputeEnergyHistogram[system])
  {
    fprintf(FilePtr,"\tEnergy histograms are written to file every %d cycles\n",WriteEnergyHistogramEvery[system]);
    fprintf(FilePtr,"\tLower limit of the energy histogram: %f\n",(double)EnergyHistogramLowerLimit[system]);
    fprintf(FilePtr,"\tUpper limit of the energy histogram: %f\n",(double)EnergyHistogramUpperLimit[system]);
    fprintf(FilePtr,"\tSize of the energy histogram: %d\n",EnergyHistogramSize[system]);
  }

  // sampling the thermodynamic factor
  fprintf(FilePtr,"Compute thermodynamic factors: %s\n",ComputeThermoDynamicFactor[system]?"yes":"no");
  if(ComputeThermoDynamicFactor[system])
    fprintf(FilePtr,"\tThermodynamic factors are written to file every %d cycles\n",WriteThermoDynamicFactorEvery[system]);

  // sampling the inter-framework spacing histogram
  fprintf(FilePtr,"Framework spacing histograms: %s\n",ComputeFrameworkSpacingHistogram[system]?"yes":"no");
  if(ComputeFrameworkSpacingHistogram[system])
  {
    fprintf(FilePtr,"\tFramework spacing histograms are written to file every %d cycles\n",WriteFrameworkSpacingHistogramEvery[system]);
    fprintf(FilePtr,"\tFramework spacing maximum range: %18.10f\n",FrameworkSpacingRange[system]);
    fprintf(FilePtr,"\tFramework spacing number of elements: %d\n",FrameworkSpacingHistogramSize[system]);
  }

  // sampling histograms of the residence times
  fprintf(FilePtr,"Residence times histograms: %s\n",ComputeResidenceTimes[system]?"yes":"no");
  if(ComputeResidenceTimes[system])
  {
    fprintf(FilePtr,"\tresidence-times histogram is written to file every %d cycles\n",WriteResidenceTimesEvery[system]);
    fprintf(FilePtr,"\tnumber of elements of the histogram: %d\n",ResidenceTimesHistogramSize[system]);
    fprintf(FilePtr,"\trange of the histogram: %lg\n",RangeResidenceTimes[system]);
  }

  // sampling histograms of the distance between 2 selected atoms
  fprintf(FilePtr,"Distance histograms: %s\n",ComputeDistanceHistograms[system]?"yes":"no");
  if(ComputeDistanceHistograms[system])
  {
    fprintf(FilePtr,"\tdistance histogram is written to file every %d cycles\n",WriteDistanceHistogramsEvery[system]);
    fprintf(FilePtr,"\tdistance histogram range: %18.10f\n",MaxRangeDistanceHistogram);
    fprintf(FilePtr,"\tdistance histogram number of elements: %d\n",NumberOfElementsDistanceHistogram);

    fprintf(FilePtr,"\tDistance histograms: %d\n",NumberOfDistanceHistogramDefinitions[system]);
    if(NumberOfDistanceHistogramDefinitions[system]>0)
    {
      for(i=0;i<NumberOfDistanceHistogramDefinitions[system];i++)
      {
        switch(DistanceHistogramDefinitions[system][i][0][0])
        {
          case FRAMEWORK:
            sprintf(buffer1,"Framework %d atom %d",DistanceHistogramDefinitions[system][i][0][1],DistanceHistogramDefinitions[system][i][0][2]);
            break;
          case ADSORBATE:
            sprintf(buffer1,"Adsorbate molecule %d atom %d",DistanceHistogramDefinitions[system][i][0][1],DistanceHistogramDefinitions[system][i][0][2]);
            break;
          case CATION:
            sprintf(buffer1,"Cation molecule %d atom %d",DistanceHistogramDefinitions[system][i][0][1],DistanceHistogramDefinitions[system][i][0][2]);
            break;
          default:
            break;
        }
        switch(DistanceHistogramDefinitions[system][i][1][0])
        {
          case FRAMEWORK:
            sprintf(buffer2,"Framework %d atom %d",DistanceHistogramDefinitions[system][i][1][1],DistanceHistogramDefinitions[system][i][1][2]);
            break;
          case ADSORBATE:
            sprintf(buffer2,"Adsorbate molecule %d atom %d",DistanceHistogramDefinitions[system][i][1][1],DistanceHistogramDefinitions[system][i][1][2]);
            break;
          case CATION:
            sprintf(buffer2,"Cation molecule %d atom %d",DistanceHistogramDefinitions[system][i][1][1],DistanceHistogramDefinitions[system][i][1][2]);
            break;
          default:
            break;
        }
        fprintf(FilePtr,"\t\tdistance histogram pair %d (%s,%s)\n",i,buffer1,buffer2);
      }
    }
  }

  // sampling histograms of the bend angle between 3 selected atoms
  fprintf(FilePtr,"Bend Angle histograms: %s\n",ComputeBendAngleHistograms[system]?"yes":"no");
  if(ComputeDistanceHistograms[system])
  {
    fprintf(FilePtr,"\tbend angle histogram is written to file every %d cycles\n",WriteBendAngleHistogramsEvery[system]);
    fprintf(FilePtr,"\tbend angle histogram range: %18.10f\n",MaxRangeBendAngleHistogram);
    fprintf(FilePtr,"\tbend angle histogram number of elements: %d\n",NumberOfElementsBendAngleHistogram);

    fprintf(FilePtr,"\tnumber of bend angles histograms: %d\n",NumberOfBendAngleHistogramDefinitions[system]);
    if(NumberOfBendAngleHistogramDefinitions[system]>0)
    {
      for(i=0;i<NumberOfBendAngleHistogramDefinitions[system];i++)
      {
        switch(BendAngleHistogramDefinitions[system][i][0][0])
        {
          case FRAMEWORK:
            sprintf(buffer1,"Framework %d atom %d",BendAngleHistogramDefinitions[system][i][0][1],BendAngleHistogramDefinitions[system][i][0][2]);
            break;
          case ADSORBATE:
            sprintf(buffer1,"Adsorbate molecule %d atom %d",BendAngleHistogramDefinitions[system][i][0][1],BendAngleHistogramDefinitions[system][i][0][2]);
            break;
          case CATION:
            sprintf(buffer1,"Cation molecule %d atom %d",BendAngleHistogramDefinitions[system][i][0][1],BendAngleHistogramDefinitions[system][i][0][2]);
            break;
          default:
            break;
        }
        switch(BendAngleHistogramDefinitions[system][i][1][0])
        {
          case FRAMEWORK:
            sprintf(buffer2,"Framework %d atom %d",BendAngleHistogramDefinitions[system][i][1][1],BendAngleHistogramDefinitions[system][i][1][2]);
            break;
          case ADSORBATE:
            sprintf(buffer2,"Adsorbate molecule %d atom %d",BendAngleHistogramDefinitions[system][i][1][1],BendAngleHistogramDefinitions[system][i][1][2]);
            break;
          case CATION:
            sprintf(buffer2,"Cation molecule %d atom %d",BendAngleHistogramDefinitions[system][i][1][1],BendAngleHistogramDefinitions[system][i][1][2]);
            break;
          default:
            break;
        }
        switch(BendAngleHistogramDefinitions[system][i][2][0])
        {
          case FRAMEWORK:
            sprintf(buffer3,"Framework %d atom %d",BendAngleHistogramDefinitions[system][i][2][1],BendAngleHistogramDefinitions[system][i][2][2]);
            break;
          case ADSORBATE:
            sprintf(buffer3,"Adsorbate molecule %d atom %d",BendAngleHistogramDefinitions[system][i][2][1],BendAngleHistogramDefinitions[system][i][2][2]);
            break;
          case CATION:
            sprintf(buffer3,"Cation molecule %d atom %d",BendAngleHistogramDefinitions[system][i][2][1],BendAngleHistogramDefinitions[system][i][2][2]);
            break;
          default:
            break;
        }
        fprintf(FilePtr,"\t\tbend angle histogram triple %d (%s,%s,%s)\n",i,buffer1,buffer2,buffer3);
      }
    }
  }

  // sampling histograms of the dihedral angle between 4 selected atoms
  fprintf(FilePtr,"Dihedral angle histograms: %s\n",ComputeDihedralAngleHistograms[system]?"yes":"no");
  if(ComputeDihedralAngleHistograms[system])
  {
    fprintf(FilePtr,"\tdihedral angle histogram is written to file every %d cycles\n",WriteDihedralAngleHistogramsEvery[system]);
    fprintf(FilePtr,"\tdihedral angle histogram range: %18.10f\n",MaxRangeDihedralAngleHistogram);
    fprintf(FilePtr,"\tdihedral angle histogram number of elements: %d\n",NumberOfElementsDihedralAngleHistogram);

    fprintf(FilePtr,"\tnumber of dihedral angle histograms: %d\n",NumberOfDihedralAngleHistogramDefinitions[system]);
    if(NumberOfDihedralAngleHistogramDefinitions[system]>0)
    {
      for(i=0;i<NumberOfDihedralAngleHistogramDefinitions[system];i++)
      {
        switch(DihedralAngleHistogramDefinitions[system][i][0][0])
        {
          case FRAMEWORK:
            sprintf(buffer1,"Framework %d atom %d",DihedralAngleHistogramDefinitions[system][i][0][1],DihedralAngleHistogramDefinitions[system][i][0][2]);
            break;
          case ADSORBATE:
            sprintf(buffer1,"Adsorbate molecule %d atom %d",DihedralAngleHistogramDefinitions[system][i][0][1],DihedralAngleHistogramDefinitions[system][i][0][2]);
            break;
          case CATION:
            sprintf(buffer1,"Cation molecule %d atom %d",DihedralAngleHistogramDefinitions[system][i][0][1],DihedralAngleHistogramDefinitions[system][i][0][2]);
            break;
          default:
            break;
        }
        switch(DihedralAngleHistogramDefinitions[system][i][1][0])
        {
          case FRAMEWORK:
            sprintf(buffer2,"Framework %d atom %d",DihedralAngleHistogramDefinitions[system][i][1][1],DihedralAngleHistogramDefinitions[system][i][1][2]);
            break;
          case ADSORBATE:
            sprintf(buffer2,"Adsorbate molecule %d atom %d",DihedralAngleHistogramDefinitions[system][i][1][1],DihedralAngleHistogramDefinitions[system][i][1][2]);
            break;
          case CATION:
            sprintf(buffer2,"Cation molecule %d atom %d",DihedralAngleHistogramDefinitions[system][i][1][1],DihedralAngleHistogramDefinitions[system][i][1][2]);
            break;
          default:
            break;
        }
        switch(DihedralAngleHistogramDefinitions[system][i][2][0])
        {
          case FRAMEWORK:
            sprintf(buffer3,"Framework %d atom %d",DihedralAngleHistogramDefinitions[system][i][2][1],DihedralAngleHistogramDefinitions[system][i][2][2]);
            break;
          case ADSORBATE:
            sprintf(buffer3,"Adsorbate molecule %d atom %d",DihedralAngleHistogramDefinitions[system][i][2][1],DihedralAngleHistogramDefinitions[system][i][2][2]);
            break;
          case CATION:
            sprintf(buffer3,"Cation molecule %d atom %d",DihedralAngleHistogramDefinitions[system][i][2][1],DihedralAngleHistogramDefinitions[system][i][2][2]);
            break;
          default:
            break;
        }
        switch(DihedralAngleHistogramDefinitions[system][i][3][0])
        {
          case FRAMEWORK:
            sprintf(buffer4,"Framework %d atom %d",DihedralAngleHistogramDefinitions[system][i][3][1],DihedralAngleHistogramDefinitions[system][i][3][2]);
            break;
          case ADSORBATE:
            sprintf(buffer4,"Adsorbate molecule %d atom %d",DihedralAngleHistogramDefinitions[system][i][3][1],DihedralAngleHistogramDefinitions[system][i][3][2]);
            break;
          case CATION:
            sprintf(buffer4,"Cation molecule %d atom %d",DihedralAngleHistogramDefinitions[system][i][3][1],DihedralAngleHistogramDefinitions[system][i][3][2]);
            break;
          default:
            break;
        }
        fprintf(FilePtr,"\t\tdihedral angle histogram quad %d (%s,%s,%s,%s)\n",i,buffer1,buffer2,buffer3,buffer4);
      }
    }
  }

  // sampling histograms of the angle between two planes (each formed by 3 chosen atoms)
  fprintf(FilePtr,"Angle between planes histograms: %s\n",ComputeAngleBetweenPlanesHistograms[system]?"yes":"no");
  if(ComputeAngleBetweenPlanesHistograms[system])
  {
    fprintf(FilePtr,"\tangle between planes histogram is written to file every %d cycles\n",WriteAngleBetweenPlanesHistogramsEvery[system]);
    fprintf(FilePtr,"\tangle between planes histogram range: %18.10f\n",MaxRangeAngleBetweenPlanesHistogram);
    fprintf(FilePtr,"\tangle between planes histogram number of elements: %d\n",NumberOfElementsAngleBetweenPlanesHistogram);

    fprintf(FilePtr,"\tnumber of angle-between-planes histograms: %d\n",NumberOfAngleBetweenPlanesHistogramDefinitions[system]);
    if(NumberOfDistanceHistogramDefinitions[system]>0)
    {
      for(i=0;i<NumberOfAngleBetweenPlanesHistogramDefinitions[system];i++)
      {
        switch(AngleBetweenPlanesHistogramDefinitions[system][i][0][0])
        {
          case FRAMEWORK:
            sprintf(buffer1,"Framework %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][0][1],AngleBetweenPlanesHistogramDefinitions[system][i][0][2]);
            break;
          case ADSORBATE:
            sprintf(buffer1,"Adsorbate molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][0][1],AngleBetweenPlanesHistogramDefinitions[system][i][0][2]);
            break;
          case CATION:
            sprintf(buffer1,"Cation molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][0][1],AngleBetweenPlanesHistogramDefinitions[system][i][0][2]);
            break;
          default:
            break;
        }
        switch(AngleBetweenPlanesHistogramDefinitions[system][i][1][0])
        {
          case FRAMEWORK:
            sprintf(buffer2,"Framework %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][1][1],AngleBetweenPlanesHistogramDefinitions[system][i][1][2]);
            break;
          case ADSORBATE:
            sprintf(buffer2,"Adsorbate molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][1][1],AngleBetweenPlanesHistogramDefinitions[system][i][1][2]);
            break;
          case CATION:
            sprintf(buffer2,"Cation molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][1][1],AngleBetweenPlanesHistogramDefinitions[system][i][1][2]);
            break;
          default:
            break;
        }
        switch(AngleBetweenPlanesHistogramDefinitions[system][i][2][0])
        {
          case FRAMEWORK:
            sprintf(buffer3,"Framework %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][2][1],AngleBetweenPlanesHistogramDefinitions[system][i][2][2]);
            break;
          case ADSORBATE:
            sprintf(buffer3,"Adsorbate molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][2][1],AngleBetweenPlanesHistogramDefinitions[system][i][2][2]);
            break;
          case CATION:
            sprintf(buffer3,"Cation molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][2][1],AngleBetweenPlanesHistogramDefinitions[system][i][2][2]);
            break;
          default:
            break;
        }
        switch(AngleBetweenPlanesHistogramDefinitions[system][i][3][0])
        {
          case FRAMEWORK:
            sprintf(buffer4,"Framework %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][3][1],AngleBetweenPlanesHistogramDefinitions[system][i][3][2]);
            break;
          case ADSORBATE:
            sprintf(buffer4,"Adsorbate molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][3][1],AngleBetweenPlanesHistogramDefinitions[system][i][3][2]);
            break;
          case CATION:
            sprintf(buffer4,"Cation molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][3][1],AngleBetweenPlanesHistogramDefinitions[system][i][3][2]);
            break;
          default:
            break;
        }
        switch(AngleBetweenPlanesHistogramDefinitions[system][i][4][0])
        {
          case FRAMEWORK:
            sprintf(buffer5,"Framework %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][4][1],AngleBetweenPlanesHistogramDefinitions[system][i][4][2]);
            break;
          case ADSORBATE:
            sprintf(buffer5,"Adsorbate molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][4][1],AngleBetweenPlanesHistogramDefinitions[system][i][4][2]);
            break;
          case CATION:
            sprintf(buffer5,"Cation molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][4][1],AngleBetweenPlanesHistogramDefinitions[system][i][4][2]);
            break;
          default:
            break;
        }
        switch(AngleBetweenPlanesHistogramDefinitions[system][i][5][0])
        {
          case FRAMEWORK:
            sprintf(buffer6,"Framework %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][5][1],AngleBetweenPlanesHistogramDefinitions[system][i][5][2]);
            break;
          case ADSORBATE:
            sprintf(buffer6,"Adsorbate molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][5][1],AngleBetweenPlanesHistogramDefinitions[system][i][5][2]);
            break;
          case CATION:
            sprintf(buffer6,"Cation molecule %d atom %d",AngleBetweenPlanesHistogramDefinitions[system][i][5][1],AngleBetweenPlanesHistogramDefinitions[system][i][5][2]);
            break;
          default:
            break;
        }
        fprintf(FilePtr,"\t\tangle between planes histogram %d (%s,%s,%s <-> %s,%s,%s)\n",i,buffer1,buffer2,buffer3,buffer4,buffer5,buffer6);
      }
    }
  }

  // sampling molecular properties (bond distance, bend angle, dihedral angle)
  fprintf(FilePtr,"Molecule properties: %s\n",ComputeMoleculeProperties[system]?"yes":"no");
  if(ComputeMoleculeProperties[system])
  {
    fprintf(FilePtr,"\tMolecular properties are written to file every %d cycles\n",WriteMoleculePropertiesEvery[system]);
    fprintf(FilePtr,"\tNumber of elements of the bond-distance histogram: %d\n",BondLengthHistogramSize[system]);
    fprintf(FilePtr,"\tNumber of elements of the bend-angle histogram: %d\n",BendAngleHistogramSize[system]);
    fprintf(FilePtr,"\tNumber of elements of the dihedral histogram: %d\n",DihedralHistogramSize[system]);
    fprintf(FilePtr,"\tRange of the bond-distance histogram: %lf\n",BondLengthRange[system]);
    fprintf(FilePtr,"\tRange of the bend-angle histogram: %lf\n",BendAngleRange[system]);
    fprintf(FilePtr,"\tRange of the dihedral histogram: %lf\n",DihedralRange[system]);
  }

  // sampling the IR spectra (spacings: 2048, 4196, 8192, 16384, 32768 points)
  fprintf(FilePtr,"Infra-red spectra: %s\n",ComputeInfraRedSpectra[system]?"yes":"no");
  if(ComputeInfraRedSpectra[system])
    fprintf(FilePtr,"\tspectra are written to file every %d cycles\n",WriteInfraRedSpectraEvery[system]);

  // sampling the mean-squared displacement using a modified order-N algorithm
  fprintf(FilePtr,"Mean-squared displacement using modified order-N algorithm: %s\n",ComputeMSDOrderN[system]?"yes":"no");
  if(ComputeMSDOrderN[system])
  {
    fprintf(FilePtr,"\tmsd is sampled every %d cycles\n",SampleMSDOrderNEvery[system]);
    fprintf(FilePtr,"\tmsd is written to file every %d cycles\n",WriteMSDOrderNEvery[system]);
    fprintf(FilePtr,"\tthe (initial) maximum number of blocks: %d\n",MaxNumberOfBlocksMSDOrderN);
    fprintf(FilePtr,"\tthe number of elements per block: %d\n",NumberOfBlockElementsMSDOrderN);
    if(ComputeIndividualMSDOrderN)
      fprintf(FilePtr,"\tthe msd's are computed per component and per molecule\n");
    else
      fprintf(FilePtr,"\tthe msd's are computed per component\n");
  }

  // sampling the velocity autocorrelation function using a modified order-N algorithm
  fprintf(FilePtr,"Velocity-autocorrelation function modified order-N algorithm: %s\n",ComputeVACFOrderN[system]?"yes":"no");
  if(ComputeVACFOrderN[system])
  {
    fprintf(FilePtr,"\tvacf is sampled every %d cycles\n",SampleVACFOrderNEvery[system]);
    fprintf(FilePtr,"\tvacf is written to file every %d cycles\n",WriteVACFOrderNEvery[system]);
    fprintf(FilePtr,"\tthe (initial) maximum number of blocks: %d\n",MaxNumberOfBlocksVACFOrderN);
    fprintf(FilePtr,"\tthe number of elements per block: %d\n",NumberOfBlockElementsVACFOrderN);
    if(ComputeIndividualVACFOrderN)
      fprintf(FilePtr,"\tthe vacf's are computed per component and per molecule\n");
    else
      fprintf(FilePtr,"\tthe vacf's are computed per component\n");
  }

  // sampling of the rotational velocity autocorrelation function using a modified order-N algorithm
  fprintf(FilePtr,"Rotational velocity-autocorrelation function modified order-N algorithm: %s\n",ComputeRVACFOrderN[system]?"yes":"no");
  if(ComputeRVACFOrderN[system])
  {
    fprintf(FilePtr,"\trvacf is sampled every %d cycles\n",SampleRVACFOrderNEvery[system]);
    fprintf(FilePtr,"\trvacf is written to file every %d cycles\n",WriteRVACFOrderNEvery[system]);
    fprintf(FilePtr,"\tthe (initial) maximum number of blocks: %d\n",MaxNumberOfBlocksRVACFOrderN);
    fprintf(FilePtr,"\tthe number of elements per block: %d\n",NumberOfBlockElementsRVACFOrderN);
    if(ComputeIndividualRVACFOrderN)
      fprintf(FilePtr,"\tthe rvacf's are computed per component and per molecule\n");
    else
      fprintf(FilePtr,"\tthe rvacf's are computed per component\n");
  }

  // sampling of the molecular orientation autocorrelation function using a modified order-N algorithm
  fprintf(FilePtr,"Molecular orientation-autocorrelation function modified order-N algorithm: %s\n",ComputeMolecularOrientationOrderN[system]?"yes":"no");
  if(ComputeMolecularOrientationOrderN[system])
  {
    fprintf(FilePtr,"\tmoacf is sampled every %d cycles\n",SampleMolecularOrientationOrderNEvery[system]);
    fprintf(FilePtr,"\tmoacf is written to file every %d cycles\n",WriteMolecularOrientationOrderNEvery[system]);
    fprintf(FilePtr,"\tthe (initial) maximum number of blocks: %d\n",MaxNumberOfBlocksMolecularOrientationOrderN);
    fprintf(FilePtr,"\tthe number of elements per block: %d\n",NumberOfBlockElementsMolecularOrientationOrderN);
    switch(MolecularOrientationType)
    {
      case END_TO_END_VECTOR:
        fprintf(FilePtr,"\tOrientation type: end-to-end vector\n");
        break;
      case MOLECULAR_VECTOR:
        fprintf(FilePtr,"\tOrientation type: fixed vector in molecular frame\n");
        fprintf(FilePtr,"\t\tOrientation group: %d\n",MolecularOrientationGroup);
        fprintf(FilePtr,"\t\tOrientation vector: %lf %lf %lf\n",MolecularOrientationVector.x,MolecularOrientationVector.y,MolecularOrientationVector.z);
        break;
      default:
        break;
    }
  }

  // sampling of the bond orientation autocorrelation function using a modified order-N algorithm
  fprintf(FilePtr,"Bond orientation-autocorrelation function modified order-N algorithm: %s\n",ComputeBondOrientationOrderN[system]?"yes":"no");
  if(ComputeBondOrientationOrderN[system])
  {
    fprintf(FilePtr,"\tboacf is sampled every %d cycles\n",SampleBondOrientationOrderNEvery[system]);
    fprintf(FilePtr,"\tboacf is written to file every %d cycles\n",WriteBondOrientationOrderNEvery[system]);
    fprintf(FilePtr,"\tthe (initial) maximum number of blocks: %d\n",MaxNumberOfBlocksBondOrientationOrderN);
    fprintf(FilePtr,"\tthe number of elements per block: %d\n",NumberOfBlockElementsBondOrientationOrderN);
    for(l=0;l<Framework[system].NumberOfFrameworks;l++)
    {
      fprintf(FilePtr,"\tthe number of sampled orientation bonds: %d\n",
          NumberOfOrientationFrameworkBonds[system][l]);
      for(i=0;i<NumberOfOrientationFrameworkBonds[system][l];i++)
      {
          fprintf(FilePtr,"\t\t framework[%d] bond[%d] => %s %s\n",
              l,i,
              OrientationFrameworkBonds[system][l][i][0],
              OrientationFrameworkBonds[system][l][i][1]);
      }
    }
  }


  // sampling the mean-square displacement function using a conventional algorithm
  fprintf(FilePtr,"Mean-squared displacement (conventional algorithm): %s\n",ComputeMSD[system]?"yes":"no");
  if(ComputeMSD[system])
  {
    fprintf(FilePtr,"\tmsd is sampled every %d cycles\n",SampleMSDEvery[system]);
    fprintf(FilePtr,"\tmsd is written to file every %d cycles\n",WriteMSDEvery[system]);
  }

  // sampling the velocity autocorrelation function using a conventional algorithm
  fprintf(FilePtr,"Velocity-autocorrelation function (conventional algorithm): %s\n",ComputeVACF[system]?"yes":"no");
  if(ComputeVACF[system])
  {
    fprintf(FilePtr,"\tvacf is sampled every %d cycles\n",SampleVACFEvery[system]);
    fprintf(FilePtr,"\tvacf is written to file every %d cycles\n",WriteVACFEvery[system]);
  }

  // sampling the 3D histograms of position (i.e. 3D free energy)
  fprintf(FilePtr,"3D density grid for adsorbates: %s\n",ComputeDensityProfile3DVTKGrid[system]?"yes":"no");
  if(ComputeDensityProfile3DVTKGrid[system])
  {
    fprintf(FilePtr,"\t3D density grids for adsorbates are written to file every %d cycles\n",WriteDensityProfile3DVTKGridEvery[system]);
    fprintf(FilePtr,"\t3D density grids consisting of %d x %d x %d points\n",DensityProfile3DVTKGridPoints.x,DensityProfile3DVTKGridPoints.y,DensityProfile3DVTKGridPoints.z);
    if(DensityAveragingTypeVTK==VTK_UNIT_CELL) fprintf(FilePtr,"\t3D density grid output as a single unit cell\n");
    else fprintf(FilePtr,"\t3D density grid made for the full simulation-cell\n");
    if(DensityAveragingTypeVTK==VTK_FULL_BOX)
    {
      if(AverageDensityOverUnitCellsVTK==TRUE)
        fprintf(FilePtr,"\t3D density grid averaged over all unit cells\n");
      else
        fprintf(FilePtr,"\t3D density grid NOT averaged over all unit cells\n");
    }
  }

  // samples the cation sites and adsorption sites
  fprintf(FilePtr,"Compute cation an/or adsorption sites: %s\n",ComputeCationAndAdsorptionSites[system]?"yes":"no");
  if(ComputeCationAndAdsorptionSites[system])
    fprintf(FilePtr,"\tCation and/or adsorption sites are written to file every %d cycles\n",WriteCationAndAdsorptionSitesEvery[system]);

  // samples initial configurations for the tranmission coefficient (dcTST)
  fprintf(FilePtr,"dcTST snapshots: %s\n",WritedcTSTSnapShotsToFile[system]?"yes":"no");
  if(WritedcTSTSnapShotsToFile[system])
    fprintf(FilePtr,"\tdcTST snapshots are written to file every %d cycles\n",WritedcTSTSnapShotsEvery[system]);

  // samples the stress and pressure
  fprintf(FilePtr,"Compute pressure and stress: %s\n",ComputeMolecularPressure[system]?"yes":"no");
  fprintf(FilePtr,"\n\n");

  fprintf(FilePtr,"VTK\n");
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"VTK fractional-range position framework atoms: [%lf,%lf] [%lf,%lf] [%lf,%lf]\n",(double)VTKFractionalFrameworkAtomsMin.x,
          (double)VTKFractionalFrameworkAtomsMax.x,(double)VTKFractionalFrameworkAtomsMin.x,(double)VTKFractionalFrameworkAtomsMax.y,
          (double)VTKFractionalFrameworkAtomsMin.z,(double)VTKFractionalFrameworkAtomsMax.z);
  fprintf(FilePtr,"VTK fractional-range position framework bonds: [%lf,%lf] [%lf,%lf] [%lf,%lf]\n",(double)VTKFractionalFrameworkBondsMin.x,
          (double)VTKFractionalFrameworkBondsMax.x,(double)VTKFractionalFrameworkBondsMin.x,(double)VTKFractionalFrameworkBondsMax.y,
          (double)VTKFractionalFrameworkBondsMin.z,(double)VTKFractionalFrameworkBondsMax.z);
  fprintf(FilePtr,"VTK fractional-range com-position adsorbate molecules: [%lf,%lf] [%lf,%lf] [%lf,%lf]\n",(double)VTKFractionalAdsorbateComMin.x,
          (double)VTKFractionalAdsorbateComMax.x,(double)VTKFractionalAdsorbateComMin.x,(double)VTKFractionalAdsorbateComMax.y,
          (double)VTKFractionalAdsorbateComMin.z,(double)VTKFractionalAdsorbateComMax.z);
  fprintf(FilePtr,"VTK fractional-range com-position cation molecules: [%lf,%lf] [%lf,%lf] [%lf,%lf]\n",(double)VTKFractionalCationComMin.x,
          (double)VTKFractionalCationComMax.x,(double)VTKFractionalCationComMin.x,(double)VTKFractionalCationComMax.y,
          (double)VTKFractionalCationComMin.z,(double)VTKFractionalCationComMax.z);
  if(FreeEnergyAveragingTypeVTK==VTK_UNIT_CELL) fprintf(FilePtr,"\t3D free energy grid output as a single unit cell\n");
  else fprintf(FilePtr,"\t3D free energy grid made for the full simulation-cell\n");
  fprintf(FilePtr,"\n\n");

  fprintf(FilePtr,"Thermo/Baro-stat NHC parameters\n");
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"External temperature: %lg [K]\n",(double)therm_baro_stats.ExternalTemperature[system]);
  fprintf(FilePtr,"Beta: %lg [energy unit]\n",(double)Beta[system]);
  if(NumberOfIsothermPressures==1)
     fprintf(FilePtr,"External Pressure: %lg [Pa]\n",(double)therm_baro_stats.ExternalPressure[system][0]*PRESSURE_CONVERSION_FACTOR);
  else
  {
    fprintf(FilePtr,"Number of isotherm points: %d\n",NumberOfIsothermPressures);
    for(i=0;i<NumberOfIsothermPressures;i++)
    {
      if(i==CurrentIsothermPressure)
        fprintf(FilePtr,"\t\t point %3d: External Pressure: %lg [Pa]\n",i,(double)therm_baro_stats.ExternalPressure[system][i]*PRESSURE_CONVERSION_FACTOR);
      else
        fprintf(FilePtr,"\t point %3d: External Pressure: %lg [Pa]\n",i,(double)therm_baro_stats.ExternalPressure[system][i]*PRESSURE_CONVERSION_FACTOR);
    }
  }
  if(therm_baro_stats.UseExternalStress)
  {
    switch(Dimension)
    {
      case 2:
        fprintf(FilePtr,"External Stress: {{%f,%f},{-,%lf}} [Pa]\n",
          (double)therm_baro_stats.ExternalStress[system].ax*PRESSURE_CONVERSION_FACTOR,
          (double)therm_baro_stats.ExternalStress[system].ay*PRESSURE_CONVERSION_FACTOR,
          (double)therm_baro_stats.ExternalStress[system].by*PRESSURE_CONVERSION_FACTOR);
        break;
      case 3:
        fprintf(FilePtr,"External Stress: {{%lf,%lf,%lf},{-,%f,%f},{-,-,%lf}} [Pa]\n",
          (double)therm_baro_stats.ExternalStress[system].ax*PRESSURE_CONVERSION_FACTOR,
          (double)therm_baro_stats.ExternalStress[system].ay*PRESSURE_CONVERSION_FACTOR,
          (double)therm_baro_stats.ExternalStress[system].az*PRESSURE_CONVERSION_FACTOR,
          (double)therm_baro_stats.ExternalStress[system].by*PRESSURE_CONVERSION_FACTOR,
          (double)therm_baro_stats.ExternalStress[system].bz*PRESSURE_CONVERSION_FACTOR,
          (double)therm_baro_stats.ExternalStress[system].cz*PRESSURE_CONVERSION_FACTOR);
        break;
    }
  }
  fprintf(FilePtr,"\n\n");


  fprintf(FilePtr,"Thermostat chain-length: %d\n",therm_baro_stats.ThermostatChainLength);
  fprintf(FilePtr,"Timescale parameter for thermostat: %lf [ps]\n",(double)therm_baro_stats.time_scale_parameter_thermostat);
  fprintf(FilePtr,"Barostat chain-length:   %d\n",therm_baro_stats.BarostatChainLength);
  fprintf(FilePtr,"Timescale parameter for barostat:   %lf [ps]\n",(double)therm_baro_stats.time_scale_parameter_barostat);
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Number of Yoshida-Suzuki decomposition steps: %d\n",therm_baro_stats.NumberOfYoshidaSuzukiSteps);
  fprintf(FilePtr,"Number of respa steps: %d\n",therm_baro_stats.NumberOfRespaSteps);
  fprintf(FilePtr,"\n\n");

  fprintf(FilePtr,"Method and settings for electrostatics\n");
  fprintf(FilePtr,"===============================================================================\n");
  fprintf(FilePtr,"Dielectric constant of the medium : %lf\n",(double)DielectricConstantOfTheMedium);
  fprintf(FilePtr,"Charge from charge-equilibration: %s\n",ChargeFromChargeEquilibration?"yes":"no");
  if(ChargeFromChargeEquilibration)
  {
    fprintf(FilePtr,"Symmetrize the charges of the framework: %s\n",SymmetrizeFrameworkCharges?"yes":"no");
  }
  switch(ChargeMethod)
  {
    case NONE:
      fprintf(FilePtr,"No electrostatics\n");
      break;
    case TRUNCATED_COULOMB:
      if(BoundaryCondition[system]==FINITE)
        fprintf(FilePtr,"Coulombic potential is used (finite system)\n");
      else
        fprintf(FilePtr,"Coulombic potential is truncated at: %f [A]\n",CutOffChargeCharge[system]);
      break;
    case SHIFTED_COULOMB:
      fprintf(FilePtr,"Coulombic potential is truncated and shifted to zero at: %f [A]\n",CutOffChargeCharge[system]);
      break;
    case SMOOTHED_COULOMB:
      fprintf(FilePtr,"Coulombic potential is smoothed using a switching function from %f to: %f [A]\n",CutOffChargeChargeSwitch[system],CutOffChargeCharge[system]);
      break;
    case WOLFS_METHOD:
      fprintf(FilePtr,"The original Wolf's method is used, potential shifted to zero at %f [A]\n",CutOffChargeCharge[system]);
      break;
    case WOLFS_METHOD_DAMPED:
      fprintf(FilePtr,"The damped Wolf's method is used, potential shifted to zero at %f [A], Alpha: %18.10f\n",CutOffChargeCharge[system],Alpha[system]);
      break;
    case WOLFS_METHOD_DAMPED_FG:
      fprintf(FilePtr,"The damped potential of Fenell and Gezelter is used, potential shifted to zero at %f [A], Alpha: %18.10f\n",CutOffChargeCharge[system],Alpha[system]);
      break;
    case EWALD:
    default:
      fprintf(FilePtr,"Ewald summation is used (exact solution of a periodic system)\n");
      if(EwaldAutomatic)
        fprintf(FilePtr,"Relative precision                : %g\n",(double)EwaldPrecision);
      fprintf(FilePtr,"Alpha convergence parameter       : %lf\n",(double)Alpha[system]);
      fprintf(FilePtr,"kvec (x,y,z)                      : %d %d %d\n",kvec[system].x,kvec[system].y,kvec[system].z);
      break;
  }
  fprintf(FilePtr,"\n\n");

  fprintf(FilePtr,"CFC-RXMC parameters\n");
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"Number of reactions: %d\n",NumberOfReactions);
  for(i=0;i<NumberOfReactions;i++)
  {
    fprintf(FilePtr,"\tReaction [%d]:",i);
    for(j=0;j<NumberOfComponents;j++)
      fprintf(FilePtr," %d",ReactantsStoichiometry[i][j]);
    fprintf(FilePtr," -->");
    for(j=0;j<NumberOfComponents;j++)
      fprintf(FilePtr," %d",ProductsStoichiometry[i][j]);
    fprintf(FilePtr,"\n");
  }
  fprintf(FilePtr,"\n\n");

  fprintf(FilePtr,"Rattle parameters\n");
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"Distance constraint type: %s\n",DistanceConstraintType==DISTANCE_R?"r-r_0":"r^2-r^2_0");
  fprintf(FilePtr,"Bend angle constraint type: %s\n",BendConstraintType==THETA?"theta-theta_0":
                   (BendConstraintType==COS_THETA?"cos(theta)-cos(theta_0)":"cos^2(theta)-cos^2(theta_0)"));
  fprintf(FilePtr,"Dihedral angle constraint type: %s\n",DihedralConstraintType==PHI?"phi-phi_0":
                   (DihedralConstraintType==COS_PHI?"cos(phi)-cos(phi_0)":"cos^2(phi)-cos^2(phi_0)"));
  fprintf(FilePtr,"Inversion-bend angle constraint type: %s\n",InversionBendConstraintType==CHI?"chi-chi_0":
                   (InversionBendConstraintType==SIN_CHI?"cos(phi)-cos(phi_0)":"cos^2(chi)-cos^2(chi_0)"));
  fprintf(FilePtr,"Out-of-plane distance constraint type: %s\n",OutOfPlaneConstraintType==DISTANCE_R?"r-r_0":"r^2-r^2_0");
  fprintf(FilePtr,"\n\n");

  fprintf(FilePtr,"Spectra parameters\n");
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"Compute normal modes: %s\n",ComputeNormalModes==TRUE?"yes":"no");
  if(ComputeNormalModes)
  {
    fprintf(FilePtr,"number of steps of the normal mode: %d\n",ModeResolution);
    fprintf(FilePtr,"from normal mode: %d to: %d\n",MinimumMode,MaximumMode);
    fprintf(FilePtr,"Correct normal modes for constraints: %s\n",CorrectNormalModesForConstraints==TRUE?"yes":"no");
  }
  fprintf(FilePtr,"\n\n");

  fprintf(FilePtr,"Minimization parameters\n");
  fprintf(FilePtr,"===========================================================================\n");

  if(MinimizationVariables==FRACTIONAL)
    fprintf(FilePtr,"Generalized coordinates are: fractional center-of-mass, elements of the orientational matrix p1,p2,p3 and strain\n");
  else
    fprintf(FilePtr,"Generalized coordinates are: Cartesian center-of-mass, elements of the orientational matrix p1,p2,p3 and strain\n");
  switch(MinimizationPotentialMethod)
  {
    case NUMERICALLY:
      fprintf(FilePtr,"Potential derivatives are evaluated: numerically\n");
      break;
    case ANALYTICALLY:
    default:
      fprintf(FilePtr,"Potential derivatives are evaluated: analytically\n");
      break;
  }
  fprintf(FilePtr,"Translation of the system is removed from the generalized Hessian: %s\n",RemoveTranslationFromHessian?"yes":"no");
  fprintf(FilePtr,"Rotation of the system is removed from the generalized Hessian: %s\n",RemoveRotationFromHessian?"yes":"no");
  fprintf(FilePtr,"Maximum step-length: %g\n",MaximumStepLengthInput);
  fprintf(FilePtr,"Convergence factor: %g\n",MinimizationConvergenceFactor);
  fprintf(FilePtr,"Maximum number of minimization steps: %d\n",MaximumNumberOfMinimizationSteps);
  fprintf(FilePtr,"Use gradients in the line-minimizations: %s\n",UseGradientInLineMinimization?"yes":"no");
  fprintf(FilePtr,"RMS gradient tolerance: %g\n",RMSGradientTolerance);
  fprintf(FilePtr,"Maximum gradient tolerance: %g\n",MaxGradientTolerance);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Distance constraints: %d\n",NumberOfDistanceConstraints[system]);
  if(NumberOfDistanceConstraints[system]>0)
  {
    for(i=0;i<NumberOfDistanceConstraints[system];i++)
    {
      switch(DistanceDefinitions[system][i][0][0])
      {
        case FRAMEWORK:
          sprintf(buffer1,"Framework %d atom %d",DistanceDefinitions[system][i][0][1],DistanceDefinitions[system][i][0][2]);
          break;
        case ADSORBATE:
          sprintf(buffer1,"Adsorbate molecule %d atom %d",DistanceDefinitions[system][i][0][1],DistanceDefinitions[system][i][0][2]);
          break;
        case CATION:
          sprintf(buffer1,"Cation molecule %d atom %d",DistanceDefinitions[system][i][0][1],DistanceDefinitions[system][i][0][2]);
          break;
        default:
          break;
      }
      switch(DistanceDefinitions[system][i][1][0])
      {
        case FRAMEWORK:
          sprintf(buffer2,"Framework %d atom %d",DistanceDefinitions[system][i][1][1],DistanceDefinitions[system][i][1][2]);
          break;
        case ADSORBATE:
          sprintf(buffer2,"Adsorbate molecule %d atom %d",DistanceDefinitions[system][i][1][1],DistanceDefinitions[system][i][1][2]);
          break;
        case CATION:
          sprintf(buffer2,"Cation molecule %d atom %d",DistanceDefinitions[system][i][1][1],DistanceDefinitions[system][i][1][2]);
          break;
        default:
          break;
      }
      fprintf(FilePtr,"\tconstraint %d (%s,%s)  constraint distance: %g [A]\n",i,buffer1,buffer2,DistanceConstraintParameter[system][i]);
    }
  }

  fprintf(FilePtr,"Angle constraints: %d\n",NumberOfAngleConstraints[system]);
  if(NumberOfAngleConstraints[system]>0)
  {
    for(i=0;i<NumberOfAngleConstraints[system];i++)
    {
      switch(AngleDefinitions[system][i][0][0])
      {
        case FRAMEWORK:
          sprintf(buffer1,"Framework %d atom %d",AngleDefinitions[system][i][0][1],AngleDefinitions[system][i][0][2]);
          break;
        case ADSORBATE:
          sprintf(buffer1,"Adsorbate molecule %d atom %d",AngleDefinitions[system][i][0][1],AngleDefinitions[system][i][0][2]);
          break;
        case CATION:
          sprintf(buffer1,"Cation molecule %d atom %d",AngleDefinitions[system][i][0][1],AngleDefinitions[system][i][0][2]);
          break;
        default:
          break;
      }
      switch(AngleDefinitions[system][i][1][0])
      {
        case FRAMEWORK:
          sprintf(buffer2,"Framework %d atom %d",AngleDefinitions[system][i][1][1],AngleDefinitions[system][i][1][2]);
          break;
        case ADSORBATE:
          sprintf(buffer2,"Adsorbate molecule %d atom %d",AngleDefinitions[system][i][1][1],AngleDefinitions[system][i][1][2]);
          break;
        case CATION:
          sprintf(buffer2,"Cation molecule %d atom %d",AngleDefinitions[system][i][1][1],AngleDefinitions[system][i][1][2]);
          break;
        default:
          break;
      }
      switch(AngleDefinitions[system][i][2][0])
      {
        case FRAMEWORK:
          sprintf(buffer3,"Framework %d atom %d",AngleDefinitions[system][i][2][1],AngleDefinitions[system][i][2][2]);
          break;
        case ADSORBATE:
          sprintf(buffer3,"Adsorbate molecule %d atom %d",AngleDefinitions[system][i][2][1],AngleDefinitions[system][i][2][2]);
          break;
        case CATION:
          sprintf(buffer3,"Cation molecule %d atom %d",AngleDefinitions[system][i][2][1],AngleDefinitions[system][i][2][2]);
          break;
        default:
          break;
      }

      fprintf(FilePtr,"\tconstraint %d (%s,%s,%s)  constraint angle: %g [degrees]\n",i,buffer1,buffer2,buffer3,AngleConstraintParameter[system][i]*RAD2DEG);
    }
  }

  fprintf(FilePtr,"Dihedral constraints: %d\n",NumberOfDihedralConstraints[system]);
  if(NumberOfDihedralConstraints[system]>0)
  {
    for(i=0;i<NumberOfDihedralConstraints[system];i++)
    {
      switch(DihedralDefinitions[system][i][0][0])
      {
        case FRAMEWORK:
          sprintf(buffer1,"Framework %d atom %d",DihedralDefinitions[system][i][0][1],DihedralDefinitions[system][i][0][2]);
          break;
        case ADSORBATE:
          sprintf(buffer1,"Adsorbate molecule %d atom %d",DihedralDefinitions[system][i][0][1],DihedralDefinitions[system][i][0][2]);
          break;
        case CATION:
          sprintf(buffer1,"Cation molecule %d atom %d",DihedralDefinitions[system][i][0][1],DihedralDefinitions[system][i][0][2]);
          break;
        default:
          break;
      }
      switch(DihedralDefinitions[system][i][1][0])
      {
        case FRAMEWORK:
          sprintf(buffer2,"Framework %d atom %d",DihedralDefinitions[system][i][1][1],DihedralDefinitions[system][i][1][2]);
          break;
        case ADSORBATE:
          sprintf(buffer2,"Adsorbate molecule %d atom %d",DihedralDefinitions[system][i][1][1],DihedralDefinitions[system][i][1][2]);
          break;
        case CATION:
          sprintf(buffer2,"Cation molecule %d atom %d",DihedralDefinitions[system][i][1][1],DihedralDefinitions[system][i][1][2]);
          break;
        default:
          break;
      }
      switch(DihedralDefinitions[system][i][2][0])
      {
        case FRAMEWORK:
          sprintf(buffer3,"Framework %d atom %d",DihedralDefinitions[system][i][2][1],DihedralDefinitions[system][i][2][2]);
          break;
        case ADSORBATE:
          sprintf(buffer3,"Adsorbate molecule %d atom %d",DihedralDefinitions[system][i][2][1],DihedralDefinitions[system][i][2][2]);
          break;
        case CATION:
          sprintf(buffer3,"Cation molecule %d atom %d",DihedralDefinitions[system][i][2][1],DihedralDefinitions[system][i][2][2]);
          break;
        default:
          break;
      }
      switch(DihedralDefinitions[system][i][3][0])
      {
        case FRAMEWORK:
          sprintf(buffer4,"Framework %d atom %d",DihedralDefinitions[system][i][3][1],DihedralDefinitions[system][i][3][2]);
          break;
        case ADSORBATE:
          sprintf(buffer4,"Adsorbate molecule %d atom %d",DihedralDefinitions[system][i][3][1],DihedralDefinitions[system][i][3][2]);
          break;
        case CATION:
          sprintf(buffer4,"Cation molecule %d atom %d",DihedralDefinitions[system][i][3][1],DihedralDefinitions[system][i][3][2]);
          break;
        default:
          break;
      }
      fprintf(FilePtr,"\tconstraint %d (%s,%s,%s,%s)  constraint angle: %g [degrees]\n",i,buffer1,buffer2,buffer3,buffer4,DihedralConstraintParameter[system][i]*RAD2DEG);
    }
  }
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Improper dihedral constraints: %d\n",NumberOfImproperDihedralConstraints[system]);
  if(NumberOfImproperDihedralConstraints[system]>0)
  {
    for(i=0;i<NumberOfImproperDihedralConstraints[system];i++)
    {
      switch(ImproperDihedralDefinitions[system][i][0][0])
      {
        case FRAMEWORK:
          sprintf(buffer1,"Framework %d atom %d",ImproperDihedralDefinitions[system][i][0][1],ImproperDihedralDefinitions[system][i][0][2]);
          break;
        case ADSORBATE:
          sprintf(buffer1,"Adsorbate molecule %d atom %d",ImproperDihedralDefinitions[system][i][0][1],ImproperDihedralDefinitions[system][i][0][2]);
          break;
        case CATION:
          sprintf(buffer1,"Cation molecule %d atom %d",ImproperDihedralDefinitions[system][i][0][1],ImproperDihedralDefinitions[system][i][0][2]);
          break;
        default:
          break;
      }
      switch(ImproperDihedralDefinitions[system][i][1][0])
      {
        case FRAMEWORK:
          sprintf(buffer2,"Framework %d atom %d",ImproperDihedralDefinitions[system][i][1][1],ImproperDihedralDefinitions[system][i][1][2]);
          break;
        case ADSORBATE:
          sprintf(buffer2,"Adsorbate molecule %d atom %d",ImproperDihedralDefinitions[system][i][1][1],ImproperDihedralDefinitions[system][i][1][2]);
          break;
        case CATION:
          sprintf(buffer2,"Cation molecule %d atom %d",ImproperDihedralDefinitions[system][i][1][1],ImproperDihedralDefinitions[system][i][1][2]);
          break;
        default:
          break;
      }
      switch(ImproperDihedralDefinitions[system][i][2][0])
      {
        case FRAMEWORK:
          sprintf(buffer3,"Framework %d atom %d",ImproperDihedralDefinitions[system][i][2][1],ImproperDihedralDefinitions[system][i][2][2]);
          break;
        case ADSORBATE:
          sprintf(buffer3,"Adsorbate molecule %d atom %d",ImproperDihedralDefinitions[system][i][2][1],ImproperDihedralDefinitions[system][i][2][2]);
          break;
        case CATION:
          sprintf(buffer3,"Cation molecule %d atom %d",ImproperDihedralDefinitions[system][i][2][1],ImproperDihedralDefinitions[system][i][2][2]);
          break;
        default:
          break;
      }
      switch(ImproperDihedralDefinitions[system][i][3][0])
      {
        case FRAMEWORK:
          sprintf(buffer4,"Framework %d atom %d",ImproperDihedralDefinitions[system][i][3][1],ImproperDihedralDefinitions[system][i][3][2]);
          break;
        case ADSORBATE:
          sprintf(buffer4,"Adsorbate molecule %d atom %d",ImproperDihedralDefinitions[system][i][3][1],ImproperDihedralDefinitions[system][i][3][2]);
          break;
        case CATION:
          sprintf(buffer4,"Cation molecule %d atom %d",ImproperDihedralDefinitions[system][i][3][1],ImproperDihedralDefinitions[system][i][3][2]);
          break;
        default:
          break;
      }
      fprintf(FilePtr,"\tconstraint %d (%s,%s,%s,%s)  constraint angle: %g [degrees]\n",i,buffer1,buffer2,buffer3,buffer4,ImproperDihedralConstraintParameter[system][i]*RAD2DEG);
    }
  }
  fprintf(FilePtr,"\n");


  fprintf(FilePtr,"Inversion-bend constraints: %d\n",NumberOfInversionBendConstraints[system]);
  if(NumberOfInversionBendConstraints[system]>0)
  {
    for(i=0;i<NumberOfInversionBendConstraints[system];i++)
    {
      switch(InversionBendDefinitions[system][i][0][0])
      {
        case FRAMEWORK:
          sprintf(buffer1,"Framework %d atom %d",InversionBendDefinitions[system][i][0][1],InversionBendDefinitions[system][i][0][2]);
          break;
        case ADSORBATE:
          sprintf(buffer1,"Adsorbate molecule %d atom %d",InversionBendDefinitions[system][i][0][1],InversionBendDefinitions[system][i][0][2]);
          break;
        case CATION:
          sprintf(buffer1,"Cation molecule %d atom %d",InversionBendDefinitions[system][i][0][1],InversionBendDefinitions[system][i][0][2]);
          break;
        default:
          break;
      }
      switch(InversionBendDefinitions[system][i][1][0])
      {
        case FRAMEWORK:
          sprintf(buffer2,"Framework %d atom %d",InversionBendDefinitions[system][i][1][1],InversionBendDefinitions[system][i][1][2]);
          break;
        case ADSORBATE:
          sprintf(buffer2,"Adsorbate molecule %d atom %d",InversionBendDefinitions[system][i][1][1],InversionBendDefinitions[system][i][1][2]);
          break;
        case CATION:
          sprintf(buffer2,"Cation molecule %d atom %d",InversionBendDefinitions[system][i][1][1],InversionBendDefinitions[system][i][1][2]);
          break;
        default:
          break;
      }
      switch(InversionBendDefinitions[system][i][2][0])
      {
        case FRAMEWORK:
          sprintf(buffer3,"Framework %d atom %d",InversionBendDefinitions[system][i][2][1],InversionBendDefinitions[system][i][2][2]);
          break;
        case ADSORBATE:
          sprintf(buffer3,"Adsorbate molecule %d atom %d",InversionBendDefinitions[system][i][2][1],InversionBendDefinitions[system][i][2][2]);
          break;
        case CATION:
          sprintf(buffer3,"Cation molecule %d atom %d",InversionBendDefinitions[system][i][2][1],InversionBendDefinitions[system][i][2][2]);
          break;
        default:
          break;
      }
      switch(InversionBendDefinitions[system][i][3][0])
      {
        case FRAMEWORK:
          sprintf(buffer4,"Framework %d atom %d",InversionBendDefinitions[system][i][3][1],InversionBendDefinitions[system][i][3][2]);
          break;
        case ADSORBATE:
          sprintf(buffer4,"Adsorbate molecule %d atom %d",InversionBendDefinitions[system][i][3][1],InversionBendDefinitions[system][i][3][2]);
          break;
        case CATION:
          sprintf(buffer4,"Cation molecule %d atom %d",InversionBendDefinitions[system][i][3][1],InversionBendDefinitions[system][i][3][2]);
          break;
        default:
          break;
      }
      fprintf(FilePtr,"\tconstraint %d (%s,%s,%s,%s)  constraint angle: %g [degrees]\n",i,buffer1,buffer2,buffer3,buffer4,InversionBendConstraintParameter[system][i]*RAD2DEG);
    }
  }
  fprintf(FilePtr,"\n");


  fprintf(FilePtr,"Out-of-plane constraints: %d\n",NumberOfOutOfPlaneDistanceConstraints[system]);
  if(NumberOfOutOfPlaneDistanceConstraints[system]>0)
  {
    for(i=0;i<NumberOfOutOfPlaneDistanceConstraints[system];i++)
    {
      switch(OutOfPlaneDistanceDefinitions[system][i][0][0])
      {
        case FRAMEWORK:
          sprintf(buffer1,"Framework %d atom %d",OutOfPlaneDistanceDefinitions[system][i][0][1],OutOfPlaneDistanceDefinitions[system][i][0][2]);
          break;
        case ADSORBATE:
          sprintf(buffer1,"Adsorbate molecule %d atom %d",OutOfPlaneDistanceDefinitions[system][i][0][1],OutOfPlaneDistanceDefinitions[system][i][0][2]);
          break;
        case CATION:
          sprintf(buffer1,"Cation molecule %d atom %d",OutOfPlaneDistanceDefinitions[system][i][0][1],OutOfPlaneDistanceDefinitions[system][i][0][2]);
          break;
        default:
          break;
      }
      switch(OutOfPlaneDistanceDefinitions[system][i][1][0])
      {
        case FRAMEWORK:
          sprintf(buffer2,"Framework %d atom %d",OutOfPlaneDistanceDefinitions[system][i][1][1],OutOfPlaneDistanceDefinitions[system][i][1][2]);
          break;
        case ADSORBATE:
          sprintf(buffer2,"Adsorbate molecule %d atom %d",OutOfPlaneDistanceDefinitions[system][i][1][1],OutOfPlaneDistanceDefinitions[system][i][1][2]);
          break;
        case CATION:
          sprintf(buffer2,"Cation molecule %d atom %d",OutOfPlaneDistanceDefinitions[system][i][1][1],OutOfPlaneDistanceDefinitions[system][i][1][2]);
          break;
        default:
          break;
      }
      switch(OutOfPlaneDistanceDefinitions[system][i][2][0])
      {
        case FRAMEWORK:
          sprintf(buffer3,"Framework %d atom %d",OutOfPlaneDistanceDefinitions[system][i][2][1],OutOfPlaneDistanceDefinitions[system][i][2][2]);
          break;
        case ADSORBATE:
          sprintf(buffer3,"Adsorbate molecule %d atom %d",OutOfPlaneDistanceDefinitions[system][i][2][1],OutOfPlaneDistanceDefinitions[system][i][2][2]);
          break;
        case CATION:
          sprintf(buffer3,"Cation molecule %d atom %d",OutOfPlaneDistanceDefinitions[system][i][2][1],OutOfPlaneDistanceDefinitions[system][i][2][2]);
          break;
        default:
          break;
      }
      switch(OutOfPlaneDistanceDefinitions[system][i][3][0])
      {
        case FRAMEWORK:
          sprintf(buffer4,"Framework %d atom %d",OutOfPlaneDistanceDefinitions[system][i][3][1],OutOfPlaneDistanceDefinitions[system][i][3][2]);
          break;
        case ADSORBATE:
          sprintf(buffer4,"Adsorbate molecule %d atom %d",OutOfPlaneDistanceDefinitions[system][i][3][1],OutOfPlaneDistanceDefinitions[system][i][3][2]);
          break;
        case CATION:
          sprintf(buffer4,"Cation molecule %d atom %d",OutOfPlaneDistanceDefinitions[system][i][3][1],OutOfPlaneDistanceDefinitions[system][i][3][2]);
          break;
        default:
          break;
      }
      fprintf(FilePtr,"\tconstraint %d (%s,%s,%s,%s)  constraint angle: %g [degrees]\n",i,buffer1,buffer2,buffer3,buffer4,OutOfPlaneDistanceConstraintParameter[system][i]*RAD2DEG);
    }
  }
  fprintf(FilePtr,"\n");



  fprintf(FilePtr,"Harmonic distance constraints: %d\n",NumberOfHarmonicDistanceConstraints[system]);
  if(NumberOfHarmonicDistanceConstraints[system]>0)
  {
    for(i=0;i<NumberOfHarmonicDistanceConstraints[system];i++)
    {
      switch(HarmonicDistanceDefinitions[system][i][0][0])
      {
        case FRAMEWORK:
          sprintf(buffer1,"Framework %d atom %d",HarmonicDistanceDefinitions[system][i][0][1],HarmonicDistanceDefinitions[system][i][0][2]);
          break;
        case ADSORBATE:
          sprintf(buffer1,"Adsorbate molecule %d atom %d",HarmonicDistanceDefinitions[system][i][0][1],HarmonicDistanceDefinitions[system][i][0][2]);
          break;
        case CATION:
          sprintf(buffer1,"Cation molecule %d atom %d",HarmonicDistanceDefinitions[system][i][0][1],HarmonicDistanceDefinitions[system][i][0][2]);
          break;
        default:
          break;
      }
      switch(HarmonicDistanceDefinitions[system][i][1][0])
      {
        case FRAMEWORK:
          sprintf(buffer2,"Framework %d atom %d",HarmonicDistanceDefinitions[system][i][1][1],HarmonicDistanceDefinitions[system][i][1][2]);
          break;
        case ADSORBATE:
          sprintf(buffer2,"Adsorbate molecule %d atom %d",HarmonicDistanceDefinitions[system][i][1][1],HarmonicDistanceDefinitions[system][i][1][2]);
          break;
        case CATION:
          sprintf(buffer2,"Cation molecule %d atom %d",HarmonicDistanceDefinitions[system][i][1][1],HarmonicDistanceDefinitions[system][i][1][2]);
          break;
        default:
          break;
      }
      fprintf(FilePtr,"\tconstraint %d (%s,%s)  bond-constant: %g [K] constraint distance: %g [A]\n",i,buffer1,buffer2,
               HarmonicDistanceConstraintParameters[system][i][0]*ENERGY_TO_KELVIN,HarmonicDistanceConstraintParameters[system][i][1]);
    }
  }

  fprintf(FilePtr,"Harmonic angle constraints: %d\n",NumberOfHarmonicAngleConstraints[system]);
  if(NumberOfHarmonicAngleConstraints[system]>0)
  {
    for(i=0;i<NumberOfHarmonicAngleConstraints[system];i++)
    {
      switch(HarmonicAngleDefinitions[system][i][0][0])
      {
        case FRAMEWORK:
          sprintf(buffer1,"Framework %d atom %d",HarmonicAngleDefinitions[system][i][0][1],HarmonicAngleDefinitions[system][i][0][2]);
          break;
        case ADSORBATE:
          sprintf(buffer1,"Adsorbate molecule %d atom %d",HarmonicAngleDefinitions[system][i][0][1],HarmonicAngleDefinitions[system][i][0][2]);
          break;
        case CATION:
          sprintf(buffer1,"Cation molecule %d atom %d",HarmonicAngleDefinitions[system][i][0][1],HarmonicAngleDefinitions[system][i][0][2]);
          break;
        default:
          break;
      }
      switch(HarmonicAngleDefinitions[system][i][1][0])
      {
        case FRAMEWORK:
          sprintf(buffer2,"Framework %d atom %d",HarmonicAngleDefinitions[system][i][1][1],HarmonicAngleDefinitions[system][i][1][2]);
          break;
        case ADSORBATE:
          sprintf(buffer2,"Adsorbate molecule %d atom %d",HarmonicAngleDefinitions[system][i][1][1],HarmonicAngleDefinitions[system][i][1][2]);
          break;
        case CATION:
          sprintf(buffer2,"Cation molecule %d atom %d",HarmonicAngleDefinitions[system][i][1][1],HarmonicAngleDefinitions[system][i][1][2]);
          break;
        default:
          break;
      }
      switch(HarmonicAngleDefinitions[system][i][2][0])
      {
        case FRAMEWORK:
          sprintf(buffer3,"Framework %d atom %d",HarmonicAngleDefinitions[system][i][2][1],HarmonicAngleDefinitions[system][i][2][2]);
          break;
        case ADSORBATE:
          sprintf(buffer3,"Adsorbate molecule %d atom %d",HarmonicAngleDefinitions[system][i][2][1],HarmonicAngleDefinitions[system][i][2][2]);
          break;
        case CATION:
          sprintf(buffer3,"Cation molecule %d atom %d",HarmonicAngleDefinitions[system][i][2][1],HarmonicAngleDefinitions[system][i][2][2]);
          break;
        default:
          break;
      }

      fprintf(FilePtr,"\tconstraint %d (%s,%s,%s)  bend-constant: %g [K] constraint angle: %g [A]\n",i,buffer1,buffer2,buffer3,
               HarmonicAngleConstraintParameters[system][i][0]*ENERGY_TO_KELVIN,HarmonicAngleConstraintParameters[system][i][1]*RAD2DEG);
    }
  }

  fprintf(FilePtr,"Harmonic dihedral constraints: %d\n",NumberOfHarmonicDihedralConstraints[system]);
  if(NumberOfHarmonicDihedralConstraints[system]>0)
  {
    for(i=0;i<NumberOfHarmonicDihedralConstraints[system];i++)
    {
      switch(HarmonicDihedralDefinitions[system][i][0][0])
      {
        case FRAMEWORK:
          sprintf(buffer1,"Framework %d atom %d",HarmonicDihedralDefinitions[system][i][0][1],HarmonicDihedralDefinitions[system][i][0][2]);
          break;
        case ADSORBATE:
          sprintf(buffer1,"Adsorbate molecule %d atom %d",HarmonicDihedralDefinitions[system][i][0][1],HarmonicDihedralDefinitions[system][i][0][2]);
          break;
        case CATION:
          sprintf(buffer1,"Cation molecule %d atom %d",HarmonicDihedralDefinitions[system][i][0][1],HarmonicDihedralDefinitions[system][i][0][2]);
          break;
        default:
          break;
      }
      switch(HarmonicDihedralDefinitions[system][i][1][0])
      {
        case FRAMEWORK:
          sprintf(buffer2,"Framework %d atom %d",HarmonicDihedralDefinitions[system][i][1][1],HarmonicDihedralDefinitions[system][i][1][2]);
          break;
        case ADSORBATE:
          sprintf(buffer2,"Adsorbate molecule %d atom %d",HarmonicDihedralDefinitions[system][i][1][1],HarmonicDihedralDefinitions[system][i][1][2]);
          break;
        case CATION:
          sprintf(buffer2,"Cation molecule %d atom %d",HarmonicDihedralDefinitions[system][i][1][1],HarmonicDihedralDefinitions[system][i][1][2]);
          break;
        default:
          break;
      }
      switch(HarmonicDihedralDefinitions[system][i][2][0])
      {
        case FRAMEWORK:
          sprintf(buffer3,"Framework %d atom %d",HarmonicDihedralDefinitions[system][i][2][1],HarmonicDihedralDefinitions[system][i][2][2]);
          break;
        case ADSORBATE:
          sprintf(buffer3,"Adsorbate molecule %d atom %d",HarmonicDihedralDefinitions[system][i][2][1],HarmonicDihedralDefinitions[system][i][2][2]);
          break;
        case CATION:
          sprintf(buffer3,"Cation molecule %d atom %d",HarmonicDihedralDefinitions[system][i][2][1],HarmonicDihedralDefinitions[system][i][2][2]);
          break;
        default:
          break;
      }
      switch(HarmonicDihedralDefinitions[system][i][3][0])
      {
        case FRAMEWORK:
          sprintf(buffer4,"Framework %d atom %d",HarmonicDihedralDefinitions[system][i][3][1],HarmonicDihedralDefinitions[system][i][3][2]);
          break;
        case ADSORBATE:
          sprintf(buffer4,"Adsorbate molecule %d atom %d",HarmonicDihedralDefinitions[system][i][3][1],HarmonicDihedralDefinitions[system][i][3][2]);
          break;
        case CATION:
          sprintf(buffer4,"Cation molecule %d atom %d",HarmonicDihedralDefinitions[system][i][3][1],HarmonicDihedralDefinitions[system][i][3][2]);
          break;
        default:
          break;
      }
      fprintf(FilePtr,"\tconstraint %d (%s,%s,%s,%s)  dihedral-constant: %g [K] constraint angle: %g [A]\n",i,buffer1,buffer2,buffer3,buffer4,
               HarmonicDihedralConstraintParameters[system][i][0]*ENERGY_TO_KELVIN,HarmonicDihedralConstraintParameters[system][i][1]*RAD2DEG);
    }
  }
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Dihedral mid-point measurements: %d\n",NumberOfTwoPointDihedralDefinitions[system]);
  if(NumberOfTwoPointDihedralDefinitions[system]>0)
  {
    for(i=0;i<NumberOfTwoPointDihedralDefinitions[system];i++)
    {
      switch(TwoPointDihedralDefinitions[system][i][0][0])
      {
        case FRAMEWORK:
          sprintf(buffer1,"Framework %d atom %d",TwoPointDihedralDefinitions[system][i][0][1],TwoPointDihedralDefinitions[system][i][0][2]);
          break;
        case ADSORBATE:
          sprintf(buffer1,"Adsorbate molecule %d atom %d",TwoPointDihedralDefinitions[system][i][0][1],TwoPointDihedralDefinitions[system][i][0][2]);
          break;
        case CATION:
          sprintf(buffer1,"Cation molecule %d atom %d",TwoPointDihedralDefinitions[system][i][0][1],TwoPointDihedralDefinitions[system][i][0][2]);
          break;
        default:
          break;
      }
      switch(TwoPointDihedralDefinitions[system][i][1][0])
      {
        case FRAMEWORK:
          sprintf(buffer2,"Framework %d atom %d",TwoPointDihedralDefinitions[system][i][1][1],TwoPointDihedralDefinitions[system][i][1][2]);
          break;
        case ADSORBATE:
          sprintf(buffer2,"Adsorbate molecule %d atom %d",TwoPointDihedralDefinitions[system][i][1][1],TwoPointDihedralDefinitions[system][i][1][2]);
          break;
        case CATION:
          sprintf(buffer2,"Cation molecule %d atom %d",TwoPointDihedralDefinitions[system][i][1][1],TwoPointDihedralDefinitions[system][i][1][2]);
          break;
        default:
          break;
      }
      switch(TwoPointDihedralDefinitions[system][i][2][0])
      {
        case FRAMEWORK:
          sprintf(buffer3,"Framework %d atom %d",TwoPointDihedralDefinitions[system][i][2][1],TwoPointDihedralDefinitions[system][i][2][2]);
          break;
        case ADSORBATE:
          sprintf(buffer3,"Adsorbate molecule %d atom %d",TwoPointDihedralDefinitions[system][i][2][1],TwoPointDihedralDefinitions[system][i][2][2]);
          break;
        case CATION:
          sprintf(buffer3,"Cation molecule %d atom %d",TwoPointDihedralDefinitions[system][i][2][1],TwoPointDihedralDefinitions[system][i][2][2]);
          break;
        default:
          break;
      }
      switch(TwoPointDihedralDefinitions[system][i][3][0])
      {
        case FRAMEWORK:
          sprintf(buffer4,"Framework %d atom %d",TwoPointDihedralDefinitions[system][i][3][1],TwoPointDihedralDefinitions[system][i][3][2]);
          break;
        case ADSORBATE:
          sprintf(buffer4,"Adsorbate molecule %d atom %d",TwoPointDihedralDefinitions[system][i][3][1],TwoPointDihedralDefinitions[system][i][3][2]);
          break;
        case CATION:
          sprintf(buffer4,"Cation molecule %d atom %d",TwoPointDihedralDefinitions[system][i][3][1],TwoPointDihedralDefinitions[system][i][3][2]);
          break;
        default:
          break;
      }
      switch(TwoPointDihedralDefinitions[system][i][4][0])
      {
        case FRAMEWORK:
          sprintf(buffer5,"Framework %d atom %d",TwoPointDihedralDefinitions[system][i][4][1],TwoPointDihedralDefinitions[system][i][4][2]);
          break;
        case ADSORBATE:
          sprintf(buffer5,"Adsorbate molecule %d atom %d",TwoPointDihedralDefinitions[system][i][4][1],TwoPointDihedralDefinitions[system][i][4][2]);
          break;
        case CATION:
          sprintf(buffer5,"Cation molecule %d atom %d",TwoPointDihedralDefinitions[system][i][4][1],TwoPointDihedralDefinitions[system][i][4][2]);
          break;
        default:
          break;
      }
      switch(TwoPointDihedralDefinitions[system][i][5][0])
      {
        case FRAMEWORK:
          sprintf(buffer6,"Framework %d atom %d",TwoPointDihedralDefinitions[system][i][5][1],TwoPointDihedralDefinitions[system][i][5][2]);
          break;
        case ADSORBATE:
          sprintf(buffer6,"Adsorbate molecule %d atom %d",TwoPointDihedralDefinitions[system][i][5][1],TwoPointDihedralDefinitions[system][i][5][2]);
          break;
        case CATION:
          sprintf(buffer6,"Cation molecule %d atom %d",TwoPointDihedralDefinitions[system][i][5][1],TwoPointDihedralDefinitions[system][i][5][2]);
          break;
        default:
          break;
      }
      fprintf(FilePtr,"\tmeasure-dihedral  %d ([%s,%s],%s,%s,[%s,%s])\n",i,buffer1,buffer2,buffer3,buffer4,buffer5,buffer6);
    }
  }
  fprintf(FilePtr,"\n");

  if(NumberOfFixedAtomTypes>0)
  {
    fprintf(FilePtr,"Fixed atom types:");
    for(i=0;i<NumberOfFixedAtomTypes;i++)
      fprintf(FilePtr," %s",FixedAtomTypes[i]);
    fprintf(FilePtr,"\n\n");
  }

  if(NumberOfActiveAtomTypes>0)
  {
    fprintf(FilePtr,"Active atom types:");
    for(i=0;i<NumberOfActiveAtomTypes;i++)
      fprintf(FilePtr," %s",ActiveAtomTypes[i]);
    fprintf(FilePtr,"\n\n");
  }

  if(Framework[system].FrameworkModel==FLEXIBLE)
  {
    for(i=0;i<Framework[system].NumberOfFrameworks;i++)
    {
      nr_fixed=nr_free=0;
      for(k=0;k<Framework[system].NumberOfAtoms[i];k++)
      {
        if(Framework[system].Atoms[i][k].Fixed.x&&Framework[system].Atoms[i][k].Fixed.y&&Framework[system].Atoms[i][k].Fixed.z)
          nr_fixed++;
        else
          nr_free++;
      }

      fprintf(FilePtr,"Number of fixed framework atoms: %d  (active: %d)\n",Framework[system].NumberOfFixedAtoms[i],Framework[system].NumberOfFreeAtoms[i]);
      fprintf(FilePtr,"fixed atom-list: ");
      for(k=0;k<Framework[system].NumberOfAtoms[i];k++)
      {
        sprintf(buffer1,"(%d, [%c%c%c])",k,Framework[system].Atoms[i][k].Fixed.x?'X':'-',
            Framework[system].Atoms[i][k].Fixed.y?'Y':'-',Framework[system].Atoms[i][k].Fixed.z?'Z':'-');
        if(Framework[system].Atoms[i][k].Fixed.x||Framework[system].Atoms[i][k].Fixed.y||Framework[system].Atoms[i][k].Fixed.z)
          fprintf(FilePtr,"%s ",buffer1);
      }
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"Number of active framework atoms: %d  (fixed: %d)\n",nr_free,nr_fixed);
      fprintf(FilePtr,"active atom-list: ");
      for(k=0;k<Framework[system].NumberOfAtoms[i];k++)
        if(Framework[system].Atoms[i][k].Fixed.x&&Framework[system].Atoms[i][k].Fixed.y&&Framework[system].Atoms[i][k].Fixed.z)
          fprintf(FilePtr,"%d ",k);
      fprintf(FilePtr,"\n");
    }
  }
  else
    fprintf(FilePtr,"All framework atoms are fixed\n");
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Fixed adsorbate atoms:  ");
  for(m=0;m<NumberOfAdsorbateMolecules[system];m++)
  {
    Type=Adsorbates[system][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(!Components[Type].Groups[l].Rigid) // free unit
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];

          sprintf(buffer1,"(%d,%d,[%c%c%c])",m,A,Adsorbates[system][m].Atoms[k].Fixed.x?'X':'-',
              Adsorbates[system][m].Atoms[k].Fixed.y?'Y':'-',Adsorbates[system][m].Atoms[k].Fixed.z?'Z':'-');

          if(Adsorbates[system][m].Atoms[k].Fixed.x||Adsorbates[system][m].Atoms[k].Fixed.y||Adsorbates[system][m].Atoms[k].Fixed.z)
            fprintf(FilePtr,"%s ",buffer1);
        }
      }
    }
  }
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Fixed adsorbate groups (center-of-mass):  ");
  for(m=0;m<NumberOfAdsorbateMolecules[system];m++)
  {
    Type=Adsorbates[system][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        sprintf(buffer1,"(%d,%d,[%c%c%c])",m,l,Adsorbates[system][m].Groups[l].FixedCenterOfMass.x?'X':'-',
            Adsorbates[system][m].Groups[l].FixedCenterOfMass.y?'Y':'-',Adsorbates[system][m].Groups[l].FixedCenterOfMass.z?'Z':'-');

        if(Adsorbates[system][m].Groups[l].FixedCenterOfMass.x||Adsorbates[system][m].Groups[l].FixedCenterOfMass.y||Adsorbates[system][m].Groups[l].FixedCenterOfMass.z)
          fprintf(FilePtr,"%s ",buffer1);
      }
    }
  }
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Fixed adsorbate groups (orientation):  ");
  for(m=0;m<NumberOfAdsorbateMolecules[system];m++)
  {
    Type=Adsorbates[system][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        sprintf(buffer1,"(%d,%d,[%c%c%c])",m,l,Adsorbates[system][m].Groups[l].FixedOrientation.x?'X':'-',
            Adsorbates[system][m].Groups[l].FixedOrientation.y?'Y':'-',Adsorbates[system][m].Groups[l].FixedOrientation.z?'Z':'-');

        if(Adsorbates[system][m].Groups[l].FixedOrientation.x||Adsorbates[system][m].Groups[l].FixedOrientation.y||Adsorbates[system][m].Groups[l].FixedOrientation.z)
          fprintf(FilePtr,"%s ",buffer1);
      }
    }
  }
  fprintf(FilePtr,"\n\n");




  fprintf(FilePtr,"Fixed cation atoms:  ");
  for(m=0;m<NumberOfCationMolecules[system];m++)
  {
    Type=Cations[system][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(!Components[Type].Groups[l].Rigid) // rigid unit
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];

          sprintf(buffer1,"(%d,%d,[%c%c%c])",m,A,Cations[system][m].Atoms[k].Fixed.x?'X':'-',
              Cations[system][m].Atoms[k].Fixed.y?'Y':'-',Cations[system][m].Atoms[k].Fixed.z?'Z':'-');

          if(Cations[system][m].Atoms[k].Fixed.x||Cations[system][m].Atoms[k].Fixed.y||Cations[system][m].Atoms[k].Fixed.z)
            fprintf(FilePtr,"%s ",buffer1);
        }
      }
    }
  }
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Fixed cation groups (center-of-mass):  ");
  for(m=0;m<NumberOfCationMolecules[system];m++)
  {
    Type=Cations[system][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        sprintf(buffer1,"(%d,%d,[%c%c%c])",m,l,Cations[system][m].Groups[l].FixedCenterOfMass.x?'X':'-',
            Cations[system][m].Groups[l].FixedCenterOfMass.y?'Y':'-',Cations[system][m].Groups[l].FixedCenterOfMass.z?'Z':'-');

        if(Cations[system][m].Groups[l].FixedCenterOfMass.x||Cations[system][m].Groups[l].FixedCenterOfMass.y||Cations[system][m].Groups[l].FixedCenterOfMass.z)
          fprintf(FilePtr,"%s ",buffer1);
      }
    }
  }
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Fixed cation groups (orientation):  ");
  for(m=0;m<NumberOfCationMolecules[system];m++)
  {
    Type=Cations[system][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        sprintf(buffer1,"(%d,%d,[%c%c%c])",m,l,Cations[system][m].Groups[l].FixedOrientation.x?'X':'-',
            Cations[system][m].Groups[l].FixedOrientation.y?'Y':'-',Cations[system][m].Groups[l].FixedOrientation.z?'Z':'-');

        if(Cations[system][m].Groups[l].FixedOrientation.x||Cations[system][m].Groups[l].FixedOrientation.y||Cations[system][m].Groups[l].FixedOrientation.z)
          fprintf(FilePtr,"%s ",buffer1);
      }
    }
  }
  fprintf(FilePtr,"\n\n\n");


  fprintf(FilePtr,"dcTST parameters\n");
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"Free energy profiles computed: %s\n",ComputeFreeEnergyProfile[system]?"yes":"no");
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].ComputeFreeEnergyProfile[system])
        fprintf(FilePtr,"\t%s-profile\n",Components[i].Name);
  }
  fprintf(FilePtr,"Free energy profiles written every %d cycles\n",WriteFreeEnergyProfileEvery[system]);
  switch(FreeEnergyMappingType[system])
  {
    case A_MAPPING:
      fprintf(FilePtr,"Free energy mapping: mapped to a-coordinate\n");
      break;
    case B_MAPPING:
      fprintf(FilePtr,"Free energy mapping: mapped to b-coordinate\n");
      break;
    case C_MAPPING:
      fprintf(FilePtr,"Free energy mapping: mapped to c-coordinate\n");
      break;
    case ABC_MAPPING:
      fprintf(FilePtr,"Free energy mapping: mapped to a,b,c-coordinates\n");
      break;
    case MAP_AB_DIAGONAL:
      fprintf(FilePtr,"Free energy mapping: mapped to the diagonal formed by a,b\n");
      break;
    case MAP_AC_DIAGONAL:
      fprintf(FilePtr,"Free energy mapping: mapped to the diagonal formed by a,c\n");
      break;
    case MAP_BC_DIAGONAL:
      fprintf(FilePtr,"Free energy mapping: mapped to the diagonal formed by b,c\n");
      break;
    case MAP_O_AB_DIAGONAL:
      fprintf(FilePtr,"Free energy mapping: mapped to the diagonal formed by O and a+b\n");
      break;
    case MAP_O_AC_DIAGONAL:
      fprintf(FilePtr,"Free energy mapping: mapped to the diagonal formed by O and a+c\n");
      break;
    case MAP_O_BC_DIAGONAL:
      fprintf(FilePtr,"Free energy mapping: mapped to the diagonal formed by O and b+c\n");
      break;
    case MAP_A_BC_DIAGONAL:
      fprintf(FilePtr,"Free energy mapping: mapped to the diagonal formed by a and b+c\n");
      break;
    case MAP_B_AC_DIAGONAL:
      fprintf(FilePtr,"Free energy mapping: mapped to the diagonal formed by b and a+c\n");
      break;
    case MAP_C_AB_DIAGONAL:
      fprintf(FilePtr,"Free energy mapping: mapped to the diagonal formed by c and a+b\n");
      break;
    case MAP_O_ABC_DIAGONAL:
      fprintf(FilePtr,"Free energy mapping: mapped to the diagonal formed by O and a+b+c\n");
      break;
  }

  fprintf(FilePtr,"BarrierPosition: %18.10lf %18.10lf %18.10lf\n",BarrierPosition[system].x,BarrierPosition[system].y,BarrierPosition[system].z);
  fprintf(FilePtr,"BarrierNormal:   %18.10lf %18.10lf %18.10lf\n",BarrierNormal[system].x,BarrierNormal[system].y,BarrierNormal[system].z);
  fprintf(FilePtr,"Start with a molecule on top of the barrier: %s\n",PutMoleculeOnBarrier[system]?"yes":"no");
  fprintf(FilePtr,"Maximum distance to barrier (e.g. distance to minumum free energy): %18.10lf [A]\n",MaxBarrierDistance[system]);
  fprintf(FilePtr,"Maximum trajectory time: %18.10lf [ps]\n",MaxBarrierTime[system]);
  fprintf(FilePtr,"Each configuration is used with %d different initial velocities\n",NumberOfVelocities[system]);
  fprintf(FilePtr,"\n\n");

  fprintf(FilePtr,"Cbmc parameters\n");
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"Biasing method: %s\n",BiasingMethod==LJ_BIASING?"using only the VDW part":
                   "using the VDW and the real part of the Ewald summation");
  fprintf(FilePtr,"Number of trial positions:                                       %d\n",NumberOfTrialPositions);
  fprintf(FilePtr,"Number of trial positions (reinsertion):                         %d\n",NumberOfTrialPositionsReinsertion);
  fprintf(FilePtr,"Number of trial positions (partial reinsertion):                 %d\n",NumberOfTrialPositionsPartialReinsertion);
  fprintf(FilePtr,"Number of trial positions (identity-change):                     %d\n",NumberOfTrialPositionsIdentityChange);
  fprintf(FilePtr,"Number of trial positions (Gibbs particle transfer):             %d\n",NumberOfTrialPositionsGibbs);
  fprintf(FilePtr,"Number of trial positions (insertion/deletion):                  %d\n",NumberOfTrialPositionsSwap);
  fprintf(FilePtr,"Number of trial positions (Widom insertion):                     %d\n",NumberOfTrialPositionsWidom);
  fprintf(FilePtr,"Number of trial positions coupled Torsion-selection:             %d\n",NumberOfTrialPositionsTorsion);
  fprintf(FilePtr,"Number of trial positions first bead:                            %d\n",NumberOfTrialPositionsForTheFirstBead);
  fprintf(FilePtr,"Number of trial positions first bead (reinsertion):              %d\n",NumberOfTrialPositionsForTheFirstBeadReinsertion);
  fprintf(FilePtr,"Number of trial positions first bead (partial reinsertion):      %d\n",NumberOfTrialPositionsForTheFirstBeadPartialReinsertion);
  fprintf(FilePtr,"Number of trial positions first bead (identity-change):          %d\n",NumberOfTrialPositionsForTheFirstBeadIdentityChange);
  fprintf(FilePtr,"Number of trial positions first bead (Gibbs particle transfer):  %d\n",NumberOfTrialPositionsForTheFirstBeadGibbs);
  fprintf(FilePtr,"Number of trial positions first bead (insertion/deletion):       %d\n",NumberOfTrialPositionsForTheFirstBeadSwap);
  fprintf(FilePtr,"Number of trial positions first bead (Widom insertion):          %d\n",NumberOfTrialPositionsForTheFirstBeadWidom);
  fprintf(FilePtr,"Number of trial moves per open bead:                             %d\n",NumberOfTrialMovesPerOpenBead);
  fprintf(FilePtr,"Target acceptance ratio small-mc scheme:                         %lf\n",(double)TargetAccRatioSmallMCScheme);
  fprintf(FilePtr,"Energy overlap criteria:                                         %lg\n",(double)EnergyOverlapCriteria);
  fprintf(FilePtr,"Minimal Rosenbluth factor:                                       %lg\n",(double)MinimumRosenbluthFactor);
  fprintf(FilePtr,"\n\n");

  if(WritePseudoAtomsToOutput)
  {
    fprintf(FilePtr,"Pseudo atoms: %d\n",NumberOfPseudoAtoms);
    fprintf(FilePtr,"===========================================================================\n");
      for(i=0;i<NumberOfPseudoAtoms;i++)
      {
        sprintf(framework_charge_string,"no");
        switch(PseudoAtoms[i].ChargeDefinitionType)
        {
          case CHARGE_ATOM_FROM_PSEUDO_ATOM_DEFINITION:
            sprintf(charge_string,"%-12.9f       ",PseudoAtoms[i].Charge1);
            if(PseudoAtoms[i].FrameworkAtom)
              sprintf(framework_charge_string,"yes (charge from pseudo-atoms file)");
            break;
          case CHARGE_ATOM_FROM_STRUCTURE_FILE:
            sprintf(charge_string,"%-12.9f   (av)",PseudoAtoms[i].Charge1);
            if(PseudoAtoms[i].FrameworkAtom)
              sprintf(framework_charge_string,"yes (charge from structure file)");
            break;
          case CHARGE_NOT_DEFINED:
            sprintf(charge_string,"%-12.9f (-na-)",PseudoAtoms[i].Charge1);
            if(PseudoAtoms[i].FrameworkAtom)
              sprintf(framework_charge_string,"yes (charge definition not found)");
            break;
          default:
            fprintf(stderr, "ERROR: undefined ChargeDefinitionType (output.c)\n");
            exit(0);
            break;
        }

        switch(PolarizationMatrix)
        {
          case ISOTROPIC:
            sprintf(polarization_string,"%-12.9f",PseudoAtoms[i].Polarization.ax*COULOMBIC_CONVERSION_FACTOR);
            break;
          case ANISOTROPIC:
            sprintf(polarization_string,"%-12.9f,%-12.9f,%-12.9f",PseudoAtoms[i].Polarization.ax*COULOMBIC_CONVERSION_FACTOR,
                PseudoAtoms[i].Polarization.by*COULOMBIC_CONVERSION_FACTOR,PseudoAtoms[i].Polarization.cz*COULOMBIC_CONVERSION_FACTOR);
            break;
          case REGULAR_UPPER_TRIANGLE:
            sprintf(polarization_string,"%-12.9f,%-12.9f,%-12.9f,%-12.9f,%-12.9f,%-12.9f",
                PseudoAtoms[i].Polarization.ax*COULOMBIC_CONVERSION_FACTOR,PseudoAtoms[i].Polarization.ay*COULOMBIC_CONVERSION_FACTOR,
                PseudoAtoms[i].Polarization.az*COULOMBIC_CONVERSION_FACTOR,PseudoAtoms[i].Polarization.by*COULOMBIC_CONVERSION_FACTOR,
                PseudoAtoms[i].Polarization.bz*COULOMBIC_CONVERSION_FACTOR,PseudoAtoms[i].Polarization.cz*COULOMBIC_CONVERSION_FACTOR);
            break;
        }

        fprintf(FilePtr,"Pseudo Atom[%4d] Name %-8s Oxidation: %-8s Element: %-4s pdb-name: %-4s Scat. Types: %3d %3d Mass=%-12.9f B-factor:%-8.3lf\n"
                      "                 Charge=%s  Polarization=%s [A^3] (considered %-s and %-s)  Interactions: %s\n"
                      "                 Anisotropic factor: %8.3lf [-] (%s), Radius: %8.3lf [A], Framework-atom: %3s\n",
           i,
           PseudoAtoms[i].Name,
           PseudoAtoms[i].OxidationStateString,
           PseudoAtoms[i].ChemicalElement,
           PseudoAtoms[i].PrintToPDBName,
           PseudoAtoms[i].ScatteringType,
           PseudoAtoms[i].AnomalousScatteringType,
           (double)PseudoAtoms[i].Mass,
           (double)PseudoAtoms[i].TemperatureFactor,
           charge_string,
           polarization_string,
           PseudoAtoms[i].HasCharges?"a charged atom":"a chargeless atom",
           PseudoAtoms[i].IsPolarizable?"polarizable":"no polarization",
           PseudoAtoms[i].Interaction?"yes":" no",
           (double)PseudoAtoms[i].AnisotropicDisplacement,
           PseudoAtoms[i].AnisotropicType==RELATIVE?"Relative":"Absolute",
           PseudoAtoms[i].Radius,
           framework_charge_string);
      }
    fprintf(FilePtr,"\n\n");
  }

  if(WriteForcefieldToOutput)
  {
    fprintf(FilePtr,"Forcefield: %s\n",ForceField);
    fprintf(FilePtr,"===========================================================================\n");
    fprintf(FilePtr,"Minimal distance: %lg\n",(double)sqrt(OverlapDistanceSquared));
    if(BoundaryCondition[system]==FINITE)
      fprintf(FilePtr,"No VDW cutOff\n");
    else
    {
      fprintf(FilePtr,"CutOff VDW : %lf (%lf)\n",(double)CutOffVDW,(double)CutOffVDWSquared);
      fprintf(FilePtr,"CutOff VDW switching on: %lf (%lf)\n",(double)CutOffVDWSwitch,(double)CutOffVDWSwitchSquared);
      fprintf(FilePtr,"CutOff charge-charge : %lf (%lf)\n",(double)CutOffChargeCharge[system],(double)CutOffChargeChargeSquared[system]);
      fprintf(FilePtr,"CutOff charge-charge switching on: %lf (%lf)\n",(double)CutOffChargeChargeSwitch[system],(double)CutOffChargeChargeSwitchSquared[system]);
      fprintf(FilePtr,"CutOff charge-bonddipole : %lf (%lf)\n",(double)CutOffChargeBondDipole,(double)CutOffChargeBondDipoleSquared);
      fprintf(FilePtr,"CutOff charge-bondipole switching on: %lf (%lf)\n",(double)CutOffChargeBondDipoleSwitch,(double)CutOffChargeBondDipoleSwitchSquared);
      fprintf(FilePtr,"CutOff bonddipole-bonddipole : %lf (%lf)\n",(double)CutOffBondDipoleBondDipole,(double)CutOffBondDipoleBondDipoleSquared);
      fprintf(FilePtr,"CutOff bonddipole-bondipole switching on: %lf (%lf)\n",(double)CutOffBondDipoleBondDipoleSwitch,(double)CutOffBondDipoleBondDipoleSwitchSquared);
    }
    if(ComputePolarization)
    {
      fprintf(FilePtr,"Polarization is computed\n");
      fprintf(FilePtr,"Framework-Framework polarization is %s\n",OmitIntraFrameworkPolarization?"ommited":"included");
      fprintf(FilePtr,"Adsorbate-Adsorbate polarization is %s\n",OmitAdsorbateAdsorbatePolarization?"ommited":"included");
      fprintf(FilePtr,"Cation-Cation polarization is %s\n",OmitCationCationPolarization?"ommited":"included");
      fprintf(FilePtr,"Adsorbate-Cation polarization is %s\n",OmitAdsorbateCationPolarization?"ommited":"included");
      if(BackPolarization)
      {
        fprintf(FilePtr,"Back polarization is computed by self-consistent iteration\n");
        fprintf(FilePtr,"%d steps are used to iterate to convergence\n",NumberOfBackPolarizationSteps);
      }
      else
        fprintf(FilePtr,"Back polarization is neglected\n");
    }
    else
      fprintf(FilePtr,"Polarization is neglected\n");
    if(OmitAdsorbateAdsorbateVDWInteractions)
      fprintf(FilePtr,"adsorbate-adsorbate VDW interactions are omitted\n");
    if(OmitAdsorbateAdsorbateCoulombInteractions)
      fprintf(FilePtr,"adsorbate-adsorbate Coulombic interactions are omitted\n");
    if(OmitAdsorbateCationVDWInteractions)
      fprintf(FilePtr,"adsorbate-cation VDW interactions are omitted\n");
    if(OmitAdsorbateCationCoulombInteractions)
      fprintf(FilePtr,"adsorbate-cation Coulombic interactions are omitted\n");
    if(OmitCationCationVDWInteractions)
      fprintf(FilePtr,"cation-cation VDW interactions are omitted\n");
    if(OmitCationCationCoulombInteractions)
      fprintf(FilePtr,"cation-cation Coulombic interactions are omitted\n");
    if(TailCorrections)
      fprintf(FilePtr,"TailCorrections are used\n");
    if(ShiftPotentials)
      fprintf(FilePtr,"All potentials are shifted to zero at the Cutoff\n");
    else
      fprintf(FilePtr,"All potentials are unshifted !!!!!!\n");
    fprintf(FilePtr,"\n");

    switch(GeneralMixingRule)
    {
      case NONE:
        fprintf(FilePtr,"No general mixing rule: all cross term specified individualy\n");
        break;
      case LORENTZ_BERTHELOT:
        fprintf(FilePtr,"General mixing rule: Lorentz-Berthelot mixing rules are used FIRST for cross terms\n");
        break;
      case JORGENSEN:
        fprintf(FilePtr,"General mixing rule: Jorgensen mixing rules are used FIRST for cross terms\n");
        break;
    }
    fprintf(FilePtr,"%d cross terms are overwritten using the individual mixing rules from the file 'force_field_mixing_rules.def'\n",
           IndividualMixingRules);
    fprintf(FilePtr,"and then %d terms are overwritten using the specific interactions from the file 'force_field.def'\n",
           IndividualInteractions);
    fprintf(FilePtr,"\n");

    fprintf(FilePtr,"The force field and all the interactions:\n");
    for(i=1;i<NumberOfPseudoAtoms;i++)
      for(j=i;j<NumberOfPseudoAtoms;j++)
      {
        switch(PotentialType[i][j])
        {
          case ZERO_POTENTIAL:
            fprintf(FilePtr,"%7s - %7s [ZERO_POTENTIAL]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name);
            break;
          case ZERO_POTENTIAL_CONTINUOUS_FRACTIONAL:
            fprintf(FilePtr,"%7s - %7s [ZERO_POTENTIAL_CONTINUOUS_FRACTIONAL]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name);
            break;
          case HARD_SPHERE:
            fprintf(FilePtr,"%7s - %7s [HARD_SPHERE] radius: %9.5lf\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]);
            break;
          case LENNARD_JONES:
            // 4*p_0*((p_1/r)^12-(p_1/r)^6)
            // ======================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            // p_2/k_B [K]    (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [LENNARD_JONES] p_0/k_B: %9.5lf [K], p_1: %7.5lf [A], shift/k_B: %12.8lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case LENNARD_JONES_SMOOTHED3:
            // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
            // ======================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            fprintf(FilePtr,"%7s - %7s [LENNARD_JONES_SMOOTHED3] p_0/k_B: %9.5lf [K], p_1: %7.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]);
            break;
          case LENNARD_JONES_SMOOTHED5:
            // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
            // ======================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            fprintf(FilePtr,"%7s - %7s [LENNARD_JONES_SMOOTHED5] p_0/k_B: %9.5lf, p_1: %7.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]);
            break;
          case LENNARD_JONES_CONTINUOUS_FRACTIONAL:
            // 4*p_0*((p_1/r)^12-(p_1/r)^6)
            // ======================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            // p_2/k_B [K]    (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [LENNARD_JONES_CONTINUOUS_FRACTIONAL] p_0/k_B: %9.5lf [K], p_1: %7.5lf [A], shift/k_B: %12.8lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case LENNARD_JONES_CONTINUOUS_FRACTIONAL_SMOOTHED3:
            // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
            // ======================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            fprintf(FilePtr,"%7s - %7s [LENNARD_JONES_CONTINUOUS_FRACTIONAL_SMOOTHED3] p_0/k_B: %9.5lf [K], p_1: %7.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]);
            break;
          case LENNARD_JONES_CONTINUOUS_FRACTIONAL_SMOOTHED5:
            // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
            // ======================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            fprintf(FilePtr,"%7s - %7s [LENNARD_JONES_CONTINUOUS_FRACTIONAL_SMOOTHED5] p_0/k_B: %9.5lf, p_1: %7.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]);
            break;
          case WCA:
            // 4*p_0*((p_1/r)^12-(p_1/r)^6)
            // ======================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            // p_2/k_B [K]    (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [WCA] p_0/k_B: %9.5lf [K], p_1: %7.5lf [A], shift/k_B: %12.8lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case FEYNMAN_HIBBS_LENNARD_JONES:
            // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
            // =============================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            // p_2     [u]    reduced mass in unified atomic mass units
            // p_3/k_B [K]    (non-zero for a shifted potential)
            // T       [K]    the temperature
            fprintf(FilePtr,"%7s - %7s [FEYNMAN_HIBBS_LENNARD_JONES] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A], p_2: %8.5lf [u] (h^2/(24p_2k_BT)=%8.5f), shift/k_B: %8.5lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2],
              (double)(FH_CONVERSION_FACTOR/PotentialParms[i][j][2]),
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3:
            // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
            // ====================================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            // p_2     [u]    reduced mass in unified atomic mass units
            // T       [K]    the temperature
            fprintf(FilePtr,"%7s - %7s [FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A], p_2: %8.5lf [u]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]);
            break;
          case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5:
            // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
            // ====================================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            // p_2     [u]    reduced mass in unified atomic mass units
            // T       [K]    the temperature
            fprintf(FilePtr,"%7s - %7s [FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A], p_2: %8.5lf [u]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]);
            break;
          case FEYNMAN_HIBBS_LENNARD_JONES2:
            // 4*p_0*((p_1/r)^12-(p_1/r)^6)+4 [p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2]*(p_2/T)
            // ======================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            // p_2     [A^2]  correction factor
            // p_3/k_B [K]    (non-zero for a shifted potential)
            // T       [K]    the temperature
            fprintf(FilePtr,"%7s - %7s [FEYNMAN_HIBBS_LENNARD_JONES2] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A], p_2: %8.5lf [A^2], shift/k_B: %8.5lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2],
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3:
            // {4*p_0*((p_1/r)^12-(p_1/r)^6)+4 [p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2]*(p_2/T)}*S(r)
            // ======================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            // p_2     [A^2]  correction factor
            // T       [K]    the temperature
            fprintf(FilePtr,"%7s - %7s [FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A], p_2: %8.5lf [A^2]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]);
            break;
          case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5:
            // {4*p_0*((p_1/r)^12-(p_1/r)^6)+4 [p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2]*(p_2/T)}*S(r)
            // ======================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            // p_2     [A^2]  correction factor
            // T       [K]    the temperature
            fprintf(FilePtr,"%7s - %7s [FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A], p_2: %8.5lf [A^2]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]);
            break;
          case LENNARD_JONES_SHIFTED_FORCE:
            // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]-[(p_1/rc)^12-(p_1/rc)^6]}+[12*(p_1/rc)^12-6*(p_1/rc)^6]*(r-rc)/rc
            // ===============================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            fprintf(FilePtr,"%7s - %7s [LENNARD_JONES_SHIFTED_FORCE] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]);
            break;
          case LENNARD_JONES_SHIFTED_FORCE2:
            // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]+[6*(p_1/rc)^12-3*(p_1/rc)^6]}*r^2/rc^2+[7*(p_1/rc)^12+4*(p_1/rc)^6]
            // =================================================================================================
            // p_0/k_B [K]    strength parameter epsilon
            // p_1     [A]    size parameter sigma
            fprintf(FilePtr,"%7s - %7s [LENNARD_JONES_SHIFTED_FORCE2] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]);
            break;
          case POTENTIAL_12_6:
            // p_0/r^12-p_1/r^6
            // ======================================================================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^6]
            // p_2/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [POTENTIAL_12_6] p_0/k_B: %8.5lf [K A^12], p_1/k_B: %8.5lf [K A^6], "
                            "shift/k_B: %8.5lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case POTENTIAL_12_6_SMOOTHED3:
            // {p_0/r^12-p_1/r^6}*S(r)
            // ======================================================================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^6]
            fprintf(FilePtr,"%7s - %7s [POTENTIAL_12_6_SMOOTHED3] p_0/k_B: %8.5lf [K A^12], p_1/k_B: %8.5lf [K A^6]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN);
            break;
          case POTENTIAL_12_6_SMOOTHED5:
            // {p_0/r^12-p_1/r^6}*S(r)
            // ======================================================================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^6]
            fprintf(FilePtr,"%7s - %7s [POTENTIAL_12_6_SMOOTHED5] p_0/k_B: %lf [K A^12], p_1/k_B: %8.5lf [K A^6]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN);
            break;
          case POTENTIAL_12_6_2_0:
            // p_0/r^12+p_1/r^6+p_2/r^2+p_3
            // ======================================================================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^6]
            // p_2/k_B [K A^2]
            // p_3/k_B [K]
            // p_4/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [POTENTIAL_12_6_2_0] p_0/k_B: %8.5lf [K A^12], p_1/k_B: %8.5lf [K A^6], "
                            "p_2/k_B: %8.5lf [K A^2], p_3/k_B: %8.5lf [K], shift/k_B: %8.5lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case POTENTIAL_12_6_2_0_SMOOTHED3:
            // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
            // ======================================================================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^6]
            // p_2/k_B [K A^2]
            // p_3/k_B [K]
            fprintf(FilePtr,"%7s - %7s [POTENTIAL_12_6_2_0_SMOOTHED3] p_0/k_B: %8.5lf [K A^12], p_1/k_B: %8.5lf [K A^6], "
                            "p_2/k_B: %8.5lf [K A^2], p_3/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN);
            break;
          case POTENTIAL_12_6_2_0_SMOOTHED5:
            // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
            // ======================================================================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^6]
            // p_2/k_B [K A^2]
            // p_3/k_B [K]
            fprintf(FilePtr,"%7s - %7s [POTENTIAL_12_6_2_0_SMOOTHED5] p_0/k_B: %8.5lf [K A^12], p_1/k_B: %8.5lf [K A^6], "
                            "p_2/k_B: %8.5lf [K A^2], p_3/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN);
            break;
          case MORSE:
            // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
            // =================================================================================
            // p_0/k_B [K]       force constant
            // p_1     [A^-1]    parameter
            // p_2     [A]       reference distance
            // p_3/k_B [K]       (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [MORSE] p_0/k_B: %9.5lf [K], p_1: %8.5lf [A^-1], p_2: %8.5lf [A], "
                            "shift/k_B: %12.8lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              PotentialParms[i][j][1],
              PotentialParms[i][j][2],
              PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case MORSE_SMOOTHED3:
            // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]*S(r)
            // =================================================================================
            // p_0/k_B [K]       force constant
            // p_1     [A^-1]    parameter
            // p_2     [A]       reference distance
            fprintf(FilePtr,"%7s - %7s [MORSE_SMOOTHED3] p_0/k_B: %9.5lf [K], p_1: %8.5lf [A^-1], p_2: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              PotentialParms[i][j][1],
              PotentialParms[i][j][2]);
            break;
          case MORSE_SMOOTHED5:
            // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]*S(r)
            // =================================================================================
            // p_0/k_B [K]       force constant
            // p_1     [A^-1]    parameter
            // p_2     [A]       reference distance
            fprintf(FilePtr,"%7s - %7s [MORSE_SMOOTHED5] p_0/k_B: %9.5lf [K], p_1: %8.5lf [A^-1], p_2: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              PotentialParms[i][j][1],
              PotentialParms[i][j][2]);
            break;
          case MORSE2:
            // p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]
            // =================================================================================
            // p_0/k_B [K]       force constant
            // p_1     [A^-1]    parameter
            // p_2     [A]       reference distance
            // p_3/k_B [K]       (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [MORSE2] p_0/k_B: %9.5lf [K], p_1: %8.5lf [A^-1], p_2: %8.5lf [A], "
                            "shift/k_B: %12.8lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              PotentialParms[i][j][1],
              PotentialParms[i][j][2],
              PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case MORSE2_SMOOTHED3:
            // {p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]}*S(r)
            // =================================================================================
            // p_0/k_B [K]       force constant
            // p_1     [A^-1]    parameter
            // p_2     [A]       reference distance
            fprintf(FilePtr,"%7s - %7s [MORSE2_SMOOTHED3] p_0/k_B: %9.5lf [K], p_1: %8.5lf [A^-1], p_2: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              PotentialParms[i][j][1],
              PotentialParms[i][j][2]);
            break;
          case MORSE2_SMOOTHED5:
            // {p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]}*S(r)
            // =================================================================================
            // p_0/k_B [K]       force constant
            // p_1     [A^-1]    parameter
            // p_2     [A]       reference distance
            fprintf(FilePtr,"%7s - %7s [MORSE2_SMOOTHED5] p_0/k_B: %9.5lf [K], p_1: %8.5lf [A^-1], p_2: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              PotentialParms[i][j][1],
              PotentialParms[i][j][2]);
            break;
          case MORSE3:
            // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))*(r/p_2-2^(1/6))])^2-1}
            // =================================================================================
            // p_0/k_B [K]       force constant
            // p_1     [A]       reference distance
            // p_2/k_B [K]       (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [MORSE3] p_0/k_B: %9.5lf [K], p_2: %8.5lf [A], "
                            "shift/k_B: %12.8lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              PotentialParms[i][j][1],
              PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case MORSE3_SMOOTHED3:
            // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))*(r/p_2-2^(1/6))])^2-1}*S(r)
            // =================================================================================
            // p_0/k_B [K]       force constant
            // p_1     [A]       reference distance
            fprintf(FilePtr,"%7s - %7s [MORSE3_SMOOTHED3] p_0/k_B: %9.5lf [K], p_2: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              PotentialParms[i][j][1]);
            break;
          case MORSE3_SMOOTHED5:
            // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))*(r/p_2-2^(1/6))])^2-1}*S(r)
            // =================================================================================
            // p_0/k_B [K]       force constant
            // p_1     [A]       reference distance
            fprintf(FilePtr,"%7s - %7s [MORSE3_SMOOTHED5] p_0/k_B: %9.5lf [K], p_2: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              PotentialParms[i][j][1]);
            break;
          case CFF_9_6:
            // p_0/r^9-p_1/r^6
            // ======================================================================================
            // p_0/k_B [K A^9]
            // p_1/k_B [K A^6]
            // p_2/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [CFF_9_6] p_0/k_B: %8.5lf [K A^9], p_1/k_B: %8.5lf [K A^6], shift/k_B: %8.5lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case CFF_9_6_SMOOTHED3:
            // {p_0/r^9-p_1/r^6}*S(r)
            // ======================================================================================
            // p_0/k_B [K A^9]
            // p_1/k_B [K A^6]
            fprintf(FilePtr,"%7s - %7s [CFF_9_6_SMOOTHED3] p_0/k_B: %8.5lf [K A^9], p_1/k_B: %8.5lf [K A^6]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN);
            break;
          case CFF_9_6_SMOOTHED5:
            // {p_0/r^9-p_1/r^6}*S(r)
            // ======================================================================================
            // p_0/k_B [K A^9]
            // p_1/k_B [K A^6]
            fprintf(FilePtr,"%7s - %7s [CFF_9_6_SMOOTHED5] p_0/k_B: %8.5lf [K A^9], p_1/k_B: %8.5lf [K A^6]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN);
            break;
          case CFF_EPS_SIGMA:
            // p_0*[2*(p_1/r)^9-3(p_1/r)^6]
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A]
            // p_2/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [CFF_EPS_SIGMA] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A], shift/k_B: %8.5lf [K], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              TailCorrection[i][j]?"yes":"no");
            break;
          case CFF_EPS_SIGMA_SMOOTHED3:
            // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A]
            fprintf(FilePtr,"%7s - %7s [CFF_EPS_SIGMA_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]);
            break;
          case CFF_EPS_SIGMA_SMOOTHED5:
            // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A]
            fprintf(FilePtr,"%7s - %7s [CFF_EPS_SIGMA_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]);
            break;
          case BUCKINGHAM:
            // p_0*exp(-p_1*r)-p_2/r^6
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^6]
            // p_3/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [BUCKINGHAM] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^6], shift/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN);
            break;
          case BUCKINGHAM_SMOOTHED3:
            // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^6]
            fprintf(FilePtr,"%7s - %7s [BUCKINGHAM_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^6]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN);
            break;
          case BUCKINGHAM_SMOOTHED5:
            // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^6]
            fprintf(FilePtr,"%7s - %7s [BUCKINGHAM_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^6]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN);
            break;
          case BUCKINGHAM2:
            // if(r<p_3) 1e10 else p_0*exp(-p_1*r)-p_2/r^6
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^6]
            // p_3/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [BUCKINGHAM2] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^6], p_3: %8.5lf [A], shift/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3],
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN);
            break;
          case BUCKINGHAM2_SMOOTHED3:
            // if(r<p_3) 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^6]
            fprintf(FilePtr,"%7s - %7s [BUCKINGHAM2_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^6], p_3: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]);
            break;
          case BUCKINGHAM2_SMOOTHED5:
            // if(r<p_3) 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^6]
            fprintf(FilePtr,"%7s - %7s [BUCKINGHAM2_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^6], p_3: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]);
            break;
          case DZUBAK2012:
            // if(r<p_3) 1e10 else p_0*exp(-p_1*r)-p_2/r^6
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^5]
            // p_3/k_B [K A^6]
            // p_4     [A]  (non-zero for a shifted potential)
            // p_5/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [DZUBAK2012] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^5], p_3/k_B: %8.5lf [K A^6], p_4: %8.5lf [A], shift/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4],
              (double)PotentialParms[i][j][5]*ENERGY_TO_KELVIN);
            break;
          case MM3_VDW:
            // sqrt(p_0^i*p_0^j)*[1.84e-5*exp(-12/P)-2.25*P^6]  if P>=3.02
            // sqrt(p_0^i*p_0^j)*192.27*P^2                     if P<3.02
            // ======================================================================================
            // p_0     [kcal/mol]
            // p_1     [A]
            // p_3     [kcal/mol]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [MM3_VDW] p_0: %8.5lf [kcal/mol], p_1: %8.5lf [A], shift: %8.5lf [kcal/mol], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KCAL_PER_MOL,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KCAL_PER_MOL,
              TailCorrection[i][j]?"yes":"no");
            break;
          case MM3_VDW_SMOOTHED3:
            // {sqrt(p_0^i*p_0^j)*[1.84e-5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
            // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                     if P<3.02
            // ======================================================================================
            // p_0     [kcal/mol]
            // p_1     [A]
            fprintf(FilePtr,"%7s - %7s [MM3_VDW_SMOOTHED3] p_0: %8.5lf [kcal/mol], p_1: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KCAL_PER_MOL,
              (double)PotentialParms[i][j][1]);
            break;
          case MM3_VDW_SMOOTHED5:
            // {sqrt(p_0^i*p_0^j)*[1.84e-5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
            // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                     if P<3.02
            // ======================================================================================
            // p_0     [kcal/mol]
            // p_1     [A]
            fprintf(FilePtr,"%7s - %7s [MM3_VDW_SMOOTHED5] p_0: %8.5lf [kcal/mol], p_1: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KCAL_PER_MOL,
              (double)PotentialParms[i][j][1]);
            break;

          // obsolete
          case MM3_HYDROGEN_VDW:
            fprintf(FilePtr,"%7s - %7s [MM3_HYDROGEN_VDW] p_0: %8.5lf [kcal/mol], p_1: %8.5lf [A], shift: %8.5lf [kcal/mol], tailcorrection: %s\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KCAL_PER_MOL,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KCAL_PER_MOL,
              TailCorrection[i][j]?"yes":"no");
            break;
          case MATSUOKA_CLEMENTI_YOSHIMINE:
            // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K]
            // p_3     [A^-1]
            // p_4/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [MATSUOKA_CLEMENTI_YOSHIMINE] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K], "
                            "p_3: %8.5lf [A^-1], shift/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3],
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN);
            break;
          case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3:
            // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K]
            // p_3     [A^-1]
            fprintf(FilePtr,"%7s - %7s [MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], "
                            "p_2/k_B: %8.5lf [K], p_3: %8.5lf [A^-1]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]);
            break;
          case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5:
            // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K]
            // p_3     [A^-1]
            fprintf(FilePtr,"%7s - %7s [MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], "
                            "p_2/k_B: %8.5lf [K], p_3: %8.5lf [A^-1]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]);
            break;
          case GENERIC:
            // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^4]
            // p_3/k_B [K A^6]
            // p_4/k_B [K A^8]
            // p_5/k_B [K A^10]
            // p_6/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [GENERIC] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^4], p_3/k_B: %8.5lf [K A^6], "
                            "p_4/k_B: %8.5lf [K A^8], p_5/k_B: %8.5lf [K A^10], shift/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][5]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][6]*ENERGY_TO_KELVIN);
            break;
          case GENERIC_SMOOTHED3:
            // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^4]
            // p_3/k_B [K A^6]
            // p_4/k_B [K A^8]
            // p_5/k_B [K A^10]
            fprintf(FilePtr,"%7s - %7s [GENERIC_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^4], p_3/k_B: %8.5lf [K A^6], "
                            "p_4/k_B: %8.5lf [K A^8], p_5/k_B: %8.5lf [K A^10]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][5]*ENERGY_TO_KELVIN);
            break;
          case GENERIC_SMOOTHED5:
            // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^4]
            // p_3/k_B [K A^6]
            // p_4/k_B [K A^8]
            // p_5/k_B [K A^10]
            fprintf(FilePtr,"%7s - %7s [GENERIC_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^4], p_3/k_B: %8.5lf [K A^6], "
                            "p_4/k_B: %8.5lf [K A^8], p_5/k_B: %8.5lf [K A^10]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][5]*ENERGY_TO_KELVIN);
            break;
          case PELLENQ_NICHOLSON:
            // p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^6]
            // p_3/k_B [K A^8]
            // p_4/k_B [K A^10]
            // p_5/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [PELLENQ_NICHOLSON] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^6], "
                            "p_3/k_B: %8.5lf [K A^8], p_4/k_B: %8.5lf [K A^10], shift/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][5]*ENERGY_TO_KELVIN);
            break;
          case PELLENQ_NICHOLSON_SMOOTHED3:
            // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^6]
            // p_3/k_B [K A^8]
            // p_4/k_B [K A^10]
            fprintf(FilePtr,"%7s - %7s [PELLENQ_NICHOLSON_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], "
                            "p_2/k_B: %8.5lf [K A^6], p_3/k_B: %8.5lf [K A^8], p_4/k_B: %8.5lf [K A^10]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN);
            break;
          case PELLENQ_NICHOLSON_SMOOTHED5:
            // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^6]
            // p_3/k_B [K A^8]
            // p_4/k_B [K A^10]
            fprintf(FilePtr,"%7s - %7s [PELLENQ_NICHOLSON_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], "
                            "p_2/k_B: %8.5lf [K A^6], p_3/k_B: %8.5lf [K A^8], p_4/k_B: %8.5lf [K A^10]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN);
            break;
          case HYDRATED_ION_WATER:
            // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^4]
            // p_3/k_B [K A^6]
            // p_4/k_B [K A^12]
            // p_5/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [HYDRATED_ION_WATER] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], p_2/k_B: %8.5lf [K A^4], "
                            "p_3/k_B: %8.5lf [K A^6], p_4/k_B: %8.5lf [K A^12], shift/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][5]*ENERGY_TO_KELVIN);
            break;
          case HYDRATED_ION_WATER_SMOOTHED3:
            // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^4]
            // p_3/k_B [K A^6]
            // p_4/k_B [K A^12]
            fprintf(FilePtr,"%7s - %7s [HYDRATED_ION_WATER_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], "
                            "p_2/k_B: %8.5lf [K A^4], p_3/k_B: %8.5lf [K A^6], p_4/k_B: %8.5lf [K A^12]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN);
            break;
          case HYDRATED_ION_WATER_SMOOTHED5:
            // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [A^-1]
            // p_2/k_B [K A^4]
            // p_3/k_B [K A^6]
            // p_4/k_B [K A^12]
            fprintf(FilePtr,"%7s - %7s [HYDRATED_ION_WATER_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1: %8.5lf [A^-1], "
                            "p_2/k_B: %8.5lf [K A^4], p_3/k_B: %8.5lf [K A^6], p_4/k_B: %8.5lf [K A^12]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN);
            break;
          case MIE:
            // p_0*[p_1/r^p_2-p_1/r^p_3]
            // ======================================================================================
            // p_0/k_B [K]
            // p_1/k_B [K A^p_2]
            // p_2     [-]
            // p_3     [-]
            // p_4/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [MIE] p_0/k_B: %8.5lf [K], p_1/k_B: %8.5lf [K A^p_2], p_2: %8.5lf [-], "
                            "p_3: %8.5lf [-], shift/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3],
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN);
            break;
          case MIE_SMOOTHED3:
            // {p_0*[p_1/r^p_2-p_1/r^p_3]}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1/k_B [K A^p_2]
            // p_2     [-]
            // p_3     [-]
            fprintf(FilePtr,"%7s - %7s [MIE_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1/k_B: %8.5lf [K A^p_2], "
                            "p_2: %8.5lf [-], p_3: %8.5lf [-]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]);
            break;
          case MIE_SMOOTHED5:
            // {p_0*[p_1/r^p_2-p_1/r^p_3]}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1/k_B [K A^p_2]
            // p_2     [-]
            // p_3     [-]
            fprintf(FilePtr,"%7s - %7s [MIE_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1/k_B: %8.5lf [K A^p_2], "
                            "p_2: %8.5lf [-], p_3: %8.5lf [-]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3]);
            break;
          case MIE_CUTOFF:
            // p_0*[p_1/r^p_2-p_1/r^p_3] if r < p_4, otherwise 0
            // ======================================================================================
            // p_0/k_B [K]
            // p_1/k_B [K A^p_2]
            // p_2     [-]
            // p_3     [-]
            // p_4     [A]
            // p_5/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [MIE_CUTOFF] p_0/k_B: %8.5lf [K], p_1/k_B: %8.5lf [K A^p_2], p_2: %8.5lf [-], "
                            "p_3: %8.5lf [-], p_4: %8.5lf [A], shift/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3],
              (double)PotentialParms[i][j][4],
              (double)PotentialParms[i][j][5]*ENERGY_TO_KELVIN);
            break;
          case MIE_SMOOTHED3_CUTOFF:
            // {p_0*[p_1/r^p_2-p_1/r^p_3]}*S(r) if r < p_4, otherwise 0
            // ======================================================================================
            // p_0/k_B [K]
            // p_1/k_B [K A^p_2]
            // p_2     [-]
            // p_3     [-]
            // p_4     [A]
            fprintf(FilePtr,"%7s - %7s [MIE_SMOOTHED3_CUTOFF] p_0/k_B: %8.5lf [K], p_1/k_B: %8.5lf [K A^p_2], "
                            "p_2: %8.5lf [-], p_3: %8.5lf [-], p_4: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3],
              (double)PotentialParms[i][j][4]);
            break;
          case MIE_SMOOTHED5_CUTOFF:
            // {p_0*[p_1/r^p_2-p_1/r^p_3]}*S(r) if r < p_4, otherwise 0
            // ======================================================================================
            // p_0/k_B [K]
            // p_1/k_B [K A^p_2]
            // p_2     [-]
            // p_3     [-]
            // p_4     [A]
            fprintf(FilePtr,"%7s - %7s [MIE_SMOOTHED5_CUTOFF] p_0/k_B: %8.5lf [K], p_1/k_B: %8.5lf [K A^p_2], "
                            "p_2: %8.5lf [-], p_3: %8.5lf [-], p_4: %8.5lf [A]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][3],
              (double)PotentialParms[i][j][4]);
            break;
          case BORN_HUGGINS_MEYER:
            // p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [-]
            // p_2     [A]
            // p_3/k_B [K A^6]
            // p_4/k_B [K A^8]
            // p_5/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [BORN_HUGGINS_MEYER] p_0/k_B: %8.5lf [K], p_1: %8.5lf [-], p_2: %8.5lf [A], "
                            "p_3/k_B: %8.5lf [K A^6], p_4/k_B: %8.5lf [K A^8], shift/k_B: %8.5lf\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2],
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][5]*ENERGY_TO_KELVIN);
            break;
          case BORN_HUGGINS_MEYER_SMOOTHED3:
            // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [-]
            // p_2     [A]
            // p_3/k_B [K A^6]
            // p_4/k_B [K A^8]
            fprintf(FilePtr,"%7s - %7s [BORN_HUGGINS_MEYER_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1: %8.5lf [-], p_2: %8.5lf [A], "
                            "p_3/k_B: %8.5lf [K A^6], p_4/k_B: %8.5lf [K A^8]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2],
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN);
            break;
          case BORN_HUGGINS_MEYER_SMOOTHED5:
            // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
            // ======================================================================================
            // p_0/k_B [K]
            // p_1     [-]
            // p_2     [A]
            // p_3/k_B [K A^6]
            // p_4/k_B [K A^8]
            fprintf(FilePtr,"%7s - %7s [BORN_HUGGINS_MEYER_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1: %8.5lf [-], p_2: %8.5lf [A], "
                            "p_3/k_B: %8.5lf [K A^6], p_4/k_B: %8.5lf [K A^8]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1],
              (double)PotentialParms[i][j][2],
              (double)PotentialParms[i][j][3]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][4]*ENERGY_TO_KELVIN);
            break;
          case HYDROGEN:
            // p_0/r^12-p_1/r^10
            // ======================================================================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^10]
            // p_2/k_B [K]  (non-zero for a shifted potential)
            fprintf(FilePtr,"%7s - %7s [HYDROGEN]  p_0/k_B: %8.5lf [K], p_1/k_B: 8.5 %lf [K], shift/k_B: %8.5lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][2]*ENERGY_TO_KELVIN);
            break;
          case HYDROGEN_SMOOTHED3:
            // {p_0/r^12-p_1/r^10}*S(r)
            // ======================================================================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^10]
            fprintf(FilePtr,"%7s - %7s [HYDROGEN_SMOOTHED3] p_0/k_B: %8.5lf [K], p_1/k_B: 8.5 %lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN);
            break;
          case HYDROGEN_SMOOTHED5:
            // {p_0/r^12-p_1/r^10}*S(r)
            // ======================================================================================
            // p_0/k_B [K A^12]
            // p_1/k_B [K A^10]
            fprintf(FilePtr,"%7s - %7s [HYDROGEN_SMOOTHED5] p_0/k_B: %8.5lf [K], p_1/k_B: 8.5 %lf [K]\n",
              PseudoAtoms[i].Name,
              PseudoAtoms[j].Name,
              (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
              (double)PotentialParms[i][j][1]*ENERGY_TO_KELVIN);
            break;
        }
      }
    fprintf(FilePtr,"\n\n");
  }

  if(WriteMoleculeDefinitionToOutput)
  {
    fprintf(FilePtr,"MoleculeDefinitions:\n");
    fprintf(FilePtr,"===========================================================================\n");

    if(NumberOfComponents==0)
       fprintf(FilePtr,"No components defined in 'simulation.input'\n");

    for(i=0;i<NumberOfComponents;i++)
    {
      if(Components[i].ExtraFrameworkMolecule)
        fprintf(FilePtr,"Component %d [%s] (Cation molecule)\n",i,Components[i].Name);
      else
        fprintf(FilePtr,"Component %d [%s] (Adsorbate molecule)\n",i,Components[i].Name);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tMoleculeDefinitions: %s\n",Components[i].MoleculeDefinition);
      fprintf(FilePtr,"\tComponent %s\n",Components[i].HasCharges?"contains (at least some) atoms which are charged":
          "contains no atoms with charge");
      fprintf(FilePtr,"\tComponent %s\n",Components[i].IsPolarizable?"contains (at least some) atoms which are polarizable":
          "contains no atoms with point dipoles (polarization)");

      NetCharge=0.0;
      for(j=0;j<Components[i].NumberOfAtoms;j++)
        NetCharge+=Components[i].Charge[j];
      fprintf(FilePtr,"\tComponent has a net charge of %lf\n",(double)NetCharge);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tIdeal chain Rosenbluth weight: %lg\n",(double)Components[i].IdealGasRosenbluthWeight[system]);
      fprintf(FilePtr,"\tIdeal chain total energy: %lf\n",(double)Components[i].IdealGasTotalEnergy[system]);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tCritical temparure [K]: %lf\n",(double)Components[i].CriticalTemperature);
      fprintf(FilePtr,"\tCritical pressure [Pa]: %lf\n",(double)Components[i].CriticalPressure);
      fprintf(FilePtr,"\tAcentric factor [-]: %lf\n",(double)Components[i].AcentricFactor);
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"\tRXMC partition factor [-]: %lf\n",(double)Components[i].PartitionFunction);
      fprintf(FilePtr,"\n");

      switch(FluidState[system])
      {
        case VAPOUR_STABLE:
          fprintf(FilePtr,"\tVapour=stable, Liquid=metastable\n");
          break;
        case LIQUID_STABLE:
          fprintf(FilePtr,"\tLiquid=stable, Vapour=metastable\n");
          break;
        case VAPOUR_LIQUID_STABLE:
          fprintf(FilePtr,"\tLiquid=stable, Vapour=stable\n");
          break;
        case SUPER_CRITICAL_FLUID:
          fprintf(FilePtr,"\tFluid is supercritical\n");
          break;
        case LIQUID:
          fprintf(FilePtr,"\tFluid is a liquid\n");
          break;
        case VAPOUR:
          fprintf(FilePtr,"\tFluid is a vapour\n");
          break;
      }
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tMolFraction:     %18.10lf [-]\n",(double)Components[i].MolFraction[system]);
      fprintf(FilePtr,"\tCompressibility: %18.10lf [-]\n",(double)Components[i].Compressibility[system]);
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"\tDensity of the bulk fluid phase: %18.10lf [kg/m^3]\n",
           (double)Components[i].BulkFluidDensity[system]);
      fprintf(FilePtr,"\n");

      if(Components[i].Swapable)
      {
        fprintf(FilePtr,"\tBinary mixture EOS parameters: ");
          for(j=0;j<NumberOfComponents;j++)
           if(Components[j].Swapable)
             fprintf(FilePtr," (%d): %lf",j,BinaryInteractionParameter[i][j]);
           else
             fprintf(FilePtr," (%d): -",j);
        fprintf(FilePtr,"\n");
      }
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tAmount of excess molecules: %18.10lf [-]\n",
                 (double)Components[i].AmountOfExcessMolecules[system]);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tConversion factor molecules/unit cell -> mol/kg: %18.10lf [-]\n",
                 (double)Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[system]);
      fprintf(FilePtr,"\tConversion factor molecules/unit cell -> gr/gr: %18.10lf [-]\n",
                 (double)Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[system]);
      fprintf(FilePtr,"\tConversion factor molecules/unit cell -> cm^3 STP/gr: %18.10lf [-]\n",
                 (double)Components[i].MOLEC_PER_UC_TO_CC_STP_G[system]);
      fprintf(FilePtr,"\tConversion factor molecules/unit cell -> cm^3 STP/cm^3: %18.10lf [-]\n",
                 (double)Components[i].MOLEC_PER_UC_TO_CC_STP_CC[system]);
      fprintf(FilePtr,"\tConversion factor mol/kg -> cm^3 STP/gr: %18.10lf [-]\n",
                 (double)Components[i].MOL_PER_KG_TO_CC_STP_G[system]);
      fprintf(FilePtr,"\tConversion factor mol/kg -> cm^3 STP/cm^3: %18.10lf [-]\n",
                 (double)Components[i].MOL_PER_KG_TO_CC_STP_CC[system]);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tPartial pressure: %24.14lf [Pa]\n",(double)Components[i].PartialPressure[system]*PRESSURE_CONVERSION_FACTOR);
      fprintf(FilePtr,"\t                  %24.14lf [Torr]\n",(double)
              Components[i].PartialPressure[system]*PRESSURE_CONVERSION_FACTOR*PA_TO_TORR);
      fprintf(FilePtr,"\t                  %24.14lf [bar]\n",(double)
              Components[i].PartialPressure[system]*PRESSURE_CONVERSION_FACTOR*PA_TO_BAR);
      fprintf(FilePtr,"\t                  %24.14lf [atm]\n",(double)
              Components[i].PartialPressure[system]*PRESSURE_CONVERSION_FACTOR*PA_TO_ATM);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tFugacity coefficient: %18.10lf [-]\n",(double)Components[i].FugacityCoefficient[system]);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tPartial fugacity: %24.14lf [Pa]\n",(double)(Components[i].FugacityCoefficient[system]*
           Components[i].PartialPressure[system]*PRESSURE_CONVERSION_FACTOR));
      fprintf(FilePtr,"\t                  %24.14lf [Torr]\n",(double)(Components[i].FugacityCoefficient[system]*
              Components[i].PartialPressure[system]*PRESSURE_CONVERSION_FACTOR*PA_TO_TORR));
      fprintf(FilePtr,"\t                  %24.14lf [bar]\n",(double)(Components[i].FugacityCoefficient[system]*
              Components[i].PartialPressure[system]*PRESSURE_CONVERSION_FACTOR*PA_TO_BAR));
      fprintf(FilePtr,"\t                  %24.14lf [atm]\n",(double)(Components[i].FugacityCoefficient[system]*
              Components[i].PartialPressure[system]*PRESSURE_CONVERSION_FACTOR*PA_TO_ATM));
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tMolecule contains %d number of atoms\n",Components[i].NumberOfAtoms);
      for(j=0;j<Components[i].NumberOfAtoms;j++)
      {
          fprintf(FilePtr,"\t\tatom: %4d  is of type: %4d [%10s] (group: %d)%s\n",
            j,Components[i].Type[j],PseudoAtoms[Components[i].Type[j]].Name,Components[i].group[j],
            Components[i].Fixed[j]?" kept fixed in minimization":"");
      }
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tMolecule contains %d chirality centers\n",Components[i].NumberOfChiralityCenters);
      for(j=0;j<Components[i].NumberOfAtoms;j++)
      {
        if(Components[i].NumberOfChiralityCenters>0)
        {
          if(Components[i].Chirality[j])
          {
             fprintf(FilePtr,"\t\tCenter around: %4d   with A=%d B=%d C=%d D=%d [%s]\n",
               Components[i].ChiralB[j],
               Components[i].ChiralA[j],
               Components[i].ChiralB[j],
               Components[i].ChiralC[j],
               Components[i].ChiralD[j],
               Components[i].ChiralityType[j]==S_CHIRAL?"S-chirality":"R-charility");
          }
        }
      }
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tMolecule contains %d number of groups\n",Components[i].NumberOfGroups);
      for(j=0;j<Components[i].NumberOfGroups;j++)
      {
        fprintf(FilePtr,"\n");
        fprintf(FilePtr,"\t\tgroup: %d containing: %d elements\n",j,Components[i].Groups[j].NumberOfGroupAtoms);
        fprintf(FilePtr,"\t\t-------------------------------------------------\n");
        switch(Components[i].Groups[j].Rigid)
        {
          case TRUE:
            switch(Components[i].Groups[j].Type)
            {
              case POINT_PARTICLE:
                fprintf(FilePtr,"\t\tthe group represents a point particle\n");
                break;
              case LINEAR_MOLECULE:
                fprintf(FilePtr,"\t\tthe group is rigid and linear\n");
                break;
              case NONLINEAR_MOLECULE:
                fprintf(FilePtr,"\t\tthe group is rigid and nonlinear\n");
                break;
            }
            break;
          case FALSE:
            fprintf(FilePtr,"\t\tthe group is modelled as flexible, i.e. no constraints\n");
            break;
        }
        fprintf(FilePtr,"\t\tMass: %lf [a.u.]\n",(double)Components[i].Groups[j].Mass);
        fprintf(FilePtr,"\t\tMass: %lf [kg/m^3]\n\n",(double)(Components[i].Groups[j].Mass/(Volume[0]*6.02214e-4)));
        fprintf(FilePtr,"\t\tRotational Degrees of freedom: %d\n",
                Components[i].Groups[j].RotationalDegreesOfFreedom);
        fprintf(FilePtr,"\t\tDiagonalized inertia vector: %18.10lf\n",(double)Components[i].Groups[j].InertiaVector.x);
        fprintf(FilePtr,"\t\t                             %18.10lf\n",(double)Components[i].Groups[j].InertiaVector.y);
        fprintf(FilePtr,"\t\t                             %18.10lf\n",(double)Components[i].Groups[j].InertiaVector.z);

        if(Components[i].Groups[j].Rigid)
        {
          fprintf(FilePtr,"\t\tnumber of atoms: %d\n",Components[i].Groups[j].NumberOfGroupAtoms);
          for(k=0;k<Components[i].Groups[j].NumberOfGroupAtoms;k++)
          {
            A=Components[i].Groups[j].Atoms[k];
            fprintf(FilePtr,"\t\t\telement: %d atom: %d [%10s] Charge: % lf Anisotropy: %lf Position: % lf % lf % lf Connectivity: %d (",
              k,
              A,
              PseudoAtoms[Components[i].Type[A]].Name,
              Components[i].Charge[A],
              PseudoAtoms[Components[i].Type[A]].AnisotropicDisplacement,
              (double)Components[i].Positions[A].x,
              (double)Components[i].Positions[A].y,
              (double)Components[i].Positions[A].z,
              Components[i].Connectivity[A]);
            for(l=0;l<Components[i].Connectivity[A];l++)
              fprintf(FilePtr,"%d ",Components[i].ConnectivityList[A][l]);
            fprintf(FilePtr,")\n");
          }

          fprintf(FilePtr,"\t\tnumber of permanent dipoles: %d\n",Components[i].Groups[j].NumberOfPermanentDipoles);
          for(k=0;k<Components[i].Groups[j].NumberOfPermanentDipoles;k++)
          {
             fprintf(FilePtr,"\t\t\tposition: %18.10lf %18.10lf %18.10lf, permanent dipole: %18.10lf %18.10lf %18.10lf\n",
                Components[i].Groups[j].PermanentDipolePositions[k].x,
                Components[i].Groups[j].PermanentDipolePositions[k].y,
                Components[i].Groups[j].PermanentDipolePositions[k].z,
                Components[i].Groups[j].PermanentDipoles[k].x*DEBYE_CONVERSION_FACTOR,
                Components[i].Groups[j].PermanentDipoles[k].y*DEBYE_CONVERSION_FACTOR,
                Components[i].Groups[j].PermanentDipoles[k].z*DEBYE_CONVERSION_FACTOR);
          }

          fprintf(FilePtr,"\t\tnumber of polarizabilities: %d\n",Components[i].Groups[j].NumberOfPolarizabilities);
          for(k=0;k<Components[i].Groups[j].NumberOfPolarizabilities;k++)
          {
             fprintf(FilePtr,"\t\t\tposition: %18.10lf %18.10lf %18.10lf, polarizability: %18.10lf %18.10lf %18.10lf\n",
                Components[i].Groups[j].PolarizabilityPositions[k].x,
                Components[i].Groups[j].PolarizabilityPositions[k].y,
                Components[i].Groups[j].PolarizabilityPositions[k].z,
                Components[i].Groups[j].Polarizabilites[k].ax*COULOMBIC_CONVERSION_FACTOR,
                Components[i].Groups[j].Polarizabilites[k].ay*COULOMBIC_CONVERSION_FACTOR,
                Components[i].Groups[j].Polarizabilites[k].az*COULOMBIC_CONVERSION_FACTOR);
             fprintf(FilePtr,"\t\t\t                                                                                    %18.10lf %18.10lf %18.10lf\n",
                Components[i].Groups[j].Polarizabilites[k].bx*COULOMBIC_CONVERSION_FACTOR,
                Components[i].Groups[j].Polarizabilites[k].by*COULOMBIC_CONVERSION_FACTOR,
                Components[i].Groups[j].Polarizabilites[k].bz*COULOMBIC_CONVERSION_FACTOR);
             fprintf(FilePtr,"\t\t\t                                                                                    %18.10lf %18.10lf %18.10lf\n",
                Components[i].Groups[j].Polarizabilites[k].cx*COULOMBIC_CONVERSION_FACTOR,
                Components[i].Groups[j].Polarizabilites[k].cy*COULOMBIC_CONVERSION_FACTOR,
                Components[i].Groups[j].Polarizabilites[k].cz*COULOMBIC_CONVERSION_FACTOR);
          }
        }
        else
        {
          for(k=0;k<Components[i].Groups[j].NumberOfGroupAtoms;k++)
          {
            A=Components[i].Groups[j].Atoms[k];
            fprintf(FilePtr,"\t\telement: %d atom: %d [%10s] Charge: % lf Anisotropy: %lf Connectivity: %d (",
              k,
              A,
              PseudoAtoms[Components[i].Type[A]].Name,
              Components[i].Charge[A],
              PseudoAtoms[Components[i].Type[A]].AnisotropicDisplacement,
              Components[i].Connectivity[A]);
            for(l=0;l<Components[i].Connectivity[A];l++)
              fprintf(FilePtr,"%d ",Components[i].ConnectivityList[A][l]);
            fprintf(FilePtr,")\n");
          }
        }
        fprintf(FilePtr,"\n");

        Dipole=ComputeDipoleMomentComponent(i,j);
        fprintf(FilePtr,"\t\tDipole:     %18.10lf [D]\n",
                (double)sqrt(SQR(Dipole.x)+SQR(Dipole.y)+SQR(Dipole.z))*DEBYE_CONVERSION_FACTOR);

        Quadrupole=ComputeQuadrupoleMomentComponent(i,j);
        EigenSystem3x3(Quadrupole,&Eigenvectors,&eigenvalues);
        fprintf(FilePtr,"\t\tQuadrupole: %18.10lf %18.10lf %18.10lf [D Angstrom]\n",
              eigenvalues.x*DEBYE_CONVERSION_FACTOR,eigenvalues.y*DEBYE_CONVERSION_FACTOR,eigenvalues.z*DEBYE_CONVERSION_FACTOR);
        fprintf(FilePtr,"\t\tQuadrupole tensor [D Angstrom]\n");
        fprintf(FilePtr,"\t\t\t\t %18.10lf %18.10lf %18.10lf\n",
              Quadrupole.ax*DEBYE_CONVERSION_FACTOR,Quadrupole.bx*DEBYE_CONVERSION_FACTOR,Quadrupole.cx*DEBYE_CONVERSION_FACTOR);
        fprintf(FilePtr,"\t\t\t\t %18.10lf %18.10lf %18.10lf\n",
              Quadrupole.ay*DEBYE_CONVERSION_FACTOR,Quadrupole.by*DEBYE_CONVERSION_FACTOR,Quadrupole.cy*DEBYE_CONVERSION_FACTOR);
        fprintf(FilePtr,"\t\t\t\t %18.10lf %18.10lf %18.10lf\n",
              Quadrupole.az*DEBYE_CONVERSION_FACTOR,Quadrupole.bz*DEBYE_CONVERSION_FACTOR,Quadrupole.cz*DEBYE_CONVERSION_FACTOR);
      }
      fprintf(FilePtr,"\n");


      fprintf(FilePtr,"\tStarting bead for growth               : %d\n",Components[i].StartingBead);
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"\tDegrees of freedom                     : %d\n",Components[i].DegreesOfFreedom);
      fprintf(FilePtr,"\tTranslational degrees of freedom       : %d\n",Components[i].TranslationalDegreesOfFreedom);
      fprintf(FilePtr,"\tRotational degrees of freedom          : %d\n",Components[i].RotationalDegreesOfFreedom);
      fprintf(FilePtr,"\tVibrational degrees of freedom         : %d\n",Components[i].VibrationalDegreesOfFreedom);
      fprintf(FilePtr,"\tConstraint degrees of freedom          : %d\n",Components[i].ConstraintDegreesOfFreedom);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tNumber of atoms                        : %d\n",Components[i].NumberOfAtoms);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tNumber of constraint bonds             : %d\n",Components[i].NumberOfConstraintBonds);
      fprintf(FilePtr,"\tNumber of constraint bends             : %d\n",Components[i].NumberOfConstraintBends);
      fprintf(FilePtr,"\tNumber of constraint inversion bends   : %d\n",Components[i].NumberOfConstraintInversionBends);
      fprintf(FilePtr,"\tNumber of constraint torsions          : %d\n",Components[i].NumberOfConstraintTorsions);
      fprintf(FilePtr,"\tNumber of constraint improper torsions : %d\n",Components[i].NumberOfConstraintImproperTorsions);
      fprintf(FilePtr,"\tNumber of constraint improper torsions : %d\n",Components[i].NumberOfConstraintOutOfPlanes);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tNumber of bonds                        : %d\n",Components[i].NumberOfBonds);
      fprintf(FilePtr,"\tNumber of Urey-Bradleys                : %d\n",Components[i].NumberOfUreyBradleys);
      fprintf(FilePtr,"\tNumber of bends                        : %d\n",Components[i].NumberOfBends);
      fprintf(FilePtr,"\tNumber of inversion bends              : %d\n",Components[i].NumberOfInversionBends);
      fprintf(FilePtr,"\tNumber of torsions                     : %d\n",Components[i].NumberOfTorsions);
      fprintf(FilePtr,"\tNumber of improper torsions            : %d\n",Components[i].NumberOfImproperTorsions);
      fprintf(FilePtr,"\tNumber of improper torsions            : %d\n",Components[i].NumberOfOutOfPlanes);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tNumber of bond/bond cross terms        : %d\n",Components[i].NumberOfBondBonds);
      fprintf(FilePtr,"\tNumber of bond/bend cross terms        : %d\n",Components[i].NumberOfBondBends);
      fprintf(FilePtr,"\tNumber of bend/bend cross terms        : %d\n",Components[i].NumberOfBendBends);
      fprintf(FilePtr,"\tNumber of stretch/torsion cross terms  : %d\n",Components[i].NumberOfBondTorsions);
      fprintf(FilePtr,"\tNumber of bend/torsion cross terms     : %d\n",Components[i].NumberOfBendTorsions);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tNumber of charges                      : %d\n",Components[i].NumberOfCharges);
      fprintf(FilePtr,"\tNumber of bond-dipoles                 : %d\n",Components[i].NumberOfBondDipoles);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tNumber of intra Van der Waals                             : %d\n",Components[i].NumberOfIntraVDW);
      fprintf(FilePtr,"\tNumber of intra charge-charge Coulomb                     : %d\n",Components[i].NumberOfIntraChargeCharge);
      fprintf(FilePtr,"\tNumber of intra charge-bonddipole Coulomb                 : %d\n",Components[i].NumberOfIntraChargeBondDipole);
      fprintf(FilePtr,"\tNumber of intra bonddipole-bonddipole Coulomb             : %d\n",Components[i].NumberOfIntraBondDipoleBondDipole);
      fprintf(FilePtr,"\n");

      // test
      fprintf(FilePtr,"\tNumber of excluded intra charge-charge Coulomb                     : %d\n",Components[i].NumberOfExcludedIntraChargeCharge);
      fprintf(FilePtr,"\tNumber of excluded intra charge-bonddipole Coulomb                 : %d\n",Components[i].NumberOfExcludedIntraChargeBondDipole);
      fprintf(FilePtr,"\tNumber of excluded intra bonddipole-bonddipole Coulomb             : %d\n",Components[i].NumberOfExcludedIntraBondDipoleBondDipole);
      fprintf(FilePtr,"\n");


      fprintf(FilePtr,"\tNumber of cbmc-config moves                               : %d\n",Components[i].NumberOfConfigMoves);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tParticle Moves:             \n");
      fprintf(FilePtr,"\t\tProbabilityTranslationMove:                  %lf\n",(double)(100.0*Components[i].ProbabilityTranslationMove));
      switch(Components[i].TranslationDirection)
      {
        case X_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      X\n");
          break;
        case Y_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      Y\n");
          break;
        case Z_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      Z\n");
          break;
        case XY_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      XY\n");
          break;
        case XZ_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      XZ\n");
          break;
        case YZ_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      YZ\n");
          break;
        case XYZ_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      XYZ\n");
          break;
        case A_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      A\n");
          break;
        case B_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      B\n");
          break;
        case C_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      C\n");
          break;
        case AB_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      AB\n");
          break;
        case AC_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      AC\n");
          break;
        case BC_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      BC\n");
          break;
        case ABC_DIR:
          fprintf(FilePtr,"\t\t\tTranslationDirection:      ABC\n");
          break;

      }
      fprintf(FilePtr,"\t\tPercentage of random translation moves:            %lf\n",(double)(100.0*Components[i].FractionOfRandomTranslationMove));
      fprintf(FilePtr,"\t\tPercentage of rotation moves:                      %lf\n",(double)(100.0*Components[i].FractionOfRotationMove));
      fprintf(FilePtr,"\t\tPercentage of random rotation moves:               %lf\n",(double)(100.0*Components[i].FractionOfRandomRotationMove));
      fprintf(FilePtr,"\t\tPercentage of partial reinsertion moves:           %lf\n",(double)(100.0*Components[i].FractionOfPartialReinsertionMove));
      fprintf(FilePtr,"\t\tPercentage of reinsertion moves:                   %lf\n",(double)(100.0*Components[i].FractionOfReinsertionMove));
      fprintf(FilePtr,"\t\tPercentage of reinsertion-in-place moves:          %lf\n",(double)(100.0*Components[i].FractionOfReinsertionInPlaceMove));
      fprintf(FilePtr,"\t\tPercentage of reinsertion-in-plane moves:          %lf\n",(double)(100.0*Components[i].FractionOfReinsertionInPlaneMove));
      fprintf(FilePtr,"\t\tPercentage of identity-change moves:               %lf\n",(double)(100.0*Components[i].FractionOfIdentityChangeMove));
      for(j=0;j<Components[i].NumberOfIdentityChanges;j++)
        fprintf(FilePtr,"\t\t\tmove %d    component %d => %d\n",j,i,Components[i].IdentityChanges[j]);
      fprintf(FilePtr,"\t\tPercentage of swap (insert/delete) moves:          %lf\n",(double)(100.0*Components[i].FractionOfSwapMove));
      fprintf(FilePtr,"\t\tPercentage of CF swap lambda moves:                %lf\n",(double)(100.0*Components[i].FractionOfCFSwapLambdaMove));
      fprintf(FilePtr,"\t\tPercentage of CB/CFMC swap lambda moves:           %lf\n",(double)(100.0*Components[i].FractionOfCBCFSwapLambdaMove));
      fprintf(FilePtr,"\t\tPercentage of Widom insertion moves:               %lf\n",(double)(100.0*Components[i].FractionOfWidomMove));
      fprintf(FilePtr,"\t\tPercentage of surface-area moves:                  %lf\n",(double)(100.0*Components[i].FractionOfSurfaceAreaMove));
      fprintf(FilePtr,"\t\tPercentage of Gibbs particle-transfer moves:       %lf\n",(double)(100.0*Components[i].FractionOfGibbsChangeMove));
      fprintf(FilePtr,"\t\tPercentage of Gibbs identity-change moves:         %lf\n",(double)(100.0*Components[i].FractionOfGibbsIdentityChangeMove));
      for(j=0;j<Components[i].NumberOfGibbsIdentityChanges;j++)
        fprintf(FilePtr,"\t\t\tmove %d    component %d => %d\n",j,i,Components[i].GibbsIdentityChanges[j]);
      fprintf(FilePtr,"\t\tPercentage of CF Gibbs lambda-transfer moves:      %lf\n",(double)(100.0*Components[i].FractionOfCFGibbsChangeMove));
      fprintf(FilePtr,"\t\tPercentage of CB/CFMC Gibbs lambda-transfer moves: %lf\n",(double)(100.0*Components[i].FractionOfCBCFGibbsChangeMove));
      fprintf(FilePtr,"\t\tPercentage of exchange frac./int. particle moves:  %lf\n",(double)(100.0*Components[i].FractionOfExchangeFractionalParticleMove));
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tSystem Moves:\n");
      fprintf(FilePtr,"\t\tPercentage of parallel-tempering moves:            %lf\n",(double)(100.0*Components[i].FractionOfParallelTemperingMove));
      fprintf(FilePtr,"\t\tPercentage of hyper-parallel-tempering moves:      %lf\n",(double)(100.0*Components[i].FractionOfHyperParallelTemperingMove));
      fprintf(FilePtr,"\t\tPercentage of parallel-mol-fraction moves:         %lf\n",(double)(100.0*Components[i].FractionOfParallelMolFractionMove));
      fprintf(FilePtr,"\t\t\t   Component A: %d B: %d\n",ParallelMolFractionComponentA,ParallelMolFractionComponentB);
      fprintf(FilePtr,"\t\tPercentage of chiral inversion moves:              %lf\n",(double)(100.0*Components[i].FractionOfChiralInversionMove));
      fprintf(FilePtr,"\t\tPercentage of Hybrid-NVE moves:                    %lf\n",(double)(100.0*Components[i].FractionOfHybridNVEMove));
      fprintf(FilePtr,"\t\tPercentage of Hybrid-NPH moves:                    %lf\n",(double)(100.0*Components[i].FractionOfHybridNPHMove));
      fprintf(FilePtr,"\t\tPercentage of Hybrid-NPHPR moves:                  %lf\n",(double)(100.0*Components[i].FractionOfHybridNPHPRMove));
      fprintf(FilePtr,"\t\tPercentage of volume-change moves:                 %lf\n",(double)(100.0*Components[i].FractionOfVolumeChangeMove));
      fprintf(FilePtr,"\t\tPercentage of box-shape-change moves:              %lf\n",(double)(100.0*Components[i].FractionOfBoxShapeChangeMove));
      fprintf(FilePtr,"\t\tPercentage of Gibbs volume-change moves:           %lf\n",(double)(100.0*Components[i].FractionOfGibbsVolumeChangeMove));
      fprintf(FilePtr,"\t\tPercentage of framework-change moves:              %lf\n",(double)(100.0*Components[i].FractionOfFrameworkChangeMove));
      fprintf(FilePtr,"\t\tPercentage of framework-shift moves:               %lf\n",(double)(100.0*Components[i].FractionOfFrameworkShiftMove));
      fprintf(FilePtr,"\t\tPercentage of reactive MC moves:                   %lf\n",(double)(100.0*Components[i].FractionOfCFCRXMCLambdaChangeMove));
      fprintf(FilePtr,"\n");


      if(Components[i].RestrictEnantionface)
      {
        fprintf(FilePtr,"\tMC Moves are restricted to leaving the enantioface unchanged : Yes\n");
        fprintf(FilePtr,"\tEnantionface: %s\n",Components[i].Enantioface==ENANTIOFACE_RE?"Re":"Si");
        switch(Components[i].EnantiofaceAtomDefinitions[0][0])
        {
          case FRAMEWORK:
            sprintf(buffer1,"Framework %d atom %d",Components[i].EnantiofaceAtomDefinitions[0][1],Components[i].EnantiofaceAtomDefinitions[0][2]);
            break;
          case ADSORBATE:
            sprintf(buffer1,"Adsorbate molecule %d atom %d",Components[i].EnantiofaceAtomDefinitions[0][1],Components[i].EnantiofaceAtomDefinitions[0][2]);
            break;
          case CATION:
            sprintf(buffer1,"Cation molecule %d atom %d",Components[i].EnantiofaceAtomDefinitions[0][1],Components[i].EnantiofaceAtomDefinitions[0][2]);
            break;
          default:
            break;
        }
        switch(Components[i].EnantiofaceAtomDefinitions[1][0])
        {
          case FRAMEWORK:
            sprintf(buffer2,"Framework %d atom %d",Components[i].EnantiofaceAtomDefinitions[1][1],Components[i].EnantiofaceAtomDefinitions[1][2]);
            break;
          case ADSORBATE:
            sprintf(buffer2,"Adsorbate molecule %d atom %d",Components[i].EnantiofaceAtomDefinitions[1][1],Components[i].EnantiofaceAtomDefinitions[1][2]);
            break;
          case CATION:
            sprintf(buffer2,"Cation molecule %d atom %d",Components[i].EnantiofaceAtomDefinitions[1][1],Components[i].EnantiofaceAtomDefinitions[1][2]);
            break;
          default:
            break;
        }
        switch(Components[i].EnantiofaceAtomDefinitions[2][0])
        {
          case FRAMEWORK:
            sprintf(buffer3,"Framework %d atom %d",Components[i].EnantiofaceAtomDefinitions[2][1],Components[i].EnantiofaceAtomDefinitions[2][2]);
            break;
          case ADSORBATE:
            sprintf(buffer3,"Adsorbate molecule %d atom %d",Components[i].EnantiofaceAtomDefinitions[2][1],Components[i].EnantiofaceAtomDefinitions[2][2]);
            break;
          case CATION:
            sprintf(buffer3,"Cation molecule %d atom %d",Components[i].EnantiofaceAtomDefinitions[2][1],Components[i].EnantiofaceAtomDefinitions[2][2]);
            break;
          default:
            break;
        }
        switch(Components[i].EnantiofaceAtomDefinitions[3][0])
        {
          case FRAMEWORK:
            sprintf(buffer4,"Framework %d atom %d",Components[i].EnantiofaceAtomDefinitions[3][1],Components[i].EnantiofaceAtomDefinitions[3][2]);
            break;
          case ADSORBATE:
            sprintf(buffer4,"Adsorbate molecule %d atom %d",Components[i].EnantiofaceAtomDefinitions[3][1],Components[i].EnantiofaceAtomDefinitions[3][2]);
            break;
          case CATION:
            sprintf(buffer4,"Cation molecule %d atom %d",Components[i].EnantiofaceAtomDefinitions[3][1],Components[i].EnantiofaceAtomDefinitions[3][2]);
            break;
          default:
            break;
        }
        switch(Components[i].EnantiofaceAtomDefinitions[4][0])
        {
          case FRAMEWORK:
            sprintf(buffer5,"Framework %d atom %d",Components[i].EnantiofaceAtomDefinitions[4][1],Components[i].EnantiofaceAtomDefinitions[4][2]);
            break;
          case ADSORBATE:
            sprintf(buffer5,"Adsorbate molecule %d atom %d",Components[i].EnantiofaceAtomDefinitions[4][1],Components[i].EnantiofaceAtomDefinitions[4][2]);
            break;
          case CATION:
            sprintf(buffer5,"Cation molecule %d atom %d",Components[i].EnantiofaceAtomDefinitions[4][1],Components[i].EnantiofaceAtomDefinitions[4][2]);
            break;
          default:
            break;
        }
        fprintf(FilePtr,"\tEnantionface atom definition: (%s,%s,%s,\n",buffer1,buffer2,buffer3);
        fprintf(FilePtr,"\t                               [%s,%s])",buffer4,buffer5);

        fprintf(FilePtr,"\n\n");
      }

      if(Components[i].RestrictMoves)
      {
        fprintf(FilePtr,"\tMoves are restricted: Yes\n\n");
        if(Components[i].RestrictMovesToBox)
        {
          fprintf(FilePtr,"\tMoves are restricted to a box: Yes\n");
          fprintf(FilePtr,"\tBox-a between :                  %lf-%lf\n",
                          (double)Components[i].BoxAxisABC_Min.x,
                          (double)Components[i].BoxAxisABC_Max.x);
          fprintf(FilePtr,"\tBox-b between :                  %lf-%lf\n",
                          (double)Components[i].BoxAxisABC_Min.y,
                          (double)Components[i].BoxAxisABC_Max.y);
          fprintf(FilePtr,"\tBox-c between :                  %lf-%lf\n\n",
                          (double)Components[i].BoxAxisABC_Min.z,
                          (double)Components[i].BoxAxisABC_Max.z);
        }
        else
          fprintf(FilePtr,"\tMoves are restricted to boxes: No\n\n");

        if(Components[i].RestrictMovesToPrisms)
        {
          for(j=0;j<MAX_NUMBER_OF_PRISMS;j++)
          {
            if(Components[i].RestrictMovesToPrism[j])
            {
              fprintf(FilePtr,"\tMoves are restricted to prism %d: Yes\n",j);
              fprintf(FilePtr,"\t\tprisms a between :                  %lf-%lf\n",
                              (double)Components[i].RestrictPrismABC_Min[j].x,
                              (double)Components[i].RestrictPrismABC_Max[j].x);
              fprintf(FilePtr,"\t\tprisms b between :                  %lf-%lf\n",
                              (double)Components[i].RestrictPrismABC_Min[j].y,
                              (double)Components[i].RestrictPrismABC_Max[j].y);
              fprintf(FilePtr,"\t\tprisms c between :                  %lf-%lf\n",
                              (double)Components[i].RestrictPrismABC_Min[j].z,
                              (double)Components[i].RestrictPrismABC_Max[j].z);
            }
          }
          fprintf(FilePtr,"\n");
        }
        else
          fprintf(FilePtr,"\tMoves are restricted to prisms: No\n\n");

        if(Components[i].RestrictMovesToCylinders)
        {
          for(j=0;j<MAX_NUMBER_OF_CYLINDERS;j++)
          {
            if(Components[i].RestrictMovesToCylinder[j])
            {
              fprintf(FilePtr,"\tMoves are restricted to cylinder %d: Yes\n",j);
              fprintf(FilePtr,"\t\tcylinder between :                  %lf-%lf\n",
                              (double)Components[i].RestrictCylinderABC_Min[j].x,
                              (double)Components[i].RestrictCylinderABC_Max[j].x);
              fprintf(FilePtr,"\t\tcylinder center :                   %lf,%lf,%lf\n",Components[i].RestrictCylinderCenter[j].x,Components[i].RestrictCylinderCenter[j].y,Components[i].RestrictCylinderCenter[j].z);
              fprintf(FilePtr,"\t\tcylinder direction :                %s\n",Components[i].RestrictCylinderDirection[j]==X_DIR?"X":(Components[i].RestrictCylinderDirection[j]==Y_DIR?"Y":"Z"));
              fprintf(FilePtr,"\t\tcylinder radius :                   %lf\n",Components[i].RestrictCylinderRadius[j]);
            }
          }
          fprintf(FilePtr,"\n");
        }
        else
          fprintf(FilePtr,"\tMoves are restricted to cylinders: No\n\n");

        if(Components[i].RestrictMovesToSpheres)
        {
          for(j=0;j<MAX_NUMBER_OF_SPHERES;j++)
          {
            if(Components[i].RestrictMovesToSphere[j])
            {
              fprintf(FilePtr,"\tMoves are restricted to sphere %d: Yes\n",j);
              fprintf(FilePtr,"\t\tsphere center :                   %lf,%lf,%lf\n",Components[i].RestrictSphereCenter[j].x,Components[i].RestrictSphereCenter[j].y,Components[i].RestrictSphereCenter[j].z);
              fprintf(FilePtr,"\t\tsphere radius :                   %lf\n",Components[i].RestrictSphereRadius[j]);
            }
          }
          fprintf(FilePtr,"\n");
        }
        else
          fprintf(FilePtr,"\tMoves are restricted to spheres: No\n\n");
      }
      else
        fprintf(FilePtr,"\tMoves are restricted: No\n");


      switch(Components[i].Biased)
      {
        case UMBRELLA:
          fprintf(FilePtr,"\tBiased sampling using: Umbrella sampling\n");
          fprintf(FilePtr,"\t  Umbrealla factor: %lf\n",
              (double)Components[i].UmbrellaFactor);
          break;
        case RUIZ_MONTERO:
          fprintf(FilePtr,"\tBiased sampling using: Ruiz-Montero\n");
          fprintf(FilePtr,"\t  Ruiz-Montero factor: %lf\n",
              (double)Components[i].RuizMonteroFactor);
          break;
        case NO_BIASING:
          fprintf(FilePtr,"\tNo biased sampling used for this component\n");
          break;
      }
      if(Components[i].Biased!=NO_BIASING)
      {
        switch(Components[i].BiasingDirection)
        {
          case A_MAPPING:
            fprintf(FilePtr,"\t\tBiasing-direction: A\n");
            break;
          case B_MAPPING:
            fprintf(FilePtr,"\t\tBiasing-direction: B\n");
            break;
          case C_MAPPING:
            fprintf(FilePtr,"\t\tBiasing-direction: C\n");
            break;
          case MAP_AB_DIAGONAL:
            fprintf(FilePtr,"\t\tBiasing-direction: AB diagonal\n");
            break;
          case MAP_AC_DIAGONAL:
            fprintf(FilePtr,"\t\tBiasing-direction: AC diagonal\n");
            break;
          case MAP_BC_DIAGONAL:
            fprintf(FilePtr,"\t\tBiasing-direction: BC diagonal\n");
            break;
        }
      }
      fprintf(FilePtr,"\n");

      // printing Bond-data
      if(Components[i].NumberOfBonds>0)
      {
        fprintf(FilePtr,"\tNumber of Bonds: %d\n",Components[i].NumberOfBonds);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfBonds;j++)
        {
           fprintf(FilePtr,"\tBond interaction %d: A=%d B=%d Type:%s\n",
             j,
             Components[i].Bonds[j].A,
             Components[i].Bonds[j].B,
             BondTypes[Components[i].BondType[j]].Name);
           if(BondTypes[Components[i].BondType[j]].nr_args>0)
           {
             switch(Components[i].BondType[j])
             {
               case HARMONIC_BOND:
                 // 0.5*p0*SQR(r-p1);
                 // ===============================================
                 // p_0/k_B [K/A^2]   force constant
                 // p_1     [A]       reference bond distance
                 fprintf(FilePtr,"\t\tHARMONIC_BOND: p_0/k_B=%-10.6f [K/A^2], p_1=%-10.6f [A]\n",
                   (double)(Components[i].BondArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BondArguments[j][1]));
                 break;
               case MORSE_BOND:
                 // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
                 // ===============================================
                 // p_0/k_B [K]       force constant
                 // p_1     [A^-1]    parameter
                 // p_2     [A]       reference bond distance
                 fprintf(FilePtr,"\t\tMORSE_BOND: p_0/k_B=%-10.6f [K], p_1=%-10.6f [A^-1], p_2=%-10.6f [A]\n",
                    (double)(Components[i].BondArguments[j][0]*ENERGY_TO_KELVIN),
                    (double)(Components[i].BondArguments[j][1]),
                    (double)(Components[i].BondArguments[j][2]));
                 break;
               case LJ_12_6_BOND:
                 // A/r_ij^12-B/r_ij^6
                 // ===============================================
                 // p_0/k_B [K A^12]
                 // p_1/k_B [K A^6]
                 fprintf(FilePtr,"\t\tLJ_12_6_BOND: p_0/k_B=%-10.6f [K A^12], p_1/k_B=%-10.6f [K A^6]\n",
                   (double)(Components[i].BondArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BondArguments[j][1]*ENERGY_TO_KELVIN));
                 break;
               case LENNARD_JONES_BOND:
                 // 4*p_0*((p_1/r)^12-(p_1/r)^6)
                 // ===============================================
                 // p_0/k_B [K]
                 // p_1     [A]
                 fprintf(FilePtr,"\t\tLENNARD_JONES_BOND: p_0/k_B=%-10.6f [K], p_1=%-10.6f [A]\n",
                   (double)(Components[i].BondArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BondArguments[j][1]));
                 break;
               case BUCKINGHAM_BOND:
                 // p_0*exp(-p_1 r)-p_2/r^6
                 // ===============================================
                 // p_0/k_B [K]
                 // p_1     [A^-1]
                 // p_2/k_B [K A^6]
                 fprintf(FilePtr,"\t\tBUCKINGHAM_BOND: p_0/k_B=%-10.6f [K], p_1=%-10.6f [A^-1], p_2/k_B=%-10.6f [K A^6]\n",
                    (double)(Components[i].BondArguments[j][0]*ENERGY_TO_KELVIN),
                    (double)(Components[i].BondArguments[j][1]),
                    (double)(Components[i].BondArguments[j][2]*ENERGY_TO_KELVIN));
                 break;
               case RESTRAINED_HARMONIC_BOND:
                 // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
                 // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
                 // ===============================================
                 // p_0/k_B [K/A^2]
                 // p_1     [A]
                 // p_2     [A]
                 fprintf(FilePtr,"\t\tRESTRAINED_HARMONIC_BOND: p_0/k_B=%10.6f [K/A^2], p_1=%10.6f [A], p_2=%10.6f [A],\n",
                    Components[i].BondArguments[j][0]*ENERGY_TO_KELVIN,
                    Components[i].BondArguments[j][1],
                    Components[i].BondArguments[j][2]);
                 break;
               case QUARTIC_BOND:
                 // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
                 // ===========================================================
                 // p_0/k_B [K/A^2]
                 // p_1     [A]
                 // p_2/k_B [K/A^3]
                 // p_3/k_B [K/A^4]
                 fprintf(FilePtr,"\t\tQUARTIC_BOND: p_0/k_B=%-10.6f [K/A^2], p_1=%-10.6f [A], p_2/k_B=%-10.6f [K/A^3], p_3/k_B=%-10.6f [K/A^4]\n",
                   (double)(Components[i].BondArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BondArguments[j][1]),
                   (double)(Components[i].BondArguments[j][2]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BondArguments[j][3]*ENERGY_TO_KELVIN));
                 break;
               case CFF_QUARTIC_BOND:
                 // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
                 // ===============================================
                 // p_0/k_B [K/A^2]
                 // p_1     [A]
                 // p_2/k_B [K/A^3]
                 // p_3/k_B [K/A^4]
                 fprintf(FilePtr,"\t\tCFF_QUARTIC_BOND: p_0/k_B=%-10.6f [K/A^2], p_1=%-10.6f [A], p_2/k_B=%-10.6f [K/A^3], p_3/k_B=%-10.6f [K/A^4]\n",
                   (double)(Components[i].BondArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BondArguments[j][1]),
                   (double)(Components[i].BondArguments[j][2]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BondArguments[j][3]*ENERGY_TO_KELVIN));
                 break;
               case MM3_BOND:
                 // (1/2)p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
                 // =================================================================
                 // p_0     [mdyne/A molecule]
                 // p_1     [A]
                 fprintf(FilePtr,"\t\tMM3_BOND: p_0=%-10.6f [mdyne/A molecule], p_1=%-10.6f [A]\n",
                   (double)(Components[i].BondArguments[j][0]/(71.94*KCAL_PER_MOL_TO_ENERGY)),
                   (double)(Components[i].BondArguments[j][1]));
                 break;
               case RIGID_BOND:
               case FIXED_BOND:
                 fprintf(FilePtr,"\t\tr_0=%-18.10f [A]\n",
                   (double)(Components[i].BondArguments[j][0]));
                 break;
             }
           }
        }
        fprintf(FilePtr,"\n");
      }

      // printing Bond-dipole data
      if(Components[i].NumberOfBondDipoles>0)
      {
        fprintf(FilePtr,"\tNumber of Bond-dipoles: %d\n",Components[i].NumberOfBondDipoles);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfBondDipoles;j++)
        {
           fprintf(FilePtr,"\tBond-dipole %d: A=%d B=%d Magnitude: %g [D]\n",
             j,
             Components[i].BondDipoles[j].A,
             Components[i].BondDipoles[j].B,
             (double)(DEBYE_CONVERSION_FACTOR*Components[i].BondDipoleMagnitude[j]));
        }
        fprintf(FilePtr,"\n");
      }

      // printing UreyBradley-data
      if(Components[i].NumberOfUreyBradleys>0)
      {
        fprintf(FilePtr,"\tNumber of UreyBradleys: %d\n",Components[i].NumberOfUreyBradleys);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfUreyBradleys;j++)
        {
           fprintf(FilePtr,"\tUreyBradley interaction %d: A=%d B=%d C=%d Type:%s\n",
             j,
             Components[i].UreyBradleys[j].A,
             Components[i].UreyBradleys[j].B,
             Components[i].UreyBradleys[j].C,
             UreyBradleyTypes[Components[i].UreyBradleyType[j]].Name);
           if(UreyBradleyTypes[Components[i].UreyBradleyType[j]].nr_args>0)
           {
             switch(Components[i].UreyBradleyType[j])
             {
               case HARMONIC_UREYBRADLEY:
                 // 0.5*p0*SQR(r-p1);
                 // ===============================================
                 // p_0/k_B [K/A^2]   force constant
                 // p_1     [A]       reference bond distance
                 fprintf(FilePtr,"\t\tHARMONIC_UREYBRADLEY: p_0/k_B=%-10.6f [K/A^2], p_1=%-10.6f [A]\n",
                   (double)(Components[i].UreyBradleyArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].UreyBradleyArguments[j][1]));
                 break;
               case MORSE_UREYBRADLEY:
                 // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
                 // ===============================================
                 // p_0/k_B [K]       force constant
                 // p_1     [A^-1]    parameter
                 // p_2     [A]       reference bond distance
                 fprintf(FilePtr,"\t\tMORSE_UREYBRADLEY: p_0/k_B=%-10.6f [K], p_1=%-10.6f [A^-1], p_2=%-10.6f [A]\n",
                    (double)(Components[i].UreyBradleyArguments[j][0]*ENERGY_TO_KELVIN),
                    (double)(Components[i].UreyBradleyArguments[j][1]),
                    (double)(Components[i].UreyBradleyArguments[j][2]));
                 break;
               case LJ_12_6_UREYBRADLEY:
                 // A/r_ij^12-B/r_ij^6
                 // ===============================================
                 // p_0/k_B [K A^12]
                 // p_1/k_B [K A^6]
                 fprintf(FilePtr,"\t\tLJ_12_6_UREYBRADLEY: p_0/k_B=%-10.6f [K A^12], p_1/k_B=%-10.6f [K A^6]\n",
                   (double)(Components[i].UreyBradleyArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].UreyBradleyArguments[j][1]*ENERGY_TO_KELVIN));
                 break;
               case LENNARD_JONES_UREYBRADLEY:
                 // 4*p_0*((p_1/r)^12-(p_1/r)^6)
                 // ===============================================
                 // p_0/k_B [K]
                 // p_1     [A]
                 fprintf(FilePtr,"\t\tLENNARD_JONES_UREYBRADLEY: p_0/k_B=%-10.6f [K], p_1=%-10.6f [A]\n",
                   (double)(Components[i].UreyBradleyArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].UreyBradleyArguments[j][1]));
                 break;
               case BUCKINGHAM_UREYBRADLEY:
                 // p_0*exp(-p_1 r)-p_2/r^6
                 // ===============================================
                 // p_0/k_B [K]
                 // p_1     [A^-1]
                 // p_2/k_B [K A^6]
                 fprintf(FilePtr,"\t\tBUCKINGHAM_UREYBRADLEY: p_0/k_B=%-10.6f [K], p_1=%-10.6f [A^-1], p_2/k_B=%-10.6f [K A^6]\n",
                    (double)(Components[i].UreyBradleyArguments[j][0]*ENERGY_TO_KELVIN),
                    (double)(Components[i].UreyBradleyArguments[j][1]),
                    (double)(Components[i].UreyBradleyArguments[j][2]*ENERGY_TO_KELVIN));
                 break;
               case RESTRAINED_HARMONIC_UREYBRADLEY:
                 // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
                 // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
                 // ===============================================
                 // p_0/k_B [K/A^2]
                 // p_1     [A]
                 // p_2     [A]
                 fprintf(FilePtr,"\t\tRESTRAINED_HARMONIC_UREYBRADLEY: p_0/k_B=%10.6f [K/A^2], p_1=%10.6f [A], p_2=%10.6f [A],\n",
                    Components[i].UreyBradleyArguments[j][0]*ENERGY_TO_KELVIN,
                    Components[i].UreyBradleyArguments[j][1],
                    Components[i].UreyBradleyArguments[j][2]);
                 break;
               case QUARTIC_UREYBRADLEY:
                 // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
                 // ===========================================================
                 // p_0/k_B [K/A^2]
                 // p_1     [A]
                 // p_2/k_B [K/A^3]
                 // p_3/k_B [K/A^4]
                 fprintf(FilePtr,"\t\tQUARTIC_UREYBRADLEY: p_0/k_B=%-10.6f [K/A^2], p_1=%-10.6f [A], p_2/k_B=%-10.6f [K/A^3], p_3/k_B=%-10.6f [K/A^4]\n",
                   (double)(Components[i].UreyBradleyArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].UreyBradleyArguments[j][1]),
                   (double)(Components[i].UreyBradleyArguments[j][2]*ENERGY_TO_KELVIN),
                   (double)(Components[i].UreyBradleyArguments[j][3]*ENERGY_TO_KELVIN));
                 break;
               case CFF_QUARTIC_UREYBRADLEY:
                 // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
                 // ===============================================
                 // p_0/k_B [K/A^2]
                 // p_1     [A]
                 // p_2/k_B [K/A^3]
                 // p_3/k_B [K/A^4]
                 fprintf(FilePtr,"\t\tCFF_QUARTIC_UREYBRADLEY: p_0/k_B=%-10.6f [K/A^2], p_1=%-10.6f [A], p_2/k_B=%-10.6f [K/A^3], p_3/k_B=%-10.6f [K/A^4]\n",
                   (double)(Components[i].UreyBradleyArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].UreyBradleyArguments[j][1]),
                   (double)(Components[i].UreyBradleyArguments[j][2]*ENERGY_TO_KELVIN),
                   (double)(Components[i].UreyBradleyArguments[j][3]*ENERGY_TO_KELVIN));
                 break;
               case MM3_UREYBRADLEY:
                 // (1/2)p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
                 // =================================================================
                 // p_0     [mdyne/A molecule]
                 // p_1     [A]
                 fprintf(FilePtr,"\t\tMM3_UREYBRADLEY: p_0=%-10.6f [mdyne/A molecule], p_1=%-10.6f [A]\n",
                   (double)(Components[i].UreyBradleyArguments[j][0]/(71.94*KCAL_PER_MOL_TO_ENERGY)),
                   (double)(Components[i].UreyBradleyArguments[j][1]));
                 break;
               case RIGID_UREYBRADLEY:
               case FIXED_UREYBRADLEY:
                 fprintf(FilePtr,"\t\tr_0=%-18.10f [A]\n",
                   (double)(Components[i].BondArguments[j][0]));
                 break;
             }
           }
        }
        fprintf(FilePtr,"\n");
      }

      // printing Bend-data
      if(Components[i].NumberOfBends>0)
      {
        fprintf(FilePtr,"\tNumber of Bends: %d\n",Components[i].NumberOfBends);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfBends;j++)
        {
           fprintf(FilePtr,"\tBend interaction %d: A=%d B=%d C=%d Type:%s\n",
             j,
             Components[i].Bends[j].A,
             Components[i].Bends[j].B,
             Components[i].Bends[j].C,
             BendTypes[Components[i].BendType[j]].Name);
           if(BendTypes[Components[i].BendType[j]].nr_args>0)
           {
             switch(Components[i].BendType[j])
             {
               case HARMONIC_BEND:
                 // (1/2)p_0*(theta-p_1)^2
                 // ===============================================
                 // p_0/k_B [K/rad^2]
                 // p_1     [degrees]
                 fprintf(FilePtr,"\t\tHARMONIC_BEND: p_0/k_B=%-10.6f [K/rad^2], p_1=%-10.6f [degrees]\n",
                   (double)(Components[i].BendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BendArguments[j][1]*RAD2DEG));
                 break;
               case QUARTIC_BEND:
                 // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
                 // ======================================================================
                 // p_0/k_B [K/rad^2]
                 // p_1     [degrees]
                 // p_2/k_B [K/rad^3]
                 // p_3/k_B [K/rad^4]
                 fprintf(FilePtr,"\t\tQUARTIC_BEND:  p_0/k_B=%-10.6f [K/rad^2], p_1=%-10.6f [degrees], p_2/k_B=%-10.6f [K/rad^3], p_3/k_B=%-10.6f [K/rad^4]\n",
                   (double)(Components[i].BendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BendArguments[j][1]*RAD2DEG),
                   (double)(Components[i].BendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BendArguments[j][0]*ENERGY_TO_KELVIN));
                 break;
               case CFF_QUARTIC_BEND:
                 // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
                 // =====================================================
                 // p_0/k_B [K/rad^2]
                 // p_1     [degrees]
                 // p_2/k_B [K/rad^3]
                 // p_3/k_B [K/rad^4]
                 fprintf(FilePtr,"\t\tCFF_QUARTIC_BEND: p_0/k_B=%-10.6f [K/rad^2], p_1=%-10.6f [degrees], p_2/k_B=%-10.6f [K/rad^3], p_3/k_B=%-10.6f [K/rad^4]\n",
                   (double)(Components[i].BendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BendArguments[j][1]*RAD2DEG),
                   (double)(Components[i].BendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BendArguments[j][0]*ENERGY_TO_KELVIN));
                 break;
               case HARMONIC_COSINE_BEND:
                 // (1/2)*p_0*(cos(theta)-cos(p_1))^2
                 // ===============================================
                 // p_0/k_B [K]
                 // p_1     [degrees]
                 fprintf(FilePtr,"\t\tHARMONIC_COSINE_BEND: p_0/k_B=%-10.6f [K], p_1=%-10.6f [degrees]\n",
                   (double)(Components[i].BendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(acos(Components[i].BendArguments[j][1])*RAD2DEG));
                 break;
               case COSINE_BEND:
                 // p_0*(1+cos(p_1*theta-p_2))
                 // ===============================================
                 // p_0/k_B [K]
                 // p_1     [-]
                 // p_2     [degrees]
                 fprintf(FilePtr,"\t\tCOSINE_BEND: p_0/k_B=%-10.6f [K], p_1=%-10.6f [-], p_2=%-10.6f [degrees]\n",
                   (double)(Components[i].BendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BendArguments[j][1]),
                   (double)(Components[i].BendArguments[j][2]*RAD2DEG));
                 break;
               case TAFIPOLSKY_BEND:
                 // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
                 // ===============================================
                 // p_0/k_B [K]
                 fprintf(FilePtr,"\t\tTAFIPOLSKY_BEND: p_0/k_B=%-10.6f [K]\n",
                   (double)(Components[i].BendArguments[j][0]*ENERGY_TO_KELVIN));
                 break;
               case MM3_BEND:
                 // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
                 // =================================================================================================
                 // p_0/k_B [mdyne A/rad^2]
                 // p_1     [degrees]
                 fprintf(FilePtr,"\t\tMM3_BEND: p_0/k_B=%-10.6f [mdyne A/rad^2], p_1=%-10.6f [degrees]\n",
                   (double)(Components[i].BendArguments[j][0]/(0.021914*KCAL_PER_MOL_TO_ENERGY)),
                   (double)(Components[i].BendArguments[j][2]*RAD2DEG));
                 break;
               case FIXED_BEND:
                 // theta_0
                 fprintf(FilePtr,"\t\ttheta_0=%-18.10f [degrees]\n",
                   (double)(Components[i].BendArguments[j][0]*RAD2DEG));
                 break;
               default:
                 fprintf(stderr, "Unknown bend-potential\n");
                 break;
             }
          }
        }
        fprintf(FilePtr,"\n");
      }

      // printing Inversion Bend-data
      if(Components[i].NumberOfInversionBends>0)
      {
        fprintf(FilePtr,"\tNumber of Inversion Bends: %d\n",Components[i].NumberOfInversionBends);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfInversionBends;j++)
        {
           fprintf(FilePtr,"\tInversion Bend interaction %d: A=%d B=%d C=%d D=%d Type:%s\n",
             j,
             Components[i].InversionBends[j].A,
             Components[i].InversionBends[j].B,
             Components[i].InversionBends[j].C,
             Components[i].InversionBends[j].D,
             InversionBendTypes[Components[i].InversionBendType[j]].Name);
           if(InversionBendTypes[Components[i].InversionBendType[j]].nr_args>0)
           {
             switch(Components[i].InversionBendType[j])
             {
               case HARMONIC_INVERSION:
                 // (1/2)*p_0*(chi-p_1)^2
                 // ===============================================
                 // p_0/k_B [K]
                 // p_1     [degrees]
                 fprintf(FilePtr,"\t\tHARMONIC_INVERSION: p_0/k_B=%-10.6f [K], p_1=%-10.6f [degrees]\n",
                    (double)(Components[i].InversionBendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].InversionBendArguments[j][1]*RAD2DEG));
                 break;
               case HARMONIC_INVERSION2:
                 // (1/2)*p_0*(chi-p_1)^2
                 // ===============================================
                 // p_0/k_B [K]
                 // p_1     [degrees]
                 fprintf(FilePtr,"\t\tHARMONIC_INVERSION2: p_0/k_B=%-10.6f [K], p_1=%-10.6f [degrees]\n",
                    (double)(Components[i].InversionBendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(Components[i].InversionBendArguments[j][1]*RAD2DEG));
                 break;
               case HARMONIC_COSINE_INVERSION:
                 // (1/2)*p_0*(cos(phi)-cos(p_1))^2
                 // ===============================================
                 // p_0/k_B [K]
                 // p_1     [degrees]
                 fprintf(FilePtr,"\t\tHARMONIC_COSINE_INVERSION: p_0/k_B=%-10.6f [K], p_1=%-10.6f [degrees]\n",
                   (double)(Components[i].InversionBendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(acos(Components[i].InversionBendArguments[j][1])*RAD2DEG));
                 break;
               case HARMONIC_COSINE_INVERSION2:
                 // (1/2)*p_0*(cos(phi)-cos(p_1))^2
                 // ===============================================
                 // p_0/k_B [K]
                 // p_1     [degrees]
                 fprintf(FilePtr,"\t\tHARMONIC_COSINE_INVERSION2: p_0/k_B=%-10.6f [K], p_1=%-10.6f [degrees]\n",
                   (double)(Components[i].InversionBendArguments[j][0]*ENERGY_TO_KELVIN),
                   (double)(acos(Components[i].InversionBendArguments[j][1])*RAD2DEG));
                 break;
               case PLANAR_INVERSION:
                 // (1/2)*p_0*(1-cos(phi))
                 // ===============================================
                 // p_0/k_B [K]
                 fprintf(FilePtr,"\t\tPLANAR_INVERSION: p_0/k_B=%-10.6f [K]\n",
                   (double)(Components[i].InversionBendArguments[j][0]*ENERGY_TO_KELVIN));
                 break;
               case MM3_INVERSION:
                 // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
                 // =================================================================================================
                 // p_0/k_B [mdyne A/rad^2]
                 // p_1     [degrees]
                 fprintf(FilePtr,"\t\tMM3_INVERSION: p_0/k_B=%-10.6f [mdyne A/rad^2], p_1=%-10.6f [degrees]\n",
                   (double)(Components[i].BendArguments[j][0]/(0.021914*KCAL_PER_MOL_TO_ENERGY)),
                   (double)(Components[i].BendArguments[j][2]*RAD2DEG));
                 break;
               case FIXED_INVERSION_BEND:
                 fprintf(FilePtr,"\t\tchi_0=%-18.10f [degrees]\n",
                   (double)(Components[i].InversionBendArguments[j][0]*RAD2DEG));
                 break;
               default:
                 fprintf(stderr, "Unknown inversion bend-potential\n");
                 break;
             }
          }
        }
        fprintf(FilePtr,"\n");
      }

      // printing Torsion-data
      if(Components[i].NumberOfTorsions>0)
      {
        fprintf(FilePtr,"\tNumber of Torsions: %d\n",Components[i].NumberOfTorsions);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfTorsions;j++)
        {
          fprintf(FilePtr,"\tTorsions interaction %d: A=%d B=%d C=%d D=%d Type:%s\n",
             j,
             Components[i].Torsions[j].A,
             Components[i].Torsions[j].B,
             Components[i].Torsions[j].C,
             Components[i].Torsions[j].D,
             TorsionTypes[Components[i].TorsionType[j]].Name);
          if(TorsionTypes[Components[i].TorsionType[j]].nr_args>0)
          {
            switch(Components[i].TorsionType[j])
            {
              case HARMONIC_DIHEDRAL:
                // (1/2)*p_0*(phi-p_1)^2
                // ===============================================
                // p_0/k_B [K/rad^2]
                // p_1     [degrees]
                fprintf(FilePtr,"\t\tHARMONIC_DIHEDRAL: p_0/k_B=%-10.6f [K/rad^2], p_2=%-10.6f [degrees]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]*RAD2DEG));
                break;
              case HARMONIC_COSINE_DIHEDRAL:
                // (1/2)*p_0*(cos(phi)-cos(p_1))^2
                // ===============================================
                // p_0/k_B [K]
                // p_1     [degrees]
                fprintf(FilePtr,"\t\tHARMONIC_COSINE_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1=%-10.6f [degrees]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(acos(Components[i].TorsionArguments[j][1])*RAD2DEG));
                break;
              case THREE_COSINE_DIHEDRAL:
                // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
                // ========================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                fprintf(FilePtr,"\t\tTHREE_COSINE_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
             case MM3_DIHEDRAL:
                // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
                // ========================================================================
                // p_0     [kcal/mol]
                // p_1     [kcal/mol]
                // p_2     [kcal/mol]
                fprintf(FilePtr,"\t\tMM3_DIHEDRAL: p_0=%-10.6f [kcal/mol], p_1=%-10.6f [kcal/mol], p_2=%-10.6f [kcal/mol]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KCAL_PER_MOL),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KCAL_PER_MOL),
                  (double)(Components[i].TorsionArguments[j][2]*ENERGY_TO_KCAL_PER_MOL));
                break;
             case CVFF_BLOCKED_DIHEDRAL:
                // 
                // ========================================================================
                // p_0     [rad]
                // p_1     [K]
                // p_2     [-]
                // p_3     [rad]
                // p_4     [rad]
                fprintf(FilePtr,"\t\tCVFF_BLOCKED_DIHEDRAL: p_0=%-10.6f [rad], p_1=%-10.6f [K], p_2=%-10.6f [-], p_3=%-10.6f [rad], p_4=%-10.6f [rad]\n",
                  (double)(Components[i].TorsionArguments[j][0]),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][2]),
                  (double)(Components[i].TorsionArguments[j][3]),
                  (double)(Components[i].TorsionArguments[j][4]));
                break;
              case CFF_DIHEDRAL:
                // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
                // ======================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                fprintf(FilePtr,"\t\tCFF_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
              case CFF_DIHEDRAL2:
                // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
                // ======================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                fprintf(FilePtr,"\t\tCFF_DIHEDRAL2: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
              case SIX_COSINE_DIHEDRAL:
                // Prod_i=0^5 p_i*cos(phi)^i
                // =========================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                // p_4/k_B [K]
                // p_5/k_B [K]
                fprintf(FilePtr,"\t\tSIX_COSINE_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K], p_4/k_B=%-10.6f [K], p_5/k_B=%-10.6f [K]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][3]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][4]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][5]*ENERGY_TO_KELVIN));
                break;
              case TRAPPE_DIHEDRAL:
                // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
                // =============================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                fprintf(FilePtr,"\t\tTRAPPE_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][3]*ENERGY_TO_KELVIN));
                break;
              case TRAPPE_DIHEDRAL_EXTENDED:
                // p_0[0]+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
                // ================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                // p_4/k_B [K]
                fprintf(FilePtr,"\t\tTRAPPE_DIHEDRAL_EXTENDED: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K], p_4/k_B=%-10.6f [K]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][3]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][4]*ENERGY_TO_KELVIN));
                break;
              case CVFF_DIHEDRAL:
                // p_0*(1+cos(p_1*phi-p_2))
                // ========================
                // p_0/k_B [K]
                // p_1     [-]
                // p_2     [degrees]
                fprintf(FilePtr,"\t\tCVFF_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1=%-10.6f [-], p_2=%-10.6f [degrees]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]),
                  (double)(Components[i].TorsionArguments[j][2]*RAD2DEG));
                break;
              case OPLS_DIHEDRAL:
                // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
                // =================================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                fprintf(FilePtr,"\t\tOPLS_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][3]*ENERGY_TO_KELVIN));
                break;
             case FOURIER_SERIES_DIHEDRAL:
                // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
                // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
                // =======================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                // p_4/k_B [K]
                // p_5/k_B [K]
                fprintf(FilePtr,"\t\tFOURIER_SERIES_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K], p_4/k_B=%-10.6f [K], p_5/k_B=%-10.6f [K]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][3]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][4]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][5]*ENERGY_TO_KELVIN));
                break;
             case FOURIER_SERIES_DIHEDRAL2:
                // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
                // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
                // =======================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                // p_4/k_B [K]
                // p_5/k_B [K]
                fprintf(FilePtr,"\t\tFOURIER_SERIES_DIHEDRAL2: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K], p_4/k_B=%-10.6f [K], p_5/k_B=%-10.6f [K]\n",
                  (double)(Components[i].TorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][3]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][4]*ENERGY_TO_KELVIN),
                  (double)(Components[i].TorsionArguments[j][5]*ENERGY_TO_KELVIN));
                break;
             case FIXED_DIHEDRAL:
                 fprintf(FilePtr,"\t\tphi_0=%-18.10f [degrees]\n",
                   (double)(Components[i].TorsionArguments[j][0]*RAD2DEG));
                break;
              default:
                fprintf(stderr, "Unknown torsion-potential\n");
                break;
            }
          }
        }
        fprintf(FilePtr,"\n");
      }

      // printing improper Torsion-data
      if(Components[i].NumberOfImproperTorsions>0)
      {
        fprintf(FilePtr,"\tNumber of Improper Torsions: %d\n",Components[i].NumberOfImproperTorsions);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfImproperTorsions;j++)
        {
          fprintf(FilePtr,"\tImproper Torsions interaction %d: A=%d B=%d C=%d D=%d Type:%s\n",
             j,
             Components[i].ImproperTorsions[j].A,
             Components[i].ImproperTorsions[j].B,
             Components[i].ImproperTorsions[j].C,
             Components[i].ImproperTorsions[j].D,
             ImproperTorsionTypes[Components[i].ImproperTorsionType[j]].Name);
          if(ImproperTorsionTypes[Components[i].ImproperTorsionType[j]].nr_args>0)
          {
            switch(Components[i].ImproperTorsionType[j])
            {
              case HARMONIC_IMPROPER_DIHEDRAL:
                // (1/2)*p_0*(phi-p_1)^2
                // ===============================================
                // p_0/k_B [K/rad^2]
                // p_1     [degrees]
                fprintf(FilePtr,"\t\tHARMONIC_IMPROPER_DIHEDRAL: p_0/k_B=%-10.6f [K/rad^2], p_2=%-10.6f [degrees]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*RAD2DEG));
                break;
              case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
                // (1/2)*p_0*(cos(phi)-cos(p_1))^2
                // ===============================================
                // p_0/k_B [K]
                // p_1     [degrees]
                fprintf(FilePtr,"\t\tHARMONIC_COSINE_IMPROPER_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1=%-10.6f [degrees]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(acos(Components[i].ImproperTorsionArguments[j][1])*RAD2DEG));
                break;
              case THREE_COSINE_IMPROPER_DIHEDRAL:
                // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
                // ========================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                fprintf(FilePtr,"\t\tTHREE_COSINE_IMPROPER_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
             case MM3_IMPROPER_DIHEDRAL:
                // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
                // ========================================================================
                // p_0     [kcal/mol]
                // p_1     [kcal/mol]
                // p_2     [kcal/mol]
                fprintf(FilePtr,"\t\tMM3_IMPROPER_DIHEDRAL: p_0=%-10.6f [kcal/mol], p_1=%-10.6f [kcal/mol], p_2=%-10.6f [kcal/mol]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KCAL_PER_MOL),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*ENERGY_TO_KCAL_PER_MOL),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*ENERGY_TO_KCAL_PER_MOL));
                break;
              case CFF_IMPROPER_DIHEDRAL:
                // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
                // ======================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                fprintf(FilePtr,"\t\tCFF_IMPROPER_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
              case CFF_IMPROPER_DIHEDRAL2:
                // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
                // ======================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                fprintf(FilePtr,"\t\tCFF_IMPROPER_DIHEDRAL2: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
              case SIX_COSINE_IMPROPER_DIHEDRAL:
                // Prod_i=0^5 p_i*cos(phi)^i
                // =========================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                // p_4/k_B [K]
                // p_5/k_B [K]
                fprintf(FilePtr,"\t\tSIX_COSINE_IMPROPER_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K], p_4/k_B=%-10.6f [K], p_5/k_B=%-10.6f [K]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][3]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][4]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][5]*ENERGY_TO_KELVIN));
                break;
              case TRAPPE_IMPROPER_DIHEDRAL:
                // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
                // =============================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                fprintf(FilePtr,"\t\tTRAPPE_IMPROPER_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][3]*ENERGY_TO_KELVIN));
                break;
              case TRAPPE_IMPROPER_DIHEDRAL_EXTENDED:
                // p_0[0]+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
                // ================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                // p_4/k_B [K]
                fprintf(FilePtr,"\t\tTRAPPE_IMPROPER_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K], p_4/k_B=%-10.6f [K]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][3]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][4]*ENERGY_TO_KELVIN));
                break;
             case CVFF_IMPROPER_DIHEDRAL:
                // p_0*(1+cos(p_1*phi-p_2))
                // ========================
                // p_0/k_B [K]
                // p_1     [-]
                // p_2     [degrees]
                fprintf(FilePtr,"\t\tCVFF_IMPROPER_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1=%-10.6f [-], p_2=%-10.6f [degrees]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*RAD2DEG));
                break;
              case OPLS_IMPROPER_DIHEDRAL:
                // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
                // =================================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                fprintf(FilePtr,"\t\tOPLS_IMPROPER_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][3]*ENERGY_TO_KELVIN));
                break;
             case FOURIER_SERIES_IMPROPER_DIHEDRAL:
                // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
                // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
                // =======================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                // p_4/k_B [K]
                // p_5/k_B [K]
                fprintf(FilePtr,"\t\tFOURIER_SERIES_IMPROPER_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K], p_4/k_B=%-10.6f [K], p_5/k_B=%-10.6f [K]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][3]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][4]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][5]*ENERGY_TO_KELVIN));
                break;
             case FOURIER_SERIES_IMPROPER_DIHEDRAL2:
                // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
                // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
                // =======================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                // p_3/k_B [K]
                // p_4/k_B [K]
                // p_5/k_B [K]
                fprintf(FilePtr,"\t\tFOURIER_SERIES_IMPROPER_DIHEDRAL2: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K], p_3/k_B=%-10.6f [K], p_4/k_B=%-10.6f [K], p_5/k_B=%-10.6f [K]\n",
                  (double)(Components[i].ImproperTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][2]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][3]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][4]*ENERGY_TO_KELVIN),
                  (double)(Components[i].ImproperTorsionArguments[j][5]*ENERGY_TO_KELVIN));
                break;
             case FIXED_IMPROPER_DIHEDRAL:
                 fprintf(FilePtr,"\t\tphi_0=%-18.10f [degrees]\n",
                   (double)(Components[i].ImproperTorsionArguments[j][0]*RAD2DEG));
                break;
              default:
                fprintf(stderr, "Unknown improper torsion-potential\n");
                break;
            }
          }
        }
        fprintf(FilePtr,"\n");
      }

      // printing out-of-plane distance data
      if(Components[i].NumberOfOutOfPlanes>0)
      {
        fprintf(FilePtr,"\tNumber of Out-of-plane distances: %d\n",Components[i].NumberOfOutOfPlanes);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfOutOfPlanes;j++)
        {
          fprintf(FilePtr,"\tOut-of-plane distance interaction %d: A=%d B=%d C=%d D=%d Type:%s\n",
             j,
             Components[i].OutOfPlanes[j].A,
             Components[i].OutOfPlanes[j].B,
             Components[i].OutOfPlanes[j].C,
             Components[i].OutOfPlanes[j].D,
             OutOfPlaneTypes[Components[i].OutOfPlaneType[j]].Name);
          if(OutOfPlaneTypes[Components[i].OutOfPlaneType[j]].nr_args>0)
          {
            switch(Components[i].OutOfPlaneType[j])
            {
              default:
                fprintf(stderr, "Unknown Out-of-plane-distance potential\n");
                break;
            }
          }
        }
        fprintf(FilePtr,"\n");
      }


      // printing Bond/Bond cross term-data
      if(Components[i].NumberOfBondBonds>0)
      {
        fprintf(FilePtr,"\tNumber of Bond/Bond cross terms: %d\n",Components[i].NumberOfBondBonds);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfBondBonds;j++)
        {
           fprintf(FilePtr,"\tBond/Bond interaction %d: A=%d B=%d C: %d Type:%s\n",
             j,
             Components[i].BondBonds[j].A,
             Components[i].BondBonds[j].B,
             Components[i].BondBonds[j].C,
             BondBondTypes[Components[i].BondBondType[j]].Name);

           if(BondBondTypes[Components[i].BondBondType[j]].nr_args>0)
           {
             switch(Components[i].BondBondType[j])
             {
               case CVFF_BOND_BOND_CROSS:
                // p_0*(rab-p_1)*(rbc-p_2)
                // =======================
                // p_0/k_B [K/A^2]
                // p_1     [A]
                // p_2     [A]
                fprintf(FilePtr,"\t\tCVFF_BOND_BOND_CROSS: p_0/k_B=%-10.6f [K/A^2], p_1=%-10.6f [A], p_2=%-10.6f [A]\n",
                  (double)(Components[i].BondBondArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BondBondArguments[j][1]),
                  (double)(Components[i].BondBondArguments[j][2]));
                 break;
              default:
                fprintf(stderr, "Unknown bond/bond-potential\n");
                break;
             }
           }
        }
        fprintf(FilePtr,"\n");
      }

      // printing Bond/Bend cross term-data
      if(Components[i].NumberOfBondBends>0)
      {
        fprintf(FilePtr,"\tNumber of Bond/Bend cross terms: %d\n",Components[i].NumberOfBondBends);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfBondBends;j++)
        {
           fprintf(FilePtr,"\tBond/Bend interaction %d: A=%d B=%d C: %d Type:%s\n",
             j,
             Components[i].BondBends[j].A,
             Components[i].BondBends[j].B,
             Components[i].BondBends[j].C,
             BondBendTypes[Components[i].BondBendType[j]].Name);
           if(BondBendTypes[Components[i].BondBendType[j]].nr_args>0)
           {
             switch(Components[i].BondBendType[j])
             {
               case CVFF_BOND_BEND_CROSS:
                 // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
                 // =========================================
                 // p_0     [degrees]
                 // p_1/k_B [K/A/rad]
                 // p_2     [A]
                 // p_3/k_B [K/A/rad]
                 // p_4     [A]
                 fprintf(FilePtr,"\t\tCVFF_BOND_BEND_CROSS: p_0=%-10.6f [degrees], p_1/k_B=%-10.6f [K/A/rad], p_2=%-10.6f [A], p_3/k_B=%-10.6f [K/A/rad], p_4=%-10.6f [A]\n",
                   (double)(Components[i].BondBendArguments[j][0]*RAD2DEG),
                   (double)(Components[i].BondBendArguments[j][1]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BondBendArguments[j][2]),
                   (double)(Components[i].BondBendArguments[j][3]*ENERGY_TO_KELVIN),
                   (double)(Components[i].BondBendArguments[j][4]));
                  break;
               case MM3_BOND_BEND_CROSS:
                 // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
                 // =====================================
                 // p_0     [mdyne/rad]
                 // p_1     [A]
                 // p_2     [A]
                 // p_3     [degrees]
                 fprintf(FilePtr,"\t\tMM3_BOND_BEND_CROSS: p_0=%-10.6f [mdyne/rad], p_1=%-10.6f [A], p_2=%-10.6f [A], p_3=%-10.6f [degrees]\n",
                   (double)(Components[i].BondBendArguments[j][0]/(2.51118*KCAL_PER_MOL_TO_ENERGY)),
                   (double)(Components[i].BondBendArguments[j][1]),
                   (double)(Components[i].BondBendArguments[j][2]),
                   (double)(Components[i].BondBendArguments[j][3]*RAD2DEG));
                 break;
               case TRUNCATED_HARMONIC:
                 // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
                 // ================================================================
                 // p_0/k_B [K/rad^2]
                 // p_1     [degrees]
                 // p_2     [A]
                fprintf(FilePtr,"\t\tTRUNCATED_HARMONIC: p_0/k_B=%-10.6f [K/rad^2], p_1=%-10.6f [degrees], p_2=%-10.6f [A]\n",
                  (double)(Components[i].BondBendArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BondBendArguments[j][1]*RAD2DEG),
                  (double)(Components[i].BondBendArguments[j][2]));
                 break;
               case SCREENED_HARMONIC:
                 // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
                 // ===============================================
                 // p_0/k_B [K/rad^2]
                 // p_1     [degrees]
                 // p_2     [A]
                 // p_3     [A]
                fprintf(FilePtr,"\t\tSCREENED_HARMONIC: p_0/k_B=%-10.6f [K/rad^2], p_1=%-10.6f [degrees], p_2=%-10.6f [A], p_3=%-10.6f [A]\n",
                  (double)(Components[i].BondBendArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BondBendArguments[j][1]*RAD2DEG),
                  (double)(Components[i].BondBendArguments[j][2]),
                  (double)(Components[i].BondBendArguments[j][3]));
                 break;
               case SCREENED_VESSAL:
                 // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
                 // ============================================================================
                 // p_0/k_B [K/rad^2]
                 // p_1     [degrees]
                 // p_2     [A]
                 // p_3     [A]
                fprintf(FilePtr,"\t\tSCREENED_VESSAL: p_0/k_B=%-10.6f [K/rad^2], p_1=%-10.6f [degrees], p_2=%-10.6f [A], p_3=%-10.6f [A]\n",
                  (double)(Components[i].BondBendArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BondBendArguments[j][1]*RAD2DEG),
                  (double)(Components[i].BondBendArguments[j][2]),
                  (double)(Components[i].BondBendArguments[j][3]));
                 break;
               case TRUNCATED_VESSAL:
                 // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
                 //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
                 // ============================================================================
                 // p_0/k_B [K/rad^(4+p_2)]
                 // p_1     [degrees]
                 // p_2     [-]
                 // p_3     [A]
                fprintf(FilePtr,"\t\tTRUNCATED_VESSAL: p_0/k_B=%-10.6f [K/rad^(4+p_2)], p_1=%-10.6f [degrees], p_2=%-10.6f [-], p_3=%-10.6f [A]\n",
                  (double)(Components[i].BondBendArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BondBendArguments[j][1]),
                  (double)(Components[i].BondBendArguments[j][2]),
                  (double)(Components[i].BondBendArguments[j][3]));
                 break;
               default:
                 fprintf(stderr, "Unknown stretch/bend-potential\n");
                 break;
             }
           }
        }
        fprintf(FilePtr,"\n");
      }


      // printing Bend/Bend cross term-data
      if(Components[i].NumberOfBendBends>0)
      {
        fprintf(FilePtr,"\tNumber of Bend/Bend cross terms: %d\n",Components[i].NumberOfBendBends);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfBendBends;j++)
        {
           fprintf(FilePtr,"\tBend/Bend interaction %d: A=%d B=%d C: %d Type:%s\n",
             j,
             Components[i].BendBends[j].A,
             Components[i].BendBends[j].B,
             Components[i].BendBends[j].C,
             BendBendTypes[Components[i].BendBendType[j]].Name);
           if(BendBendTypes[Components[i].BendBendType[j]].nr_args>0)
           {
             switch(Components[i].BendBendType[j])
             {
               case CVFF_BEND_BEND_CROSS:
                 // p_0*(Theta1-p_1)*(Theta2-p_2)
                 // ===================================
                 // p_0/k_B [K/rad^2)]
                 // p_1     [degrees]
                 // p_2     [degrees]
                fprintf(FilePtr,"\t\tCVFF_BEND_BEND_CROSS: p_0/k_B=%-18.10f [K/rad^2], p_1=%-10.6f [degrees], p_2=%-10.6f [degrees]\n",
                  (double)(Components[i].BendBendArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendBendArguments[j][1]*RAD2DEG),
                  (double)(Components[i].BendBendArguments[j][1]*RAD2DEG));
                 break;
               case MM3_BEND_BEND_CROSS:
                 // -p_0*(Theta1-p_1)*(Theta2-p_2)
                 // ===================================
                 // p_0     [mdyne A/rad^2]
                 // p_1     [degrees]
                 // p_2     [degrees]
                fprintf(FilePtr,"\t\tMM3_BEND_BEND_CROSS: p_0/k_B=%-10.6f [mdyne A/rad^2], p_1=%-10.6f [degrees], p_2=%-10.6f [degrees]\n",
                  (double)(Components[i].BendBendArguments[j][0]/(0.02191418*KCAL_PER_MOL_TO_ENERGY)),
                  (double)(Components[i].BendBendArguments[j][1]*RAD2DEG),
                  (double)(Components[i].BendBendArguments[j][1]*RAD2DEG));
                 break;
               default:
                 fprintf(stderr, "Unknown bend/bend-potential\n");
                 break;
              }
           }
        }
        fprintf(FilePtr,"\n");
      }

      // printing Bond/Torsion cross term-data
      if(Components[i].NumberOfBondTorsions>0)
      {
        fprintf(FilePtr,"\tNumber of Bond/Torsion cross terms: %d\n",Components[i].NumberOfBondTorsions);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfBondTorsions;j++)
        {
          fprintf(FilePtr,"\tBond/Torsion interaction %d: A=%d B=%d C: %d D:%d Type:%s\n",
             j,
             Components[i].BondTorsions[j].A,
             Components[i].BondTorsions[j].B,
             Components[i].BondTorsions[j].C,
             Components[i].BondTorsions[j].D,
             BondTorsionTypes[Components[i].BondTorsionType[j]].Name);
          if(BondTorsionTypes[Components[i].BondTorsionType[j]].nr_args>0)
          {
            switch(Components[i].BondTorsionType[j])
            {
              case MM3_BOND_TORSION_CROSS:
                // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
                // =====================================================================================
                // p_0     [kcal/A mole]
                // p_1     [kcal/A mole]
                // p_2     [kcal/A mole]
                // p_3     [A]
                fprintf(FilePtr,"\t\tMM3_BOND_TORSION_CROSS: p_0=%-10.6f [kcal/mole], p_1=%-10.6f [kcal/mole], p_2=%-10.6f [kcal/mole], p_3=%-10.6f [A]\n",
                  (double)(Components[i].BondTorsionArguments[j][0]*ENERGY_TO_KCAL_PER_MOL),
                  (double)(Components[i].BondTorsionArguments[j][1]*ENERGY_TO_KCAL_PER_MOL),
                  (double)(Components[i].BondTorsionArguments[j][2]*ENERGY_TO_KCAL_PER_MOL),
                  (double)(Components[i].BondTorsionArguments[j][3]));
                break;
                break;
              default:
                fprintf(stderr, "Unknown stretch/torsion-potential\n");
                break;
            }
          }
        }
        fprintf(FilePtr,"\n");
      }

      // printing Bend/Torsion cross term-data
      if(Components[i].NumberOfBendTorsions>0)
      {
        fprintf(FilePtr,"\tNumber of Bend/Torsion cross terms: %d\n",Components[i].NumberOfBendTorsions);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfBendTorsions;j++)
        {
          fprintf(FilePtr,"\tBend/Torsion interaction %d: A=%d B=%d C: %d D:%d Type:%s\n",
             j,
             Components[i].BendTorsions[j].A,
             Components[i].BendTorsions[j].B,
             Components[i].BendTorsions[j].C,
             Components[i].BendTorsions[j].D,
             BendTorsionTypes[Components[i].BendTorsionType[j]].Name);
          if(BendTorsionTypes[Components[i].BendTorsionType[j]].nr_args>0)
          {
            switch(Components[i].BendTorsionType[j])
            {
              case CFF_BEND_TORSION_CROSS:
                // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
                // =====================================================================================
                // p_0/k_B [K/rad^3]
                // p_1     [degrees]
                // p_2     [degrees]
                fprintf(FilePtr,"\t\tCFF_BEND_TORSION_CROSS: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [degrees], p_2/k_B=%-10.6f [degrees]\n",
                  (double)(Components[i].BendTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
              case CVFF_BEND_TORSION_CROSS:
                // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
                // =====================================================================================
                // p_0/k_B [K/rad^3]
                // p_1     [degrees]
                // p_2     [degrees]
                // 0.5*parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi
                fprintf(FilePtr,"\t\tCVFF_BEND_TORSION_CROSS: p_0/k_B=%-10.6f [K/rad^3], p_1=%-10.6f [degrees], p_2=%-10.6f [degrees]\n",
                  (double)(Components[i].BendTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][1]*RAD2DEG),
                  (double)(Components[i].BendTorsionArguments[j][2]*RAD2DEG));
                break;
              case SMOOTHED_DIHEDRAL:
                // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2))]*S(Theta2)
                // ======================================================================================
                // p_0/k_B [K/rad^2]
                // p_1     [-]
                // p_2     [degrees]
                fprintf(FilePtr,"\t\tSMOOTHED_DIHEDRAL: p_0/k_B=%-10.6f [K/rad^2], p_1=%-10.6f [-], p_2/k_B=%-10.6f [degrees]\n",
                  (double)(Components[i].BendTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][1]),
                  (double)(Components[i].BendTorsionArguments[j][2]*RAD2DEG));
                break;
              case SMOOTHED_THREE_COSINE_DIHEDRAL:
                // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
                // ======================================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                fprintf(FilePtr,"\t\tSMOOTHED_THREE_COSINE_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K]\n",
                  (double)(Components[i].BendTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
              case NICHOLAS_DIHEDRAL:
                // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
                // ======================================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                fprintf(FilePtr,"\t\tNICHOLAS_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K]\n",
                  (double)(Components[i].BendTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
              case SMOOTHED_CFF_DIHEDRAL:
                // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
                // ======================================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                fprintf(FilePtr,"\t\tSMOOTHED_CFF_DIHEDRAL: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K]\n",
                  (double)(Components[i].BendTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
              case SMOOTHED_CFF_DIHEDRAL2:
                // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
                // ======================================================================================
                // p_0/k_B [K]
                // p_1/k_B [K]
                // p_2/k_B [K]
                fprintf(FilePtr,"\t\tSMOOTHED_CFF_DIHEDRAL2: p_0/k_B=%-10.6f [K], p_1/k_B=%-10.6f [K], p_2/k_B=%-10.6f [K]\n",
                  (double)(Components[i].BendTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][1]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][2]*ENERGY_TO_KELVIN));
                break;
              case SMOOTHED_CFF_BEND_TORSION_CROSS:
                // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
                // ======================================================================================
                // p_0/k_B [K/rad^3]
                // p_1     [degrees]
                // p_2     [degrees]
                fprintf(FilePtr,"\t\tSMOOTHED_CFF_BEND_TORSION_CROSS: p_0/k_B=%-10.6f [K], p_1=%-10.6f [degrees], p_2=%-10.6f [degrees]\n",
                  (double)(Components[i].BendTorsionArguments[j][0]*ENERGY_TO_KELVIN),
                  (double)(Components[i].BendTorsionArguments[j][1]*RAD2DEG),
                  (double)(Components[i].BendTorsionArguments[j][2]*RAD2DEG));
                break;
              default:
                fprintf(stderr, "Unknown bend/torsion-potential\n");
                break;
            }
          }
        }
        fprintf(FilePtr,"\n");
      }


      if(Components[i].NumberOfIntraVDW>0)
      {
        fprintf(FilePtr,"\tnumber of intra VDW interactions: %d\n",Components[i].NumberOfIntraVDW);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfIntraVDW;j++)
          fprintf(FilePtr,"\tinteraction %d: A=%d B=%d  scaling factor: %g\n",
                j,
                Components[i].IntraVDW[j].A,
                Components[i].IntraVDW[j].B,
                Components[i].IntraVDWScaling[j]);
        fprintf(FilePtr,"\n");
      }

      if(Components[i].NumberOfIntraChargeCharge>0)
      {
        fprintf(FilePtr,"\tnumber of intra charge-charge Coulomb interactions: %d\n",Components[i].NumberOfIntraChargeCharge);
        fprintf(FilePtr,"\t-----------------------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfIntraChargeCharge;j++)
          fprintf(FilePtr,"\tinteraction %d: A=%d B=%d  scaling factor: %g\n",
                j,
                Components[i].IntraChargeCharge[j].A,
                Components[i].IntraChargeCharge[j].B,
                Components[i].IntraChargeChargeScaling[j]);
        fprintf(FilePtr,"\n");
      }

      if(Components[i].NumberOfIntraChargeBondDipole>0)
      {
        fprintf(FilePtr,"\tnumber of intra charge-bonddipole Coulomb interactions: %d\n",Components[i].NumberOfIntraChargeBondDipole);
        fprintf(FilePtr,"\t-------------------------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfIntraChargeBondDipole;j++)
        {
          A=Components[i].IntraChargeBondDipole[j].A;
          B=Components[i].IntraChargeBondDipole[j].B;
          fprintf(FilePtr,"\tinteraction %d: A=%d B=%d  [atoms: (%d) - (%d,%d)]\n",
                j,
                A,
                B,
                A,
                Components[i].BondDipoles[B].A,
                Components[i].BondDipoles[B].B);
        }
        fprintf(FilePtr,"\n");
      }

      if(Components[i].NumberOfIntraBondDipoleBondDipole>0)
      {
        fprintf(FilePtr,"\tnumber of intra bonddipole-bonddipole Coulomb interactions: %d\n",Components[i].NumberOfIntraBondDipoleBondDipole);
        fprintf(FilePtr,"\t-----------------------------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfIntraBondDipoleBondDipole;j++)
        {
          A=Components[i].IntraBondDipoleBondDipole[j].A;
          B=Components[i].IntraBondDipoleBondDipole[j].B;
          fprintf(FilePtr,"\tinteraction %d: A=%d B=%d  [atoms: (%d,%d) - (%d,%d)]\n",
                j,
                A,
                B,
                Components[i].BondDipoles[A].A,
                Components[i].BondDipoles[A].B,
                Components[i].BondDipoles[B].A,
                Components[i].BondDipoles[B].B);
        }
        fprintf(FilePtr,"\n");
      }

      if(Components[i].NumberOfConfigMoves>0)
      {
        fprintf(FilePtr,"\tnumber of cbmc-config changes: %d\n",
                Components[i].NumberOfConfigMoves);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfConfigMoves;j++)
        {
          fprintf(FilePtr,"\tnr fixed %d: ",Components[i].NumberOfUnchangedAtomsConfig[j]);
          for(k=0;k<Components[i].NumberOfUnchangedAtomsConfig[j];k++)
            fprintf(FilePtr,"%d ",Components[i].UnchangedAtomsConfig[j][k]);
          fprintf(FilePtr,"\n");
        }
      }
      fprintf(FilePtr,"\n");

      if(Components[i].NumberOfIdentityChanges>0)
      {
        fprintf(FilePtr,"\tnumber of identity changes: %d\n",
                Components[i].NumberOfIdentityChanges);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfIdentityChanges;j++)
          fprintf(FilePtr,"\t%d (%s) to %d (%s)\n",i,Components[i].Name,
          Components[i].IdentityChanges[j],Components[Components[i].IdentityChanges[j]].Name);
        fprintf(FilePtr,"\n");
      }
      if(Components[i].NumberOfIdentityConfigMoves>0)
      {
        fprintf(FilePtr,"\tnumber of identity-config changes: %d\n",
                Components[i].NumberOfIdentityConfigMoves);
        fprintf(FilePtr,"\t--------------------------------------------\n");
        for(j=0;j<Components[i].NumberOfIdentityConfigMoves;j++)
        {
          fprintf(FilePtr,"\tnr fixed %d: ",Components[i].NumberOfUnchangedAtomsIdentityConfig[j]);
          for(k=0;k<Components[i].NumberOfUnchangedAtomsIdentityConfig[j];k++)
            fprintf(FilePtr,"%d ",Components[i].UnchangedAtomsIdentityConfig[j][k]);
          fprintf(FilePtr,"\n");
        }
      }
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"\tNumber of pockets blocked in a unitcell: %d\n",
        Components[i].NumberOfBlockCenters[system]);
      if(Components[i].BlockPockets[system])
      {
        fprintf(FilePtr,"\tPockets are blocked for this component\n");
        fprintf(FilePtr,"\tBlock-pockets Filename: %s.block\n",
          Components[i].BlockPocketsFilename[system]);
        for(j=0;j<Components[i].NumberOfBlockCenters[system];j++)
          fprintf(FilePtr,"\t\tBlock center: %lf %lf %lf  distance: %lf\n",
            (double)Components[i].BlockCenters[system][j].x,
            (double)Components[i].BlockCenters[system][j].y,
            (double)Components[i].BlockCenters[system][j].z,
            (double)Components[i].BlockDistance[system][j]);
      }
      else
        fprintf(FilePtr,"\t\tPockets are NOT blocked for this component\n");

      fprintf(FilePtr,"\n");

    }
    fprintf(FilePtr,"\n\n");

    if(Framework[system].FrameworkModel!=NONE)
    {
      fprintf(FilePtr,"Framework Status\n");
      fprintf(FilePtr,"===========================================================================\n");

      if(Lowenstein[system])
        fprintf(FilePtr,"Lowenstein's rule obeyed by framework\n");
      else
        fprintf(FilePtr,"WARNING: Lowenstein's rule NOT obeyed by framework !!!!!!!!!!\n");

      // print information for a flexible framework
      // for a rigid framework, print the bond-dipoles
      switch(Framework[system].FrameworkModel)
      {
        case FLEXIBLE:
          fprintf(FilePtr,"\tFramework is modelled as: flexible using the model of %s\n",Framework[system].FrameworkDefinitions);
          fprintf(FilePtr,"\t==================================================================================\n");
          fprintf(FilePtr,"\n");

          if(Framework[system].NumberOfCoreShellDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Core/Shell definitions: %d\n",
                    Framework[system].NumberOfCoreShellDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfCoreShellDefinitions;i++)
            {
              fprintf(FilePtr,"\tCore/shell pair: core [%s] and shell [%s] (%d pairs)\n",
                   PseudoAtoms[Framework[system].CoreShellDefinitions[i].A].Name,
                   PseudoAtoms[Framework[system].CoreShellDefinitions[i].B].Name,
                   Framework[system].NumberOfCoreShellsPerType[i]);
            }
            fprintf(FilePtr,"\n");
          }

          if(Framework[system].NumberOfBondsDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Bonds definitions: %d\n",
                    Framework[system].NumberOfBondsDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfBondsDefinitions;i++)
            {
               fprintf(FilePtr,"\tBond interaction between [%s-%s] using bond-type: %s (%d bonds)\n",
                 PseudoAtoms[Framework[system].BondDefinitions[i].A].Name,
                 PseudoAtoms[Framework[system].BondDefinitions[i].B].Name,
                 BondTypes[Framework[system].BondDefinitionType[i]].Name,
                 Framework[system].NumberOfBondsPerType[i]);
                 fprintf(FilePtr,"\t  args: %d   ",BondTypes[Framework[system].BondDefinitionType[i]].nr_args);
                 for(k=0;k<BondTypes[Framework[system].BondDefinitionType[i]].nr_args;k++)
                   fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].BondArgumentDefinitions[i][k]);
                 fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }

          if(Framework[system].NumberOfBondDipoleDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of BondDipole definitions: %d\n",
                    Framework[system].NumberOfBondDipoleDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfBondDipoleDefinitions;i++)
            {
               fprintf(FilePtr,"\tBond dipole interaction for [%s-%s] with magnitude: %lf [D] (%d bonds)\n",
                 PseudoAtoms[Framework[system].BondDipoleDefinitions[i].A].Name,
                 PseudoAtoms[Framework[system].BondDipoleDefinitions[i].B].Name,
                 (double)(Framework[system].BondDipoleArgumentDefinition[i]*DEBYE_CONVERSION_FACTOR),
                 Framework[system].NumberOfBondDipolesPerType[i]);
            }
            fprintf(FilePtr,"\n");
          }


          if(Framework[system].NumberOfUreyBradleyDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Urey-Bradley definitions: %d\n",
                    Framework[system].NumberOfUreyBradleyDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfUreyBradleyDefinitions;i++)
            {
               fprintf(FilePtr,"\tUrey-Bradley interaction between [%s-(%s)-%s] using Urey-Bradley-type:%s (%d Urey-Bradleys)\n",
                 PseudoAtoms[Framework[system].UreyBradleyDefinitions[i].A].Name,
                 PseudoAtoms[Framework[system].UreyBradleyDefinitions[i].B].Name,
                 PseudoAtoms[Framework[system].UreyBradleyDefinitions[i].C].Name,
                 BondTypes[Framework[system].UreyBradleyDefinitionType[i]].Name,
                 Framework[system].NumberOfUreyBradleysPerType[i]);
                 fprintf(FilePtr,"\t  args: %d   ",BondTypes[Framework[system].UreyBradleyDefinitionType[i]].nr_args);
                 for(k=0;k<BondTypes[Framework[system].UreyBradleyDefinitionType[i]].nr_args;k++)
                   fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].UreyBradleyArgumentDefinitions[i][k]);
                 fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }


          if(Framework[system].NumberOfBendDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Bend definitions: %d\n",
                    Framework[system].NumberOfBendDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfBendDefinitions;i++)
            {
              fprintf(FilePtr,"\tBend interaction between [%s-%s-%s] using bend-type: %s (%d bends)\n",
                PseudoAtoms[Framework[system].BendDefinitions[i].A].Name,
                PseudoAtoms[Framework[system].BendDefinitions[i].B].Name,
                PseudoAtoms[Framework[system].BendDefinitions[i].C].Name,
                BendTypes[Framework[system].BendDefinitionType[i]].Name,
                Framework[system].NumberOfBendsPerType[i]);
                fprintf(FilePtr,"\t  args: %d   ",BendTypes[Framework[system].BendDefinitionType[i]].nr_args);
                for(k=0;k<BendTypes[Framework[system].BendDefinitionType[i]].nr_args;k++)
                  fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].BendArgumentDefinitions[i][k]);
                fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }

          if(Framework[system].NumberOfInversionBendDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Inversion bend definitions: %d\n",
                    Framework[system].NumberOfInversionBendDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfInversionBendDefinitions;i++)
            {
               fprintf(FilePtr,"\tInversion bend interaction between [%s-%s-%s-%s] using inversion-bend type:%s (%d inversion-bends)\n",
                 PseudoAtoms[Framework[system].InversionBendDefinitions[i].A].Name,
                 PseudoAtoms[Framework[system].InversionBendDefinitions[i].B].Name,
                 PseudoAtoms[Framework[system].InversionBendDefinitions[i].C].Name,
                 PseudoAtoms[Framework[system].InversionBendDefinitions[i].D].Name,
                 InversionBendTypes[Framework[system].InversionBendDefinitionType[i]].Name,
                 Framework[system].NumberOfInversionBendsPerType[i]);
                 fprintf(FilePtr,"\t  args: %d   ",InversionBendTypes[Framework[system].InversionBendDefinitionType[i]].nr_args);
                 for(k=0;k<InversionBendTypes[Framework[system].InversionBendDefinitionType[i]].nr_args;k++)
                   fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].InversionBendArgumentDefinitions[i][k]);
                 fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }

          if(Framework[system].NumberOfTorsionDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Torsion definitions: %d\n",
                    Framework[system].NumberOfTorsionDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfTorsionDefinitions;i++)
            {
               fprintf(FilePtr,"\tTorsion interaction between [%s-%s-%s-%s] using torsion-type:%s (%d torsions)\n",
                 PseudoAtoms[Framework[system].TorsionDefinitions[i].A].Name,
                 PseudoAtoms[Framework[system].TorsionDefinitions[i].B].Name,
                 PseudoAtoms[Framework[system].TorsionDefinitions[i].C].Name,
                 PseudoAtoms[Framework[system].TorsionDefinitions[i].D].Name,
                 TorsionTypes[Framework[system].TorsionDefinitionType[i]].Name,
                 Framework[system].NumberOfTorsionsPerType[i]);
                 fprintf(FilePtr,"\t  args: %d   ",TorsionTypes[Framework[system].TorsionDefinitionType[i]].nr_args);
                 for(k=0;k<TorsionTypes[Framework[system].TorsionDefinitionType[i]].nr_args;k++)
                   fprintf(FilePtr,"arg[%d]=%-12.6lf ",k,(double)Framework[system].TorsionArgumentDefinitions[i][k]);
                 fprintf(FilePtr," VDW-scaling: %-12.6lf  q-scaling: %-12.6lf\n",Framework[system].TorsionArgumentDefinitions[i][6],
                   Framework[system].TorsionArgumentDefinitions[i][7]);
            }
            fprintf(FilePtr,"\n");
          }

          if(Framework[system].NumberOfImproperTorsionDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Improper Torsion definitions: %d\n",
                    Framework[system].NumberOfImproperTorsionDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfImproperTorsionDefinitions;i++)
            {
               fprintf(FilePtr,"\tImproper Torsion interaction between [%s-%s-%s-%s] using improper torsion-type:%s (%d improper torsions)\n",
                 PseudoAtoms[Framework[system].ImproperTorsionDefinitions[i].A].Name,
                 PseudoAtoms[Framework[system].ImproperTorsionDefinitions[i].B].Name,
                 PseudoAtoms[Framework[system].ImproperTorsionDefinitions[i].C].Name,
                 PseudoAtoms[Framework[system].ImproperTorsionDefinitions[i].D].Name,
                 ImproperTorsionTypes[Framework[system].ImproperTorsionDefinitionType[i]].Name,
                 Framework[system].NumberOfImproperTorsionsPerType[i]);
                 fprintf(FilePtr,"\t  args: %d   ",ImproperTorsionTypes[Framework[system].ImproperTorsionDefinitionType[i]].nr_args);
                 for(k=0;k<ImproperTorsionTypes[Framework[system].ImproperTorsionDefinitionType[i]].nr_args;k++)
                   fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].ImproperTorsionArgumentDefinitions[i][k]);
                 fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }

          if(Framework[system].NumberOfOutOfPlaneDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Out-of-plane distance definitions: %d\n",
                    Framework[system].NumberOfOutOfPlaneDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfOutOfPlaneDefinitions;i++)
            {
               fprintf(FilePtr,"\tOut-of-plane interaction between [%s-%s-%s-%s] using out-of-plane-type:%s (%d out-of-planes)\n",
                 PseudoAtoms[Framework[system].OutOfPlaneDefinitions[i].A].Name,
                 PseudoAtoms[Framework[system].OutOfPlaneDefinitions[i].B].Name,
                 PseudoAtoms[Framework[system].OutOfPlaneDefinitions[i].C].Name,
                 PseudoAtoms[Framework[system].OutOfPlaneDefinitions[i].D].Name,
                 OutOfPlaneTypes[Framework[system].OutOfPlaneDefinitionType[i]].Name,
                 Framework[system].NumberOfOutOfPlanesPerType[i]);
                 fprintf(FilePtr,"\t  args: %d   ",OutOfPlaneTypes[Framework[system].OutOfPlaneDefinitionType[i]].nr_args);
                 for(k=0;k<OutOfPlaneTypes[Framework[system].OutOfPlaneDefinitionType[i]].nr_args;k++)
                   fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].OutOfPlaneArgumentDefinitions[i][k]);
                 fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }


          if(Framework[system].NumberOfBondBondDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Bond-bond cross-term definitions: %d\n",
                    Framework[system].NumberOfBondBondDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfBondBondDefinitions;i++)
            {
              fprintf(FilePtr,"\tBond-bond cross-term interaction between [%s-%s-%s] using bond/bond-type:%s (%d bond-bonds\n",
                PseudoAtoms[Framework[system].BondBondDefinitions[i].A].Name,
                PseudoAtoms[Framework[system].BondBondDefinitions[i].B].Name,
                PseudoAtoms[Framework[system].BondBondDefinitions[i].C].Name,
                BondBondTypes[Framework[system].BondBondDefinitionType[i]].Name,
                Framework[system].NumberOfBondBondsPerType[i]);
                fprintf(FilePtr,"\t  args: %d   ",BondBondTypes[Framework[system].BondBondDefinitionType[i]].nr_args);
                for(k=0;k<BondBondTypes[Framework[system].BondBondDefinitionType[i]].nr_args;k++)
                  fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].BondBondArgumentDefinitions[i][k]);
                fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }

          if(Framework[system].NumberOfBondBendDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Bond/Bend cross-term definitions: %d\n",
                    Framework[system].NumberOfBondBendDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfBondBendDefinitions;i++)
            {
              fprintf(FilePtr,"\tBond-bend cross-term interaction between [%s-%s-%s] using bond/bend-type:%s (%d bond/bends)\n",
                PseudoAtoms[Framework[system].BondBendDefinitions[i].A].Name,
                PseudoAtoms[Framework[system].BondBendDefinitions[i].B].Name,
                PseudoAtoms[Framework[system].BondBendDefinitions[i].C].Name,
                BondBendTypes[Framework[system].BondBendDefinitionType[i]].Name,
                Framework[system].NumberOfBondBendsPerType[i]);
                fprintf(FilePtr,"\t  args: %d   ",BondBendTypes[Framework[system].BondBendDefinitionType[i]].nr_args);
                for(k=0;k<BondBendTypes[Framework[system].BondBendDefinitionType[i]].nr_args;k++)
                  fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].BondBendArgumentDefinitions[i][k]);
                fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }

          if(Framework[system].NumberOfBendBendDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Bend/Bend cross-term definitions: %d\n",
                    Framework[system].NumberOfBendBendDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfBendBendDefinitions;i++)
            {
              fprintf(FilePtr,"\tBend-bend cross-term interaction between [%s-%s-%s-%s] using bend/bend-type:%s (%d bend-bends)\n",
                PseudoAtoms[Framework[system].BendBendDefinitions[i].A].Name,
                PseudoAtoms[Framework[system].BendBendDefinitions[i].B].Name,
                PseudoAtoms[Framework[system].BendBendDefinitions[i].C].Name,
                PseudoAtoms[Framework[system].BendBendDefinitions[i].D].Name,
                BendBendTypes[Framework[system].BendBendDefinitionType[i]].Name,
                Framework[system].NumberOfBendBendsPerType[i]);
                fprintf(FilePtr,"\t  args: %d   ",BendBendTypes[Framework[system].BendBendDefinitionType[i]].nr_args);
                for(k=0;k<BendBendTypes[Framework[system].BendBendDefinitionType[i]].nr_args;k++)
                  fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].BendBendArgumentDefinitions[i][k]);
                fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }

          if(Framework[system].NumberOfBondTorsionDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Bond/Torsion cross-term definitions: %d\n",
                    Framework[system].NumberOfBondTorsionDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfBondTorsionDefinitions;i++)
            {
              fprintf(FilePtr,"\tBond/Torsion cross-term interaction between [%s-%s-%s-%s] using stretch/torsion-type:%s (%d stretch-torsions)\n",
                PseudoAtoms[Framework[system].BondTorsionDefinitions[i].A].Name,
                PseudoAtoms[Framework[system].BondTorsionDefinitions[i].B].Name,
                PseudoAtoms[Framework[system].BondTorsionDefinitions[i].C].Name,
                PseudoAtoms[Framework[system].BondTorsionDefinitions[i].D].Name,
                BondTorsionTypes[Framework[system].BondTorsionDefinitionType[i]].Name,
                Framework[system].NumberOfBondTorsionsPerType[i]);
                fprintf(FilePtr,"\t  args: %d   ",BondTorsionTypes[Framework[system].BondTorsionDefinitionType[i]].nr_args);
                for(k=0;k<BondTorsionTypes[Framework[system].BondTorsionDefinitionType[i]].nr_args;k++)
                  fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].BondTorsionArgumentDefinitions[i][k]);
                fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }

          if(Framework[system].NumberOfBendTorsionDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of Bend/Torsion cross-term definitions: %d\n",
                    Framework[system].NumberOfBendTorsionDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfBendTorsionDefinitions;i++)
            {
              fprintf(FilePtr,"\tBend/Torsion cross-term interaction between [%s-%s-%s-%s] using bend/torsion-type:%s (%d bend-torsions)\n",
                PseudoAtoms[Framework[system].BendTorsionDefinitions[i].A].Name,
                PseudoAtoms[Framework[system].BendTorsionDefinitions[i].B].Name,
                PseudoAtoms[Framework[system].BendTorsionDefinitions[i].C].Name,
                PseudoAtoms[Framework[system].BendTorsionDefinitions[i].D].Name,
                BendTorsionTypes[Framework[system].BendTorsionDefinitionType[i]].Name,
                Framework[system].NumberOfBendTorsionsPerType[i]);
                nr_args=BendTorsionTypes[Framework[system].BendTorsionDefinitionType[i]].nr_args;
                fprintf(FilePtr,"\t  args: %d   ",nr_args);
                for(k=0;k<nr_args;k++)
                  fprintf(FilePtr,"arg[%d]=%lf ",k,(double)Framework[system].BendTorsionArgumentDefinitions[i][k]);
                fprintf(FilePtr,"\n");
            }
            fprintf(FilePtr,"\n");
          }

          for(CurrentFramework=0;CurrentFramework<Framework[system].NumberOfFrameworks;CurrentFramework++)
          {
            fprintf(FilePtr,"\n");
            fprintf(FilePtr,"\tNumber of core/shells:                           %d\n",Framework[system].NumberOfCoreShells[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of bonds:                                 %d\n",Framework[system].NumberOfBonds[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of bond-dipoles:                          %d\n",Framework[system].NumberOfBondDipoles[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of bends:                                 %d\n",Framework[system].NumberOfBends[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of Urey-Bradleys:                         %d\n",Framework[system].NumberOfUreyBradleys[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of Inversion bends:                       %d\n",Framework[system].NumberOfInversionBends[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of Torsions:                              %d\n",Framework[system].NumberOfTorsions[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of Improper Torsions:                     %d\n",Framework[system].NumberOfImproperTorsions[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of Out-of-planes:                         %d\n",Framework[system].NumberOfOutOfPlanes[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of bond-bond cross terms:                 %d\n",Framework[system].NumberOfBondBonds[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of bond-bend cross terms:                 %d\n",Framework[system].NumberOfBondBends[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of bend-bend cross terms:                 %d\n",Framework[system].NumberOfBendBends[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of stretch-torsion cross terms:           %d\n",Framework[system].NumberOfBondTorsions[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of bend-torsion cross terms:              %d\n",Framework[system].NumberOfBendTorsions[CurrentFramework]);
            fprintf(FilePtr,"\n");

            fprintf(FilePtr,"\tNumber of charges:                               %d\n",Framework[system].NumberOfCharges[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of bonddipoles:                           %d\n",Framework[system].NumberOfBondDipoles[CurrentFramework]);
            fprintf(FilePtr,"\n");

            fprintf(FilePtr,"\tNumber of Intra VDW:                             %d\n",Framework[system].NumberOfIntraVDW[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of Intra Coulomb charge-charge:           %d\n",Framework[system].NumberOfIntraCharges[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of Intra Coulomb charge-bonddipole:       %d\n",Framework[system].NumberOfIntraChargeBondDipoles[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of Intra Coulomb bonddipoles:             %d\n",Framework[system].NumberOfIntraBondDipoles[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of excluded Intra VDW:                    %d\n",Framework[system].NumberOfExcludedIntraVDW[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of excluded intra charge-charge:          %d\n",Framework[system].NumberOfExcludedIntraChargeCharge[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of excluded intra charge-bonddipole:      %d\n",Framework[system].NumberOfExcludedIntraChargeBondDipole[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of excluded intra bonddipole-bonddipole:  %d\n",Framework[system].NumberOfExcludedIntraBondDipoleBondDipole[CurrentFramework]);
          }
          fprintf(FilePtr,"\n\n");
          fprintf(FilePtr,"\tNumber of 1-2 intra-framework pairs: %d\n",Framework[system].NumberOfIntra12Interactions);
          fprintf(FilePtr,"\tNumber of 1-3 intra-framework pairs: %d\n",Framework[system].NumberOfIntra13Interactions);
          fprintf(FilePtr,"\tNumber of 1-4 intra-framework pairs: %d\n",Framework[system].NumberOfIntra14Interactions);
          fprintf(FilePtr,"\tNumber of 1-2,1-3 intra-framework pairs: %d\n",Framework[system].NumberOfIntra123Interactions);
          fprintf(FilePtr,"\tNumber of 1-2,1-3,1-4 intra-framework pairs: %d\n",Framework[system].NumberOfIntra1234Interactions);
          if(RemoveBondNeighboursFromLongRangeInteraction)
            fprintf(FilePtr,"\tBond interactions are removed (excluded from VDW and Coulomb interactions)\n");
          else
            fprintf(FilePtr,"\tBond interactions are *included* in the VDW and Coulomb interactions\n");
          if(RemoveBendNeighboursFromLongRangeInteraction)
            fprintf(FilePtr,"\tBend interactions are removed (excluded from VDW and Coulomb interactions)\n");
          else
            fprintf(FilePtr,"\tBend interactions are *included* in the VDW and Coulomb interactions\n");
          if(RemoveTorsionNeighboursFromLongRangeInteraction)
            fprintf(FilePtr,"\tTorsion interactions are removed (excluded from VDW and Coulomb interactions)\n");
          else
            fprintf(FilePtr,"\tTorsion interactions are *included* in the VDW and Coulomb interactions\n");
          fprintf(FilePtr,"\tRemove12NeighboursFromVDWInteraction %s\n",Remove12NeighboursFromVDWInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove13NeighboursFromVDWInteraction %s\n",Remove13NeighboursFromVDWInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove14NeighboursFromVDWInteraction %s\n",Remove14NeighboursFromVDWInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove12NeighboursFromChargeChargeInteraction %s\n",Remove12NeighboursFromChargeChargeInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove13NeighboursFromChargeChargeInteraction %s\n",Remove13NeighboursFromChargeChargeInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove14NeighboursFromChargeChargeInteraction %s\n",Remove14NeighboursFromChargeChargeInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove11NeighboursFromChargeBondDipoleInteraction %s\n",Remove11NeighboursFromChargeBondDipoleInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove12NeighboursFromChargeBondDipoleInteraction %s\n",Remove12NeighboursFromChargeBondDipoleInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove13NeighboursFromChargeBondDipoleInteraction %s\n",Remove13NeighboursFromChargeBondDipoleInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove14NeighboursFromChargeBondDipoleInteraction %s\n",Remove14NeighboursFromChargeBondDipoleInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove12NeighboursFromBondDipoleBondDipoleInteraction %s\n",Remove12NeighboursFromBondDipoleBondDipoleInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove13NeighboursFromBondDipoleBondDipoleInteraction %s\n",Remove13NeighboursFromBondDipoleBondDipoleInteraction?"yes":"no");
          fprintf(FilePtr,"\tRemove14NeighboursFromBondDipoleBondDipoleInteraction %s\n",Remove14NeighboursFromBondDipoleBondDipoleInteraction?"yes":"no");
          fprintf(FilePtr,"\tInternalFrameworkLennardJonesInteractions %s\n",InternalFrameworkLennardJonesInteractions?"yes":"no");

          fprintf(FilePtr,"\n");
          break;
        default:
          fprintf(FilePtr,"\tFramework is modelled as: rigid\n\n");

          if(Framework[system].NumberOfBondDipoleDefinitions>0)
          {
            fprintf(FilePtr,"\tNumber of BondDipole definitions: %d\n",
                    Framework[system].NumberOfBondDipoleDefinitions);
            fprintf(FilePtr,"\t--------------------------------------------\n");
            for(i=0;i<Framework[system].NumberOfBondDipoleDefinitions;i++)
            {
               fprintf(FilePtr,"\t\tBond dipole interaction for [%s-%s] with magnitude: %lf [D] (%d bonds)\n",
                 PseudoAtoms[Framework[system].BondDipoleDefinitions[i].A].Name,
                 PseudoAtoms[Framework[system].BondDipoleDefinitions[i].B].Name,
                 (double)(Framework[system].BondDipoleArgumentDefinition[i]*DEBYE_CONVERSION_FACTOR),
                 Framework[system].NumberOfBondDipolesPerType[i]);
            }
            fprintf(FilePtr,"\n");
          }
          for(CurrentFramework=0;CurrentFramework<Framework[system].NumberOfFrameworks;CurrentFramework++)
          {
            fprintf(FilePtr,"\tNumber of charges:                               %d\n",Framework[system].NumberOfCharges[CurrentFramework]);
            fprintf(FilePtr,"\tNumber of bonddipoles:                           %d\n",Framework[system].NumberOfBondDipoles[CurrentFramework]);
            fprintf(FilePtr,"\n");
          }
          break;
      }
      fprintf(FilePtr,"\n");
    }
  }


  fprintf(FilePtr,"System Properties\n");
  fprintf(FilePtr,"===========================================================================\n");

  if(BoundaryCondition[system]==FINITE)
  {
    fprintf(FilePtr,"Finite size system: no boundary conditions are applied\n");
    fprintf(FilePtr,"\n");
  }
  else
  {
    if(UseCellLists[system])
    {
      fprintf(FilePtr,"Energy and force calculations use the cell-list method: %d %d %d cells\n\n",
       NumberOfCellListCells[system].x,NumberOfCellListCells[system].y,NumberOfCellListCells[system].z);
    }

  switch(Dimension)
  {
    case 2:
      fprintf(FilePtr,"Unit cell size: %lf %lf\n",
              (double)UnitCellSize[system].x,(double)UnitCellSize[system].y);
      fprintf(FilePtr,"Cell angles (radians)  alpha: %lf\n",(double)AlphaAngle[system]);
      fprintf(FilePtr,"Cell angles (degrees)  alpha: %lf\n",(double)AlphaAngle[system]*180.0/M_PI);
      fprintf(FilePtr,"Number of unitcells [a]: %d\n",NumberOfUnitCells[system].x);
      fprintf(FilePtr,"Number of unitcells [b]: %d\n",NumberOfUnitCells[system].y);
      fprintf(FilePtr,"\n");

      switch(BoundaryCondition[system])
      {
        case NONE:
          fprintf(FilePtr,"No periodic boundaries\n");
          break;
        case CUBIC:
          fprintf(FilePtr,"CUBIC Boundary conditions: alpha=90 degrees, a=b=c\n");
          break;
        case RECTANGULAR:
          fprintf(FilePtr,"RECTANGULAR Boundary conditions: alpha=90 degrees\n");
          break;
        case TRICLINIC:
          fprintf(FilePtr,"TRICLINIC Boundary conditions: alpha!=90 or beta!=90 or gamma!=90\n");
      }
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"Cartesian axis A is collinear with crystallographic axis a\n");
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"lengths of cell vectors:\n");
      fprintf(FilePtr,"%9.5lf %9.5lf\n",(double)BoxProperties[system].ax,(double)BoxProperties[system].ay);
      fprintf(FilePtr,"cosines of cell angles:\n");
      fprintf(FilePtr,"%9.5lf\n",(double)BoxProperties[system].bx);
      fprintf(FilePtr,"perpendicular cell widths:\n");
      fprintf(FilePtr,"%9.5lf %9.5lf\n",(double)BoxProperties[system].cx,(double)BoxProperties[system].cy);
      fprintf(FilePtr,"area of the cell: %18.12lf (A^2)\n",(double)Volume[system]);
      fprintf(FilePtr,"\n");


      fprintf(FilePtr,"Orthogonalization matrix Box\n");
      fprintf(FilePtr,"Transforms fractional coordinates abc into orthonormal Cartesian coordinates xyz\n");
      fprintf(FilePtr,"Deorthogonalization matrix InverseBox\n");
      fprintf(FilePtr,"Transforms orthonormal Cartesian coordinates xyz into fractional coordinates xyz\n");
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"Box[%d]:\n",system);
      fprintf(FilePtr,"\t%18.12lf %18.12lf\n",(double)Box[system].ax,(double)Box[system].bx);
      fprintf(FilePtr,"\t%18.12lf %18.12lf\n",(double)Box[system].ay,(double)Box[system].by);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"Inverse box[%d]:\n",system);
      fprintf(FilePtr,"\t%18.12lf %18.12lf\n",(double)InverseBox[system].ax,(double)InverseBox[system].bx);
      fprintf(FilePtr,"\t%18.12lf %18.12lf\n",(double)InverseBox[system].ay,(double)InverseBox[system].by);
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"Unitcell box[%d]:\n",system);
      fprintf(FilePtr,"\t%18.12lf %18.12lf\n",(double)UnitCellBox[system].ax,(double)UnitCellBox[system].bx);
      fprintf(FilePtr,"\t%18.12lf %18.12lf\n",(double)UnitCellBox[system].ay,(double)UnitCellBox[system].by);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"Unitcell inverse box[%d]:\n",system);
      fprintf(FilePtr,"\t%18.12lf %18.12lf\n",(double)InverseUnitCellBox[system].ax,(double)InverseUnitCellBox[system].bx);
      fprintf(FilePtr,"\t%18.12lf %18.12lf\n",(double)InverseUnitCellBox[system].ay,(double)InverseUnitCellBox[system].by);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"lengths of cell vectors (inverse box):\n");
      fprintf(FilePtr,"%9.5lf %9.5lf\n",(double)InverseBoxProperties[system].ax,(double)InverseBoxProperties[system].ay);
      fprintf(FilePtr,"cosine of cell angles (inverse box):\n");
      fprintf(FilePtr,"%9.5lf\n",(double)InverseBoxProperties[system].bx);
      fprintf(FilePtr,"perpendicular cell widths (inverse):\n");
      fprintf(FilePtr,"%9.5lf %9.5lf\n",(double)InverseBoxProperties[system].cx,(double)InverseBoxProperties[system].cy);
      fprintf(FilePtr,"area of the cell: %18.12lf (A^2)\n",(double)Volume[system]);
      break;
    case 3:
      fprintf(FilePtr,"Unit cell size: %lf %lf %lf\n",
              (double)UnitCellSize[system].x,(double)UnitCellSize[system].y,(double)UnitCellSize[system].z);
      fprintf(FilePtr,"Cell angles (radians)  alpha: %lf beta: %lf gamma: %lf\n",
              (double)AlphaAngle[system],(double)BetaAngle[system],(double)GammaAngle[system]);
      fprintf(FilePtr,"Cell angles (degrees)  alpha: %lf beta: %lf gamma: %lf\n",
              (double)AlphaAngle[system]*180.0/M_PI,(double)BetaAngle[system]*180.0/M_PI,(double)GammaAngle[system]*180.0/M_PI);
      fprintf(FilePtr,"Number of unitcells [a]: %d\n",NumberOfUnitCells[system].x);
      fprintf(FilePtr,"Number of unitcells [b]: %d\n",NumberOfUnitCells[system].y);
      fprintf(FilePtr,"Number of unitcells [c]: %d\n",NumberOfUnitCells[system].z);
      fprintf(FilePtr,"\n");

      switch(BoundaryCondition[system])
      {
        case NONE:
          fprintf(FilePtr,"No periodic boundaries\n");
          break;
        case CUBIC:
          fprintf(FilePtr,"CUBIC Boundary conditions: alpha=beta=gamma=90 degrees, a=b=c\n");
          break;
        case RECTANGULAR:
          fprintf(FilePtr,"RECTANGULAR Boundary conditions: alpha=beta=gamma=90 degrees\n");
          break;
        case TRICLINIC:
          fprintf(FilePtr,"TRICLINIC Boundary conditions: alpha!=90 or beta!=90 or gamma!=90\n");
      }
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"Cartesian axis A is collinear with crystallographic axis a\n");
      fprintf(FilePtr,"Cartesian axis B is collinear with (axb)xA\n");
      fprintf(FilePtr,"Cartesian axis C is collinear with (axb)\n");
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"lengths of cell vectors:\n");
      fprintf(FilePtr,"%9.5lf %9.5lf %9.5lf\n",(double)BoxProperties[system].ax,
         (double)BoxProperties[system].ay,(double)BoxProperties[system].az);
      fprintf(FilePtr,"cosines of cell angles:\n");
      fprintf(FilePtr,"%9.5lf %9.5lf %9.5lf\n",(double)BoxProperties[system].bx,
         (double)BoxProperties[system].by,(double)BoxProperties[system].bz);
      fprintf(FilePtr,"perpendicular cell widths:\n");
      fprintf(FilePtr,"%9.5lf %9.5lf %9.5lf\n",(double)BoxProperties[system].cx,
         (double)BoxProperties[system].cy,(double)BoxProperties[system].cz);
      fprintf(FilePtr,"volume of the cell: %18.12lf (A^3)\n",(double)Volume[system]);
      fprintf(FilePtr,"\n");


      fprintf(FilePtr,"Orthogonalization matrix Box\n");
      fprintf(FilePtr,"Transforms fractional coordinates abc into orthonormal Cartesian coordinates xyz\n");
      fprintf(FilePtr,"Deorthogonalization matrix InverseBox\n");
      fprintf(FilePtr,"Transforms orthonormal Cartesian coordinates xyz into fractional coordinates xyz\n");
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"Box[%d]:\n",system);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)Box[system].ax,(double)Box[system].bx,(double)Box[system].cx);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)Box[system].ay,(double)Box[system].by,(double)Box[system].cy);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)Box[system].az,(double)Box[system].bz,(double)Box[system].cz);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"Inverse box[%d]:\n",system);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)InverseBox[system].ax,
         (double)InverseBox[system].bx,(double)InverseBox[system].cx);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)InverseBox[system].ay,
          (double)InverseBox[system].by,(double)InverseBox[system].cy);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)InverseBox[system].az,
          (double)InverseBox[system].bz,(double)InverseBox[system].cz);
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"Unitcell box[%d]:\n",system);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)UnitCellBox[system].ax,
         (double)UnitCellBox[system].bx,(double)UnitCellBox[system].cx);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)UnitCellBox[system].ay,
         (double)UnitCellBox[system].by,(double)UnitCellBox[system].cy);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)UnitCellBox[system].az,
         (double)UnitCellBox[system].bz,(double)UnitCellBox[system].cz);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"Unitcell inverse box[%d]:\n",system);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)InverseUnitCellBox[system].ax,
         (double)InverseUnitCellBox[system].bx,(double)InverseUnitCellBox[system].cx);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)InverseUnitCellBox[system].ay,
         (double)InverseUnitCellBox[system].by,(double)InverseUnitCellBox[system].cy);
      fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)InverseUnitCellBox[system].az,
         (double)InverseUnitCellBox[system].bz,(double)InverseUnitCellBox[system].cz);
      fprintf(FilePtr,"\n");

      fprintf(FilePtr,"lengths of cell vectors (inverse box):\n");
      fprintf(FilePtr,"%9.5lf %9.5lf %9.5lf\n",(double)InverseBoxProperties[system].ax,
         (double)InverseBoxProperties[system].ay,(double)InverseBoxProperties[system].az);
      fprintf(FilePtr,"cosines of cell angles (inverse box):\n");
      fprintf(FilePtr,"%9.5lf %9.5lf %9.5lf\n",(double)InverseBoxProperties[system].bx,
         (double)InverseBoxProperties[system].by,(double)InverseBoxProperties[system].bz);
      fprintf(FilePtr,"perpendicular cell widths (inverse):\n");
      fprintf(FilePtr,"%9.5lf %9.5lf %9.5lf\n",(double)InverseBoxProperties[system].cx,
         (double)InverseBoxProperties[system].cy,(double)InverseBoxProperties[system].cz);
      fprintf(FilePtr,"volume of the cell: %18.12lf (A^3)\n",(double)Volume[system]);
      break;
  }
  fprintf(FilePtr,"\n");

  if(UseReplicas[system])
  {
    fprintf(FilePtr,"Number of Replica Cells: %d %d %d\n",NumberOfReplicaCells[system].x,NumberOfReplicaCells[system].y,NumberOfReplicaCells[system].z);
    fprintf(FilePtr,"Replica Box:\n");
    fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)ReplicaBox[system].ax,(double)ReplicaBox[system].bx,(double)ReplicaBox[system].cx);
    fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)ReplicaBox[system].ay,(double)ReplicaBox[system].by,(double)ReplicaBox[system].cy);
    fprintf(FilePtr,"\t%18.12lf %18.12lf %18.12lf\n",(double)ReplicaBox[system].az,(double)ReplicaBox[system].bz,(double)ReplicaBox[system].cz);
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Shift of the Replica cells:\n");
    ncell=0;
    for(k1=0;k1<NumberOfReplicaCells[system].x;k1++)
      for(k2=0;k2<NumberOfReplicaCells[system].y;k2++)
        for(k3=0;k3<NumberOfReplicaCells[system].z;k3++)
        {
          fprintf(FilePtr,"\t[cell: %2d] %18.12f %18.12f %18.12f\n",ncell,ReplicaShift[ncell].x,ReplicaShift[ncell].y,ReplicaShift[ncell].z);
          ncell++;
        }
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"No replicas are used\n");
  }

  if(Framework[system].FrameworkModel==FLEXIBLE)
    fprintf(FilePtr,"Framework is simulated as 'flexible'\n");
  else
    fprintf(FilePtr,"Framework is simulated as 'rigid'\n");

  fprintf(FilePtr,"Number of framework atoms: %d\n",Framework[system].TotalNumberOfAtoms);
  fprintf(FilePtr,"Number of framework atoms in the unit cell: %d\n",Framework[system].TotalNumberOfUnitCellAtoms);
  fprintf(FilePtr,"Framework Mass: %18.12lf [g/mol]\n",(double)Framework[system].FrameworkMass);
  fprintf(FilePtr,"Framework Density: %18.12lf [kg/m^3] %18.13lf [cm^3/g]\n",
       (double)(1e-3*Framework[system].FrameworkDensity),
        (double)(1e6*Volume[system]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT)/Framework[system].FrameworkMass);
  fprintf(FilePtr,"Helium void fraction: %13.8lf\n",(double)HeliumVoidFraction[system]);
  fprintf(FilePtr,"Available pore volume: %13.8lf [A^3] %13.8lf [cm^3/g]\n",(double)(HeliumVoidFraction[system]*Volume[system]),
        HeliumVoidFraction[system]*(1e6*Volume[system]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT)/Framework[system].FrameworkMass);
  fprintf(FilePtr,"Conversion factor from molecule/unit cell -> kmol/m^3: %lg, kmol/m^3 accesible pore volume: %lg\n",
     1.0/(1e3*Volume[system]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT),
     1.0/(HeliumVoidFraction[system]*(1e3*Volume[system]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT)));

  if(NumberOfModificationRules>0)
  {
    fprintf(FilePtr,"The framework atom types are modified at reading the input, %d rules\n",NumberOfModificationRules);
    for(i=0;i<NumberOfModificationRules;i++)
    {
      switch(ModificationRuleType[i])
      {
        case MODIFY_FRAMEWORKATOM_CONNECTED_TO:
          fprintf(FilePtr,"\trule %d: modify (%s) to (%s) when connected to (%s) and (%s)\n",i,
            ModifyFrameworkAtoms[i][0],ModifyFrameworkAtoms[i][1],ModifyFrameworkAtoms[i][2],ModifyFrameworkAtoms[i][3]);
          break;
        case MODIFY_FRAMEWORKATOM_DIMER:
          fprintf(FilePtr,"\trule %d: modify (%s-%s) to (%s-%s)\n",i,
            ModifyFrameworkAtoms[i][0],ModifyFrameworkAtoms[i][1],ModifyFrameworkAtoms[i][2],ModifyFrameworkAtoms[i][3]);
          break;
        case MODIFY_FRAMEWORKATOM_TRIPLE:
          fprintf(FilePtr,"\trule %d: modify (%s-%s-%s) to (%s-%s-%s)\n",i,
            ModifyFrameworkAtoms[i][0],ModifyFrameworkAtoms[i][1],ModifyFrameworkAtoms[i][2],
            ModifyFrameworkAtoms[i][3],ModifyFrameworkAtoms[i][4],ModifyFrameworkAtoms[i][5]);
          break;
        case MODIFY_FRAMEWORKATOM_PLANAR:
          fprintf(FilePtr,"\trule %d: modify (%s-%s-%s-%s-%s) to (%s-%s-%s-%s-%s)\n",i,
            ModifyFrameworkAtoms[i][0],ModifyFrameworkAtoms[i][1],ModifyFrameworkAtoms[i][2],ModifyFrameworkAtoms[i][3],ModifyFrameworkAtoms[i][4],
            ModifyFrameworkAtoms[i][5],ModifyFrameworkAtoms[i][6],ModifyFrameworkAtoms[i][7],ModifyFrameworkAtoms[i][8],ModifyFrameworkAtoms[i][9]);
          break;
      }
    }
  }
  if(NumberOfForbiddenConnectivityRules>0)
  {
    fprintf(FilePtr,"Forbidden rules for framework atom types, %d rules\n",NumberOfForbiddenConnectivityRules);
    for(i=0;i<NumberOfForbiddenConnectivityRules;i++)
      fprintf(FilePtr,"\trule %d: forbidden framework sequence (%s-%s-%s)\n",i,
        ForbiddenConnectivityAtoms[i][0],ForbiddenConnectivityAtoms[i][1],ForbiddenConnectivityAtoms[i][2]);
  }

  if(NumberOfSubstitutions>0)
  {
    fprintf(FilePtr,"Numnber of Atom-substitutions: %d (%d fixed, then %d random)\n",NumberOfSubstitutions,NumberOfSingleSubstitutionRules,NumberOfRandomSubstitutions);
    for(i=0;i<NumberOfSubstitutions;i++)
      fprintf(FilePtr,"\tsubstitued the %d-th atom of type %s by type %s\n",ListOfAtomSubstitutions[i][0],
        PseudoAtoms[ListOfAtomSubstitutions[i][1]].Name,PseudoAtoms[ListOfAtomSubstitutions[i][2]].Name);
  }

  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Number Of Frameworks (per system): %d\n",Framework[system].NumberOfFrameworks);
  fprintf(FilePtr,"----------------------------------------------------------------------\n");
  for(CurrentFramework=0;CurrentFramework<Framework[system].NumberOfFrameworks;CurrentFramework++)
  {
    fprintf(FilePtr,"Framework name: %s\n",Framework[system].Name[CurrentFramework]);
    if(Framework[system].NumberOfCitations[CurrentFramework]>0)
    {
       for(i=0;i<Framework[system].NumberOfCitations[CurrentFramework];i++)
       {
         fprintf(FilePtr,"Citation:\n");
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationId)>0)
           fprintf(FilePtr,"\tID:                 %s\n",Framework[system].CitationInformation[CurrentFramework][i].CitationId);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationAuthorName)>0)
           fprintf(FilePtr,"\tauthor name:        '%s'\n",Framework[system].CitationInformation[CurrentFramework][i].CitationAuthorName);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationCoordinateLinkage)>0)
           fprintf(FilePtr,"\tcoordinate linkage: %s\n",Framework[system].CitationInformation[CurrentFramework][i].CitationCoordinateLinkage);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationTitle)>0)
           fprintf(FilePtr,"\ttitle:              '%s'\n",Framework[system].CitationInformation[CurrentFramework][i].CitationTitle);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationCountry)>0)
           fprintf(FilePtr,"\tcountry:            '%s'\n",Framework[system].CitationInformation[CurrentFramework][i].CitationCountry);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationJournalAbbrev)>0)
           fprintf(FilePtr,"\tjournal abbrev.:    '%s'\n",Framework[system].CitationInformation[CurrentFramework][i].CitationJournalAbbrev);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationJournalVolume)>0)
           fprintf(FilePtr,"\tjournal volume:     %s\n",Framework[system].CitationInformation[CurrentFramework][i].CitationJournalVolume);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationJournalIssue)>0)
           fprintf(FilePtr,"\tjournal issue:      %s\n",Framework[system].CitationInformation[CurrentFramework][i].CitationJournalIssue);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationPageFirst)>0)
           fprintf(FilePtr,"\tfirst page:         %s\n",Framework[system].CitationInformation[CurrentFramework][i].CitationPageFirst);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationPageLast)>0)
           fprintf(FilePtr,"\tlast page:          %s\n",Framework[system].CitationInformation[CurrentFramework][i].CitationPageLast);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationYear)>0)
           fprintf(FilePtr,"\tyear:               %s\n",Framework[system].CitationInformation[CurrentFramework][i].CitationYear);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationJournalID_ASTM)>0)
           fprintf(FilePtr,"\tjournal ID ASTM:    %s\n",Framework[system].CitationInformation[CurrentFramework][i].CitationJournalID_ASTM);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationJournalID_ISSN)>0)
           fprintf(FilePtr,"\tjournal ID ISSN:    %s\n",Framework[system].CitationInformation[CurrentFramework][i].CitationJournalID_ISSN);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationBookTitle)>0)
           fprintf(FilePtr,"\tbook title:         '%s'\n",Framework[system].CitationInformation[CurrentFramework][i].CitationBookTitle);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationBookPublisher)>0)
           fprintf(FilePtr,"\tbook publisher:     '%s'\n",Framework[system].CitationInformation[CurrentFramework][i].CitationBookPublisher);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationBookID_ISBN)>0)
           fprintf(FilePtr,"\tbook ID ISBN:       %s\n",Framework[system].CitationInformation[CurrentFramework][i].CitationBookID_ISBN);
         if(strlen(Framework[system].CitationInformation[CurrentFramework][i].CitationSpecialDetails)>0)
           fprintf(FilePtr,"\tspecial details:    '%s'\n",Framework[system].CitationInformation[CurrentFramework][i].CitationSpecialDetails);
       }
    }
    fprintf(FilePtr,"Space group: %d\n",SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].Number);
    fprintf(FilePtr,"\tIdentifier: %d\n",Framework[system].SpaceGroupIdentifier[CurrentFramework]);
    fprintf(FilePtr,"\tshort international Hermann-Mauguin symbol: %s\n",
          SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].ShortInternationalHermannMauguinSpaceGroupSymbol);
    fprintf(FilePtr,"\tlong international Hermann-Mauguin symbol: %s\n",
          SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].LongInternationalHermannMauguinSpaceGroupSymbol);
    fprintf(FilePtr,"\tHall symbol: %s\n",
          SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].HallSpaceGroupSymbol);

    fprintf(FilePtr,"\tNumber of lattice translations: %d ",SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].NumberOfLatticeTranslations);
    if(SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].NumberOfLatticeTranslations>0)
    {
      fprintf(FilePtr,"[ ");
      for(i=0;i<SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].NumberOfLatticeTranslations;i++)
        fprintf(FilePtr,"(%g,%g,%g) ",SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].LatticeTranslation[i].x,
                                      SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].LatticeTranslation[i].y,
                                      SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].LatticeTranslation[i].z);
      fprintf(FilePtr,"]");
    }
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"\tacentric/centric: %s\n",(SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].Centric)?"centric":"acentric");
    fprintf(FilePtr,"\tchiral: %s\n",(SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].Chiral)?"yes":"no");
    fprintf(FilePtr,"\tenantiomorphic: %s\n",(SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].Enantiomorphic)?"yes":"no");
    fprintf(FilePtr,"\tnumber of operators: %d\n",SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].NumberOfOperators);

    for(j=0;j<SpaceGroupData[Framework[system].SpaceGroupIdentifier[CurrentFramework]].NumberOfOperators;j++)
      fprintf(FilePtr,"\t\t'%s'\n",StringSpaceGroupElements[Framework[system].SpaceGroupIdentifier[CurrentFramework]][j]);

    if(Framework[system].FrameworkModels[CurrentFramework]==FLEXIBLE)
      fprintf(FilePtr,"Framework is simulated as 'flexible'\n");
    else
      fprintf(FilePtr,"Framework is simulated as 'rigid'\n");

    fprintf(FilePtr,"Shift: %lf %lf %lf\n",
                       Framework[system].ShiftUnitCell[CurrentFramework].x,
                       Framework[system].ShiftUnitCell[CurrentFramework].y,
                       Framework[system].ShiftUnitCell[CurrentFramework].z);
    fprintf(FilePtr,"Number of framework atoms: %d\n",Framework[system].NumberOfAtoms[CurrentFramework]);
    fprintf(FilePtr,"Number of asymmetric atoms: %d\n",Framework[system].NumberOfAsymmetricAtoms[CurrentFramework]);
    fprintf(FilePtr,"Number of free framework atoms: %d\n",Framework[system].NumberOfFreeAtoms[CurrentFramework]);
    fprintf(FilePtr,"Number of fixed framework atoms: %d\n",Framework[system].NumberOfFixedAtoms[CurrentFramework]);
    fprintf(FilePtr,"Number of framework atoms in the unit cell: %d\n",Framework[system].NumberOfUnitCellAtoms[CurrentFramework]);
    fprintf(FilePtr,"Framework Mass: %18.12lf [g/mol]\n",(double)Framework[system].FrameworkMassPerComponent[CurrentFramework]);
    fprintf(FilePtr,"Framework Density: %18.12lf [kg/m^3]\n",(double)(1e-3*
       Framework[system].FrameworkDensityPerComponent[CurrentFramework]));

    charge=0.0;
    smallest_charge=DBL_MAX;
    largest_charge=-DBL_MAX;
    for(i=0;i<Framework[system].NumberOfAtoms[CurrentFramework];i++)
    {
      Type=Framework[system].Atoms[CurrentFramework][i].Type;
      charge+=Framework[system].Atoms[CurrentFramework][i].Charge;
      if(Framework[system].Atoms[CurrentFramework][i].Charge>largest_charge) largest_charge=Framework[system].Atoms[CurrentFramework][i].Charge;
      if(Framework[system].Atoms[CurrentFramework][i].Charge<smallest_charge) smallest_charge=Framework[system].Atoms[CurrentFramework][i].Charge;
    }

    if(CorrectNetChargeOnPseudoAtom>=0)
      fprintf(FilePtr,"Charge corrected on PseudoAtom %d [%s] by subtracting: %lf\n",
        CorrectNetChargeOnPseudoAtom,PseudoAtoms[CorrectNetChargeOnPseudoAtom].Name,(double)CorrectNetChargeOnPseudoAtomValue);

    fprintf(FilePtr,"Framework has net charge: %lf\n",(double)charge);
    fprintf(FilePtr,"         largest charge : %lf\n",(double)largest_charge);
    fprintf(FilePtr,"         smallest charge: %lf\n",(double)smallest_charge);
    fprintf(FilePtr,"\n");
  }

  if(Framework[system].NumberOfAsymmetricIons>0)
  {
    fprintf(FilePtr,"Number of crystallographi (ion-) site types in the asymmetric unit cell: %d\n",Framework[system].NumberOfAsymmetricIons);
    for(i=0;i<Framework[system].NumberOfAsymmetricIons;i++)
      fprintf(FilePtr,"\t%4d asymmetric fractional position: %f %f %f   number of sites in full unit cell: %d\n",
          i,
          Framework[system].IonsAsymmetric[i].Position.x,
          Framework[system].IonsAsymmetric[i].Position.y,
          Framework[system].IonsAsymmetric[i].Position.z,
          crystallographic_stats[system].NumberOfCationSites[i]);
    fprintf(FilePtr,"\n");

    fprintf(FilePtr,"\tTotal number of crystallographic (ion-) sites in the unit cell: %d\n",Framework[system].NumberOfIons);
    for(i=0;i<Framework[system].NumberOfIons;i++)
    {
      fprintf(FilePtr,"\t%4d position: %12.8f %12.8f %12.8f type: %d\n",
          i,
          Framework[system].Ions[i].Position.x,
          Framework[system].Ions[i].Position.y,
          Framework[system].Ions[i].Position.z,
          Framework[system].Ions[i].AssymetricType);
    }
    fprintf(FilePtr,"\n");
  }

  switch(Framework[system].FrameworkModel)
  {
    case FULL:
      fprintf(FilePtr,"Using FULL Host-guest interaction calculation (for testing purposes)\n");
      break;
    case RIGID:
      fprintf(FilePtr,"Using RIGID Host-guest interactions\n");
      break;
    case GRID:
      fprintf(FilePtr,"Using an interpolation grid, testing the grid for each pseudo atom:\n");

      // make sure the order of random-numbers is the same (used when restarting from a binary restart)
      StoreRandomNumberStatus();
      CurrentSystem=system;
      TestForceGrid(FilePtr);
      RetrieveRandomNumberStatus();
      fprintf(FilePtr,"\n\n");

      if(NumberOfGridSeeds>0)
      {
        fprintf(FilePtr,"Seeds to block the VDW-grids\n");
        fprintf(FilePtr,"===========================================================================\n");
        for(i=0;i<NumberOfGridSeeds;i++)
          fprintf(FilePtr,"Grid seed (%d): %g %g %g\n",i,GridSeeds[i].x,GridSeeds[i].y,GridSeeds[i].z);
        fprintf(FilePtr,"\n\n");
      }
      break;
  }


  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Current Atom Status\n");
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"Number of framework atoms        : %d\n",Framework[system].TotalNumberOfAtoms);
  fprintf(FilePtr,"Number of cations molecules      : %d\n",NumberOfCationMolecules[system]);
  fprintf(FilePtr,"Number of adsorbate molecules    : %d\n",NumberOfAdsorbateMolecules[system]);

  for(i=0;i<NumberOfComponents;i++)
    fprintf(FilePtr,"Component %4d : %4d molecules\n",i,Components[i].NumberOfMolecules[system]);

  for(i=0;i<NumberOfPseudoAtoms;i++)
    fprintf(FilePtr,"Pseudo Atoms %4d [%8s]: %4d atoms\n",
       i,
       PseudoAtoms[i].Name,
       NumberOfPseudoAtomsType[system][i]);

  if(ComputePowderDiffractionPattern)
  {
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Powder diffraction options\n");
    fprintf(FilePtr,"===========================================================================\n");

    switch(Diffraction.Type)
    {
      case XRAY_DIFFRACTION:
        fprintf(FilePtr,"X-ray Diffraction\n");
        switch(Diffraction.RadiationType)
        {
          case CHROMIUM_RADIATION:
            fprintf(FilePtr,"Radiation type: Chromium\n");
            break;
          case IRON_RADIATION:
            fprintf(FilePtr,"Radiation type: Iron\n");
            break;
          case COPPER_RADIATION:
            fprintf(FilePtr,"Radiation type: Copper\n");
            break;
          case MOLYBDENUM_RADIATION:
            fprintf(FilePtr,"Radiation type: Molybdenum\n");
            break;
          case SILVER_RADIATION:
            fprintf(FilePtr,"Radiation type: Silver\n");
            break;
          case SYNCHROTRON_RADIATION:
            fprintf(FilePtr,"Radiation type: Synchrotron\n");
            break;
        }
        switch(Diffraction.lambda_type)
        {
          case DIFFRACTION_SINGLE:
            fprintf(FilePtr,"single wavelength: %f\n",Diffraction.lambda);
            break;
          case DIFFRACTION_DOUBLET:
            fprintf(FilePtr,"wavelength doublet: %f and %f\n",Diffraction.lambda,Diffraction.lambda2);
            break;
        }
        break;
      case ELECTRON_DIFFRACTION:
        fprintf(FilePtr,"Electron diffraction\n");
        fprintf(FilePtr,"Wavelength: %f\n",Diffraction.lambda);
        break;
      case NEUTRON_DIFFRACTION:
        fprintf(FilePtr,"Neutron diffraction\n");
        fprintf(FilePtr,"Wavelength: %f\n",Diffraction.lambda);
        break;
    }
    fprintf(FilePtr,"Spectrum obtained from 2theta=%f to %f using steps of %f (%d points)\n",
      Diffraction.two_theta_min,Diffraction.two_theta_max,Diffraction.two_theta_step,Diffraction.n);
    switch(Diffraction.PeakShape)
    {
      case DIFFRACTION_GAUSSIAN:
        fprintf(FilePtr,"Peaks are assumned to have Gaussian shapes\n");
        break;
      case DIFFRACTION_LORENTZIAN:
        fprintf(FilePtr,"Peaks are assumned to have Lorentzian shapes\n");
        break;
      case DIFFRACTION_PSEUDO_VOIGT:
        fprintf(FilePtr,"Peaks are assumned to have Pseudo-Voigt shapes\n");
        break;
    }
    fprintf(FilePtr,"Shape parameters: %f %f %f\n",Diffraction.u,Diffraction.v,Diffraction.w);
  }



  CurrentSystem=system;
  if(!ContinueAfterCrash)
  {
    InitializesEnergiesCurrentSystem();
    CalculateEnergy();
  }
  PrintEnergyStatus(FilePtr,"initial full energy");

  fprintf(FilePtr,"\n\n\n\n\n\n");
  fprintf(FilePtr,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  fprintf(FilePtr,"Starting simulation\n");
  fprintf(FilePtr,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  fprintf(FilePtr,"\n");
  fflush(FilePtr);
}

void PrintPostSimulationStatus(void)
{
  FILE *FilePtr;
  time_t curtime;
  struct tm *loctime;
  char buffer[256];

  /* Get the current time.  */
  curtime = time (NULL);

  /* Convert it to local time representation.  */
  loctime = localtime (&curtime);



  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    FilePtr=OutputFilePtr[CurrentSystem];

    fprintf(FilePtr,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    fprintf(FilePtr,"Finishing simulation\n");
    fprintf(FilePtr,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    fprintf(FilePtr,"\n\n");

    PrintEnergyStatus(FilePtr,"running energy");

    fprintf(FilePtr,"Monte-Carlo moves statistics\n");
    fprintf(FilePtr,"===========================================================================\n");
    fprintf(FilePtr,"\n");
    PrintSmallMCAddStatistics(FilePtr);
    PrintTranslationStatistics(FilePtr);
    PrintRandomTranslationStatistics(FilePtr);
    PrintRotationStatistics(FilePtr);
    PrintRandomRotationStatistics(FilePtr);
    PrintSwapAddStatistics(FilePtr);
    PrintSwapRemoveStatistics(FilePtr);
    PrintReinsertionStatistics(FilePtr);
    PrintReinsertionInPlaneStatistics(FilePtr);
    PrintReinsertionInPlaceStatistics(FilePtr);
    PrintPartialReinsertionStatistics(FilePtr);
    PrintIdentityChangeStatistics(FilePtr);
    PrintParallelTemperingStatistics(FilePtr);
    PrintHyperParallelTemperingStatistics(FilePtr);
    PrintParallelMolFractionStatistics(FilePtr);
    PrintChiralInversionStatistics(FilePtr);
    PrintVolumeChangeStatistics(FilePtr);
    PrintBoxShapeChangeStatistics(FilePtr);
    PrintFrameworkStatistics(FilePtr);
    PrintFrameworkShiftStatistics(FilePtr);
    PrintHybridNVEStatistics(FilePtr);
    PrintHybridNPHStatistics(FilePtr);
    PrintHybridNPHPRStatistics(FilePtr);
    PrintGibbsVolumeChangeStatistics(FilePtr);
    PrintGibbsSwapStatistics(FilePtr);
    PrintGibbsIdentityChangeStatistics(FilePtr);
    PrintCFSwapLambdaStatistics(FilePtr);
    PrintCBCFSwapLambdaStatistics(FilePtr);
    PrintCFGibbsLambdaStatistics(FilePtr);
    PrintCBCFGibbsLambdaStatistics(FilePtr);
    PrintRXMCStatistics(FilePtr);
    PrintExchangeFractionalParticleStatistics(FilePtr);
    fprintf(FilePtr,"\n\n");
    PrintCPUStatistics(FilePtr);
    fprintf(FilePtr,"\n\n");

    UHostHostRunning=UHostHost[CurrentSystem];
    UHostHostVDWRunning=UHostHostVDW[CurrentSystem];
    UHostHostCoulombRunning=UHostHostCoulomb[CurrentSystem];
    UHostHostChargeChargeRealRunning=UHostHostChargeChargeReal[CurrentSystem];
    UHostHostChargeChargeFourierRunning=UHostHostChargeChargeFourier[CurrentSystem];
    UHostHostChargeBondDipoleRealRunning=UHostHostChargeBondDipoleReal[CurrentSystem];
    UHostHostChargeBondDipoleFourierRunning=UHostHostChargeBondDipoleFourier[CurrentSystem];
    UHostHostBondDipoleBondDipoleRealRunning=UHostHostBondDipoleBondDipoleReal[CurrentSystem];
    UHostHostBondDipoleBondDipoleFourierRunning=UHostHostBondDipoleBondDipoleFourier[CurrentSystem];

    UHostAdsorbateRunning=UHostAdsorbate[CurrentSystem];
    UHostAdsorbateVDWRunning=UHostAdsorbateVDW[CurrentSystem];
    UHostAdsorbateCoulombRunning=UHostAdsorbateCoulomb[CurrentSystem];
    UHostAdsorbateChargeChargeRealRunning=UHostAdsorbateChargeChargeReal[CurrentSystem];
    UHostAdsorbateChargeChargeFourierRunning=UHostAdsorbateChargeChargeFourier[CurrentSystem];
    UHostAdsorbateChargeBondDipoleRealRunning=UHostAdsorbateChargeBondDipoleReal[CurrentSystem];
    UHostAdsorbateChargeBondDipoleFourierRunning=UHostAdsorbateChargeBondDipoleFourier[CurrentSystem];
    UHostAdsorbateBondDipoleBondDipoleRealRunning=UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
    UHostAdsorbateBondDipoleBondDipoleFourierRunning=UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];

    UHostCationRunning=UHostCation[CurrentSystem];
    UHostCationVDWRunning=UHostCationVDW[CurrentSystem];
    UHostCationCoulombRunning=UHostCationCoulomb[CurrentSystem];
    UHostCationChargeChargeRealRunning=UHostCationChargeChargeReal[CurrentSystem];
    UHostCationChargeChargeFourierRunning=UHostCationChargeChargeFourier[CurrentSystem];
    UHostCationChargeBondDipoleRealRunning=UHostCationChargeBondDipoleReal[CurrentSystem];
    UHostCationChargeBondDipoleFourierRunning=UHostCationChargeBondDipoleFourier[CurrentSystem];
    UHostCationBondDipoleBondDipoleRealRunning=UHostCationBondDipoleBondDipoleReal[CurrentSystem];
    UHostCationBondDipoleBondDipoleFourierRunning=UHostCationBondDipoleBondDipoleFourier[CurrentSystem];

    UAdsorbateAdsorbateRunning=UAdsorbateAdsorbate[CurrentSystem];
    UAdsorbateAdsorbateVDWRunning=UAdsorbateAdsorbateVDW[CurrentSystem];
    UAdsorbateAdsorbateCoulombRunning=UAdsorbateAdsorbateCoulomb[CurrentSystem];
    UAdsorbateAdsorbateChargeChargeRealRunning=UAdsorbateAdsorbateChargeChargeReal[CurrentSystem];
    UAdsorbateAdsorbateChargeChargeFourierRunning=UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem];
    UAdsorbateAdsorbateChargeBondDipoleRealRunning=UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem];
    UAdsorbateAdsorbateChargeBondDipoleFourierRunning=UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem];
    UAdsorbateAdsorbateBondDipoleBondDipoleRealRunning=UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
    UAdsorbateAdsorbateBondDipoleBondDipoleFourierRunning=UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];

    UCationCationRunning=UCationCation[CurrentSystem];
    UCationCationVDWRunning=UCationCationVDW[CurrentSystem];
    UCationCationCoulombRunning=UCationCationCoulomb[CurrentSystem];
    UCationCationChargeChargeRealRunning=UCationCationChargeChargeReal[CurrentSystem];
    UCationCationChargeChargeFourierRunning=UCationCationChargeChargeFourier[CurrentSystem];
    UCationCationChargeBondDipoleRealRunning=UCationCationChargeBondDipoleReal[CurrentSystem];
    UCationCationChargeBondDipoleFourierRunning=UCationCationChargeBondDipoleFourier[CurrentSystem];
    UCationCationBondDipoleBondDipoleRealRunning=UCationCationBondDipoleBondDipoleReal[CurrentSystem];
    UCationCationBondDipoleBondDipoleFourierRunning=UCationCationBondDipoleBondDipoleFourier[CurrentSystem];

    UAdsorbateCationRunning=UAdsorbateCation[CurrentSystem];
    UAdsorbateCationVDWRunning=UAdsorbateCationVDW[CurrentSystem];
    UAdsorbateCationCoulombRunning=UAdsorbateCationCoulomb[CurrentSystem];
    UAdsorbateCationChargeChargeRealRunning=UAdsorbateCationChargeChargeReal[CurrentSystem];
    UAdsorbateCationChargeChargeFourierRunning=UAdsorbateCationChargeChargeFourier[CurrentSystem];
    UAdsorbateCationChargeBondDipoleRealRunning=UAdsorbateCationChargeBondDipoleReal[CurrentSystem];
    UAdsorbateCationChargeBondDipoleFourierRunning=UAdsorbateCationChargeBondDipoleFourier[CurrentSystem];
    UAdsorbateCationBondDipoleBondDipoleRealRunning=UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem];
    UAdsorbateCationBondDipoleBondDipoleFourierRunning=UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem];

    UHostBondRunning=UHostBond[CurrentSystem];
    UHostUreyBradleyRunning=UHostUreyBradley[CurrentSystem];
    UHostBendRunning=UHostBend[CurrentSystem];
    UHostInversionBendRunning=UHostInversionBend[CurrentSystem];
    UHostTorsionRunning=UHostTorsion[CurrentSystem];
    UHostImproperTorsionRunning=UHostImproperTorsion[CurrentSystem];
    UHostOutOfPlaneRunning=UHostOutOfPlane[CurrentSystem];
    UHostBondBondRunning=UHostBondBond[CurrentSystem];
    UHostBondBendRunning=UHostBondBend[CurrentSystem];
    UHostBendBendRunning=UHostBendBend[CurrentSystem];
    UHostBondTorsionRunning=UHostBondTorsion[CurrentSystem];
    UHostBendTorsionRunning=UHostBendTorsion[CurrentSystem];

    UAdsorbateBondRunning=UAdsorbateBond[CurrentSystem];
    UAdsorbateUreyBradleyRunning=UAdsorbateUreyBradley[CurrentSystem];
    UAdsorbateBendRunning=UAdsorbateBend[CurrentSystem];
    UAdsorbateInversionBendRunning=UAdsorbateInversionBend[CurrentSystem];
    UAdsorbateTorsionRunning=UAdsorbateTorsion[CurrentSystem];
    UAdsorbateImproperTorsionRunning=UAdsorbateImproperTorsion[CurrentSystem];
    UAdsorbateOutOfPlaneRunning=UAdsorbateOutOfPlane[CurrentSystem];
    UAdsorbateBondBondRunning=UAdsorbateBondBond[CurrentSystem];
    UAdsorbateBondBendRunning=UAdsorbateBondBend[CurrentSystem];
    UAdsorbateBendBendRunning=UAdsorbateBendBend[CurrentSystem];
    UAdsorbateBondTorsionRunning=UAdsorbateBondTorsion[CurrentSystem];
    UAdsorbateBendTorsionRunning=UAdsorbateBendTorsion[CurrentSystem];
    UAdsorbateIntraVDWRunning=UAdsorbateIntraVDW[CurrentSystem];
    UAdsorbateIntraChargeChargeCoulombRunning=UAdsorbateIntraChargeCharge[CurrentSystem];
    UAdsorbateIntraChargeBondDipoleCoulombRunning=UAdsorbateIntraChargeBondDipole[CurrentSystem];
    UAdsorbateIntraBondDipoleBondDipoleCoulombRunning=UAdsorbateIntraBondDipoleBondDipole[CurrentSystem];

    UCationBondRunning=UCationBond[CurrentSystem];
    UCationUreyBradleyRunning=UCationUreyBradley[CurrentSystem];
    UCationBendRunning=UCationBend[CurrentSystem];
    UCationInversionBendRunning=UCationInversionBend[CurrentSystem];
    UCationTorsionRunning=UCationTorsion[CurrentSystem];
    UCationImproperTorsionRunning=UCationImproperTorsion[CurrentSystem];
    UCationOutOfPlaneRunning=UCationOutOfPlane[CurrentSystem];
    UCationBondBondRunning=UCationBondBond[CurrentSystem];
    UCationBondBendRunning=UCationBondBend[CurrentSystem];
    UCationBendBendRunning=UCationBendBend[CurrentSystem];
    UCationBondTorsionRunning=UCationBondTorsion[CurrentSystem];
    UCationBendTorsionRunning=UCationBendTorsion[CurrentSystem];
    UCationIntraVDWRunning=UCationIntraVDW[CurrentSystem];
    UCationIntraChargeChargeCoulombRunning=UCationIntraChargeCharge[CurrentSystem];
    UCationIntraChargeBondDipoleCoulombRunning=UCationIntraChargeBondDipole[CurrentSystem];
    UCationIntraBondDipoleBondDipoleCoulombRunning=UCationIntraBondDipoleBondDipole[CurrentSystem];

    UTailCorrectionRunning=UTailCorrection[CurrentSystem];

    UDistanceConstraintsRunning=UDistanceConstraints[CurrentSystem];
    UAngleConstraintsRunning=UAngleConstraints[CurrentSystem];
    UDihedralConstraintsRunning=UDihedralConstraints[CurrentSystem];
    UInversionBendConstraintsRunning=UInversionBendConstraints[CurrentSystem];
    UOutOfPlaneDistanceConstraintsRunning=UOutOfPlaneDistanceConstraints[CurrentSystem];
    UExclusionConstraintsRunning=UExclusionConstraints[CurrentSystem];

    UHostPolarizationRunning=UHostPolarization[CurrentSystem];
    UAdsorbatePolarizationRunning=UAdsorbatePolarization[CurrentSystem];
    UCationPolarizationRunning=UCationPolarization[CurrentSystem];

    UHostBackPolarizationRunning=UHostBackPolarization[CurrentSystem];
    UAdsorbateBackPolarizationRunning=UAdsorbateBackPolarization[CurrentSystem];
    UCationBackPolarizationRunning=UCationBackPolarization[CurrentSystem];

    UTotalRunning=UTotal[CurrentSystem];

    InitializesEnergiesCurrentSystem();


    // switch off cell-list for a 'really' full computation of the energy
    UseCellLists[CurrentSystem]=FALSE;

    InitializeForces();
    CalculateForce();

    PrintEnergyStatus(FilePtr,"full final energy");
    PrintEnergyDriftStatus(FilePtr);

    PrintAverageTotalSystemEnergiesMC(FilePtr);

    fprintf(FilePtr,"\nSimulation finished,  %d warnings\n",NumberOfWarnings[CurrentSystem]);
    PrintWarningStatus();
    fprintf(FilePtr,"\n\n");

    /* Print out the date and time in the standard format.  */
    fprintf(FilePtr,"%s",asctime(loctime));

    /* Print it out in a nice format.  */
    strftime (buffer, 256, "Simulation finished on %A, %B %d.", loctime);
    fprintf(FilePtr,"%s\n",buffer);
    strftime (buffer, 256, "The end time was %I:%M %p.", loctime);
    fprintf(FilePtr,"%s\n\n",buffer);

  }
}
void PrintEnergyStatusToStdOut()
{
  printf("\n");
  printf("\n");
  printf("Current Energy Status\n");
  printf("===========================================================================\n");
  printf("\n");

  printf("Internal energy:\n");
  printf("Host stretch energy:                                    %18.8lf\n",(double)UHostBond[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host UreyBradley energy:                                %18.8lf\n",(double)UHostUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host bend energy:                                       %18.8lf\n",(double)UHostBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host inversion-bend energy:                             %18.8lf\n",(double)UHostInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host torsion energy:                                    %18.8lf\n",(double)UHostTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host improper torsion energy:                           %18.8lf\n",(double)UHostImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host out-of-plane energy:                               %18.8lf\n",(double)UHostOutOfPlane[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host stretch/stretch energy:                            %18.8lf\n",(double)UHostBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host bend/bend energy:                                  %18.8lf\n",(double)UHostBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host stretch/bend energy:                               %18.8lf\n",(double)UHostBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host stretch/torsion energy:                            %18.8lf\n",(double)UHostBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Host bend/torsion energy:                               %18.8lf\n",(double)UHostBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");

  printf("Adsorbate stretch energy:                               %18.8lf\n",(double)UAdsorbateBond[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate UreyBradley energy:                           %18.8lf\n",(double)UAdsorbateUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate bend energy:                                  %18.8lf\n",(double)UAdsorbateBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate inversion-bend energy:                        %18.8lf\n",(double)UAdsorbateInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate torsion energy:                               %18.8lf\n",(double)UAdsorbateTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate improper torsion energy:                      %18.8lf\n",(double)UAdsorbateImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate out-of-plane energy:                          %18.8lf\n",(double)UAdsorbateOutOfPlane[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate stretch/stretch energy:                       %18.8lf\n",(double)UAdsorbateBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate bend/bend energy:                             %18.8lf\n",(double)UAdsorbateBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate stretch/bend energy:                          %18.8lf\n",(double)UAdsorbateBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate stretch/torsion energy:                       %18.8lf\n",(double)UAdsorbateBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate bend/torsion energy:                          %18.8lf\n",(double)UAdsorbateBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate intra VDW energy:                             %18.8lf\n",(double)UAdsorbateIntraVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate intra charge-charge Coulomb energy:           %18.8lf\n",(double)UAdsorbateIntraChargeCharge[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate intra charge-bonddipole Coulomb energy:       %18.8lf\n",(double)UAdsorbateIntraChargeBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Adsorbate intra bonddipole-bonddipole Coulomb energy:   %18.8lf\n",(double)UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");

  printf("Cation stretch energy:                                  %18.8lf\n",(double)UCationBond[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation UreyBradley energy:                              %18.8lf\n",(double)UCationUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation bend energy:                                     %18.8lf\n",(double)UCationBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation inversion-bend energy:                           %18.8lf\n",(double)UCationInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation torsion energy:                                  %18.8lf\n",(double)UCationTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation improper torsion energy:                         %18.8lf\n",(double)UCationImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation out-of-plane energy:                             %18.8lf\n",(double)UCationOutOfPlane[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation stretch/stretch energy:                          %18.8lf\n",(double)UCationBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation bend/bend energy:                                %18.8lf\n",(double)UCationBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation stretch/bend energy:                             %18.8lf\n",(double)UCationBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation stretch/torsion energy:                          %18.8lf\n",(double)UCationBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation bend/torsion energy:                             %18.8lf\n",(double)UCationBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation intra VDW energy:                                %18.8lf\n",(double)UCationIntraVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation intra charge-charge Coulomb energy:              %18.8lf\n",(double)UCationIntraChargeCharge[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation intra charge-bonddipole Coulomb energy:          %18.8lf\n",(double)UCationIntraChargeBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Cation intra bonddipole-bonddipole Coulomb energy:      %18.8lf\n",(double)UCationIntraBondDipoleBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");

  printf("Host/Host energy:                                     %18.8lf\n",(double)UHostHost[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Host VDW energy:                                 %18.8lf\n",(double)UHostHostVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Host Coulomb energy:                             %18.8lf\n",(double)UHostHostCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Host charge-charge Real energy:                  %18.8lf\n",(double)UHostHostChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Host charge-charge Fourier energy:               %18.8lf\n",(double)UHostHostChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Host charge-bonddipole Real energy:              %18.8lf\n",(double)UHostHostChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Host charge-bonddipole Fourier energy:           %18.8lf\n",(double)UHostHostChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Host bondipole-bonddipole Real energy:           %18.8lf\n",(double)UHostHostBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Host bondipole-bonddipole Fourier energy:        %18.8lf\n",(double)UHostHostBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");

  printf("Host/Adsorbate energy:                                %18.8lf\n",(double)UHostAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Adsorbate VDW energy:                            %18.8lf\n",(double)UHostAdsorbateVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Adsorbate Coulomb energy:                        %18.8lf\n",(double)UHostAdsorbateCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Adsorbate charge-charge Real energy:             %18.8lf\n",(double)UHostAdsorbateChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Adsorbate charge-charge Fourier energy:          %18.8lf\n",(double)UHostAdsorbateChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Adsorbate charge-bonddipole Real energy:         %18.8lf\n",(double)UHostAdsorbateChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Adsorbate charge-bonddipole Fourier energy:      %18.8lf\n",(double)UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Adsorbate bondipole-bonddipole Real energy:      %18.8lf\n",(double)UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Adsorbate bondipole-bonddipole Fourier energy:   %18.8lf\n",(double)UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");

  printf("Host/Cation energy:                                   %18.8lf\n",(double)UHostCation[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Cation VDW energy:                               %18.8lf\n",(double)UHostCationVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Cation Coulomb energy:                           %18.8lf\n",(double)UHostCationCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Cation charge-charge Real energy:                %18.8lf\n",(double)UHostCationChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Cation charge-charge Fourier energy:             %18.8lf\n",(double)UHostCationChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Cation charge-bonddipole Real energy:            %18.8lf\n",(double)UHostCationChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Cation charge-bonddipole Fourier energy:         %18.8lf\n",(double)UHostCationChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Cation bondipole-bonddipole Real energy:         %18.8lf\n",(double)UHostCationBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost/Cation bondipole-bonddipole Fourier energy:      %18.8lf\n",(double)UHostCationBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");

  printf("Adsorbate/Adsorbate energy:                                %18.8lf\n",(double)UAdsorbateAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Adsorbate VDW energy:                            %18.8lf\n",(double)UAdsorbateAdsorbateVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Adsorbate Coulomb energy:                        %18.8lf\n",(double)UAdsorbateAdsorbateCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Adsorbate charge-charge Real energy:             %18.8lf\n",(double)UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Adsorbate charge-charge Fourier energy:          %18.8lf\n",(double)UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Adsorbate charge-bonddipole Real energy:         %18.8lf\n",(double)UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Adsorbate charge-bonddipole Fourier energy:      %18.8lf\n",(double)UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Adsorbate bondipole-bonddipole Real energy:      %18.8lf\n",(double)UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Adsorbate bondipole-bonddipole Fourier energy:   %18.8lf\n",(double)UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");

  printf("Adsorbate/Cation energy:                                   %18.8lf\n",(double)UAdsorbateCation[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Cation VDW energy:                               %18.8lf\n",(double)UAdsorbateCationVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Cation Coulomb energy:                           %18.8lf\n",(double)UAdsorbateCationCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Cation charge-charge Real energy:                %18.8lf\n",(double)UAdsorbateCationChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Cation charge-charge Fourier energy:             %18.8lf\n",(double)UAdsorbateCationChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Cation charge-bonddipole Real energy:            %18.8lf\n",(double)UAdsorbateCationChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Cation charge-bonddipole Fourier energy:         %18.8lf\n",(double)UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Cation bondipole-bonddipole Real energy:         %18.8lf\n",(double)UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate/Cation bondipole-bonddipole Fourier energy:      %18.8lf\n",(double)UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");

  printf("Cation/Cation energy:                                   %18.8lf\n",(double)UCationCation[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tCation/Cation VDW energy:                               %18.8lf\n",(double)UCationCationVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tCation/Cation Coulomb energy:                           %18.8lf\n",(double)UCationCationCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tCation/Cation charge-charge Real energy:                %18.8lf\n",(double)UCationCationChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tCation/Cation charge-charge Fourier energy:             %18.8lf\n",(double)UCationCationChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tCation/Cation charge-bonddipole Real energy:            %18.8lf\n",(double)UCationCationChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tCation/Cation charge-bonddipole Fourier energy:         %18.8lf\n",(double)UCationCationChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tCation/Cation bondipole-bonddipole Real energy:         %18.8lf\n",(double)UCationCationBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tCation/Cation bondipole-bonddipole Fourier energy:      %18.8lf\n",(double)UCationCationBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");


  printf("Polarization energy:\n");
  printf("\tHost polarization energy:                             %18.8lf\n",(double)UHostPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate polarization energy:                        %18.8lf\n",(double)UAdsorbatePolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tCation polarization energy:                           %18.8lf\n",(double)UCationPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tHost back-polarization energy:                             %18.8lf\n",(double)UHostBackPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tAdsorbate back-polarization energy:                        %18.8lf\n",(double)UAdsorbateBackPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tCation back-polarization energy:                           %18.8lf\n",(double)UCationBackPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");

  printf("Tail-correction energy:                               %18.8lf\n",(double)UTailCorrection[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");
  printf("Distance constraints energy:                          %18.8lf\n",(double)UDistanceConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Angle constraints energy:                             %18.8lf\n",(double)UAngleConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Dihedral constraints energy:                          %18.8lf\n",(double)UDihedralConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Inversion-bend constraints energy:                    %18.8lf\n",(double)UInversionBendConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Out-of-plane distance constraints energy:             %18.8lf\n",(double)UOutOfPlaneDistanceConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("Exclusion constraints energy:                         %18.8lf\n",(double)UExclusionConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\n");

  printf("===================================================================\n");
  printf("Total energy: %18.12lf\n",(double)UTotal[CurrentSystem]*ENERGY_TO_KELVIN);
  printf("\tTotal Van der Waals: %lf\n",(double)
     (UHostHostVDW[CurrentSystem]+UHostAdsorbateVDW[CurrentSystem]+UHostCationVDW[CurrentSystem]+
      UAdsorbateAdsorbateVDW[CurrentSystem]+UCationCationVDW[CurrentSystem]+UAdsorbateCationVDW[CurrentSystem])*ENERGY_TO_KELVIN);
  printf("\tTotal Coulomb: %lf\n\n",(double)
     (UHostHostCoulomb[CurrentSystem]+UHostCationCoulomb[CurrentSystem]+UHostAdsorbateCoulomb[CurrentSystem]+
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+UCationCationCoulomb[CurrentSystem]+UAdsorbateCationCoulomb[CurrentSystem]+
      UAdsorbateIntraChargeCharge[CurrentSystem]+UAdsorbateIntraChargeBondDipole[CurrentSystem]+
      UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+UCationIntraChargeCharge[CurrentSystem]+
      UCationIntraChargeBondDipole[CurrentSystem]+UCationIntraBondDipoleBondDipole[CurrentSystem])*ENERGY_TO_KELVIN);
  printf("\tTotal Polarization: %lf\n\n",(double)(UHostPolarization[CurrentSystem]+
      UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem]+UHostBackPolarization[CurrentSystem]+
      UAdsorbateBackPolarization[CurrentSystem]+UCationBackPolarization[CurrentSystem])*ENERGY_TO_KELVIN);
  printf("\n\n\n");
}


void PrintEnergyStatus(FILE *FilePtr,char *string)
{
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Current (%s) Energy Status\n",string);
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Internal energy:\n");
  fprintf(FilePtr,"Host stretch energy:                                    %18.8lf\n",(double)UHostBond[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host UreyBradley energy:                                %18.8lf\n",(double)UHostUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host bend energy:                                       %18.8lf\n",(double)UHostBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host inversion-bend energy:                             %18.8lf\n",(double)UHostInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host torsion energy:                                    %18.8lf\n",(double)UHostTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host improper torsion energy:                           %18.8lf\n",(double)UHostImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host out-of-plane energy:                               %18.8lf\n",(double)UHostOutOfPlane[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host stretch/stretch energy:                            %18.8lf\n",(double)UHostBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host bend/bend energy:                                  %18.8lf\n",(double)UHostBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host stretch/bend energy:                               %18.8lf\n",(double)UHostBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host stretch/torsion energy:                            %18.8lf\n",(double)UHostBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host bend/torsion energy:                               %18.8lf\n",(double)UHostBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Adsorbate stretch energy:                               %18.8lf\n",(double)UAdsorbateBond[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate UreyBradley energy:                           %18.8lf\n",(double)UAdsorbateUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate bend energy:                                  %18.8lf\n",(double)UAdsorbateBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate inversion-bend energy:                        %18.8lf\n",(double)UAdsorbateInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate torsion energy:                               %18.8lf\n",(double)UAdsorbateTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate improper torsion energy:                      %18.8lf\n",(double)UAdsorbateImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate out-of-plane energy:                          %18.8lf\n",(double)UAdsorbateOutOfPlane[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate stretch/stretch energy:                       %18.8lf\n",(double)UAdsorbateBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate bend/bend energy:                             %18.8lf\n",(double)UAdsorbateBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate stretch/bend energy:                          %18.8lf\n",(double)UAdsorbateBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate stretch/torsion energy:                       %18.8lf\n",(double)UAdsorbateBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate bend/torsion energy:                          %18.8lf\n",(double)UAdsorbateBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate intra VDW energy:                             %18.8lf\n",(double)UAdsorbateIntraVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate intra charge-charge Coulomb energy:           %18.8lf\n",(double)UAdsorbateIntraChargeCharge[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate intra charge-bonddipole Coulomb energy:       %18.8lf\n",(double)UAdsorbateIntraChargeBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate intra bonddipole-bonddipole Coulomb energy:   %18.8lf\n",(double)UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Cation stretch energy:                                  %18.8lf\n",(double)UCationBond[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation UreyBradley energy:                              %18.8lf\n",(double)UCationUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation bend energy:                                     %18.8lf\n",(double)UCationBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation inversion-bend energy:                           %18.8lf\n",(double)UCationInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation torsion energy:                                  %18.8lf\n",(double)UCationTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation improper torsion energy:                         %18.8lf\n",(double)UCationImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation out-of-plane energy:                             %18.8lf\n",(double)UCationOutOfPlane[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation stretch/stretch energy:                          %18.8lf\n",(double)UCationBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation bend/bend energy:                                %18.8lf\n",(double)UCationBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation stretch/bend energy:                             %18.8lf\n",(double)UCationBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation stretch/torsion energy:                          %18.8lf\n",(double)UCationBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation bend/torsion energy:                             %18.8lf\n",(double)UCationBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation intra VDW energy:                                %18.8lf\n",(double)UCationIntraVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation intra charge-charge Coulomb energy:              %18.8lf\n",(double)UCationIntraChargeCharge[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation intra charge-bonddipole Coulomb energy:          %18.8lf\n",(double)UCationIntraChargeBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation intra bonddipole-bonddipole Coulomb energy:      %18.8lf\n",(double)UCationIntraBondDipoleBondDipole[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Host/Host energy:                                     %18.8lf\n",(double)UHostHost[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Host VDW energy:                                 %18.8lf\n",(double)UHostHostVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Host Coulomb energy:                             %18.8lf\n",(double)UHostHostCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Host charge-charge Real energy:                  %18.8lf\n",(double)UHostHostChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Host charge-charge Fourier energy:               %18.8lf\n",(double)UHostHostChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Host charge-bonddipole Real energy:              %18.8lf\n",(double)UHostHostChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Host charge-bonddipole Fourier energy:           %18.8lf\n",(double)UHostHostChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Host bondipole-bonddipole Real energy:           %18.8lf\n",(double)UHostHostBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Host bondipole-bonddipole Fourier energy:        %18.8lf\n",(double)UHostHostBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Host/Adsorbate energy:                                %18.8lf\n",(double)UHostAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Adsorbate VDW energy:                            %18.8lf\n",(double)UHostAdsorbateVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Adsorbate Coulomb energy:                        %18.8lf\n",(double)UHostAdsorbateCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Adsorbate charge-charge Real energy:             %18.8lf\n",(double)UHostAdsorbateChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Adsorbate charge-charge Fourier energy:          %18.8lf\n",(double)UHostAdsorbateChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Adsorbate charge-bonddipole Real energy:         %18.8lf\n",(double)UHostAdsorbateChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Adsorbate charge-bonddipole Fourier energy:      %18.8lf\n",(double)UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Adsorbate bondipole-bonddipole Real energy:      %18.8lf\n",(double)UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Adsorbate bondipole-bonddipole Fourier energy:   %18.8lf\n",(double)UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Host/Cation energy:                                   %18.8lf\n",(double)UHostCation[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Cation VDW energy:                               %18.8lf\n",(double)UHostCationVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Cation Coulomb energy:                           %18.8lf\n",(double)UHostCationCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Cation charge-charge Real energy:                %18.8lf\n",(double)UHostCationChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Cation charge-charge Fourier energy:             %18.8lf\n",(double)UHostCationChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Cation charge-bonddipole Real energy:            %18.8lf\n",(double)UHostCationChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Cation charge-bonddipole Fourier energy:         %18.8lf\n",(double)UHostCationChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Cation bondipole-bonddipole Real energy:         %18.8lf\n",(double)UHostCationBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Cation bondipole-bonddipole Fourier energy:      %18.8lf\n",(double)UHostCationBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Adsorbate/Adsorbate energy:                                %18.8lf\n",(double)UAdsorbateAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Adsorbate VDW energy:                            %18.8lf\n",(double)UAdsorbateAdsorbateVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Adsorbate Coulomb energy:                        %18.8lf\n",(double)UAdsorbateAdsorbateCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Adsorbate charge-charge Real energy:             %18.8lf\n",(double)UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Adsorbate charge-charge Fourier energy:          %18.8lf\n",(double)UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Adsorbate charge-bonddipole Real energy:         %18.8lf\n",(double)UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Adsorbate charge-bonddipole Fourier energy:      %18.8lf\n",(double)UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Adsorbate bondipole-bonddipole Real energy:      %18.8lf\n",(double)UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Adsorbate bondipole-bonddipole Fourier energy:   %18.8lf\n",(double)UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Adsorbate/Cation energy:                                   %18.8lf\n",(double)UAdsorbateCation[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Cation VDW energy:                               %18.8lf\n",(double)UAdsorbateCationVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Cation Coulomb energy:                           %18.8lf\n",(double)UAdsorbateCationCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Cation charge-charge Real energy:                %18.8lf\n",(double)UAdsorbateCationChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Cation charge-charge Fourier energy:             %18.8lf\n",(double)UAdsorbateCationChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Cation charge-bonddipole Real energy:            %18.8lf\n",(double)UAdsorbateCationChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Cation charge-bonddipole Fourier energy:         %18.8lf\n",(double)UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Cation bondipole-bonddipole Real energy:         %18.8lf\n",(double)UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Cation bondipole-bonddipole Fourier energy:      %18.8lf\n",(double)UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Cation/Cation energy:                                   %18.8lf\n",(double)UCationCation[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation/Cation VDW energy:                               %18.8lf\n",(double)UCationCationVDW[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation/Cation Coulomb energy:                           %18.8lf\n",(double)UCationCationCoulomb[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation/Cation charge-charge Real energy:                %18.8lf\n",(double)UCationCationChargeChargeReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation/Cation charge-charge Fourier energy:             %18.8lf\n",(double)UCationCationChargeChargeFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation/Cation charge-bonddipole Real energy:            %18.8lf\n",(double)UCationCationChargeBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation/Cation charge-bonddipole Fourier energy:         %18.8lf\n",(double)UCationCationChargeBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation/Cation bondipole-bonddipole Real energy:         %18.8lf\n",(double)UCationCationBondDipoleBondDipoleReal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation/Cation bondipole-bonddipole Fourier energy:      %18.8lf\n",(double)UCationCationBondDipoleBondDipoleFourier[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");


  fprintf(FilePtr,"Polarization energy:\n");
  fprintf(FilePtr,"\tHost polarization energy:                             %18.8lf\n",(double)UHostPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate polarization energy:                        %18.8lf\n",(double)UAdsorbatePolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation polarization energy:                           %18.8lf\n",(double)UCationPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost back-polarization energy:                             %18.8lf\n",(double)UHostBackPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate back-polarization energy:                        %18.8lf\n",(double)UAdsorbateBackPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation back-polarization energy:                           %18.8lf\n",(double)UCationBackPolarization[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Tail-correction energy:                               %18.8lf\n",(double)UTailCorrection[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Distance constraints energy:                          %18.8lf\n",(double)UDistanceConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Angle constraints energy:                             %18.8lf\n",(double)UAngleConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Dihedral constraints energy:                          %18.8lf\n",(double)UDihedralConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Inversion-bend constraints energy:                    %18.8lf\n",(double)UInversionBendConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Out-of-plane distance constraints energy:             %18.8lf\n",(double)UOutOfPlaneDistanceConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Exclusion constraints energy:                         %18.8lf\n",(double)UExclusionConstraints[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"===================================================================\n");
  fprintf(FilePtr,"Total energy: %18.12lf\n",(double)UTotal[CurrentSystem]*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tTotal Van der Waals: %lf\n",(double)
      (UHostHostVDW[CurrentSystem]+UHostAdsorbateVDW[CurrentSystem]+UHostCationVDW[CurrentSystem]+
       UAdsorbateAdsorbateVDW[CurrentSystem]+UCationCationVDW[CurrentSystem]+UAdsorbateCationVDW[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tTotal Coulomb: %lf\n\n",(double)
      (UHostHostCoulomb[CurrentSystem]+UHostCationCoulomb[CurrentSystem]+UHostAdsorbateCoulomb[CurrentSystem]+
       UAdsorbateAdsorbateCoulomb[CurrentSystem]+UCationCationCoulomb[CurrentSystem]+UAdsorbateCationCoulomb[CurrentSystem]+
       UAdsorbateIntraChargeCharge[CurrentSystem]+UAdsorbateIntraChargeBondDipole[CurrentSystem]+
       UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+UCationIntraChargeCharge[CurrentSystem]+
       UCationIntraChargeBondDipole[CurrentSystem]+UCationIntraBondDipoleBondDipole[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tTotal Polarization: %lf\n\n",(double)(UHostPolarization[CurrentSystem]+
       UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem]+UHostBackPolarization[CurrentSystem]+
       UAdsorbateBackPolarization[CurrentSystem]+UCationBackPolarization[CurrentSystem])*ENERGY_TO_KELVIN);
}

void PrintEnergyDriftStatus(FILE *FilePtr)
{
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Energy-drift status\n");
  fprintf(FilePtr,"===========================================================================\n");
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Internal energy:\n");
  fprintf(FilePtr,"Host stretch energy-drift:                                           %lg\n",
                   (double)(UHostBondRunning-UHostBond[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host UreyBradley energy-drift:                                       %lg\n",
                   (double)(UHostUreyBradleyRunning-UHostUreyBradley[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host bend energy-drift:                                              %lg\n",
                   (double)(UHostBendRunning-UHostBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host inversion-bend energy-drift:                                    %lg\n",
                   (double)(UHostInversionBendRunning-UHostInversionBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host torsion energy-drift:                                           %lg\n",
                   (double)(UHostTorsionRunning-UHostTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host torsion improper energy-drift:                                  %lg\n",
                   (double)(UHostImproperTorsionRunning-UHostImproperTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host out-of-plane energy-drift:                                      %lg\n",
                   (double)(UHostOutOfPlaneRunning-UHostOutOfPlane[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host stretch/stretch energy-drift:                                   %lg\n",
                   (double)(UHostBondBondRunning-UHostBondBond[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host stretch/bend energy-drift:                                      %lg\n",
                   (double)(UHostBondBendRunning-UHostBondBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host bend/bend energy-drift:                                         %lg\n",
                   (double)(UHostBendBendRunning-UHostBendBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host stretch/torsion energy-drift:                                   %lg\n",
                   (double)(UHostBondTorsionRunning-UHostBondTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Host bend/torsion energy-drift:                                      %lg\n",
                   (double)(UHostBendTorsionRunning-UHostBendTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Adsorbate stretch energy-drift:                                      %lg\n",
                   (double)(UAdsorbateBondRunning-UAdsorbateBond[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate UreyBradley energy-drift:                                  %lg\n",
                   (double)(UAdsorbateUreyBradleyRunning-UAdsorbateUreyBradley[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate bend energy-drift:                                         %lg\n",
                   (double)(UAdsorbateBendRunning-UAdsorbateBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate inversion-bend energy-drift:                               %lg\n",
                   (double)(UAdsorbateInversionBendRunning-UAdsorbateInversionBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate torsion energy-drift:                                      %lg\n",
                   (double)(UAdsorbateTorsionRunning-UAdsorbateTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate improper torsion energy-drift:                             %lg\n",
                   (double)(UAdsorbateImproperTorsionRunning-UAdsorbateImproperTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate out-of-plane energy-drift:                                 %lg\n",
                   (double)(UAdsorbateOutOfPlaneRunning-UAdsorbateOutOfPlane[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate stretch/stretch energy-drift:                              %lg\n",
                   (double)(UAdsorbateBondBondRunning-UAdsorbateBondBond[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate stretch/bend energy-drift:                                 %lg\n",
                   (double)(UAdsorbateBondBendRunning-UAdsorbateBondBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate bend/bend energy-drift:                                    %lg\n",
                   (double)(UAdsorbateBendBendRunning-UAdsorbateBendBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate stretch/torsion energy-drift:                              %lg\n",
                   (double)(UAdsorbateBondTorsionRunning-UAdsorbateBondTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate bend/torsion energy-drift:                                 %lg\n",
                   (double)(UAdsorbateBendTorsionRunning-UAdsorbateBendTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate intra VDW energy-drift:                                    %lg\n",
                   (double)(UAdsorbateIntraVDWRunning-UAdsorbateIntraVDW[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate intra charge-charge Coulomb energy-drift:                  %lg\n",
                   (double)(UAdsorbateIntraChargeChargeCoulombRunning-UAdsorbateIntraChargeCharge[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate intra charge-bonddipole Coulomb energy-drift:              %lg\n",
                   (double)(UAdsorbateIntraChargeBondDipoleCoulombRunning-UAdsorbateIntraChargeBondDipole[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Adsorbate intra bonddipole-bonddipole Coulomb energy-drift:          %lg\n",
                   (double)(UAdsorbateIntraBondDipoleBondDipoleCoulombRunning-UAdsorbateIntraBondDipoleBondDipole[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Cation stretch energy-drift:                                         %lg\n",
                   (double)(UCationBondRunning-UCationBond[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation UreyBradley energy-drift:                                     %lg\n",
                   (double)(UCationUreyBradleyRunning-UCationUreyBradley[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation bend energy-drift:                                            %lg\n",
                   (double)(UCationBendRunning-UCationBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation inversion-bend energy-drift:                                  %lg\n",
                   (double)(UCationInversionBendRunning-UCationInversionBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation torsion energy-drift:                                         %lg\n",
                   (double)(UCationTorsionRunning-UCationTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation improper torsion energy-drift:                                %lg\n",
                   (double)(UCationImproperTorsionRunning-UCationImproperTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation out-of-plane energy-drift:                                    %lg\n",
                   (double)(UCationOutOfPlaneRunning-UCationOutOfPlane[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation stretch/stretch energy-drift:                                 %lg\n",
                   (double)(UCationBondBondRunning-UCationBondBond[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation stretch/bend energy-drift:                                    %lg\n",
                   (double)(UCationBondBendRunning-UCationBondBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation bend/bend energy-drift:                                       %lg\n",
                   (double)(UCationBendBendRunning-UCationBendBend[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation stretch/torsion energy-drift:                                 %lg\n",
                   (double)(UCationBondTorsionRunning-UCationBondTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation bend/torsion energy-drift:                                    %lg\n",
                   (double)(UCationBendTorsionRunning-UCationBendTorsion[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation intra VDW energy-drift:                                       %lg\n",
                   (double)(UCationIntraVDWRunning-UCationIntraVDW[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation intra Coulomb charge-charge energy-drift:                     %lg\n",
                   (double)(UCationIntraChargeChargeCoulombRunning-UCationIntraChargeCharge[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation intra Coulomb charge-bonddipole energy-drift:                 %lg\n",
                   (double)(UCationIntraChargeBondDipoleCoulombRunning-UCationIntraChargeBondDipole[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Cation intra Coulomb bonddipole-bonddipole energy-drift:             %lg\n",
                   (double)(UCationIntraBondDipoleBondDipoleCoulombRunning-UCationIntraBondDipoleBondDipole[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Host/Host energy-drift:                                              %lg\n",
                   (double)(UHostHostRunning-UHostHost[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Host VDW energy-drift:                                        %lg\n",
                   (double)(UHostHostVDWRunning-UHostHostVDW[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Host Coulomb energy-drift:                                    %lg\n",
                   (double)(UHostHostCoulombRunning-UHostHostCoulomb[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Host Real charge-charge energy-drift:                       %lg\n",
                   (double)(UHostHostChargeChargeRealRunning-UHostHostChargeChargeReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Host Fourier charge-charge energy-drift:                    %lg\n",
                   (double)(UHostHostChargeChargeFourierRunning-UHostHostChargeChargeFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Host Real charge-bonddipole energy-drift:                   %lg\n",
                   (double)(UHostHostChargeBondDipoleRealRunning-UHostHostChargeBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Host Fourier charge-bonddipole energy-drift:                %lg\n",
                   (double)(UHostHostChargeBondDipoleFourierRunning-UHostHostChargeBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Host Real bonddipole-bonddipole energy-drift:               %lg\n",
                   (double)(UHostHostBondDipoleBondDipoleRealRunning-UHostHostBondDipoleBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Host Fourier bonddipole-bonddipole energy-drift:            %lg\n",
                   (double)(UHostHostBondDipoleBondDipoleFourierRunning-UHostHostBondDipoleBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);

  fprintf(FilePtr,"Host/Adsorbate energy-drift:                                         %lg\n",
                   (double)(UHostAdsorbateRunning-UHostAdsorbate[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Adsorbate VDW energy-drift:                                   %lg\n",
                   (double)(UHostAdsorbateVDWRunning-UHostAdsorbateVDW[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Adsorbate Coulomb energy-drift:                               %lg\n",
                   (double)(UHostAdsorbateCoulombRunning-UHostAdsorbateCoulomb[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Adsorbate Real charge-charge energy-drift:                  %lg\n",
                   (double)(UHostAdsorbateChargeChargeRealRunning-UHostAdsorbateChargeChargeReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Adsorbate Fourier charge-charge energy-drift:               %lg\n",
                   (double)(UHostAdsorbateChargeChargeFourierRunning-UHostAdsorbateChargeChargeFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Adsorbate Real charge-bonddipole energy-drift:              %lg\n",
                   (double)(UHostAdsorbateChargeBondDipoleRealRunning-UHostAdsorbateChargeBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Adsorbate Fourier charge-bonddipole energy-drift:           %lg\n",
                   (double)(UHostAdsorbateChargeBondDipoleFourierRunning-UHostAdsorbateChargeBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Adsorbate Real bonddipole-bonddipole energy-drift:          %lg\n",
                   (double)(UHostAdsorbateBondDipoleBondDipoleRealRunning-UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Adsorbate Fourier bonddipole-bonddipole energy-drift:       %lg\n",
                   (double)(UHostAdsorbateBondDipoleBondDipoleFourierRunning-UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);

  fprintf(FilePtr,"Host/Cation energy-drift:                                            %lg\n",
                   (double)(UHostCationRunning-UHostCation[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Cation VDW energy-drift:                                      %lg\n",
                   (double)(UHostCationVDWRunning-UHostCationVDW[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost/Cation Coulomb energy-drift:                                  %lg\n",
                   (double)(UHostCationCoulombRunning-UHostCationCoulomb[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Cation Real charge-charge energy-drift:                     %lg\n",
                   (double)(UHostCationChargeChargeRealRunning-UHostCationChargeChargeReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Cation Fourier charge-charge energy-drift:                  %lg\n",
                   (double)(UHostCationChargeChargeFourierRunning-UHostCationChargeChargeFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Cation Real charge-bonddipole energy-drift:                 %lg\n",
                   (double)(UHostCationChargeBondDipoleRealRunning-UHostCationChargeBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Cation Fourier charge-bonddipole energy-drift:              %lg\n",
                   (double)(UHostCationChargeBondDipoleFourierRunning-UHostCationChargeBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Cation Real bonddipole-bonddipole energy-drift:             %lg\n",
                   (double)(UHostCationBondDipoleBondDipoleRealRunning-UHostCationBondDipoleBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tHost/Cation Fourier bonddipole-bonddipole energy-drift:          %lg\n",
                   (double)(UHostCationBondDipoleBondDipoleFourierRunning-UHostCationBondDipoleBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);

  fprintf(FilePtr,"Adsorbate/Adsorbate energy-drift:                                     %lg\n",
                   (double)(UAdsorbateAdsorbateRunning-UAdsorbateAdsorbate[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Adsorbate VDW energy-drift:                               %lg\n",
                   (double)(UAdsorbateAdsorbateVDWRunning-UAdsorbateAdsorbateVDW[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Adsorbate Coulomb energy-drift:                           %lg\n",
                   (double)(UAdsorbateAdsorbateCoulombRunning-UAdsorbateAdsorbateCoulomb[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Adsorbate Real charge-charge energy-drift:              %lg\n",
                   (double)(UAdsorbateAdsorbateChargeChargeRealRunning-UAdsorbateAdsorbateChargeChargeReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Adsorbate Fourier charge-charge energy-drift:           %lg\n",
                   (double)(UAdsorbateAdsorbateChargeChargeFourierRunning-UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Adsorbate Real charge-bonddipole energy-drift:          %lg\n",
                   (double)(UAdsorbateAdsorbateChargeBondDipoleRealRunning-UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Adsorbate Fourier charge-bonddipole energy-drift:       %lg\n",
                   (double)(UAdsorbateAdsorbateChargeBondDipoleFourierRunning-UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Adsorbate Real bonddipole-bonddipole energy-drift:      %lg\n",
                   (double)(UAdsorbateAdsorbateBondDipoleBondDipoleRealRunning-UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Adsorbate Fourier bonddipole-bonddipole energy-drift:   %lg\n",
                   (double)(UAdsorbateAdsorbateBondDipoleBondDipoleFourierRunning-UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);

  fprintf(FilePtr,"Cation/Cation energy-drift:                                           %lg\n",
                   (double)(UCationCationRunning-UCationCation[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation/Cation VDW energy-drift:                                     %lg\n",
                   (double)(UCationCationVDWRunning-UCationCationVDW[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation/Cation Coulomb energy-drift:                                 %lg\n",
                   (double)(UCationCationCoulombRunning-UCationCationCoulomb[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tCation/Cation Real charge-charge energy-drift:                    %lg\n",
                   (double)(UCationCationChargeChargeRealRunning-UCationCationChargeChargeReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tCation/Cation Fourier charge-charge energy-drift:                 %lg\n",
                   (double)(UCationCationChargeChargeFourierRunning-UCationCationChargeChargeFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tCation/Cation Real charge-bonddipole energy-drift:                %lg\n",
                   (double)(UCationCationChargeBondDipoleRealRunning-UCationCationChargeBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tCation/Cation Fourier charge-bonddipole energy-drift:             %lg\n",
                   (double)(UCationCationChargeBondDipoleFourierRunning-UCationCationChargeBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tCation/Cation Real bonddipole-bonddipole energy-drift:            %lg\n",
                   (double)(UCationCationBondDipoleBondDipoleRealRunning-UCationCationBondDipoleBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tCation/Cation Fourier bonddipole-bonddipole energy-drift:         %lg\n",
                   (double)(UCationCationBondDipoleBondDipoleFourierRunning-UCationCationBondDipoleBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);

  fprintf(FilePtr,"Adsorbate/Cation energy-drift:                                        %lg\n",
                   (double)(UAdsorbateCationRunning-UAdsorbateCation[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Cation VDW energy-drift:                                  %lg\n",
                   (double)(UAdsorbateCationVDWRunning-UAdsorbateCationVDW[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate/Cation Coulomb energy-drift:                              %lg\n",
                   (double)(UAdsorbateCationCoulombRunning-UAdsorbateCationCoulomb[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Cation Real charge-charge energy-drift:                 %lg\n",
                   (double)(UAdsorbateCationChargeChargeRealRunning-UAdsorbateCationChargeChargeReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Cation Fourier charge-charge energy-drift:              %lg\n",
                   (double)(UAdsorbateCationChargeChargeFourierRunning-UAdsorbateCationChargeChargeFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Cation Real charge-bonddipole energy-drift:             %lg\n",
                   (double)(UAdsorbateCationChargeBondDipoleRealRunning-UAdsorbateCationChargeBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Cation Fourier charge-bonddipole energy-drift:          %lg\n",
                   (double)(UAdsorbateCationChargeBondDipoleFourierRunning-UAdsorbateCationChargeBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Cation Real bonddipole-bonddipole energy-drift:         %lg\n",
                   (double)(UAdsorbateCationBondDipoleBondDipoleRealRunning-UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\t\tAdsorbate/Cation Fourier bonddipole-bonddipole energy-drift:      %lg\n",
                   (double)(UAdsorbateCationBondDipoleBondDipoleFourierRunning-UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Polarization energy-drift:\n");
  fprintf(FilePtr,"\tHost polarization energy-drift:                %lg\n",
                   (double)(UHostPolarizationRunning-UHostPolarization[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate polarization energy-drift:           %lg\n",
                   (double)(UAdsorbatePolarizationRunning-UAdsorbatePolarization[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation polarization energy-drift:              %lg\n",
                   (double)(UCationPolarizationRunning-UCationPolarization[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tHost back-polarization energy-drift:                %lg\n",
                   (double)(UHostBackPolarizationRunning-UHostBackPolarization[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tAdsorbate back-polarization energy-drift:           %lg\n",
                   (double)(UAdsorbateBackPolarizationRunning-UAdsorbateBackPolarization[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCation back-polarization energy-drift:              %lg\n",
                   (double)(UCationBackPolarizationRunning-UCationBackPolarization[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Tail-correction energy-drift:                  %lg\n",
                   (double)(UTailCorrectionRunning-UTailCorrection[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Distance constraints energy-drift:                  %lg\n",
                   (double)(UDistanceConstraintsRunning-UDistanceConstraints[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Angle constraints energy-drift:                     %lg\n",
                   (double)(UAngleConstraintsRunning-UAngleConstraints[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Dihedral constraints energy-drift:                  %lg\n",
                   (double)(UDihedralConstraintsRunning-UDihedralConstraints[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Inversion-bend constraints energy-drift:                    %lg\n",
                   (double)(UInversionBendConstraintsRunning-UInversionBendConstraints[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Out-of-plane distance constraints energy-drift:                    %lg\n",
                   (double)(UOutOfPlaneDistanceConstraintsRunning-UOutOfPlaneDistanceConstraints[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"Exclusion constraints energy-drift:                 %lg\n",
                   (double)(UExclusionConstraintsRunning-UExclusionConstraints[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"===================================================================\n");
  fprintf(FilePtr,"Total energy-drift: %lg\n",(double)(UTotalRunning-UTotal[CurrentSystem])*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\n\n");

  if(fabs(UTotalRunning-UTotal[CurrentSystem])*ENERGY_TO_KELVIN>1e-2)
  {
    if(NumberOfWarnings[NumberOfWarnings[CurrentSystem]]<MAX_NUMBER_OF_WARNINGS)
    {
      Warnings[CurrentSystem][NumberOfWarnings[CurrentSystem]]=ENERGY_DRIFT;
      NumberOfWarnings[CurrentSystem]++;
    }
  }
}

void PrintRestartFile(void)
{
  int i,j,k,l;
  FILE *FilePtrOut;
  char buffer[1024];
  int index;
  int ncell;
  int FractionalMolecule;
  VECTOR flexible_drift;
  VECTOR com;
  REAL shift;

  if (STREAM)
      return;

  mkdir("Restart",S_IRWXU);
  sprintf(buffer,"Restart/System_%d",CurrentSystem);
  mkdir(buffer,S_IRWXU);

  sprintf(buffer,"Restart/System_%d/restart_%s_%d.%d.%d_%lf_%lg%s",
          CurrentSystem,
          Framework[CurrentSystem].Name[0],
          NumberOfUnitCells[CurrentSystem].x,
          NumberOfUnitCells[CurrentSystem].y,
          NumberOfUnitCells[CurrentSystem].z,
          (double)therm_baro_stats.ExternalTemperature[CurrentSystem],
          (double)therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR,
          FileNameAppend);

  FilePtrOut=fopen(buffer,"w");

  fprintf(FilePtrOut,"Cell info:\n");
  fprintf(FilePtrOut,"========================================================================\n");
  fprintf(FilePtrOut,"number-of-unit-cells: %d %d %d\n",NumberOfUnitCells[CurrentSystem].x,NumberOfUnitCells[CurrentSystem].y,NumberOfUnitCells[CurrentSystem].z);
  fprintf(FilePtrOut,"unit-cell-vector-a: %18.12f %18.12f %18.12f\n",(double)UnitCellBox[CurrentSystem].ax,(double)UnitCellBox[CurrentSystem].ay,(double)UnitCellBox[CurrentSystem].az);
  fprintf(FilePtrOut,"unit-cell-vector-b: %18.12f %18.12f %18.12f\n",(double)UnitCellBox[CurrentSystem].bx,(double)UnitCellBox[CurrentSystem].by,(double)UnitCellBox[CurrentSystem].bz);
  fprintf(FilePtrOut,"unit-cell-vector-c: %18.12f %18.12f %18.12f\n",(double)UnitCellBox[CurrentSystem].cx,(double)UnitCellBox[CurrentSystem].cy,(double)UnitCellBox[CurrentSystem].cz);
  fprintf(FilePtrOut,"\n");
  fprintf(FilePtrOut,"cell-vector-a: %18.12f %18.12f %18.12f\n",(double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].az);
  fprintf(FilePtrOut,"cell-vector-b: %18.12f %18.12f %18.12f\n",(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].bz);
  fprintf(FilePtrOut,"cell-vector-c: %18.12f %18.12f %18.12f\n",(double)Box[CurrentSystem].cx,(double)Box[CurrentSystem].cy,(double)Box[CurrentSystem].cz);
  fprintf(FilePtrOut,"\n");

  fprintf(FilePtrOut,"cell-lengths: %18.12f %18.12f %18.12f\n",UnitCellSize[CurrentSystem].x,UnitCellSize[CurrentSystem].y,UnitCellSize[CurrentSystem].z);
  fprintf(FilePtrOut,"cell-angles: %18.12f %18.12f %18.12f\n",AlphaAngle[CurrentSystem]*RAD2DEG,BetaAngle[CurrentSystem]*RAD2DEG,GammaAngle[CurrentSystem]*RAD2DEG);
  fprintf(FilePtrOut,"\n\n");

  flexible_drift.x=flexible_drift.y=flexible_drift.z=0.0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    com=GetFrameworkCenterOfMass();
    flexible_drift.x=com.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x;
    flexible_drift.y=com.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y;
    flexible_drift.z=com.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z;
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    fprintf(FilePtrOut,"Framework:\n");
    fprintf(FilePtrOut,"========================================================================\n");
    fprintf(FilePtrOut,"Initial-framework-center-of-mass: %18.12f %18.12f %18.12f\n",(double)Framework[CurrentSystem].IntialCenterOfMassPosition.x,
      (double)Framework[CurrentSystem].IntialCenterOfMassPosition.y,(double)Framework[CurrentSystem].IntialCenterOfMassPosition.z);
    for(j=0;j<Framework[CurrentSystem].NumberOfFrameworks;j++)
    {
      fprintf(FilePtrOut,"Maximum-translation-change framework %d: %lf\n", j,
         (double)FrameworkMaximumTranslation[CurrentSystem][j]);
      fprintf(FilePtrOut,"Maximum-translation-shift-change framework %d: %lf %lf %lf\n", j,
         (double)FrameworkMaximumShiftTranslation[CurrentSystem][j].x,(double)FrameworkMaximumShiftTranslation[CurrentSystem][j].y,(double)FrameworkMaximumShiftTranslation[CurrentSystem][j].z);
    }
    fprintf(FilePtrOut,"\n\n");
  }

  fprintf(FilePtrOut,"Maximum changes for MC-moves:\n");
  fprintf(FilePtrOut,"========================================================================\n");
  fprintf(FilePtrOut,"Maximum-volume-change: %lf\n",(double)MaximumVolumeChange[CurrentSystem]);
  fprintf(FilePtrOut,"Maximum-Gibbs-volume-change: %lf\n",(double)MaximumGibbsVolumeChange[CurrentSystem]);
  fprintf(FilePtrOut,"Maximum-box-shape-change: %lf %lf %lf, %lf %lf %lf, %lf %lf %lf\n",
       (double)MaximumBoxShapeChange[CurrentSystem].ax,(double)MaximumBoxShapeChange[CurrentSystem].bx,(double)MaximumBoxShapeChange[CurrentSystem].cx,
       (double)MaximumBoxShapeChange[CurrentSystem].ay,(double)MaximumBoxShapeChange[CurrentSystem].by,(double)MaximumBoxShapeChange[CurrentSystem].cy,
       (double)MaximumBoxShapeChange[CurrentSystem].az,(double)MaximumBoxShapeChange[CurrentSystem].bz,(double)MaximumBoxShapeChange[CurrentSystem].cz);
  fprintf(FilePtrOut,"\n\n");

  fprintf(FilePtrOut,"Acceptance targets for MC-moves:\n");
  fprintf(FilePtrOut,"========================================================================\n");
  fprintf(FilePtrOut,"Target-volume-change: %lf\n",(double)TargetAccRatioVolumeChange);
  fprintf(FilePtrOut,"Target-box-shape-change: %lf\n",(double)TargetAccRatioBoxShapeChange);
  fprintf(FilePtrOut,"Target-Gibbs-volume-change: %lf\n",(double)TargetAccRatioGibbsVolumeChange);
  fprintf(FilePtrOut,"\n\n");


  


  fprintf(FilePtrOut,"Components: %d (Adsorbates %d, Cations %d)\n",NumberOfComponents,
     NumberOfAdsorbateMolecules[CurrentSystem],NumberOfCationMolecules[CurrentSystem]);
  fprintf(FilePtrOut,"========================================================================\n");
  for(l=0;l<NumberOfComponents;l++)
  {
    // write CB/CFMC biasing data
    FractionalMolecule=Components[l].FractionalMolecule[CurrentSystem];
    if(FractionalMolecule>=0)
    {
      fprintf(FilePtrOut,"Component %d (%s)\n",l,Components[l].Name);
      fprintf(FilePtrOut,"\tFractional-molecule-id component %d: %d\n",l,Components[l].FractionalMolecule[CurrentSystem]);

      fprintf(FilePtrOut,"\tLambda-factors component %d: ",l);
      for(i=0;i<Components[l].NumberOfAtoms;i++)
      {
        if(Components[l].ExtraFrameworkMolecule)
          fprintf(FilePtrOut," %lf",Cations[CurrentSystem][FractionalMolecule].Atoms[i].CFVDWScalingParameter);
        else
          fprintf(FilePtrOut," %lf",Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFVDWScalingParameter);
      }
      fprintf(FilePtrOut,"\n");

      fprintf(FilePtrOut,"\tNumber-of-biasing-factors component %d: %d\n",l,Components[l].CFLambdaHistogramSize);
      fprintf(FilePtrOut,"\tBiasing-factors component %d: ",l);
      for(i=0;i<Components[l].CFLambdaHistogramSize;i++)
        fprintf(FilePtrOut," %lf",Components[l].CFBiasingFactors[CurrentSystem][i]);
      fprintf(FilePtrOut,"\n");

      fprintf(FilePtrOut,"\tMaximum-CF-Lambda-change component %d: %lf\n",l,MaximumCFLambdaChange[CurrentSystem][l]);
      fprintf(FilePtrOut,"\tMaximum-CBCF-Lambda-change component %d: %lf\n",l,MaximumCBCFLambdaChange[CurrentSystem][l]);
      fprintf(FilePtrOut,"\n");
    }
    fprintf(FilePtrOut,"\tMaximum-translation-change component %d: %lf,%lf,%lf\n",l,
       (double)MaximumTranslation[CurrentSystem][l].x,(double)MaximumTranslation[CurrentSystem][l].y,(double)MaximumTranslation[CurrentSystem][l].z);
    fprintf(FilePtrOut,"\tMaximum-translation-in-plane-change component %d: %lf,%lf,%lf\n",l,
       (double)MaximumTranslationInPlane[CurrentSystem][l].x,(double)MaximumTranslationInPlane[CurrentSystem][l].y,(double)MaximumTranslationInPlane[CurrentSystem][l].z);
    fprintf(FilePtrOut,"\tMaximum-rotation-change component %d: %lf %lf %lf\n",l,(double)MaximumRotation[CurrentSystem][l].x,
            (double)MaximumRotation[CurrentSystem][l].y,(double)MaximumRotation[CurrentSystem][l].z);
  }
  fprintf(FilePtrOut,"\n");

  // write CFRXMC biasing data
  fprintf(FilePtrOut,"Reactions: %d\n",NumberOfReactions);
  if(NumberOfReactions>0)
  {
    fprintf(FilePtrOut,"========================================================================\n");
    fprintf(FilePtrOut,"Number-of-biasing-factors-reaction: %d\n",RXMCLambdaHistogramSize);
    for(i=0;i<NumberOfReactions;i++)
    {
      fprintf(FilePtrOut,"Reaction %d\n",i);
      fprintf(FilePtrOut,"\tLambda-factor reaction %d: %lf\n",i,CFCRXMCLambda[CurrentSystem][i]);
      fprintf(FilePtrOut,"\tMaximum-Lambda-reaction-change reaction %d: %lf\n",i,MaximumReactionLambdaChange[CurrentSystem][i]);
      fprintf(FilePtrOut,"\tFractional-product-molecules reaction %d:",i);
      for(j=0;j<NumberOfComponents;j++)
      {
        if(ProductsStoichiometry[i][j]>0)
        {
          for(k=0;k<ProductsStoichiometry[i][j];k++)
            fprintf(FilePtrOut," %d",Components[j].ProductFractionalMolecules[CurrentSystem][i][k]);
        }
      }
      fprintf(FilePtrOut,"\n");
      fprintf(FilePtrOut,"\tFractional-reactant-molecules reaction %d:",i);
      for(j=0;j<NumberOfComponents;j++)
      {
        if(ReactantsStoichiometry[i][j]>0)
        {
          for(k=0;k<ReactantsStoichiometry[i][j];k++)
            fprintf(FilePtrOut," %d",Components[j].ReactantFractionalMolecules[CurrentSystem][i][k]);
        }
      }
      fprintf(FilePtrOut,"\n");

      shift=RXMCBiasingFactors[CurrentSystem][i][0];
      fprintf(FilePtrOut,"\tReaction-biasing-factors reaction %d: ",i);
      for(k=0;k<RXMCLambdaHistogramSize;k++)
        fprintf(FilePtrOut,"%lf ",RXMCBiasingFactors[CurrentSystem][i][k]-shift);
      fprintf(FilePtrOut,"\n");
    }
    fprintf(FilePtrOut,"\n");
  }
  fprintf(FilePtrOut,"\n");

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    fprintf(FilePtrOut,"Framework atomic positions\n");
    fprintf(FilePtrOut,"========================================================================\n");

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        fprintf(FilePtrOut,"Framework-atom-position: %d %d %18.12f %18.12f %18.12f\n",
          CurrentFramework,j,
          (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Position.x-flexible_drift.x,
          (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Position.y-flexible_drift.y,
          (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Position.z-flexible_drift.z);
    }
    fprintf(FilePtrOut,"\n");

    fprintf(FilePtrOut,"Framework atomic velocities\n");
    fprintf(FilePtrOut,"========================================================================\n");
    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        fprintf(FilePtrOut,"Framework-atom-velocity: %d %d %18.12f %18.12f %18.12f\n",
          CurrentFramework,j,
          (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Velocity.x,
          (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Velocity.y,
          (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Velocity.z);
    }
    fprintf(FilePtrOut,"\n");

    fprintf(FilePtrOut,"Framework atomic forces\n");
    fprintf(FilePtrOut,"========================================================================\n");
    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        fprintf(FilePtrOut,"Framework-atom-force:    %d %d %18.12f %18.12f %18.12f\n",
          CurrentFramework,j,
          (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Force.x,
          (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Force.y,
          (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Force.z);
    }
    fprintf(FilePtrOut,"\n");

    fprintf(FilePtrOut,"Framework atomic charges\n");
    fprintf(FilePtrOut,"========================================================================\n");
    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        fprintf(FilePtrOut,"Framework-atom-charge:   %d %d %18.12f\n",
          CurrentFramework,j,
          (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Charge);
    }
    fprintf(FilePtrOut,"\n");

    fprintf(FilePtrOut,"Framework atomic fixed/free\n");
    fprintf(FilePtrOut,"========================================================================\n");
    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        fprintf(FilePtrOut,"Framework-atom-fixed:    %d %d       %d %d %d\n",
          CurrentFramework,j,
          Framework[CurrentSystem].Atoms[CurrentFramework][j].Fixed.x,
          Framework[CurrentSystem].Atoms[CurrentFramework][j].Fixed.y,
          Framework[CurrentSystem].Atoms[CurrentFramework][j].Fixed.z);
    }
    fprintf(FilePtrOut,"\n");
  }

  for(j=0;j<NumberOfComponents;j++)
  {
    fprintf(FilePtrOut,"Component: %-5d %s    %d molecules of %s\n",
      j,
      Components[j].ExtraFrameworkMolecule?"Cation   ":"Adsorbate",
      Components[j].NumberOfMolecules[CurrentSystem],
      Components[j].Name);
    fprintf(FilePtrOut,"------------------------------------------------------------------------\n");

    if(Components[j].ExtraFrameworkMolecule)
    {
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        if(Cations[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Cation-atom-position: %d %d %18.12f %18.12f %18.12f\n",
              k,l,
              (double)Cations[CurrentSystem][k].Atoms[l].Position.x-flexible_drift.x,
              (double)Cations[CurrentSystem][k].Atoms[l].Position.y-flexible_drift.y,
              (double)Cations[CurrentSystem][k].Atoms[l].Position.z-flexible_drift.z);
        }
      }
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        if(Cations[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Cation-atom-velocity: %d %d %18.12f %18.12f %18.12f\n",
              k,l,
              (double)Cations[CurrentSystem][k].Atoms[l].Velocity.x,
              (double)Cations[CurrentSystem][k].Atoms[l].Velocity.y,
              (double)Cations[CurrentSystem][k].Atoms[l].Velocity.z);
        }
      }
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        if(Cations[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Cation-atom-force:    %d %d %18.12f %18.12f %18.12f\n",
              k,l,
              (double)Cations[CurrentSystem][k].Atoms[l].Force.x,
              (double)Cations[CurrentSystem][k].Atoms[l].Force.y,
              (double)Cations[CurrentSystem][k].Atoms[l].Force.z);
        }
      }
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        if(Cations[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Cation-atom-charge:   %d %d %18.12f\n",
              k,l,
              (double)Cations[CurrentSystem][k].Atoms[l].Charge);
        }
      }
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        if(Cations[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Cation-atom-scaling:  %d %d %18.12f\n",
              k,l,
              (double)Cations[CurrentSystem][k].Atoms[l].CFVDWScalingParameter);
        }
      }
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        if(Cations[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Cation-atom-fixed:    %d %d       %d %d %d\n",
              k,l,
              Cations[CurrentSystem][k].Atoms[l].Fixed.x,
              Cations[CurrentSystem][k].Atoms[l].Fixed.y,
              Cations[CurrentSystem][k].Atoms[l].Fixed.z);
        }
      }
    }
    else
    {
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        if(Adsorbates[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Adsorbate-atom-position: %d %d %18.12f %18.12f %18.12f\n",
              k,l,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].Position.x-flexible_drift.x,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].Position.y-flexible_drift.y,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].Position.z-flexible_drift.z);
        }
      }
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        if(Adsorbates[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Adsorbate-atom-velocity: %d %d %18.12f %18.12f %18.12f\n",
              k,l,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].Velocity.x,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].Velocity.y,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].Velocity.z);
        }
      }
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        if(Adsorbates[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Adsorbate-atom-force:    %d %d %18.12f %18.12f %18.12f\n",
              k,l,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].Force.x,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].Force.y,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].Force.z);
        }
      }
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        if(Adsorbates[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Adsorbate-atom-charge:   %d %d %18.12f\n",
              k,l,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].Charge);
        }
      }
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        if(Adsorbates[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Adsorbate-atom-scaling:  %d %d %18.12f\n",
              k,l,
              (double)Adsorbates[CurrentSystem][k].Atoms[l].CFVDWScalingParameter);
        }
      }
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        if(Adsorbates[CurrentSystem][k].Type==j)
        {
          for(l=0;l<Components[j].NumberOfAtoms;l++)
            fprintf(FilePtrOut,"Adsorbate-atom-fixed:    %d %d       %d %d %d\n",
              k,l,
              Adsorbates[CurrentSystem][k].Atoms[l].Fixed.x,
              Adsorbates[CurrentSystem][k].Atoms[l].Fixed.y,
              Adsorbates[CurrentSystem][k].Atoms[l].Fixed.z);
        }
      }
    }

    fprintf(FilePtrOut,"\n");
  }
  fprintf(FilePtrOut,"\n");
  fclose(FilePtrOut);

  if(UseReplicas[CurrentSystem])
  {
    sprintf(buffer,"Restart/System_%d/restart_replicas_%s_%d.%d.%d_%lf_%lg%s",
            CurrentSystem,
            Framework[CurrentSystem].Name[0],1,1,1,
            (double)therm_baro_stats.ExternalTemperature[CurrentSystem],
            (double)therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR,
            FileNameAppend);

    FilePtrOut=fopen(buffer,"w");

    fprintf(FilePtrOut,"Cell info:\n");
    fprintf(FilePtrOut,"========================================================================\n");
    fprintf(FilePtrOut,"number-of-unit-cells: %d %d %d\n",NumberOfUnitCells[CurrentSystem].x,NumberOfUnitCells[CurrentSystem].y,NumberOfUnitCells[CurrentSystem].z);
    fprintf(FilePtrOut,"unit-cell-vector-a: %18.12f %18.12f %18.12f\n",(double)UnitCellBox[CurrentSystem].ax,(double)UnitCellBox[CurrentSystem].ay,(double)UnitCellBox[CurrentSystem].az);
    fprintf(FilePtrOut,"unit-cell-vector-b: %18.12f %18.12f %18.12f\n",(double)UnitCellBox[CurrentSystem].bx,(double)UnitCellBox[CurrentSystem].by,(double)UnitCellBox[CurrentSystem].bz);
    fprintf(FilePtrOut,"unit-cell-vector-c: %18.12f %18.12f %18.12f\n",(double)UnitCellBox[CurrentSystem].cx,(double)UnitCellBox[CurrentSystem].cy,(double)UnitCellBox[CurrentSystem].cz);
    fprintf(FilePtrOut,"\n");
    fprintf(FilePtrOut,"cell-vector-a: %18.12f %18.12f %18.12f\n",(double)ReplicaBox[CurrentSystem].ax,(double)ReplicaBox[CurrentSystem].ay,(double)ReplicaBox[CurrentSystem].az);
    fprintf(FilePtrOut,"cell-vector-b: %18.12f %18.12f %18.12f\n",(double)ReplicaBox[CurrentSystem].bx,(double)ReplicaBox[CurrentSystem].by,(double)ReplicaBox[CurrentSystem].bz);
    fprintf(FilePtrOut,"cell-vector-c: %18.12f %18.12f %18.12f\n",(double)ReplicaBox[CurrentSystem].cx,(double)ReplicaBox[CurrentSystem].cy,(double)ReplicaBox[CurrentSystem].cz);
    fprintf(FilePtrOut,"\n");

    fprintf(FilePtrOut,"cell-lengths: %18.12f %18.12f %18.12f\n",UnitCellSize[CurrentSystem].x,UnitCellSize[CurrentSystem].y,UnitCellSize[CurrentSystem].z);
    fprintf(FilePtrOut,"cell-angles: %18.12f %18.12f %18.12f\n",AlphaAngle[CurrentSystem]*RAD2DEG,BetaAngle[CurrentSystem]*RAD2DEG,GammaAngle[CurrentSystem]*RAD2DEG);
    fprintf(FilePtrOut,"\n\n");

    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      fprintf(FilePtrOut,"Framework:\n");
      fprintf(FilePtrOut,"========================================================================\n");
      fprintf(FilePtrOut,"Initial-framework-center-of-mass: %18.12f %18.12f %18.12f\n",(double)Framework[CurrentSystem].IntialCenterOfMassPosition.x,
        (double)Framework[CurrentSystem].IntialCenterOfMassPosition.y,(double)Framework[CurrentSystem].IntialCenterOfMassPosition.z);
      for(j=0;j<Framework[CurrentSystem].NumberOfFrameworks;j++)
      {
        fprintf(FilePtrOut,"Maximum-translation-change framework %d: %lf\n", j,
           (double)FrameworkMaximumTranslation[CurrentSystem][j]);
        fprintf(FilePtrOut,"Maximum-translation-shift-change framework %d: %lf %lf %lf\n", j,
           (double)FrameworkMaximumShiftTranslation[CurrentSystem][j].x,(double)FrameworkMaximumShiftTranslation[CurrentSystem][j].y,(double)FrameworkMaximumShiftTranslation[CurrentSystem][j].z);
      }
      fprintf(FilePtrOut,"\n\n");
    }

    fprintf(FilePtrOut,"Maximum changes for MC-moves:\n");
    fprintf(FilePtrOut,"========================================================================\n");
    fprintf(FilePtrOut,"Maximum-volume-change: %lf\n",(double)MaximumVolumeChange[CurrentSystem]);
    fprintf(FilePtrOut,"Maximum-Gibbs-volume-change: %lf\n",(double)MaximumGibbsVolumeChange[CurrentSystem]);
    fprintf(FilePtrOut,"Maximum-box-shape-change: %lf %lf %lf, %lf %lf %lf, %lf %lf %lf\n",
         (double)MaximumBoxShapeChange[CurrentSystem].ax,(double)MaximumBoxShapeChange[CurrentSystem].bx,(double)MaximumBoxShapeChange[CurrentSystem].cx,
         (double)MaximumBoxShapeChange[CurrentSystem].ay,(double)MaximumBoxShapeChange[CurrentSystem].by,(double)MaximumBoxShapeChange[CurrentSystem].cy,
         (double)MaximumBoxShapeChange[CurrentSystem].az,(double)MaximumBoxShapeChange[CurrentSystem].bz,(double)MaximumBoxShapeChange[CurrentSystem].cz);
    fprintf(FilePtrOut,"\n\n");

    fprintf(FilePtrOut,"Acceptance targets for MC-moves:\n");
    fprintf(FilePtrOut,"========================================================================\n");
    fprintf(FilePtrOut,"Target-volume-change: %lf\n",(double)TargetAccRatioVolumeChange);
    fprintf(FilePtrOut,"Target-box-shape-change: %lf\n",(double)TargetAccRatioBoxShapeChange);
    fprintf(FilePtrOut,"Target-Gibbs-volume-change: %lf\n",(double)TargetAccRatioGibbsVolumeChange);
    fprintf(FilePtrOut,"\n\n");

    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      index=0;
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        {
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            fprintf(FilePtrOut,"Framework-atom-position: %d %d %18.12f %18.12f %18.12f\n",
              index,CurrentFramework,
              (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Position.x+ReplicaShift[ncell].x-flexible_drift.x,
              (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Position.y+ReplicaShift[ncell].y-flexible_drift.y,
              (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Position.z+ReplicaShift[ncell].z-flexible_drift.z);
            fprintf(FilePtrOut,"Framework-atom-velocity: %d %d %18.12f %18.12f %18.12f\n",
              index,CurrentFramework,
              (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Velocity.x,
              (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Velocity.y,
              (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Velocity.z);
            fprintf(FilePtrOut,"Framework-atom-force:    %d %d %18.12f %18.12f %18.12f\n",
              index,CurrentFramework,
              (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Force.x,
              (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Force.y,
              (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Force.z);
            fprintf(FilePtrOut,"Framework-atom-charge:   %d %d %18.12f\n",
              index,CurrentFramework,
              (double)Framework[CurrentSystem].Atoms[CurrentFramework][j].Charge);
            fprintf(FilePtrOut,"Framework-atom-fixed:    %d %d       %d %d %d\n",
              index,CurrentFramework,
              Framework[CurrentSystem].Atoms[CurrentFramework][j].Fixed.x,
              Framework[CurrentSystem].Atoms[CurrentFramework][j].Fixed.y,
              Framework[CurrentSystem].Atoms[CurrentFramework][j].Fixed.z);
          }
          index++;
        }
      }
      fprintf(FilePtrOut,"\n");
    }

    index=0;
    for(j=0;j<NumberOfComponents;j++)
    {
      fprintf(FilePtrOut,"Component: %-5d %s    %d molecules of %s\n",
        j,
        Components[j].ExtraFrameworkMolecule?"Cation   ":"Adsorbate",
        TotalNumberOfReplicaCells[CurrentSystem]*Components[j].NumberOfMolecules[CurrentSystem],
        Components[j].Name);
      fprintf(FilePtrOut,"------------------------------------------------------------------------\n");

      if(Components[j].ExtraFrameworkMolecule)
      {
        index=0;
        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
          {
            if(Cations[CurrentSystem][k].Type==j)
            {
              for(l=0;l<Components[j].NumberOfAtoms;l++)
              {
                fprintf(FilePtrOut,"Cation-atom-position: %d %d %18.12f %18.12f %18.12f\n",
                  index,l,
                  (double)Cations[CurrentSystem][k].Atoms[l].Position.x+ReplicaShift[ncell].x-flexible_drift.x,
                  (double)Cations[CurrentSystem][k].Atoms[l].Position.y+ReplicaShift[ncell].y-flexible_drift.y,
                  (double)Cations[CurrentSystem][k].Atoms[l].Position.z+ReplicaShift[ncell].z-flexible_drift.z);
                fprintf(FilePtrOut,"Cation-atom-velocity: %d %d %18.12f %18.12f %18.12f\n",
                  index,l,
                  (double)Cations[CurrentSystem][k].Atoms[l].Velocity.x,
                  (double)Cations[CurrentSystem][k].Atoms[l].Velocity.y,
                  (double)Cations[CurrentSystem][k].Atoms[l].Velocity.z);
                fprintf(FilePtrOut,"Cation-atom-force:    %d %d %18.12f %18.12f %18.12f\n",
                  index,l,
                  (double)Cations[CurrentSystem][k].Atoms[l].Force.x,
                  (double)Cations[CurrentSystem][k].Atoms[l].Force.y,
                  (double)Cations[CurrentSystem][k].Atoms[l].Force.z);
                fprintf(FilePtrOut,"Cation-atom-charge:   %d %d %18.12f\n",
                  index,l,
                  (double)Cations[CurrentSystem][k].Atoms[l].Charge);
                fprintf(FilePtrOut,"Cation-atom-scaling:  %d %d %18.12f\n",
                  index,l,
                  (double)Cations[CurrentSystem][k].Atoms[l].CFVDWScalingParameter);
                fprintf(FilePtrOut,"Cation-atom-fixed:    %d %d       %d %d %d\n",
                  index,l,
                  Cations[CurrentSystem][k].Atoms[l].Fixed.x,
                  Cations[CurrentSystem][k].Atoms[l].Fixed.y,
                  Cations[CurrentSystem][k].Atoms[l].Fixed.z);
              }
              index++;
            }
          }
        }
      }
      else
      {
        index=0;
        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
          {
            if(Adsorbates[CurrentSystem][k].Type==j)
            {
              for(l=0;l<Components[j].NumberOfAtoms;l++)
              {
                fprintf(FilePtrOut,"Adsorbate-atom-position: %d %d %18.12f %18.12f %18.12f\n",
                  index,l,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].Position.x+ReplicaShift[ncell].x-flexible_drift.x,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].Position.y+ReplicaShift[ncell].y-flexible_drift.y,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].Position.z+ReplicaShift[ncell].z-flexible_drift.z);
                fprintf(FilePtrOut,"Adsorbate-atom-velocity: %d %d %18.12f %18.12f %18.12f\n",
                  index,l,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].Velocity.x,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].Velocity.y,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].Velocity.z);
                fprintf(FilePtrOut,"Adsorbate-atom-force:    %d %d %18.12f %18.12f %18.12f\n",
                  index,l,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].Force.x,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].Force.y,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].Force.z);
                fprintf(FilePtrOut,"Adsorbate-atom-charge:   %d %d %18.12f\n",
                  index,l,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].Charge);
                fprintf(FilePtrOut,"Adsorbate-atom-scaling:  %d %d %18.12f\n",
                  index,l,
                  (double)Adsorbates[CurrentSystem][k].Atoms[l].CFVDWScalingParameter);
                fprintf(FilePtrOut,"Adsorbate-atom-fixed:    %d %d       %d %d %d\n",
                  index,l,
                  Adsorbates[CurrentSystem][k].Atoms[l].Fixed.x,
                  Adsorbates[CurrentSystem][k].Atoms[l].Fixed.y,
                  Adsorbates[CurrentSystem][k].Atoms[l].Fixed.z);
              }
              index++;
            }
          }
        }
      }
      fprintf(FilePtrOut,"\n");
    }
    fprintf(FilePtrOut,"\n");
    fclose(FilePtrOut);
  }
}


void WriteBinaryRestartFiles(void)
{
  FILE *FilePtr;
  char buffer[1024];

  if (STREAM)
    return;

  fprintf(stderr, "Writing Crash-file!: %lld\n",CurrentCycle);

  mkdir("CrashRestart",S_IRWXU);

  sprintf(buffer,"CrashRestart/binary_restart%s.dat",FileNameAppend);

  FilePtr=fopen(buffer,"w");

  WriteRestartConstants(FilePtr);
  WriteRestartSimulation(FilePtr);
  WriteRestartWarnings(FilePtr);
  WriteRestartPseudoAtoms(FilePtr);
  WriteRestartComponent(FilePtr);
  WriteRestartMolecules(FilePtr);
  WriteRestartFramework(FilePtr);
  WriteRestartCBMC(FilePtr);
  WriteRestartEwald(FilePtr);
  WriteRestartStatistics(FilePtr);
  WriteRestartMcMoves(FilePtr);
  WriteRestartSample(FilePtr);
  WriteRestartThermoBarostats(FilePtr);
  WriteRestartEquationOfState(FilePtr);
  WriteRestartGrids(FilePtr);
  WriteRestartMinimization(FilePtr);
  WriteRestartUtils(FilePtr);
  WriteRestartMovies(FilePtr);
  WriteRestartOutput(FilePtr);

  fclose(FilePtr);
}

void ReadRestartOutput(FILE* FilePtr)
{
  int i;
  fpos_t pos;
  char buffer[1024],buffer2[256];
  REAL Check;
  char test_byte;

  if (STREAM)
    return;

  // open output-file for systems
  mkdir("Output",S_IRWXU);
  for(i=0;i<NumberOfSystems;i++)
  {
    sprintf(buffer,"Output/System_%d",i);
    mkdir(buffer,S_IRWXU);
  }


  OutputFilePtr=(FILE**)calloc(NumberOfSystems,sizeof(FILE*));
  FILE_CONTENTS = (char**)malloc(NumberOfSystems * sizeof(char*));
  FILE_SIZES=(size_t*)malloc(NumberOfSystems * sizeof(size_t));

  for(i=0;i<NumberOfSystems;i++)
  {
    // read the current position into the output file (in bytes)
    fread(&pos,1,sizeof(fpos_t),FilePtr);


    // construct outputfilename
    sprintf(buffer,"Output/System_%d/output_%s_%d.%d.%d_%lf_%lg%s",
            i,
            Framework[i].Name[0],
            NumberOfUnitCells[i].x,
            NumberOfUnitCells[i].y,
            NumberOfUnitCells[i].z,
            (double)therm_baro_stats.ExternalTemperature[i],
            (double)(therm_baro_stats.ExternalPressure[i][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
            FileNameAppend);

    // limit length of file-name
    strncpy(buffer2,buffer,250);
    sprintf(buffer2,"%s.data",buffer2);

    // check if the file exist
    if( access(buffer2,F_OK )==0) 
    {
       // output-file exist
    #if defined (__APPLE__)
      if(((OutputFilePtr[i]=fopen(buffer2,"r+"))!=NULL)||(pos<0))
    #else
      if(((OutputFilePtr[i]=fopen(buffer2,"r+"))!=NULL)||(pos.__pos<0))
    #endif
      {
        // get the size of file in bytes
        fseek(OutputFilePtr[i], 0L, SEEK_END);
        long sz = ftell(OutputFilePtr[i]);

        // rewind to beginning
        fseek(OutputFilePtr[i], 0L, SEEK_SET);

        // try to reposition to the saved state
        fsetpos(OutputFilePtr[i],&pos);

        if(sz>0)
        {
          fseek(OutputFilePtr[i],-1L, SEEK_CUR);              // move the position 1 byte back
          fread(&test_byte,1, sizeof(char),OutputFilePtr[i]); // do an fread() of 1 byte to check if beyond the file's end
        }
        else
          clearerr(OutputFilePtr[i]); // file-size 0
        
      #if defined (__APPLE__)
        if(feof(OutputFilePtr[i])||(pos>sz))
      #else
        if(feof(OutputFilePtr[i])||(pos.__pos>sz))
      #endif
        {
          // the position is beyond the file' end
        #if defined (__APPLE__)
          fprintf(stderr,"Failed to Reposition output-file at %lld ( beyond file-length %ld)\n",pos,sz);
        #else
          fprintf(stderr,"Failed to Reposition output-file at %ld ( beyond file-length %ld)\n",pos.__pos,sz);
        #endif
          fclose(OutputFilePtr[i]);                          // close file
          OutputFilePtr[i]=fopen(buffer2,"w");      // create new file by reopening as "w"
          PrintPreSimulationStatusCurrentSystem(i); // print the pre-simulation status again
        }
        else
        {
          fsetpos(OutputFilePtr[i],&pos);
        #if defined (__APPLE__)
          fprintf(stderr,"Succesfully repositioned output-file of size %ld to position %lld\n",sz,pos);
        #else
          fprintf(stderr,"Succesfully repositioned output-file of size %ld to position %ld\n",sz,pos.__pos);
        #endif
        }
      }
      else
      {
        fprintf(stderr,"Failed to reposition output-file\n");
        fclose(FilePtr);
        OutputFilePtr[i]=fopen(buffer2,"w");        // create new file by reopening as "w"
        PrintPreSimulationStatusCurrentSystem(i);   // print the pre-simulation status again
      }
    }
    else // the file does not exist
    {
      fprintf(stderr,"File does not exist, recreating output-file\n");
      OutputFilePtr[i]=fopen(buffer2,"w");      // create new file
      PrintPreSimulationStatusCurrentSystem(i); // print the pre-simulation status again
    }
  }
  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartOutput)\n");
    ContinueAfterCrash=FALSE;
  }
}

void WriteRestartOutput(FILE* FilePtr)
{
  int i;
  fpos_t pos;
  REAL Check;

  for(i=0;i<NumberOfSystems;i++)
  {
    fflush(OutputFilePtr[i]);
    fgetpos(OutputFilePtr[i],&pos);
    fwrite(&pos,1,sizeof(fpos_t),FilePtr);
  }

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}
