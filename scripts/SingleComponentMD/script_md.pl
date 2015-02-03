#!/usr/bin/perl -w
use File::Copy;
use File::Path;

# framework-information
@framework = ("ZnBDC");  # list of structures
@HeliumVoidFraction = ("0.516026");  # list of helium-voidfractions for the structures
@UnitCells = ("4 2 2"); # list of sizes of the unit cells for the structures
@Forcefield=("GenericMOFs");

# temperature-information
@temperature = (313.0); # list of temperatures
@molecule = ("hexane","3-methylpentane","22-dimethylbutane"); # list of molecules
@mol_fraction = (1.0,1.0,1.0);

# loading-information
$loading_start = 1;
$loading_end = 80;
$number_of_loading_points = 32; # 32 point equally spaced

# simulation-information
$SimulationType="MolecularDynamics";
$NumberOfCycles="100000000000";
$NumberOfInitializationCycles="10000";
$NumberOfEquilibrationCycles="100000";
$PrintEvery="10000";
$RestartFile="no";

# system and queuing information
$divide_into_batches="yes"; # combine serial run in larger blocks
$batches = 8; # combine into an 8-core job
$queue = "mof1"; # the queue type
$job_name = "MD"; # name of the job

for($i=0;$i<$number_of_loading_points;$i++)
{
  $loading[$i]=($loading_end-$loading_start)*($i/($number_of_loading_points-1))+$loading_start;
}

# get cluster-name
chomp($cluster = `hostname -s`);

# empty 'submit-file'
open(DATw4, ">submit") || die("Could not open file!");
close(DATw4);

$index_batches=0;
$count=0;

$index_framework=0;
foreach (@framework)
{
  $dir_fr=$framework[$index_framework];
  mkpath($dir_fr);

  $index_temperature=0;
  foreach (@temperature)
  {
    $dir_temp="$dir_fr/$temperature[$index_temperature]K";
    mkpath($dir_temp);

    $index_molecule=0;
    foreach (@molecule)
    {
      $dir_mol="$dir_temp/$molecule[$index_molecule]";
      mkpath($dir_mol);

      $index_loading=0;
      foreach (@loading)
      {
        $dir_loading="$dir_mol/$loading[$index_loading]_mol_per_cell";
        mkpath($dir_loading);

        open(DATw1, ">$dir_loading/run") || die("Could not open file!");
        printf DATw1 "#! /bin/sh -f\n";
        printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
        printf DATw1 "\${RASPA_DIR}/bin/simulate \$1";
        close(DATw1);
        chmod 0755, "$dir_loading/run";

        if(($cluster eq "login3") || ($cluster eq "login4"))
        {
          open(DATw1, ">$dir_loading/bsub.job") || die("Could not open file!");
          printf DATw1 "#!/bin/bash\n";
          printf DATw1 "#PBS -S /bin/bash\n";
          printf DATw1 "#PBS -l nodes=1\n";
          printf DATw1 "#PBS -N $job_name\n";
          printf DATw1 "#PBS -l walltime=120:00:00\n";
          printf DATw1 "#PBS -o pbs.out\n";
          printf DATw1 "#PBS -e pbs.err\n";
          printf DATw1 "#PBS -r n\n";
          printf DATw1 "#PBS -V\n";
          printf DATw1 "#PBS -mba\n";
          printf DATw1 "\n";
          printf DATw1 "cd \${PBS_O_WORKDIR}\n";
          printf DATw1 "\n";
          printf DATw1 "echo \$JOB_ID > jobid\n";
          printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
          printf DATw1 "\${RASPA_DIR}/bin//simulate \$1\n";
          close(DATw1);
          chmod 0755, "$dir_loading/bsub.job";
        }
        elsif($cluster eq "carbon")
        {
          open(DATw1, ">$dir_loading/bsub.job") || die("Could not open file!");
          printf DATw1 "#!/bin/bash\n";
          printf DATw1 "# script for Grid Engine\n";
          printf DATw1 "#\$ -S /bin/bash\n";
          printf DATw1 "#\$ -q $queue\n";
          printf DATw1 "#\$ -N $job_name\n";
          printf DATw1 "#\$ -V\n";
          printf DATw1 "#\$ -cwd\n";
          printf DATw1 "#\$ -notify\n";
          printf DATw1 "\n";
          printf DATw1 "echo \$JOB_ID > jobid\n";
          printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
          printf DATw1 "\${RASPA_DIR}/bin/simulate \$1\n";
          close(DATw1);
          chmod 0755, "$dir_loading/bsub.job";
        }
        elsif($cluster eq "kraken")
        {
          open(DATw1, ">$dir_loading/bsub.job") || die("Could not open file!");
          printf DATw1 "#!/bin/bash\n";
          printf DATw1 "# script for Grid Engine\n";
          printf DATw1 "#\$ -S /bin/bash\n";
          printf DATw1 "#\$ -q $queue\n";
          printf DATw1 "#\$ -N $job_name\n";
          printf DATw1 "#\$ -V\n";
          printf DATw1 "#\$ -cwd\n";
          printf DATw1 "#\$ -notify\n";
          printf DATw1 "\n";
          printf DATw1 "# set path for intel compiler libraries\n";
          printf DATw1 "source /opt/intel/cce/10.1.018/bin/iccvars.sh\n";
          printf DATw1 "\n";
          printf DATw1 "echo \$JOB_ID > jobid\n";
          printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
          printf DATw1 "\${RASPA_DIR}/bin/simulate \$1\n";
          close(DATw1);
          chmod 0755, "$dir_loading/bsub.job";
        }

        open(DATw3, ">$dir_loading/simulation.input") || die("Could not open file!");
        print DATw3 "SimulationType                $SimulationType\n";
        print DATw3 "NumberOfCycles                $NumberOfCycles\n";
        print DATw3 "NumberOfInitializationCycles  $NumberOfInitializationCycles\n";
        print DATw3 "NumberOfEquilibrationCycles   $NumberOfEquilibrationCycles\n";
        print DATw3 "PrintEvery                    $PrintEvery\n";
        print DATw3 "RestartFile                   $RestartFile\n";
        print DATw3 "\n";
        print DATw3 "ContinueAfterCrash no\n";
        print DATw3 "WriteBinaryRestartFileEvery 5000\n";
        print DATw3 "\n";
        print DATw3 "Ensemble                      NVT\n";
        print DATw3 "\n";
        print DATw3 "ChargeMethod                  Ewald\n";
        print DATw3 "Forcefield                    $Forcefield[$index_framework]\n";
        print DATw3 "EwaldPrecision                1e-6\n";
        print DATw3 "TimeStep                      0.0005\n";
        print DATw3 "\n";
        print DATw3 "Framework           0\n";
        print DATw3 "FrameworkName       $framework[$index_framework]\n";
        print DATw3 "UnitCells           $UnitCells[$index_framework]\n";
        print DATw3 "HeliumVoidFraction  $HeliumVoidFraction[$index_framework]\n";
        print DATw3 "ExternalTemperature $temperature[$index_temperature]\n";
        print DATw3 "ExternalPressure    0\n";
        print DATw3 "ComputeMSD          yes\n";
        print DATw3 "WriteMSDEvery       5000\n";
        print DATw3 "\n";

        print DATw3 "\n";
        print DATw3 "Component 0 MoleculeName              $molecule[$index_molecule]\n";
        print DATw3 "            StartingBead              0\n";
        print DATw3 "            MoleculeDefinition        TraPPE\n";
        print DATw3 "            FugacityCoefficient       1.0\n";
        print DATw3 "            TranslationProbability    1.0\n";
        print DATw3 "            RotationProbability       1.0\n";
        print DATw3 "            ReinsertionProbability    1.0\n";
        print DATw3 "            CBMCProbability           1.0\n";
        print DATw3 "            SwapProbability           0.0\n";
        printf DATw3 "            CreateNumberOfMolecules   %d\n",$loading[$index_loading]*$mol_fraction[$index_molecule];
        close(DATw3);

        if($divide_into_batches eq "yes")
        {
          if($count%$batches==0)
          {
            $index_batches=$index_batches+1;

            if(($cluster eq "login3") || ($cluster eq "login4"))
            {
              open(DATw4, ">submit_$index_batches") || die("Could not open file!");
              print DATw4 "\#\!/bin/bash\n";
              print DATw4 "#PBS -S /bin/bash\n";
              print DATw4 "#PBS -l nodes=1\n";
              print DATw4 "#PBS -N $job_name\n";
              print DATw4 "#PBS -l walltime=120:00:00\n";
              print DATw4 "#PBS -o pbs.out\n";
              print DATw4 "#PBS -e pbs.err\n";
              print DATw4 "#PBS -r n\n";
              print DATw4 "#PBS -V\n";
              print DATw4 "#PBS -mba\n";
              print DATw4 "\n";
              print DATw4 "cd \$\{PBS_O_WORKDIR\}\n";
              print DATw4 "\n";
              close(DATw4);
              chmod 0755, "submit_$index_batches"
            }
            elsif($cluster eq "carbon")
            {
              open(DATw4, ">submit_$index_batches") || die("Could not open file!");
              printf DATw4 "#!/bin/bash\n";
              printf DATw4 "# script for Grid Engine\n";
              printf DATw4 "#\$ -S /bin/bash\n";
              printf DATw4 "#\$ -q $queue\n";
              printf DATw4 "#\$ -N $job_name\n";
              printf DATw4 "#\$ -V\n";
              printf DATw4 "#\$ -cwd\n";
              printf DATw4 "#\$ -pe orte $batches\n";
              printf DATw4 "#\$ -notify\n";
              printf DATw4 "\n";
              close(DATw4);
              chmod 0755, "$dir_loading/bsub.job";
            }
            elsif($cluster eq "kraken")
            {
              open(DATw4, ">submit_$index_batches") || die("Could not open file!");
              printf DATw4 "#!/bin/bash\n";
              printf DATw4 "# script for Grid Engine\n";
              printf DATw4 "#\$ -S /bin/bash\n";
              printf DATw4 "#\$ -q $queue\n";
              printf DATw4 "#\$ -N $job_name\n";
              printf DATw4 "#\$ -V\n";
              printf DATw4 "#\$ -cwd\n";
              printf DATw4 "#\$ -pe orte $batches\n";
              printf DATw4 "#\$ -notify\n";
              printf DATw4 "\n";
              printf DATw4 "# set path for intel compiler libraries\n";
              printf DATw4 "source /opt/intel/cce/10.1.018/bin/iccvars.sh\n";
              printf DATw4 "\n";
              close(DATw4);
              chmod 0755, "$dir_loading/bsub.job";
            }

            open(DATw5, ">>submit") || die("Could not open file!");
            print DATw5 "qsub ./submit_$index_batches\n";
            close(DATw5);
          }

          open(DATw5, ">>submit_$index_batches") || die("Could not open file!");
          print DATw5 "cd \"$dir_loading\"\n";
          print DATw5 ".\/run > out 2\>\&1\&\n";
          print DATw5 "cd -\n";

          if(($count+1)%$batches==0) {print DATw5 "wait\n";}
          close(DATw5);

          $count=$count+1;
        }
        else
        {
          open(DATw4, ">>submit") || die("Could not open file!");
          print DATw4 "cd \"$dir_loading\"\n";
          print DATw4 "qsub bsub.job\n";
          print DATw4 "cd -\n";
          close(DATw4);
        }

        $index_loading=$index_loading+1;
      }
      $index_molecule=$index_molecule+1;
    }
    $index_temperature=$index_temperature+1;
  }
  $index_framework=$index_framework+1;
}

# put a "wait" on the last submit which (probably) is not a full batch
if($divide_into_batches eq "yes")
{
  open(DATw5, ">>submit_$index_batches") || die("Could not open file!");
  if(($count)%$batches!=0) {print DATw5 "wait\n";}
  close(DATw5);
}

chmod 0755, "submit"

