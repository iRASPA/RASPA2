#!/usr/bin/perl -w
use File::Copy;
use File::Path;

# framework-information
@Forcefield=("GarciaPerez2012");
@CutOffVDW=(12.0);

# temperature-information
@temperature = (140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460); # list of temperatures
@molecule = ("heptane"); # list of molecules
@idealgas = ([1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0]); # list of IG Rosenbluth weights for each temperature

# simulation-information
$SimulationType="MonteCarlo";
$NumberOfCycles="100000";
$NumberOfInitializationCycles="50000";
$PrintEvery="5000";
$RestartFile="no";

# system and queuing information
$divide_into_batches="no"; # combine serial run in larger blocks
$batches = 8; # combine into an 8-core job
$queue = "serial,jaipur,bio-serial"; # the queue type
$job_name = "MC_GIBBS-pentane"; # name of the job
@file_list = (); # list of files copied to all the directories

# get cluster-name
chomp($cluster = `hostname -s`);

# empty 'submit-file'
open(DATw4, ">submit") || die("Could not open file!");
close(DATw4);

$index_batches=0;
$count=0;

$index_molecule=0;
foreach (@molecule)
{
  $index_temperature=0;
  foreach (@temperature)
  {
    $dir_press="$molecule[$index_molecule]/$temperature[$index_temperature]";
    mkpath($dir_press);

    open(DATw1, ">$dir_press/run") || die("Could not open file!");
    printf DATw1 "#! /bin/sh -f\n";
    printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
    printf DATw1 "\${RASPA_DIR}/bin/simulate \$1\n";
    close(DATw1);
    chmod 0755, "$dir_press/run";

    for($i = 0; $i < scalar @file_list; $i++)
    {
      copy("$file_list[$i]","$dir_press/$file_list[$i]") or die "Copy failed: $!";
    }

    if(($cluster eq "login3") || ($cluster eq "login4"))
    {
      open(DATw1, ">$dir_press/bsub.job") || die("Could not open file!");
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
      printf DATw1 "echo \$PBS_JOBID > jobid\n";
      printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
      printf DATw1 "\${RASPA_DIR}/bin//simulate \$1\n";
      close(DATw1);
      chmod 0755, "$dir_press/bsub.job";
    }
    elsif($cluster eq "carbon")
    {
      open(DATw1, ">$dir_press/bsub.job") || die("Could not open file!");
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
      chmod 0755, "$dir_press/bsub.job";
    }
    elsif($cluster eq "kraken")
    {
      open(DATw1, ">$dir_press/bsub.job") || die("Could not open file!");
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
      chmod 0755, "$dir_press/bsub.job";
    }
    else
    {
      die("Unknown cluster!");
    }

    open(DATw3, ">$dir_press/simulation.input") || die("Could not open file!");
    print DATw3 "SimulationType                $SimulationType\n";
    print DATw3 "NumberOfCycles                $NumberOfCycles\n";
    print DATw3 "NumberOfInitializationCycles  $NumberOfInitializationCycles\n";
    print DATw3 "PrintEvery                    $PrintEvery\n";
    print DATw3 "RestartFile                   $RestartFile\n";
    print DATw3 "\n";
    print DATw3 "ContinueAfterCrash no\n";
    print DATw3 "WriteBinaryRestartFileEvery 5000\n";
    print DATw3 "\n";
    print DATw3 "ChargeMethod                  Ewald\n";
    print DATw3 "Forcefield                    $Forcefield[$index_molecule]\n";
    print DATw3 "CutOffVDW                     $CutOffVDW[$index_molecule]\n";
    print DATw3 "\n";
    print DATw3 "Box 0\n";
    print DATw3 "BoxLengths 35 35 35\n";
    print DATw3 "BoxAngles 90 90 90\n";
    print DATw3 "ExternalTemperature $temperature[$index_temperature]\n";
    print DATw3 "\n";
    print DATw3 "Box 1\n";
    print DATw3 "BoxLengths 70 70 70\n";
    print DATw3 "BoxAngles 90 90 90\n";
    print DATw3 "ExternalTemperature $temperature[$index_temperature]\n";
    print DATw3 "\n";
    print DATw3 "GibbsVolumeChangeProbability 1.0\n";
    print DATw3 "\n";

    print DATw3 "Component 0 MoleculeName              $molecule[$index_molecule]\n";
    print DATw3 "            StartingBead              0\n";
    print DATw3 "            MoleculeDefinition        TraPPE\n";
    print DATw3 "            IdealGasRosenbluthWeight  $idealgas[$index_temperature][$index_molecule]\n";
    print DATw3 "            FugacityCoefficient       1.0\n";
    print DATw3 "            TranslationProbability    1.0\n";
    print DATw3 "            RotationProbability       1.0\n";
    print DATw3 "            ReinsertionProbability    1.0\n";
    print DATw3 "            CBMCProbability           1.0\n";
    print DATw3 "            GibbsSwapProbability      1.0\n";
    print DATw3 "            CreateNumberOfMolecules   200 100\n";
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
          printf DATw4 "\n";
          close(DATw4);
          chmod 0755, "$dir_press/bsub.job";
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
          printf DATw4 "\n";
          printf DATw4 "# set path for intel compiler libraries\n";
          printf DATw4 "source /opt/intel/cce/10.1.018/bin/iccvars.sh\n";
          printf DATw4 "\n";
          close(DATw4);
          chmod 0755, "$dir_press/bsub.job";
        }
        else
        {
          die("Unknown cluster!");
        }

        open(DATw5, ">>submit") || die("Could not open file!");
        print DATw5 "qsub ./submit_$index_batches\n";
        close(DATw5);
      }

      open(DATw5, ">>submit_$index_batches") || die("Could not open file!");
      print DATw5 "cd \"$dir_press\"\n";
      print DATw5 ".\/run > out 2\>\&1\&\n";
      print DATw5 "cd -\n";

      if(($count+1)%$batches==0) {print DATw5 "wait\n";}
      close(DATw5);

      $count=$count+1;
    }
    else
    {
      open(DATw4, ">>submit") || die("Could not open file!");
      print DATw4 "cd \"$dir_press\"\n";
      print DATw4 "qsub bsub.job\n";
      print DATw4 "cd -\n";
      close(DATw4);
    }
    $index_temperature=$index_temperature+1;
  }
  $index_molecule=$index_molecule+1;
}

# put a "wait" on the last submit which (probably) is not a full batch
if($divide_into_batches eq "yes")
{
  open(DATw5, ">>submit_$index_batches") || die("Could not open file!");
  if(($count)%$batches!=0) {print DATw5 "wait\n";}
  close(DATw5);
}

chmod 0755, "submit"

