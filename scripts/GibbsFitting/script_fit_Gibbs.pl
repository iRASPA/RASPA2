#!/usr/bin/perl -w
use File::Copy;
use File::Path;

# framework-information
@Forcefield=("GarciaPerez2012");

@epsilon=([27,79],[28,80],[32,85],[40,90],[56,53]);
@sigma=([2.76,3.0],[2.78,30.2],[2.80,3.05],[2.82,3.1],[3.96,3.0]);

@FitAtoms=("C_co2","O_co2");

# temperature-information
@temperature = (120,180,240,280); # list of temperatures
@molecule = ("CO2"); # list of molecules
@idealgas = ([1.0],[1.0],[1.0],[1.0]); # list of IG Rosenbluth weights for each temperature

# simulation-information
$SimulationType="MonteCarlo";
$NumberOfCycles="75000";
$NumberOfInitializationCycles="40000";
$PrintEvery="5000";
$RestartFile="no";

# system and queuing information
$divide_into_batches="yes"; # combine serial run in larger blocks
$batches = 8; # combine into an 8-core job
$queue = "default"; # the queue type
$job_name = "MC_GIBBS_FIT"; # name of the job
@file_list = ("force_field_mixing_rules.def","pseudo_atoms.def"); # list of files copied to all the directories

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
    $index_sigma=0;
    foreach (@sigma)
    {
      $sigma_list=join("-",@{$sigma[$index_sigma]});
      $index_epsilon=0;
      foreach (@epsilon)
      {
        $epsilon_list=join("-",@{$epsilon[$index_epsilon]});
        $dir_press="$molecule[$index_molecule]/$temperature[$index_temperature]K/$sigma_list/$epsilon_list";
        mkpath($dir_press);

        open(DATw1, ">$dir_press/run") || die("Could not open file!");
        printf DATw1 "#! /bin/sh -f\n";
        printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
        printf DATw1 "\${RASPA_DIR}/bin/simulate \$1\n";
        close(DATw1);
        chmod 0755, "$dir_press/run";

        # copy files to the dir
        for($i = 0; $i < scalar @file_list; $i++)
        {
          copy("$file_list[$i]","$dir_press/$file_list[$i]") or die "Copy failed: $!";
        }

        # read in the forcefield file
        open(FILE,"$dir_press/force_field_mixing_rules.def") || die("Cannot Open File");
        my(@fcont) = <FILE>;
        close FILE;
       
        # replace atom LJ paramters by the trial ones 
        open(FOUT,">$dir_press/force_field_mixing_rules.def") || die("Cannot Open File");
        foreach $line (@fcont)
        {
          for($i = 0; $i < scalar @FitAtoms; $i++)
          {
            if($line =~ /$FitAtoms[$i]/)
            {
                $line = "$FitAtoms[$i]        Lennard-jones     $epsilon[$index_epsilon][$i]       $sigma[$index_sigma][$i]\n";
            }
          }
          print FOUT $line;
        }
        close FOUT;

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
          #die("Unknown cluster!");
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
        print DATw3 "\n";
        print DATw3 "Box 0\n";
        print DATw3 "BoxLengths 30 30 30\n";
        print DATw3 "BoxAngles 90 90 90\n";
        print DATw3 "ExternalTemperature $temperature[$index_temperature]\n";
        print DATw3 "\n";
        print DATw3 "Box 1\n";
        print DATw3 "BoxLengths 50 50 50\n";
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
        print DATw3 "            CreateNumberOfMolecules   400 100\n";
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
              #die("Unknown cluster!");
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
        $index_epsilon=$index_epsilon+1;
      }
      $index_sigma=$index_sigma+1;
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

