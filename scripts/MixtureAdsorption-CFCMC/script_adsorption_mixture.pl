#!/usr/bin/perl -w
use File::Copy;
use File::Path;

# framework-information
@Framework = ("ITQ-29");  # list of structures
@BlockPockets = ("yes");
@BlockPocketsFilename = ("ITQ-29");  # list of blocking pockets filenames
@HeliumVoidFraction = ("0.401506");  # list of helium-voidfractions for the structures
@UnitCells = ("2 2 2"); # list of sizes of the unit cells for the structures
@Forcefield=("ExampleZeolitesForceField");
@UseChargesFromCIFFile=("yes");
@RemoveAtomNumberCodeFromLabel=("yes");
@CutOff=("11.8");

# temperature-information
@Temperature = (300.0); # list of temperatures
@Molecule = ("CO2","methane"); # list of molecules
@MoleculeDefinition = ("ExampleDefinitions","ExampleDefinitions");
@Idealgas = ([1.0,1.0]); # list of IG Rosenbluth weights for each temperature

# pressure-information
$pressure_type = "fugacity";     # "pressure" or "fugacity"
$pressure_start = 1e2;          # lowest pressure
$pressure_end = 1e16;            # highest pressure
$number_of_pressure_points = 36; # number of points equally spaced in log-scale or linear scale
$pressure_unit="Pa";             # "bar", "kPa", or "Pa"
$pressure_scale="log";           # "log" or "linear"

# simulation-information
$SimulationType="MonteCarlo";
$NumberOfCycles="500000";
$NumberOfInitializationCycles="100000";
$NumberOfEquilibrationCycles="100000";
$PrintEvery="5000";
$RestartFile="no";

# system and queuing information
$queueing_system="gridengine";  # "gridengine" or "torque"
$divide_into_batches="no"; # combine serial run in larger blocks
$batches = 8; # combine into an 8-core job
$queue = "mof4,mof5,VASP"; # the queue type
$job_name = "Mixture"; # name of the job
@file_list = (); # list of files copied to all the directories

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

if($pressure_scale eq "log")
{
  for($i = 0; $i < $number_of_pressure_points; $i++)
  {
    $pressure[$i]=10**(((log10($pressure_end)-log10($pressure_start))*($i/($number_of_pressure_points-1)))+log10($pressure_start));
  }
}
else
{
  for($i = 0; $i < $number_of_pressure_points; $i++)
  {
    $pressure[$i]=($pressure_end-$pressure_start)*($i/($number_of_pressure_points-1))+$pressure_start;
  }
}

# empty 'submit-file'
open(DATw4, ">submit") || die("Could not open file!");
close(DATw4);

$index_batches=0;
$count=0;

$index_framework=0;
foreach (@Framework)
{
  $dir_fr=$Framework[$index_framework];
  mkpath($dir_fr);

  $index_temperature=0;
  foreach (@Temperature)
  {
    # concat the molecule names with "-" between them
    $molecule_list = join("-",@Molecule);

    $dir_temp="$dir_fr/$Temperature[$index_temperature]K/$molecule_list";
    mkpath($dir_temp);

    $index_pressure=0;
    foreach (@pressure)
    {
      $dir_press="$dir_temp/$pressure[$index_pressure]$pressure_unit";
      mkpath($dir_press);

      open(DATw1, ">$dir_press/run") || die("Could not open file!");
      printf DATw1 "#! /bin/sh -f\n";
      printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
      printf DATw1 "\${RASPA_DIR}/bin/simulate \$1";
      close(DATw1);
      chmod 0755, "$dir_press/run";

      # copy files to the dir
      for($i = 0; $i < scalar @file_list; $i++)
      {
        copy("$file_list[$i]","$dir_press/$file_list[$i]") or die "Copy failed: $!";
      }

      if($queueing_system eq "torque")
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
        printf DATw1 "echo \$JOB_ID > jobid\n";
        printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
        printf DATw1 "\${RASPA_DIR}/bin//simulate \$1\n";
        close(DATw1);
        chmod 0755, "$dir_press/bsub.job";
      }
      elsif($queueing_system eq "gridengine")
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
      else
      {
        die("Unknown queuing system!");
      }

      
      open(DATw3, ">$dir_press/simulation.input") || die("Could not open file!");
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
      print DATw3 "ChargeMethod                  Ewald\n";
      print DATw3 "Forcefield                    $Forcefield[$index_framework]\n";
      print DATw3 "CutOffVDW                     $CutOff[$index_framework]\n";
      print DATw3 "RemoveAtomNumberCodeFromLabel $RemoveAtomNumberCodeFromLabel[$index_framework]\n";
      print DATw3 "\n";
      print DATw3 "Framework             0\n";
      print DATw3 "FrameworkName         $Framework[$index_framework]\n";
      print DATw3 "UseChargesFromCIFFile $UseChargesFromCIFFile[$index_framework]\n";
      print DATw3 "UnitCells             $UnitCells[$index_framework]\n";
      print DATw3 "HeliumVoidFraction    $HeliumVoidFraction[$index_framework]\n";
      print DATw3 "ExternalTemperature   $Temperature[$index_temperature]\n";
      if($pressure_unit eq "bar") {printf DATw3 "ExternalPressure      %g\n", $pressure[$index_pressure]*1e5}
      elsif($pressure_unit eq "kPa") {printf DATw3 "ExternalPressure      %g\n", $pressure[$index_pressure]*1e3}
      elsif($pressure_unit eq "Pa")  {printf DATw3 "ExternalPressure      %g\n", $pressure[$index_pressure]}
      print DATw3 "\n";

      $index_molecule=0;
      foreach (@Molecule)
      {
        print DATw3 "\n";
        print DATw3 "Component $index_molecule MoleculeName              $Molecule[$index_molecule]\n";
        print DATw3 "            StartingBead              0\n";
        print DATw3 "            MoleculeDefinition        $MoleculeDefinition[$index_molecule]\n";
        print DATw3 "            BlockPockets              $BlockPockets[$index_framework]\n";
        print DATw3 "            BlockPocketsFilename      $BlockPocketsFilename[$index_framework]\n";
        print DATw3 "            IdealGasRosenbluthWeight  $Idealgas[$index_temperature][$index_molecule]\n";
        if($pressure_type eq "fugacity") {printf DATw3 "            FugacityCoefficient       1.0\n";}
        print DATw3 "            TranslationProbability    1.0\n";
        print DATw3 "            RotationProbability       1.0\n";
        print DATw3 "            ReinsertionProbability    1.0\n";
        print DATw3 "            CBMCProbability           1.0\n";
        print DATw3 "            IdentityChangeProbability 1.0\n";
        print DATw3 "              NumberOfIdentityChanges ",scalar @Molecule,"\n";
        print DATw3 "              IdentityChangesList    ";
        for($i = 0; $i < scalar @Molecule; $i++) {print DATw3 " $i"};
        print DATw3 "\n";
        print DATw3 "            CFSwapLambdaProbability   1.0\n";
        print DATw3 "            CreateNumberOfMolecules   0\n";
        $index_molecule=$index_molecule+1;
      }
      close(DATw3);

      if($divide_into_batches eq "yes")
      {
        if($count%$batches==0)
        {
          $index_batches=$index_batches+1;

          if($queueing_system eq "torque")
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
          elsif($queueing_system eq "gridengine")
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
            chmod 0755, "$dir_press/bsub.job";
          }
          else
          {
            die("Unknown queueing system!");
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


      $index_pressure=$index_pressure+1;
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

