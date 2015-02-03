#!/usr/bin/perl -w
use File::Copy;
use File::Path;
use File::Find;

# temperature-information
@temperature = (120,180,240,280,300,320); # list of temperatures
@molecule = ("CO2"); # list of molecules

@epsilon=([27,79],[28,80]);
@sigma=([2.80,3.05],[2.78,3.0]);

# create the results directory
mkpath("Results");

open(DATw1, ">Results/Results.dat") || die("Could not open file!");

# loop over frameworks, temperatures, molecules
$index_molecule=0;
foreach (@molecule)
{
  $index_sigma=0;
  foreach (@sigma)
  {
    $sigma_list=join("-",@{$sigma[$index_sigma]});
    $index_epsilon=0;
    foreach (@epsilon)
    {
      $epsilon_list=join("-",@{$epsilon[$index_epsilon]});
      open(DATw2, ">Results/Results.dat-$molecule[$index_molecule]-$sigma_list-$epsilon_list") || die("Could not open file!");
      $index_temperature=0;
      foreach (@temperature)
      {
        $dir_molecule="$molecule[$index_molecule]/$temperature[$index_temperature]K/$sigma_list/$epsilon_list";

        #search for all RASPA output-files
        @files=();
        find (\&search, "$dir_molecule");
        sub search {push (@files,$File::Find::name) if(-f and /output*/);}

        foreach $file (@files)
        {
          $density=`gawk '/Average Density:/{getline; getline; getline; getline; getline; getline; getline; getline; ++c;if(c==1) print \$2" "\$5" "}' '$file'`;

          printf DATw1 "$temperature[$index_temperature] $density $sigma_list $epsilon_list";
          printf DATw2 "$temperature[$index_temperature] $density\n";
        } 
        $index_temperature=$index_temperature+1;
      }
      close(DATw2);
      $index_epsilon=$index_epsilon+1;
    }
    $index_sigma=$index_sigma+1;
  }
  $index_molecule=$index_molecule+1;
}
close(DATw1);

