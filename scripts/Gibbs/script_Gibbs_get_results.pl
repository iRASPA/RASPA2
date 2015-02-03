#!/usr/bin/perl -w
use File::Copy;
use File::Path;
use File::Find;


# temperature-information
@temperature = (100,110,120,130,140,150,160,165,170,175); # list of temperatures
@molecule = ("methane"); # list of molecules

# create the results directory
mkpath("Results");

open(DATw1, ">Results/Results.dat") || die("Could not open file!");

# loop over frameworks, temperatures, molecules
$index_molecule=0;
foreach (@molecule)
{
  $index_temperature=0;
  foreach (@temperature)
  {
    $dir_molecule="$molecule[$index_molecule]/$temperature[$index_temperature]";

    open(DATw2, ">Results/Results.dat-$molecule[$index_molecule]-$temperature[$index_temperature]K") || die("Could not open file!");

    #search for all RASPA output-files
    @files=();
    find (\&search, "$dir_molecule");
    sub search {push (@files,$File::Find::name) if(-f and /output*/);}

    foreach $file (@files)
    {
      $density=`gawk '/Average Density:/{getline; getline; getline; getline; getline; getline; getline; getline; ++c;if(c==1) print \$2" "\$5" "}' '$file'`;

      printf DATw1 "$temperature[$index_temperature] $density";
      printf DATw2 "$density\n";
    } 
    close(DATw2);

    $index_temperature=$index_temperature+1;
  }
  $index_molecule=$index_molecule+1;
}
close(DATw1);

