#!/usr/bin/perl -w
use File::Copy;
use File::Path;
use File::Find;

# list of frameworks
@framework = ("fe2bdp3-minusb");

# temperature-information
@temperature = (433.0); # list of temperatures
@molecule = ("hexane","2-methylpentane","3-methylpentane","22-dimethylbutane","23-dimethylbutane","heptane","2-methylhexane","3-methylhexane","22-dimethylpentane","23-dimethylpentane"); # list of molecules

# create the results directory
mkpath("Results");

# loop over frameworks, temperatures, molecules
$index_framework=0;
foreach (@framework)
{
  $index_temperature=0;
  foreach (@temperature)
  {
    $index_molecule=0;
    foreach (@molecule)
    {
      $dir_molecule="$framework[$index_framework]/$temperature[$index_temperature]K/$molecule[$index_molecule]";

      open(DATw2, ">data");

      #search for all RASPA output-files
      @files=();
      find (\&search, "$dir_molecule");
      sub search {push (@files,$File::Find::name) if(-f and /output*/);}

      foreach $file (@files)
      {
        $fugacity=`gawk '/Partial fugacity:/{++c;if(c==1) printf \$3" "}' '$file'`;
        $fugacitycoefficient=`gawk '/Fugacity coefficient:/{++c;if(c==1) printf \$3" "}' '$file'`;
        $press=`gawk '/Partial pressure:/{++c;if(c==1) printf \$3" "}' '$file'`;
        $absloading=`gawk '/Average loading absolute \\[molecules/{printf pr1" "pr2}{pr1=\$2;pr2=\$4" "}' '$file'`;
        $absloading_molec_uc=`gawk '/Average loading absolute \\[molecules/{++c;if(c==1) printf \$6" "\$8" "}' '$file'`;
        $absloading_mol_kg=`gawk '/Average loading absolute \\[mol\\/kg/{printf \$6" "\$8" "}' '$file'`;
        $absloading_mg_g=`gawk '/Average loading absolute \\[mil/{++c;if(c==1) printf \$6" "\$8" "}' '$file'`;
        $excloading=`gawk '/Average loading excess \\[molecules/{printf pr1" "pr2}{pr1=\$2;pr2=\$4}' '$file'`;
        $excloading_molec_uc=`gawk '/Average loading excess \\[molecules/{++c;if(c==1) printf \$6" "\$8" "}' '$file'`;
        $excloading_mol_kg=`gawk '/Average loading excess \\[mol\\/kg/{++c;if(c==1) printf \$6" "\$8" "}' '$file'`;
        $excloading_mg_g=`gawk '/Average loading excess \\[mil/{++c;if(c==1) printf \$6" "\$8" "}' '$file'`;
        $heats=`gawk '/Heat of desorption:/{getline; getline; getline; getline; getline; getline; getline; getline; ++c;if(c==1) print \$2" "\$4" "}' '$file'`;

        printf DATw2 "$fugacity $fugacitycoefficient $press $absloading $absloading_molec_uc $absloading_mol_kg $absloading_mg_g $excloading $excloading_molec_uc $excloading_mol_kg $excloading_mg_g $heats\n";
      } 
      close(DATw2);

      `sort -n data > data_sorted`;
      `sed '/^\$/d' data_sorted > data_sorted_no_empty_lines`;

      open(DATw2, ">Results/Results.dat-$framework[$index_framework]-$temperature[$index_temperature]K-$molecule[$index_molecule]") || die("Could not open file!");
      printf DATw2 "# $framework[$index_framework], $temperature[$index_temperature]K, $molecule[$index_molecule]\n";
      printf DATw2 "#  1: pressure [Pa]\n";
      printf DATw2 "#  2: fugacity coefficient [-]\n";
      printf DATw2 "#  3: fugacity [Pa]\n";
      printf DATw2 "#  4: (absolute) loading [molecules per total cell]\n";
      printf DATw2 "#  5: (absolute) error [molecules per total cell]\n";
      printf DATw2 "#  6: (absolute) loading [molecules per unit cell]\n";
      printf DATw2 "#  7: (absolute) error [molecules per unit cell]\n";
      printf DATw2 "#  8: (absolute) loading [mol/kg]\n";
      printf DATw2 "#  9: (absolute) error [mol/kg]\n";
      printf DATw2 "# 10: (absolute) loading [mg/g]\n";
      printf DATw2 "# 11: (absolute) error [mg/g]\n";
      printf DATw2 "# 12: (excess) loading [molecules per total cell]\n";
      printf DATw2 "# 13: (excess) error [molecules per total cell]\n";
      printf DATw2 "# 14: (excess) loading [molecules per unit cell]\n";
      printf DATw2 "# 15: (excess) error [molecules per unit cell]\n";
      printf DATw2 "# 16: (excess) loading [mol/kg]\n";
      printf DATw2 "# 17: (excess) error [mol/kg]\n";
      printf DATw2 "# 18: (excess) loading [mg/g]\n";
      printf DATw2 "# 19: (excess) error [mg/g]\n";
      printf DATw2 "# 20: heat of desorption [K]\n";
      printf DATw2 "# 21: heat of desorption error [K]\n";
      printf DATw2 "\n";
      close(DATw2);

      `cat data_sorted_no_empty_lines >> 'Results/Results.dat-$framework[$index_framework]-$temperature[$index_temperature]K-$molecule[$index_molecule]'`;

      $index_molecule=$index_molecule+1;
    }
    $index_temperature=$index_temperature+1;
  }
  $index_framework=$index_framework+1;
}

unlink("data","data_sorted","data_sorted_no_empty_lines");
